/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
 * Department of Computer Science.                                                 *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

#include "parallelcoordinatesaxescreator.h"

#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"

#include "modules/ensembleanalysis/utils/ensemblehash.h"
#include "modules/ensembleanalysis/utils/utils.h"

#include <random>


namespace voreen {

ParallelCoordinatesAxesCreator::ParallelCoordinatesAxesCreator()
    : AsyncComputeProcessor<ParallelCoordinatesAxesCreatorInput, ParallelCoordinatesAxesCreatorOutput>()
    , ensembleport_(Port::INPORT, "port_ensemble", "Ensemble Input", false, Processor::VALID )
    , seedMask_(Port::INPORT, "port_volume", "Volume Mask Input", false, Processor::VALID )
    , axesport_(Port::OUTPORT, "port_axes", "Parallel Coordinates Axes" )
    , propertyMembers_("property_members", "Selected Members", Processor::VALID )
    , propertyFields_("property_fields", "Selected Fields", Processor::VALID )
    , propertySpatialSampleCount_("property_spatial_sample_count", "Spatial Sample Count", 32768, 1, 4194304, Processor::VALID )
    , propertyTemporalSampleCount_("property_temporal_sample_count", "Temporal Sample Count", 1, 1, std::numeric_limits<int>::max(), Processor::VALID )
    , propertySampleRegion_("sampleRegion", "Sample Region")
    , propertySeedTime_("property_seed_time", "Current Random Seed", static_cast<int>(time(0)), std::numeric_limits<int>::min(), std::numeric_limits<int>::max())
    , propertyAggregateMembers_("property_aggregate_members", "Aggregate Members", false, Processor::VALID )
    , propertyFileDialog_("property_file_dialog", "File Output", "Select File...", "", "Voreen Parallel Coordinates (*.vpc)", FileDialogProperty::FileMode::SAVE_FILE, Processor::VALID )
    , propertySaveButton_("property_save", "Save File", Processor::VALID )
{
    // --- Initialize Ports --- //
    this->addPort(ensembleport_ );
    ON_CHANGE(ensembleport_, ParallelCoordinatesAxesCreator, adjustToEnsemble);
    this->addPort(seedMask_ );
    this->addPort(axesport_ );

    // --- Initialize Properties --- //
    this->addProperty(propertyMembers_ );
    ON_CHANGE(propertyMembers_, ParallelCoordinatesAxesCreator, adjustToSelection);
    this->addProperty(propertyFields_ );
    this->addProperty(propertySpatialSampleCount_);
    this->addProperty(propertyTemporalSampleCount_ );
    //this->addProperty(propertySampleRegion_); // Only allow ensemble bounds to sync with ParallelCoordinatesVoxelSelection
    propertySampleRegion_.addOption("bounds", "Ensemble Bounds"); // Will be selected.
    propertySampleRegion_.addOption("common", "Common Bounds");
    this->addProperty(propertySeedTime_);
    this->addProperty(propertyAggregateMembers_ );
    this->addProperty(propertyFileDialog_ );
    this->addProperty(propertySaveButton_ );
    propertySaveButton_.onChange(LambdaFunctionCallback([this] {
        if(axesport_.hasData() && !propertyFileDialog_.get().empty() )
            axesport_.getData()->serialize(propertyFileDialog_.get() );
    } ) );
}

Processor* ParallelCoordinatesAxesCreator::create() const {
    return new ParallelCoordinatesAxesCreator();
}

bool ParallelCoordinatesAxesCreator::isReady() const {
    if(!ensembleport_.isReady()) {
        setNotReadyErrorMessage("No input");
        return false;
    }

    // Note: volumeport is optional!

    if(!axesport_.isReady()) {
        setNotReadyErrorMessage("Outport not connected");
        return false;
    }

    return true;
}

void ParallelCoordinatesAxesCreator::adjustToEnsemble() {
    bool hasData = ensembleport_.hasData();
    propertyMembers_.setReadOnlyFlag(!hasData);
    propertyFields_.setReadOnlyFlag(!hasData);

    if(!hasData) {
        return;
    }

    const EnsembleDataset& ensemble = *ensembleport_.getData();

    std::string hash = EnsembleHash(ensemble).getHash();
    if(hash != hash_) {
        propertyMembers_.blockCallbacks(true );
        propertyMembers_.reset();
        for( const auto& member : ensemble.getMembers() ) {
            propertyMembers_.addRow(member.getName(), member.getColor());
        }
        propertyMembers_.blockCallbacks(false );
        propertyMembers_.invalidate();

        hash_ = hash;
    }
}

void ParallelCoordinatesAxesCreator::adjustToSelection() {
    auto ensemble = EnsembleDataset();
    for( const auto member : propertyMembers_.get() )
        ensemble.addMember(ensembleport_.getData()->getMembers()[member]);

    auto fields = std::vector<std::pair<std::string, int>>();
    for( const auto& fieldName : ensemble.getCommonFieldNames()) {
        size_t numChannels = ensemble.getNumChannels(fieldName);
        for(size_t channel=0; channel<numChannels; channel++) {
            fields.emplace_back(std::pair<std::string, int>(fieldName, channel));
        }
    }

    if(fields != fields_) {

        propertyFields_.blockCallbacks(true);
        propertyFields_.reset();

        axesLabels_.clear();
        for( const auto& fieldName : ensemble.getCommonFieldNames()) {
            size_t numChannels = ensemble.getNumChannels(fieldName);
            if(numChannels > 1) {
                for(size_t channel=0; channel<numChannels; channel++) {
                    // We add the channel index to the field name (using 1-based counting) to identify the channels.
                    std::string axisLabel = fieldName + " (" + std::to_string(channel + 1) + ')';
                    axesLabels_.emplace_back(axisLabel);
                    propertyFields_.addRow(axisLabel);
                }
            }
            else {
                axesLabels_.emplace_back(fieldName);
                propertyFields_.addRow(fieldName);
            }
        }

        propertyFields_.blockCallbacks(false);
        propertyFields_.invalidate();

        propertyTemporalSampleCount_.setMaxValue(ensemble.getMaxNumTimeSteps());
        //propertyTemporalSampleCount_.set(ensemble.getMaxNumTimeSteps());

        fields_ = fields;
    }
}

ParallelCoordinatesAxesCreatorInput ParallelCoordinatesAxesCreator::prepareComputeInput() {

    const EnsembleDataset* ensemble = ensembleport_.getThreadSafeData();
    if(!ensemble) {
        throw InvalidInputException("No input", InvalidInputException::S_ERROR);
    }

    tgt::Bounds bounds;
    if (propertySampleRegion_.get() == "bounds") {
        bounds = ensemble->getBounds();
    }
    else if (propertySampleRegion_.get() == "common") {
        bounds = ensemble->getCommonBounds();
    }
    else {
        throw InvalidInputException("Unknown sample region", InvalidInputException::S_ERROR);
    }

    if(!bounds.isDefined()) {
        throw InvalidInputException("Bounding box is empty", InvalidInputException::S_ERROR);
    }

    // Set up random generator.
    std::function<float()> rnd(
            std::bind(std::uniform_real_distribution<float>(0.0f, 1.0f), std::mt19937(propertySeedTime_.get())));


    const VolumeBase* seedMask = seedMask_.getThreadSafeData();
    auto numSeedPoints = static_cast<size_t>(propertySpatialSampleCount_.get());
    std::vector<tgt::vec3> seedPoints;
    seedPoints.reserve(numSeedPoints);
    if (seedMask) {
        bounds.intersectVolume(seedMask->getBoundingBox().getBoundingBox());
        if(!bounds.isDefined()) {
            throw InvalidInputException("Seed Mask does not overlap with ensemble ROI", InvalidInputException::S_ERROR);
        }

        VolumeRAMRepresentationLock seedMaskLock(seedMask);
        tgt::mat4 seedMaskVoxelToWorldMatrix = seedMask->getVoxelToWorldMatrix();
        tgt::svec3 dim = seedMaskLock->getDimensions();
        for(size_t z=0; z < dim.z; z++) {
            for(size_t y=0; y < dim.y; y++) {
                for(size_t x=0; x < dim.x; x++) {
                    if(seedMaskLock->getVoxelNormalized(x, y, z) != 0.0f) {
                        tgt::vec3 pos = seedMaskVoxelToWorldMatrix * tgt::vec3(x, y, z);
                        if(bounds.containsPoint(pos)) {
                            seedPoints.emplace_back(pos);
                        }
                    }
                }
            }
        }

        if (seedPoints.empty()) {
            throw InvalidInputException("No seed points found in ROI", InvalidInputException::S_ERROR);
        }

        tgt::shuffle(seedPoints.begin(), seedPoints.end(), std::mt19937(propertySeedTime_.get()));
        seedPoints.resize(std::min(seedPoints.size(), numSeedPoints));

        LINFO("Restricting seed points to volume mask using " << seedPoints.size() << " seeds");
    }
    else {
        // Without a seed mask, we uniformly sample the whole space enclosed by the roi.
        for (size_t k = 0; k<numSeedPoints; k++) {
            // Since argument evaluation order is unspecified in c++, we need to ensure the order manually.
            float x = rnd();
            float y = rnd();
            float z = rnd();

            tgt::vec3 seedPoint(x, y, z);
            seedPoint = bounds.getLLF() + seedPoint * bounds.diagonal();
            seedPoints.push_back(seedPoint);
        }
    }

    tgtAssert(!seedPoints.empty(), "no seed points found");
    if (seedPoints.empty()) {
        throw InvalidInputException("No seed points found", InvalidInputException::S_ERROR);
    }

    // --- Gather values --- //
    std::unique_ptr<EnsembleDataset> selectedEnsemble(new EnsembleDataset());
    for( const auto member : propertyMembers_.get() )
        selectedEnsemble->addMember(ensembleport_.getData()->getMembers()[member]);

    const auto temporalSampleCount = propertyTemporalSampleCount_.get();
    const auto& members = ensemble->getMembers();

    if(members.empty()) {
        throw InvalidInputException("No members selected", InvalidInputException::S_ERROR);
    }

    auto fields = std::vector<std::pair<std::string, int>>();
    for( const auto field : propertyFields_.get() )
        fields.push_back( fields_[field] );

    auto axesLabels = std::vector<std::string>();
    for( const auto field : propertyFields_.get() )
        axesLabels.push_back(axesLabels_[field] );

    if(fields.empty()) {
        throw InvalidInputException("No fields selected", InvalidInputException::S_ERROR);
    }

    // --- Gather ranges --- //
    auto ranges = std::vector<tgt::vec2>( fields.size() );
    for( size_t i = 0; i < fields.size(); ++i ) {
        ranges[i] = ensemble->getValueRange( fields[i].first );
    }

    return ComputeInput {
        std::move(selectedEnsemble),
        EnsembleHash(*ensembleport_.getData()).getHash(), // Original hash.
        temporalSampleCount,
        propertyAggregateMembers_.get(),
        std::move(seedPoints),
        std::move(fields),
        std::move(axesLabels),
        std::move(ranges)
    };
}

ParallelCoordinatesAxesCreatorOutput ParallelCoordinatesAxesCreator::compute(ComputeInput input, ProgressReporter& progressReporter) const {

    auto ensemble = std::move(input.ensemble);
    const auto& members = ensemble->getMembers();
    size_t temporalSampleCount = input.temporalSampleCount;
    auto timeInterval = ensemble->getCommonTimeInterval();
    auto seedPoints = std::move(input.seedPoints);
    auto fields = std::move(input.fields);
    auto axesLabels = std::move(input.axesLabels);
    auto ranges = std::move(input.ranges);
    std::string ensembleHash = std::move(input.ensembleHash);

    if( input.aggregate ) {
        auto memberNames = std::vector<std::string> { "Aggregated" };
        auto values = std::vector<float>(members.size() * temporalSampleCount * fields.size() * seedPoints.size() );
        auto* currentValue = values.data();
        for(size_t i = 0; i < temporalSampleCount; ++i, currentValue += members.size() * fields.size() * seedPoints.size() )
        {
            for( size_t j = 0; j < fields.size(); ++j )
            {
                auto* dst = currentValue + j;
                for(const auto& member : members)
                {
                    auto t = mapRange(i, tgt::svec2(0, temporalSampleCount), timeInterval);
                    const auto& timestep = member.getTimeSteps()[member.getTimeStep(t)];
                    const auto& fieldName = fields[j].first;
                    const auto& channel = fields[j].second;

                    LDEBUG("[ParallelCoordinatesAxesCreator] Collecting voxels: Member=" << member.getName() << ", Timestep=" << j << ", Field=" << fieldName);

                    const auto* volumeHandle = timestep.getVolume(fieldName);
                    tgt::Bounds bounds = volumeHandle->getBoundingBox().getBoundingBox();
                    tgt::mat4 worldToVoxel = volumeHandle->getWorldToVoxelMatrix();
                    RealWorldMapping rwm = volumeHandle->getRealWorldMapping();
                    const auto volume = VolumeRAMRepresentationLock( volumeHandle );

                    for( const tgt::vec3& seedPoint : seedPoints )
                    {
                        if(bounds.containsPoint(seedPoint)) {
                            *dst = rwm.normalizedToRealWorld(volume->getVoxelNormalized( worldToVoxel * seedPoint, channel ));
                        }
                        else {
                            *dst = 0.0f; // TODO: 0 might not be a good choice depending on the dataset.
                        }
                        dst += fields.size();
                    }
                }
            }
            progressReporter.setProgress(1.0f * i / temporalSampleCount);
        }

        return ComputeOutput {
            std::unique_ptr<ParallelCoordinatesAxes>(
                    new ParallelCoordinatesAxes(
                            ensembleHash,
                            std::move(memberNames),
                            std::move(fields),
                            std::move(axesLabels),
                            std::move(ranges),
                            std::move(values),
                            temporalSampleCount,
                            seedPoints.size() * members.size()
                        )
                    )
        };
    }
    else {
        auto memberNames = std::vector<std::string>( members.size() );
        auto values = std::vector<float>(members.size() * temporalSampleCount * fields.size() * seedPoints.size() );

        for( size_t i = 0; i < members.size(); ++i )
        {
            const auto& member = members[i];
            memberNames[i] = member.getName();

            auto* currentValue = values.data() + i * (temporalSampleCount * fields.size() * seedPoints.size() );
            SubtaskProgressReporter memberProgress(progressReporter, tgt::vec2(i, i+1) / tgt::vec2(members.size()));
            for(size_t j = 0; j < temporalSampleCount; ++j, currentValue += fields.size() * seedPoints.size() )
            {
                auto t = mapRange(j, tgt::svec2(0, temporalSampleCount), timeInterval);
                const auto& timestep = member.getTimeSteps()[member.getTimeStep(t)];
                for( size_t k = 0; k < fields.size(); ++k )
                {
                    const auto& fieldName = fields[k].first;
                    const auto& channel = fields[k].second;
                    LDEBUG("[ParallelCoordinatesAxesCreator] Collecting voxels: Member=" << member.getName() << ", Timestep=" << j << ", Field=" << fieldName);

                    const auto* volumeHandle = timestep.getVolume(fieldName);
                    tgt::Bounds bounds = volumeHandle->getBoundingBox().getBoundingBox();
                    tgt::mat4 worldToVoxel = volumeHandle->getWorldToVoxelMatrix();
                    RealWorldMapping rwm = volumeHandle->getRealWorldMapping();
                    const auto volume = VolumeRAMRepresentationLock( volumeHandle );

                    auto* dst = currentValue + k;
                    for( const tgt::vec3& seedPoint : seedPoints )
                    {
                        if(bounds.containsPoint(seedPoint)) {
                            *dst = rwm.normalizedToRealWorld(volume->getVoxelNormalized( worldToVoxel * seedPoint, channel ));
                        }
                        else {
                            *dst = 0.0f; // TODO: 0 might not be a good choice depending on the dataset.
                        }
                        dst += fields.size();
                    }
                    memberProgress.setProgress(1.0f * j / temporalSampleCount);
                }
            }
        }

        return ComputeOutput {
            std::unique_ptr<ParallelCoordinatesAxes>(
                    new ParallelCoordinatesAxes(
                            ensembleHash,
                            std::move(memberNames),
                            std::move(fields),
                            std::move(axesLabels),
                            std::move(ranges),
                            std::move(values),
                            temporalSampleCount,
                            seedPoints.size()
                        )
                    )
        };
    }
}

void ParallelCoordinatesAxesCreator::processComputeOutput(ComputeOutput output) {
    axesport_.setData(output.axes.release());
}

std::vector<std::reference_wrapper<Port>> ParallelCoordinatesAxesCreator::getCriticalPorts() {
    auto criticalPorts = AsyncComputeProcessor<ComputeInput, ComputeOutput>::getCriticalPorts();
    criticalPorts.erase(std::remove_if(criticalPorts.begin(), criticalPorts.end(), [this] (const std::reference_wrapper<Port>& port){
        return port.get().getID() == seedMask_.getID();
    }), criticalPorts.end());
    return criticalPorts;
}

void ParallelCoordinatesAxesCreator::serialize(Serializer& s) const {
    AsyncComputeProcessor::serialize(s);
    s.serialize("hash", hash_);
    s.serialize("fields", fields_);
    s.serialize("axesLabels", axesLabels_);
}

void ParallelCoordinatesAxesCreator::deserialize(Deserializer& s) {
    AsyncComputeProcessor::deserialize(s);
    s.optionalDeserialize("hash", hash_, std::string(""));
    s.optionalDeserialize("fields", fields_, decltype(fields_)());
    s.optionalDeserialize("axesLabels", axesLabels_, decltype(axesLabels_)());

    if(fields_.size() != axesLabels_.size() || fields_.size() != propertyFields_.getNumRows()) {
        LWARNINGC("voreen.EnsembleAnalysis.ParallelCoordinatesAxesCreator",
                "Unexpected attribute mismatch after deserialization. Resetting selection.");

        hash_.clear();
        fields_.clear();
        axesLabels_.clear();

        propertyMembers_.reset();
        propertyFields_.reset();
    }
}

}