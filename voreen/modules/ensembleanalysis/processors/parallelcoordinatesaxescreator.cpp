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

#include <random>
#include <tuple>

namespace voreen {

ParallelCoordinatesAxesCreator::ParallelCoordinatesAxesCreator()
    : AsyncComputeProcessor<ParallelCoordianesAxesCreatorInput, ParallelCoordianesAxesCreatorOutput>()
    , ensembleport_(Port::INPORT, "port_ensemble", "Ensemble Input", false, Processor::VALID )
    , volumeport_(Port::INPORT, "port_volume", "Volume Mask Input", false, Processor::VALID )
    , axesport_(Port::OUTPORT, "port_axes", "Parallel Coordinates Axes" )
    , propertyRuns_("property_runs", "Selected Runs", Processor::VALID )
    , propertyFields_("property_fields", "Selected Fields", Processor::VALID )
    , propertySampleCount_("property_sample_count", "Sample Count", 1, 1, std::numeric_limits<int>::max(), Processor::VALID )
    , propertyAggregateRuns_("property_aggregate_runs", "Aggregate Runs", false, Processor::VALID )
    , propertyFileDialog_("property_file_dialog", "File Output", "Select File...", "", "Parallel Coordinates (*.pc)", FileDialogProperty::FileMode::SAVE_FILE, Processor::VALID )
    , propertySaveButton_("property_save", "Save File", Processor::VALID )
{
    // --- Initialize Ports --- //
    this->addPort(ensembleport_ );
    this->addPort(volumeport_ );
    this->addPort(axesport_ );

    // --- Initialize Properties --- //
    this->addProperty(propertyRuns_ );
    this->addProperty(propertyFields_ );
    this->addProperty(propertySampleCount_ );
    this->addProperty(propertyAggregateRuns_ );
    this->addProperty(propertyFileDialog_ );
    this->addProperty(propertySaveButton_ );

    // --- Initialize Callbacks --- //
    ensembleport_.onNewData(LambdaFunctionCallback([this] {
        propertyRuns_.blockCallbacks(true );
        propertyRuns_.reset();
        if( ensembleport_.hasData() )
            for( const auto& run : ensembleport_.getData()->getRuns() )
                propertyRuns_.addRow(run.getName(), run.getColor() );
        propertyRuns_.blockCallbacks(false );
        propertyRuns_.invalidate();
    } ) );
    volumeport_.onChange(MemberFunctionCallback<ParallelCoordinatesAxesCreator>(this, &ParallelCoordinatesAxesCreator::updateValidVoxels ) );
    propertyRuns_.onChange(LambdaFunctionCallback([this]{
        auto ensemble = EnsembleDataset();
        for( const auto run : propertyRuns_.get() )
            ensemble.addRun(ensembleport_.getData()->getRuns()[run] );

        propertyFields_.blockCallbacks(true );
        propertyFields_.reset();
        for( const auto& field : ensemble.getCommonFieldNames() )
            propertyFields_.addRow(field);
        propertyFields_.blockCallbacks(false );
        propertyFields_.invalidate();
    } ) );
    propertyFields_.onChange(MemberFunctionCallback<ParallelCoordinatesAxesCreator>(this, &ParallelCoordinatesAxesCreator::updateValidVoxels ) );
    propertySaveButton_.onChange(LambdaFunctionCallback([this] {
        if(axesport_.hasData() && !propertyFileDialog_.get().empty() )
            axesport_.getData()->serialize(propertyFileDialog_.get() );
    } ) );
}

Processor* ParallelCoordinatesAxesCreator::create() const {
    return new ParallelCoordinatesAxesCreator();
}
std::string ParallelCoordinatesAxesCreator::getClassName() const {
    return "ParallelCoordinatesAxesCreator";
}
std::string ParallelCoordinatesAxesCreator::getCategory() const {
    return "ParallelCoordinates";
}
bool ParallelCoordinatesAxesCreator::isReady() const {
    return ensembleport_.isReady() && axesport_.isReady();
}

void ParallelCoordinatesAxesCreator::updateValidVoxels() {
    if( volumeport_.hasData() ) {
        const auto volume = volumeport_.getData()->getRepresentation<VolumeRAM>();

        validVoxels_.clear();
        validVoxels_.reserve(volume->getNumVoxels() );
        for( size_t i = 0; i < volume->getNumVoxels(); ++i )
            if( volume->getVoxelNormalized( i ) )
                validVoxels_.push_back(static_cast<int32_t>( i ) );
        validVoxels_.shrink_to_fit();
    }
    else {
        if( ensembleport_.hasData() )
        {
            auto ensemble = EnsembleDataset();
            for( const auto run : propertyRuns_.get() )
                ensemble.addRun(ensembleport_.getData()->getRuns()[run] );

            if(ensemble.getRuns().size() && propertyFields_.get().size() ) {
                const auto field = ensemble.getCommonFieldNames()[propertyFields_.get().front()];
                const auto voxelCount = ensemble.getRuns().front().getTimeSteps().front().getVolume(field)->getNumVoxels();

                validVoxels_.resize(voxelCount );
                std::iota(validVoxels_.begin(), validVoxels_.end(), 0 );
            }
        }
        else validVoxels_.clear();
    }
    propertySampleCount_.setMaxValue(static_cast<int>( validVoxels_.size() ) );
}

ParallelCoordianesAxesCreatorInput ParallelCoordinatesAxesCreator::prepareComputeInput() {

    if(!ensembleport_.hasData()) {
        throw InvalidInputException("No input", InvalidInputException::S_ERROR);
    }

    // --- Shuffle valid voxels --- //
    const auto sampleCount = propertySampleCount_.get();
    auto validVoxels = validVoxels_;
    std::shuffle( validVoxels.begin(), validVoxels.end(), std::mt19937( std::random_device()( ) ) );
    validVoxels.resize( sampleCount );

    // --- Gather values --- //
    auto ensemble = EnsembleDataset();
    for( const auto run : propertyRuns_.get() )
        ensemble.addRun(ensembleport_.getData()->getRuns()[run] );

    auto fieldNames = std::vector<std::string>();
    for( const auto field : propertyFields_.get() )
        fieldNames.push_back( ensemble.getCommonFieldNames()[field] );

    const auto numTimesteps = ensemble.getMinNumTimeSteps();
    const auto& runs = ensemble.getRuns();

    if( runs.empty() || fieldNames.empty() || numTimesteps == 0 ) {
        throw InvalidInputException("Empty input", InvalidInputException::S_ERROR);
    }

    // --- Gather ranges --- //
    auto ranges = std::vector<std::pair<float, float>>( fieldNames.size() );
    for( size_t i = 0; i < fieldNames.size(); ++i ) {
        const auto range = ensemble.getValueRange( fieldNames[i] );
        ranges[i] = std::make_pair( range.x, range.y );
    }

    return ComputeInput {
        std::move(ensemble),
        sampleCount,
        numTimesteps,
        std::move(validVoxels),
        std::move(fieldNames),
        std::move(ranges)
    };
}

ParallelCoordianesAxesCreatorOutput ParallelCoordinatesAxesCreator::compute(ComputeInput input, ProgressReporter& progressReporter) const {

    const auto& runs = input.ensemble.getRuns();
    size_t sampleCount = input.sampleCount;
    size_t numTimesteps = input.numTimeSteps;
    auto validVoxels = std::move(input.validVoxels);
    auto fieldNames = std::move(input.fieldNames);
    auto ranges = std::move(input.ranges);

    if( propertyAggregateRuns_.get() ) {
        auto runNames = std::vector<std::string> { "Aggregated" };
        auto values = std::vector<float>( runs.size() * numTimesteps * fieldNames.size() * sampleCount );
        auto currentValue = values.data();
        for( size_t i = 0; i < numTimesteps; ++i, currentValue += runs.size() * fieldNames.size() * sampleCount )
        {
            for( size_t j = 0; j < fieldNames.size(); ++j )
            {
                auto dst = currentValue + j;
                for( size_t k = 0; k < runs.size(); ++k )
                {
                    const auto& run = runs[k];
                    const auto& timestep = run.getTimeSteps()[i];
                    const auto& fieldName = fieldNames[j];

                    std::cout << "[ParallelCoordinatesAxesCreator] Collecting voxels: Run=" << run.getName() << ", Timestep=" << j << ", Field=" << fieldName << std::endl;

                    const auto volumeHandle = timestep.getVolume(fieldName);
                    const auto volume = VolumeRAMRepresentationLock( volumeHandle );

                    for( const auto index : validVoxels )
                    {
                        *dst = static_cast<float>( volume->getVoxelNormalized( index ) );
                        dst += fieldNames.size();
                    }
                }
            }
            progressReporter.setProgress(1.0f * i / numTimesteps);
        }

        return ComputeOutput {
            std::unique_ptr<ParallelCoordinatesAxes>(new ParallelCoordinatesAxes(std::move(runNames), std::move(fieldNames), std::move(ranges), std::move(values), numTimesteps, sampleCount * runs.size() ) )
        };
    }
    else {
        auto runNames = std::vector<std::string>( runs.size() );
        auto values = std::vector<float>( runs.size() * numTimesteps * fieldNames.size() * sampleCount );

        for( size_t i = 0; i < runs.size(); ++i )
        {
            const auto& run = runs[i];
            runNames[i] = run.getName();

            auto currentValue = values.data() + i * ( numTimesteps * fieldNames.size() * sampleCount );
            SubtaskProgressReporter runProgress(progressReporter, tgt::vec2(i, i+1) / tgt::vec2(runs.size()));
            for( size_t j = 0; j < numTimesteps; ++j, currentValue += fieldNames.size() * sampleCount )
            {
                const auto& timestep = run.getTimeSteps()[j];
                for( size_t k = 0; k < fieldNames.size(); ++k )
                {
                    const auto& fieldName = fieldNames[k];
                    std::cout << "[ParallelCoordinatesAxesCreator] Collecting voxels: Run=" << run.getName() << ", Timestep=" << j << ", Field=" << fieldName << std::endl;

                    const auto volumeHandle = timestep.getVolume(fieldName);
                    const auto volume = VolumeRAMRepresentationLock( volumeHandle );
                    const auto dim = volume->getDimensions();

                    auto dst = currentValue + k;
                    for( const auto index : validVoxels )
                    {
                        if( dim.z == 1 ) *dst = static_cast<float>( volume->getVoxelNormalized( index % ( dim.x * dim.y ) ) );
                        else *dst = static_cast<float>( volume->getVoxelNormalized( index ) );
                        dst += fieldNames.size();
                    }
                    runProgress.setProgress(1.0f * j / numTimesteps);
                }
            }
        }

        return ComputeOutput {
            std::unique_ptr<ParallelCoordinatesAxes>(new ParallelCoordinatesAxes(std::move(runNames ), std::move(fieldNames ), std::move(ranges ), std::move(values ), numTimesteps, sampleCount ) )
        };
    }
}

void ParallelCoordinatesAxesCreator::processComputeOutput(ComputeOutput output) {
    axesport_.setData(output.axes.release());
}

std::vector<std::reference_wrapper<Port>> ParallelCoordinatesAxesCreator::getCriticalPorts() {
    auto criticalPorts = AsyncComputeProcessor<ComputeInput, ComputeOutput>::getCriticalPorts();
    criticalPorts.erase(std::remove_if(criticalPorts.begin(), criticalPorts.end(), [this] (const std::reference_wrapper<Port>& port){
        return port.get().getID() == volumeport_.getID();
    }), criticalPorts.end());
    return criticalPorts;
}

}