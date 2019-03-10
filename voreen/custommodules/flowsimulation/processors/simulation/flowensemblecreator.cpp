/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "flowensemblecreator.h"

#include "voreen/core/ports/conditions/portconditionvolumetype.h"

#include "modules/hdf5/io/hdf5volumereader.h"
#include "modules/hdf5/io/hdf5volumewriter.h"
#include "modules/hdf5/io/hdf5filevolume.h"

// Include actual Volume Filters.
#include "../../flowfeatures/magnitudefeature.h"
//#include "../../flowfeatures/curlfeature.h"       // TODO: implement
//#include "../../flowfeatures/divergencefeature.h" // TODO: implement
//#include "../../flowfeatures/helicityfeature.h"   // TODO: implement

namespace voreen {

const std::string FlowEnsembleCreator::loggerCat_("voreen.flowsimulation.FlowEnsembleCreator");

FlowEnsembleCreator::FlowEnsembleCreator()
    : AsyncComputeProcessor<FlowEnsembleCreatorInput, FlowEnsembleCreatorOutput>()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , deleteOriginalData_("deleteOriginalData", "Delete original Data", false)
    , simulationResultPath_("simulationResultPath", "Simulation Result Path", "Path", "", "", FileDialogProperty::DIRECTORY, Processor::INVALID_RESULT, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
    , ensembleOutputPath_("ensembleOutputPath", "Ensemble Output Path", "Path", "", "", FileDialogProperty::DIRECTORY, Processor::INVALID_RESULT, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
    , outputVolumeDeflateLevel_("outputVolumeDeflateLevel", "Deflate Level", 1, 0, 9, Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_DEFAULT)
    , featureList_("featureList", "Feature List", false)
{
    addPort(inport_);
    inport_.addCondition(new PortConditionVolumeList(new PortConditionVolumeType3xFloat()));

    addProperty(featureList_);
    featureList_.setGroupID("feature");
    featureList_.setDuplicationAllowed(false);
    setPropertyGroupGuiName("feature", "Feature");

    // Add features.
    //addFeature(new MagnitudeFeature()); // Will always be added implicitly.
    //addFeature(new CurlFeature())

    // Technical stuff.
    addProperty(deleteOriginalData_);
    deleteOriginalData_.setGroupID("output");
    addProperty(ensembleOutputPath_);
    ensembleOutputPath_.setGroupID("output");
    addProperty(outputVolumeDeflateLevel_);
    outputVolumeDeflateLevel_.setGroupID("output");
    setPropertyGroupGuiName("output", "Output");
}

FlowEnsembleCreator::~FlowEnsembleCreator() {
}

bool FlowEnsembleCreator::isReady() const {
    if(!isInitialized()) {
        setNotReadyErrorMessage("Not initialized.");
        return false;
    }
    if(!inport_.isReady()) {
        setNotReadyErrorMessage("Inport not ready.");
        return false;
    }
    return true;
}

Processor* FlowEnsembleCreator::create() const {
    return new FlowEnsembleCreator();
}

FlowEnsembleCreatorInput FlowEnsembleCreator::prepareComputeInput() {

    if(!inport_.hasData()) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    if(featureList_.getInstances().empty()) {
        throw InvalidInputException("No filter selected", InvalidInputException::S_ERROR);
    }

    const VolumeList* volumeList = inport_.getThreadSafeData();

    return FlowEnsembleCreatorInput{
            simulationResultPath_.get(),
            ensembleOutputPath_.get(),
            volumeList
    };
}
FlowEnsembleCreatorOutput FlowEnsembleCreator::compute(FlowEnsembleCreatorInput input, ProgressReporter& progressReporter) const {

/*
    std::vector<std::string> runs = tgt::FileSystem::listSubDirectories(input.simulationResultPath, true);
    float progressPerRun = 1.0f / runs.size();
    for(const std::string& run : runs) {
        std::string runPath = input.simulationResultPath + "/" + run;
        std::vector<std::string> fileNames = tgt::FileSystem::readDirectory(runPath, true, false);

    }

    for(size_t i=0; i<input.measuredData->size(); i++) {

        const VolumeBase& volume = *input.measuredData->at(i);

        VolumeFilterStackBuilder builder(volume);
        std::string baseType = volume.getBaseType();
        for (const InteractiveListProperty::Instance &instance : featureList_.getInstances()) {

            // Create new instance of selected feature aka. filter.
            VolumeFilter *filter = flowFeatures_[instance.itemId_]->create();

            // Base type of output volume is determined by last filter output type.
            baseType = filter->getSliceBaseType();

            builder.addLayer(std::unique_ptr<VolumeFilter>(filter));
        }

        std::unique_ptr<SliceReader> sliceReader = builder.build(0);

        const std::string volumeFilePath = outputVolumeFilePath_.get();
        const std::string volumeLocation = HDF5VolumeWriter::VOLUME_DATASET_NAME;
        const tgt::svec3 dim = volume.getDimensions();

        if (volumeFilePath.empty()) {
            throw InvalidInputException("No volume file path specified!", InvalidInputException::S_ERROR);
        }

        std::unique_ptr<HDF5FileVolume> outputVolume = nullptr;
        try {
            outputVolume = std::unique_ptr<HDF5FileVolume>(
                    HDF5FileVolume::createVolume(volumeFilePath, volumeLocation, baseType, dim, 1, true,
                                                 outputVolumeDeflateLevel_.get(), tgt::svec3(dim.x, dim.y, 1), false));
        } catch (tgt::IOException e) {
            throw InvalidInputException("Could not create output volume.", InvalidInputException::S_ERROR);
        }

        outputVolume->writeSpacing(volume.getSpacing());
        outputVolume->writeOffset(volume.getOffset());
        outputVolume->writeRealWorldMapping(volume.getRealWorldMapping());

        writeSlicesToHDF5File(*sliceReader, *outputVolume, &progressReporter);
    }
    */

    return {};
}
void FlowEnsembleCreator::processComputeOutput(FlowEnsembleCreatorOutput output) {
}

// private methods
//

void FlowEnsembleCreator::addFeature(FlowFeature* feature) {
    featureList_.addItem(feature->getName());
    flowFeatures_.push_back(std::unique_ptr<FlowFeature>(feature));
}

}   // namespace
