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

#include "voreen/core/ports/conditions/portconditionvolumelist.h"

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
    , spatialResolution_("spatialResolution", "Spatial Resolution", 128, 32, 1024)
    , spacing_("spacing", "Spacing", tgt::vec3::one, tgt::vec3::zero, tgt::vec3(1000.0f))
    , offset_("offset", "Offset", tgt::vec3::zero, -tgt::vec3(1000.0f), tgt::vec3(1000.0f))
    , simulationResultPath_("simulationResultPath", "Simulation Result Path", "Path", "", "", FileDialogProperty::DIRECTORY, Processor::INVALID_RESULT, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
    , ensembleOutputPath_("ensembleOutputPath", "Ensemble Output Path", "Path", "", "", FileDialogProperty::DIRECTORY, Processor::INVALID_RESULT, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
    , outputVolumeDeflateLevel_("outputVolumeDeflateLevel", "Deflate Level", 1, 0, 9, Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_DEFAULT)
    , featureList_("featureList", "Feature List", false)
{
    addPort(inport_);
    inport_.addCondition(new PortConditionVolumeListEnsemble());
    inport_.addCondition(new PortConditionVolumeListAdapter(new PortConditionVolumeTypeFloat()));

    addProperty(featureList_);
    featureList_.setGroupID("feature");
    featureList_.setDuplicationAllowed(false);
    setPropertyGroupGuiName("feature", "Feature");

    // Add features.
    //addFeature(new MagnitudeFeature()); // Will always be added implicitly.
    //addFeature(new CurlFeature())

    // Input Properties.
    addProperty(spatialResolution_);
    spatialResolution_.setGroupID("input");
    addProperty(spacing_);
    spacing_.setGroupID("input");
    addProperty(offset_);
    offset_.setGroupID("input");
    setPropertyGroupGuiName("input", "Input");

    // Technical stuff.
    addProperty(deleteOriginalData_);
    deleteOriginalData_.setGroupID("output");
    addProperty(simulationResultPath_);
    simulationResultPath_.setGroupID("output");
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

    // TODO: make more dynamic!
    std::map<std::string, size_t> CHANNEL_COUNT;
    CHANNEL_COUNT["pressure"] = 1;
    CHANNEL_COUNT["velocity"] = 3;

    std::vector<std::string> runs = tgt::FileSystem::listSubDirectories(input.simulationResultPath, true);
    float progressPerRun = 1.0f / (runs.size() + input.measuredData->size());
    int maxIteration = 0; // We need to know the max iteration for the measured data.
    for(size_t i=0; i<runs.size(); i++) {
        const std::string &run = runs[i];
        std::string runPath = input.simulationResultPath + "/" + run;
        std::vector<std::string> fileNames = tgt::FileSystem::readDirectory(runPath, true, false);

        std::map<int, std::map<std::string, std::vector<std::string>>> timeSteps;

        // Only look for .raw files.
        for (const std::string &file : fileNames) {

            size_t underscore = file.find_first_of('_');
            size_t dot = file.find_first_of('.');

            if (underscore != std::string::npos && dot != std::string::npos &&
                file.substr(dot) == ".raw") {

                std::string property = file.substr(0, underscore);
                int iteration = std::atoi(file.substr(underscore + 1, dot).c_str());
                maxIteration = std::max(maxIteration, iteration);

                // HACK: count somehow different the amount of channels.
                for (size_t channel = 0; channel < CHANNEL_COUNT[property]; channel++) {
                    timeSteps[iteration][property].push_back(file);
                }
            }
        }

        SubtaskProgressReporter timeStepProgressReporter(progressReporter, tgt::vec2(i, i+1)*progressPerRun);
        size_t timeStepIdx = 0;
        float progressPerTimeStep = 1.0f / timeSteps.size();
        for(const auto& pair : timeSteps) {
            int iteration = pair.first;
            const std::map<std::string, std::vector<std::string>> properties = pair.second;

            const std::string volumeFilePath = input.ensembleOutputPath + "/" + run + "/" + std::to_string(iteration) + ".h5";
            const std::string volumeLocation = HDF5VolumeWriter::VOLUME_DATASET_NAME;
            const tgt::svec3 dim(spatialResolution_.get());

            std::unique_ptr<HDF5FileVolume> outputVolume = nullptr;
            try {
                outputVolume = std::unique_ptr<HDF5FileVolume>(
                        HDF5FileVolume::createVolume(volumeFilePath, volumeLocation, "float", dim, properties.size(), true,
                                                     outputVolumeDeflateLevel_.get(), tgt::svec3(dim.x, dim.y, 1), false));
            } catch (const tgt::IOException& e) {
                LERROR("Could not create output volume.");
                continue;
            }

            outputVolume->writeSpacing(spacing_.get());
            outputVolume->writeOffset(offset_.get());
            //outputVolume->writeRealWorldMapping(volume.getRealWorldMapping());

            size_t propertyIdx = 0;
            for(const auto& property : properties) {

                const std::string& name = property.first;
                const std::vector<std::string>& files = property.second;
                if(files.size() != 1) {
                    LWARNING("Multiple properties per property!");
                }
                std::string propertyFilename = files.front();
                std::fstream propertyFile(propertyFilename, std::ios::in | std::ios::binary);
                if(!propertyFile) {
                    LERROR("Can not read T=" << std::to_string(iteration) << ", Property=" << name);
                    continue;
                }

                // Initialize volume.
                std::unique_ptr<VolumeRAM> volumeData = nullptr;
                switch(CHANNEL_COUNT[name]) {
                case 1:
                    volumeData.reset(new VolumeRAM_Float(dim, true));
                    break;
                case 3:
                    volumeData.reset(new VolumeRAM_3xFloat(dim, true));
                    break;
                default:
                    break;
                }
                if(!volumeData) {
                    LERROR("Unkown number of channels for property " << name);
                    continue;
                }

                // Read volume data.
                propertyFile.read(reinterpret_cast<char*>(volumeData->getData()), volumeData->getNumBytes());
                if(!propertyFile) {
                    LERROR("Can not read data of T=" << std::to_string(iteration) << ", Property=" << name);
                    continue;
                }

                Volume* volume = new Volume(volumeData.release(), spacing_.get(), offset_.get());
                if(volume->getNumChannels() > 1) {

                    // If we have more than one channel (namely 3), apply flow features.
                    for (const InteractiveListProperty::Instance& instance : featureList_.getInstances()) {

                        VolumeFilterStackBuilder builder(*volume);
                        std::string baseType = volume->getBaseType();

                        // Create new instance of selected feature aka. filter.
                        VolumeFilter *filter = flowFeatures_[instance.itemId_]->create();

                        // Base type of output volume is determined by last filter output type.
                        baseType = filter->getSliceBaseType();

                        builder.addLayer(std::unique_ptr<VolumeFilter>(filter));

                        std::unique_ptr<SliceReader> sliceReader = builder.build(0);
                        sliceReader->seek(0);
                        for(size_t z = 0; z < dim.z; ++z) {
                            if(z > 0) {
                                sliceReader->advance();
                            }
                            outputVolume->writeSlices(sliceReader->getCurrentSlice(), z, 0, 1);
                        }

                    }
                }
                else {
                    outputVolume->writeVolume(volume->getRepresentation<VolumeRAM>(), propertyIdx, 1);
                }

                // TODO: use name somehow as meta data?
                propertyIdx++;
                //outputVolume->writeMetaData();
                //Volume* volume = new Volume(data, spacing_.get(), offset_.get());
                //volume->getMetaDataContainer().addMetaData("name", new StringMetaData(prop));
            }
            timeStepProgressReporter.setProgress(timeStepIdx * progressPerTimeStep);
            timeStepIdx++;
        }
        progressReporter.setProgress(i*progressPerRun);
    }
/*
    // Handle measured data.
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

        //const std::string volumeFilePath = outputVolumeFilePath_.get();
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
    // Delete old data, if desired.
    if(deleteOriginalData_.get()) {
        tgt::FileSystem::deleteDirectoryRecursive(simulationResultPath_.get());
    }
}

// private methods
//

void FlowEnsembleCreator::adjustPropertiesToInput() {

    const VolumeBase* reference = nullptr;
    if(inport_.hasData() && !inport_.getData()->empty()) {
        reference = inport_.getData()->first();
    }

    setPropertyGroupVisible("input", !reference);
    if(reference) {
        spatialResolution_.set(tgt::max(reference->getDimensions()));
        spacing_.set(reference->getSpacing());
        offset_.set(reference->getOffset());
    }
}

void FlowEnsembleCreator::addFeature(FlowFeature* feature) {
    featureList_.addItem(feature->getName());
    flowFeatures_.push_back(std::unique_ptr<FlowFeature>(feature));
}

}   // namespace
