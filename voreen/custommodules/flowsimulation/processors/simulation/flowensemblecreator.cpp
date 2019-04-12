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

#include "voreen/core/datastructures/volume/volumedecorator.h"

#include "modules/core/io/vvdvolumewriter.h"

namespace voreen {

const std::string FlowEnsembleCreator::loggerCat_("voreen.flowsimulation.FlowEnsembleCreator");

FlowEnsembleCreator::FlowEnsembleCreator()
    : AsyncComputeProcessor<FlowEnsembleCreatorInput, FlowEnsembleCreatorOutput>()
    , magnitudeInport_(Port::INPORT, "magnitudeInput", "Measured Magnitude Input (Optional)", false)
    , velocityInport_(Port::INPORT, "velocityInput", "Measured Velocity Input (Optional)", false)
    , deleteOriginalData_("deleteOriginalData", "Delete original Data", false)
    , simulationResultPath_("simulationResultPath", "Simulation Result Path", "Path", "", "", FileDialogProperty::DIRECTORY, Processor::INVALID_RESULT, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
    , ensembleOutputPath_("ensembleOutputPath", "Ensemble Output Path", "Path", "", "", FileDialogProperty::DIRECTORY, Processor::INVALID_RESULT, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
    , measuredDataPath_("measuredDataPath", "Measured Data Path (Optional)", "Path", "", "", FileDialogProperty::DIRECTORY, Processor::INVALID_RESULT, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
    , measuredDataName_("measuredDataName", "Measured Data Name", "4d_pc_mri")
    , simulationTime_("simulationTime", "Simulation Time (Link with FlowCharacteristics)", 1.0f, 0.0f, 20.0f)
    , temporalResolution_("temporalResolution", "Temporal Resolution (Link with FlowCharacteristics", 3.1f, 0.1f, 100.0f)
{
    addPort(magnitudeInport_);
    magnitudeInport_.addCondition(new PortConditionVolumeListEnsemble());
    magnitudeInport_.addCondition(new PortConditionVolumeListAdapter(new PortConditionVolumeTypeFloat()));
    addPort(velocityInport_);
    velocityInport_.addCondition(new PortConditionVolumeListEnsemble());
    velocityInport_.addCondition(new PortConditionVolumeListAdapter(new PortConditionVolumeType3xFloat()));

    // Technical stuff.
    addProperty(deleteOriginalData_);
    deleteOriginalData_.setGroupID("output");
    addProperty(simulationResultPath_);
    simulationResultPath_.setGroupID("output");
    addProperty(ensembleOutputPath_);
    ensembleOutputPath_.setGroupID("output");
    addProperty(measuredDataPath_);
    measuredDataPath_.setGroupID("output");
    //addProperty(measuredDataName_); // TODO: implement functionality
    measuredDataName_.setGroupID("output");
    setPropertyGroupGuiName("output", "Output");

    addProperty(simulationTime_);
    simulationTime_.setVisibleFlag(false);
    addProperty(temporalResolution_);
    temporalResolution_.setVisibleFlag(false);
}

FlowEnsembleCreator::~FlowEnsembleCreator() {
}

bool FlowEnsembleCreator::isReady() const {
    if(!isInitialized()) {
        setNotReadyErrorMessage("Not initialized.");
        return false;
    }
    return true;
}

Processor* FlowEnsembleCreator::create() const {
    return new FlowEnsembleCreator();
}

FlowEnsembleCreatorInput FlowEnsembleCreator::prepareComputeInput() {
    return FlowEnsembleCreatorInput{
            simulationResultPath_.get(),
            measuredDataPath_.get(),
            ensembleOutputPath_.get(),
            measuredDataName_.get(),
            magnitudeInport_.getThreadSafeData(),
            velocityInport_.getThreadSafeData(),
            deleteOriginalData_.get()
    };
}

FlowEnsembleCreatorOutput FlowEnsembleCreator::compute(FlowEnsembleCreatorInput input, ProgressReporter& progressReporter) const {

    // Map Old -> New filename.
    std::map<std::string, std::string> files;

    // Gather all files.
    std::vector<std::string> ensembles = tgt::FileSystem::listSubDirectories(input.simulationResultPath);
    for(size_t i=0; i<ensembles.size(); i++) {
        const std::string& ensemble = ensembles[i];
        std::string ensemblePath = input.simulationResultPath + "/" + ensemble;
        std::vector<std::string> runs = tgt::FileSystem::listSubDirectories(ensemblePath);
        for (size_t j = 0; j < runs.size(); j++) {
            const std::string& run = runs[j];
            std::string runPath = ensemblePath + "/" + run;
            std::vector<std::string> fileNames = tgt::FileSystem::readDirectory(runPath, true, false);

            // Only look for .raw files.
            for (const std::string& file : fileNames) {
                std::string oldPath = runPath + "/" + file;

                size_t underscore = file.find_last_of('_');
                size_t dot = file.find_last_of('.');

                if (underscore != std::string::npos && dot != std::string::npos &&
                        (file.substr(dot) == ".vvd" || file.substr(dot) == ".raw")) {

                    std::string channel = file.substr(0, underscore);
                    std::string iteration = file.substr(underscore + 1);

                    std::string newPath = input.ensembleOutputPath + "/";
                    newPath += channel   + "/";
                    newPath += ensemble  + "/";
                    //newPath += run       + "/t_"; // VVD file contains filename that would have to be changed.
                    newPath += run  + "/" + channel + "_";
                    newPath += iteration;
                    files[oldPath] = newPath;
                }
            }
        }

        // Handle measured data.
        writeMeasuredData(input.magnitudeData, input.ensembleOutputPath, ensemble, "magnitude", input.measuredDataName);
        writeMeasuredData(input.velocityData,  input.ensembleOutputPath, ensemble, "velocity",  input.measuredDataName);
    }

    // Copy / Move.
    size_t fileIdx = 0;
    float progressPerFile = 1.0f / files.size();
    for(const auto& file : files) {

        // Create output directory.
        tgt::FileSystem::createDirectoryRecursive(tgt::FileSystem::dirName(file.second));

        if(input.deleteOriginalData) {
            tgt::FileSystem::renameFile(file.first, file.second, false);
        }
        else {
            tgt::FileSystem::copyFile(file.first, file.second);
        }

        fileIdx++;
        progressReporter.setProgress(fileIdx*progressPerFile);
    }

    return {};
}

void FlowEnsembleCreator::processComputeOutput(FlowEnsembleCreatorOutput output) {
    // Delete old data, if desired.
    if(deleteOriginalData_.get()) {
        tgt::FileSystem::deleteDirectoryRecursive(simulationResultPath_.get());
    }
}

void FlowEnsembleCreator::writeMeasuredData(const VolumeList* measuredData,
                                            const std::string& ensembleOutputPath,
                                            const std::string& ensemble,
                                            const std::string& channel,
                                            const std::string& measuredDataName) const {

    if (measuredData && !measuredData->empty()) {

        std::string runPath = ensembleOutputPath + "/" + channel + "/" + ensemble + "/" + measuredDataName;
        tgt::FileSystem::createDirectoryRecursive(runPath);

        if (measuredData->size() == 1) {
            const VolumeBase *original = measuredData->first();
            VvdVolumeWriter().write(runPath + "/t0.vvd", original);
            std::unique_ptr<VolumeBase> duplicate(new VolumeDecoratorReplaceTimestep(original, 1.0f));
            VvdVolumeWriter().write(runPath + "/t1.vvd", duplicate.get());
        } else {
            for (size_t j = 0; j < measuredData->size(); j++) {
                float timeStep = j * temporalResolution_.get();
                std::unique_ptr<VolumeBase> original(
                        new VolumeDecoratorReplaceTimestep(measuredData->at(j), timeStep));
                VvdVolumeWriter().write(runPath + "/t" + std::to_string(j) + ".vvd", original.get());
            }
        }
    }
}

}   // namespace
