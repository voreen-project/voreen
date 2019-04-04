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
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , deleteOriginalData_("deleteOriginalData", "Delete original Data", false)
    , simulationResultPath_("simulationResultPath", "Simulation Result Path", "Path", "", "", FileDialogProperty::DIRECTORY, Processor::INVALID_RESULT, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
    , ensembleOutputPath_("ensembleOutputPath", "Ensemble Output Path", "Path", "", "", FileDialogProperty::DIRECTORY, Processor::INVALID_RESULT, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
    , setting_("setting", "Setting", "4d_pc_mri")
    , simulationTime_("simulationTime", "Simulation Time (Link with FlowCharacteristics)", 1.0f, 0.0f, 20.0f)
    , temporalResolution_("temporalResolution", "Temporal Resolution (Link with FlowCharacteristics", 3.1f, 0.1f, 100.0f)
{
    addPort(inport_);
    inport_.addCondition(new PortConditionVolumeListEnsemble());
    PortConditionLogicalOr* typeCondition = new PortConditionLogicalOr();
    typeCondition->addLinkedCondition(new PortConditionVolumeTypeFloat());
    typeCondition->addLinkedCondition(new PortConditionVolumeType3xFloat());
    inport_.addCondition(new PortConditionVolumeListAdapter(typeCondition));

    // Technical stuff.
    addProperty(deleteOriginalData_);
    deleteOriginalData_.setGroupID("output");
    addProperty(simulationResultPath_);
    simulationResultPath_.setGroupID("output");
    addProperty(ensembleOutputPath_);
    ensembleOutputPath_.setGroupID("output");
    addProperty(setting_);
    setting_.setGroupID("output");
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

    const VolumeList* volumeList = inport_.getThreadSafeData();

    return FlowEnsembleCreatorInput{
            simulationResultPath_.get(),
            ensembleOutputPath_.get(),
            volumeList,
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
        if(input.measuredData && !input.measuredData->empty()) {

            size_t numChannels = input.measuredData->first()->getNumChannels();
            if (numChannels != 1 && numChannels != 3) {
                continue;
            }

            std::string setting = setting_.get();
            std::string channel = numChannels == 1 ? "magnitude" : "velocity";
            std::string runPath = input.ensembleOutputPath + "/" + channel + "/" + ensemble + "/" + setting;

            tgt::FileSystem::createDirectoryRecursive(runPath);

            if (input.measuredData->size() == 1) {
                const VolumeBase* original = input.measuredData->first();
                VvdVolumeWriter().write(runPath + "/t0.vvd", original);
                std::unique_ptr<VolumeBase> duplicate(new VolumeDecoratorReplaceTimestep(original, 1.0f));
                VvdVolumeWriter().write(runPath + "/t1.vvd", duplicate.get());
            } else {
                for (size_t j = 0; j < input.measuredData->size(); j++) {
                    float timeStep = j * temporalResolution_.get();
                    std::unique_ptr<VolumeBase> original(
                            new VolumeDecoratorReplaceTimestep(input.measuredData->at(j), timeStep));
                    VvdVolumeWriter().write(runPath + "/t" + std::to_string(j) + ".vvd", original.get());
                }
            }
        }
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

}   // namespace
