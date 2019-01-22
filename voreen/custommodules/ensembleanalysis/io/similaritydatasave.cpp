/***********************************************************************************
*                                                                                 *
* Voreen - The Volume Rendering Engine                                            *
*                                                                                 *
* Copyright (C) 2005-2017 University of Muenster, Germany.                        *
* Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "similaritydatasave.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/volume/volumedisk.h"

#include "../datastructures/similaritydata.h"

#include "tgt/filesystem.h"
#include "tgt/vector.h"

namespace voreen {

const std::string SimilarityDataSave::loggerCat_("voreen.viscontest2018.SimilarityDataSave");

SimilarityDataSave::SimilarityDataSave()
    : Processor()
    // ports
    , ensembleInport_(Port::INPORT, "ensembleInport", "Ensemble Input", false)
    // properties
    , filenameProp_("filenameprop", "Save File as", "Select file...", VoreenApplication::app()->getUserDataPath(), "similarity data (*.sd)", FileDialogProperty::DIRECTORY, Processor::INVALID_PATH)
    , saveButton_("saveButton", "Save")
    , temporalResolution_("temporalResolution", "Temporal Resolution", 0.2f, 0.01f, 10000.0f)
    // members
    , saveSimilarityData_(false)
{
    addPort(ensembleInport_);
    ensembleInport_.onChange(MemberFunctionCallback<SimilarityDataSave>(this, &SimilarityDataSave::adjustToEnsemble));

    saveButton_.onChange(MemberFunctionCallback<SimilarityDataSave>(this, &SimilarityDataSave::saveSimilarityData));

    addProperty(filenameProp_);
    addProperty(saveButton_);
    addProperty(temporalResolution_);
}

SimilarityDataSave::~SimilarityDataSave() {
}

bool SimilarityDataSave::isReady() const {
    if (!isInitialized())
        return false;
    // only streamlineport must be ready
    if (!ensembleInport_.isReady())
        return false;

    return true;
}

void SimilarityDataSave::invalidate(int inv) {
    Processor::invalidate(inv);

    if (inv == Processor::INVALID_PATH && isInitialized()) {
        saveSimilarityData_ = true;
    }
}

void SimilarityDataSave::process() {
    if (saveSimilarityData_){
        saveSimilarityData();
        saveSimilarityData_ = false;
    }
}

//---------------------------------------------------------------------------------
//      Callbacks
//---------------------------------------------------------------------------------
void SimilarityDataSave::saveSimilarityData() {
    if (!isInitialized())
        return;
    if (!ensembleInport_.hasData()) {
        LWARNING("no input ensemble");
        return;
    }
    if (filenameProp_.get().empty()) {
        LWARNING("no filename specified");
        return;
    }
    if(!tgt::FileSystem::dirExists(filenameProp_.get())) {
        LWARNING("path is no directory");
        return;
    }

    const EnsembleDataset* ensemble = ensembleInport_.getData();
    float temporalResolution = temporalResolution_.get();
    const int numTimeSteps = static_cast<int>(ensemble->getEndTime() / temporalResolution);
    float progressIncrement = temporalResolution / (ensemble->getCommonChannels().size() * ensemble->getMaxTotalDuration());
    setProgress(0.0f);

    // Create a folder for each scalar.
    std::map<std::string, std::vector<std::vector<float>>> data;
    for(const std::string& channel : ensemble->getCommonChannels()) {
        std::string dir = filenameProp_.get() + "/" + channel;
        if (!tgt::FileSystem::dirExists(dir) && !tgt::FileSystem::createDirectory(dir)) {
            LERROR("Could not create directory: " + dir);
            return;
        }

        std::vector<std::vector<float>>& field = data[channel];

        field.resize(tgt::hmul(ensemble->getDimensions()));

        // Each timestep
        for(int t = 0; t < numTimeSteps; t++) {

            // Reset voxel values.
            for (size_t i = 0; i < field.size(); i++) {
                field[i].resize(ensemble->getRuns().size());

                // Initialize each value with 0, since there might not be a suitable timestep of each run.
                for (size_t j = 0; j < field[i].size(); j++) {
                    field[i][j] = 0.0f;
                }
            }

            tgt::svec3 pos;
            const tgt::vec3& dim = ensemble->getDimensions();
            for (size_t r = 0; r < ensemble->getRuns().size(); r++) {
                for (pos.z = 0; pos.z < dim.z; pos.z++) {
                    for (pos.y = 0; pos.y < dim.y; pos.y++) {
                        for (pos.x = 0; pos.x < dim.x; pos.x++) {

                            const EnsembleDataset::Run& run = ensemble->getRuns()[r];

                            // TODO: Optimize: do not seach linearly from the beginning each time entering this loop
                            for (size_t ts = 0; ts < run.timeSteps_.size(); ts++) {
                                if (t*temporalResolution >= run.timeSteps_[ts].time_ && t*temporalResolution < run.timeSteps_[ts].time_ + run.timeSteps_[ts].duration_) {
                                    const VolumeBase* volume = run.timeSteps_[ts].channels_.at(channel);
                                    const VolumeRAM_Float* volumeData = dynamic_cast<const VolumeRAM_Float*>(volume->getRepresentation<VolumeRAM>());
                                    tgtAssert(volumeData, "Could not create VolumeRAM representation");

                                    size_t idx = VolumeRAM_Float::calcPos(dim, pos);
                                    field[idx][r] = ensemble->pickSample(volumeData, volume->getSpacing(), pos);
                                    // We found a suitable time step.
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            // Write file.
            try {

                std::string name = itos(t, static_cast<int>(std::to_string(numTimeSteps).length()));
                std::string path = dir + "/" + name + ".sd";

                std::ofstream outFile;
                outFile.open(path);

                JsonSerializer serializer;
                Serializer s(serializer);
                s.serialize("rawSimilarityData", field);
                s.serialize("dimensions", tgt::vec3(ensemble->getDimensions()));
                serializer.write(outFile, false);

                outFile.close();
            }
            catch(tgt::Exception& e) {
                LERROR(e.what());
            }

            setProgress(getProgress() + progressIncrement);
        }
    }

    setProgress(1.0f);
}

void SimilarityDataSave::adjustToEnsemble() {

    temporalResolution_.setReadOnlyFlag(true);

    if(!ensembleInport_.hasData())
        return;

    temporalResolution_.setMinValue(0.01f);
    temporalResolution_.setMaxValue(ensembleInport_.getData()->getMaxTimeStepDuration());
    temporalResolution_.set(ensembleInport_.getData()->getMinTimeStepDuration());

    temporalResolution_.setReadOnlyFlag(false);
}

}   // namespace
