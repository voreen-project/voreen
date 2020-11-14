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

#include "ensembledatasource.h"

#include "voreen/core/datastructures/volume/volumelist.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/io/volumereader.h"
#include "voreen/core/io/volumeserializer.h"

#include "../datastructures/ensembledataset.h"
#include "../utils/ensemblehash.h"

namespace voreen {

const std::string EnsembleDataSource::SCALAR_FIELD_NAME = "Scalar";
const std::string EnsembleDataSource::NAME_FIELD_NAME = "name";
const std::string EnsembleDataSource::SIMULATED_TIME_NAME = "simulated_time";
const std::string EnsembleDataSource::MEMBER_NAME = "member_name";
const std::string EnsembleDataSource::FALLBACK_FIELD_NAME = "unnamed";
const std::string EnsembleDataSource::loggerCat_("voreen.ensembleanalysis.EnsembleDataSource");

EnsembleDataSource::EnsembleDataSource()
    : Processor()
    , outport_(Port::OUTPORT, "ensembledataset", "EnsembleDataset Output", false)
    , ensemblePath_("ensemblepath", "Ensemble Path", "Select Ensemble root folder", "", "", FileDialogProperty::DIRECTORY, Processor::INVALID_PATH)
    , loadingStrategy_("loadingStrategy", "Loading Strategy", Processor::VALID)
    , loadDatasetButton_("loadDataset", "Load Dataset")
    , memberProgress_("memberProgress", "Members loaded")
    , timeStepProgress_("timeStepProgress", "Time Steps loaded")
    , loadedMembers_("loadedMembers", "Loaded Members", 5)
    , printEnsemble_("printEnsemble", "Print Ensemble", "Print Ensemble", "", "HTML (*.html)", FileDialogProperty::SAVE_FILE)
    , colorMap_("colorMap", "Color Map")
    , overrideTime_("overrideTime", "Override Time", false, Processor::VALID, Property::LOD_ADVANCED)
    , hash_("hash", "Hash", "", Processor::VALID, Property::LOD_DEBUG)
    , clearCache_("clearCache", "Clear Cache", Processor::VALID, Property::LOD_ADVANCED)
{
    addPort(outport_);
    addProperty(ensemblePath_);
    addProperty(loadingStrategy_);
    loadingStrategy_.addOption("manual", "Manual");
    loadingStrategy_.addOption("full", "Full");
    loadingStrategy_.addOption("lazy", "Lazy");
    addProperty(loadDatasetButton_);
    ON_CHANGE(loadDatasetButton_, EnsembleDataSource, buildEnsembleDataset);
    addProperty(memberProgress_);
    addProgressBar(&memberProgress_);
    addProperty(timeStepProgress_);

    addProperty(loadedMembers_);
    loadedMembers_.setColumnLabel(0, "Name");
    loadedMembers_.setColumnLabel(1, "Num Time Steps");
    loadedMembers_.setColumnLabel(2, "Start Time");
    loadedMembers_.setColumnLabel(3, "End Time");
    loadedMembers_.setColumnLabel(4, "Duration");

    addProperty(printEnsemble_);
    ON_CHANGE(printEnsemble_, EnsembleDataSource, printEnsembleDataset);

    addProperty(colorMap_);
    std::vector<tgt::Color> colors;
    colors.push_back(tgt::Color(0.0f, 0.0f, 1.0f, 1.0f));
    colors.push_back(tgt::Color(1.0f, 0.0f, 0.0f, 1.0f));
    colors.push_back(tgt::Color(1.0f, 1.0f, 0.0f, 1.0f));
    colorMap_.set(ColorMap::createFromVector(colors));

    addProperty(overrideTime_);
    addProperty(hash_);
    hash_.setEditable(false);

    addProperty(clearCache_);
    ON_CHANGE_LAMBDA(clearCache_, [this] {
        tgt::FileSystem::deleteDirectoryRecursive(getCachePath());
    });
}

Processor* EnsembleDataSource::create() const {
    return new EnsembleDataSource();
}

void EnsembleDataSource::process() {

    // Reload whole ensemble, if file watching was enabled and some file changed.
    if(invalidationLevel_ >= INVALID_PATH && ensemblePath_.isFileWatchEnabled()) {
        buildEnsembleDataset();
    }

    // Just set the data, because connecting another port would require to reload the data otherwise.
    // This also enables file watching.
    outport_.setData(output_.get(), false);
}

void EnsembleDataSource::deinitialize() {
    clearEnsembleDataset();
    Processor::deinitialize();
}

void EnsembleDataSource::deserialize(Deserializer& s) {
    Processor::deserialize(s);

    if (loadingStrategy_.get() == "full") {
        buildEnsembleDataset();
    }
    else if(loadingStrategy_.get() == "lazy" && VoreenApplication::app()->useCaching()) {

        std::string filename = getCachePath() + "/" + hash_.get() + ".ensemble";

        std::ifstream inFile;
        inFile.open(filename);

        if(inFile.good()) {
            output_.reset(new EnsembleDataset());
            try {
                JsonDeserializer d;
                d.read(inFile, true);
                Deserializer deserializer(d);
                deserializer.deserialize("ensemble", *output_);
            }
            catch (tgt::Exception& e) {
                LWARNING("Loading ensemble from cache failed: " << e.what());
                output_.reset();
            }
            inFile.close();
        }
        else if(!ensemblePath_.get().empty()) {
            LWARNING("Could not read cache file, rebuilding cache...");
            buildEnsembleDataset();
        }
    }
}

void EnsembleDataSource::clearEnsembleDataset() {
    outport_.clear();
    output_.reset(); // Important: clear the output before deleting volumes!
    volumes_.clear();
    setProgress(0.0f);
    timeStepProgress_.setProgress(0.0f);
    loadedMembers_.reset();
    hash_.reset();
}

void EnsembleDataSource::buildEnsembleDataset() {

    // Delete old data.
    clearEnsembleDataset();

    if(ensemblePath_.get().empty())
        return;

    std::unique_ptr<EnsembleDataset> dataset(new EnsembleDataset());

    std::vector<std::string> members = tgt::FileSystem::listSubDirectories(ensemblePath_.get(), true);
    float progressPerMember = 1.0f / members.size();

    VolumeSerializerPopulator populator;
    ColorMap::InterpolationIterator colorIter = colorMap_.get().getInterpolationIterator(members.size());

    for(const std::string& member : members) {
        std::string memberPath = ensemblePath_.get() + "/" + member;
        std::vector<std::string> fileNames = tgt::FileSystem::readDirectory(memberPath, true, false);

        timeStepProgress_.setProgress(0.0f);
        float progressPerTimeStep = 1.0f / fileNames.size();

        std::vector<TimeStep> timeSteps;
        for (const std::string& fileName : fileNames) {

            // Skip raw files. They belong to VVD files or can't be read anyway.
            std::string ext = tgt::FileSystem::fileExtension(fileName, true);
            if (ext.empty() || ext == "raw") {
                continue;
            }

            std::string url = memberPath + "/" + fileName;

            VolumeReader* reader = EnsembleDataset::getVolumeReader(url);
            if (!reader) {
                LERROR("No suitable reader found for " << fileName);
                continue;
            }

            std::map<std::string, const VolumeBase*> volumeData;
            float time = 0.0f;
            float duration = 0.0f;
            bool timeIsSet = false;

            std::vector<VolumeURL> subURLs = reader->listVolumes(url);
            for(const VolumeURL& subURL : subURLs) {
                std::unique_ptr<VolumeBase> volumeHandle;
                try {
                    volumeHandle.reset(reader->read(subURL));
                }
                catch (tgt::IOException& e) {
                    LERROR("Error reading " << subURL.getURL() << ": " << e.what());
                }

                if(!volumeHandle)
                    break;

                float currentTime = 0.0f;
                if(!overrideTime_.get()) {
                    if (volumeHandle->hasMetaData(VolumeBase::META_DATA_NAME_TIMESTEP)) {
                        currentTime = volumeHandle->getTimestep();
                    } else if (volumeHandle->hasMetaData(SIMULATED_TIME_NAME)) {
                        currentTime = volumeHandle->getMetaDataValue<FloatMetaData>(SIMULATED_TIME_NAME, 0.0f);
                    }
                    else {
                        currentTime = 1.0f * timeSteps.size();
                        LWARNING("Actual time information not found for time step " << timeSteps.size() << " of member " << member);
                    }
                }
                else {
                    currentTime = 1.0f * timeSteps.size();
                }

                if (!timeIsSet) {
                    time = currentTime;
                    timeIsSet = true;
                }
                else if (currentTime != time) {
                    LWARNING("Time stamp not equal channel-wise for t=" << timeSteps.size() << " of member " << member);
                }

                std::string fieldName;
                if(volumeHandle->hasMetaData(NAME_FIELD_NAME)) {
                    fieldName = volumeHandle->getMetaData(NAME_FIELD_NAME)->toString();
                }
                else if(volumeHandle->hasMetaData(SCALAR_FIELD_NAME)) {
                    fieldName = volumeHandle->getMetaData(SCALAR_FIELD_NAME)->toString();
                }
                else {
                    fieldName = FALLBACK_FIELD_NAME;
                }

                // Add additional information gained reading the file structure.
                Volume* volume = dynamic_cast<Volume*>(volumeHandle.get());
                tgtAssert(volume, "volumeHandle must be volume");
                volume->setMetaDataValue<StringMetaData>(MEMBER_NAME, member);

                volumeData[fieldName] = volumeHandle.get();

                // Ownership remains.
                volumes_.push_back(std::move(volumeHandle));
            }

            // Calculate duration the current timeStep is valid.
            // Note that the last time step has a duration of 0.
            if (!timeSteps.empty())
                duration = time - timeSteps.back().getTime();

            timeSteps.emplace_back(TimeStep{volumeData, time, duration});

            // Update progress bar.
            timeStepProgress_.setProgress(std::min(timeStepProgress_.getProgress() + progressPerTimeStep, 1.0f));
        }


        auto timeStepCompare = [] (const TimeStep& t1, const TimeStep& t2) {
            return t1.getTime() < t2.getTime();
        };
        std::stable_sort(timeSteps.begin(), timeSteps.end(), timeStepCompare);

        std::vector<TimeStep> timeStepsVec(std::make_move_iterator(timeSteps.begin()), std::make_move_iterator(timeSteps.end()));

        // Update overview table.
        std::vector<std::string> row(5);
        row[0] = member; // Name
        row[1] = std::to_string(timeSteps.size()); // Num Time Steps
        if (!timeSteps.empty()) {
            row[2] = std::to_string(timeStepsVec.front().getTime()); // Start time
            row[3] = std::to_string(timeStepsVec.back().getTime()); // End time
            row[4] = std::to_string(timeStepsVec.back().getTime() - timeStepsVec.front().getTime()); // Duration
        }
        else {
            row[2] = row[3] = row[4] = "N/A";
        }
        loadedMembers_.addRow(row);

        // Update dataset.
        tgt::Color color = *colorIter;
        dataset->addMember({member, color.xyz(), timeStepsVec});
        ++colorIter;

        // Update progress bar.
        setProgress(getProgress() + progressPerMember);
    }

    hash_.set(EnsembleHash(*dataset).getHash());
    output_ = std::move(dataset);

    if(loadingStrategy_.get() == "lazy" && VoreenApplication::app()->useCaching()) {
        tgt::FileSystem::createDirectoryRecursive(getCachePath());
        std::string filename = getCachePath() + "/" + hash_.get() + ".ensemble";

        std::ofstream outFile;
        outFile.open(filename);
        if(outFile.good()) {
            try {
                JsonSerializer s;
                Serializer serializer(s);
                serializer.serialize("ensemble", *output_);
                s.write(outFile, false, true);
            }
            catch (tgt::Exception& e) {
                LWARNING("Storing ensemble meta data to cache failed: " << e.what());
            }
            outFile.close();
        }
        else {
            LWARNING("Could not write cache file");
        }
    }

    timeStepProgress_.setProgress(1.0f);
    setProgress(1.0f);
}

void EnsembleDataSource::printEnsembleDataset() {

    if(!output_) {
        LWARNING("No ensemble loaded");
        return;
    }

    std::fstream file(printEnsemble_.get(), std::ios::out);
    file << output_->toHTML();
    if (!file.good()) {
        LERROR("Could not write " << printEnsemble_.get() << " file");
    }
}

} // namespace
