/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "voreen/core/io/volumereader.h"

#include "../utils/ensemblehash.h"
#include "../utils/utils.h"

namespace voreen {

// Deprecated names:
static const char* SCALAR_FIELD_NAME = "Scalar";
static const char* NAME_FIELD_NAME = "name";
static const char* SIMULATED_TIME_NAME = "simulated_time";

const std::string EnsembleDataSource::loggerCat_("voreen.ensembleanalysis.EnsembleDataSource");

EnsembleDataSource::EnsembleDataSource()
    : Processor()
    , outport_(Port::OUTPORT, "ensembledataset", "EnsembleDataset Output", false)
    , ensemblePath_("ensemblepath", "Ensemble Path", "Select Ensemble root folder", "", "", FileDialogProperty::DIRECTORY, Processor::VALID)
    , loadingStrategy_("loadingStrategy", "Loading Strategy", Processor::VALID)
    , loadDatasetButton_("loadDataset", "Load Dataset")
    , memberProgress_("memberProgress", "Members loaded")
    , timeStepProgress_("timeStepProgress", "Time Steps loaded")
    , loadedMembers_("loadedMembers", "Loaded Members", 5)
    , printEnsemble_("printEnsemble", "Print Ensemble", "Print Ensemble", "", "HTML (*.html)", FileDialogProperty::SAVE_FILE, Processor::VALID)
    , colorMap_("colorMap", "Color Map", ColorMap::createPET(), Processor::VALID)
    , overrideTime_("overrideTime", "Override Time", false, Processor::VALID, Property::LOD_ADVANCED)
    , overrideFieldName_("overrideFieldName", "Override Field Name", false, Processor::VALID, Property::LOD_ADVANCED)
    , showProgressDialog_("showProgressDialog", "Show Progress Dialog", true, Processor::VALID, Property::LOD_DEBUG)
    , hash_("hash", "Hash", "", Processor::VALID, Property::LOD_DEBUG)
    , clearCache_("clearCache", "Clear Cache", Processor::VALID, Property::LOD_ADVANCED)
{
    addPort(outport_);
    addProperty(ensemblePath_);
    ON_CHANGE_LAMBDA(ensemblePath_, [this] {
        if(ensemblePath_.isFileWatchEnabled()) {
            buildEnsembleDataset();
        }
    });
    addProperty(loadingStrategy_);
    loadingStrategy_.addOption("manual", "Manual");
    loadingStrategy_.addOption("full", "Full");
    loadingStrategy_.addOption("cached", "Cached");
    addProperty(loadDatasetButton_);
    ON_CHANGE(loadDatasetButton_, EnsembleDataSource, loadEnsembleDataset);
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
    // This is an alternative color map used for the deep water asteroid impact ensemble.
    colors.emplace_back(tgt::Color(230,  25,  75, 255)/255.0f);
    colors.emplace_back(tgt::Color( 60, 180,  75, 255)/255.0f);
    colors.emplace_back(tgt::Color(255, 225,  25, 255)/255.0f);
    colors.emplace_back(tgt::Color(  0, 130, 200, 255)/255.0f);
    colors.emplace_back(tgt::Color(245, 130,  48, 255)/255.0f);
    colors.emplace_back(tgt::Color(145,  30, 180, 255)/255.0f);
    colors.emplace_back(tgt::Color( 70, 240, 240, 255)/255.0f);
    colors.emplace_back(tgt::Color( 0, 0, 0, 255)/255.0f); // Needs to be added since colormap iterators are implemented weirdly.
    colorMap_.set(ColorMap::createFromVector(colors));

    addProperty(overrideTime_);
    addProperty(overrideFieldName_);
    addProperty(showProgressDialog_);
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
    // Just set the data, because connecting another port would require to reload the data otherwise.
    outport_.setData(output_.get(), false);
}

void EnsembleDataSource::deinitialize() {
    clearEnsembleDataset();
    Processor::deinitialize();
}

void EnsembleDataSource::deserialize(Deserializer& s) {
    Processor::deserialize(s);
    if(loadingStrategy_.get() != "manual") {
        loadEnsembleDataset();
    }
}

void EnsembleDataSource::clearEnsembleDataset() {
    outport_.clear();
    output_.reset(); // Important: clear the output before deleting volumes!
    volumes_.clear();
    setProgress(0.0f);
    timeStepProgress_.setProgress(0.0f);
    loadedMembers_.reset();
    hash_.set("");
}

void EnsembleDataSource::buildEnsembleDataset() {

    // Delete old data.
    clearEnsembleDataset();

    std::unique_ptr<EnsembleDataset> dataset(new EnsembleDataset());

    std::vector<std::string> members = tgt::FileSystem::listSubDirectories(ensemblePath_.get(), true);
    float progressPerMember = 1.0f / members.size();

    EnsembleVolumeReaderPopulator populator;
    ColorMap::InterpolationIterator colorIter = colorMap_.get().getInterpolationIterator(members.size());

    for(size_t i=0; i<members.size(); i++) {

        const std::string& member = members[i];

        std::unique_ptr<ProgressBar> progressDialog;
        if (showProgressDialog_.get() && VoreenApplication::app()) {
            progressDialog.reset(VoreenApplication::app()->createProgressDialog());
        }
        if (progressDialog) {
            progressDialog->setTitle("Loading Member " + std::to_string(i+1) + "/" + std::to_string(members.size()));
            progressDialog->setProgressMessage("Loading " + member + " ...");
            progressDialog->show();
            progressDialog->setProgress(0.f);
            progressDialog->forceUpdate();
        }

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

            VolumeReader* reader = populator.getVolumeReader(url);
            if (!reader) {
                LERROR("No suitable reader found for " << fileName);
                continue;
            }

            std::map<std::string, const VolumeBase*> volumeData;
            float time = 0.0f;
            float duration = 0.0f;
            bool timeIsSet = false;

            std::vector<VolumeURL> subURLs = reader->listVolumes(url);
            for(size_t k = 0; k<subURLs.size(); k++) {
                const VolumeURL& subURL = subURLs[k];

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
                    if (volumeHandle->hasMetaData(SIMULATED_TIME_NAME)) { // deprecated
                        currentTime = volumeHandle->getMetaDataValue<FloatMetaData>(SIMULATED_TIME_NAME, 0.0f);
                    }
                    else if (volumeHandle->hasMetaData(VolumeBase::META_DATA_NAME_TIMESTEP)) {
                        currentTime = volumeHandle->getTimestep();
                    }
                    else {
                        currentTime = 1.0f * timeSteps.size();
                        LWARNING("No time step information found for time step " << timeSteps.size() << " of member " << member);
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
                    LWARNING("Time stamp not equal field-wise for " << subURL.getURL());
                }

                std::string fieldName;
                if(!overrideFieldName_.get()) {
                    if (volumeHandle->hasMetaData(NAME_FIELD_NAME)) { // deprecated
                        fieldName = volumeHandle->getMetaData(NAME_FIELD_NAME)->toString();
                    } else if (volumeHandle->hasMetaData(SCALAR_FIELD_NAME)) { // deprecated
                        fieldName = volumeHandle->getMetaData(SCALAR_FIELD_NAME)->toString();
                    } else if (volumeHandle->hasMetaData(VolumeBase::META_DATA_NAME_MODALITY)) {
                        fieldName = volumeHandle->getModality().getName();
                    } else {
                        fieldName = "field_" + std::to_string(k);
                        LWARNING("No field name information found for field " << k << " in " << subURL.getURL());
                    }
                }
                else {
                    fieldName = "field_" + std::to_string(k);
                }

                volumeData[fieldName] = volumeHandle.get();

                // Ownership remains.
                volumes_.push_back(std::move(volumeHandle));
            }

            // Calculate duration of the current timeStep.
            // Note that the last time step has a duration of 0.
            // TODO: duration might change after sorting step.
            if (!timeSteps.empty()) {
                duration = time - timeSteps.back().getTime();
            }

            timeSteps.emplace_back(TimeStep{volumeData, time, duration});

            // Update progress bar.
            float progress = std::min(timeStepProgress_.getProgress() + progressPerTimeStep, 1.0f);
            timeStepProgress_.setProgress(progress);
            if(progressDialog) {
                progressDialog->setProgress(progress);
            }
        }

        auto timeStepCompare = [] (const TimeStep& t1, const TimeStep& t2) {
            return t1.getTime() < t2.getTime();
        };
        std::stable_sort(timeSteps.begin(), timeSteps.end(), timeStepCompare);

        // Set color.
        tgt::Color color = *colorIter;
        dataset->addMember({member, color.xyz(), timeSteps});
        ++colorIter;

        // Update progress bar.
        setProgress(getProgress() + progressPerMember);
        if(progressDialog) {
            progressDialog->hide();
        }
    }

    hash_.set(EnsembleHash(*dataset).getHash());
    output_ = std::move(dataset);
    updateTable();

    if(loadingStrategy_.get() == "cached" && VoreenApplication::app()->useCaching()) {
        tgt::FileSystem::createDirectoryRecursive(getCachePath());

        std::ofstream outFile;
        outFile.open(getEnsembleCachePath());
        if(outFile.good()) {
            try {
                JsonSerializer s;
                Serializer serializer(s);
                serializer.serialize("hash", hash_.get());
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

void EnsembleDataSource::loadEnsembleDataset() {

    if(ensemblePath_.get().empty()) {
        LWARNING("No path specified");
        return;
    }

    if(loadingStrategy_.get() == "cached" && VoreenApplication::app()->useCaching()) {

        std::ifstream inFile;
        inFile.open(getEnsembleCachePath());

        if(inFile.good()) {

            std::unique_ptr<EnsembleDataset> ensemble(new EnsembleDataset);
            try {
                std::string hash;

                JsonDeserializer d;
                d.read(inFile, true);
                Deserializer deserializer(d);
                deserializer.deserialize("hash", hash);
                deserializer.deserialize("ensemble", *ensemble);

                clearEnsembleDataset();
                if(hash == EnsembleHash(*ensemble).getHash()) {
                    hash_.set(hash);

                    output_ = std::move(ensemble);
                    updateTable();
                    setProgress(1.0f);
                    timeStepProgress_.setProgress(1.0f);
                    return; // Done.
                }

                LWARNING("Hash mismatch, rebuilding cache...");
            }
            catch (tgt::Exception& e) {
                LWARNING("Loading ensemble from cache failed, rebuilding cache...");
            }
        }
        else if(!ensemblePath_.get().empty()) {
            LWARNING("Could not read cache file, rebuilding cache...");
        }
    }

    // Rebuild ensemble as fallback.
    buildEnsembleDataset();
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

void EnsembleDataSource::updateTable() {
    loadedMembers_.reset();

    if(!output_) {
        return;
    }

    for(const EnsembleMember& member : output_->getMembers()) {
        std::vector<std::string> row(5);
        row[0] = member.getName(); // Name
        row[1] = std::to_string(member.getTimeSteps().size()); // Num Time Steps
        if (!member.getTimeSteps().empty()) {
            row[2] = std::to_string(member.getTimeSteps().front().getTime()); // Start time
            row[3] = std::to_string(member.getTimeSteps().back().getTime()); // End time
            row[4] = std::to_string(member.getTimeSteps().back().getTime() - member.getTimeSteps().front().getTime()); // Duration
        } else {
            row[2] = row[3] = row[4] = "N/A";
        }
        loadedMembers_.addRow(row);
    }
}

std::string EnsembleDataSource::getEnsembleCachePath() const {
    return getCachePath() + "/" + VoreenHash::getHash(ensemblePath_.get()) + ".ensemble";
}

} // namespace
