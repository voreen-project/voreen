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
const std::string EnsembleDataSource::RUN_NAME = "run_name";
const std::string EnsembleDataSource::FALLBACK_FIELD_NAME = "unnamed";
const std::string EnsembleDataSource::loggerCat_("voreen.ensembleanalysis.EnsembleDataSource");

EnsembleDataSource::EnsembleDataSource()
    : Processor()
    , outport_(Port::OUTPORT, "ensembledataset", "EnsembleDataset Output", false)
    , ensemblePath_("ensemblepath", "Ensemble Path", "Select Ensemble root folder", "", "", FileDialogProperty::DIRECTORY)
    , loadDatasetButton_("loadDataset", "Load Dataset")
    , runProgress_("runProgress", "Runs loaded")
    , timeStepProgress_("timeStepProgress", "Time Steps loaded")
    , loadedRuns_("loadedRuns", "Loaded Runs", 5)
    , colorMap_("colorMap", "Color Map")
    , hash_("hash", "Hash", "", Processor::VALID, Property::LOD_DEBUG)
{
    addPort(outport_);
    addProperty(ensemblePath_);
    addProperty(loadDatasetButton_);
    addProperty(runProgress_);
    addProgressBar(&runProgress_);
    addProperty(timeStepProgress_);

    addProperty(loadedRuns_);
    loadedRuns_.setColumnLabel(0, "Name");
    loadedRuns_.setColumnLabel(1, "Num Time Steps");
    loadedRuns_.setColumnLabel(2, "Start Time");
    loadedRuns_.setColumnLabel(3, "End Time");
    loadedRuns_.setColumnLabel(4, "Duration");

    addProperty(colorMap_);
    std::vector<tgt::Color> colors;
    colors.push_back(tgt::Color(0.0f, 0.0f, 1.0f, 1.0f));
    colors.push_back(tgt::Color(1.0f, 0.0f, 0.0f, 1.0f));
    colors.push_back(tgt::Color(1.0f, 1.0f, 0.0f, 1.0f));
    colorMap_.set(ColorMap::createFromVector(colors));

    addProperty(hash_);
    hash_.setEditable(false);

    loadDatasetButton_.onChange(MemberFunctionCallback<EnsembleDataSource>(this, &EnsembleDataSource::buildEnsembleDataset));
}

EnsembleDataSource::~EnsembleDataSource() {
}

Processor* EnsembleDataSource::create() const {
    return new EnsembleDataSource();
}

void EnsembleDataSource::process() {
    // Nothing to do.
}

void EnsembleDataSource::initialize() {
    Processor::initialize();
    //buildEnsembleDataset();
}

void EnsembleDataSource::deinitialize() {
    clearEnsembleDataset();
    Processor::deinitialize();
}

void EnsembleDataSource::clearEnsembleDataset() {
    outport_.clear(); // Important: clear the output before deleting volumes!
    volumes_.clear();
    setProgress(0.0f);
    timeStepProgress_.setProgress(0.0f);
    loadedRuns_.reset();
    hash_.reset();
}

void EnsembleDataSource::buildEnsembleDataset() {

    // Delete old data.
    clearEnsembleDataset();

    if(ensemblePath_.get().empty())
        return;

    loadDatasetButton_.setReadOnlyFlag(true);

    EnsembleDataset* dataset = new EnsembleDataset();

    std::vector<std::string> runs = tgt::FileSystem::listSubDirectories(ensemblePath_.get(), true);
    float progressPerRun = 1.0f / runs.size();

    VolumeSerializerPopulator populator;
    ColorMap::InterpolationIterator colorIter = colorMap_.get().getInterpolationIterator(runs.size());

    for(const std::string& run : runs) {
        std::string runPath = ensemblePath_.get() + "/" + run;
        std::vector<std::string> fileNames = tgt::FileSystem::readDirectory(runPath, true, false);

        timeStepProgress_.setProgress(0.0f);
        float progressPerTimeStep = 1.0f / fileNames.size();

        std::vector<EnsembleDataset::TimeStep> timeSteps;
        for(const std::string& fileName : fileNames) {

            // Skip raw files. They belong to VVD files or can't be read anyway.
            if(tgt::FileSystem::fileExtension(fileName, true) == "raw") {
                continue;
            }

            std::string url = runPath + "/" + fileName;
            std::vector<VolumeReader*> readers = populator.getVolumeSerializer()->getReaders(url);
            if(readers.empty()) {
                LERROR("No valid volume reader found for " << url);
                break;
            }

            VolumeReader* reader = readers.front();
            tgtAssert(reader, "Reader was null");

            EnsembleDataset::TimeStep timeStep;
            timeStep.path_ = fileName;
            timeStep.time_ = 0.0f;
            timeStep.duration_ = 0.0f;
            bool timeIsSet = false;

            const std::vector<VolumeURL>& subURLs = reader->listVolumes(url);
            for(const VolumeURL& subURL : subURLs) {
                std::unique_ptr<VolumeBase> volumeHandle(reader->read(subURL));
                if(!volumeHandle)
                    break;

                float time = volumeHandle->getTimestep();
                if(volumeHandle->hasMetaData(SIMULATED_TIME_NAME)) {
                    time = dynamic_cast<const FloatMetaData*>(volumeHandle->getMetaData(SIMULATED_TIME_NAME))->getValue();
                }

                if (!timeIsSet) {
                    timeStep.time_ = time;
                    timeIsSet = true;
                }
                else if (timeStep.time_ != time)
                    LWARNING("Meta data '" << SIMULATED_TIME_NAME << "' not equal channel-wise in t=" << timeSteps.size() << " of run" << run);

                std::string name;
                if(volumeHandle->hasMetaData(NAME_FIELD_NAME)) {
                    name = volumeHandle->getMetaData(NAME_FIELD_NAME)->toString();
                }
                else if(volumeHandle->hasMetaData(SCALAR_FIELD_NAME)) {
                    name = volumeHandle->getMetaData(SCALAR_FIELD_NAME)->toString();
                }
                else {
                    name = FALLBACK_FIELD_NAME;
                }

                // Add additional information gained reading the file structure.
                Volume* volume = dynamic_cast<Volume*>(volumeHandle.get());
                tgtAssert(volume, "volumeHandle must be volume");
                volume->getMetaDataContainer().addMetaData(RUN_NAME, new StringMetaData(run));

                timeStep.fieldNames_[name] = volumeHandle.get();

                // Ownership remains.
                volumes_.push_back(std::move(volumeHandle));
            }

            // Calculate duration the current timeStep is valid.
            // Note that the last time step has a duration of 0.
            if (!timeSteps.empty())
                timeSteps.back().duration_ = timeStep.time_ - timeSteps.back().time_;

            timeSteps.push_back(timeStep);

            // Update progress bar.
            timeStepProgress_.setProgress(std::min(timeStepProgress_.getProgress() + progressPerTimeStep, 1.0f));
        }

        // Update overview table.
        std::vector<std::string> row(5);
        row[0] = run; // Name
        row[1] = std::to_string(timeSteps.size()); // Num Time Steps
        if (!timeSteps.empty()) {
            row[2] = std::to_string(timeSteps.front().time_); // Start time
            row[3] = std::to_string(timeSteps.back().time_); // End time
            row[4] = std::to_string(timeSteps.back().time_ - timeSteps.front().time_); // Duration
        }
        else {
            row[2] = row[3] = row[4] = "N/A";
        }
        loadedRuns_.addRow(row);//addRow(run, color);

        // Update dataset.
        tgt::Color color = *colorIter;
        dataset->addRun(EnsembleDataset::Run{ run, color.xyz(), timeSteps });
        ++colorIter;

        // Update progress bar.
        setProgress(getProgress() + progressPerRun);
    }

    hash_.set(EnsembleHash(*dataset).getHash());
    outport_.setData(dataset, true);

    timeStepProgress_.setProgress(1.0f);
    setProgress(1.0f);
    loadDatasetButton_.setReadOnlyFlag(false);
}

} // namespace
