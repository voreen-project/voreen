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
#include "../utils/colorpool.h"

namespace voreen {

const std::string EnsembleDataSource::SCALAR_FIELD_NAME = "Scalar";
const std::string EnsembleDataSource::SIMULATED_TIME_NAME = "simulated_time";
const std::string EnsembleDataSource::loggerCat_("voreen.viscontest2018.EnsembleDataSource");

EnsembleDataSource::EnsembleDataSource()
    : Processor()
    , ensemblePath_("ensemblepath", "Ensemble Path", "Select Ensemble root folder", "", "", FileDialogProperty::DIRECTORY)
    , loadDatasetButton_("loadDataset", "Load Dataset")
    , runProgress_("runProgress", "Runs loaded")
    , timeStepProgress_("timeStepProgress", "Time Steps loaded")
    , outport_(Port::OUTPORT, "ensembledataset", "EnsembleDataset Output", false)
    , loadedRuns_("loadedRuns", "Loaded Runs", 5)
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

    loadDatasetButton_.onChange(MemberFunctionCallback<EnsembleDataSource>(this, &EnsembleDataSource::buildEnsembleDataset));
}

EnsembleDataSource::~EnsembleDataSource() {
    clearEnsembleDataset();
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

void EnsembleDataSource::clearEnsembleDataset() {
    if(outport_.hasData()) {
        for(const VolumeBase* volume : outport_.getData()->getVolumes())
            delete volume;
        outport_.setData(nullptr);
    }
    setProgress(0.0f);
    timeStepProgress_.setProgress(0.0f);
    loadedRuns_.reset();
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

    for(const std::string& run : runs) {
        std::string runPath = ensemblePath_.get() + "/" + run;
        std::vector<std::string> fileNames = tgt::FileSystem::readDirectory(runPath, true, false);

        timeStepProgress_.setProgress(0.0f);
        float progressPerTimeStep = 1.0f / fileNames.size();

        std::vector<EnsembleDataset::TimeStep> timeSteps;
        for(const std::string& fileName : fileNames) {
            std::string url = runPath + "/" + fileName;
            std::vector<VolumeReader*> readers = populator_.getVolumeSerializer()->getReaders(url);
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
            bool timeSet = false;

            const std::vector<VolumeURL>& subURLs = reader->listVolumes(url);
            for(const VolumeURL& subURL : subURLs) {
                VolumeBase* volumeHandle = reader->read(subURL);
                if(!volumeHandle)
                    break;

                if (!volumeHandle->hasDerivedData<VolumeMinMax>())
                    LWARNING("Volume does not contain min max information - needs to be calculated");

                const FloatMetaData* time = dynamic_cast<const FloatMetaData*>(volumeHandle->getMetaData(SIMULATED_TIME_NAME));
                if(!time) {
                    delete volumeHandle;
                    LERROR("Meta data '" << SIMULATED_TIME_NAME << "' not present for " << subURL.getPath());
                    break;
                }

                if (!timeSet) {
                    timeStep.time_ = time->getValue();
                    timeSet = true;
                }
                else if (timeStep.time_ != time->getValue())
                    LWARNING("Meta data '" << SIMULATED_TIME_NAME << "' not equal channel-wise in t=" << timeSteps.size() << " of run" << run);

                const MetaDataBase* scalar = volumeHandle->getMetaData(SCALAR_FIELD_NAME);
                if(!scalar) {
                    delete volumeHandle;
                    LERROR("Meta data '" << SCALAR_FIELD_NAME << "' not present for " << subURL.getPath());
                    break;
                }

                // Add additional information gained reading the file structure.
                Volume* volume = dynamic_cast<Volume*>(volumeHandle);
                tgtAssert(volume, "volumeHandle must be volume");
                volume->getMetaDataContainer().addMetaData("run_name", new StringMetaData(run));

                timeStep.channels_[scalar->toString()] = volumeHandle;
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
        tgt::vec3 color = ColorPool::getDistinctColor(dataset->getRuns().size());
        dataset->addRun(EnsembleDataset::Run{ run, color, timeSteps });

        // Update progress bar.
        setProgress(getProgress() + progressPerRun);
    }

    outport_.setData(dataset, true);

    timeStepProgress_.setProgress(1.0f);
    setProgress(1.0f);
    loadDatasetButton_.setReadOnlyFlag(false);
}

} // namespace
