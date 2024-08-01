/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "timeserieslistcreator.h"

#include <limits>
#include <random>

namespace voreen {
    const std::string TimeSeriesListCreator::loggerCat_("voreen.scivis2021.TimeSeriesListCreator");

    TimeSeriesListCreator::TimeSeriesListCreator()
        : AsyncComputeProcessor<TimeSeriesListCreatorInput, TimeSeriesListCreatorOutput>()
        , ensembleInport_(Port::INPORT, "ensemble.in", "Ensemble Inport", false)
        , maskInport_(Port::INPORT, "maskInport", "Mask Input (optional)")
        , timeSeriesOutport_(Port::OUTPORT, "time.series.list.out", "TimeSeries Output", false)
        , numberOfSamples_("sample.number", "Number of samples", 1, 0, 1000)
        , sampleOnSlice_("sampleOnSlice", "Sample on slice")
        , depth_("depth", "Depth", 660, 0, 2886)
        , fields_("fields", "Included fields")
    {
        addPort(ensembleInport_);
        addPort(maskInport_);
        addPort(timeSeriesOutport_);

        ensembleInport_.onChange(MemberFunctionCallback<TimeSeriesListCreator>(this, &TimeSeriesListCreator::onInportChange));

        addProperty(numberOfSamples_);
        addProperty(sampleOnSlice_);
        ON_CHANGE_LAMBDA(sampleOnSlice_, [this] {
            depth_.setVisibleFlag(sampleOnSlice_.get());
        });
        addProperty(depth_);
        depth_.setVisibleFlag(sampleOnSlice_.get());
        addProperty(fields_);
    }
    TimeSeriesListCreator::~TimeSeriesListCreator() {

    }

    Processor* TimeSeriesListCreator::create() const {
        return new TimeSeriesListCreator();
    }

    bool TimeSeriesListCreator::isReady() const {
        if (!isInitialized()) {
            setNotReadyErrorMessage("Not initialized");
            return false;
        }

        if (!ensembleInport_.isReady()) {
            setNotReadyErrorMessage("No input");
            return false;
        }

        //mask is optional

        return true;
    }

    TimeSeriesListCreator::ComputeInput TimeSeriesListCreator::prepareComputeInput() {
        auto ensemble = ensembleInport_.getThreadSafeData();
        if (!ensemble) {
            throw InvalidInputException("No input", InvalidInputException::S_WARNING);
        }
        if (ensemble->getMembers().empty()) {
            throw InvalidInputException("Need at least a single run", InvalidInputException::S_ERROR);
        }
        if (fields_.getSelectedRowIndices().size() == 0) {
            throw InvalidInputException("Need at least one selected field", InvalidInputException::S_ERROR);
        }

        timeSeriesOutport_.setData(nullptr);

        float depth = -1;
        if(sampleOnSlice_.get())
        {
            depth = depth_.get();
        }

        return ComputeInput {
            std::move(ensemble),
            numberOfSamples_.get(),
            depth,
            fields_.getSelectedRowIndices()
        };
    }
    TimeSeriesListCreator::ComputeOutput TimeSeriesListCreator::compute(ComputeInput input, ProgressReporter& progressReporter) const {
        auto ensemble = std::move(input.ensemble);
        progressReporter.setProgress(0.01);

        // get some ensemble attributes
        const size_t numMembers = ensemble->getMembers().size();
        std::vector<std::string> fieldNames;
        for(int i : input.fieldIndices)
            fieldNames.push_back(ensemble->getMembers()[0].getTimeSteps().at(0).getFieldNames()[i]);
        //size_t numFields = fieldNames.size();
        float startTime = ensemble->getStartTime();
        float endTime = ensemble->getEndTime();
        std::vector<unsigned int> channelsPerField;
        for (auto& fieldName : fieldNames) {
            channelsPerField.push_back(ensemble->getNumChannels(fieldName));
        }

        // get the number of voxels, assuming all volumes have the same number of voxels
        auto first = ensemble->getMembers()[0].getTimeSteps().at(0).getVolume(fieldNames.at(0));
        tgt::svec3 dimensions = first->getDimensions();
        size_t numVoxels = dimensions.x * dimensions.y * dimensions.z;

        // create time series list
        std::unique_ptr<TimeSeriesList> timeSeriesList(new TimeSeriesList(startTime, endTime,
            fieldNames, channelsPerField));
        timeSeriesList->setVolumeDimensions(dimensions);

        // Set fields
        for (auto& fieldName : fieldNames) {
            timeSeriesList->setRange(fieldName, ensemble->getValueRange(fieldName));
        }

        // for now we assume that we have only one member
        auto member = ensemble->getMembers()[0];
        auto& ensembleTimeSteps = member.getTimeSteps();
        size_t numTimeSteps = ensembleTimeSteps.size();

        // random uniform distribution
        std::default_random_engine generator;
        std::uniform_int_distribution<int> xDistribution(0, dimensions.x - 1);
        std::uniform_int_distribution<int> yDistribution(0, dimensions.y - 1);
        std::uniform_int_distribution<int> zDistribution(0, dimensions.z - 1);
        progressReporter.setProgress(0.03);

        const VolumeBase* seedMask = maskInport_.getThreadSafeData();

        // Create Time Series
        size_t numSamples = input.numberOfSamples;
        for (size_t s = 0; s < numSamples; ++s) {
            // random voxel position
            int rx = xDistribution(generator);
            int ry = yDistribution(generator);
            int rz = zDistribution(generator);
            tgt::svec3 pos(rx, ry, rz);

            if(input.depth > 0) {
                float thickness = ensemble->getBounds().getURB().y - ensemble->getBounds().getLLF().y;
                // Assume same spacing for all members
                float vox =  (thickness - input.depth)/first->getSpacing().y;
                pos.y = int(std::round(vox));
                std::cout << pos.y << std::endl;
            }

            if(seedMask) {
                VolumeRAMRepresentationLock seedMaskLock(seedMask);
                int count = 0;
                while(seedMaskLock->getVoxelNormalized(rx, ry, rz) == 0.0f && count < 1000) {
                    rx = xDistribution(generator);
                    ry = yDistribution(generator);
                    rz = zDistribution(generator);
                    pos = tgt::svec3(rx, ry, rz);
                    count++;
                }
            }

            // create first time series
            TimeSeries timeSeries(pos);

            float currentStart = std::numeric_limits<float>::max();
            float currentEnd = -std::numeric_limits<float>::max();
            timeSeries.setStartTime(currentStart);
            timeSeries.setEndTime(currentEnd);
            timeSeriesList->addTimeSeries(timeSeries);
        }
        progressReporter.setProgress(0.05);
//#pragma omp parallel for
        for (long t = 0; t < static_cast<long>(numTimeSteps); ++t) {
            auto& timeStep = ensembleTimeSteps[t];
//            for(auto& fieldName : fieldNames) {
//                const voreen::VolumeBase* vol = timeStep.getVolume(fieldName);
//                VolumeRAMRepresentationLock lock(vol);
//                voreen::RealWorldMapping rwm = vol->getRealWorldMapping();
//
//                for (size_t s = 0; s < numSamples; ++s) {
//                    TimeSeries* timeSeries = timeSeriesList->getTimeSeriesEditable(s);
//                    timeSeries->setStartTime(std::min(timeSeriesList->getTimeSeries(s).getStartTime(), timeStep.getTime()));
//                    timeSeries->setEndTime(std::max(timeSeriesList->getTimeSeries(s).getEndTime(), timeStep.getTime()));
//                    //                TimeSeriesStep timeSeriesStep(timeStep.getTime(), fieldValues);
////                timeSeries->addTimeSeriesStep(timeSeriesStep);
//                }
//            }
            for (size_t s = 0; s < numSamples; ++s) {
                TimeSeries* timeSeries = timeSeriesList->getTimeSeriesEditable(s);
                timeSeries->setStartTime(std::min(timeSeriesList->getTimeSeries(s).getStartTime(), timeStep.getTime()));
                timeSeries->setEndTime(std::max(timeSeriesList->getTimeSeries(s).getEndTime(), timeStep.getTime()));

                // get fields for current time step
                std::vector<float> fieldValues;
                for (auto& fieldName : fieldNames) {
                    const voreen::VolumeBase* vol = timeStep.getVolume(fieldName);
                    VolumeRAMRepresentationLock lock(vol);
                    voreen::RealWorldMapping rwm = vol->getRealWorldMapping();

                    size_t numChannels = vol->getNumChannels();
                    for (size_t channel = 0; channel < numChannels; ++channel) {
                        float value = rwm.normalizedToRealWorld(lock->getVoxelNormalized(timeSeries->getPosition(), channel));
                        tgt::vec2 range = timeSeriesList->getRange(fieldName);
                        value = (value - range.x)/(range.y-range.x);
                        fieldValues.push_back(value);
                    }
                }
                // add time series step
                TimeSeriesStep timeSeriesStep(timeStep.getTime(), fieldValues);
                timeSeries->addTimeSeriesStep(timeSeriesStep);
            }
            progressReporter.setProgress(0.05+0.95f * t / numTimeSteps);
        }

        return ComputeOutput{
            std::move(timeSeriesList)
        };
    }
    void TimeSeriesListCreator::processComputeOutput(ComputeOutput output) {
        timeSeriesOutport_.setData(output.timeSeriesList.release());
    }

    void TimeSeriesListCreator::onInportChange() {
        if (!ensembleInport_.isReady())
            return;

        // grab ensemble
        auto ensemble = ensembleInport_.getThreadSafeData();
        if (ensemble->getMembers().size() == 0)
            return;

        // grab first timestep
        auto member = ensemble->getMembers()[0];
        auto& ensembleTimeSteps = member.getTimeSteps();
        size_t numTimeSteps = ensembleTimeSteps.size();
        if (numTimeSteps == 0)
            return;

        // grab field names
        auto& timeStep = ensembleTimeSteps[0];
        auto fieldNames = timeStep.getFieldNames();
        if (fieldNames.size() == 0)
            return;

        fields_.reset();
        for(auto& fieldName : fieldNames) {
            fields_.addRow(fieldName);
        }

        // grab volume for first field
        const voreen::VolumeBase* vol = timeStep.getVolume(fieldNames[0]);
        if (vol == nullptr)
            return;

        // get the number of voxels
        tgt::svec3 dimensions = vol->getDimensions();
        size_t numberOfVoxels = dimensions.x * dimensions.y * dimensions.z;
        numberOfSamples_.setMaxValue(numberOfVoxels);
    }

} // namespace
