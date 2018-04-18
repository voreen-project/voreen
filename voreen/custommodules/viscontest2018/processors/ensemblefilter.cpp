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

#include "ensemblefilter.h"

#include "voreen/core/datastructures/callback/lambdacallback.h"

#include "voreen/core/properties/boundingboxproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"

#include "../properties/stringlistproperty.h"

namespace voreen {

//----------------------------------------
// Filter : Run
//----------------------------------------

class FilterRun : public Filter {
public:
    FilterRun()
        : runs_("runs", "Selected Runs")
    {
        runs_.setDescription("Selects multiple runs from the ensemble data.");
    }

    Property& getProperty() {
        return runs_;
    }

    EnsembleDataset* applyFilter(const EnsembleDataset* ensemble) {
        
        EnsembleDataset* dataset = new EnsembleDataset();
        dataset->setRoi(ensemble->getRoi());

        for (int row : runs_.getSelectedRowIndices()) {
            EnsembleDataset::Run run = ensemble->getRuns()[row];
            dataset->addRun(run);
        }

        return dataset;
    }

    void adjustToEnsemble(const EnsembleDataset* ensemble) {
        
        // Reset range.
        runs_.reset();

        if (ensemble) {
            // Adjust range to data.
            std::vector<int> selectedRunIndices;
            for (const EnsembleDataset::Run& run : ensemble->getRuns()) {
                runs_.addRow(run.name_, run.color_);
                selectedRunIndices.push_back(selectedRunIndices.size());
            }
            runs_.setSelectedRowIndices(selectedRunIndices);
        }
    }

private:
    StringListProperty runs_;
};

//----------------------------------------
// Filter : Time Step
//----------------------------------------
class FilterTimeStep : public Filter {
public:
    FilterTimeStep()
        : timeSteps_("timeSteps", "Selected Time Steps", tgt::ivec2(-1, -1), -1, -1)
    {
        timeSteps_.setDescription("Selects a range from time steps from the ensemble data.");
    }

    Property& getProperty() {
        return timeSteps_;
    }

    EnsembleDataset* applyFilter(const EnsembleDataset* ensemble) {
        
        EnsembleDataset* dataset = new EnsembleDataset();
        dataset->setRoi(ensemble->getRoi());

        for (const EnsembleDataset::Run& run : ensemble->getRuns()) {
            if (run.timeSteps_.empty())
                continue;

            std::vector<EnsembleDataset::TimeStep> timeSteps;
            int max = std::min<int>(run.timeSteps_.size() - 1, timeSteps_.get().y);
            for (int i = timeSteps_.get().x; i <= max; i++) {
                timeSteps.push_back(run.timeSteps_[i]);
            }

            dataset->addRun(EnsembleDataset::Run{ run.name_, run.color_, timeSteps });
        }

        return dataset;
    }

    void adjustToEnsemble(const EnsembleDataset* ensemble) {
        // Reset range.
        timeSteps_.setMinValue(-1);
        timeSteps_.setMaxValue(-1);
        
        // Adjust range to data.
        if (ensemble && ensemble->getMaxNumTimeSteps() > 0) {
            timeSteps_.setMinValue(0);
            timeSteps_.setMaxValue(ensemble->getMaxNumTimeSteps() - 1);
            timeSteps_.set(tgt::ivec2(0, ensemble->getMaxNumTimeSteps() - 1));
        }
    }

private:

    IntIntervalProperty timeSteps_;
};

//----------------------------------------
// Filter : Time Interval
//----------------------------------------

class FilterTimeInterval : public Filter {
public:
    FilterTimeInterval()
        : timeInterval_("timeInterval", "Selected Time Interval", tgt::vec2(0.0f, 0.0f), 0.0f, 0.0f)
    {
        timeInterval_.setDescription("Selects all time steps within the configured interval from the ensemble data.");
    }

    Property& getProperty() {
        return timeInterval_;
    }

    EnsembleDataset* applyFilter(const EnsembleDataset* ensemble) {

        EnsembleDataset* dataset = new EnsembleDataset();
        dataset->setRoi(ensemble->getRoi());

        for (const EnsembleDataset::Run& run : ensemble->getRuns()) {
            std::vector<EnsembleDataset::TimeStep> timeSteps;
            for (size_t i = 0; i < run.timeSteps_.size(); i++) {
                if (run.timeSteps_[i].time_ > timeInterval_.get().y)
                    break;
                if (run.timeSteps_[i].time_ >= timeInterval_.get().x)
                    timeSteps.push_back(run.timeSteps_[i]);
            }

            dataset->addRun(EnsembleDataset::Run{ run.name_, run.color_, timeSteps });
        }

        return dataset;
    }

    void adjustToEnsemble(const EnsembleDataset* ensemble) {
        if (!ensemble) {
            // Reset range.
            timeInterval_.setMinValue(0.0f);
            timeInterval_.setMaxValue(0.0f);
        }
        else {
            // Adjust range to data.
            timeInterval_.setMinValue(ensemble->getStartTime());
            timeInterval_.setMaxValue(ensemble->getEndTime());
            timeInterval_.set(tgt::vec2(ensemble->getStartTime(), ensemble->getEndTime()));
        }
    }

private:

    FloatIntervalProperty timeInterval_;
};

//----------------------------------------
// Filter : Channel
//----------------------------------------

class FilterChannel : public Filter {
public:
    FilterChannel()
        : channels_("channel", "Selected Channel")
    {
        channels_.setDescription("Selects a single channel from the ensemble data.");
    }

    Property& getProperty() {
        return channels_;
    }

    EnsembleDataset* applyFilter(const EnsembleDataset* ensemble) {

        EnsembleDataset* dataset = new EnsembleDataset();
        dataset->setRoi(ensemble->getRoi());

        for (const EnsembleDataset::Run& run : ensemble->getRuns()) {
            std::vector<EnsembleDataset::TimeStep> timesteps;
            for (const EnsembleDataset::TimeStep& timestep : run.timeSteps_) {
                EnsembleDataset::TimeStep filteredTimeStep = timestep;

                std::map<std::string, const VolumeBase*> filteredChannels;
                filteredChannels[channels_.getValue()] = timestep.channels_.at(channels_.getValue());
                filteredTimeStep.channels_ = filteredChannels;

                timesteps.push_back(filteredTimeStep);
            }
            dataset->addRun(EnsembleDataset::Run{ run.name_, run.color_, timesteps });
        }

        return dataset;
    }

    void adjustToEnsemble(const EnsembleDataset* ensemble) {
        
        channels_.setOptions(std::deque<Option<std::string>>());

        if (ensemble) {
            for (const std::string& channel : ensemble->getCommonChannels()) {
                channels_.addOption(channel, channel, channel);
            }
        }
    }

private:

    OptionProperty<std::string> channels_;
};

//----------------------------------------
// Filter : ROI
//----------------------------------------

class FilterROI : public Filter {
public:
    FilterROI()
        : regionOfInterest_("roi", "Region of interest")
    {
        regionOfInterest_.setDescription("Modifies the region of interest (ROI) of the ensemble dataset.");
    }

    Property& getProperty() {
        return regionOfInterest_;
    }

    EnsembleDataset* applyFilter(const EnsembleDataset* ensemble) {
        // Clone input data and adjust roi.
        EnsembleDataset* dataset = new EnsembleDataset(ensemble);
        dataset->setRoi(regionOfInterest_.get());
        return dataset;
    }

    void adjustToEnsemble(const EnsembleDataset* ensemble) {
        if (!ensemble) {
            regionOfInterest_.setReadOnlyFlag(true);
            regionOfInterest_.setMinValue(tgt::ivec3::zero);
            regionOfInterest_.setMaxValue(tgt::ivec3::zero);
        }
        else {
            regionOfInterest_.setMinValue(ensemble->getRoi().getLLF());
            regionOfInterest_.setMaxValue(ensemble->getRoi().getURB());
            regionOfInterest_.set(ensemble->getRoi());
            regionOfInterest_.setReadOnlyFlag(false);
        }
    }

private:

    IntBoundingBoxProperty regionOfInterest_;
};

//----------------------------------------
// EnsembleFilter
//----------------------------------------

EnsembleFilter::EnsembleFilter()
    : Processor()
    , ensembleInport_(Port::INPORT, "ensembledatastructurein", "Ensemble Datastructure Input", false)
    , ensembleOutport_(Port::OUTPORT, "ensembledatastructureout", "Ensemble Datastructure Output", false)
    , needsProcess_(false)
{
    addPort(ensembleInport_);
    ON_CHANGE(ensembleInport_, EnsembleFilter, adjustToEnsemble);
    addPort(ensembleOutport_);

    addFilter(new FilterRun());
    //addFilter(new FilterTimeStep());
    addFilter(new FilterTimeInterval());
    addFilter(new FilterChannel());
    addFilter(new FilterROI());
}

EnsembleFilter::~EnsembleFilter() {
    for (Filter* filter : filters_)
        delete filter;
}

Processor* EnsembleFilter::create() const {
    return new EnsembleFilter();
}

void EnsembleFilter::process() {
    if(needsProcess_) {
        applyFilter();
        needsProcess_ = false;
    }
}

void EnsembleFilter::invalidate(int inv) {
    Processor::invalidate(inv);

    if (inv == Processor::INVALID_RESULT && isInitialized())
        needsProcess_ = true;
}

void EnsembleFilter::addFilter(Filter* filter) {
    filters_.push_back(filter);
    addProperty(filter->getProperty());
    ON_CHANGE_LAMBDA(filter->getProperty(), [this] { invalidate(Processor::INVALID_RESULT); });
}

void EnsembleFilter::adjustToEnsemble() {
    for (Filter* filter : filters_)
        filter->adjustToEnsemble(ensembleInport_.getData());
}

void EnsembleFilter::applyFilter() {

    ensembleOutport_.setData(nullptr);

    if (ensembleInport_.hasData()) {
        std::unique_ptr<EnsembleDataset> ensemble(new EnsembleDataset(ensembleInport_.getData()));
        for (Filter* filter : filters_) {
            ensemble.reset(filter->applyFilter(ensemble.get()));
        }

        ensembleOutport_.setData(ensemble.release(), true);
    }
}

} // namespace
