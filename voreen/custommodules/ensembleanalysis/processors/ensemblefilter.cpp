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
#include "../utils/ensemblehash.h"

namespace voreen {

class Filter {
public:
    virtual Property& getProperty() = 0;
    virtual EnsembleDataset* applyFilter(const EnsembleDataset& ensemble) = 0;
    virtual void adjustToEnsemble(const EnsembleDataset* ensemble) = 0;
};

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

    EnsembleDataset* applyFilter(const EnsembleDataset& ensemble) {
        
        EnsembleDataset* dataset = new EnsembleDataset();

        for (int row : runs_.getSelectedRowIndices()) {
            EnsembleDataset::Run run = ensemble.getRuns()[row];
            dataset->addRun(run);
        }

        dataset->setRoi(ensemble.getRoi());

        return dataset;
    }

    void adjustToEnsemble(const EnsembleDataset* ensemble) {
        
        // Reset range.
        runs_.reset();

        if (ensemble) {
            // Adjust range to data.
            std::vector<int> selectedRunIndices;
            for (const EnsembleDataset::Run& run : ensemble->getRuns()) {
                runs_.addRow(run.name_, ensemble->getColor(selectedRunIndices.size()));
                selectedRunIndices.push_back(static_cast<int>(selectedRunIndices.size()));
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

    EnsembleDataset* applyFilter(const EnsembleDataset& ensemble) {
        
        EnsembleDataset* dataset = new EnsembleDataset();

        for (const EnsembleDataset::Run& run : ensemble.getRuns()) {
            if (run.timeSteps_.empty())
                continue;

            std::vector<EnsembleDataset::TimeStep> timeSteps;
            int max = std::min(static_cast<int>(run.timeSteps_.size()) - 1, timeSteps_.get().y);
            for (int i = timeSteps_.get().x; i <= max; i++) {
                timeSteps.push_back(run.timeSteps_[i]);
            }

            dataset->addRun(EnsembleDataset::Run{ run.name_, run.color_, timeSteps });
        }

        dataset->setRoi(ensemble.getRoi());

        return dataset;
    }

    void adjustToEnsemble(const EnsembleDataset* ensemble) {
        // Reset range.
        timeSteps_.setMinValue(-1);
        timeSteps_.setMaxValue(-1);
        
        // Adjust range to data.
        if (ensemble && ensemble->getMaxNumTimeSteps() > 0) {
            timeSteps_.setMinValue(0);
            timeSteps_.setMaxValue(static_cast<int>(ensemble->getMaxNumTimeSteps()) - 1);
            timeSteps_.set(tgt::ivec2(0, static_cast<int>(ensemble->getMaxNumTimeSteps()) - 1));
        }
    }

private:

    IntIntervalProperty timeSteps_;
};

//----------------------------------------
// Filter : Remove first Time Step
//----------------------------------------
class FilterRemoveFirstTimeStep : public Filter {
public:
    FilterRemoveFirstTimeStep()
            : enableRemoveFirstTimeStep_("enableRemoveFirstTimeStep", "Remove first Time Step", false)
    {
        enableRemoveFirstTimeStep_.setDescription("Removes the first time step of each run.");
    }

    Property& getProperty() {
        return enableRemoveFirstTimeStep_;
    }

    EnsembleDataset* applyFilter(const EnsembleDataset& ensemble) {

        // Clone input, if not enabled.
        if(!enableRemoveFirstTimeStep_.get()) {
            return new EnsembleDataset(ensemble);
        }

        EnsembleDataset* dataset = new EnsembleDataset();

        for (const EnsembleDataset::Run& run : ensemble.getRuns()) {
            if (run.timeSteps_.empty())
                continue;

            std::vector<EnsembleDataset::TimeStep> timeSteps;
            for (size_t i = 1; i < run.timeSteps_.size(); i++) {
                timeSteps.push_back(run.timeSteps_[i]);
            }

            dataset->addRun(EnsembleDataset::Run{ run.name_, run.color_, timeSteps });
        }

        dataset->setRoi(ensemble.getRoi());

        return dataset;
    }

    void adjustToEnsemble(const EnsembleDataset* ensemble) {
    }

private:

    BoolProperty enableRemoveFirstTimeStep_;
};

//----------------------------------------
// Filter : Select last Time Step
//----------------------------------------
class FilterSelectLastTimeStep : public Filter {
public:
    FilterSelectLastTimeStep()
        : enableSelectLastTimeStep_("enableSelectLastTimeStep", "Select last Time Step", false)
    {
        enableSelectLastTimeStep_.setDescription("Selects only the last time step of each run.");
    }

    Property& getProperty() {
        return enableSelectLastTimeStep_;
    }

    EnsembleDataset* applyFilter(const EnsembleDataset& ensemble) {

        // Clone input, if not enabled.
        if(!enableSelectLastTimeStep_.get()) {
            return new EnsembleDataset(ensemble);
        }

        EnsembleDataset* dataset = new EnsembleDataset();

        for (const EnsembleDataset::Run& run : ensemble.getRuns()) {
            if (run.timeSteps_.empty())
                continue;

            std::vector<EnsembleDataset::TimeStep> timeSteps;
            timeSteps.push_back(run.timeSteps_.back());

            dataset->addRun(EnsembleDataset::Run{ run.name_, run.color_, timeSteps });
        }

        dataset->setRoi(ensemble.getRoi());

        return dataset;
    }

    void adjustToEnsemble(const EnsembleDataset* ensemble) {
    }

private:

    BoolProperty enableSelectLastTimeStep_;
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

    EnsembleDataset* applyFilter(const EnsembleDataset& ensemble) {

        EnsembleDataset* dataset = new EnsembleDataset();

        for (const EnsembleDataset::Run& run : ensemble.getRuns()) {
            std::vector<EnsembleDataset::TimeStep> timeSteps;
            for (size_t i = 0; i < run.timeSteps_.size(); i++) {
                if (run.timeSteps_[i].time_ > timeInterval_.get().y)
                    break;
                if (run.timeSteps_[i].time_ >= timeInterval_.get().x)
                    timeSteps.push_back(run.timeSteps_[i]);
            }

            dataset->addRun(EnsembleDataset::Run{ run.name_, run.color_, timeSteps });
        }

        dataset->setRoi(ensemble.getRoi());

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
// Filter : Field
//----------------------------------------

class FilterField : public Filter {
public:
    FilterField()
        : fields_("channel", "Selected Field", Processor::INVALID_RESULT, true)
    {
        fields_.setDescription("Selects a single field from the ensemble data."
                                 "<br>"
                                 "(*) Marks common fields across all runs."
        );
    }

    Property& getProperty() {
        return fields_;
    }

    EnsembleDataset* applyFilter(const EnsembleDataset& ensemble) {

        EnsembleDataset* dataset = new EnsembleDataset();

        for (const EnsembleDataset::Run& run : ensemble.getRuns()) {
            std::vector<EnsembleDataset::TimeStep> timesteps;
            for (const EnsembleDataset::TimeStep& timestep : run.timeSteps_) {
                // Only add time step, if selected field is available.
                if(timestep.fieldNames_.count(fields_.getValue()) > 0) {
                    EnsembleDataset::TimeStep filteredTimeStep = timestep;

                    std::map<std::string, const VolumeBase*> filteredFields;
                    filteredFields[fields_.getValue()] = timestep.fieldNames_.at(fields_.getValue());
                    filteredTimeStep.fieldNames_ = filteredFields;

                    timesteps.push_back(filteredTimeStep);
                }
            }
            dataset->addRun(EnsembleDataset::Run{ run.name_, run.color_, timesteps });
        }

        dataset->setRoi(ensemble.getRoi());

        return dataset;
    }

    void adjustToEnsemble(const EnsembleDataset* ensemble) {
        
        fields_.setOptions(std::deque<Option<std::string>>());

        if (ensemble) {
            const std::set<std::string> common(ensemble->getCommonFieldNames().begin(), ensemble->getCommonFieldNames().end());
            for (const std::string& fieldName : ensemble->getUniqueFieldNames()) {
                bool isCommon = common.find(fieldName) != common.end();
                fields_.addOption(fieldName, fieldName + (isCommon ? " (*)" : ""), fieldName);
            }
        }
    }

private:

    OptionProperty<std::string> fields_;
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

    EnsembleDataset* applyFilter(const EnsembleDataset& ensemble) {
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

    FloatBoundingBoxProperty regionOfInterest_;
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
    //addFilter(new FilterTimeStep()); // replaced by FilterTimeInterval
    addFilter(new FilterTimeInterval());
    addFilter(new FilterRemoveFirstTimeStep());
    addFilter(new FilterSelectLastTimeStep());
    addFilter(new FilterField());
    addFilter(new FilterROI());
}

EnsembleFilter::~EnsembleFilter() {
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
    filters_.push_back(std::unique_ptr<Filter>(filter));
    addProperty(filter->getProperty());
    ON_CHANGE_LAMBDA(filter->getProperty(), [this] { invalidate(Processor::INVALID_RESULT); });
}

void EnsembleFilter::adjustToEnsemble() {
    
    ensembleOutport_.clear();
    
    if(ensembleInport_.hasData()) {
        std::string hash = EnsembleHash(*ensembleInport_.getData()).getHash();
        if(hash != hash_) {
            for (auto& filter : filters_)
                filter->adjustToEnsemble(ensembleInport_.getData());

            hash_ = hash;
        }
    }
}

void EnsembleFilter::applyFilter() {

    ensembleOutport_.clear();

    if (ensembleInport_.hasData()) {
        std::unique_ptr<EnsembleDataset> ensemble(new EnsembleDataset(*ensembleInport_.getData()));
        for (auto& filter : filters_) {
            ensemble.reset(filter->applyFilter(*ensemble));
        }

        ensembleOutport_.setData(ensemble.release(), true);
    }
}

void EnsembleFilter::serialize(Serializer& s) const {
    Processor::serialize(s);
    s.serialize("hash", hash_);
}

void EnsembleFilter::deserialize(Deserializer& s) {
    Processor::deserialize(s);
    s.optionalDeserialize("hash", hash_, std::string(""));
}

} // namespace
