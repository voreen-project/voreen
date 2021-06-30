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

#include "ensemblefilter.h"

#include "voreen/core/datastructures/callback/lambdacallback.h"

#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/string/stringlistproperty.h"

#include "../utils/ensemblehash.h"

namespace voreen {

class Filter {
public:
    virtual ~Filter() {}
    virtual std::vector<Property*> getProperties() = 0;
    virtual EnsembleDataset* applyFilter(const EnsembleDataset& ensemble) = 0;
    virtual void adjustToEnsemble(const EnsembleDataset* ensemble) = 0;
    virtual bool isActive() const = 0;
};

//----------------------------------------
// Filter : Member
//----------------------------------------

class FilterMember : public Filter {
public:
    FilterMember()
        : members_("members", "Selected Members")
    {
        members_.setDescription("Selects multiple members from the ensemble data.");
    }

    std::vector<Property*> getProperties() {
        return std::vector<Property*>(1, &members_);
    }

    EnsembleDataset* applyFilter(const EnsembleDataset& ensemble) {
        
        EnsembleDataset* dataset = new EnsembleDataset();

        for (int row : members_.getSelectedRowIndices()) {
            EnsembleMember member = ensemble.getMembers()[row];
            dataset->addMember(member);
        }

        return dataset;
    }

    void adjustToEnsemble(const EnsembleDataset* ensemble) {
        
        // Reset range.
        members_.reset();

        if (ensemble) {
            // Adjust range to data.
            std::vector<int> selectedMemberIndices;
            for (const EnsembleMember& member : ensemble->getMembers()) {
                members_.addRow(member.getName(), member.getColor());
                selectedMemberIndices.push_back(static_cast<int>(selectedMemberIndices.size()));
            }
            members_.setSelectedRowIndices(selectedMemberIndices);
        }
    }

    bool isActive() const {
        return static_cast<size_t>(members_.getNumRows()) != members_.getSelectedRowIndices().size();
    }

private:
    StringListProperty members_;
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

    std::vector<Property*> getProperties() {
        return std::vector<Property*>(1, &timeSteps_);
    }

    EnsembleDataset* applyFilter(const EnsembleDataset& ensemble) {
        
        EnsembleDataset* dataset = new EnsembleDataset();

        for (const EnsembleMember& member : ensemble.getMembers()) {
            std::vector<TimeStep> timeSteps;
            int max = std::min(static_cast<int>(member.getTimeSteps().size()) - 1, timeSteps_.get().y);
            for (int i = timeSteps_.get().x; i <= max; i++) {
                timeSteps.push_back(member.getTimeSteps()[i]);
            }

            if(!timeSteps.empty()) {
                dataset->addMember(EnsembleMember{member.getName(), member.getColor(), timeSteps});
            }
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
            timeSteps_.setMaxValue(static_cast<int>(ensemble->getMaxNumTimeSteps()) - 1);
            timeSteps_.set(tgt::ivec2(0, static_cast<int>(ensemble->getMaxNumTimeSteps()) - 1));
        }
    }

    bool isActive() const {
        return timeSteps_.get() != tgt::ivec2(timeSteps_.getMinValue(), timeSteps_.getMaxValue());
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
        , keepIfOnlyTimeStep_("keepIfOnlyTimeStep", "Keep if only single Time Step", true)
    {
        ON_CHANGE_LAMBDA(enableRemoveFirstTimeStep_, [this] {
            keepIfOnlyTimeStep_.setVisibleFlag(enableRemoveFirstTimeStep_.get());
        });
        enableRemoveFirstTimeStep_.setDescription("Removes the first time step of each member.");
        keepIfOnlyTimeStep_.setDescription("Keep Time Step, if member only has a single one.");
        enableRemoveFirstTimeStep_.invalidate();
    }

    std::vector<Property*> getProperties() {
        std::vector<Property*> properties;
        properties.push_back(&enableRemoveFirstTimeStep_);
        properties.push_back(&keepIfOnlyTimeStep_);
        return properties;
    }

    EnsembleDataset* applyFilter(const EnsembleDataset& ensemble) {

        EnsembleDataset* dataset = new EnsembleDataset();

        for (const EnsembleMember& member : ensemble.getMembers()) {
            std::vector<TimeStep> timeSteps;

            // If the member only contains a single time step, we keep it
            if(member.getTimeSteps().size() == 1 && keepIfOnlyTimeStep_.get()) {
                timeSteps.push_back(member.getTimeSteps().front());
            }

            for (size_t i = 1; i < member.getTimeSteps().size(); i++) {
                timeSteps.push_back(member.getTimeSteps()[i]);
            }

            if(!timeSteps.empty()) {
                dataset->addMember(EnsembleMember{member.getName(), member.getColor(), timeSteps});
            }
        }

        return dataset;
    }

    void adjustToEnsemble(const EnsembleDataset* ensemble) {
    }

    bool isActive() const {
        return enableRemoveFirstTimeStep_.get();
    }

private:

    BoolProperty enableRemoveFirstTimeStep_;
    BoolProperty keepIfOnlyTimeStep_;
};

//----------------------------------------
// Filter : Select last Time Step
//----------------------------------------
class FilterSelectLastTimeStep : public Filter {
public:
    FilterSelectLastTimeStep()
        : enableSelectLastTimeStep_("enableSelectLastTimeStep", "Select last Time Step", false)
    {
        enableSelectLastTimeStep_.setDescription("Selects only the last time step of each member.");
    }

    std::vector<Property*> getProperties() {
        return std::vector<Property*>(1, &enableSelectLastTimeStep_);
    }

    EnsembleDataset* applyFilter(const EnsembleDataset& ensemble) {

        EnsembleDataset* dataset = new EnsembleDataset();

        for (const EnsembleMember& member : ensemble.getMembers()) {
            if (member.getTimeSteps().empty())
                continue;

            std::vector<TimeStep> timeSteps;
            timeSteps.push_back(member.getTimeSteps().back());

            dataset->addMember(EnsembleMember{member.getName(), member.getColor(), timeSteps});
        }

        return dataset;
    }

    void adjustToEnsemble(const EnsembleDataset* ensemble) {
    }

    bool isActive() const {
        return enableSelectLastTimeStep_.get();
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

    std::vector<Property*> getProperties() {
        return std::vector<Property*>(1, &timeInterval_);
    }

    EnsembleDataset* applyFilter(const EnsembleDataset& ensemble) {

        EnsembleDataset* dataset = new EnsembleDataset();
        float epsilon = ensemble.getMinTimeStepDuration();

        for (const EnsembleMember& member : ensemble.getMembers()) {
            std::vector<TimeStep> timeSteps;
            for (const auto& timeStep : member.getTimeSteps()) {
                if (timeStep.getTime() > timeInterval_.get().y + epsilon)
                    break;
                if (timeStep.getTime() >= timeInterval_.get().x)
                    timeSteps.push_back(timeStep);
            }

            if(!timeSteps.empty()) {
                dataset->addMember(EnsembleMember{member.getName(), member.getColor(), timeSteps});
            }
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

    bool isActive() const {
        return timeInterval_.get() != tgt::vec2(timeInterval_.getMinValue(), timeInterval_.getMaxValue());
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
                                 "(*) Marks common fields across all members."
        );
    }

    std::vector<Property*> getProperties() {
        return std::vector<Property*>(1, &fields_);
    }

    EnsembleDataset* applyFilter(const EnsembleDataset& ensemble) {

        EnsembleDataset* dataset = new EnsembleDataset();

        for (const EnsembleMember& member : ensemble.getMembers()) {
            std::vector<TimeStep> timeSteps;
            for (const TimeStep& timeStep : member.getTimeSteps()) {

                std::vector<std::string> fieldNames;
                fieldNames.push_back(fields_.getValue());

                // Only add time step, if selected field is available.
                TimeStep filtered = timeStep.createSubset(fieldNames);
                if(!filtered.getFieldNames().empty()) {
                    timeSteps.push_back(filtered);
                }
            }
            if(!timeSteps.empty()) {
                dataset->addMember(EnsembleMember{member.getName(), member.getColor(), timeSteps});
            }
        }

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

    bool isActive() const {
        return fields_.getOptions().size() > 1;
    }

private:

    OptionProperty<std::string> fields_;
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

    addFilter(new FilterMember());
    //addFilter(new FilterTimeStep()); // replaced by FilterTimeInterval
    addFilter(new FilterTimeInterval());
    addFilter(new FilterRemoveFirstTimeStep());
    addFilter(new FilterSelectLastTimeStep());
    addFilter(new FilterField());
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

    if (inv == Processor::INVALID_RESULT)
        needsProcess_ = true;
}

void EnsembleFilter::addFilter(Filter* filter) {
    filters_.emplace_back(std::unique_ptr<Filter>(filter));
    for(Property* property : filter->getProperties()) {
        addProperty(property);
    }
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
            if(filter->isActive()) {
                ensemble.reset(filter->applyFilter(*ensemble));
            }
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
