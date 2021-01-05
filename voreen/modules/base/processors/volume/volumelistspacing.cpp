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

#include "volumelistspacing.h"

#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"

namespace voreen {

const std::string VolumeListSpacing::loggerCat_("voreen.base.VolumeListSpacing");

VolumeListSpacing::VolumeListSpacing()
    : VolumeProcessor()
    , inport_(Port::INPORT, "volumelist.input", "VolumeList Input")
    , outport_(Port::OUTPORT, "volumelist.output", "VolumeList Output", false)
    , enableProcessing_("enabled", "Enable", false)
    , mode_("mode", "Mode")
    , uniformSpacing_("uniformSpacing", "Uniform Spacing", false)
    , spacingX_("spacingX", "Spacing X (mm)", 1.f, 0.000001f, 1000.f)
    , spacingY_("spacingY", "Spacing Y (mm)", 1.f, 0.000001f, 1000.f)
    , spacingZ_("spacingZ", "Spacing Z (mm)", 1.f, 0.000001f, 1000.f)
    , reset_("reset", "Reset Spacing")
    , currentlySelected_("currentlySelected", "Currently Displayed Volume", 0, 0, std::numeric_limits<int>::max(), Processor::VALID, NumericProperty<int>::DYNAMIC)
    , spacingDisplay_("spacingDisplay", "Resulting Spacing (mm)", tgt::vec3(1.0f), tgt::vec3(0.0f), tgt::vec3(1000.f), Processor::VALID)
{
    addPort(inport_);
    inport_.Observable<PortObserver>::addObserver(this);
    addPort(outport_);

    enableProcessing_.onChange(MemberFunctionCallback<VolumeListSpacing>(this, &VolumeListSpacing::adjustPropertyVisibility));
    addProperty(enableProcessing_);

    mode_.addOption("replace", "Replace");
    mode_.addOption("scale", "Scale");
    addProperty(mode_);

    addProperty(uniformSpacing_);

    spacingX_.setNumDecimals(6);
    spacingX_.setStepping(0.001f);
    spacingY_.setNumDecimals(6);
    spacingY_.setStepping(0.001f);
    spacingZ_.setNumDecimals(6);
    spacingZ_.setStepping(0.001f);
    addProperty(spacingX_);
    addProperty(spacingY_);
    addProperty(spacingZ_);
    reset_.onChange(MemberFunctionCallback<VolumeListSpacing>(this, &VolumeListSpacing::resetSpacing));
    addProperty(reset_);

    currentlySelected_.setGroupID("currentDisplay");
    spacingDisplay_.setGroupID("currentDisplay");
    addProperty(currentlySelected_);
    spacingDisplay_.setReadOnlyFlag(true);
    spacingDisplay_.setNumDecimals(6);
    addProperty(spacingDisplay_);
    setPropertyGroupGuiName("currentDisplay", "Currently Selected Volume Spacing");

    spacingX_.onChange(
        MemberFunctionOneParameterCallback<VolumeListSpacing, int>(this, &VolumeListSpacing::spacingChanged, 0));
    spacingY_.onChange(
        MemberFunctionOneParameterCallback<VolumeListSpacing, int>(this, &VolumeListSpacing::spacingChanged, 1));
    spacingZ_.onChange(
        MemberFunctionOneParameterCallback<VolumeListSpacing, int>(this, &VolumeListSpacing::spacingChanged, 2));

    uniformSpacing_.onChange(MemberFunctionCallback<VolumeListSpacing>(this, &VolumeListSpacing::uniformScalingChanged));

    currentlySelected_.onChange(MemberFunctionCallback<VolumeListSpacing>(this, &VolumeListSpacing::updateCurrentlySelected));
}

VolumeListSpacing::~VolumeListSpacing() {
    // We need to remove the observer, since outport is destructed before inport (stack hierarchy)
    // and the destruction of inport will clear the (already destructed) outport otherwise.
    inport_.Observable<PortObserver>::removeObserver(this);
}

Processor* VolumeListSpacing::create() const {
    return new VolumeListSpacing();
}

void VolumeListSpacing::initialize() {
    VolumeProcessor::initialize();
    adjustPropertyVisibility();
}

void VolumeListSpacing::deinitialize() {
    clearOutput();
    VolumeProcessor::deinitialize();
}

void VolumeListSpacing::process() {
    // clear the old data
    clearOutput();

    // get input data
    const VolumeList* inputList = inport_.getData();

    // disabled or empty input list: just propagate input list
    if (!enableProcessing_.get() || !inputList || inputList->empty()) {
        outport_.setData(inputList, false);
        updateCurrentlySelected();  // update the displayed spacing in case the input has changed and the processor is disabled
        return;
    }

    // process the list
    VolumeList* outputList = new VolumeList();
    tgt::vec3 spacing(spacingX_.get(), spacingY_.get(), spacingZ_.get());
    bool scale = mode_.isSelected("scale");

    for (size_t i = 0; i < inputList->size(); ++i) {
        const VolumeBase* inputVolume = inputList->at(i);
        tgt::vec3 tmpSpacing = spacing;
        if (scale)
            tmpSpacing *= inputVolume->getSpacing();
        VolumeBase* outputVolume =
            new VolumeDecoratorReplaceSpacing(inputVolume, tmpSpacing);

        outputList->add(outputVolume);

        decorators_.push_back(std::unique_ptr<VolumeBase>(outputVolume));
    }
    outport_.setData(outputList, true);

    // update the spacing display
    updateCurrentlySelected();
}

void VolumeListSpacing::clearOutput() {
    outport_.clear();
    decorators_.clear();
}

void VolumeListSpacing::spacingChanged(int dim) {

    if (!uniformSpacing_.get())
        return;

    if (dim == 0) {
        float xScale = spacingX_.get();
        spacingY_.set(xScale);
        spacingZ_.set(xScale);
    }
    else if (dim == 1) {
        float yScale = spacingY_.get();
        spacingX_.set(yScale);
        spacingZ_.set(yScale);
    }
    else if (dim == 2) {
        float zScale = spacingZ_.get();
        spacingX_.set(zScale);
        spacingY_.set(zScale);
    }
}

void VolumeListSpacing::uniformScalingChanged() {
    if (uniformSpacing_.get())
        spacingChanged(0);
}

void VolumeListSpacing::adjustPropertiesToInput() {
    const VolumeList* inputList = inport_.getData();
    int max = ((inputList != 0) ? static_cast<int>(inputList->size() - 1) : 0);
    currentlySelected_.setMaxValue(std::max(0,max)); // the spacing is automatically adjusted afterwards
    updateCurrentlySelected();
}

void VolumeListSpacing::updateCurrentlySelected() {

    if (!inport_.getData() || inport_.getData()->empty())
        return;

    const VolumeBase* v = inport_.getData()->at(currentlySelected_.get());

    tgt::vec3 curSpacing = tgt::vec3(spacingX_.get(), spacingY_.get(), spacingZ_.get());
    if (!enableProcessing_.get()) {
        spacingDisplay_.set(v->getSpacing());
    }
    else {
        if (mode_.isSelected("replace")) {
            spacingDisplay_.set(curSpacing);
        }
        else if (mode_.isSelected("scale")) {
            spacingDisplay_.set(curSpacing * v->getSpacing());
        }
    }
}

void VolumeListSpacing::resetSpacing() {
    mode_.selectByKey("scale");
    spacingX_.set(1.f);
    spacingY_.set(1.f);
    spacingZ_.set(1.f);
}

void VolumeListSpacing::adjustPropertyVisibility() {
    bool enabled = enableProcessing_.get();
    mode_.setReadOnlyFlag(!enabled);
    uniformSpacing_.setReadOnlyFlag(!enabled);
    spacingX_.setReadOnlyFlag(!enabled);
    spacingY_.setReadOnlyFlag(!enabled);
    spacingZ_.setReadOnlyFlag(!enabled);
    reset_.setReadOnlyFlag(!enabled);
}

void VolumeListSpacing::dataWillChange(const Port* source) {
    clearOutput();
}

}   // namespace
