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

#include "volumelistoffset.h"

#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"

namespace voreen {

const std::string VolumeListOffset::loggerCat_("voreen.base.VolumeListOffset");

VolumeListOffset::VolumeListOffset()
    : VolumeProcessor()
    , inport_(Port::INPORT, "volumelist.input", "VolumeList Input")
    , outport_(Port::OUTPORT, "volumelist.output", "VolumeList Output", false)
    , enableProcessing_("enabled", "Enable", false)
    , offset_("offset", "Set Offset", tgt::vec3(0.0f), tgt::vec3(-10000000.0), tgt::vec3(10000000.0))
    , reset_("reset", "Reset Offset")
    , currentlySelected_("currentlySelected", "Currently Displayed Volume", 0, 0, std::numeric_limits<int>::max(), Processor::VALID, NumericProperty<int>::DYNAMIC)
    , offsetDisplay_("offsetDisplay", "Resulting Offset (mm)", tgt::vec3(0.0f), tgt::vec3(-10000000.0), tgt::vec3(10000000.f), Processor::VALID)
{
    addPort(inport_);
    inport_.Observable<PortObserver>::addObserver(this);
    addPort(outport_);

    enableProcessing_.onChange(MemberFunctionCallback<VolumeListOffset>(this, &VolumeListOffset::adjustPropertyVisibility));
    addProperty(enableProcessing_);

    offset_.setNumDecimals(6);
    offset_.setStepping(tgt::vec3(0.001f));
    addProperty(offset_);

    reset_.onChange(MemberFunctionCallback<VolumeListOffset>(this, &VolumeListOffset::resetOffset));
    addProperty(reset_);

    currentlySelected_.setGroupID("currentDisplay");
    offsetDisplay_.setGroupID("currentDisplay");
    addProperty(currentlySelected_);
    offsetDisplay_.setReadOnlyFlag(true);
    offsetDisplay_.setNumDecimals(6);
    addProperty(offsetDisplay_);
    setPropertyGroupGuiName("currentDisplay", "Currently Selected Volume Offset");

    currentlySelected_.onChange(MemberFunctionCallback<VolumeListOffset>(this, &VolumeListOffset::updateCurrentlySelected));
}

VolumeListOffset::~VolumeListOffset() {
    // We need to remove the observer, since outport is destructed before inport (stack hierarchy)
    // and the destruction of inport will clear the (already destructed) outport otherwise.
    inport_.Observable<PortObserver>::removeObserver(this);
}

Processor* VolumeListOffset::create() const {
    return new VolumeListOffset();
}

void VolumeListOffset::initialize() {
    VolumeProcessor::initialize();
    adjustPropertyVisibility();
}

void VolumeListOffset::deinitialize() {
    clearOutput();
    VolumeProcessor::deinitialize();
}

void VolumeListOffset::process() {
    // clear the old data
    clearOutput();

    // get input data
    const VolumeList* inputList = inport_.getData();

    // disabled or empty input list: just propagate input list
    if (!enableProcessing_.get() || !inputList || inputList->empty()) {
        outport_.setData(inputList, false);
        updateCurrentlySelected();  // update the displayed offset in case the input has changed and the processor is disabled
        return;
    }

    // process the list
    VolumeList* outputList = new VolumeList();

    for (size_t i = 0; i < inputList->size(); ++i) {
        const VolumeBase* inputVolume = inputList->at(i);

        VolumeBase* outputVolume =
            new VolumeDecoratorReplaceOffset(inputVolume, offset_.get());

        outputList->add(outputVolume);

        decorators_.push_back(std::unique_ptr<VolumeBase>(outputVolume));
    }
    outport_.setData(outputList, true);

    // update the offset display
    updateCurrentlySelected();
}

void VolumeListOffset::clearOutput() {
    outport_.clear();
    decorators_.clear();
}

void VolumeListOffset::adjustPropertiesToInput() {
    const VolumeList* inputList = inport_.getData();
    int max = ((inputList != 0) ? static_cast<int>(inputList->size() - 1) : 0);
    currentlySelected_.setMaxValue(std::max(0,max)); // the offset is automatically adjusted afterwards
    updateCurrentlySelected();
}

void VolumeListOffset::updateCurrentlySelected() {

    if (!inport_.getData() || inport_.getData()->empty())
        return;

    const VolumeBase* v = inport_.getData()->at(currentlySelected_.get());

    if (!enableProcessing_.get()) {
        offsetDisplay_.set(v->getOffset());
    }
    else {
        offsetDisplay_.set(offset_.get());
    }
}

void VolumeListOffset::resetOffset() {
    offset_.set(tgt::vec3::zero);
}

void VolumeListOffset::adjustPropertyVisibility() {
    bool enabled = enableProcessing_.get();
    offset_.setReadOnlyFlag(!enabled);
    reset_.setReadOnlyFlag(!enabled);
}

void VolumeListOffset::dataWillChange(const Port* source) {
    clearOutput();
}

}   // namespace
