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

#include "volumelistrealworldmapping.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"

namespace voreen {

const std::string VolumeListRealWorldMapping::loggerCat_("voreen.base.VolumeListRealWorldMapping");

VolumeListRealWorldMapping::VolumeListRealWorldMapping()
    : VolumeProcessor()
    , inport_(Port::INPORT, "volumelist.input", "Volume Input")
    , outport_(Port::OUTPORT, "volumelist.output", "Volume Output", false)
    , enableProcessing_("enabled", "Enable", false)
    , replace_("replace", "Replace", true)
    , mode_("mode_", "Mode")
    , realWorldRange_("realWorldRange", "Real world range", tgt::vec2(0.0f, 1.0f), 0, 100000)
    , scale_("scale", "Scale", 1.0f, 0.0f, 10000.0f)
    , offset_("offset", "Offset", 0.0f, -100000.0f, 100000.0f)
{
    addPort(inport_);
    inport_.Observable<PortObserver>::addObserver(this);
    addPort(outport_);
    addProperty(enableProcessing_);

    addProperty(replace_);

    addProperty(mode_);
        mode_.addOption("normalize", "Normalized", Mode::NORMALIZED);
        mode_.addOption("denormalize", "Denormalized", Mode::DENORMALIZED);
        mode_.addOption("scaleOffset", "Scale and Offset", Mode::SCALE_OFFSET);
        mode_.addOption("interval", "Interval", Mode::INTERVAL);
    addProperty(scale_);
    scale_.setNumDecimals(4);
    addProperty(offset_);
    offset_.setNumDecimals(4);

    addProperty(realWorldRange_);

    updatePropertyVisibility();
    ON_CHANGE(mode_, VolumeListRealWorldMapping, updatePropertyVisibility);

    ON_CHANGE(scale_, VolumeListRealWorldMapping, updatePropertyValues);
    ON_CHANGE(offset_, VolumeListRealWorldMapping, updatePropertyValues);
    ON_CHANGE(realWorldRange_, VolumeListRealWorldMapping, updatePropertyValues);

}

VolumeListRealWorldMapping::~VolumeListRealWorldMapping() {
    // We need to remove the observer, since outport is destructed before inport (stack hierarchy)
    // and the destruction of inport will clear the (already destructed) outport otherwise.
    inport_.Observable<PortObserver>::removeObserver(this);
}

Processor* VolumeListRealWorldMapping::create() const {
    return new VolumeListRealWorldMapping();
}

void VolumeListRealWorldMapping::initialize() {
    VolumeProcessor::initialize();
}

void VolumeListRealWorldMapping::deinitialize() {
    clearOutput();
    VolumeProcessor::deinitialize();
}

void VolumeListRealWorldMapping::process() {
    clearOutput();

    const VolumeList* inputList = inport_.getData();

    if (!enableProcessing_.get() || !inputList || inputList->empty()) {
        outport_.setData(inputList, false);
        return;
    }

    VolumeList* outputList = new VolumeList();

    for (size_t i = 0; i < inputList->size(); ++i) {
        const VolumeBase* inputVolume = inputList->at(i);

        RealWorldMapping rwm = getRealWorldMapping(inputVolume);
        if(!replace_.get()) {
            rwm = RealWorldMapping::combine(inputVolume->getRealWorldMapping(), rwm);
        }

        VolumeBase* outputVolume = new VolumeDecoratorReplaceRealWorldMapping(inputVolume, rwm);
        outputList->add(outputVolume);

        decorators_.emplace_back(std::unique_ptr<VolumeBase>(outputVolume));
    }
    outport_.setData(outputList, true);
}


void VolumeListRealWorldMapping::reset() {
    realWorldRange_.set(tgt::vec2(0.0f, 1.0f));
}

void VolumeListRealWorldMapping::clearOutput() {
    outport_.clear();
    decorators_.clear();
}

void VolumeListRealWorldMapping::updatePropertyVisibility(){
    Mode mode = mode_.getValue();
    bool scaleOffsetVisible = mode == SCALE_OFFSET;
    bool intervalVisible = mode == INTERVAL;

    scale_.setVisibleFlag(scaleOffsetVisible);
    offset_.setVisibleFlag(scaleOffsetVisible);

    realWorldRange_.setVisibleFlag(intervalVisible);
}

void VolumeListRealWorldMapping::updatePropertyValues(){
    switch (mode_.getValue())
    {
    case INTERVAL:
    {
        tgt::vec2 range = realWorldRange_.get();
        scale_.set(range.y-range.x);
        offset_.set(range.x);
    }
    break;
    case SCALE_OFFSET:
    {
        float scale = scale_.get();
        float offset = offset_.get();
        realWorldRange_.set(tgt::vec2(offset, offset+scale));
    }
    break;
    default:
        break;
    }
}

RealWorldMapping VolumeListRealWorldMapping::getRealWorldMapping(const VolumeBase* volume){
    switch (mode_.getValue())
    {
    case INTERVAL:
    {
        tgt::vec2 range = realWorldRange_.get();
        return RealWorldMapping(range.y-range.x, range.x, "");
    }
    case SCALE_OFFSET:
    {
        float scale = scale_.get();
        float offset = offset_.get();
        return RealWorldMapping(scale, offset, "");
    }
    case NORMALIZED:
        return RealWorldMapping(1.0f, 0.0f, "");
    case DENORMALIZED:
    {
        return RealWorldMapping::createDenormalizingMapping(volume->getBaseType());
    }
    default:
        break;
    }
    return RealWorldMapping(1, 0.0f, "");
}

void VolumeListRealWorldMapping::dataWillChange(const Port* source) {
    clearOutput();
}

}   // namespace
