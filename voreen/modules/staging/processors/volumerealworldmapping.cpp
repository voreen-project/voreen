/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "volumerealworldmapping.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"

namespace voreen {

const std::string VolumeRealWorldMapping::loggerCat_("voreen.base.VolumeRVMNormalization");

VolumeRealWorldMapping::VolumeRealWorldMapping()
    : VolumeProcessor()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false)
    , enableProcessing_("enabled", "Enable", false)
    , mode_("mode_", "Mode")
    , realWorldRange_("realWorldRange", "Real world range", tgt::vec2(0.0f, 1.0f), 0, 100000)
    , scale_("scale", "Scale", 1.0f, 0.0f, 100000.0f)
    , offset_("offset", "Offset", 0.0f, 0.0f, 100000.0f)
{
    addPort(inport_);
    addPort(outport_);
    addProperty(enableProcessing_);

    addProperty(mode_);
        mode_.addOption("normalize", "Normalized", Mode::NORMALIZED);
        mode_.addOption("denormalize", "Denormalized", Mode::DENORMALIZED);
        mode_.addOption("scaleOffset", "Scale and Offset", Mode::SCALE_OFFSET);
        mode_.addOption("interval", "Interval", Mode::INTERVAL);
    addProperty(scale_);
    addProperty(offset_);

    addProperty(realWorldRange_);

    updatePropertyVisibility();
    ON_CHANGE(mode_, VolumeRealWorldMapping, updatePropertyVisibility);

    ON_CHANGE(scale_, VolumeRealWorldMapping, updatePropertyValues);
    ON_CHANGE(offset_, VolumeRealWorldMapping, updatePropertyValues);
    ON_CHANGE(realWorldRange_, VolumeRealWorldMapping, updatePropertyValues);

}

Processor* VolumeRealWorldMapping::create() const {
    return new VolumeRealWorldMapping();
}

void VolumeRealWorldMapping::initialize() {
    VolumeProcessor::initialize();
}

void VolumeRealWorldMapping::process() {
    const VolumeBase* inputVolume = inport_.getData();

    if (!enableProcessing_.get()) {
        outport_.setData(inputVolume, false);
        return;
    }

    RealWorldMapping rvm = getRealWorldMapping();

    VolumeBase* outputVolume = new VolumeDecoratorReplaceRealWorldMapping(inputVolume, rvm);
    outport_.setData(outputVolume);
}


void VolumeRealWorldMapping::reset() {
    realWorldRange_.set(tgt::vec2(0.0f, 1.0f));
}

void VolumeRealWorldMapping::updatePropertyVisibility(){
    Mode mode = mode_.getValue();
    bool scaleOffsetVisible = mode == SCALE_OFFSET;
    bool intervalVisible = mode == INTERVAL;

    scale_.setVisibleFlag(scaleOffsetVisible);
    offset_.setVisibleFlag(scaleOffsetVisible);

    realWorldRange_.setVisibleFlag(intervalVisible);
}

void VolumeRealWorldMapping::updatePropertyValues(){
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

RealWorldMapping VolumeRealWorldMapping::getRealWorldMapping(){
    switch (mode_.getValue())
    {
    case INTERVAL:
        {
            tgt::vec2 range = realWorldRange_.get();
            return RealWorldMapping(range.y-range.x, range.x, "");
        }
        break;
    case SCALE_OFFSET:
        {
            float scale = scale_.get();
            float offset = offset_.get();
            return RealWorldMapping(scale, offset, "");
        }
        break;
    case NORMALIZED:
        return RealWorldMapping(1.0f, 0.0f, "");
        break;
    case DENORMALIZED:
        {
            const VolumeBase* inputVolume = inport_.getData();
            if (!inputVolume) break;
            RealWorldMapping rvm = inputVolume->getRealWorldMapping();
            return rvm.createDenormalizingMapping(inputVolume->getBaseType());
        }
    default:
        break;
    }
    return RealWorldMapping(1, 0.0f, "");
}

}   // namespace
