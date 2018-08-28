/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "intensitymeasure.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"

#include "tgt/textureunit.h"

using tgt::MouseEvent;
using tgt::TextureUnit;
using tgt::ivec3;
using tgt::svec3;
using tgt::vec3;
using std::string;

namespace voreen {

IntensityMeasure::IntensityMeasure()
    : ImageProcessor("image/distance")
    , fhpInport_(Port::INPORT, "fhp", "First-hit-points Input")
    , outport_(Port::OUTPORT, "output", "Text Output")
    , volumePort0_(Port::INPORT, "volume.input0", "Volume0 Input")
    , volumePort1_(Port::INPORT, "volume.input1", "Volume1 Input")
    , volumePort2_(Port::INPORT, "volume.input2", "Volume2 Input")
    , volumePort3_(Port::INPORT, "volume.input3", "Volume3 Input")
    //, requirePress_("requirePress", "Require Mouse Press?", false)
    , useVolume0_("useVolume0", "Use volume 0", true)
    , useVolume1_("useVolume1", "Use volume 1", false)
    , useVolume2_("useVolume2", "Use volume 2", false)
    , useVolume3_("useVolume3", "Use volume 3", false)
    , numDigits0_("numDigits0", "Number of digits 0", 2, 1, 12)
    , numDigits1_("numDigits1", "Number of digits 1", 2, 1, 12)
    , numDigits2_("numDigits2", "Number of digits 2", 2, 1, 12)
    , numDigits3_("numDigits3", "Number of digits 3", 2, 1, 12)
{
    addPort(fhpInport_);
    addPort(volumePort0_);
    addPort(volumePort1_);
    addPort(volumePort2_);
    addPort(volumePort3_);
    addPort(outport_);

    //addProperty(requirePress_);
    addProperty(useVolume0_);
    addProperty(numDigits0_);
    addProperty(useVolume1_);
    addProperty(numDigits1_);
    addProperty(useVolume2_);
    addProperty(numDigits2_);
    addProperty(useVolume3_);
    addProperty(numDigits3_);

    mouseEventProp_ = new EventProperty<IntensityMeasure>(
        "mouseEvent.measure", "Intensity measure", this,
        &IntensityMeasure::measure,
        tgt::MouseEvent::MOUSE_BUTTON_ALL,
        tgt::MouseEvent::MOTION | tgt::MouseEvent::PRESSED | tgt::MouseEvent::RELEASED,
        tgt::Event::ALT, true);
    //addEventProperty(mouseEventProp_);
    showText_ = false;
    mousePos2D_ = tgt::ivec2(0, 0);
    mousePos3D_ = tgt::Color(0.0f);

    font_ = new tgt::Font(VoreenApplication::app()->getFontPath("VeraMono.ttf"), 11);
}

IntensityMeasure::~IntensityMeasure() {
    delete font_;
    delete mouseEventProp_;
}

Processor* IntensityMeasure::create() const {
    return new IntensityMeasure();
}

bool IntensityMeasure::isReady() const {
    if (!isInitialized() || !fhpInport_.isReady() || !outport_.isReady())
        return false;

    if (!volumePort0_.isReady())
        return false;

    return true;
}

tgt::ivec2 IntensityMeasure::cropToViewport(tgt::ivec2 mousePos, tgt::ivec2 viewport) {
    tgt::ivec2 result = mousePos;
    tgt::ivec2 size = viewport;

    if (result.x < 0)
        result.x = 0;
    else if (result.x > size.x-1)
        result.x = size.x-1;

    if (result.y < 0)
        result.y = 0;
    else if (result.y > size.y-1)
        result.y = size.y-1;

    return result;
}

void IntensityMeasure::measure(tgt::MouseEvent* e) {
    if (((e->action() & tgt::MouseEvent::PRESSED) && !showText_) ||
        ((e->action() & tgt::MouseEvent::MOTION) && showText_) //||
        //!requirePress_.get()
        ) {
        mousePos2D_ = cropToViewport(tgt::ivec2(e->coord().x, e->viewport().y-e->coord().y), e->viewport());

        tgt::Color curColor = fhpInport_.getRenderTarget()->getColorAtPos(mousePos2D_);

        if (curColor.a > 0.5f) {
            if (!showText_)
                showText_ = true;
            mousePos3D_ = curColor;
            invalidate();
            e->accept();
        }
    }

    if (e->action() & tgt::MouseEvent::RELEASED) {
        if (showText_) {
            showText_ = false;
            mousePos2D_ = tgt::ivec2(0, 0);
            mousePos3D_ = tgt::Color(0.0f);
            invalidate();
            e->accept();
        }
        outport_.setData("");
    }
}

void IntensityMeasure::onEvent(tgt::Event* e) {
    MouseEvent* me = dynamic_cast<MouseEvent*>(e);
    if(me) {
        if(me->action() & MouseEvent::EXIT) {
            showText_ = false;
            outport_.setData("");
        }
        //else if(me->action() & MouseEvent::ENTER)
            //showText_ = true;
        else if(me->action() & MouseEvent::MOTION) {
            mousePos2D_ = cropToViewport(tgt::ivec2(me->coord().x, me->viewport().y-me->coord().y), me->viewport());

            tgt::Color curColor = fhpInport_.getRenderTarget()->getColorAtPos(mousePos2D_);

            if (curColor.a > 0.5f) {
                showText_ = true;
                mousePos3D_ = curColor;
                invalidate();
                //e->accept();
            }
            else {
                showText_ = false;
                outport_.setData("");
            }
        }

    }
    RenderProcessor::onEvent(e);
}

bool IntensityMeasure::measureIntensity(tgt::vec3 samplePos, const VolumeBase* vh, float& intensity, std::string& unit) {
    tgt::vec3 voxelSamplePos = vh->getWorldToVoxelMatrix() * samplePos;

    if ((min(voxelSamplePos) >= 0.0) && tgt::hand(lessThan(voxelSamplePos, vec3(vh->getDimensions())))) {
        const VolumeRAM* volume = vh->getRepresentation<VolumeRAM>();

        intensity = volume->getVoxelNormalizedLinear(voxelSamplePos);

        RealWorldMapping rwm = vh->getRealWorldMapping();
        intensity = rwm.normalizedToRealWorld(intensity);
        unit = rwm.getUnit();
        return true;
    }
    else
        return false;
}

std::string generateLabel(float intensity, int numDigits, const std::string& unit) {
    std::string label = ftos(intensity);
    if(label.find(".") != std::string::npos)
        label = label.substr(0, label.find(".") + numDigits + 1);
    label += " " + unit;
    label += " ";
    return label;
}

void IntensityMeasure::process() {
    if (showText_) {
        // construct the string for display
        std::string label;
        tgt::vec3 samplePos = static_cast<tgt::vec3>(mousePos3D_.xyz());

        if (useVolume0_.get() && volumePort0_.isReady()) {
            float intensity = 0.0f;
            string unit = "";
            if(measureIntensity(samplePos, volumePort0_.getData(), intensity, unit)) {
                std::stringstream strstr;
                //strstr << worldToVoxel(samplePos, volumePort0_.getData()) << " ";
                label += strstr.str();
                label += generateLabel(intensity, numDigits0_.get(), unit);
                label += "\n";
            }
        }
        if (useVolume1_.get() && volumePort1_.isReady()) {
            float intensity = 0.0f;
            string unit = "";
            if(measureIntensity(samplePos, volumePort1_.getData(), intensity, unit)) {
                std::stringstream strstr;
                //strstr << worldToVoxel(samplePos, volumePort1_.getData()) << " ";
                label += strstr.str();
                label += generateLabel(intensity, numDigits1_.get(), unit);
                label += "\n";
            }
        }
        if (useVolume2_.get() && volumePort2_.isReady()) {
            float intensity = 0.0f;
            string unit = "";
            if(measureIntensity(samplePos, volumePort2_.getData(), intensity, unit)) {
                std::stringstream strstr;
                //strstr << worldToVoxel(samplePos, volumePort2_.getData()) << " ";
                label += strstr.str();
                label += generateLabel(intensity, numDigits2_.get(), unit);
                label += "\n";
            }
        }
        if (useVolume3_.get() && volumePort3_.isReady()) {
            float intensity = 0.0f;
            string unit = "";
            if(measureIntensity(samplePos, volumePort3_.getData(), intensity, unit)) {
                std::stringstream strstr;
                //strstr << worldToVoxel(samplePos, volumePort3_.getData()) << " ";
                label += strstr.str();
                label += generateLabel(intensity, numDigits3_.get(), unit);
                label += "\n";
            }
        }

        outport_.setData(label);
    }
    else {
        outport_.setData("");
    }

    LGL_ERROR;
}

} // namespace voreen
