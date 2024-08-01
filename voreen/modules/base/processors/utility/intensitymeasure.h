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

#ifndef VRN_INTENSITYMEASURE_H
#define VRN_INTENSITYMEASURE_H

#include "voreen/core/processors/imageprocessor.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/utils/stringutils.h"
#include "tgt/font.h"
#include "tgt/glmath.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/textport.h"

namespace voreen {

class VRN_CORE_API IntensityMeasure : public ImageProcessor {
public:
    IntensityMeasure();
    ~IntensityMeasure();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "IntensityMeasure"; }
    virtual std::string getCategory() const  { return "Utility";          }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE;  }
    virtual bool isUtility() const           { return true; }

    virtual bool isReady() const;

    void measure(tgt::MouseEvent* e);
    virtual void onEvent(tgt::Event* e);

protected:
    virtual void setDescriptions() {
        setDescription("Allows to interactively measure intensities.");
    }

    void process();

private:
    tgt::ivec2 cropToViewport(tgt::ivec2 mousePos, tgt::ivec2 viewport);
    bool measureIntensity(tgt::vec3 samplePos, const VolumeBase* vh, float& intensity, std::string& unit);
    tgt::ivec3 worldToVoxel(tgt::vec3 samplePos, const VolumeBase* vh);

    RenderPort fhpInport_;
    TextPort outport_;
    VolumePort volumePort0_;
    VolumePort volumePort1_;
    VolumePort volumePort2_;
    VolumePort volumePort3_;

    EventProperty<IntensityMeasure>* mouseEventProp_;
    //BoolProperty requirePress_;
    BoolProperty useVolume0_;
    BoolProperty useVolume1_;
    BoolProperty useVolume2_;
    BoolProperty useVolume3_;
    IntProperty numDigits0_;
    IntProperty numDigits1_;
    IntProperty numDigits2_;
    IntProperty numDigits3_;

    tgt::ivec2 mousePos2D_;
    tgt::Color mousePos3D_;
    bool showText_;

    tgt::Font* font_;
};

} // namespace

#endif // VRN_INTENSITYMEASURE_H
