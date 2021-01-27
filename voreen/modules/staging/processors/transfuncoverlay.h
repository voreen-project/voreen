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

#ifndef VRN_TRANSFUNCOVERLAY_H
#define VRN_TRANSFUNCOVERLAY_H

#include "voreen/core/processors/imageprocessor.h"

#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/fontproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "tgt/texturemanager.h"
#include "tgt/gpucapabilities.h"

namespace voreen {

/**
 * Overlays a transfer function on top of the rendering.
 */
class TransFuncOverlay : public ImageProcessor {
public:
    TransFuncOverlay();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "TransFuncOverlay"; }
    virtual std::string getCategory() const   { return "Image Processing"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_TESTING;  }

    bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Provides an overlay that renders a transfer function.");
    }

    virtual void beforeProcess();
    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);


private:
    void onChangeUsePixelCoordinates();
    void onChangeDeriveRangeFromTransFunc();

    /*
     * Determine overlay dimensions and bottom-left in float pixel coords.
     */
    void getCurrentOverlayPixelCoordinates(tgt::ivec2& bl, tgt::ivec2& dim) const;

    RenderPort imageInport_;
    RenderPort privatePort_;
    RenderPort outport_;

    FontProperty fontProp_;
    TransFunc1DKeysProperty transferFunc_;
    BoolProperty renderPreIntegrationTable_;
    BoolProperty renderOverlay_;
    OptionProperty<bool> usePixelCoordinates_;
    IntVec2Property overlayBottomLeft_;             ///< pixel coordinates
    IntVec2Property overlayDimensions_;             ///< pixel coordinates
    FloatVec2Property overlayBottomLeftRelative_;   ///< normalized coordinates
    FloatVec2Property overlayDimensionsRelative_;   ///< normalized coordinates
    FloatProperty overlayOpacity_;
    ColorProperty fontColor_;
    BoolProperty renderRangeProp_;                  ///< render the tf range only if true
    FloatVec2Property valueRangeProp_;              ///< value range for transfunc domain
    BoolProperty deriveRangeFromTransFunc_;
    IntOptionProperty rangePrecisionProp_;          ///< number decimals in the overlay
    StringProperty tfUnit_;
    BoolProperty renderBorder_;
    FloatProperty borderWidth_;
    ColorProperty borderColor_;

    tgt::Shader* copyShader_;
};

} // namespace

#endif
