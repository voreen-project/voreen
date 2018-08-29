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

#ifndef VRN_COLORCODINGABILITY_H
#define VRN_COLORCODINGABILITY_H

#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"

namespace voreen {

class ColorCodingAbility {
public:
    enum ColorCodingMode {
        COLOR_CODING_MODE_MAGNITUDE,
        COLOR_CODING_MODE_DIRECTION,
        COLOR_CODING_MODE_MONOCHROME
    };

public:
    ColorCodingAbility();
    ~ColorCodingAbility();

    ColorCodingMode getColorCodingMode() const { return colorCodingModeProp_.getValue(); }
    OptionProperty<ColorCodingMode>& getColorCodingModeProp() const { return colorCodingModeProp_; }
    IntOptionProperty& getColorTableProp() const { return colorTableProp_; }
    ColorProperty& getColorProp() const { return colorProp_; }
    TransFunc1DKeysProperty& getTransFuncProp() const { return transFuncProp_; }
    std::string getShaderDefines() const;

    void onColorCodingModeChange();

private:
    mutable OptionProperty<ColorCodingMode> colorCodingModeProp_;
    mutable IntOptionProperty colorTableProp_;
    mutable ColorProperty colorProp_;
    mutable TransFunc1DKeysProperty transFuncProp_;
};

}   // namespace

#endif
