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

#include "colorcodingability.h"

#include <string>
#include <vector>

namespace voreen {

ColorCodingAbility::ColorCodingAbility()
    : colorCodingModeProp_("colorCoding", "color coding: "),
      colorTableProp_("colorTableProp_", "color table: "),
      colorProp_("colorProp_", "color: ", tgt::vec4(1.0f)),
      transFuncProp_("transFuncProp_", "Color Texture")
{
    colorCodingModeProp_.addOption("magnitude", "magnitude", COLOR_CODING_MODE_MAGNITUDE);
    colorCodingModeProp_.addOption("direction", "direction", COLOR_CODING_MODE_DIRECTION);
    colorCodingModeProp_.addOption("monochrome", "monochrome", COLOR_CODING_MODE_MONOCHROME);
    colorCodingModeProp_.onChange(
        MemberFunctionCallback<ColorCodingAbility>(this, &ColorCodingAbility::onColorCodingModeChange));

    colorTableProp_.addOption("rainbow", "rainbow", 0);
    colorTableProp_.addOption("hot metal", "hot metal", 1);
    colorTableProp_.addOption("texture", "texture", 2);
}

ColorCodingAbility::~ColorCodingAbility() {
}

std::string ColorCodingAbility::getShaderDefines() const {
    std::ostringstream oss;
    oss << "#define COLOR_MODE " << static_cast<int>(colorCodingModeProp_.getValue()) << "\n";
    if (colorCodingModeProp_.getValue() == COLOR_CODING_MODE_MAGNITUDE)
        oss << "#define COLOR_TABLE " << colorTableProp_.getValue() << "\n";
    return oss.str();
}

void ColorCodingAbility::onColorCodingModeChange() {
    switch (colorCodingModeProp_.getValue()) {
        default:
        case COLOR_CODING_MODE_MAGNITUDE:
            colorProp_.setVisibleFlag(false);
            colorTableProp_.setVisibleFlag(true);
            transFuncProp_.setVisibleFlag(true);
            break;
        case COLOR_CODING_MODE_DIRECTION:
            colorProp_.setVisibleFlag(false);
            colorTableProp_.setVisibleFlag(false);
            transFuncProp_.setVisibleFlag(false);
            break;
        case COLOR_CODING_MODE_MONOCHROME:
            colorProp_.setVisibleFlag(true);
            colorTableProp_.setVisibleFlag(false);
            transFuncProp_.setVisibleFlag(false);
            break;
    }
}

}   // namespace
