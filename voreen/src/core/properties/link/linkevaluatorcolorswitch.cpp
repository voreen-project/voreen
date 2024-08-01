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

#include "voreen/core/properties/link/linkevaluatorcolorswitch.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/colorproperty.h"

namespace voreen {

void LinkEvaluatorColorSwitchBool::eval(Property* src, Property* dst) {
    if (dynamic_cast<BoolProperty*>(src)){
        ColorSwitchProperty* colorprop = dynamic_cast<ColorSwitchProperty*>(dst);
        BoolProperty* boolprop = dynamic_cast<BoolProperty*>(src);
        colorprop->setUseActiveColor(boolprop->get());

    }else{
        ColorSwitchProperty* colorprop = dynamic_cast<ColorSwitchProperty*>(src);
        BoolProperty* boolprop = dynamic_cast<BoolProperty*>(dst);
        boolprop->set(colorprop->getUseActiveColor());
    }
}

bool LinkEvaluatorColorSwitchBool::arePropertiesLinkable(const Property* src, const Property* dst) const
{
    if (dynamic_cast<const BoolProperty*>(src) && dynamic_cast<const ColorSwitchProperty*>(dst))
        return true;

    if (dynamic_cast<const BoolProperty*>(dst) && dynamic_cast<const ColorSwitchProperty*>(src))
        return true;
    return false;
}


/*void LinkEvaluatorColorSwitchColor::eval(Property* src, Property* dst) {

}*/

/*bool LinkEvaluatorColorSwitchColor::arePropertiesLinkable(const Property* p1, const Property* p2) const
{

}*/

} // namespace
