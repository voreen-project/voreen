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

#include "touchtableboolwidget.h"
#include "../touchtableoverlay.h"


namespace voreen {

TouchTableBoolWidget::TouchTableBoolWidget()
    : TouchTableWidget()
    , boolProp_("boolProp", "Bool Property", false)
{
    addProperty(boolProp_);
    boolProp_.onChange(MemberFunctionCallback<TouchTableBoolWidget>(this, &TouchTableBoolWidget::updateWidget));
}


Processor* TouchTableBoolWidget::create() const {
    return new TouchTableBoolWidget();
}

std::string TouchTableBoolWidget::getClassName() const {
    return "TouchTableBoolWidget";
}

void TouchTableBoolWidget::setDescriptions() {
    setDescription("Widget for use with TouchScreenOverlay. Pressing and releasing it triggers a bool property that may be linked to other processors in the\
            Voreen network to trigger functionality.");
    boolProp_.setDescription("Link this with the property to like to manipulate with this widget.");
}

void TouchTableBoolWidget::pressed() {
    //trigger bool property and set widget active / inactive depending on its state
    boolProp_.set(!boolProp_.get());
    //isActive_.set(boolProp_.get()); done by callback
}

//----------------------------------------------------------------------------------------------
//  Callbacks
//----------------------------------------------------------------------------------------------
void TouchTableBoolWidget::updateWidget() {
    isActive_.set(boolProp_.get());
}

} // namespace
