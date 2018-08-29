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

#include "touchtablemovewidget.h"
#include "../touchtableoverlay.h"


namespace voreen {

TouchTableMoveWidget::TouchTableMoveWidget()
    : TouchTableWidget()
{
    symbolFile_.setVisibleFlag(false);
}


Processor* TouchTableMoveWidget::create() const {
    return new TouchTableMoveWidget();
}

std::string TouchTableMoveWidget::getClassName() const {
    return "TouchTableMoveWidget";
}

void TouchTableMoveWidget::setDescriptions() {
    setDescription("Widget for use with TouchScreenOverlay. Pressing and releasing it switches the movement mode of the overlay processor, \
alternating between changing the camera position and changing the camera orientation.");
}

void TouchTableMoveWidget::initialize() {
    Processor::initialize();

    //set default texture
    TexMgr.addPath(VoreenApplication::app()->getModulePath("touchtable") + "/textures");
    symbolTex_=TexMgr.load("moveAction.png");
    symbolFile_.set(VoreenApplication::app()->getModulePath("touchtable") + "/textures" + "/moveAction.png");
}

void TouchTableMoveWidget::pressed() {
    //trigger movement state and set widget active / inactive depending on its state
    if (overlay_) {
        overlay_->setCameraMovement(!overlay_->getCameraMovement());
        isActive_.set(overlay_->getCameraMovement());
    }
}

void TouchTableMoveWidget::setOverlay(TouchTableOverlay* overlay) {
    overlay_ = overlay;
    if (overlay)
        isActive_.set(overlay->getCameraMovement());
    else
        isActive_.set(false);
}

} // namespace
