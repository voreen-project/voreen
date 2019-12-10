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

#include "pblinkcontrol.h"

namespace voreen {

const std::string PBLinkControl::loggerCat_("voreen.bovenkamp.PBLinkControl");

PBLinkControl::PBLinkControl()
    : Processor()
    , renderMagnitudeOverlay_("renderMagnitudeOverlay","Render Magnitude Overlay",false)
    , renderColorCodingOverlay_("renderColorCodingOverlay","Render Color Coding Overlay",false)
    , colorProp_("colorProp","Color:")
{
    //overlay
    renderMagnitudeOverlay_.setVisibleFlag(false);
    renderColorCodingOverlay_.setVisibleFlag(false);
    addProperty(renderMagnitudeOverlay_);
    addProperty(renderColorCodingOverlay_);

    colorProp_.addOption("velocity" , "Veocity" , StreamlineRenderer3D::COLOR_VELOCITY);
    colorProp_.addOption("direction" , "Direction" , StreamlineRenderer3D::COLOR_DIRECTION);
    addProperty(colorProp_);

    colorProp_.onChange(MemberFunctionCallback<PBLinkControl>(this, &PBLinkControl::updateOverlayCallback));



}

Processor* PBLinkControl::create() const {
    return new PBLinkControl();
}

void PBLinkControl::process() { }

void PBLinkControl::updateOverlayCallback() {
    if (colorProp_.get() == "velocity")
        renderMagnitudeOverlay_.set(true);
    else
        renderMagnitudeOverlay_.set(false);

    if (colorProp_.get() == "direction")
        renderColorCodingOverlay_.set(true);
    else
        renderColorCodingOverlay_.set(false);
}

} // namespace voreen

