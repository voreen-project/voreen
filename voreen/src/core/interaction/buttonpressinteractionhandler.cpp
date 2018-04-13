/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "voreen/core/interaction/buttonpressinteractionhandler.h"

using tgt::Event;

namespace voreen {

voreen::ButtonPressInteractionHandler::ButtonPressInteractionHandler() :
    InteractionHandler("dummy", "dummy")
{ }

ButtonPressInteractionHandler::ButtonPressInteractionHandler(const std::string& id, const std::string& name,
        ButtonProperty* buttonProp, tgt::KeyEvent::KeyCode keyCode, tgt::Event::Modifier modifier, bool sharing, bool enabled)
    : InteractionHandler(id, name, keyCode, modifier, sharing, enabled)
{
    tgtAssert(buttonProp, "No button property");
    buttonProp_ = buttonProp;
}

void ButtonPressInteractionHandler::onEvent(tgt::Event* eve) {
    tgtAssert(buttonProp_, "No button property");

    if (tgt::KeyEvent* keyEve = dynamic_cast<tgt::KeyEvent*>(eve)) {
        if(!keyEve->pressed())
            buttonProp_->clicked();

        eve->accept();
    }
}

} // namespace
