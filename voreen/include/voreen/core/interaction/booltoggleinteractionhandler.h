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

#ifndef VRN_BOOLTOGGLEINTERACTIONHANDLER_H
#define VRN_BOOLTOGGLEINTERACTIONHANDLER_H

#include "voreen/core/interaction/interactionhandler.h"
#include "voreen/core/properties/boolproperty.h"

#include "tgt/event/keyevent.h"

namespace voreen {

class VRN_CORE_API BoolToggleInteractionHandler: public InteractionHandler {

public:
    /// Default constructor needed for serialization. Do not call it directly!
    BoolToggleInteractionHandler();

    /**
     * Constructor.
     *
     * @param id Identifier that must be unique across all interaction handlers
     *  of a processor. Must not be empty.
     * @param name the string that is to be displayed in the GUI
     * @param boolProp numeric property the handler will operate on. Must not be null.
     * @param The modifier that has to be active for allowing the handler to events
     * @param sharing @see InteractionHandler
     * @param enabled @see InteractionHandler
     */
    BoolToggleInteractionHandler(const std::string& id, const std::string& name, BoolProperty* boolProp, tgt::KeyEvent::KeyCode keyCode,
        tgt::Event::Modifier modifier = tgt::Event::MODIFIER_NONE, bool sharing = false, bool enabled = true);

    virtual ~BoolToggleInteractionHandler() {}

    virtual std::string getClassName() const   { return "BoolToggleInteractionHandler";     }
    virtual InteractionHandler* create() const { return new BoolToggleInteractionHandler(); }

private:
    virtual void onEvent(tgt::Event* e);

    BoolProperty* boolProp_;
};

} // namespace

#endif
