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

#ifndef VRN_SLICECAMERAINTERACTIONHANDLER_H
#define VRN_SLICECAMERAINTERACTIONHANDLER_H

#include "voreen/core/interaction/interactionhandler.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/optionproperty.h"

#include "tgt/event/eventhandler.h"

#include <bitset>

namespace voreen {

class CameraProperty;

class VRN_CORE_API SliceCameraInteractionHandler : public InteractionHandler {
    friend class EventProperty<SliceCameraInteractionHandler>;

public:
    /// Default constructor needed for serialization. Do not call it directly.
    SliceCameraInteractionHandler();

    /**
     * Constructor.
     *
     * @param id Identifier that must be unique across all interaction handlers
     *  of a processor. Must not be empty.
     * @param guiText the string that is to be displayed in the GUI
     * @param cameraProp camera property the handler will operate on. Must not be null.
     */
    SliceCameraInteractionHandler(const std::string& id, const std::string& guiText, CameraProperty* cameraProp,
        bool sharing = false, bool enabled = true);

    virtual ~SliceCameraInteractionHandler();

    virtual std::string getClassName() const   { return "SliceCameraInteractionHandler";     }
    virtual InteractionHandler* create() const { return new SliceCameraInteractionHandler(); }

private:
    /// @see InteractionHandler::onEvent
    virtual void onEvent(tgt::Event* e);

    // functions called by the event properties
    void zoomEvent(tgt::MouseEvent* e);
    void shiftEvent(tgt::MouseEvent* e);
    void handleMultitouch(tgt::TouchEvent* e);

    // camera property the handler operates on
    CameraProperty* cameraProp_;

    // trackball properties
    EventProperty<SliceCameraInteractionHandler>* zoomEvent_;
    EventProperty<SliceCameraInteractionHandler>* shiftEvent_;
    EventProperty<SliceCameraInteractionHandler>* wheelZoomEvent_;
    EventProperty<SliceCameraInteractionHandler>* multiTouchEvent_;

    // Stores if a mousebutton (LEFT, MIDDLE, RIGHT) has been pressed but not released yet.
    std::bitset<3> pressedMouseButtons_;

    tgt::ivec2 lastMousePos_;
    float lastDistance_;
};

} // namespace

#endif
