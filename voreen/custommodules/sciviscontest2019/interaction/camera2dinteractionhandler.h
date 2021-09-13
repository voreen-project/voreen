/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2014 University of Muenster, Germany.                        *
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

#ifndef VRN_CAMERA2DINTERACTIONHANDLER_H
#define VRN_CAMERA2DINTERACTIONHANDLER_H

#include "voreen/core/interaction/interactionhandler.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/vectorproperty.h"

#include "tgt/event/eventhandler.h"

#include <bitset>

namespace voreen {

class CameraProperty;

/**
 * Class that implements Interaction with a 2D camera defined by a center and a window width
 */
class VRN_CORE_API Camera2DInteractionHandler : public InteractionHandler {
    friend class EventProperty<Camera2DInteractionHandler>;

public:
    /// Default constructor needed for serialization. Do not call it directly.
    Camera2DInteractionHandler();

    /**
     * Constructor.
     *
     * @param id Identifier that must be unique across all interaction handlers
     *  of a processor. Must not be empty.
     * @param guiText the string that is to be displayed in the GUI
     * @param center Center of the camera view rectangle, in world coordinates.
     * @param width Width of the camera view rectangle, in world coordinates.
     * @param outportWidth Width of outport. Used to determine exact mouse positions
     *
     */
    Camera2DInteractionHandler(const std::string& id, const std::string& guiText, FloatVec2Property* center, FloatProperty* width, const size_t* outPortWidth, bool sharing = false, bool enabled = true);

    virtual ~Camera2DInteractionHandler();

    virtual std::string getClassName() const   { return "Camera2DInteractionHandler";     }
    virtual InteractionHandler* create() const { return new Camera2DInteractionHandler(); }

private:
    /// @see InteractionHandler::onEvent
    virtual void onEvent(tgt::Event* e);

    // functions called by the event properties
    void zoomEvent(tgt::MouseEvent* e);
    void shiftEvent(tgt::MouseEvent* e);
    void multiTouchEvent(tgt::TouchEvent* e);

    /// Center of the camera view rectangle, in world coordinates.
    FloatVec2Property* center_;
    /// Width of the camera view rectangle, in world coordinates.
    FloatProperty* width_;
    /// Width Width of outport. Used to determine exact mouse positions
    const size_t* outPortWidth_;

    // Event properties
    EventProperty<Camera2DInteractionHandler>* multiTouchEvent_;
    EventProperty<Camera2DInteractionHandler>* zoomEvent_;
    EventProperty<Camera2DInteractionHandler>* shiftEvent_;
    EventProperty<Camera2DInteractionHandler>* wheelZoomEvent_;

    // Stores if a mousebutton (LEFT, MIDDLE, RIGHT) has been pressed but not released yet.
    std::bitset<3> pressedMouseButtons_;

    /// Last mouse position (e.g. to shift the camera)
    tgt::ivec2 lastMousePos_;
    /// Last distance between to points. Used for 2 finger zoom
    float lastDistance_;
};

} // namespace

#endif
