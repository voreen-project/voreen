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

#ifndef VRN_CAMERAINTERACTIONHANDLER_H
#define VRN_CAMERAINTERACTIONHANDLER_H

#include "voreen/core/datastructures/geometry/meshlistgeometry.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"
#include "voreen/core/interaction/interactionhandler.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/optionproperty.h"

#include "tgt/event/eventhandler.h"

#include <bitset>

namespace tgt {
    class Timer;
}

namespace voreen {

class CameraProperty;
class TrackballNavigation;
class FirstPersonNavigation;

/**
 * Interaction handler for camera interactions providing
 * two navigation modes: trackball and first-person.
 *
 * It is possible to use multiple camera interaction handlers
 * in a single processor, but then you have to make sure that
 * they have unique names.
 *
 * @see TrackballNavigation
 * @see FirstPersonNavigation
 */
class VRN_CORE_API CameraInteractionHandler : public InteractionHandler {

    friend class EventProperty<CameraInteractionHandler>;

public:
    /**
     * Used to specify the current camera interaction behaviour.
     */
    enum NavigationMetaphor {
        TRACKBALL,
        FIRST_PERSON
    };

    /// Default constructor needed for serialization. Do not call it directly.
    CameraInteractionHandler();

    /**
     * Constructor.
     *
     * @param id Identifier that must be unique across all interaction handlers
     *  of a processor. Must not be empty.
     * @param guiText the string that is to be displayed in the GUI
     * @param cameraProp camera property the handler will operate on. Must not be null.
     */
    CameraInteractionHandler(const std::string& id, const std::string& guiText, CameraProperty* cameraProp,
        bool sharing = false, bool enabled = true);

    virtual ~CameraInteractionHandler();

    virtual std::string getClassName() const   { return "CameraInteractionHandler";     }
    virtual InteractionHandler* create() const { return new CameraInteractionHandler(); }

    /**
     * Determines the the navigation metaphor to be used
     * for interpreting the user input.
     */
    void setNavigationMode(NavigationMetaphor behaviour);

    /**
     * Returns the currently used navigation mode.
     */
    NavigationMetaphor getNavigationMode() const;

    /**
     * @see InteractionHandler::setVisible
     */
    void setVisible(bool state);

    /**
     * Should the trackball center shift with the camera?
     */
    void adjustCenterShift();

    /**
     * After shifting, reset the camera focus to the center of the trackball.
     */
    void resetCameraFocusToTrackballCenter();

private:

    /// @see InteractionHandler::onEvent
    virtual void onEvent(tgt::Event* e);

    // functions called by the event properties
    void rotateEvent(tgt::MouseEvent* e);
    void zoomEvent(tgt::MouseEvent* e);
    void shiftEvent(tgt::MouseEvent* e);
    void keyEvent(tgt::KeyEvent* e);
    void handleMultitouch(tgt::TouchEvent* e);

    // Writes the states of the fp event properties to the interaction handler.
    void updateFpKeySettings();

    /**
     * Adjusts the event property widget's visibility.
     * Called after the navigation mode has changed.
     */
    void adjustWidgetStates();

    // camera property the handler operates on
    CameraProperty* cameraProp_;

    // trackball <-> first-person navigation
    OptionProperty<NavigationMetaphor> navigationMetaphor_;

    // trackball properties
    EventProperty<CameraInteractionHandler>* rotateEvent_;
    EventProperty<CameraInteractionHandler>* contRotateEvent_;
    EventProperty<CameraInteractionHandler>* zoomEvent_;
    EventProperty<CameraInteractionHandler>* shiftEvent_;
    EventProperty<CameraInteractionHandler>* wheelZoomEvent_;
    EventProperty<CameraInteractionHandler>* multiTouchEvent_;

    // first-person properties
    EventProperty<CameraInteractionHandler>* moveForwardEvent_;
    EventProperty<CameraInteractionHandler>* moveLeftEvent_;
    EventProperty<CameraInteractionHandler>* moveBackwardEvent_;
    EventProperty<CameraInteractionHandler>* moveRightEvent_;
    EventProperty<CameraInteractionHandler>* moveUpEvent_;
    EventProperty<CameraInteractionHandler>* moveDownEvent_;

    // navigation objects
    TrackballNavigation* tbNavi_;
    FirstPersonNavigation* fpNavi_;

    // Needed for FirstPersonNavigation for continuous movement
    tgt::Timer* motionTimer_;
    tgt::EventHandler* timerEventHandler_;

    // Stores if a mousebutton (LEFT, MIDDLE, RIGHT) has been pressed but not released yet.
    std::bitset<3> pressedMouseButtons_;
};

} // namespace

#endif // VRN_INTERACTIONHANDLER_H
