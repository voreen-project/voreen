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

#ifndef VRN_TOUCHTABLECAMERAWIDGET_H
#define VRN_TOUCHTABLECAMERAWIDGET_H


#include "../touchtableoverlay.h"
#include "touchtablemenuwidgetscrollable.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "modules/base/processors/image/orientationoverlay.h"


#include "../touchtablecontrolelement.h"
#include "../touchtablemenuframe.h"

#include "tgt/filesystem.h"
#include "tgt/quaternion.h"

#include "tgt/timer.h"
#include "voreen/core/voreenapplication.h"
#include "tgt/event/timeevent.h"

namespace voreen {

    /// standard orientations for the camera
    enum CameraOrientation {
        AXIAL, CORONAL, SAGITTAL, AXIAL_INVERSE, CORONAL_INVERSE, SAGITTAL_INVERSE
    };

    /**
    * This widget configures settings of a linked camera property (eg. standard perspectives).
    * Camera positions and orientations may be bookmarked.
    * Additionally, rendering of an orientation cube or coordinate axes is provided through linking with an OrientationOverlay processor.
    */
    class TouchTableCameraWidget : public TouchTableMenuWidgetScrollable {

        ///struct for displaying messages when loading / saving a preset
        struct Message {

            Message() {
                message_ = "";
                display_ = false;
            }

            std::string message_;
            bool display_;
        };

    public:

        TouchTableCameraWidget();

        virtual ~TouchTableCameraWidget();

        Processor* create() const {
            return new TouchTableCameraWidget();
        }

        virtual bool isReady() const {
            return port_.isConnected();
        }

        virtual std::string getClassName() const {
            return "TouchTableCameraWidget";
        }

        /**
         * Handles touch points within the camera overlay menu.
         *
         * @param tp a std::deque containing all the touch points that should be handled by the widget.
         */
        virtual void handleTouchPoints(const std::deque<tgt::TouchPoint>& tp);

        /**
         * Renders the control elements inside the camera overlay menu.
         */
        virtual void renderComponents();

        /**
         * Handles timer events for selecting default camera orientations, continuous rotation and display messages.
         */
        virtual void timerEvent(tgt::TimeEvent* te);

        ///load camera preset
        virtual void handleScrollableMenuSelection(std::string file);

        /**
         * @see TouchTableMenuWidgetScrollable
         */
        virtual void updateScrollableMenuPosition();

        /**
         * @see TouchTableMenuWidgetScrollable
         */
        virtual void updateScrollableMenuContent();

    protected:

        void setDescriptions(){
            setDescription("Camera Widget for managing standard camera orientations, continuous rotation and camera presets.\n \
                    Link bool properties for rendering orientation box and axes with OrientationOverlay processor and camera with raycaster / proxy geometry.");
        }

        virtual void initialize();

        virtual void deinitialize();

        /**
         * Is called if properties change to update the control element attributes.
         * @see TouchTableMenuWidgetScrollable
         */
        virtual void updateComponents();

        /**
         * Save the current camera property to a file.
         * The path is set by the file dialog property, the file name is constructed from the current date and time.
         */
        virtual void saveCamera();

        /**
         * Changes the camera orientation using 10 animation steps while keeping the trackball distance.
         */
        virtual void orientationChanged(CameraOrientation o);

        /**
         * Checks if continuous rotation is enabled and starts / stops the corresponding timer.
         */
        virtual void setRotationTimerState();

        CameraProperty cam_;                                ///< camera property (to be linked with raycaster and other processors in the voreen network)

        BoolProperty renderOrientationBox_;                 ///< determines if a textured cube is rendered for camera orientation
        BoolProperty renderOrientationAxes_;                ///< determines if the coordinate axes are rendered for camera orientation
        BoolProperty enableOrientationProp_;                                        ///< used to disable orientation overlay (must be linked)
        OptionProperty<OrientationOverlay::OrientationType> orientationTypeProp_;   ///< determining the current used overlay

        FileDialogProperty presetDirectory_;                ///< property for selecting a directory in which to search for camera presets

        bool rotateX_;                                      ///< is continuous rotation around x axis enabled?
        bool rotateY_;                                      ///< is continuous rotation around y axis enabled?
        bool rotateZ_;                                      ///< is continuous rotation around z axis enabled?

        //control elements that control standard camera positions
        TouchTableControlElement coronal_;
        TouchTableControlElement sagittal_;
        TouchTableControlElement axial_;
        TouchTableControlElement coronalReverse_;
        TouchTableControlElement sagittalReverse_;
        TouchTableControlElement axialReverse_;

        //control elements for orientation rendering
        TouchTableControlElement orientationBox_;
        TouchTableControlElement orientationAxes_;

        //control elements for rotation around the axes
        TouchTableControlElement rotateXElem_, rotateYElem_, rotateZElem_;

        TouchTableControlElement loadPreset_;                ///< control element for loading presets
        TouchTableControlElement savePreset_;                ///< control element for savong the current camera as a preset

        //textures for control elements
        tgt::Texture* axialTex_, * axialInvTex_;
        tgt::Texture* sagittalTex_, * sagittalInvTex_;
        tgt::Texture* coronalTex_, * coronalInvTex_;

        tgt::Texture* orientationBoxTex_, * orientationAxesTex_;

        tgt::Texture* rotationXTex_, * rotationYTex_, * rotationZTex_;

        tgt::Texture* presetTex_, * saveTex_;

        tgt::Timer* orientationTimer_;                  ///< Timer object for setting default camera orientations
        tgt::Timer* rotationTimer_;                     ///< Timer object for rotating the camera around the coordinate axes
        tgt::Timer* messageTimer_;                      ///< Timer object for displaying save / load messages
        tgt::EventHandler eventHandler_;                ///< A local eventhandle which is added to the timer

        CameraOrientation currentRotation_;             ///< current orientation for animation
        tgt::quat oldQuat_;                             ///< old quaternion for rotation animation

        std::vector<std::pair<std::string, tgt::Texture*> > cameraPresetFiles_;     ///< contains all camera preset file names availabe in the preset menu

        Message fileMessage_;                           ///< message displayed for load / save actions

        static const std::string loggerCat_;
    };

} // namespace

#endif // VRN_TOUCHTABLECAMERAWIDGET_H
