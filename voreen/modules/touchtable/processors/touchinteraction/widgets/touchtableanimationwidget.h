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

#ifndef VRN_TOUCHTABLESANIMATIONWIDGET_H
#define VRN_TOUCHTABLESANIMATIONWIDGET_H

#include "../touchtableoverlay.h"
#include "touchtablemenuwidget.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/fontproperty.h"
#include "voreen/core/animation/animatedprocessor.h"
#include "voreen/core/animation/animation.h"
#include "voreen/core/animation/templatepropertytimeline.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"

#include "../touchtablecontrolelement.h"
#include "../touchtablemenuframe.h"
#include "../touchtablescrollabletimelinemenu.h"

#include "../touchtableslider.h"

namespace voreen {

    /**
    * TouchTableAnimationWidget inherits from TouchTableMenuWidget. When active through this widget
    * animations can be managed. States of the Clipping, Transfunc Widget, the camera etc. can be recorded
    * and be animated through a certain interpolation function between the recorded states
    */
    class TouchTableAnimationWidget : public TouchTableMenuWidget {

    public:

        TouchTableAnimationWidget();

        ~TouchTableAnimationWidget();

        Processor* create() const {
            return new TouchTableAnimationWidget();
        }

        virtual bool isReady() const {
            return port_.isConnected();
        }

        virtual std::string getClassName() const {
            return "TouchTableAnimationWidget";
        }

        /**
         * Handles touch points within the animation overlay menu
         *
         * @param tp a std::deque containing all the touch points that should be handled by the widget.
         * @param overlay the TouchTableOverlay processor which has been catching the touch point events before.
         */
        virtual void handleTouchPoints(const std::deque<tgt::TouchPoint>& tp);

    protected:

        void setDescriptions(){
            setDescription("manages animations of changes between recorded states");
        }

        /**
         * Renders the animation overlay menu, the control elements, the ti
         */
        virtual void renderComponents();

        /**
        * Renders the current Time
        */
        virtual void renderTime();

        /**
        * Renders the scrollable timeline with content
        */
        virtual void renderScrollableTimeline();

        virtual void initialize();

        virtual void deinitialize();

         /**
         * Is called if properties change to update the control element attributes.
         */
        virtual void updateComponents();

        //update functions called, when properties change
        /**
        * updates the timeline GUI
        */
        virtual void updateTimeline();


        /**
        * set the selection handler of the timelineMenu. Here: When touchpoint in timeline, propertystate for given time is being rendered
        *
        *@param time: time at which the propertystate is being rendered
        *@param bool: true if a key value has been hit, false otherwise
        */
        virtual void handleSelection(boost::tuple<float, bool, float>);

        /**
        * updates the slider, according to cursor position in the timeline
        */
        virtual void updateTimelineSlider();

        /**
        * deletes a keyvalue at the given time.
        *
        *@param time= time at which keyvalue shall be deleted
        */
        virtual void deleteKeyValue(float time);

        /**
        * handlers TimerEvents called by the animationTimer
        *
        *@param te TimeEvent to be handled
        */
        virtual void timerEvent(tgt::TimeEvent* te);

        /**
        * sets up an interface for a new Animation
        */
        virtual void setUpNewAnimationMode();

        /**
        * sets up the new animation Object with all its propertyTimelines
        */
        virtual void setUpNewAnimation();

        /**
        * renders all text displayed, when new animatoin is being set up
        */
        virtual void renderNewAnimationText();

        //properties
        FontProperty fontProperty_;
        IntProperty timelineDuration_;
        IntProperty timelineStepSize_;
        FloatProperty fps_;
        FloatProperty windFactor_;
        //animatedProperties
        TransFunc1DKeysProperty transfuncProp_;
        CameraProperty cameraProp_;
        BoolProperty toggleClipping0Prop_, toggleClipping1Prop_, toggleClipping2Prop_;
        FloatVec3Property planeNormal0Prop_, planeNormal1Prop_, planeNormal2Prop_;
        FloatProperty planePosition0Prop_, planePosition1Prop_, planePosition2Prop_;

        //textures for control elements
        tgt::Texture* resetTex_, /* * undoTex_,* redoTex_,*/ * jumpEndTex_, * jumpBeginningTex_, * rewindTex_, * fastforwardTex_, * pauseTex_, *playTex_,
            * recordTex_, * cameraTex_, * clippingTex_, * transfuncTex_, *deleteTex_, /**exportVideoTex_, */ *confirmTex_;

        //control elements
        TouchTableControlElement resetElem_, /*undoElem_, redoElem_,*/ jumpEndElem_, jumpBeginningElem_, rewindElem_, fastforwardElem_, pauseElem_,
            playElem_, recordElem_, cameraElem_, clippingElem_, transfuncElem_, deleteElem_, /*exportVideoElem_, */ confirmElem_;

        bool newAnimationMode_;                                ///< set when resetElem has been hit, to set up a new animation

        TouchTableSlider durationSlider_;                    ///< slider to set up the duration of the new Animation

        std::vector<int> durationSliderIds_;

        //menu with current time
         TouchTableMenuFrame timeMenu_;
        std::string currentTime_;

        //mapping of touchpoints and control elements
        std::map<int, TouchTableControlElement*> tpControlElementMapping_;

        //control element vector containing all control elements
        std::vector<TouchTableControlElement*> controlElements_;

        //scrollable menu containing the timeline
        TouchTableScrollableTimelineMenu timelineMenu_;

        //propertytimelines managing the animation of the animated properties
        PropertyTimelineTransFunc1DKeys* transfuncTimeline_;                //< pointer to timeline for transfer function
        PropertyTimelineCamera* cameraTimeline_;                    //< pointer to timeline for camera
        //clipping timelines
        std::vector<PropertyTimelineFloat*> clippingPlanePositionTimelines_;
        std::vector<PropertyTimelineVec3*> clippingPlaneNormalTimelines_;
        std::vector<PropertyTimelineBool*> clippingToggleTimelines_;

        //timer and eventHanlder to call/handle animation events
        tgt::Timer* playTimer_;
        tgt::Timer* fastforwardTimer_;
        tgt::Timer* rewindTimer_;
        tgt::EventHandler eventHandler_;


        //animation class organising the animations
        Animation* animation_;

        ProcessorNetwork* network_;
        NetworkEvaluator* networkEval_;

        //processors camera, transfer function or clipping processors to be animated
        std::vector<AnimatedProcessor*> animatedProcessors_;

        static const std::string loggerCat_;
    };

} // namespace

#endif // VRN_TOUCHTABLEANIMATONWIDGET_H
