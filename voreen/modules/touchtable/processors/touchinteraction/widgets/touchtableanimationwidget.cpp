/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "touchtableanimationwidget.h"
#include <time.h>
#include "voreen/core/voreenapplication.h"
#include "voreen/core/network/processornetwork.h"
#include "voreen/core/network/networkevaluator.h"
#include "voreen/core/animation/serializationfactories.h"
#include "voreen/core/animation/templatepropertytimeline.h"


namespace voreen {


const std::string TouchTableAnimationWidget::loggerCat_("voreen.touchtable.TouchTableAnimationWidget");


TouchTableAnimationWidget::TouchTableAnimationWidget()
    : TouchTableMenuWidget("animation.png")
    , resetElem_(25, tgt::ivec2(10 + 20, 175),0, 1.0f, tgt::vec4(0.0f), tgt::vec4(0.0f), 0, false, false)
    //, undoElem_(25, tgt::ivec2(10 + 20, 175),0, 1.0f, tgt::vec4(0.0f), tgt::vec4(0.0f), 1, false, false)
    //, redoElem_(25, tgt::ivec2(10 + 20, 175),0, 1.0f, tgt::vec4(0.0f), tgt::vec4(0.0f), 2, false, false)
    , jumpEndElem_(25, tgt::ivec2(10 + 20, 175),0, 1.0f, tgt::vec4(0.0f), tgt::vec4(0.0f), 3, false, false)
    , jumpBeginningElem_(25, tgt::ivec2(10 + 20, 175),0, 1.0f, tgt::vec4(0.0f), tgt::vec4(0.0f), 4, false, false)
    , rewindElem_(25, tgt::ivec2(10 + 20, 175),0, 1.0f, tgt::vec4(0.0f), tgt::vec4(0.0f), 5, false, true)
    , fastforwardElem_(25, tgt::ivec2(10 + 20, 175),0, 1.0f, tgt::vec4(0.0f), tgt::vec4(0.0f), 6, false, true)
    , pauseElem_(25, tgt::ivec2(10 + 20, 175),0, 1.0f, tgt::vec4(0.0f), tgt::vec4(0.0f), 7, false, true)
    , playElem_(25, tgt::ivec2(10 + 20, 175),0, 1.0f, tgt::vec4(0.0f), tgt::vec4(0.0f), 8, false, true)
    , recordElem_(25, tgt::ivec2(10 + 20, 175),0, 1.0f, tgt::vec4(0.0f), tgt::vec4(0.0f), 9, false, false)
    , cameraElem_(25, tgt::ivec2(10 + 20, 175),0, 1.0f, tgt::vec4(0.0f), tgt::vec4(0.0f), 10, false, true)
    , clippingElem_(25, tgt::ivec2(10 + 20, 175),0, 1.0f, tgt::vec4(0.0f), tgt::vec4(0.0f), 11, false, true)
    , transfuncElem_(25, tgt::ivec2(10 + 20, 175),0, 1.0f, tgt::vec4(0.0f), tgt::vec4(0.0f), 12, false, true)
    , deleteElem_(25, tgt::ivec2(10 + 20, 175),0, 1.0f, tgt::vec4(0.0f), tgt::vec4(0.0f), 13, false, false)
    //, exportVideoElem_(25, tgt::ivec2(10 + 20, 175),0, 1.0f, tgt::vec4(0.0f), tgt::vec4(0.0f), 13, false, false)
    , confirmElem_(25, tgt::ivec2(10 + 20, 175),0, 1.0f, tgt::vec4(0.0f), tgt::vec4(0.0f), 13, false, false)
    , timeMenu_(tgt::ivec2(0), tgt::vec4(1.0f, 1.0f, 1.0f, 0.4f))
    , fontProperty_("fontProperty", "Font Property")
    , timelineDuration_("timelineDuration", "Timeline Duration", 900, 900, 15000)
    , timelineStepSize_("timelineStepSize", "Timeline Stepsize", 60, 10, 100)
    , fps_("fps", "FPS", 20, 1, 60)
    , windFactor_("windFactor", "Wind Factor", 4, 1, 20)
    , timelineMenu_(getPosition(), tgt::vec4(0.1f, 0.1f, 0.1f, 1.f), 500 ,40)
    , currentTime_("00:00")
    , transfuncProp_("transfunc", "Transfer Function")
    , cameraProp_("camera", "Camera")
    , toggleClipping0Prop_("toggleClipping0",  "Toggle Clipping 0")
    , toggleClipping1Prop_("toggleClipping1", "Toggle Clipping 1")
    , toggleClipping2Prop_("toggleClipping2", "Toggle Clipping 2")
    , planeNormal0Prop_("planeNormal0" , "Plane Normal 0",tgt::vec3(0.f), tgt::vec3(-1.1f), tgt::vec3(1.1f))
    , planeNormal1Prop_("planeNormal1" , "Plane Normal 1",tgt::vec3(0.f), tgt::vec3(-1.1f), tgt::vec3(1.1f))
    , planeNormal2Prop_("planeNormal2" , "Plane Normal 2",tgt::vec3(0.f), tgt::vec3(-1.1f), tgt::vec3(1.1f))
    , planePosition0Prop_("planePosition0", "Plane Position 0", 0.f, -FLT_MAX, FLT_MAX)
    , planePosition1Prop_("planePosition1", "Plane Position 1", 0.f, -FLT_MAX, FLT_MAX)
    , planePosition2Prop_("planePosition2", "Plane Position 2", 0.f, -FLT_MAX, FLT_MAX)
    , playTimer_(0)
    , fastforwardTimer_(0)
    , rewindTimer_(0)
    , eventHandler_()
    , animation_(0)
    , network_(0)
    , networkEval_(0)
    , newAnimationMode_(true)
    , durationSlider_(tgt::ivec2(100,100), 10, 10, HORIZONTAL, 25, 25, 0)
{

    //set timer
    eventHandler_.addListenerToBack(this);
    playTimer_ = VoreenApplication::app()->createTimer(&eventHandler_);
    fastforwardTimer_ = VoreenApplication::app()->createTimer(&eventHandler_);
    rewindTimer_ = VoreenApplication::app()->createTimer(&eventHandler_);

    //save elements in vector for easier access, when modifications are the same for all control elements
    controlElements_.push_back(&transfuncElem_);
    controlElements_.push_back(&clippingElem_);
    controlElements_.push_back(&cameraElem_);
    controlElements_.push_back(&recordElem_);
    controlElements_.push_back(&playElem_);
    controlElements_.push_back(&pauseElem_);
    controlElements_.push_back(&fastforwardElem_);
    controlElements_.push_back(&rewindElem_);
    controlElements_.push_back(&jumpBeginningElem_);
    controlElements_.push_back(&jumpEndElem_);
    //controlElements_.push_back(&redoElem_);
    //controlElements_.push_back(&undoElem_);
    controlElements_.push_back(&resetElem_);
    //controlElements_.push_back(&exportVideoElem_);

    //add Properties
    addProperty(fontProperty_);
    addProperty(timelineDuration_);
    addProperty(timelineStepSize_);
    addProperty(fps_);
    addProperty(windFactor_);
    //add animated Properties
    addProperty(transfuncProp_); transfuncProp_.setVisibleFlag(false);
    addProperty(cameraProp_); cameraProp_.setVisibleFlag(false);
    addProperty(toggleClipping0Prop_); toggleClipping0Prop_.setVisibleFlag(false);
    addProperty(toggleClipping1Prop_); toggleClipping1Prop_.setVisibleFlag(false);
    addProperty(toggleClipping2Prop_); toggleClipping2Prop_.setVisibleFlag(false);
    addProperty(planeNormal0Prop_); planeNormal0Prop_.setVisibleFlag(false);
    addProperty(planeNormal1Prop_); planeNormal1Prop_.setVisibleFlag(false);
    addProperty(planeNormal2Prop_); planeNormal2Prop_.setVisibleFlag(false);
    addProperty(planePosition0Prop_); planePosition0Prop_.setVisibleFlag(false);
    addProperty(planePosition1Prop_); planePosition1Prop_.setVisibleFlag(false);
    addProperty(planePosition2Prop_); planePosition2Prop_.setVisibleFlag(false);

    //on change
    timelineDuration_.onChange(MemberFunctionCallback<TouchTableAnimationWidget>(this, &TouchTableAnimationWidget::updateTimeline));
    timelineStepSize_.onChange(MemberFunctionCallback<TouchTableAnimationWidget>(this, &TouchTableAnimationWidget::updateTimeline));
}

TouchTableAnimationWidget::~TouchTableAnimationWidget(){
    delete playTimer_;
    delete fastforwardTimer_;
    delete rewindTimer_;
    delete animation_;
}

void TouchTableAnimationWidget::initialize() {
    TouchTableMenuWidget::initialize();

    //load textures for control elements
    resetTex_ = TexMgr.load("newAnimation.png");
    //undoTex_ = TexMgr.load("undo.png");
    //redoTex_ = TexMgr.load("redo.png");
    jumpEndTex_ = TexMgr.load("jumpEnd.png");
    jumpBeginningTex_ = TexMgr.load("jumpBeginning.png");
    rewindTex_ = TexMgr.load("rewind.png");
    fastforwardTex_ = TexMgr.load("fastforward.png");
    pauseTex_ = TexMgr.load("pause.png");
    playTex_ = TexMgr.load("play.png");
    recordTex_ = TexMgr.load("record.png");
    cameraTex_= TexMgr.load("camera.png");
    clippingTex_ = TexMgr.load("Scalpel-icon.png");
    transfuncTex_ = TexMgr.load("tfwidget.png");
    deleteTex_ = TexMgr.load("delete.png");
    //exportVideoTex_ = TexMgr.load("exportVideo.png");
    confirmTex_ = TexMgr.load("apply.png");

    //set textures
    resetElem_.setSymbolTexture(resetTex_);
    //undoElem_.setSymbolTexture(undoTex_);
    //redoElem_.setSymbolTexture(redoTex_);
    jumpEndElem_.setSymbolTexture(jumpEndTex_);
    jumpBeginningElem_.setSymbolTexture(jumpBeginningTex_);
    rewindElem_.setSymbolTexture(rewindTex_);
    fastforwardElem_.setSymbolTexture(fastforwardTex_);
    pauseElem_.setSymbolTexture(pauseTex_);
    playElem_.setSymbolTexture(playTex_);
    recordElem_.setSymbolTexture(recordTex_);
    cameraElem_.setSymbolTexture(cameraTex_);
    clippingElem_.setSymbolTexture(clippingTex_);
    transfuncElem_.setSymbolTexture(transfuncTex_);
    deleteElem_.setSymbolTexture(deleteTex_);
    //exportVideoElem_.setSymbolTexture(exportVideoTex_);
    confirmElem_.setSymbolTexture(confirmTex_);

    //set font attributes
    fontProperty_.get()->setTextAlignment(tgt::Font::Center);
    fontProperty_.get()->setVerticalTextAlignment(tgt::Font::Middle);
    fontProperty_.get()->setFontSize(25);
    //fontProperty_.get()->setLineWidth(static_cast<float>(timeMenu_.getUR().x - timeMenu_.getLL().x - 20));
    fontProperty_.get()->setLineWidth(30.f);


    /*OLD SOLUTION; CREATING PROPERTYTIMELINES VIA LINKED PROPERTIES
    //initialize property timelines
    transfuncTimeline_ = PropertyTimelineFactory::getInstance()->createTimeline(&transfuncProp_);
    cameraTimeline_ = PropertyTimelineFactory::getInstance()->createTimeline(&cameraProp_);
    clippingTimelines_.push_back(PropertyTimelineFactory::getInstance()->createTimeline(&toggleClipping0Prop_));
    clippingTimelines_.push_back(PropertyTimelineFactory::getInstance()->createTimeline(&toggleClipping1Prop_));
    clippingTimelines_.push_back(PropertyTimelineFactory::getInstance()->createTimeline(&toggleClipping2Prop_));
    clippingTimelines_.push_back(PropertyTimelineFactory::getInstance()->createTimeline(&planeNormal0Prop_));
    clippingTimelines_.push_back(PropertyTimelineFactory::getInstance()->createTimeline(&planeNormal1Prop_));
    clippingTimelines_.push_back(PropertyTimelineFactory::getInstance()->createTimeline(&planeNormal2Prop_));
    clippingTimelines_.push_back(PropertyTimelineFactory::getInstance()->createTimeline(&planePosition0Prop_));
    clippingTimelines_.push_back(PropertyTimelineFactory::getInstance()->createTimeline(&planePosition1Prop_));
    clippingTimelines_.push_back(PropertyTimelineFactory::getInstance()->createTimeline(&planePosition2Prop_));
    //add all timelines in one vector for common operations
    propertyTimelines_.push_back(transfuncTimeline_);
    propertyTimelines_.push_back(cameraTimeline_);
    for(std::vector<PropertyTimeline*>::iterator timelineIter = clippingTimelines_.begin(); timelineIter != clippingTimelines_.end(); ++timelineIter){
        propertyTimelines_.push_back(*timelineIter);
    }*/


    //default settings
    fontProperty_.get()->setFontSize(21);
    fontProperty_.get()->setLineWidth(80);
    fontProperty_.get()->setVerticalTextAlignment(tgt::Font::Middle);
    fontProperty_.get()->setTextAlignment(tgt::Font::Left);
    timelineStepSize_.set(60);
    menuDimensions_.set(tgt::ivec2(924,236));
    cameraElem_.setActive(true);
    transfuncElem_.setActive(true);
    clippingElem_.setActive(true);
    timelineDuration_.set(900);

    //set selection handler for timelineMenu
    timelineMenu_.setSelectionHandler(std::bind1st(std::mem_fun(&TouchTableAnimationWidget::handleSelection), this));

}

void TouchTableAnimationWidget::deinitialize() {

    playTimer_->stop();
    fastforwardTimer_->stop();
    rewindTimer_->stop();

    //clean up

    //dispose Textures
    TexMgr.dispose(resetTex_); resetTex_ = 0;
    //TexMgr.dispose(undoTex_); undoTex_ = 0;
    //TexMgr.dispose(redoTex_); reduTex = 0;
    TexMgr.dispose(confirmTex_); confirmTex_ = 0;
    TexMgr.dispose(jumpEndTex_); jumpEndTex_ = 0;
    TexMgr.dispose(jumpBeginningTex_); jumpBeginningTex_ = 0;
    TexMgr.dispose(rewindTex_); rewindTex_ = 0;
    TexMgr.dispose(fastforwardTex_); fastforwardTex_ = 0;
    TexMgr.dispose(pauseTex_); pauseTex_ = 0;
    TexMgr.dispose(playTex_); playTex_ = 0;
    TexMgr.dispose(recordTex_); recordTex_ = 0;
    TexMgr.dispose(cameraTex_); cameraTex_ = 0;
    TexMgr.dispose(clippingTex_); clippingTex_ = 0;
    TexMgr.dispose(transfuncTex_); transfuncTex_ = 0;
    TexMgr.dispose(deleteTex_); deleteTex_= 0;

    TouchTableMenuWidget::deinitialize();
}

void TouchTableAnimationWidget::handleTouchPoints(const std::deque<tgt::TouchPoint>& tp) {

    std::deque<tgt::TouchPoint> pointsForTimelineMenu;

    for (std::deque<tgt::TouchPoint>::const_iterator currentTp = tp.begin(); currentTp != tp.end(); ++currentTp) {

        tgt::vec2 tpPos = convertToMenuCoordinates(currentTp->pos());

        //touchpoints in timelineMenu are handeld by the same
        if(timelineMenu_.contains(tpPos) && !newAnimationMode_){
            pointsForTimelineMenu.push_back(*currentTp);
            continue;
        }

        if (currentTp->state() == tgt::TouchPoint::TouchPointPressed) {
            //if control element has been hit add to mapping with according touchpoint
            for(std::vector<TouchTableControlElement*>::iterator controlElementIter = controlElements_.begin(); controlElementIter != controlElements_.end(); ++controlElementIter){
                TouchTableControlElement* currentElem = *controlElementIter;
                if(hitsControlElement(tpPos, *currentElem) && !newAnimationMode_ ){
                    tpControlElementMapping_.insert(std::make_pair(currentTp->id(), currentElem));
                    continue;
                }
            }
            //if deleteElem hit
            if(hitsControlElement(tpPos, deleteElem_)&& !newAnimationMode_){
                tpControlElementMapping_.insert(std::make_pair(currentTp->id(), &deleteElem_));
                continue;
            }
            // hits confirm element
            if(hitsControlElement(tpPos, confirmElem_) && newAnimationMode_){
                tpControlElementMapping_.insert(std::make_pair(currentTp->id(), &confirmElem_));
                continue;
            }
            // hits slider
            if(durationSlider_.checkIndicatorHit(tpPos)&& newAnimationMode_){
                durationSliderIds_.push_back(currentTp->id());
                continue;
            }


        }

        if (currentTp->state() == tgt::TouchPoint::TouchPointMoved) {
            std::vector<int>::iterator i = std::find(durationSliderIds_.begin(), durationSliderIds_.end(), currentTp->id());

            if (i != durationSliderIds_.end() && newAnimationMode_) {
                durationSlider_.updateIndicatorPosition(tpPos);
                timelineDuration_.set(static_cast<int>(durationSlider_.getIndicatorPos() * (timelineDuration_.getMaxValue() - timelineDuration_.getMinValue()) + timelineDuration_.getMinValue()));
                invalidate();
                continue;
            }

        }

        if(currentTp->state() == tgt::TouchPoint::TouchPointReleased){
            //check whether the slider has been released and erase tp from mapping
            std::vector<int>::iterator i = std::find(durationSliderIds_.begin(), durationSliderIds_.end(), currentTp->id());
            if(i != durationSliderIds_.end() && newAnimationMode_){
                durationSliderIds_.erase(i);
                continue;
            }

            //check if touchpoint in mapping, else do nothing and jump to next jp in deque
            std::map<int,TouchTableControlElement*>::iterator mappingIter = tpControlElementMapping_.find(currentTp->id());
            if( mappingIter== tpControlElementMapping_.end() && !newAnimationMode_)
                continue;

            //play element hit
            if(hitsControlElement(tpPos, playElem_)&& !newAnimationMode_){
                fastforwardTimer_->stop(); rewindTimer_->stop();
                playTimer_->start(static_cast<int>(1000.0 / animation_->getFPS()), 0);
                tpControlElementMapping_.erase(mappingIter);
                continue;
            }

            //fastforward element hit
            if(hitsControlElement(tpPos, fastforwardElem_)&& !newAnimationMode_){
                playTimer_->stop(); rewindTimer_->stop();
                fastforwardTimer_->start(static_cast<int>(1000.0 / animation_->getFPS()), 0);
                tpControlElementMapping_.erase(mappingIter);
                continue;
            }

            //rewind element hit
            if(hitsControlElement(tpPos, rewindElem_)&& !newAnimationMode_){
                playTimer_->stop(); fastforwardTimer_->stop();
                rewindTimer_->start(static_cast<int>(1000.0 / animation_->getFPS()), 0);
                tpControlElementMapping_.erase(mappingIter);
                continue;
            }

            //pause element hit
            if(hitsControlElement(tpPos, pauseElem_)&& !newAnimationMode_){
                playTimer_->stop(); fastforwardTimer_->stop(); rewindTimer_->stop();
                tpControlElementMapping_.erase(mappingIter);
                continue;
            }

            //jump to Beginning element hit
            if(hitsControlElement(tpPos, jumpBeginningElem_)&& !newAnimationMode_){
                //set gui
                timelineMenu_.getTimeline().setCurrentPos(1- timelineMenu_.getTimeline().getTimeOffset());
                timelineMenu_.setIndicatorPosition(0);
                animation_->renderAt(0.0);
                tpControlElementMapping_.erase(mappingIter);
                invalidate();
                continue;
            }

            //jump to End element hit
            if(hitsControlElement(tpPos, jumpEndElem_)&& !newAnimationMode_){
                //set gui
                timelineMenu_.getTimeline().setCurrentPos(timelineDuration_.get()- timelineMenu_.getTimeline().getTimeOffset()-1);
                timelineMenu_.setIndicatorPosition(1);
                animation_->renderAt(timelineMenu_.getTimeline().getTimeForPos(timelineDuration_.get()));
                tpControlElementMapping_.erase(mappingIter);
                invalidate();
                continue;
            }

            //record element hit
            if(hitsControlElement(tpPos, recordElem_)&& !newAnimationMode_){
                //add new keyValue to PropertyTimelines
                /************GUI****************************/
                bool keyValueAdded = timelineMenu_.addKeyValue(); //add keyValue to the timeline GUI, if no keyvalue already exists at current position
                /**************internal logic***************/
                if(keyValueAdded){
                    float currentTime = timelineMenu_.getTimeline().getCurrentTime();
                    transfuncTimeline_->setCurrentSettingAsKeyvalue(currentTime, true);
                    cameraTimeline_->setCurrentSettingAsKeyvalue(currentTime, false);
                    for(std::vector<PropertyTimelineBool*>::iterator i = clippingToggleTimelines_.begin(); i != clippingToggleTimelines_.end(); ++i){
                        (*i)->setCurrentSettingAsKeyvalue(currentTime,false);
                    }
                    for(std::vector<PropertyTimelineFloat*>::iterator i = clippingPlanePositionTimelines_.begin(); i != clippingPlanePositionTimelines_.end(); ++i){
                        (*i)->setCurrentSettingAsKeyvalue(currentTime, false);
                    }
                    clippingPlanePositionTimelines_.clear();
                    for(std::vector<PropertyTimelineVec3*>::iterator i = clippingPlaneNormalTimelines_.begin(); i != clippingPlaneNormalTimelines_.end(); ++i){
                        (*i)->setCurrentSettingAsKeyvalue(currentTime, false);
                    }
                    //animation_->setActualNetworkAsKeyvalues(timelineMenu_.getTimeline().getCurrentTime());
                }

                tpControlElementMapping_.erase(mappingIter);
                invalidate();
                continue;
            }

            //delete elemet hit
            if(hitsControlElement(tpPos, deleteElem_)&& !newAnimationMode_){
                /**************internal logic***************/
                KeyValue* keyValue = timelineMenu_.getTimeline().getSelectedKeyValue();
                if(keyValue)
                    deleteKeyValue(keyValue->time_);
                /************GUI****************************/
                timelineMenu_.getTimeline().eraseSelectedKeyValue();

                tpControlElementMapping_.erase(mappingIter);
                invalidate();
                continue;
            }

            //camera element hit
            if(hitsControlElement(tpPos, cameraElem_)&& !newAnimationMode_){
                cameraTimeline_->setActiveOnRendering(!cameraElem_.isActive());
                cameraElem_.setActive(!cameraElem_.isActive());
                tpControlElementMapping_.erase(mappingIter);
                invalidate();
                continue;
            }

            //transfunc element hit
            if(hitsControlElement(tpPos, transfuncElem_)&& !newAnimationMode_){
                transfuncTimeline_->setActiveOnRendering(!transfuncElem_.isActive());
                transfuncElem_.setActive(!transfuncElem_.isActive());
                tpControlElementMapping_.erase(mappingIter);
                invalidate();
                continue;
            }

            //clipping element hit
            if(hitsControlElement(tpPos, clippingElem_)&& !newAnimationMode_){
                for(std::vector<PropertyTimelineBool*>::iterator i = clippingToggleTimelines_.begin(); i!= clippingToggleTimelines_.end(); ++i){
                    (*i)->setActiveOnRendering(!clippingElem_.isActive());
                }for(std::vector<PropertyTimelineFloat*>::iterator i = clippingPlanePositionTimelines_.begin(); i!= clippingPlanePositionTimelines_.end(); ++i){
                    (*i)->setActiveOnRendering(!clippingElem_.isActive());
                }for(std::vector<PropertyTimelineVec3*>::iterator i = clippingPlaneNormalTimelines_.begin(); i!= clippingPlaneNormalTimelines_.end(); ++i){
                    (*i)->setActiveOnRendering(!clippingElem_.isActive());
                }
                clippingElem_.setActive(!clippingElem_.isActive());
                tpControlElementMapping_.erase(mappingIter);
                invalidate();
                continue;
            }


            //reset element hit
            if(hitsControlElement(tpPos, resetElem_)&& !newAnimationMode_){
                setUpNewAnimationMode();
                invalidate();
                continue;
            }

            //confirm element hit
            if(hitsControlElement(tpPos, confirmElem_) && newAnimationMode_){
                newAnimationMode_ = false;
                timelineDuration_.set(static_cast<int>(durationSlider_.getIndicatorPos() * timelineDuration_.getMaxValue()));
                //TODO set up the animation here
                setUpNewAnimation();
                tpControlElementMapping_.erase(mappingIter);
                invalidate();
                continue;
            }

        }

    }//end for

    if(!pointsForTimelineMenu.empty()){
        timelineMenu_.handleTouchPoints(convertTouchPointsToMenuCoordinates(pointsForTimelineMenu));
        invalidate();
    }
}

void TouchTableAnimationWidget::setUpNewAnimation(){
    //****************************Setting up Animation******************************************************************************************/
    //getting the animatedProcessors from network
    //WARNING: cast from "const ProcessorNetwork*" to "ProcessorNetwork"
    networkEval_ = VoreenApplication::app()->getNetworkEvaluator();
    network_  = (ProcessorNetwork*) networkEval_->getProcessorNetwork() ;
    animation_ =new Animation(network_);
    animation_->setNetwork(network_);
    animation_->setDuration(timelineMenu_.getTimeline().getTimeForPos(timelineDuration_.get()));
    animation_->setFPS(fps_.get());
    animation_->setUndoSteps(1);
    animatedProcessors_ = animation_->getAnimatedProcessors();
    //add timelines to those processors that are being animated.
    for(std::vector<AnimatedProcessor*>::iterator animatedProcIter = animatedProcessors_.begin(); animatedProcIter != animatedProcessors_.end(); ++animatedProcIter){
        AnimatedProcessor* currentAnimatedProc = *animatedProcIter;
        std::string currentName = currentAnimatedProc->getProcessorName();

        if(currentName == "TouchTableClippingWidget" ){
            Processor* proc = currentAnimatedProc->getCorrespondingProcessor();
            std::vector<Property*> props = proc->getProperties();
            for(std::vector<Property*>::iterator propIter = props.begin(); propIter != props.end(); ++propIter){
                Property* currentProp = *propIter;
                std::string propName = currentProp->getID();
                if(propName == "clipplane0toggle" || propName == "clipplane1toggle" || propName == "clipplane2toggle"){
                    currentAnimatedProc->addTimeline(currentProp);

                }
                if(propName == "clipplane0planeNormal" || propName == "clipplane1planeNormal" || propName == "clipplane2planeNormal"){
                    currentAnimatedProc->addTimeline(currentProp);
                }
                if(propName == "clipplane0planeposition" || propName == "clipplane1planeposition" || propName == "clipplane2planeposition"){
                    currentAnimatedProc->addTimeline(currentProp);
                }

            }
        }
        if(currentName == "TouchTableTransFuncWidget"){
            Processor* proc = currentAnimatedProc->getCorrespondingProcessor();
            std::vector<Property*> props = proc->getProperties();
            for(std::vector<Property*>::iterator propIter = props.begin(); propIter != props.end(); ++propIter){
                Property* currentProp = *propIter;
                std::string propName = currentProp->getID();
                if(propName == "transfunc"){
                    currentAnimatedProc->addTimeline(currentProp);
                    transfuncTimeline_ = (PropertyTimelineTransFunc1DKeys*) currentAnimatedProc->getPropertyTimelines().back();
                    transfuncTimeline_->registerUndoObserver(animation_);
                    //transfuncTimeline_->addObserver(animation_);
                }
            }
        }
        if(currentName == "SingleVolumeRaycaster"){
            Processor* proc = currentAnimatedProc->getCorrespondingProcessor();
            std::vector<Property*> props = proc->getProperties();
            for(std::vector<Property*>::iterator propIter = props.begin(); propIter != props.end(); ++propIter){
                Property* currentProp = *propIter;
                std::string propName = currentProp->getID();
                if(propName == "camera"){
                    currentAnimatedProc->addTimeline(currentProp);
                    cameraTimeline_  = (PropertyTimelineCamera*) currentAnimatedProc->getPropertyTimelines().back();
                    cameraTimeline_->registerUndoObserver(animation_);
                    //cameraTimeline_->addObserver(animation_);
                }
            }
        }
    }
    for(std::vector<AnimatedProcessor*>::iterator i = animatedProcessors_.begin(); i!=animatedProcessors_.end(); ++i){
        std::vector<PropertyTimeline*> temp = (*i)->getPropertyTimelines();
        for(std::vector<PropertyTimeline*>::iterator it = temp.begin(); it != temp.end(); ++ it){
            std::string propName = (*it)->getProperty()->getID();
            if(propName == "clipplane0toggle" || propName == "clipplane1toggle" || propName == "clipplane2toggle"){
                clippingToggleTimelines_.push_back((PropertyTimelineBool*) *it);
                (*it)->registerUndoObserver(animation_);
            }

            if(propName == "clipplane0planeNormal" || propName == "clipplane1planeNormal" || propName == "clipplane2planeNormal"){
                clippingPlaneNormalTimelines_.push_back((PropertyTimelineVec3*) *it);
                (*it)->registerUndoObserver(animation_);
            }

            if(propName == "clipplane0planeposition" || propName == "clipplane1planeposition" || propName == "clipplane2planeposition"){
                clippingPlanePositionTimelines_.push_back((PropertyTimelineFloat*) *it);
                (*it)->registerUndoObserver(animation_);
            }
        }
    }
    /*************************Setting up Animation****************************************************************/

}

void TouchTableAnimationWidget::setUpNewAnimationMode(){
    newAnimationMode_ = true;
    //set up slider
    durationSlider_.setIndicatorWidth(700);
    durationSlider_.setIndicatorPos(0.5);
    timelineDuration_.set(timelineDuration_.getMaxValue() / 2);

    //delete all Animation Settings, delete all existing timelines
    cameraTimeline_ = 0;
    transfuncTimeline_ = 0 ;
    for(std::vector<PropertyTimelineBool*>::iterator i = clippingToggleTimelines_.begin(); i != clippingToggleTimelines_.end(); ++i){
        (*i) = 0;
    }
    clippingToggleTimelines_.clear();
    for(std::vector<PropertyTimelineFloat*>::iterator i = clippingPlanePositionTimelines_.begin(); i != clippingPlanePositionTimelines_.end(); ++i){
        (*i) = 0;
    }
    clippingPlanePositionTimelines_.clear();
    for(std::vector<PropertyTimelineVec3*>::iterator i = clippingPlaneNormalTimelines_.begin(); i != clippingPlaneNormalTimelines_.end(); ++i){
        (*i) = 0;
    }
    clippingPlaneNormalTimelines_.clear();
    for(std::vector<AnimatedProcessor*>::iterator procIter = animatedProcessors_.begin(); procIter != animatedProcessors_.end(); ++procIter){
        std::vector<PropertyTimeline*> propTimelines = (*procIter)->getPropertyTimelines();
        for(std::vector<PropertyTimeline*>::iterator timelineIter = propTimelines.begin(); timelineIter != propTimelines.end(); ++timelineIter){
            (*procIter)->removeTimeline((*timelineIter)->getProperty());
        }
    }
    animatedProcessors_.clear();
    if(animation_)
        delete animation_;
    animation_ = 0;

    //reset GUI
    timelineMenu_.getTimeline().getKeyValues().clear();

}

void TouchTableAnimationWidget::deleteKeyValue(float time){
    std::map<float, PropertyKeyValue<TransFunc1DKeys*>*> const transfuncMap = transfuncTimeline_->getKeyValues();
    for(std::map<float, PropertyKeyValue<TransFunc1DKeys*>*>::const_iterator i= transfuncMap.begin(); i != transfuncMap.end(); ++i){
        if(std::abs(i->first - time) <0.1)
            transfuncTimeline_->deleteKeyValue(i->second);
    }

    std::map<float, PropertyKeyValue<tgt::Camera>*> const cameraMap = cameraTimeline_->getKeyValues();
    for(std::map<float, PropertyKeyValue<tgt::Camera>*>::const_iterator i= cameraMap.begin(); i != cameraMap.end(); ++i){
        if(std::abs(i->first - time) <0.1)
            cameraTimeline_->deleteKeyValue(i->second);
    }

    for(std::vector<PropertyTimelineBool*>::iterator i = clippingToggleTimelines_.begin(); i!= clippingToggleTimelines_.end(); ++i){
        std::map<float, PropertyKeyValue<bool>*> const toggleMap = (*i)->getKeyValues();
        for(std::map<float, PropertyKeyValue<bool>*>::const_iterator it = toggleMap.begin(); it!= toggleMap.end(); ++it){
            if(std::abs(it->first - time) < 0.1)
                (*i)->deleteKeyValue(it->second);
        }
    }
    for(std::vector<PropertyTimelineVec3*>::iterator i = clippingPlaneNormalTimelines_.begin(); i!= clippingPlaneNormalTimelines_.end(); ++i){
        std::map<float, PropertyKeyValue<tgt::vec3>*> const planeNormalMap = (*i)->getKeyValues();
        for(std::map<float, PropertyKeyValue<tgt::vec3>*>::const_iterator it = planeNormalMap.begin(); it!= planeNormalMap.end(); ++it){
            if(std::abs(it->first - time) < 0.1)
                (*i)->deleteKeyValue(it->second);
        }
    }
    for(std::vector<PropertyTimelineFloat*>::iterator i = clippingPlanePositionTimelines_.begin(); i!= clippingPlanePositionTimelines_.end(); ++i){
        std::map<float, PropertyKeyValue<float>*> const planePositionMap = (*i)->getKeyValues();
        for(std::map<float, PropertyKeyValue<float>*>::const_iterator it = planePositionMap.begin(); it!= planePositionMap.end(); ++it){
            if(std::abs(it->first - time) < 0.1)
                (*i)->deleteKeyValue(it->second);
        }
    }
}

void TouchTableAnimationWidget::handleSelection(boost::tuple<float,bool,float> t){
    float newTime = boost::get<0>(t);
    float oldTime = boost::get<2>(t);
    bool keyValueHit = boost::get<1>(t);
    //if key Value has been hit, change position of all propertyKeyValues according to the hit keyValue in the gui
    if(keyValueHit){
        std::map<float, PropertyKeyValue<TransFunc1DKeys*>*> const transfuncMap = transfuncTimeline_->getKeyValues();
        for(std::map<float, PropertyKeyValue<TransFunc1DKeys*>*>::const_iterator i= transfuncMap.begin(); i != transfuncMap.end(); ++i){
            if(std::abs(i->first - oldTime) < 0.1)
                transfuncTimeline_->changeTimeOfKeyValue(newTime, i->second);
        }

        std::map<float, PropertyKeyValue<tgt::Camera>*> const cameraMap = cameraTimeline_->getKeyValues();
        for(std::map<float, PropertyKeyValue<tgt::Camera>*>::const_iterator i= cameraMap.begin(); i != cameraMap.end(); ++i){
            if(std::abs(i->first - oldTime) <0.1)
            cameraTimeline_->changeTimeOfKeyValue(newTime, i->second);
        }

        for(std::vector<PropertyTimelineBool*>::iterator i = clippingToggleTimelines_.begin(); i!= clippingToggleTimelines_.end(); ++i){
            std::map<float, PropertyKeyValue<bool>*> const toggleMap = (*i)->getKeyValues();
            for(std::map<float, PropertyKeyValue<bool>*>::const_iterator it = toggleMap.begin(); it!= toggleMap.end(); ++it){
                if(std::abs(it->first - oldTime) < 0.1)
                    (*i)->changeTimeOfKeyValue(newTime, it->second);
            }
        }
        for(std::vector<PropertyTimelineVec3*>::iterator i = clippingPlaneNormalTimelines_.begin(); i!= clippingPlaneNormalTimelines_.end(); ++i){
            std::map<float, PropertyKeyValue<tgt::vec3>*> const planeNormalMap = (*i)->getKeyValues();
            for(std::map<float, PropertyKeyValue<tgt::vec3>*>::const_iterator it = planeNormalMap.begin(); it!= planeNormalMap.end(); ++it){
                if(std::abs(it->first - oldTime) < 0.1)
                    (*i)->changeTimeOfKeyValue(newTime, it->second);
            }
        }
        for(std::vector<PropertyTimelineFloat*>::iterator i = clippingPlanePositionTimelines_.begin(); i!= clippingPlanePositionTimelines_.end(); ++i){
            std::map<float, PropertyKeyValue<float>*> const planePositionMap = (*i)->getKeyValues();
            for(std::map<float, PropertyKeyValue<float>*>::const_iterator it = planePositionMap.begin(); it!= planePositionMap.end(); ++it){
                if(std::abs(it->first - oldTime) < 0.1)
                    (*i)->changeTimeOfKeyValue(newTime, it->second);
            }
        }
    }else{
        animation_->renderAt(newTime);
    }
}

void TouchTableAnimationWidget::timerEvent(tgt::TimeEvent* te){
    tgt::Timer* currentTimer = te->getTimer();
    if(currentTimer == playTimer_){
        if(timelineMenu_.getTimeline().getCurrentPos() >= timelineDuration_.get()){//stop when reached the end
            timelineMenu_.getTimeline().setCurrentPos(timelineDuration_.get()-1-timelineMenu_.getTimeline().getTimeOffset());
            currentTimer->stop();
        }
        else{
            timelineMenu_.getTimeline().setCurrentTime(timelineMenu_.getTimeline().getCurrentTime() + (1.0f/animation_->getFPS()));
            animation_->renderAt(timelineMenu_.getTimeline().getCurrentTime());
        }


    }
    if(currentTimer == rewindTimer_){
        if(timelineMenu_.getTimeline().getCurrentPos() <= 1){//stop when reached the end
            timelineMenu_.getTimeline().setCurrentPos(1-timelineMenu_.getTimeline().getTimeOffset());
            currentTimer->stop();
        }else{
            timelineMenu_.getTimeline().setCurrentTime(timelineMenu_.getTimeline().getCurrentTime() - 1.0f/animation_->getFPS()*windFactor_.get());
            animation_->renderAt(timelineMenu_.getTimeline().getCurrentTime());
        }

    }
    if(currentTimer == fastforwardTimer_){
        if(timelineMenu_.getTimeline().getCurrentPos() >= timelineDuration_.get()-1){//stop when reached the end
            timelineMenu_.getTimeline().setCurrentPos(timelineDuration_.get()-1-timelineMenu_.getTimeline().getTimeOffset());
            currentTimer->stop();
        }else{
            timelineMenu_.getTimeline().setCurrentTime(timelineMenu_.getTimeline().getCurrentTime() + 1.0f/animation_->getFPS()*windFactor_.get());
            animation_->renderAt(timelineMenu_.getTimeline().getCurrentTime());
        }

    }
    //if cursor runs out of visible space of contentMenu, reset slider position
    updateTimelineSlider();
    invalidate();
}

void TouchTableAnimationWidget::renderScrollableTimeline(){
    tgt::ivec2 offset = menu_.getLL();
    overlay_->renderMenuFrame(timelineMenu_, true, offset);
    overlay_->renderMenuFrame(timelineMenu_.getContentMenu(), true, offset);
    overlay_->renderSlider(&timelineMenu_.getSlider(), offset);

    //render timeline
    // set OpenGL status: depth func and blending
    glDepthFunc(GL_ALWAYS);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    //set transformation to use pixel coordinates
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadMatrix(tgt::mat4::createOrtho(0, overlay_->getOutportSize().x, 0, overlay_->getOutportSize().y, -1, 1));
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();


        //if menu of the mode widget is inverted and this should be taken into account: set transformation
        if (isMenuInverted()) {
            float xt = - (static_cast<float>(menu_.getLL().x) + static_cast<float>(menuDimensions_.get().x) / 2.f);
            float yt = - (static_cast<float>(menu_.getLL().y) + static_cast<float>(menuDimensions_.get().y) / 2.f);

            MatStack.translate(-xt,-yt, 0.f);
            MatStack.rotate(180.f, 0.f, 0.f, 1.f);
            MatStack.translate(xt,yt, 0.f);
        }

    glLineWidth(4.f);
    glColor3f(0.f,0.f,0.f);
     TouchTableMenuFrame& contentMenu = timelineMenu_.getContentMenu();
    int yPos = offset.y + contentMenu.getLL().y + (contentMenu.getUR().y - contentMenu.getLL().y) /2;
    //render horizontal line
    glBegin(GL_LINES);
        glVertex2f(static_cast<float>(offset.x+ contentMenu.getLL().x),static_cast<float>(yPos));
        glVertex2f(static_cast<float>(offset.x + contentMenu.getUR().x), static_cast<float>(yPos));
    glEnd();

    //render vertical lines
    glLineWidth(2.f);
    int stepSize = timelineMenu_.getTimeline().getStepSize();
    int duration = timelineMenu_.getTimeline().getDuration();
    int timeOffset = timelineMenu_.getTimeline().getTimeOffset();
    for(int xPos= stepSize - (timeOffset % stepSize) ; xPos<= duration && xPos< (contentMenu.getUR().x - contentMenu.getLL().x); xPos+= stepSize){
        glBegin(GL_LINES);
        glVertex2f(static_cast<float>(contentMenu.getLL().x + offset.x+ xPos),static_cast<float>(yPos +10));
        glVertex2f(static_cast<float>(contentMenu.getLL().x + offset.x + xPos), static_cast<float>(yPos-10));
        glEnd();
    }

    //render current position if in visible space
    glLineWidth(3.f);
    glColor3f(0.f,1.f,0.f);
    int currentPos = timelineMenu_.getTimeline().getCurrentPos();
    if(currentPos > timeOffset && currentPos < timeOffset + timelineMenu_.getContentMenu().getUR().x - timelineMenu_.getContentMenu().getLL().x){
        glBegin(GL_LINES);
        glVertex2f(static_cast<float>(contentMenu.getLL().x + offset.x + currentPos - timeOffset), static_cast<float>(offset.y + contentMenu.getUR().y));
        glVertex2f(static_cast<float>(contentMenu.getLL().x + offset.x + currentPos - timeOffset), static_cast<float>(offset.y + contentMenu.getLL().y));
        glEnd();
    }

    //render keyValues
    std::vector<KeyValue> keyValues = timelineMenu_.getTimeline().getKeyValues();
    for(std::vector<KeyValue>::iterator keyValueIter = keyValues.begin(); keyValueIter != keyValues.end(); ++keyValueIter){
        tgt::vec2 posLL( static_cast<float>(contentMenu.getLL().x + offset.x + keyValueIter->pos_ - keyValueIter->width_/2 - timeOffset), static_cast<float>(yPos - keyValueIter->height_/2));
        tgt::vec2 posUR( static_cast<float>(contentMenu.getLL().x + offset.x + keyValueIter->pos_ + keyValueIter->width_/2 - timeOffset), static_cast<float>(yPos + keyValueIter->height_/2));
        tgt::vec2 llBound = contentMenu.getLL() + offset;
        tgt::vec2 urBound = contentMenu.getUR() + offset;
        //if(posUR.x <= contentMenu.getUR().x + offset.x && posLL.x >= contentMenu.getLL().x + offset.x)
            if(keyValueIter->isSelected_)
                overlay_->renderColoredQuad(tgt::clamp(posLL,llBound,urBound),  tgt::clamp(posUR, llBound, urBound),tgt::vec4(0.1f,0.1f,0.3f,1.f));
            else
                overlay_->renderColoredQuad(tgt::clamp(posLL,llBound,urBound),  tgt::clamp(posUR, llBound, urBound), tgt::vec4(0.1f,0.1f,0.3f,0.7));
    }

    // set back OpenGL status
    glLineWidth(1.f);
    glColor3f(1.f,1.f,1.f);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glDisable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ZERO);

    //set transformation to default
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();

}

void TouchTableAnimationWidget::renderTime(){
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadMatrix(tgt::mat4::createOrtho(0, overlay_->getOutportSize().x, 0, overlay_->getOutportSize().y, -1, 1));
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();


        //if menu of the mode widget is inverted and this should be taken into account: set transformation
        if (isMenuInverted()) {
            float xt = - (static_cast<float>(getMenuLL().x) + static_cast<float>(menuDimensions_.get().x) / 2.f);
            float yt = - (static_cast<float>(getMenuLL().y) + static_cast<float>(menuDimensions_.get().y) / 2.f);

            MatStack.translate(-xt,-yt, 0.f);
            MatStack.rotate(180.f, 0.f, 0.f, 1.f);
            MatStack.translate(xt,yt, 0.f);
        }

    glDisable(GL_DEPTH_TEST);
    glColor4f(0.f,0.f,0.f,1.f);
    //fontProperty_.get()->setLineWidth(static_cast<float>(timeMenu_.getUR().x - timeMenu_.getLL().x - 20));

    //render time in timeMenu
    float yPos = static_cast<float>(timeMenu_.getUR().y - (timeMenu_.getUR().y- timeMenu_.getLL().y - fontProperty_.get()->getLineHeight()*2 + 5)/2);
    tgt::vec3 position(static_cast<float>(timeMenu_.getLL().x + 10), yPos , 0.f);
    fontProperty_.get()->render(position, overlay_->getOutportSize(), "Time:");
    position = tgt::vec3(static_cast<float>(timeMenu_.getLL().x + 10), yPos - 30, 0.f);
    fontProperty_.get()->render(position, overlay_->getOutportSize(), timelineMenu_.getTimeline().getCurrentTimeAsString());

    //render time in timeline
    int oldFontSize = fontProperty_.get()->getFontSize();
    fontProperty_.get()->setFontSize(12);
    fontProperty_.get()->setVerticalTextAlignment(tgt::Font::Top);
    int stepSize = timelineMenu_.getTimeline().getStepSize();
    int duration = timelineMenu_.getTimeline().getDuration();
    int timeOffset = timelineMenu_.getTimeline().getTimeOffset();
    int textSize = (int) fontProperty_.get()->getSize(tgt::vec3(0.f,0.f,0.f), overlay_->getOutportSize(), "00:00").x;
    for(int xPos= stepSize - (timeOffset % stepSize) - textSize/2 ; xPos<= duration && xPos< (timelineMenu_.getContentMenu().getUR().x - timelineMenu_.getContentMenu().getLL().x) - textSize; xPos+= stepSize){
        if(xPos < 0 ) continue;
        yPos = (float) menu_.getLL().y + timelineMenu_.getContentMenu().getLL().y + (timelineMenu_.getContentMenu().getUR().y - timelineMenu_.getContentMenu().getLL().y) /2;
        tgt::vec3 pos(static_cast<float>(timelineMenu_.getContentMenu().getLL().x + menu_.getLL().x + xPos),
                        static_cast<float>(yPos + 11), 0.f);
        std::string time = timelineMenu_.getTimeline().getTimeForPosAsString(xPos+textSize);

        fontProperty_.get()->render(pos, overlay_->getOutportSize(), time);
    }

    fontProperty_.get()->setFontSize(oldFontSize);
    fontProperty_.get()->setVerticalTextAlignment(tgt::Font::Middle);


    glEnable(GL_DEPTH_TEST);
    glColor4f(1.f,1.f,1.f,1.f);
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();

    LGL_ERROR;
}

void TouchTableAnimationWidget::renderNewAnimationText(){
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadMatrix(tgt::mat4::createOrtho(0, overlay_->getOutportSize().x, 0, overlay_->getOutportSize().y, -1, 1));
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();


        //if menu of the mode widget is inverted and this should be taken into account: set transformation
        if (isMenuInverted()) {
            float xt = - (static_cast<float>(getMenuLL().x) + static_cast<float>(menuDimensions_.get().x) / 2.f);
            float yt = - (static_cast<float>(getMenuLL().y) + static_cast<float>(menuDimensions_.get().y) / 2.f);

            MatStack.translate(-xt,-yt, 0.f);
            MatStack.rotate(180.f, 0.f, 0.f, 1.f);
            MatStack.translate(xt,yt, 0.f);
        }

    glDisable(GL_DEPTH_TEST);
    glColor4f(0.f,0.f,0.f,1.f);

    //save current settings of the fontProperty
    int lineWidth = static_cast<int>(fontProperty_.get()->getLineWidth());

    //render text
    fontProperty_.get()->setLineWidth(static_cast<float>(durationSlider_.getBarLength()));
    std::string text= "Laenge der Animation: ";
    tgt::vec3 textPos(static_cast<float>(menu_.getLL().x) + durationSlider_.getBarLL().x, static_cast<float>(menu_.getLL().y) + durationSlider_.getBarUR().y + 40, 0.f);

    fontProperty_.get()->render(textPos, overlay_->getOutportSize(), text);

    //render duration
    fontProperty_.get()->setLineWidth(200);
    int duration = timelineDuration_.get();
    std::string durationText = timelineMenu_.getTimeline().getTimeForPosAsString(duration) + " Minuten";
    tgt::vec3 durationPos(static_cast<float>(menu_.getLL().x) + durationSlider_.getBarLength()/2 - fontProperty_.get()->getLineWidth()/2 , static_cast<float>(menu_.getLL().y) + durationSlider_.getBarLL().y- fontProperty_.get()->getFontSize() - 10, 0.f);

    fontProperty_.get()->render(durationPos, overlay_->getOutportSize(), durationText);

    //restore old settings of fontproperty
    fontProperty_.get()->setLineWidth(static_cast<float>(lineWidth));

    glEnable(GL_DEPTH_TEST);
    glColor4f(1.f,1.f,1.f,1.f);
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();

    LGL_ERROR;
}

void TouchTableAnimationWidget::renderComponents() {

    if(newAnimationMode_){
        //render controlElement
        overlay_->renderControlElement(&confirmElem_, menu_.getLL());

        //render slider and text according to slider
        overlay_->renderSlider(&durationSlider_, menu_.getLL());
        renderNewAnimationText();

    }else{
        //render control elements
        for(std::vector<TouchTableControlElement*>::iterator controlElementIter = controlElements_.begin(); controlElementIter != controlElements_.end(); ++controlElementIter){
            overlay_->renderControlElement((*controlElementIter), menu_.getLL());
        }
        //if a keyValue is selected show delete control element
        if(timelineMenu_.getTimeline().anyKeyValueSelected()){
            overlay_->renderControlElement((&deleteElem_),menu_.getLL());
        }else{
            overlay_->renderControlElement(&deleteElem_, menu_.getLL(), 0.7f);
        }
        //render timeline
        renderScrollableTimeline();

        //render timeMenu and time
        overlay_->renderMenuFrame(timeMenu_);
        renderTime();
    }


}

void TouchTableAnimationWidget::updateComponents(){

    //update control element radius
    for(std::vector<TouchTableControlElement*>::iterator controlElementIter = controlElements_.begin(); controlElementIter != controlElements_.end(); ++controlElementIter){
        (*controlElementIter)->setRadius(controlElementRadius_.get());
    }

    //update positions of control elements
    int recordSetSize = 7;
    int smallSpace = 5;  //space between control elements in one set
    int yPos = menuDimensions_.get().y - smallSpace - controlElementRadius_.get();
    int controlElement = 2*controlElementRadius_.get(); //control element diameter
    //recordSet (play, pause ...)
    int recordsetStartPos = (menuDimensions_.get().x - ( 2*recordSetSize*controlElementRadius_.get() + 10*smallSpace)) / 2 - (5*smallSpace - 3*controlElementRadius_.get())/2;
    jumpBeginningElem_.setPosition(tgt::ivec2(recordsetStartPos + controlElementRadius_.get() + smallSpace, yPos));
    rewindElem_.setPosition(tgt::ivec2(recordsetStartPos + 3*controlElementRadius_.get() + 2*smallSpace,yPos));
    pauseElem_.setPosition(tgt::ivec2(recordsetStartPos + 5*controlElementRadius_.get() + 3*smallSpace,yPos));
    playElem_.setPosition(tgt::ivec2(recordsetStartPos + 7*controlElementRadius_.get() + 4*smallSpace,yPos));
    fastforwardElem_.setPosition(tgt::ivec2(recordsetStartPos + 9*controlElementRadius_.get() + 5*smallSpace,yPos));
    jumpEndElem_.setPosition(tgt::ivec2(recordsetStartPos + 11*controlElementRadius_.get() + 6*smallSpace,yPos));
    recordElem_.setPosition(tgt::ivec2(recordsetStartPos + 13*controlElementRadius_.get() + 10*smallSpace,yPos));
    //undo/redo Set
    //redoElem_.setPosition(tgt::ivec2(recordsetStartPos - 4*smallSpace - controlElementRadius_.get(), yPos));
    //undoElem_.setPosition(tgt::ivec2(recordsetStartPos - 5*smallSpace - 3*controlElementRadius_.get(), yPos));
    //resetAnimation
    resetElem_.setPosition(tgt::ivec2(smallSpace + controlElementRadius_.get(), yPos));
    //exportVideoElem_.setPosition(tgt::ivec2(2*smallSpace + 3*controlElementRadius_.get(), yPos));
    //widgetSet (camera, clipping ...)
    cameraElem_.setPosition(tgt::ivec2(menuDimensions_.get().x - smallSpace - controlElementRadius_.get(), yPos));
    clippingElem_.setPosition(tgt::ivec2(menuDimensions_.get().x - 2*smallSpace - 3*controlElementRadius_.get(), yPos));
    transfuncElem_.setPosition(tgt::ivec2(menuDimensions_.get().x - 3* smallSpace - 5* controlElementRadius_.get(), yPos));

    //update timeMenu
    timeMenu_.setLL(menu_.getLL() + tgt::ivec2(2*smallSpace, 2*controlElementRadius_.get() + 2*smallSpace));
    timeMenu_.setUR(menu_.getLL() + tgt::ivec2(10 + menuDimensions_.get().x / 10,menuDimensions_.get().y - 3*smallSpace - 2* controlElementRadius_.get()));

    //delete element
    deleteElem_.setPosition(tgt::ivec2((timeMenu_.getUR().x - timeMenu_.getLL().x)/2, 5 + controlElementRadius_.get()));

    //update timeline
    updateTimeline();

    //set up elements for newAnimation Mode
    confirmElem_.setPosition(tgt::ivec2(10+ controlElementRadius_.get() ,10 + controlElementRadius_.get()));
    confirmElem_.setRadius(controlElementRadius_.get());
    durationSlider_.setBarLength((menu_.getUR().x - menu_.getLL().x) * 8/10);
    durationSlider_.setIndicatorWidth(25); durationSlider_.setIndicatorHeight(25);

}

void TouchTableAnimationWidget::updateTimeline(){
    int smallSpace = 5;  //space between control elements in one set

    timelineMenu_.setLL(tgt::ivec2(timeMenu_.getUR().x - menu_.getLL().x + 10, 10));
    timelineMenu_.setUR(tgt::ivec2(menuDimensions_.get().x - 10,menuDimensions_.get().y - 3*smallSpace - 2* controlElementRadius_.get()));
    timelineMenu_.getTimeline().setDuration(timelineDuration_.get());
    timelineMenu_.getTimeline().setStepSize(timelineStepSize_.get());
}

void TouchTableAnimationWidget::updateTimelineSlider(){
    int currentPos = timelineMenu_.getTimeline().getCurrentPos();
    int contentMenuLength = (timelineMenu_.getContentMenu().getUR().x - timelineMenu_.getContentMenu().getLL().x);
    tgt::ivec2 visibleMenuSpace = tgt::ivec2(timelineMenu_.getTimeline().getTimeOffset(), timelineMenu_.getTimeline().getTimeOffset() +contentMenuLength);
    if(currentPos >= visibleMenuSpace.y){
        float newIndicatorPos = currentPos >= (timelineDuration_.get() - contentMenuLength) ? 1 : static_cast<float>(currentPos) / static_cast<float>(timelineDuration_.get() - contentMenuLength);
        timelineMenu_.setIndicatorPosition(newIndicatorPos);
    }else if(currentPos <= visibleMenuSpace.x){
        float newIndicatorPos = currentPos <= contentMenuLength ? 0 : static_cast<float>(currentPos- contentMenuLength) / static_cast<float>(timelineDuration_.get() - contentMenuLength);
        timelineMenu_.setIndicatorPosition(newIndicatorPos);
    }
}


} // namespace
