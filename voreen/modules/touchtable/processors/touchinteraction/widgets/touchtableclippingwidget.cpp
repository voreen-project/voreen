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

#include "iostream"
#include "touchtableclippingwidget.h"

#include "voreen/core/voreenapplication.h"

namespace voreen {


    const std::string TouchTableClippingWidget::loggerCat_ = "voreen.touchtable.TouchTableClippingWidget";

    TouchTableClippingWidget::TouchTableClippingWidget()
        : TouchTableMenuWidget("Scalpel-icon.png")
        , maxPlanes_("maxplanes", "Maximum Number of Planes", 3, 1, 3)
        , contrelemBoundingBox_(25, tgt::ivec2(25 +5, 25 +5),0,1.f,
        tgt::vec4(1.f,0.f,0.f,0.8f),tgt::vec4(0.5f,0.f,0.f,0.f),4,false, true)
        , contrelemAdd_(25, tgt::ivec2(3 * 25 + 10, 25 + 5), 0, 1.f, tgt::vec4(0.5f, 1.f, 0.f, 0.8f),
        tgt::vec4(0.5f, 0.f, 0.f, 0.f), 5 , false, true)
        , contrelemAddXY_(25, tgt::ivec2(3 * 25 + 10, 25 + 5), 0, 1.f, tgt::vec4(0.5f, 1.f, 0.f, 0.8f),
        tgt::vec4(0.5f, 0.f, 0.f, 0.f), 6 , false, false)
        , contrelemAddXZ_(25, tgt::ivec2(3 * 25 + 10, 25 + 5), 0, 1.f, tgt::vec4(0.5f, 1.f, 0.f, 0.8f),
        tgt::vec4(0.5f, 0.f, 0.f, 0.f), 7 , false, false)
        , contrelemAddYZ_(25, tgt::ivec2(3 * 25 + 10, 25 + 5), 0, 1.f, tgt::vec4(0.5f, 1.f, 0.f, 0.8f),
        tgt::vec4(0.5f, 0.f, 0.f, 0.f), 8 , false, false)
        , boundingBoxOn_("boundingBoxOn", "BoundingBox on" , false)
        , cameraProp_("cameraProp", "Camera Property")
        , exclusiveModeSet_(false)
        , freeClippingLineParams_(freeClippingParameters(-1,tgt::ivec2(0,0),tgt::ivec2(0,0)))
        , sliderLength_(200)
        , sliderRange_("sliderrange", "Slider Range", tgt::vec2(-1.f, 1.f), tgt::vec2(-FLT_MAX), tgt::vec2(FLT_MAX), Processor::VALID)
    {

        addProperty(maxPlanes_);

        //initialize and add the properties and push them into the vector propertySets
        for(int i = 0; i < 3; ++i) {
            ClippingPropertySet* clipPropSet = new ClippingPropertySet(i);
            addProperty(clipPropSet->color_);
            addProperty(clipPropSet->renderBoolProperty_);
            addProperty(clipPropSet->toggleBoolProperty_);
            addProperty(clipPropSet->planeNormalProperty_);
            addProperty(clipPropSet->planePositionProperty_);
            addProperty(clipPropSet->invertProperty_);
            addProperty(clipPropSet->renderManipulation_);
            propertySets_.push_back(clipPropSet);
        }

        addProperty(boundingBoxOn_);
        addProperty(cameraProp_);

        //initialize control elements
        //clippingElements
        tgt::ivec2 startPos = tgt::ivec2(0, 50);
        for (int i = ((int) propertySets_.size())-1; i >=0; --i) {
            inactiveClippingSets_.push_back(ClippingSet(tgt::ivec2(getPosition()),
                startPos + tgt::ivec2(0, i * (15 + 2 * controlElementRadius_.get())),propertySets_.at(i)->color_.get(), controlElementRadius_.get(), propertySets_.at(i), i, sliderLength_,i));
        }

        //onChange
        isActive_.onChange(MemberFunctionCallback<TouchTableClippingWidget>(this, &TouchTableClippingWidget::deactivateWidget));

        std::vector<ClippingPropertySet*>::iterator propIter;
        for (propIter = propertySets_.begin(); propIter != propertySets_.end(); ++propIter)
        {
            ClippingPropertySet* currentPropertySet = *propIter;
            currentPropertySet->color_.onChange(MemberFunctionCallback<TouchTableClippingWidget>(this,&TouchTableClippingWidget::updateClippingSetColor));
            currentPropertySet->planePositionProperty_.onChange(MemberFunctionCallback<TouchTableClippingWidget>(this,&TouchTableClippingWidget::updateSliderPosition));
        }

        addProperty(sliderRange_);
        sliderRange_.setVisibleFlag(false);
    }

    /*
    * TouchTableClippingWidget desctructor
    */
    TouchTableClippingWidget::~TouchTableClippingWidget(){
        std::vector<ClippingPropertySet*>::iterator propertySetIter;
        for(propertySetIter = propertySets_.begin(); propertySetIter != propertySets_.end(); ++propertySetIter){
            delete *propertySetIter;
        }
    }

    /*
    * textures are set to their according control elements and inital values are set
    */
    void TouchTableClippingWidget::initialize() {
        TouchTableMenuWidget::initialize();

        //load textures for control elements
        symbolTexClipping_= TexMgr.load("eye.png");
        symbolTexBounding_=TexMgr.load("bounding.png");
        symbolTexRenderGeometry_=TexMgr.load("renderGeometry.png");
        symbolTexRenderManipulation_ = TexMgr.load("renderManipulation.png");
        symbolTexInvert_=TexMgr.load("invert.png");
        symbolTexDelete_=TexMgr.load("delete.png");
        symbolTexFreeClipping_= TexMgr.load("razorblade.png");
        symbolTexAdd_ = TexMgr.load("add.png");
        symbolTexAddXY_ = TexMgr.load("z_up.png");
        symbolTexAddXZ_ = TexMgr.load("y_up.png");
        symbolTexAddYZ_ = TexMgr.load("x_up.png");

        setSymbolTexClippingSet(symbolTexClipping_,symbolTexRenderGeometry_, symbolTexInvert_, symbolTexFreeClipping_, symbolTexDelete_, symbolTexRenderManipulation_);
        contrelemBoundingBox_.setSymbolTexture(symbolTexBounding_);
        contrelemAdd_.setSymbolTexture(symbolTexAdd_);
        contrelemAddXY_.setSymbolTexture(symbolTexAddXY_);
        contrelemAddXZ_.setSymbolTexture(symbolTexAddXZ_);
        contrelemAddYZ_.setSymbolTexture(symbolTexAddYZ_);

        /*//set properties false on default
        std::vector<ClippingPropertySet*>::iterator propIter;
        for(propIter = propertySets_.begin(); propIter!=propertySets_.end(); ++propIter){
            ClippingPropertySet* currentPropertySet = *propIter;
            currentPropertySet->toggleBoolProperty_.set(false);
            currentPropertySet->renderBoolProperty_.set(false);
        }*/

        for(std::vector<ClippingSet>::iterator clippingIter = inactiveClippingSets_.begin(); clippingIter != inactiveClippingSets_.end(); ){
            if(clippingIter->propertySet_->toggleBoolProperty_.get()){
                clippingIter->toggleClippingPlane_.setActive(true);
                clippingIter->renderClippingPlaneGeometry_.setActive(true);
                clippingSets_.push_back(*clippingIter);
                clippingIter = inactiveClippingSets_.erase(clippingIter);
            }else{
                ++clippingIter;
            }
        }

    }

    void TouchTableClippingWidget::deinitialize() {

        //clean up
        TexMgr.dispose(symbolTexClipping_);
        TexMgr.dispose(symbolTexBounding_);
        TexMgr.dispose(symbolTexRenderGeometry_);
        TexMgr.dispose(symbolTexRenderManipulation_);
        TexMgr.dispose(symbolTexInvert_);
        TexMgr.dispose(symbolTexDelete_);
        TexMgr.dispose(symbolTexFreeClipping_);
        TexMgr.dispose(symbolTexAdd_);
        TexMgr.dispose(symbolTexAddXY_);
        TexMgr.dispose(symbolTexAddXZ_);
        TexMgr.dispose(symbolTexAddYZ_);

        TouchTableMenuWidget::deinitialize();
    }

    void TouchTableClippingWidget::setSymbolTexClippingSet(tgt::Texture* texClipping, tgt::Texture* texRenderGeometry, tgt::Texture* texInvert, tgt::Texture* texFreeClipping, tgt::Texture* texDelete, tgt::Texture* texRenderManipulation){
        std::vector<ClippingSet>::iterator clippingIter;
        for(clippingIter = clippingSets_.begin(); clippingIter != clippingSets_.end(); ++clippingIter){
            clippingIter->toggleClippingPlane_.setSymbolTexture(texClipping);

            clippingIter->renderClippingPlaneGeometry_.setSymbolTexture(texRenderGeometry);
            clippingIter->renderManipulation_.setSymbolTexture(texRenderManipulation);
            clippingIter->invertClippingPlane_.setSymbolTexture(texInvert);
            clippingIter->freeClipping_.setSymbolTexture(texFreeClipping);
            clippingIter->deletePlane_.setSymbolTexture(texDelete);

        }
        for(clippingIter = inactiveClippingSets_.begin(); clippingIter != inactiveClippingSets_.end(); ++clippingIter){
            clippingIter->toggleClippingPlane_.setSymbolTexture(texClipping);
            clippingIter->renderClippingPlaneGeometry_.setSymbolTexture(texRenderGeometry);
            clippingIter->renderManipulation_.setSymbolTexture(texRenderManipulation);
            clippingIter->invertClippingPlane_.setSymbolTexture(texInvert);
            clippingIter->freeClipping_.setSymbolTexture(texFreeClipping);
            clippingIter->deletePlane_.setSymbolTexture(texDelete);

        }
    }

    void TouchTableClippingWidget::deactivateWidget(){
        if(isActive_.get()){
            std::vector<ClippingSet>& clippingSets = clippingSets_;
            std::vector<ClippingSet>::iterator clippingSetIter;
            for(clippingSetIter = clippingSets.begin(); clippingSetIter != clippingSets.end(); ++clippingSetIter){
                clippingSetIter->propertySet_->renderBoolProperty_.set(clippingSetIter->renderClippingPlaneGeometry_.isActive());
                clippingSetIter->propertySet_->renderManipulation_.set(clippingSetIter->renderManipulation_.isActive());
            }
        }else{
            std::vector<ClippingPropertySet*>::iterator propertySetIter;
            for(propertySetIter = propertySets_.begin(); propertySetIter != propertySets_.end(); ++propertySetIter){
                ClippingPropertySet* currentPropSet = *propertySetIter;
                currentPropSet->renderBoolProperty_.set(false);
                currentPropSet->renderManipulation_.set(false);

            }
        }
    }

    /*void TouchTableClippingWidget::setMenuSize(){
    int yLength = (((int)clippingSets_.size())+1) * (2 * controlElementRadius_+10) + 40;
    int xLength = 10 * controlElementRadius_ + 65 + sliderLength_.get();
    menu_.setAnchor(getPosition());
    //menu_->setLL(tgt::ivec2(menu_->getAnchor().x + widgetRadius_,menu_->getAnchor().y+ widgetRadius_));
    menu_.setUR(tgt::ivec2(menu_.getLL().x +xLength, menu_.getLL().y + yLength));
    }*/

    void TouchTableClippingWidget::handleTouchPoints(const std::deque<tgt::TouchPoint>& tp) {

        //iterate over all touchpoints
        std::deque<tgt::TouchPoint>::const_iterator tpIter;
        for(tpIter = tp.begin(); tpIter != tp.end(); ++tpIter){
            //convert tp position, so it's relative to the menu coordinates
            tgt::ivec2 tpPos = convertToMenuCoordinates(tpIter->pos());

            //normal Mode <--> exclusiveMode
            if(!exclusiveModeSet_){
                //get touchpoint state
                tgt::TouchPoint::State state = tpIter->state();
                if(state == tgt::TouchPoint::TouchPointPressed)
                    handleTouchPointPressed(*tpIter);
                else if(state == tgt::TouchPoint::TouchPointReleased)
                    handleTouchPointReleased(tpPos, tpIter->id());
                else if(state= tgt::TouchPoint::TouchPointMoved)
                    handleTouchPointMoved(*tpIter);
            }else{
                tgt::TouchPoint::State state = tpIter->state();
                if(state == tgt::TouchPoint::TouchPointPressed)
                    handleExclusiveTouchPointPressed(*tpIter);
                else if(state == tgt::TouchPoint::TouchPointMoved){
                    handleExclusiveTouchPointMoved(*tpIter);
                    render();

                }
                else if(state == tgt::TouchPoint::TouchPointReleased)
                    handleExclusiveTouchPointReleased(*tpIter);
            }


        }
    }

    void TouchTableClippingWidget::handleExclusiveTouchPointPressed(const tgt::TouchPoint& tp){
        float yOutportSize = static_cast<float>(overlay_->getOutportSize().y);
        //save the first touchpoint in a map and do not handle any other touchpoints until released
        if(freeClippingLineParams_.id_ == -1){
            tgt::ivec2 startingPos = tgt::ivec2((int) tp.pos().x, (int) (yOutportSize - tp.pos().y));
            freeClippingLineParams_.id_ = tp.id();
            freeClippingLineParams_.startPos_= startingPos;
            freeClippingLineParams_.currentPos_ = startingPos;
        }
    }

    void TouchTableClippingWidget::handleExclusiveTouchPointMoved(const tgt::TouchPoint& tp){
        float yOutportSize = static_cast<float>(overlay_->getOutportSize().y);

        if(freeClippingLineParams_.id_ == tp.id()){
            freeClippingLineParams_.currentPos_ = tgt::ivec2((int) tp.pos().x, (int) (yOutportSize - tp.pos().y)) ;

        }
    }

    void TouchTableClippingWidget::handleExclusiveTouchPointReleased(const tgt::TouchPoint& tp){

        //cancel drawing clipping line by hitting add button again
        tgt::ivec2 tpPos = convertToMenuCoordinates(tp.pos());

        if(hitsControlElement(tpPos, freeClippingLineParams_.activeClippingSet_->freeClipping_) && freeClippingLineParams_.activeClippingSet_->freeClipping_.isActive()){
                freeClippingLineParams_.id_ = -1;
                //reset exclusive mode
                overlay_->setExclusiveMode(false);
                exclusiveModeSet_ = false;
                freeClippingLineParams_.activeClippingSet_->freeClipping_.setActive(false);
                //reset pointer to active clipping set in freeClippingLineParams
                freeClippingLineParams_.activeClippingSet_= 0;
                //ealry return because clipping was canceld
                return;
        }

        if(hitsControlElement(tpPos, contrelemAdd_) && contrelemAdd_.isActive()){

            //reset pointer to active clipping set in freeClippingLineParams
            freeClippingLineParams_.activeClippingSet_->toggleClippingPlane_.setActive(false);
            freeClippingLineParams_.activeClippingSet_->renderClippingPlaneGeometry_.setActive(false);
            freeClippingLineParams_.activeClippingSet_= 0;
            freeClippingLineParams_.id_ = -1;
            //reset exclusive mode
            overlay_->setExclusiveMode(false);
            exclusiveModeSet_ = false;
            contrelemAdd_.setActive(false);

        }else{
            if(freeClippingLineParams_.id_== tp.id()){
                freeClippingLineParams_.id_ = -1;
                overlay_->setExclusiveMode(false);
                exclusiveModeSet_=false;
                freeClippingLineParams_.activeClippingSet_->toggleClippingPlane_.setActive(true);
                freeClippingLineParams_.activeClippingSet_->propertySet_->toggleBoolProperty_.set(true);
                freeClippingLineParams_.activeClippingSet_->renderClippingPlaneGeometry_.setActive(true);
                freeClippingLineParams_.activeClippingSet_->propertySet_->renderBoolProperty_.set(true);
                for(std::vector<ClippingSet>::iterator i = clippingSets_.begin(); i != clippingSets_.end(); ++i){
                    i->propertySet_->renderManipulation_.set(false);
                    i->renderManipulation_.setActive(false);
                }
                freeClippingLineParams_.activeClippingSet_->renderManipulation_.setActive(true);
                freeClippingLineParams_.activeClippingSet_->propertySet_->renderManipulation_.set(true);

                setClippingPropsFreeClipping(freeClippingLineParams_.currentPos_ ,freeClippingLineParams_.startPos_,freeClippingLineParams_.activeClippingSet_->propertySet_->planePositionProperty_,freeClippingLineParams_.activeClippingSet_->propertySet_->planeNormalProperty_);

                freeClippingLineParams_.activeClippingSet_->freeClipping_.setActive(false);
                //if add button was pressed, swap activated clipping set from inactive clipping sets to active clipping sets
                if(contrelemAdd_.isActive()){
                    inactiveClippingSets_.pop_back();
                    clippingSets_.push_back(*freeClippingLineParams_.activeClippingSet_);
                    contrelemAdd_.setActive(false);
                }
                freeClippingLineParams_.activeClippingSet_= 0;

                updateMenuCoordinates();
                updateSliderPosition();

            }
        }
    }

    void TouchTableClippingWidget::handleTouchPointPressed(const tgt::TouchPoint tp){
        //correct the position of the TouchPoint
        tgt::ivec2 tpPos = convertToMenuCoordinates(tp.pos());

        //check if one of the controlelements has been hit
        std::vector<ClippingSet>::iterator clipIter;
        for(clipIter = clippingSets_.begin(); clipIter != clippingSets_.end(); ++clipIter){
            //slider indicator has been hit

            if(clipIter->slider_.checkIndicatorHit(tpPos)){
                touchPointSliderMapping_.insert(std::make_pair(tp.id(), &(*clipIter)));
            }
        }
    }

    void TouchTableClippingWidget::handleTouchPointMoved(const tgt::TouchPoint tp){
        //correct the position of the TouchPoint
        tgt::ivec2 tpPos = convertToMenuCoordinates(tp.pos());

        std::map<int, ClippingSet*>::iterator iter;
        iter = touchPointSliderMapping_.find(tp.id());

        if (iter != touchPointSliderMapping_.end()) {
            ClippingSet* clippingSet = iter->second;
            TouchTableSlider& slider = clippingSet->slider_;
            //update position of slider
            slider.updateIndicatorPosition(tpPos);
            invalidate();
            //update according properties in the clipping set
            //float sliderPosition = slider.getIndicatorPos() *2 -1;

            updateSliderRange();

            float sliderPosition = slider.getIndicatorPos() * (sliderRange_.get().y - sliderRange_.get().x) + sliderRange_.get().x;

            clippingSet->propertySet_->planePositionProperty_.set(sliderPosition);


        }
    }

    void TouchTableClippingWidget::handleTouchPointReleased(tgt::ivec2 tpPos, int tpId){
        //check if one of the controlelements has been hit
        std::vector<ClippingSet>::iterator clipIter;
        for(clipIter = clippingSets_.begin(); clipIter != clippingSets_.end(); ++clipIter){
            //check whether the slider has been released and erase tp from mapping
            if(touchPointSliderMapping_.find(tpId) != touchPointSliderMapping_.end()){
                touchPointSliderMapping_.erase(touchPointSliderMapping_.find(tpId));
                break;
            }
            //Toggle Control Element Hit
            TouchTableControlElement* currentControlElemIter = &(clipIter->toggleClippingPlane_);
            //if(tgt::length(tpPos - currentControlElemIter->getPos()) <= currentControlElemIter->getRadius()){
            if(hitsControlElement(tpPos, *currentControlElemIter)){
                bool active = currentControlElemIter->isActive();
                clipIter->propertySet_->toggleBoolProperty_.set(!active);
                currentControlElemIter->setActive(!active);
                //set associated render geometry control element accordingly
                clipIter->renderClippingPlaneGeometry_.setActive(!active);
                clipIter->propertySet_->renderBoolProperty_.set(!active);
                clipIter->renderManipulation_.setActive(false);
                clipIter->propertySet_->renderManipulation_.set(false);
                continue; //no other element from this set could have been hit
            }
            //rendering Control Element Hit
            currentControlElemIter = &(clipIter->renderClippingPlaneGeometry_);
            if(hitsControlElement(tpPos, *currentControlElemIter)){
                bool active = currentControlElemIter->isActive();
                //((TouchTableClippingWidget*) widget_)->setControlElementBoolOn(!active, currentControlElemIter->getId()); //cast widget_ to ClippingWidget throughout the whole class?!
                clipIter->propertySet_->renderBoolProperty_.set(!active);
                currentControlElemIter->setActive(!active);
                continue; //no other element from this set could have been hit
            }
            //render manipulation control element hit
            currentControlElemIter = &(clipIter->renderManipulation_);
            if(hitsControlElement(tpPos, *currentControlElemIter)){
                bool active = currentControlElemIter->isActive();
                //deactivate rendering of manipulation in every other clipping set
                for(std::vector<ClippingSet>::iterator i = clippingSets_.begin(); i != clippingSets_.end(); ++i){
                    i->renderManipulation_.setActive(false);
                    i->propertySet_->renderManipulation_.set(false);
                }
                clipIter->propertySet_->renderManipulation_.set(!active);
                currentControlElemIter->setActive(!active);
                invalidate();
                continue;
            }
            //invert Control Element Hit
            currentControlElemIter = &(clipIter->invertClippingPlane_);
            if(hitsControlElement(tpPos, *currentControlElemIter)){
                bool active = currentControlElemIter->isActive();
                clipIter->propertySet_->invertProperty_.set(!active);
                //clipIter->propertySet_->planeNormalProperty_.set((-1.f) * clipIter->propertySet_->planeNormalProperty_.get());
                //clipIter->propertySet_->planePositionProperty_.set((-1.f) * clipIter->propertySet_->planePositionProperty_.get());
                currentControlElemIter->setActive(!active);
            }
            //free clipping element hit
            currentControlElemIter = &(clipIter->freeClipping_);
            if(hitsControlElement(tpPos, *currentControlElemIter)){
                currentControlElemIter->setActive(true);
                //((TouchTableClippingWidget*) widget_)->setControlElementBoolOn(!active, currentControlElemIter->getId()); //cast widget_ to ClippingWidget throughout the whole class?!
                freeClippingLineParams_.activeClippingSet_ = &(*clipIter);
                overlay_->setExclusiveMode(true);
                exclusiveModeSet_ = true;
            }
            //delete plane element hit
            currentControlElemIter = &(clipIter->deletePlane_);
            if(hitsControlElement(tpPos, *currentControlElemIter)){
                if(clippingSets_.size() !=0){
                    //deactivate all functions of the object being deleted
                    clipIter->toggleClippingPlane_.setActive(false);
                    clipIter->propertySet_->toggleBoolProperty_.set(false);
                    clipIter->renderClippingPlaneGeometry_.setActive(false);
                    clipIter->propertySet_->renderBoolProperty_.set(false);
                    clipIter->invertClippingPlane_.setActive(false);
                    clipIter->propertySet_->renderManipulation_.set(false);
                    //clipIter->propertySet_->invertBoolProperty_.set(false);
                    //push deleted clipping set into vector with inactive sets
                    inactiveClippingSets_.push_back(*clipIter);
                    //delete clipping set from active clipping sets
                    clippingSets_.erase(clipIter);
                    //update positions of active clipping sets
                    std::vector<ClippingSet>::iterator newClippingSetsIter;
                    int posInMenu =0;
                    for(newClippingSetsIter = clippingSets_.begin(); newClippingSetsIter != clippingSets_.end(); ++newClippingSetsIter){
                        updateClippingSetPosition(*newClippingSetsIter, posInMenu++);
                    }
                    updateMenuCoordinates();
                    render();
                    break;
                }
            }

        }
        //add element hit
        if(hitsControlElement(tpPos, contrelemAdd_)){
            if ((inactiveClippingSets_.size() != 0) && (clippingSets_.size() < maxPlanes_.get())) {
                contrelemAdd_.setActive(true);
                //take a clipping set from the inactive ones and add it to the active ones
                ClippingSet* activatedClippingSet = &inactiveClippingSets_.back();
                //render geometry and toggle are activated by default
                activatedClippingSet->toggleClippingPlane_.setActive(true);
                activatedClippingSet->renderClippingPlaneGeometry_.setActive(true);

                updateClippingSetPosition(*activatedClippingSet, ((int)clippingSets_.size()));

                //add activated clipping set to free clipping parameters
                freeClippingLineParams_.activeClippingSet_ = activatedClippingSet;
                //set exclusive mode
                overlay_->setExclusiveMode(true);
                exclusiveModeSet_ = true;
            }
        }
        //addXY element hit
        if(hitsControlElement(tpPos, contrelemAddXY_)){
            if ((inactiveClippingSets_.size() != 0) && (clippingSets_.size() < maxPlanes_.get())) {
                ClippingSet& activatedClippingSet = inactiveClippingSets_.back();
                //activate clipping in activated clipping set
                activatedClippingSet.propertySet_->toggleBoolProperty_.set(true);
                activatedClippingSet.toggleClippingPlane_.setActive(true);
                activatedClippingSet.propertySet_->renderBoolProperty_.set(true);
                activatedClippingSet.renderClippingPlaneGeometry_.setActive(true);
                for(std::vector<ClippingSet>::iterator i = clippingSets_.begin(); i != clippingSets_.end(); ++i){
                    i->propertySet_->renderManipulation_.set(false);
                    i->renderManipulation_.setActive(false);
                }
                activatedClippingSet.renderManipulation_.setActive(true);
                activatedClippingSet.propertySet_->renderManipulation_.set(true);
                //set parameters (planePosition and planeNormal)
                setClippingPropsFreeClipping(tgt::ivec2(5,overlay_->getOutportSize().y/2), tgt::ivec2(overlay_->getOutportSize().x-5,overlay_->getOutportSize().y/2), activatedClippingSet.propertySet_->planePositionProperty_, activatedClippingSet.propertySet_->planeNormalProperty_);
                activatedClippingSet.propertySet_->planePositionProperty_.set(static_cast<float>(0.5 * (sliderRange_.get().y - sliderRange_.get().x) + sliderRange_.get().x));
                //move activated clipping set from inactive sets to active ones
                inactiveClippingSets_.pop_back();
                clippingSets_.push_back(activatedClippingSet);
                updateMenuCoordinates();
                invalidate();
            }
        }
        //addXZ element hit
        if(hitsControlElement(tpPos, contrelemAddXZ_)){
            if ((inactiveClippingSets_.size() != 0) && (clippingSets_.size() < maxPlanes_.get())) {
                ClippingSet& activatedClippingSet = inactiveClippingSets_.back();
                //activate clipping in activated clipping set
                activatedClippingSet.propertySet_->toggleBoolProperty_.set(true);
                activatedClippingSet.toggleClippingPlane_.setActive(true);
                activatedClippingSet.propertySet_->renderBoolProperty_.set(true);
                activatedClippingSet.renderClippingPlaneGeometry_.setActive(true);
                for(std::vector<ClippingSet>::iterator i = clippingSets_.begin(); i != clippingSets_.end(); ++i){
                    i->propertySet_->renderManipulation_.set(false);
                    i->renderManipulation_.setActive(false);
                }
                activatedClippingSet.renderManipulation_.setActive(true);
                activatedClippingSet.propertySet_->renderManipulation_.set(true);
                //set parameters (planePosition and planeNormal)
                setClippingPropsFreeClipping(tgt::ivec2(overlay_->getOutportSize().x/2,5), tgt::ivec2(overlay_->getOutportSize().x/2,overlay_->getOutportSize().y-5), activatedClippingSet.propertySet_->planePositionProperty_, activatedClippingSet.propertySet_->planeNormalProperty_);
                activatedClippingSet.propertySet_->planePositionProperty_.set(static_cast<const float>(0.5 * (sliderRange_.get().y - sliderRange_.get().x) + sliderRange_.get().x));
                //move activated clipping set from inactive sets to active ones
                inactiveClippingSets_.pop_back();
                clippingSets_.push_back(activatedClippingSet);
                updateMenuCoordinates();
                invalidate();
            }
        }
        //addYZ element hit
        if(hitsControlElement(tpPos, contrelemAddYZ_)){
            if ((inactiveClippingSets_.size() != 0) && (clippingSets_.size() < maxPlanes_.get())) {
                ClippingSet& activatedClippingSet = inactiveClippingSets_.back();
                //activate clipping in activated clipping set
                activatedClippingSet.propertySet_->toggleBoolProperty_.set(true);
                activatedClippingSet.toggleClippingPlane_.setActive(true);
                activatedClippingSet.propertySet_->renderBoolProperty_.set(true);
                activatedClippingSet.renderClippingPlaneGeometry_.setActive(true);
                for(std::vector<ClippingSet>::iterator i = clippingSets_.begin(); i != clippingSets_.end(); ++i){
                    i->propertySet_->renderManipulation_.set(false);
                    i->renderManipulation_.setActive(false);
                }
                activatedClippingSet.renderManipulation_.setActive(true);
                activatedClippingSet.propertySet_->renderManipulation_.set(true);
                //set parameters (planePosition and planeNormal)
                activatedClippingSet.propertySet_->planeNormalProperty_.set(cameraProp_.get().getLook());
                activatedClippingSet.propertySet_->planePositionProperty_.set(static_cast<const float>(0.5 * (sliderRange_.get().y - sliderRange_.get().x) + sliderRange_.get().x));
                //move activated clipping set from inactive sets to active ones
                inactiveClippingSets_.pop_back();
                clippingSets_.push_back(activatedClippingSet);
                updateMenuCoordinates();
                invalidate();
            }
        }
        //bounding box element hit
        if(hitsControlElement(tpPos, contrelemBoundingBox_)){
            bool active = contrelemBoundingBox_.isActive();
            contrelemBoundingBox_.setActive(!active);
            boundingBoxOn_.set(!active);
        }

    }

    void TouchTableClippingWidget::setClippingPropsFreeClipping(tgt::ivec2 startPos,tgt::ivec2 endPos, FloatProperty& planePosition, FloatVec3Property& plainNormal){
        const TouchtableModule* module = dynamic_cast<const TouchtableModule*>(VoreenApplication::app()->getModule("touchtable"));
        bool transform = false;
        if (module)
           transform = module->usePerspectiveTransformClipping();

        const tgt::Camera& camera = cameraProp_.get();

        //normalize coordinates
        tgt::vec2 nEndPos = tgt::vec2(endPos) / tgt::vec2(overlay_->getOutportSize());
        nEndPos = nEndPos * 2.f - tgt::vec2(1.f);
        tgt::vec2 nStartPos = tgt::vec2(startPos) / tgt::vec2(overlay_->getOutportSize());
        nStartPos = nStartPos * 2.f - tgt::vec2(1.f);

        //get inverse view and projection matrices
        tgt::mat4 projectionMatrixInverse;
        camera.getProjectionMatrix(overlay_->getOutportSize()).invert(projectionMatrixInverse);
        tgt::mat4 viewMatrixInverse = camera.getViewMatrixInverse();

        //get near/far plane and look vector of camera
        float nearDist = camera.getNearDist();
        float farDist = camera.getFarDist();
        tgt::vec3 look = camera.getLook();

        //compute plane normal and position
        tgt::vec3 newPlaneNormal;
        float planePos;

        if (transform) {

            tgt::vec4 startWorld = viewMatrixInverse *  (projectionMatrixInverse *  tgt::vec4(nStartPos, -1.f, 1.f));
            tgt::vec4 endWorld = viewMatrixInverse *  (projectionMatrixInverse *  tgt::vec4(nEndPos, -1.f, 1.f));

            tgt::vec4 endBackWorld = viewMatrixInverse * (projectionMatrixInverse * tgt::vec4(nEndPos, 1.f, 1.f));

            startWorld *= 1.f/startWorld.w;
            endWorld *= 1.f/endWorld.w;
            endBackWorld *= 1.f/endBackWorld.w;

            tgt::vec4 lineOne = endWorld - startWorld;
            tgt::vec4 lineTwo = endBackWorld - endWorld;

            newPlaneNormal = tgt::vec3(tgt::normalize(tgt::cross(lineOne.xyz(), lineTwo.xyz())));

            tgt::vec4 point = (startWorld + endWorld + endBackWorld) / 3.f;
            planePos = tgt::dot(point.xyz(), newPlaneNormal);

            //planePos = tgt::clamp(planePos, -1.01f, 1.01f);
            planePos = tgt::clamp(planePos, sliderRange_.get().x, sliderRange_.get().y);
        }
        else {
            tgt::vec4 startWorld = viewMatrixInverse *  (projectionMatrixInverse *  tgt::vec4(nStartPos, nearDist, 1.f));
            tgt::vec4 endWorld = viewMatrixInverse *  (projectionMatrixInverse *  tgt::vec4(nEndPos, nearDist, 1.f));


            /*startWorld *= 1.f/startWorld.w;
            endWorld *= 1.f/endWorld.w;*/

            tgt::vec4 lineWorld = endWorld - startWorld;
            newPlaneNormal = tgt::vec3(tgt::normalize(tgt::cross(tgt::vec3(lineWorld.x, lineWorld.y,lineWorld.z),look)));

            tgt::vec4 point = tgt::vec4(0.5) * (startWorld + endWorld);
            planePos = tgt::dot(point.xyz(), newPlaneNormal);
        }

        plainNormal.set(newPlaneNormal);

        planePosition.set(planePos);
    }

    void TouchTableClippingWidget::renderClippingLine(tgt::ivec2 currentPos, tgt::ivec2 startingPos, tgt::vec4 color){

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

    glLineWidth(4.f);
    glColor3f(color.r, color.g, color.b);

    glBegin(GL_LINES);
        glVertex2f(static_cast<float>(startingPos.x), static_cast<float>(startingPos.y));
        glVertex2f(static_cast<float>(currentPos.x), static_cast<float>(currentPos.y));
    glEnd();

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

    void TouchTableClippingWidget::updateClippingSetPosition(ClippingSet& clippingSet, int posInMenu){
        //update control element positions
        //update the position of the menu the clipping set is located in
        tgt::ivec2 posLL =menu_.getLL() +  tgt::ivec2(5,posInMenu * (2 * controlElementRadius_.get() +15)) + tgt::ivec2(0,2*controlElementRadius_. get() + 10);
        clippingSet.menu_.setLL(posLL);
        clippingSet.menu_.setUR(tgt::ivec2(menu_.getUR().x - 5 ,posLL.y + 10+2*controlElementRadius_.get()));
        //int yPos = (overlay->getWidgetRadius() + 15 + controlElementRadius_) + posInMenu * (15 + 2 * controlElementRadius_);
        int yPos = (3* controlElementRadius_.get() +15) + posInMenu * (15 + 2 * controlElementRadius_.get());

        clippingSet.toggleClippingPlane_.setRadius(controlElementRadius_.get());
        clippingSet.toggleClippingPlane_.setPosition(tgt::ivec2(5 + controlElementRadius_.get(), yPos));

        clippingSet.renderClippingPlaneGeometry_.setRadius(controlElementRadius_.get());
        clippingSet.renderClippingPlaneGeometry_.setPosition(tgt::ivec2(10+3*controlElementRadius_.get(), yPos));

        clippingSet.renderManipulation_.setRadius(controlElementRadius_.get());
        clippingSet.renderManipulation_.setPosition(tgt::ivec2(15+5*controlElementRadius_.get(), yPos));

        clippingSet.invertClippingPlane_.setRadius(controlElementRadius_.get());
        clippingSet.invertClippingPlane_.setPosition(tgt::ivec2(20 + 7*controlElementRadius_.get(), yPos));

        clippingSet.freeClipping_.setRadius(controlElementRadius_.get());
        clippingSet.freeClipping_.setPosition(tgt::ivec2(25+9*controlElementRadius_.get(), yPos));

        clippingSet.deletePlane_.setRadius(controlElementRadius_.get());
        clippingSet.deletePlane_.setPosition(tgt::ivec2(30+11*controlElementRadius_.get(), yPos));

        sliderLength_ = clippingSet.menu_.getUR().x - clippingSet.menu_.getLL().x - (40) - 12 * controlElementRadius_.get();
        clippingSet.slider_.setIndicatorHeight(controlElementRadius_.get() + controlElementRadius_.get()/2);
        clippingSet.slider_.setIndicatorWidth(controlElementRadius_.get());
        clippingSet.slider_.setBarOrigin(tgt::ivec2(35+12*controlElementRadius_.get(), yPos));
        clippingSet.slider_.setBarWidth(controlElementRadius_.get() *3 / 4);
        clippingSet.slider_.setBarLength(sliderLength_);

    }

    void TouchTableClippingWidget::updateRenderControlElements(){
        std::vector<ClippingSet>::iterator clippingIter;
        for(clippingIter = clippingSets_.begin(); clippingIter != clippingSets_.end(); ++clippingIter){
            clippingIter->renderClippingPlaneGeometry_.setActive(clippingIter->propertySet_->renderBoolProperty_.get());
        }
    }

    void TouchTableClippingWidget::updateToggleControlElements(){
        std::vector<ClippingSet>::iterator clippingIter;
        for(clippingIter = clippingSets_.begin(); clippingIter != clippingSets_.end(); ++clippingIter){
            clippingIter->toggleClippingPlane_.setActive(clippingIter->propertySet_->toggleBoolProperty_.get());
        }

    }

    void TouchTableClippingWidget::updateSliderRange() {

        //get scene bounds, check every corner and compute slider range
        tgt::vec3 llf = cameraProp_.getSceneBounds().getLLF();
        tgt::vec3 urb = cameraProp_.getSceneBounds().getURB();
        float length = std::max(tgt::length(llf), tgt::length(urb));

        tgt::vec3 llb = tgt::vec3(llf.x, llf.y, urb.z);
        tgt::vec3 urf = tgt::vec3(urb.x, urb.y, llf.z);
        length = std::max(length, std::max(tgt::length(llb), tgt::length(urf)));

        tgt::vec3 lrf = tgt::vec3(llf.x, urb.y, llf.z);
        tgt::vec3 ulb = tgt::vec3(urb.x, llf.y, urb.z);
        length = std::max(length, std::max(tgt::length(lrf), tgt::length(ulb)));

        tgt::vec3 lrb = tgt::vec3(llf.x, urb.y, urb.z);
        tgt::vec3 ulf = tgt::vec3(urb.x, llf.y, llf.z);
        length = std::max(length, std::max(tgt::length(lrb), tgt::length(ulf)));

        if (length == 0.f)
            length = 1.f;

        sliderRange_.set(tgt::vec2(-length, length));
    }

    void TouchTableClippingWidget::updateSliderPosition() {

        updateSliderRange();

        std::vector<ClippingSet>::iterator i;
        for(i = clippingSets_.begin(); i != clippingSets_.end(); ++i){
            //update position of each slider in every clippingSet
            //float indicatorPos = (clippingIter->propertySet_->planePositionProperty_.get() +1)/2;

            float indicatorPos = (i->propertySet_->planePositionProperty_.get() - sliderRange_.get().x) / (sliderRange_.get().y - sliderRange_.get().x);
            i->slider_.setIndicatorPos(indicatorPos);
        }
        invalidate();
    }

    void TouchTableClippingWidget::updateClippingSetColor(){
        std::vector<ClippingSet>::iterator clippingIter;
        for(clippingIter = clippingSets_.begin(); clippingIter != clippingSets_.end(); ++clippingIter){
            //set color of control elements
            clippingIter->invertClippingPlane_.setColorMod(clippingIter->propertySet_->color_.get());
            clippingIter->toggleClippingPlane_.setColorMod(clippingIter->propertySet_->color_.get());
            clippingIter->renderClippingPlaneGeometry_.setColorMod(clippingIter->propertySet_->color_.get());
            clippingIter->freeClipping_.setColorMod(clippingIter->propertySet_->color_.get());
            clippingIter->deletePlane_.setColorMod(clippingIter->propertySet_->color_.get());
            clippingIter->renderManipulation_.setColorMod(clippingIter->propertySet_->color_.get());

            //set color of menus of clippingsets according to their color property
            clippingIter->menu_.setColor(clippingIter->propertySet_->color_.get() - tgt::vec4(0.f,0.f,0.f,0.3f));
        }
        clippingIter;
        for(clippingIter = inactiveClippingSets_.begin(); clippingIter != inactiveClippingSets_.end(); ++clippingIter){
            //set color of control elements
            clippingIter->invertClippingPlane_.setColorMod(clippingIter->propertySet_->color_.get());
            clippingIter->toggleClippingPlane_.setColorMod(clippingIter->propertySet_->color_.get());
            clippingIter->renderClippingPlaneGeometry_.setColorMod(clippingIter->propertySet_->color_.get());
            clippingIter->freeClipping_.setColorMod(clippingIter->propertySet_->color_.get());
            clippingIter->deletePlane_.setColorMod(clippingIter->propertySet_->color_.get());
            clippingIter->renderManipulation_.setColorMod(clippingIter->propertySet_->color_.get());

            //set color of menus of clippingsets according to their color property
            clippingIter->menu_.setColor(clippingIter->propertySet_->color_.get() - tgt::vec4(0.f,0.f,0.f,0.3f));
        }
    }

    void TouchTableClippingWidget::updateMenuCoordinates() {
        //update y size of menu dimensions
        int clippingSetFaktor = clippingSets_.size() <=1 ? 1 : static_cast<int>(clippingSets_.size());
        int yLength = (static_cast<int>(clippingSetFaktor)) * (2 * controlElementRadius_.get() + 10) + 2*controlElementRadius_.get() + 30;
        if (yLength != menuDimensions_.get().y) {
            int xLength = menuDimensions_.get().x;
            menuDimensions_.set(tgt::ivec2(xLength, yLength));
            menuDimensions_.invalidate();
        }

        TouchTableMenuWidget::updateMenuCoordinates();
    }

    void TouchTableClippingWidget::updateComponents() {

        //update boundingbox control element
        contrelemBoundingBox_.setPosition(tgt::ivec2(controlElementRadius_.get() +5, controlElementRadius_.get() +5));
        contrelemBoundingBox_.setRadius(controlElementRadius_.get());

        //update add control elements
        contrelemAdd_.setPosition(tgt::ivec2(3* controlElementRadius_.get() +10, controlElementRadius_.get() +5));
        contrelemAdd_.setRadius(controlElementRadius_.get());

        contrelemAddXY_.setPosition( contrelemAdd_.getPosition() + tgt::ivec2(controlElementRadius_.get()*2 + 5, 0));
        contrelemAddXY_.setRadius(controlElementRadius_.get());

        contrelemAddXZ_.setPosition( contrelemAddXY_.getPosition() + tgt::ivec2(controlElementRadius_.get()*2 + 5, 0));
        contrelemAddXZ_.setRadius(controlElementRadius_.get());

        contrelemAddYZ_.setPosition(contrelemAddXZ_.getPosition() + tgt::ivec2(controlElementRadius_.get()*2 + 5, 0));
        contrelemAddYZ_.setRadius(controlElementRadius_.get());

        //update clipping sets
        std::vector<ClippingSet>::iterator clippingIter;
        int posInMenu=0;
        for(clippingIter = clippingSets_.begin(); clippingIter != clippingSets_.end(); ++clippingIter){

            updateClippingSetPosition(*clippingIter, posInMenu);

            posInMenu++;
        }

    }

    void TouchTableClippingWidget::renderComponents() {

        updateClippingSetColor();

        if(exclusiveModeSet_){ //gray out control elements that connot be used
            std::vector<ClippingSet>::iterator clippingIter;
            float grayOutFactor = 0.8f;
            bool grayoutAdd = false;
            //render clipping sets
            for(clippingIter = clippingSets_.begin(); clippingIter != clippingSets_.end(); ++clippingIter){
                overlay_->renderMenuFrame(clippingIter->menu_);
                overlay_->renderControlElement(&(clippingIter->toggleClippingPlane_), menu_.getLL(),grayOutFactor);
                overlay_->renderControlElement(&(clippingIter->renderClippingPlaneGeometry_), menu_.getLL(),grayOutFactor);
                overlay_->renderControlElement(&(clippingIter->invertClippingPlane_), menu_.getLL(),grayOutFactor);
                overlay_->renderControlElement(&(clippingIter->freeClipping_), menu_.getLL(),grayOutFactor);
                overlay_->renderControlElement(&(clippingIter->deletePlane_), menu_.getLL(),grayOutFactor);
                overlay_->renderControlElement(&(clippingIter->renderManipulation_), menu_.getLL(), grayOutFactor);
                overlay_->renderSlider(&(clippingIter->slider_), menu_.getLL(), grayOutFactor);
                //if activated clipping set has already been added to the active clipping sets, the free clipping has been called by the activated clipping set and the
                //free clipping control elements must not be grayed out. Add button needs to be grayed out
                if(clippingIter->freeClipping_.getId() == freeClippingLineParams_.activeClippingSet_->freeClipping_.getId()){
                    overlay_->renderControlElement(&(clippingIter->freeClipping_), menu_.getLL());
                    grayoutAdd = true;
                }
            }
            //render bounding box control element
            overlay_->renderControlElement(&contrelemBoundingBox_, menu_.getLL(), grayOutFactor);


            //render add control element
            if(!inactiveClippingSets_.empty() && (clippingSets_.size() < maxPlanes_.get())){
                if(grayoutAdd){
                    overlay_->renderControlElement(&contrelemAdd_, menu_.getLL(), grayOutFactor);
                }
                else{
                    overlay_->renderControlElement(&contrelemAdd_, menu_.getLL());
                }
                overlay_->renderControlElement(&contrelemAddXY_, menu_.getLL(), grayOutFactor);
                overlay_->renderControlElement(&contrelemAddXZ_, menu_.getLL(), grayOutFactor);
                overlay_->renderControlElement(&contrelemAddYZ_, menu_.getLL(), grayOutFactor);
            }

        }else{ // normal rendering of all control elements
            std::vector<ClippingSet>::iterator clippingIter;
            //render clipping sets
            for(clippingIter = clippingSets_.begin(); clippingIter != clippingSets_.end(); ++clippingIter){
                overlay_->renderMenuFrame(clippingIter->menu_);
                overlay_->renderControlElement(&(clippingIter->toggleClippingPlane_), menu_.getLL());
                overlay_->renderControlElement(&(clippingIter->renderClippingPlaneGeometry_), menu_.getLL());
                overlay_->renderControlElement(&(clippingIter->invertClippingPlane_), menu_.getLL());
                overlay_->renderControlElement(&(clippingIter->freeClipping_), menu_.getLL());
                overlay_->renderControlElement(&(clippingIter->deletePlane_), menu_.getLL());
                overlay_->renderControlElement(&(clippingIter->renderManipulation_), menu_.getLL());
                overlay_->renderSlider(&(clippingIter->slider_), menu_.getLL());

            }
            //render bounding box control element
            overlay_->renderControlElement(&contrelemBoundingBox_, menu_.getLL());
            //render add control element
            if(!inactiveClippingSets_.empty() && (clippingSets_.size() < maxPlanes_.get())){
                overlay_->renderControlElement(&contrelemAdd_, menu_.getLL());
                overlay_->renderControlElement(&contrelemAddXY_, menu_.getLL());
                overlay_->renderControlElement(&contrelemAddXZ_, menu_.getLL());
                overlay_->renderControlElement(&contrelemAddYZ_, menu_.getLL());
            }
        }


        if(exclusiveModeSet_ && (freeClippingLineParams_.id_ != -1)) {
            renderClippingLine(freeClippingLineParams_.startPos_, freeClippingLineParams_.currentPos_, freeClippingLineParams_.activeClippingSet_->freeClipping_.getColorMod());
            invalidate();
        }
    }

} // namespace*/
