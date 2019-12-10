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

#include "touchtablelightsourcewidget.h"

#include "tgt/material.h"

namespace voreen {

    const std::string TouchTableLightSourceWidget::loggerCat_("voreen.touchtable.TouchTableLightSourceWidget");

    TouchTableLightSourceWidget::TouchTableLightSourceWidget()
        : TouchTableMenuWidget("lightbulb.png")
        , lightDistance_(tgt::ivec2(30,20),100,15,HORIZONTAL,25,20)
        , lightPosition_("lightSourcePosition","lightSourcePosition", tgt::vec4(1.0f))
        , renderLightSource_("renderlight", "Show Light Widget", false)
        , camera_("camera", "Camera", tgt::Camera( tgt::vec3(1.f,1.f,1.f), tgt::vec3(0.f,0.f,0.f), tgt::vec3(0.f, 1.f, 0.f)))
        , shininess_("shininess","shininess",14.f, 3.f, 28.f) ///< Property because it handles min/max comfortably
        , diffuseLight_("diffuseLight","diffuse Light",tgt::vec4(0.3f,0.3f,0.3f,1.0f))
        , specularLight_("specularLight","specular Light", tgt::vec4(0.3f,0.3f,0.3f,1.0f))
        , ambientLight_("ambientLight", "ambient Light", tgt::vec4(0.3f,0.3f,0.3f,1.0f))
        , followCamera_(25, tgt::ivec2(30,205), 0, 1.0f, tgt::vec4(0.0f), tgt::vec4(1.0f), 0, true, true)
        , invertLight_(25, tgt::ivec2(30,140), 0, 1.0f, tgt::vec4(0.0f), tgt::vec4(1.0f), 0, true, false)
        , lightElem_(25, tgt::ivec2(30, 75))
        , lightTex_(0)
    {

        quadUR_ = menu_.getUR() - controlElementRadius_.get();
        quadLL_.y = menu_.getLL().y + controlElementRadius_.get()*2;
        quadLL_.x = quadUR_.y - quadLL_.y;

        addProperty(diffuseLight_);
        addProperty(specularLight_);
        addProperty(ambientLight_);

        addProperty(camera_);
        addProperty(lightPosition_);

        addProperty(renderLightSource_);
        renderLightSource_.onChange(MemberFunctionCallback<TouchTableLightSourceWidget>(this, &TouchTableLightSourceWidget::updateComponents));

        lightPosition_.setCamera(&camera_);
        lightDistance_.setIndicatorPos(0.5f);
        camera_.onChange(MemberFunctionCallback<TouchTableLightSourceWidget>(this, &TouchTableLightSourceWidget::adjustLightToCamera));
    }

    void TouchTableLightSourceWidget::handleTouchPoints(const std::deque<tgt::TouchPoint>& tp){
        for(std::deque<tgt::TouchPoint>::const_iterator currentTp = tp.begin(); currentTp != tp.end(); ++currentTp){
            tgt::vec2 tpPos = convertToMenuCoordinates(currentTp->pos());

            if(currentTp->state() == tgt::TouchPoint::TouchPointPressed){
                if(lightDistance_.checkIndicatorHit(tpPos)){
                    sliderID_ = currentTp->id();
                    continue;
                }
                if(hitsQuad(tpPos) ){
                    quadID_ = currentTp->id();
                    continue;
                }
                if(hitsControlElement(tpPos, followCamera_)){
                    followCamera_.setActive(!followCamera_.isActive());
                    lightPosition_.setFollowCam(followCamera_.isActive());
                    continue;
                }
                if(hitsControlElement(tpPos, invertLight_)){
                    invertLight_.setActive(!invertLight_.isActive());
                    lightPosition_.set(-lightPosition_.get());
                    continue;
                }
            }
            if(currentTp->state() == tgt::TouchPoint::TouchPointMoved){
                if(currentTp->id() == sliderID_){
                    lightDistance_.updateIndicatorPosition(tpPos);
                }
                if(currentTp->id() == quadID_){
                    updateLightSourcePosition(tpPos);
                }
            }
            if(currentTp->state() == tgt::TouchPoint::TouchPointReleased){
                if(currentTp->id() == sliderID_){
                    sliderID_ = -1;
                    continue;
                }
                if(currentTp->id() == quadID_){
                    quadID_ = -1;
                    continue;
                }

                if (hitsControlElement(tpPos, lightElem_))
                    renderLightSource_.set(!renderLightSource_.get());
            }
        }
        updateProperties();
        invalidate();
    }

    void TouchTableLightSourceWidget::initialize() {
        TouchTableMenuWidget::initialize();

        sliderID_ = -1;
        quadID_ = -1;

        followCameraTex_ = TexMgr.load("camera.png");
        followCamera_.setSymbolTexture(followCameraTex_);

        invertLightTex_ = TexMgr.load("invert.png");
        invertLight_.setSymbolTexture(invertLightTex_);

        lightTex_ = TexMgr.load("light.png");
        lightElem_.setSymbolTexture(lightTex_);

        menuDimensions_.set(tgt::ivec2(271,256));

        oldCameraNormal_ = tgt::normalize(camera_.get().getPosition() - camera_.get().getFocus());

        lightDistance_.setIndicatorHeight(controlElementRadius_.get() + controlElementRadius_.get()/2);
        lightDistance_.setIndicatorWidth(controlElementRadius_.get());
        lightDistance_.setBarOrigin(tgt::ivec2(controlElementRadius_.get() + controlElementRadius_.get()/2, controlElementRadius_.get()*3/4+5));
        lightDistance_.setBarWidth(controlElementRadius_.get()*3/4);
        lightDistance_.setBarLength(menuDimensions_.get().x -4 * controlElementRadius_.get());
    }

    void TouchTableLightSourceWidget::deinitialize() {
        TexMgr.dispose(followCameraTex_); followCameraTex_ = 0;
        TexMgr.dispose(invertLightTex_); invertLightTex_ = 0;

        TexMgr.dispose(lightTex_); lightTex_ = 0;
        lightElem_.setSymbolTexture(0);

        TouchTableMenuWidget::deinitialize();
    }

    void TouchTableLightSourceWidget::renderComponents(){
        float radius = (quadUR_.y - quadLL_.y)/2.0f - 5.0f;
        overlay_->renderSlider(&lightDistance_, menu_.getLL());
        overlay_->renderColoredQuad(quadLL_,quadUR_,tgt::vec4(0.f, 0.f, 0.f, 0.8f));

        overlay_->renderSphere(radius, shininess_.get(), lightPosition_.get() ,diffuseLight_.get(), specularLight_.get(), ambientLight_.get(), quadUR_- static_cast<int>(radius) - 5);
        overlay_->renderControlElement(&followCamera_, menu_.getLL());
        overlay_->renderControlElement(&invertLight_, menu_.getLL());
        overlay_->renderControlElement(&lightElem_, menu_.getLL());
    }

    void TouchTableLightSourceWidget::process(){
        TouchTableMenuWidget::process();
    }

    void TouchTableLightSourceWidget::updateComponents(){
        lightDistance_.setBarLength(menu_.getUR().x - menu_.getLL().x - (controlElementRadius_.get()*3)/2);
        quadUR_ = menu_.getUR() - controlElementRadius_.get();
        quadLL_.y = menu_.getLL().y + controlElementRadius_.get()*2;
        quadLL_.x = menu_.getLL().x + controlElementRadius_.get()*2 + 15;

        lightElem_.setActive(renderLightSource_.get());
    }

    void TouchTableLightSourceWidget::updateProperties(){
        float t = lightPosition_.getMaxDist();
        tgt::vec4 x = tgt::normalize(lightPosition_.getLightPos())*lightPosition_.getMaxDist();
        lightPosition_.setLightPos(x/10.0f+lightDistance_.getIndicatorPos()*(x-x/10.0f));
        shininess_.set(28.f - 25.0f*lightDistance_.getIndicatorPos());
    }

    bool TouchTableLightSourceWidget::hitsQuad(const tgt::ivec2& pos){
        tgt::ivec2 ll(controlElementRadius_.get()*2+10);
        tgt::ivec2 ur = ll+ 10+ quadUR_.y - quadLL_.y;
        if(pos.x <= ur.x && pos.x >= ll.x && pos.y <= ur.y && pos.y >= ll.y){
            return true;
        }else{
            return false;
        }
    }

    void TouchTableLightSourceWidget::adjustLightToCamera(){
        if(!followCamera_.isActive()){
            tgt::vec3 cameranormal = tgt::normalize(camera_.get().getPosition() - camera_.get().getFocus());
            tgt::vec3 rotaxis = tgt::cross(cameranormal,oldCameraNormal_);
            tgt::vec3 lightPos = tgt::quat::rotate(lightPosition_.getLightPos().xyz(),-acos(tgt::dot(cameranormal,oldCameraNormal_))/180.0f*tgt::PIf,rotaxis);
            lightPosition_.set(tgt::vec4(lightPos.x,lightPos.y,lightPos.z,0.0f));
            invertLight_.setActive(lightPos.z < 0.0f);
        }
    }

    void TouchTableLightSourceWidget::updateLightSourcePosition(const tgt::ivec2& pos){
        float radius = (quadUR_.y - quadLL_.y)/2.0f - 5.0f;
        tgt::ivec2 center = tgt::ivec2((int)radius) + controlElementRadius_.get()*2+10;
        //tgt::ivec2 center = menu_.getUR() - menu_.getLL() - controlElementRadius_.get() - (int)(radius) - 5;
        tgt::vec2 relPos = tgt::vec2(static_cast<float>(pos.x),static_cast<float>(pos.y));
        relPos.x -= center.x;
        relPos.y -= center.y;
        if(isMenuInverted()){
            relPos.x *= -1;
            relPos.y *= -1;
        }
        if(tgt::length(relPos) > radius){
            relPos = tgt::normalize(relPos) * (radius - 0.1f);
        }
        float z = std::sqrt(radius*radius - relPos.x*relPos.x - relPos.y*relPos.y);
        tgt::vec4 lightpos;;
        if(invertLight_.isActive()){
            lightpos = tgt::vec4(-relPos.x,-relPos.y,-z,0.0f);
        }else{
            lightpos = tgt::vec4(relPos.x,relPos.y,z,0.0f);
        }

        oldCameraNormal_ = -tgt::normalize(camera_.get().getFocus() - camera_.get().getPosition());
        tgt::vec3 zAxis = tgt::vec3(0.0f,0.0f,-1.0f);
        tgt::vec3 rotationAxis = tgt::cross(zAxis,oldCameraNormal_);
        float rotationAngle = acos(tgt::dot(oldCameraNormal_,zAxis))/180.0f*tgt::PIf;
        tgt::vec3 lightPos = tgt::quat::rotate(lightpos.xyz(), rotationAngle, rotationAxis);
        lightpos.x = lightPos.x;
        lightpos.y = lightPos.y;
        lightpos.z = lightPos.z;

         lightPosition_.setLightPos(lightpos);
    }

} // namespace
