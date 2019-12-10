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

#include "touchtablebodypartswidget.h"
#include "tgt/texturereader.h"
#include "voreen/core/voreenapplication.h"

namespace voreen{

    const std::string TouchTableBodyPartsWidget::loggerCat_("voreen.touchtable.bodypartswidget");

    TouchTableBodyPartsWidget::TouchTableBodyPartsWidget()
        : TouchTableMenuWidgetScrollable("bodyparts.png")
        , inPort_(Port::INPORT,"BodyPartsInPort","BodyParts Inport")
        , restoreElem_(25, tgt::ivec2(10+20,175), 0, 1.0f, tgt::vec4(0.0f), tgt::vec4(0.0f), 0, false, true)
        , bodypartsMenuHeight_("BodyPartsWidgetMenuHeight","Menu Height",200,100,1000)
        , bodypartsMenuWidth_("BodyPartsWidgetMenuWidth","Menu Width",200,100,1000)
        , restoreSignal_("RestoreSignal","Restore Signal",-1,-INT_MAX,INT_MAX)
        , bodyPartsMenu_(tgt::ivec2(0,0), tgt::vec4(0.2f,0.f,0.8f,0.5f), fontProp_.get())
        , preview_(tgt::ivec2(0),0)
    {
        addProperty(bodypartsMenuHeight_);
        addProperty(bodypartsMenuWidth_);
        addProperty(restoreSignal_);

        restoreSignal_.setVisibleFlag(false);
        selected_ = 1;
        sliderPos_ = 0.0f;

        addPort(inPort_);
    }

    TouchTableBodyPartsWidget::~TouchTableBodyPartsWidget(){

    }

    void TouchTableBodyPartsWidget::handleTouchPoints(const std::deque<tgt::TouchPoint>& tp){
        std::deque<tgt::TouchPoint> tpForBodyPartsMenu;

        for(std::deque<tgt::TouchPoint>::const_iterator currentTp = tp.begin(); currentTp != tp.end(); ++currentTp){
            tgt::vec2 tpPos = convertToMenuCoordinates(currentTp->pos());

            if(currentTp->state() == tgt::TouchPoint::TouchPointReleased){
                if(hitsControlElement(tpPos, restoreElem_)){
                    for(std::vector<TouchTablePreviewIcon>::iterator bodypartsIter = bodyParts_.begin(); bodypartsIter != bodyParts_.end(); ++bodypartsIter){
                        if(bodypartsIter->isSelected()){
                            sliderPos_ = scrollableMenu_.getSlider().getIndicatorPos();
                            restoreSignal_.set(static_cast<int>(bodypartsIter-bodyParts_.begin()));
                            preview_.previewTexture_ = 0;
                            restoreSignal_.set(-1);
                            if(++bodypartsIter == bodyParts_.end()){
                                selected_--;
                            }
                            invalidate();
                            break;
                        }
                    }
                }
            }
        }
    }

    void TouchTableBodyPartsWidget::renderComponents(){
        overlay_->renderControlElement(&restoreElem_, menu_.getLL());
        int bpMenuWidth = bodyPartsMenu_.getUR().x - bodyPartsMenu_.getLL().x;
        int menuheight = menu_.getUR().y - menu_.getLL().y;

        overlay_->renderColoredQuad(preview_.posLL_,preview_.posUR_, tgt::vec4(0.f, 0.f, 0.f, 0.8f));
        if(preview_.previewTexture_){
            overlay_->renderTexturedQuad(preview_.posLL_,preview_.posUR_,preview_.previewTexture_);
        }
    }

    void TouchTableBodyPartsWidget::initialize() {
        TouchTableMenuWidget::initialize();

        restoreTex_ = TexMgr.load("undo.png");
        restoreElem_.setSymbolTexture(restoreTex_);
        defaultTex_ = TexMgr.load("nosymbol.png");

        restoreSignal_.set(-1);

        bodyParts_.clear();

        fontProp_.get()->setFontSize(14);
        menuDimensions_.set(tgt::ivec2(814,263));
        bodypartsMenuHeight_.set(175);
        bodypartsMenuWidth_.set(515);

        updateBodyPartsMenu();

        scrollableMenu_.setFont(fontProp_.get());
        scrollableMenu_.setRowHeight(40);
        scrollableMenuIsOpen_.set(true);

        int height = scrollableMenu_.getUR().y-scrollableMenu_.getLL().y;
        int width = scrollableMenu_.getUR().x - scrollableMenu_.getLL().x;
        preview_.posLL_ = tgt::ivec2(menu_.getLL().x + 20 + width, menu_.getLL().y + scrollableMenu_.getLL().y);
        preview_.posUR_ = tgt::ivec2(preview_.posLL_.x + height *14/9, preview_.posLL_.y + height);
    }

    void TouchTableBodyPartsWidget::deinitialize() {
        TexMgr.dispose(restoreTex_); restoreTex_ = 0;
        TouchTableMenuWidget::deinitialize();
    }

    void TouchTableBodyPartsWidget::updateComponents() {
        int yPos = menu_.getUR().y - menu_.getLL().y - 10 - controlElementRadius_.get();
        restoreElem_.setPosition(tgt::ivec2(10+controlElementRadius_.get(),yPos));
        updateScrollableMenuPosition();

        int height = scrollableMenu_.getUR().y-scrollableMenu_.getLL().y;
        int width = scrollableMenu_.getUR().x - scrollableMenu_.getLL().x;
        preview_.posLL_ = tgt::ivec2(menu_.getLL().x + 20 + width, menu_.getLL().y + scrollableMenu_.getLL().y);
        preview_.posUR_ = tgt::ivec2(preview_.posLL_.x + height *14/9, preview_.posLL_.y + height);
    }

    void TouchTableBodyPartsWidget::updateBodyPartsMenu(){

        std::vector<std::pair<std::string, tgt::Texture*> > menuContent;
        for(std::vector<TouchTablePreviewIcon>::iterator bodyPartsIter = bodyParts_.begin() ; bodyPartsIter != bodyParts_.end() ; ++bodyPartsIter){
            menuContent.push_back(std::make_pair(bodyPartsIter->getCaption(), bodyPartsIter->getTexture()));
        }
        scrollableMenu_.setContent(menuContent);
    }

    void TouchTableBodyPartsWidget::updateScrollableMenuPosition(){
        scrollableMenu_.setAnchor(getPosition());
        int yPos = menu_.getUR().y - menu_.getLL().y - 2*controlElementRadius_.get() - 25;
        scrollableMenu_.setLL(tgt::ivec2(10, yPos - bodypartsMenuHeight_.get()));
        scrollableMenu_.setUR(tgt::ivec2(10+bodypartsMenuWidth_.get(), yPos));

        scrollableMenu_.placeContentMenu();
        scrollableMenu_.placeSlider();

        /*scrollableMenu_.setTextColumnWidth((scrollableMenu_.getUR().x - scrollableMenu_.getLL().x) * 3/4);
        if(!bodyParts_.empty()){
            float previewRatio = static_cast<float>(bodyParts_.back().getTexture()->getDimensions().y) / static_cast<float>(bodyParts_.back().getTexture()->getDimensions().x);
            scrollableMenu_.setRowHeight(static_cast<int>(scrollableMenu_.getTextColumnWidth() * 1/4 * previewRatio));
        }*/
        if(!bodyParts_.empty()){
        float previewRatio = static_cast<float>(bodyParts_.back().getTexture()->getDimensions().x) / static_cast<float>(bodyParts_.back().getTexture()->getDimensions().y);
        scrollableMenu_.setRowHeight((scrollableMenu_.getUR().y - scrollableMenu_.getLL().y - 20)/2);
        scrollableMenu_.setTextColumnWidth((scrollableMenu_.getUR().x - scrollableMenu_.getLL().x) - scrollableMenu_.getRowHeight()*16/9);
        }

        fontProp_.get()->setLineWidth(static_cast<float>(scrollableMenu_.getTextColumnWidth()));
        //fontProp_.invalidate();
        fontProp_.get()->setFontSize(14);
    }

    void TouchTableBodyPartsWidget::handleScrollableMenuSelection(std::string filename){
        bool gotSelection_ = false;
        for(std::vector<TouchTablePreviewIcon>::iterator bpIter = bodyParts_.begin() ; bpIter != bodyParts_.end(); ++bpIter){
            if(bpIter->getCaption() == filename){
                gotSelection_ = true;
                selected_ = static_cast<int>(bpIter - bodyParts_.begin());
                //select new hit snapshot and show it in bigger preview
                bpIter->setIsSelected(true);
                preview_.previewTexture_= bpIter->getTexture();
                invalidate();
            }else{
                //deselect current snapshot in Preview if necessary
                bpIter->setIsSelected(false);
            }
        }
        if(!gotSelection_){
            preview_.previewTexture_ = 0;
        }
        invalidate();
    }

    void TouchTableBodyPartsWidget::process(){
        TouchTableMenuWidgetScrollable::process();
        if(inPort_.hasChanged()){
            bodyParts_.clear();
            std::string selectedDescription = "no Port Data";
            if(inPort_.hasData()){
                for(size_t i = 0; i < inPort_.getData()->size(); ++i){
                    bodyParts_.push_back(TouchTablePreviewIcon(tgt::ivec2(0,0), 0, 0, inPort_.getData()->at(i).first, "", inPort_.getData()->at(i).second, false));
                }
                if(selected_ >= 0 && selected_ < bodyParts_.size()){
                    selectedDescription = bodyParts_[selected_].getCaption();
                }else{
                    selectedDescription = bodyParts_.back().getCaption();
                }
                sliderPos_ = std::min(1.0f,std::max(0.0f,static_cast<float>(selected_-1)/static_cast<float>(inPort_.getData()->size())));
            }else{
                bodyParts_.push_back(TouchTablePreviewIcon(tgt::ivec2(0,0), 0, 0, defaultTex_, "", selectedDescription, true));
                sliderPos_ = 0.0f;
                selected_ = 0;
                selectedDescription = "";
            }

            updateBodyPartsMenu();
            scrollableMenu_.getSlider().setIndicatorPos(sliderPos_);
            scrollableMenu_.setContentOffset(std::max(0,std::min(selected_,static_cast<int>(bodyParts_.size()-2))));
            handleScrollableMenuSelection(selectedDescription);

        }
    }

} //namespace voreen
