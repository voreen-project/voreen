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

#include "voreen/qt/networkeditor/graphicitems/nwebasegraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/tooltips/tooltipbasegraphicsitem.h"
#include "voreen/qt/networkeditor/styles/nwestyle_base.h"

#include <QGraphicsDropShadowEffect>
#include <QGraphicsSceneHoverEvent>
#include <QAction>
#include <QMenu>

namespace voreen {

NWEBaseGraphicsItem::NWEBaseGraphicsItem(NetworkEditor* nwe)
    : QGraphicsObject(nullptr)
    , networkEditor_(nwe)
    , toolTipItem_(nullptr)
    , shadowEffect_(nullptr)
    , paintHasBeenInitialized_(false)
    , customContextMenu_(new QMenu())
{
    tgtAssert(nwe,"NWEBaseGraphicsItem: No NetworkEditor passed!");
}

NWEBaseGraphicsItem::~NWEBaseGraphicsItem() {
    if(customContextMenu_) {
        foreach(QAction* act, customContextMenu_->actions()) {
            delete act;
        }
        delete customContextMenu_;
    }
}

//---------------------------------------------------------------------------------------------------------------
//                  getter and setter
//---------------------------------------------------------------------------------------------------------------
NetworkEditor* NWEBaseGraphicsItem::getNetworkEditor() const{
    return networkEditor_;
}

NetworkEditorLayer NWEBaseGraphicsItem::currentLayer() const {
    return networkEditor_->getCurrentLayer();
}

NetworkEditorCursorMode NWEBaseGraphicsItem::currentCursorMode() const {
    return networkEditor_->getCurrentCursorMode();
}

NWEStyle_Base* NWEBaseGraphicsItem::currentStyle() const{
    return networkEditor_->getCurrentStyle();
}

bool NWEBaseGraphicsItem::currentToolTipMode() const{
    return networkEditor_->getCurrentToolTipMode();
}

void NWEBaseGraphicsItem::prepareGeometryChange(){
    QGraphicsItem::prepareGeometryChange();
}

void NWEBaseGraphicsItem::setToolTipGraphicsItem(ToolTipBaseGraphicsItem* tooltip) {
    toolTipItem_ = tooltip;
}

ToolTipBaseGraphicsItem* NWEBaseGraphicsItem::getToolTipGraphicsItem() const {
    return toolTipItem_;
}

bool NWEBaseGraphicsItem::isPaintInitialized() const {
    return paintHasBeenInitialized_;
}

void NWEBaseGraphicsItem::resetPaintInitialization() {
    paintHasBeenInitialized_ = false;
}

void NWEBaseGraphicsItem::enableShadows(bool enable) {
    if (enable) {
        shadowEffect_ = new QGraphicsDropShadowEffect;
        shadowEffect_->setXOffset(0);
        shadowEffect_->setYOffset(0);
        shadowEffect_->setBlurRadius(15);
        shadowEffect_->setColor(QColor(0, 0, 0, 180));
    }
    else {
        shadowEffect_ = nullptr;
    }
    setGraphicsEffect(shadowEffect_);
}

//---------------------------------------------------------------------------------------------------------------
//                  paint
//---------------------------------------------------------------------------------------------------------------
void NWEBaseGraphicsItem::paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget){
    if(!paintHasBeenInitialized_) {
        initializePaintSettings();
        layoutChildItems();
        paintHasBeenInitialized_ = true;
    }
    NWEItemSettings setting = beforePaint(painter,option);
    mainPaint(painter,option,widget,setting);
}

NWEItemSettings NWEBaseGraphicsItem::beforePaint(QPainter* painter, const QStyleOptionGraphicsItem* option){
    return currentStyle()->getItemSettings(this,option);
}

//---------------------------------------------------------------------------------------------------------------
//                  context menu
//---------------------------------------------------------------------------------------------------------------
void NWEBaseGraphicsItem::addActionsToContextMenu(QMenu* menu) {
    if(acceptHoverEvents()) //y?
        menu->addActions(customContextMenu_->actions());
}

void NWEBaseGraphicsItem::setContextMenuActions() {
    foreach(QAction* a, customContextMenu_->actions())
        delete a;
    customContextMenu_->clear();
}
//---------------------------------------------------------------------------------------------------------------
//                  events
//---------------------------------------------------------------------------------------------------------------
void NWEBaseGraphicsItem::hoverEnterEvent (QGraphicsSceneHoverEvent* event){
    if(toolTipItem_ && currentToolTipMode()){
        QPointF p = event->scenePos();
        toolTipItem_->setToolTipTimerTriggertMousePosition(p);
        toolTipItem_->startTimer();
    }
    if(shadowEffect_) {
        shadowEffect_->setBlurRadius(30);
        shadowEffect_->setColor(QColor(0, 0, 0, 200));
    }
    QGraphicsItem::hoverEnterEvent(event);
}

void NWEBaseGraphicsItem::mouseMoveEvent (QGraphicsSceneMouseEvent* event){
    if(!currentToolTipMode() && toolTipItem_) {
        if(toolTipItem_->isVisible()){
            toolTipItem_->setVisible(false);
        } else {
            toolTipItem_->stopTimer();
        }
    }
    if(isSelected() && shadowEffect_) {
        //TODO: still we have heavy ghosting effects since shadows seem to be somewhat broken.
        //scene()->update();
    }
    QGraphicsItem::mouseMoveEvent(event);
}

void NWEBaseGraphicsItem::hoverLeaveEvent (QGraphicsSceneHoverEvent* event){
    if(toolTipItem_){
        toolTipItem_->stopTimer();
        toolTipItem_->setVisible(false);
    }
    if(shadowEffect_) {
        shadowEffect_->setBlurRadius(15);
        shadowEffect_->setColor(QColor(0, 0, 0, 150));
    }
    QGraphicsItem::hoverLeaveEvent(event);
}


} // namespace voreen
