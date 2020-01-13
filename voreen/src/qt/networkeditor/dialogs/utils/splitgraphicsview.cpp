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

//header
#include "voreen/qt/networkeditor/dialogs/utils/splitgraphicsview.h"

//general
#include "voreen/qt/networkeditor/editor_settings.h"
#include "tgt/assert.h"

//core
#include "voreen/core/properties/property.h"

//graphicsitem
#include "voreen/qt/networkeditor/graphicitems/core/propertygraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/propertylistgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/propertyownergraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/connections/propertylinkarrowgraphicsitem.h"

//Qt
#include <QWheelEvent>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QGraphicsItem>

namespace voreen {

SplitGraphicsView::SplitGraphicsView(QWidget* parent)
    : QGraphicsView(parent)
    , buddyView_(0)
    , globalViewAlignment_(GVA_UNKNOWN)
    , propertyList_(0)
    , currentDragedArrow_(0)
{
    setScene(new QGraphicsScene(this));
    setMouseTracking(true);
    setBackgroundBrush(QBrush(Qt::darkGray));
    setRenderHint(QPainter::Antialiasing);
    setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setSceneRect(-50000,-50000,100000,100000);
    setViewportUpdateMode(QGraphicsView::NoViewportUpdate); //two automatic updated graphicsviews get stuck
    connect(scene(),SIGNAL(selectionChanged()),this,SLOT(updateBuddySelection()));
}

SplitGraphicsView::~SplitGraphicsView() {
}

//---------------------------------------------------------------------------------------------------------------
//                  getter and setter
//---------------------------------------------------------------------------------------------------------------
SplitGraphicsView* SplitGraphicsView::getBuddyView() const{
    return buddyView_;
}

void SplitGraphicsView::setBuddyView(SplitGraphicsView* view) {
    //test, if alignment is left or right
    tgtAssert(globalViewAlignment_ == GVA_RIGHT || globalViewAlignment_ == GVA_LEFT,"No buddy view allowed in this global view alignment!!!");
    //test, if buddy alignment fits to this alignment
    tgtAssert((globalViewAlignment_ == GVA_RIGHT && view->getGlobalViewAlignment() == GVA_LEFT) ||
            (globalViewAlignment_ == GVA_LEFT && view->getGlobalViewAlignment() == GVA_RIGHT),"Buddy view has wrong global view alignment!!!");

    buddyView_ = view;
}

SplitGraphicsView::GlobalViewAlignment SplitGraphicsView::getGlobalViewAlignment() const{
    return globalViewAlignment_;
}

void SplitGraphicsView::setGlobalViewAlignment(GlobalViewAlignment gva) {
    tgtAssert(gva != GVA_UNKNOWN, "GVA_UNKNOWN can not be set!!!");
    switch(gva) {
    case GVA_CENTER:
        setMinimumSize(400, 200);
        break;
    case GVA_LEFT:
        setMinimumSize(200, 200);
        break;
    case GVA_RIGHT:
        setMinimumSize(200, 200);
        break;
    default:
        tgtAssert(false,"shoudn't get here");
        break;
    }
    globalViewAlignment_ = gva;
}

PropertyListGraphicsItem* SplitGraphicsView::getPropertyList() const{
    return propertyList_;
}

void SplitGraphicsView::setPropertyList(PropertyListGraphicsItem* list, PropertyOwnerGraphicsItem* destOwner) {
    tgtAssert(getGlobalViewAlignment() != GVA_UNKNOWN,"GlobalViewAlignment must be set befor calling setPropertyList!!!");
    //set list
    propertyList_ = list;
    //update to new scene
    propertyList_->setScene(scene());
    //prepare list design
    list->setHeaderMode(PropertyListGraphicsItem::BIG_HEADER);
    if(list->currentLayer() == NetworkEditorLayerCameraLinking) {
        list->setPropertyVisibleMode(PropertyListGraphicsItem::SHOW_CAMERA_PROPERTIES);
    } else {
        list->setPropertyVisibleMode(PropertyListGraphicsItem::SHOW_ALL_PROPERTIES);
    }
    list->setPropertyVisibleModifier(PropertyListGraphicsItem::HIDE_NO_PROPERTIES);
    list->setLinkArrowMode(PropertyListGraphicsItem::SHOW_FRIEND_LIST_ARROWS, destOwner->getPropertyList());
    //set position
    updatePropertyListPosition(true);
    //set visible
    list->setVisible(true);
    //connect all propertylinkarrows with the hover slots
    foreach(PropertyGraphicsItem* propertyItem, list->getPropertyItems()) {
        //if(Property)
    }
}

void SplitGraphicsView::updatePropertyListPosition(bool resetVerticalPos) {
    if(resetVerticalPos)
        propertyList_->setPos(0.0,0.0);

    qreal viewportTop;
    qreal viewportBottom;
    QRectF boundingBox = propertyList_->boundingRect();

    switch(getGlobalViewAlignment()){
    case SplitGraphicsView::GVA_CENTER:
        centerOn(0.0,0.0);
        viewportTop = mapToScene(0.0,viewport()->rect().top()).y();
        viewportBottom = mapToScene(0.0,viewport()->rect().bottom()).y();
        if (viewportBottom > propertyList_->pos().y()+boundingBox.height())
            propertyList_->setPos(0.0,viewportBottom-boundingBox.height());
        if(viewportTop < propertyList_->pos().y())
            propertyList_->setPosition(-boundingBox.width()/2.0, viewportTop);
        else
            propertyList_->setPosition(-boundingBox.width()/2.0, propertyList_->pos().y());
        break;
    case SplitGraphicsView::GVA_LEFT:
        centerOn(-size().width()/2.f, 0.0);
        viewportTop = mapToScene(0.0,viewport()->rect().top()).y();
        viewportBottom = mapToScene(0.0,viewport()->rect().bottom()).y();
        if (viewportBottom > propertyList_->pos().y()+boundingBox.height())
            propertyList_->setPos(0.0,viewportBottom-boundingBox.height());
        if(viewportTop < propertyList_->pos().y())
            propertyList_->setPosition(-size().width()/2.0 - boundingBox.width()/2.0, viewportTop);
        else
            propertyList_->setPosition(-size().width()/2.0 - boundingBox.width()/2.0, propertyList_->pos().y());
        break;
    case SplitGraphicsView::GVA_RIGHT:
        centerOn(size().width()/2.f, 0.0);
        viewportTop = mapToScene(0.0,viewport()->rect().top()).y();
        viewportBottom = mapToScene(0.0,viewport()->rect().bottom()).y();
        if (viewportBottom > propertyList_->pos().y()+boundingBox.height())
            propertyList_->setPos(0.0,viewportBottom-boundingBox.height());
        if(viewportTop < propertyList_->pos().y())
            propertyList_->setPosition(size().width()/2.0 - boundingBox.width()/2.0, viewportTop);
        else
            propertyList_->setPosition(size().width()/2.0 - boundingBox.width()/2.0, propertyList_->pos().y());
        break;
    default:
        tgtAssert(false, "Should not get here");
        break;
    }
}

//---------------------------------------------------------------------------------------------------------------
//                  budy list functions
//---------------------------------------------------------------------------------------------------------------
void SplitGraphicsView::updateViewport() {
    viewport()->update();
    if(buddyView_){
        buddyView_->viewport()->update();
    }
}

void SplitGraphicsView::updateBuddySelection() {
    if(buddyView_){
        buddyView_->scene()->blockSignals(true);
        buddyView_->scene()->clearSelection();
        QList<QGraphicsItem*> list = scene()->selectedItems();
        foreach(QGraphicsItem* item, list) {
            if(item->type() == UserTypesPropertyLinkArrowGraphicsItem){
                PropertyLinkArrowGraphicsItem* arrow = dynamic_cast<PropertyLinkArrowGraphicsItem*>(item);
                if(arrow->getSourceItem()->scene() == scene()) {
                    foreach(PropertyLinkArrowGraphicsItem* pa, arrow->getDestinationItem()->getDestinationLinkList()){
                        if(pa->getSourceItem() == arrow->getSourceItem()) {
                            pa->setSelected(true);
                            pa->setZValue(ZValuesSelectedPropertyLinkArrowGraphicsItem);
                            break;
                        }
                    }
                } else {
                    foreach(PropertyLinkArrowGraphicsItem* pa, arrow->getSourceItem()->getSourceLinkList()){
                        if(pa->getDestinationItem() == arrow->getDestinationItem()) {
                            pa->setSelected(true);
                            pa->setZValue(ZValuesSelectedPropertyLinkArrowGraphicsItem);
                            break;
                        }
                    }
                }
                break;
            }
        }
        buddyView_->scene()->blockSignals(false);
    }
}

/*void SplitGraphicsView::setBuddyHoverArrow(PropertyLinkArrowGraphicsItem* pItem) {
    @note is done by the PropertyLinkArrowItems themselves.
}

void SplitGraphicsView::unsetBuddyHoverArrow(PropertyLinkArrowGraphicsItem* pItem) {
    @note is done by the PropertyLinkArrowItems themselves.
}*/

//---------------------------------------------------------------------------------------------------------------
//                  current draged arrow functions
//---------------------------------------------------------------------------------------------------------------
PropertyLinkArrowGraphicsItem* SplitGraphicsView::getCurrentDragedArrow() {
    return currentDragedArrow_;
}

void SplitGraphicsView::initCurrentDragedArrow(PropertyGraphicsItem* item, QPointF mousePos) {
    tgtAssert(!currentDragedArrow_,"currentDragedArrow is already benn initilized");
    currentDragedArrow_ = new PropertyLinkArrowGraphicsItem(item);
    currentDragedArrow_->setDestinationPoint(mousePos);
    currentDragedArrow_->setSourceArrowHeadDirection(ConnectionBaseGraphicsItem::ArrowHeadDirectionAny);
    currentDragedArrow_->setIsHovered(true);
    currentDragedArrow_->setVisible(true);
    if(buddyView_ && !buddyView_->getCurrentDragedArrow()){
        buddyView_->initCurrentDragedArrow(item,mousePos);
        buddyView_->getCurrentDragedArrow()->scene()->removeItem(buddyView_->getCurrentDragedArrow());
        buddyView_->scene()->addItem(buddyView_->getCurrentDragedArrow());
    }
}

void SplitGraphicsView::deleteCurrentDragedArrow() {
    delete currentDragedArrow_;
    currentDragedArrow_ = 0;
    if(buddyView_ && buddyView_->getCurrentDragedArrow())
        buddyView_->deleteCurrentDragedArrow();
}

//---------------------------------------------------------------------------------------------------------------
//                  events
//---------------------------------------------------------------------------------------------------------------
void SplitGraphicsView::wheelEvent(QWheelEvent* event) {
    QPointF newPos = propertyList_->pos();
    newPos.setY(newPos.y() + event->delta()/2.0);

    QRectF listBoundingRect = propertyList_->boundingRect();

    qreal visibleAreaTop = mapToScene(0,rect().top()).y();
    qreal visibleAreaBottom = mapToScene(0,rect().bottom()).y();

    if(listBoundingRect.height() <= rect().height()) return;

    if(listBoundingRect.height() <= (visibleAreaBottom - visibleAreaTop)) {// list fits in view
        if(newPos.y() < visibleAreaTop)
            newPos.setY(visibleAreaTop);
        else if(newPos.y() + listBoundingRect.height() > visibleAreaBottom)
            newPos.setY(visibleAreaBottom - listBoundingRect.height());
    } else { //list doesn't fit in view
        if(newPos.y() > visibleAreaTop)
            newPos.setY(visibleAreaTop);
        else if(newPos.y() + listBoundingRect.height() < visibleAreaBottom)
            newPos.setY(visibleAreaBottom - listBoundingRect.height());
    }

    propertyList_->setPosition(newPos);
    updateViewport();
}

void SplitGraphicsView::mousePressEvent(QMouseEvent* event) {
    QGraphicsView::mousePressEvent(event);
    if (event->buttons() & Qt::LeftButton){
       QList<QGraphicsItem*> list = scene()->items(mapToScene(event->pos()));
       foreach(QGraphicsItem* item, list){
           if(item->type() == UserTypesPropertyGraphicsItem){
               initCurrentDragedArrow(qgraphicsitem_cast<PropertyGraphicsItem*>(item), mapToScene(event->pos()));
               break;
           }
       }
    } else {
        deleteCurrentDragedArrow();
    }
    updateViewport();
}

void SplitGraphicsView::mouseMoveEvent(QMouseEvent* event) {
    if(currentDragedArrow_){
        currentDragedArrow_->setDestinationPoint(mapToScene(event->pos()));
        currentDragedArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_DEFAULT);

        if(buddyView_){
            buddyView_->getCurrentDragedArrow()->setDestinationPoint(mapToScene(event->pos()));
            buddyView_->getCurrentDragedArrow()->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_DEFAULT);
            QList<QGraphicsItem*> iList = buddyView_->scene()->items(mapToScene(event->pos()));
            foreach(QGraphicsItem* item, iList) {
                if((item->type() == UserTypesPropertyGraphicsItem)){
                    const Property* dstProp = dynamic_cast<const PropertyGraphicsItem*>(item)->getProperty();
                    if (currentDragedArrow_->getSourceItem()->getProperty()->getClassName() == dstProp->getClassName()){
                        currentDragedArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_YES);
                        buddyView_->getCurrentDragedArrow()->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_YES);
                    } else if (currentDragedArrow_->getSourceItem()->getProperty()->isLinkableWith(dstProp)){
                        currentDragedArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_MAYBE);
                        buddyView_->getCurrentDragedArrow()->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_MAYBE);
                    } else {
                        currentDragedArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_NO);
                        buddyView_->getCurrentDragedArrow()->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_NO);
                    }
                break;
                }
            }
        } else {
            QList<QGraphicsItem*> iList = scene()->items(mapToScene(event->pos()));
            foreach(QGraphicsItem* item, iList) {
                if((item != currentDragedArrow_->getSourceItem()) && (item->type() == UserTypesPropertyGraphicsItem)){
                    const Property* dstProp = dynamic_cast<const PropertyGraphicsItem*>(item)->getProperty();
                    if (currentDragedArrow_->getSourceItem()->getProperty()->getClassName() == dstProp->getClassName())
                        currentDragedArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_YES);
                    else if (currentDragedArrow_->getSourceItem()->getProperty()->isLinkableWith(dstProp))
                        currentDragedArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_MAYBE);
                    else
                        currentDragedArrow_->setColorConnectableMode(ConnectionBaseGraphicsItem::CCM_NO);
                break;
                }
            }
        }
    }
    QGraphicsView::mouseMoveEvent(event);
    updateViewport();
}

void SplitGraphicsView::mouseReleaseEvent(QMouseEvent* event) {
    if (currentDragedArrow_) {
        if(buddyView_) {
            QList<QGraphicsItem*> iList = buddyView_->scene()->items(mapToScene(event->pos()));
            foreach(QGraphicsItem* item, iList) {
                if(item->type() == UserTypesPropertyGraphicsItem){
                    PropertyGraphicsItem* dstPropItem = dynamic_cast<PropertyGraphicsItem*>(item);
                    const Property* dstProp = dstPropItem->getProperty();
                    if (currentDragedArrow_->getSourceItem()->getProperty()->isLinkableWith(dstProp))
                        emit linkProperties(currentDragedArrow_->getSourceItem(),dstPropItem);
                    break;
                }
            }
        } else {
            QList<QGraphicsItem*> iList = scene()->items(mapToScene(event->pos()));
            foreach(QGraphicsItem* item, iList) {
                if((item != currentDragedArrow_->getSourceItem()) && (item->type() == UserTypesPropertyGraphicsItem)){
                    PropertyGraphicsItem* dstPropItem = dynamic_cast<PropertyGraphicsItem*>(item);
                    const Property* dstProp = dstPropItem->getProperty();
                    if (currentDragedArrow_->getSourceItem()->getProperty()->isLinkableWith(dstProp))
                        emit linkProperties(currentDragedArrow_->getSourceItem(),dstPropItem);
                    break;
                }
            }
        }
        deleteCurrentDragedArrow();
    }
    QGraphicsView::mouseReleaseEvent(event);
    updateViewport();
}

void SplitGraphicsView::keyPressEvent(QKeyEvent* event) {
    QGraphicsView::keyPressEvent(event);
}

void SplitGraphicsView::resizeEvent(QResizeEvent* event) {
    QGraphicsView::resizeEvent(event);
    updatePropertyListPosition(false);
    updateViewport();
}

} // namespace voreen
