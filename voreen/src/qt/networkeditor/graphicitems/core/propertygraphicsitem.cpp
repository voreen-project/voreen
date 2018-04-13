/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
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

#include "voreen/qt/networkeditor/graphicitems/core/propertygraphicsitem.h"
#include "voreen/qt/networkeditor/styles/nwestyle_base.h"

//core
#include "voreen/core/properties/property.h"
//gi
#include "voreen/qt/networkeditor/graphicitems/core/propertylistgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/propertyownergraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/connections/propertylinkarrowgraphicsitem.h"

#include <QGraphicsSceneMouseEvent>

namespace voreen {

PropertyGraphicsItem::PropertyGraphicsItem(const Property* prop, PropertyListGraphicsItem* parent)
    : NWEBaseGraphicsItem(parent->getNetworkEditor())
    , property_(prop)
    , propertyListItem_(parent)
    , propertyLabelItem_(0)
    , propertyListLinkedTo_(0)
    , showLinkArrows_(true)
    , boundingWidth_(0.f)
{
    tgtAssert(prop != 0, "passed null pointer");
    tgtAssert(parent != 0, "passed null pointer");
    //set parent
    setParent(parent);
    setParentItem(parent);
    //set flags
    setFlag(ItemIsSelectable);  //needed?
    setFlag(ItemIsPanel);       //needed?
    setAcceptHoverEvents(true);
    setZValue(ZValuesPropertyGraphicsItem);
    //set label
    propertyLabelItem_ = new QGraphicsTextItem(this);
    propertyLabelItem_->setParentItem(this);
    //because of the dynamic behavior of the PropertyList, we have to set label text now
    initializePaintSettings();
    paintHasBeenInitialized_ = true;
}

PropertyGraphicsItem::~PropertyGraphicsItem() {
    deleteChildItems();
}

//---------------------------------------------------------------------------------------------------------------
//                  nwebasegraphicsitem functions
//---------------------------------------------------------------------------------------------------------------
QRectF PropertyGraphicsItem::boundingRect() const {
    return currentStyle()->PropertyGI_boundingRect(this);
}

QPainterPath PropertyGraphicsItem::shape() const{
    return currentStyle()->PropertyGI_shape(this);
}

void PropertyGraphicsItem::initializePaintSettings() {
    currentStyle()->PropertyGI_initializePaintSettings(this);
    resetSize(propertyLabelItem_->boundingRect().width());
    getPropertyListItem()->findNewListWidth();
}

void PropertyGraphicsItem::mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    currentStyle()->PropertyGI_paint(this, painter, option, widget,setting);
}

void PropertyGraphicsItem::deleteChildItems() {
    foreach (PropertyLinkArrowGraphicsItem* pl, sourceLinkList_){
        removeGraphicalLink(pl->getDestinationItem());
    }
    foreach (PropertyLinkArrowGraphicsItem* pl, destinationLinkList_){
        pl->getSourceItem()->removeGraphicalLink(this);
    }
}

void PropertyGraphicsItem::prepareGeometryChange() {
    NWEBaseGraphicsItem::prepareGeometryChange();
    foreach (PropertyLinkArrowGraphicsItem* pl, sourceLinkList_){
        pl->prepareGeometryChange();
    }
    foreach (PropertyLinkArrowGraphicsItem* pl, destinationLinkList_){
        pl->prepareGeometryChange();
    }
}

//---------------------------------------------------------------------------------------------------------------
//                  getter and setter
//---------------------------------------------------------------------------------------------------------------
const Property* PropertyGraphicsItem::getProperty() const {
    return property_;
}

PropertyListGraphicsItem* PropertyGraphicsItem::getPropertyListItem() const {
    return propertyListItem_;
}

QGraphicsTextItem* PropertyGraphicsItem::getPropertyLabelItem() {
    return propertyLabelItem_;
}

const QGraphicsTextItem* PropertyGraphicsItem::getPropertyLabelItem() const{
    return propertyLabelItem_;
}

PropertyOwnerGraphicsItem* PropertyGraphicsItem::getPropertyOwnerItem() const {
    return propertyListItem_->getParentPropertyOwnerItem();
}

qreal PropertyGraphicsItem::getBoundingWidth() const {
    return boundingWidth_;
}

void PropertyGraphicsItem::resetSize(qreal width) {
    boundingWidth_ = width;
}

QList<PropertyLinkArrowGraphicsItem*> PropertyGraphicsItem::getSourceLinkList() {
    return sourceLinkList_;
}

QList<PropertyLinkArrowGraphicsItem*> PropertyGraphicsItem::getDestinationLinkList() {
    return destinationLinkList_;
}

QPointF PropertyGraphicsItem::getLeftDockingPoint() const {
    return mapToScene(QPointF(boundingRect().x(), boundingRect().height() / 2.f));
}

QPointF PropertyGraphicsItem::getRightDockingPoint() const {
    return mapToScene(QPointF(boundingRect().width(), boundingRect().height() / 2.f));
}

void PropertyGraphicsItem::setShowLinkArrows(bool b) {
    showLinkArrows_ = b;
}

PropertyListGraphicsItem* PropertyGraphicsItem::getPropertyListLinkedTo() const {
    return propertyListLinkedTo_;
}

void PropertyGraphicsItem::setPropertyListLinkedTo(PropertyListGraphicsItem* list) {
    propertyListLinkedTo_ = list;
}

void PropertyGraphicsItem::setVisible(bool visible){
    QGraphicsItem::setVisible(visible);
    if(!showLinkArrows_){
        foreach(PropertyLinkArrowGraphicsItem* item, sourceLinkList_) {
            item->setVisible(false);
        }
        foreach (PropertyLinkArrowGraphicsItem* item, destinationLinkList_) {
            item->setVisible(false);
        }
    } else
    if(visible) {
        if(!propertyListLinkedTo_){
            foreach(PropertyLinkArrowGraphicsItem* item, sourceLinkList_) {
                item->setVisible(true);
            }
            foreach (PropertyLinkArrowGraphicsItem* item, destinationLinkList_) {
                item->setVisible(true);
            }
        } else {
            QList<PropertyGraphicsItem*> list = propertyListLinkedTo_->getAllPropertyItems();
            foreach(PropertyLinkArrowGraphicsItem* item, sourceLinkList_) {
                if(list.contains(item->getDestinationItem()))
                    item->setVisible(true);
                else
                    item->setVisible(false);
            }
            foreach (PropertyLinkArrowGraphicsItem* item, destinationLinkList_) {
                if(list.contains(item->getSourceItem()))
                    item->setVisible(true);
                else
                    item->setVisible(false);
            }
        }
    } else {
        foreach(PropertyLinkArrowGraphicsItem* item, sourceLinkList_) {
            item->setVisible(false);
        }
        foreach (PropertyLinkArrowGraphicsItem* item, destinationLinkList_) {
            item->setVisible(false);
        }
        if(!propertyListLinkedTo_){
            foreach(PropertyLinkArrowGraphicsItem* item, sourceLinkList_) {
                if(item->getDestinationItem()->isVisible())
                    item->getDestinationItem()->setVisible(true);
            }
            foreach (PropertyLinkArrowGraphicsItem* item, destinationLinkList_) {
                if(item->getSourceItem()->isVisible())
                    item->getSourceItem()->setVisible(true);
            }
        } else {

        }
    }
}

//---------------------------------------------------------------------------------------------------------------
//                  linking functions
//---------------------------------------------------------------------------------------------------------------
PropertyLinkArrowGraphicsItem* PropertyGraphicsItem::addGraphicalLink(PropertyGraphicsItem* propItem) {
    //if the link already exists, return the existing one
    foreach(PropertyLinkArrowGraphicsItem* arrow, getSourceLinkList()){
        if(arrow->getDestinationItem() == propItem)
            return arrow;
    }
    //create and add link
    PropertyLinkArrowGraphicsItem* link = new PropertyLinkArrowGraphicsItem(this);
    link->setDestinationItem(propItem);
    sourceLinkList_.append(link);
    //make links bidirectional if the other link already exists
    foreach(PropertyLinkArrowGraphicsItem* item, propItem->getSourceLinkList()){
        if(item->getDestinationItem() == this){
            item->setSourceArrowHeadDirection(ConnectionBaseGraphicsItem::ArrowHeadDirectionAny);
            link->setSourceArrowHeadDirection(ConnectionBaseGraphicsItem::ArrowHeadDirectionAny);
            //set incomming link bidirectional
            foreach(PropertyLinkArrowGraphicsItem* item, getDestinationLinkList()){
                if(item->getSourceItem() == propItem){
                    item->setSourceArrowHeadDirection(ConnectionBaseGraphicsItem::ArrowHeadDirectionAny);
                    break;
                }
            }
            break;
        }
    }

    //add graphical link to destination item
    PropertyLinkArrowGraphicsItem* tmpArrow = propItem->addDestinationLink(this);
    //set corresponding arrow pointer
    link->setCorrespondingPropertyLinkArrowItem(tmpArrow);
    tmpArrow->setCorrespondingPropertyLinkArrowItem(link);

    setVisible(isVisible());// needed?
    return link;
}

PropertyLinkArrowGraphicsItem* PropertyGraphicsItem::addDestinationLink(PropertyGraphicsItem* propItem) {
    PropertyLinkArrowGraphicsItem* link = new PropertyLinkArrowGraphicsItem(propItem);
    //change parent and scene
    link->setParent(this);
    link->scene()->removeItem(link);
    scene()->addItem(link);
    link->setDestinationItem(this);
    destinationLinkList_.append(link);
    //make new link bidirectinal
    foreach(PropertyLinkArrowGraphicsItem* item, getSourceLinkList()){
        if(item->getDestinationItem() == propItem){
            link->setSourceArrowHeadDirection(ConnectionBaseGraphicsItem::ArrowHeadDirectionAny);
            break;
        }
    }
    return link;
}

void PropertyGraphicsItem::removeGraphicalLink(PropertyGraphicsItem* propItem) {
    //remove link
    foreach (PropertyLinkArrowGraphicsItem* linkItem, sourceLinkList_) {
        if(linkItem->getDestinationItem() == propItem){
            sourceLinkList_.removeOne(linkItem);
            propItem->removeDestinationLink(this);
            delete linkItem;
            break;
        }
    }
    //remove bidirectinal link if existed
    foreach (PropertyLinkArrowGraphicsItem* linkItem, destinationLinkList_) {
        if(linkItem->getSourceItem() == propItem){
            linkItem->setSourceArrowHeadDirection(ConnectionBaseGraphicsItem::ArrowHeadDirectionNone);
            return;
        }
    }
}

void PropertyGraphicsItem::removeDestinationLink(PropertyGraphicsItem* propItem) {
    //remove link
    foreach (PropertyLinkArrowGraphicsItem* linkItem, destinationLinkList_) {
        if(linkItem->getSourceItem() == propItem){
            destinationLinkList_.removeOne(linkItem);
            delete linkItem;
            break;
        }
    }
    //remove bidirectianl link
    foreach (PropertyLinkArrowGraphicsItem* linkItem, sourceLinkList_) {
        if(linkItem->getDestinationItem() == propItem){
            linkItem->setSourceArrowHeadDirection(ConnectionBaseGraphicsItem::ArrowHeadDirectionNone);
            return;
        }
    }
}



} // namespace
