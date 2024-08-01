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

#include "voreen/qt/networkeditor/graphicitems/connections/propertylinkarrowgraphicsitem.h"
#include "voreen/qt/networkeditor/styles/nwestyle_base.h"

//core

//gi
#include "voreen/qt/networkeditor/graphicitems/core/propertygraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/portownergraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/portgraphicsitem.h"

#include <QGraphicsSceneEvent>

namespace voreen {

PropertyLinkArrowGraphicsItem::PropertyLinkArrowGraphicsItem(PropertyGraphicsItem* sourceItem, PropertyGraphicsItem* destinationItem)
    : ConnectionBaseGraphicsItem(sourceItem, destinationItem)
    , correspondingPropertyLinkArrowItem_(0)
{
    sourceItem->scene()->addItem(this);

    setZValue(ZValuesPropertyLinkArrowGraphicsItem);
    destinationArrowHeadDirection_ = ArrowHeadDirectionRight;
    setFlag(ItemIsSelectable);
    setAcceptHoverEvents(true);
    setFlag(ItemIsPanel,true);

    hide();
}

PropertyLinkArrowGraphicsItem::~PropertyLinkArrowGraphicsItem() {}

//---------------------------------------------------------------------------------------------------------------
//                  getter and setter
//---------------------------------------------------------------------------------------------------------------
PropertyLinkArrowGraphicsItem* PropertyLinkArrowGraphicsItem::getCorrespondingPropertyLinkArrowItem() const {
    return correspondingPropertyLinkArrowItem_;
}

void PropertyLinkArrowGraphicsItem::setCorrespondingPropertyLinkArrowItem(PropertyLinkArrowGraphicsItem* arrow) {
    correspondingPropertyLinkArrowItem_ = arrow;
}

PropertyGraphicsItem* PropertyLinkArrowGraphicsItem::getSourceItem() const {
    if (sourceItem_) {
        PropertyGraphicsItem* i = static_cast<PropertyGraphicsItem*>(sourceItem_);
        tgtAssert(i, "sourceItem was of wrong type");
        return i;
    } else
        return 0;
}

PropertyGraphicsItem* PropertyLinkArrowGraphicsItem::getDestinationItem() const {
    if (destinationItem_) {
        PropertyGraphicsItem* i = static_cast<PropertyGraphicsItem*>(destinationItem_);
        tgtAssert(i, "destinationItem was of wrong type");
        return i;
    }
    else
        return 0;
}

QPointF PropertyLinkArrowGraphicsItem::getSourcePoint() const {
    if(!getSourceItem())
        return mapFromScene(sourcePoint_);
    if(getSourceItem()->isVisible()){
        DockPositions dp = getDockPositions();
        if(dp == DockPositionsLeftLeft || dp == DockPositionsLeftRight)
            return mapFromScene(getSourceItem()->getLeftDockingPoint());
        else
            return mapFromScene(getSourceItem()->getRightDockingPoint());
    } else {
        if(dynamic_cast<PortOwnerGraphicsItem*>(getSourceItem()->getPropertyOwnerItem()))
            return mapFromScene(dynamic_cast<PortOwnerGraphicsItem*>(getSourceItem()->getPropertyOwnerItem())->getPropertyListButton().dockingPoint());
        else
            return mapFromScene(dynamic_cast<PortGraphicsItem*>(getSourceItem()->getPropertyOwnerItem())->getPortOwner()->getPropertyListButton().dockingPoint());
    }
}

QPointF PropertyLinkArrowGraphicsItem::getDestinationPoint() const {
    if(!getDestinationItem())
        return mapFromScene(destinationPoint_);
    if(getDestinationItem()->isVisible()){
        DockPositions dp = getDockPositions();
        if(dp == DockPositionsLeftLeft || dp == DockPositionsRightLeft)
            return mapFromScene(getDestinationItem()->getLeftDockingPoint());
        else
            return mapFromScene(getDestinationItem()->getRightDockingPoint());
    } else {
        if(dynamic_cast<PortOwnerGraphicsItem*>(getDestinationItem()->getPropertyOwnerItem()))
            return mapFromScene(dynamic_cast<PortOwnerGraphicsItem*>(getDestinationItem()->getPropertyOwnerItem())->getPropertyListButton().dockingPoint());
        else
            return mapFromScene(dynamic_cast<PortGraphicsItem*>(getDestinationItem()->getPropertyOwnerItem())->getPortOwner()->getPropertyListButton().dockingPoint());
    }
}

PropertyLinkArrowGraphicsItem::DockPositions PropertyLinkArrowGraphicsItem::getDockPositions() const {
    //scene coordinates of source and destination items
    QPointF sourcePropLeftPoint, sourcePropMiddlePoint, sourcePropRightPoint;
    QPointF destinationPropLeftPoint, destinationPropMiddlePoint, destinationPropRightPoint;
    //middle is referencing to invisible propertie using the propertylistbutton as docking point
    bool sourceMiddle = false;
    bool destinationMiddle = false;

    // handle source item
    if(!getSourceItem()) { //has no source. Connection items always have a source.
        tgtAssert(false, "Shouldn't get here!!!");
        sourcePropMiddlePoint = sourcePoint_;
        sourceMiddle = true;
    } else // arrow has a source
    if(getSourceItem()->isVisible()) { //if is visible, set property docking points
        sourcePropLeftPoint = getSourceItem()->getLeftDockingPoint();
        sourcePropRightPoint = getSourceItem()->getRightDockingPoint();
    } else { //if invisible, use property list button of port owner
        if(dynamic_cast<PortOwnerGraphicsItem*>(getSourceItem()->getPropertyOwnerItem()))
            sourcePropMiddlePoint = dynamic_cast<PortOwnerGraphicsItem*>(getSourceItem()->getPropertyOwnerItem())->getPropertyListButton().dockingPoint();
        else
            sourcePropMiddlePoint = dynamic_cast<PortGraphicsItem*>(getSourceItem()->getPropertyOwnerItem())->getPortOwner()->getPropertyListButton().dockingPoint();
        sourceMiddle = true;
    }

    //handle destination item
    if(!getDestinationItem()) { //if link is drawn and has no destination item
        destinationPropMiddlePoint = destinationPoint_;
        destinationMiddle = true;
    } else
    if(getDestinationItem()->isVisible()){ //if is visible, set property docking points
        destinationPropLeftPoint = getDestinationItem()->getLeftDockingPoint();
        destinationPropRightPoint = getDestinationItem()->getRightDockingPoint();
    } else { //if invisible, use property list button of port owner
        if(dynamic_cast<PortOwnerGraphicsItem*>(getDestinationItem()->getPropertyOwnerItem()))
            destinationPropMiddlePoint = dynamic_cast<PortOwnerGraphicsItem*>(getDestinationItem()->getPropertyOwnerItem())->getPropertyListButton().dockingPoint();
        else
            destinationPropMiddlePoint = dynamic_cast<PortGraphicsItem*>(getDestinationItem()->getPropertyOwnerItem())->getPortOwner()->getPropertyListButton().dockingPoint();
        destinationMiddle = true;
    }

    //handle connection from visible property to visible property
    if(!destinationMiddle && !sourceMiddle){
        qreal leftLeft = QLineF(sourcePropLeftPoint, destinationPropLeftPoint).length();
        qreal leftRight = QLineF(sourcePropLeftPoint, destinationPropRightPoint).length();
        qreal rightLeft = QLineF(sourcePropRightPoint, destinationPropLeftPoint).length();
        qreal rightRight = QLineF(sourcePropRightPoint, destinationPropRightPoint).length();

        if ((leftLeft <= leftRight) && (leftLeft <= rightLeft) && (leftLeft <= rightRight)) {
            return DockPositionsLeftLeft;
        }
        else if ((leftRight <= leftLeft) && (leftRight <= rightLeft) && (leftRight <= rightRight)) {
            return DockPositionsLeftRight;
        }
        else if ((rightLeft <= leftLeft) && (rightLeft <= leftRight) && (rightLeft <= rightRight)) {
            return DockPositionsRightLeft;
        }
        else {
            return DockPositionsRightRight;
        }
    } else //destination item does not exist or is invisible
        if (destinationMiddle) {
            qreal leftMiddle = QLineF(sourcePropLeftPoint, destinationPropMiddlePoint).length();
            qreal rightMiddle = QLineF(sourcePropRightPoint, destinationPropMiddlePoint).length();

            //check if we are in the PropertyLinkEditor and want to link to the same property owner
            if(getSourceItem()->getPropertyListLinkedTo() == getSourceItem()->getPropertyListItem()) {
                if (leftMiddle < rightMiddle)
                    return DockPositionsLeftLeft;
                else
                    return DockPositionsRightRight;
            } else { //link between diffrernt property owners
                if (leftMiddle < rightMiddle)
                    return DockPositionsLeftRight;
                else
                    return DockPositionsRightLeft;
            }
        } else { //source item does not exist or is invisible
            qreal middleLeft = QLineF(sourcePropMiddlePoint, destinationPropLeftPoint).length();
            qreal middleRight = QLineF(sourcePropMiddlePoint, destinationPropRightPoint).length();

            if (middleLeft < middleRight)
                return DockPositionsRightLeft;
            else
                return DockPositionsLeftRight;
        }

}

//---------------------------------------------------------------------------------------------------------------
//                  nwebasegraphicsitem functions
//---------------------------------------------------------------------------------------------------------------
QRectF PropertyLinkArrowGraphicsItem::boundingRect() const {
    return currentStyle()->PropertyLinkArrowGI_boundingRect(this);
}

QPainterPath PropertyLinkArrowGraphicsItem::shape() const {
    return currentStyle()->PropertyLinkArrowGI_shape(this);
}

void PropertyLinkArrowGraphicsItem::initializePaintSettings(){
    currentStyle()->PropertyLinkArrowGI_initializePaintSettings(this);
}

void PropertyLinkArrowGraphicsItem::mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    currentStyle()->PropertyLinkArrowGI_paint(this, painter, option, widget, setting);
}

//---------------------------------------------------------------------------------------------------------------
//                  events
//---------------------------------------------------------------------------------------------------------------
void PropertyLinkArrowGraphicsItem::mousePressEvent(QGraphicsSceneMouseEvent* event) {
    ConnectionBaseGraphicsItem::mousePressEvent(event);
}

void PropertyLinkArrowGraphicsItem::mouseMoveEvent(QGraphicsSceneMouseEvent* event) {
    ConnectionBaseGraphicsItem::mouseMoveEvent(event);
}

void PropertyLinkArrowGraphicsItem::mouseReleaseEvent(QGraphicsSceneMouseEvent* event) {
    ConnectionBaseGraphicsItem::mouseReleaseEvent(event);
}

void PropertyLinkArrowGraphicsItem::hoverEnterEvent(QGraphicsSceneHoverEvent* event) {
    ConnectionBaseGraphicsItem::hoverEnterEvent(event);
    //emit hoverEnter(this); //not needed yet
    //set corresponding arrow hovered. Used in PropertyLinkDialog to synchronize the SplitViews.
    if(getCorrespondingPropertyLinkArrowItem()){
        //set zValue to fix update delay in SplitView
        getCorrespondingPropertyLinkArrowItem()->setZValue(100);
        getCorrespondingPropertyLinkArrowItem()->setIsHovered(true);
    }
}

void PropertyLinkArrowGraphicsItem::hoverLeaveEvent(QGraphicsSceneHoverEvent* event) {
    ConnectionBaseGraphicsItem::hoverLeaveEvent(event);
    //emit hoverLeave(this); //not needed yet
    //set corresponding arrow no more hovered. Used in PropertyLinkDialog to synchronize the SplitViews.
    if(getCorrespondingPropertyLinkArrowItem())
        getCorrespondingPropertyLinkArrowItem()->setIsHovered(false);
}

} // voreen namespace

