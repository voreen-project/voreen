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

#include "voreen/qt/networkeditor/graphicitems/connections/portsizelinkarrowgraphicsitem.h"
#include "voreen/qt/networkeditor/editor_settings.h"

#include "tgt/assert.h"

//style
#include "voreen/qt/networkeditor/styles/nwestyle_base.h"

//core
#include "voreen/core/ports/renderport.h"

//gi
#include "voreen/qt/networkeditor/graphicitems/core/portgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/propertygraphicsitem.h"

#include <QPainterPath>
#include <QGraphicsSceneMouseEvent>


namespace voreen {

PortSizeLinkArrowGraphicsItem::PortSizeLinkArrowGraphicsItem(PropertyGraphicsItem* sourceItem, PortGraphicsItem* sourcePort)
    : ConnectionBaseGraphicsItem(sourceItem, 0)
    , sourcePort_(sourcePort), destinationPort_(0)
    , clickedPoint_(), arrowDraged_(false), dragedArrowMoved_(false), currentArrowInitialized_(false)
{
    tgtAssert((sourceItem->getProperty()->getClassName() == "RenderSizeOriginProperty") || (sourceItem->getProperty()->getClassName() == "RenderSizeReceiveProperty"),"wrong property class");
    tgtAssert(sourcePort, "no source port");
    destinationArrowHeadDirection_ = ArrowHeadDirectionDown;

    setParent(sourcePort);

    setFlag(ItemIsSelectable);
    setZValue(ZValuesPortSizeLinkArrowGraphicsItem);
    setAcceptHoverEvents(true);
    setFlag(ItemIsPanel,true);

    sourcePort->scene()->addItem(this);

    updateNWELayerAndCursor();
}

PortSizeLinkArrowGraphicsItem::~PortSizeLinkArrowGraphicsItem() {
}

//---------------------------------------------------------------------------------------------------------------
//                  getter and setter
//---------------------------------------------------------------------------------------------------------------
QPointF PortSizeLinkArrowGraphicsItem::getSourcePoint() const {
    QPointF s = mapFromItem(sourcePort_,0,0);

    s.rx() += sourcePort_->boundingRect().width() / 2.0;
    s.ry() += sourcePort_->boundingRect().height() / 2.0;

    return s;
}

QPointF PortSizeLinkArrowGraphicsItem::getDestinationPoint() const {
    if (!destinationItem_ || !destinationPort_)
        return destinationPoint_;

    QPointF d = mapFromItem(destinationPort_,0,0);

    d.rx() += destinationPort_->boundingRect().width() / 2.0;
    d.ry() += destinationPort_->boundingRect().height() / 2.0;

    return d;
}

PropertyGraphicsItem* PortSizeLinkArrowGraphicsItem::getSourceItem() const {
    PropertyGraphicsItem* i = static_cast<PropertyGraphicsItem*>(sourceItem_);
    tgtAssert(i, "sourceItem was of wrong type");
    return i;
}

PropertyGraphicsItem* PortSizeLinkArrowGraphicsItem::getDestinationItem() const {
    if (destinationItem_) {
        PropertyGraphicsItem* i = static_cast<PropertyGraphicsItem*>(destinationItem_);
        tgtAssert(i, "destinationItem was of wrong type");
        return i;
    }
    else
        return 0;
}

PortGraphicsItem* PortSizeLinkArrowGraphicsItem::getSourcePort() const {
    return sourcePort_;
}

PortGraphicsItem* PortSizeLinkArrowGraphicsItem::getDestinationPort() const{
    return destinationPort_;
}

void PortSizeLinkArrowGraphicsItem::setDestinationPort(PortGraphicsItem* dst) {
    destinationPort_ = dst;
}

//---------------------------------------------------------------------------------------------------------------
//                  nwebasegraphicsitem functions
//---------------------------------------------------------------------------------------------------------------
QRectF PortSizeLinkArrowGraphicsItem::boundingRect() const {
    return currentStyle()->PortSizeLinkArrowGI_boundingRect(this);
}

QPainterPath PortSizeLinkArrowGraphicsItem::shape() const {
    return currentStyle()->PortSizeLinkArrowGI_shape(this);
}

void PortSizeLinkArrowGraphicsItem::initializePaintSettings() {
    currentStyle()->PortSizeLinkArrowGI_initializePaintSettings(this);
}

void PortSizeLinkArrowGraphicsItem::mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    currentStyle()->PortSizeLinkArrowGI_paint(this, painter, option, widget, setting);
}

void PortSizeLinkArrowGraphicsItem::updateNWELayerAndCursor() {
    switch(currentCursorMode()) {
    case NetworkEditorCursorSelectMode:
        switch (currentLayer()) {
        case NetworkEditorLayerDataFlow:
        case NetworkEditorLayerGeneralLinking:
        case NetworkEditorLayerCameraLinking:
            setVisible(false);
            break;
        case NetworkEditorLayerPortSizeLinking:
            setVisible(true);
            setFlag(ItemIsSelectable, true);
            setAcceptHoverEvents(true);
            break;
        default:
            tgtAssert(false, "shouldn't get here");
            break;
        }
        break;
    case NetworkEditorCursorMoveMode:
        switch (currentLayer()) {
        case NetworkEditorLayerDataFlow:
        case NetworkEditorLayerGeneralLinking:
        case NetworkEditorLayerCameraLinking:
            setVisible(false);
            break;
        case NetworkEditorLayerPortSizeLinking:
            setVisible(true);
            setFlag(ItemIsSelectable, false);
            setAcceptHoverEvents(false);
            break;
        default:
            tgtAssert(false, "shouldn't get here");
            break;
        }
        break;
    default:
        tgtAssert(false,"should not get here");
        break;
    }
}

//---------------------------------------------------------------------------------------------------------------
//                  events
//---------------------------------------------------------------------------------------------------------------
void PortSizeLinkArrowGraphicsItem::mousePressEvent(QGraphicsSceneMouseEvent* event) {
    if(event->buttons() == Qt::LeftButton){
        clickedPoint_ = event->scenePos();
        arrowDraged_ = true;
        currentArrowInitialized_ = false;
        dragedArrowMoved_ = false;
    }
    QGraphicsItem::mousePressEvent(event);
}

void PortSizeLinkArrowGraphicsItem::mouseMoveEvent(QGraphicsSceneMouseEvent* event) {
    if(arrowDraged_){
        if (!dragedArrowMoved_) {
            QPointF dis = event->scenePos()- clickedPoint_;
            if (sqrt(dis.x()*dis.x() + dis.y()+dis.y()) > paEventRadius){
                dragedArrowMoved_ = true;
            }
        }

        if (dragedArrowMoved_) {
            if(!currentArrowInitialized_){
                getSourcePort()->setCurrentPortSizeLinkArrow(this);
                currentArrowInitialized_ = true;
            }
            getSourcePort()->mouseMoveEvent(event);
        }
    } else
        QGraphicsItem::mouseMoveEvent(event);
}

void PortSizeLinkArrowGraphicsItem::mouseReleaseEvent(QGraphicsSceneMouseEvent* event) {
    if(dragedArrowMoved_)
        getSourcePort()->mouseReleaseEvent(event);
    else
        QGraphicsItem::mouseReleaseEvent(event);
}

} // namespace
