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

#include "voreen/qt/networkeditor/styles/nwestyle_classic.h"
#include "voreen/qt/networkeditor/editor_settings.h"

#include "voreen/qt/voreenapplicationqt.h"
#include "voreen/core/utils/stringutils.h"
#include <sstream>

//gi
#include "voreen/qt/networkeditor/graphicitems/core/portgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/processorgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/propertygraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/propertylistgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/connections/portarrowgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/connections/portsizelinkarrowgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/connections/portownerlinkarrowgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/connections/propertylinkarrowgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/utils/progressbargraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/utils/propertylistbuttongraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/utils/widgettogglebuttongraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/tooltips/tooltipportgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/tooltips/tooltipprocessorgraphicsitem.h"

//core
#include "voreen/core/ports/port.h"
#include "voreen/core/ports/coprocessorport.h"
#include "voreen/core/ports/renderport.h"
#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/processors/processorwidget.h"

//qt
#include <QPainter>
#include <QPainterPath>
#include <QStyle>
#include <QStyleOption>
#include <QtSvg/QSvgRenderer>

namespace voreen{

/*********************************************************************
 *                       General Color Defines
 ********************************************************************/
//editor
const QColor NWEStyle_Base::NWEStyle_NWEButtonBackgroundColor = QColor(50,50,50,100);
const QBrush NWEStyle_Base::NWEStyle_NWEBackgroundBrush =  QBrush(Qt::darkGray);

//general graphicsitem
const QColor NWEStyle_Base::NWEStyle_SelectionColor = QColor(0,0,237,255);
const QColor NWEStyle_Base::NWEStyle_HoverColor = QColor(0,0,97,255);

const QColor NWEStyle_Base::NWEStyle_ConnectionYes = Qt::green;
const QColor NWEStyle_Base::NWEStyle_ConnectionMaybe = Qt::yellow;
const QColor NWEStyle_Base::NWEStyle_ConnectionNo = Qt::red;

QSvgRenderer NWEStyle_Base::NWEStyle_Error1Renderer;
QSvgRenderer NWEStyle_Base::NWEStyle_Error2Renderer;
const QString NWEStyle_Base::NWEStyle_Error1SVGPath = QString(":/qt/icons/yellow-warning.svg");
const QString NWEStyle_Base::NWEStyle_Error2SVGPath = QString(":/qt/icons/red-warning.svg");

//aggregation
const QColor NWEStyle_Base::NWEStyle_AggregationColor1 = QColor(50,50,50,255);
//port
//processor
const QColor NWEStyle_Base::NWEStyle_ProcessorColor1 = QColor(50,50,50,255);
//portarrow
const QColor NWEStyle_Base::NWEStyle_PortArrowColor = Qt::black;
//portownerlinkarrow
const QColor NWEStyle_Base::NWEStyle_PortOwnerLinkArrowColor = Qt::lightGray;
//propertylinkarrow
const QColor NWEStyle_Base::NWEStyle_PropertyLinkArrowColor = Qt::lightGray;
//portsizelinkarrow
const QColor NWEStyle_Base::NWEStyle_PortSizeLinkArrowColor = Qt::black;
//tool tips
const QColor NWEStyle_Base::NWEStyle_ToolTipBackgroundColor = Qt::white;
//text boxes
const QColor NWEStyle_Base::NWEStyle_TextBoxBaseMainColor = QColor(255,234,168,255);

/*********************************************************************
 *                       General Functions
 ********************************************************************/
NWEStyle_Classic::NWEStyle_Classic(NetworkEditor* networkeditor)
    : NWEStyle_Base(networkeditor)
    {}

NWEStyle_Classic::~NWEStyle_Classic()
    {}

QPolygonF createArrowHead(ConnectionBaseGraphicsItem::ArrowHeadDirection direction, const QPointF& basePoint) {
    double angle = -tgt::PI / 2.0;

    QPointF arrowP1;
    QPointF arrowP2;
    switch (direction) {
    case ConnectionBaseGraphicsItem::ArrowHeadDirectionNone:
        return QPolygonF();
    case ConnectionBaseGraphicsItem::ArrowHeadDirectionDown:
        arrowP1 = basePoint + QPointF(sin(angle - tgt::PI / 3.0) * cbArrowHeadSize          , cos(angle - tgt::PI / 3.0) * cbArrowHeadSize);
        arrowP2 = basePoint + QPointF(sin(angle - tgt::PI + tgt::PI / 3.0) * cbArrowHeadSize, cos(angle - tgt::PI + tgt::PI / 3.0) * cbArrowHeadSize);
        break;
    case ConnectionBaseGraphicsItem::ArrowHeadDirectionUp:
        arrowP1 = basePoint - QPointF(sin(angle - tgt::PI / 3.0) * cbArrowHeadSize          , cos(angle - tgt::PI / 3.0) * cbArrowHeadSize);
        arrowP2 = basePoint - QPointF(sin(angle - tgt::PI + tgt::PI / 3.0) * cbArrowHeadSize, cos(angle - tgt::PI + tgt::PI / 3.0) * cbArrowHeadSize);
        break;
    case ConnectionBaseGraphicsItem::ArrowHeadDirectionRight:
        arrowP1 = basePoint + QPointF(cos(angle - tgt::PI / 3.0) * cbArrowHeadSize          , sin(angle - tgt::PI / 3.0) * cbArrowHeadSize);
        arrowP2 = basePoint + QPointF(cos(angle - tgt::PI + tgt::PI / 3.0) * cbArrowHeadSize, sin(angle - tgt::PI + tgt::PI / 3.0) * cbArrowHeadSize);
        break;
    case ConnectionBaseGraphicsItem::ArrowHeadDirectionLeft:
        arrowP1 = basePoint - QPointF(cos(angle - tgt::PI / 3.0) * cbArrowHeadSize          , sin(angle - tgt::PI / 3.0) * cbArrowHeadSize);
        arrowP2 = basePoint - QPointF(cos(angle - tgt::PI + tgt::PI / 3.0) * cbArrowHeadSize, sin(angle - tgt::PI + tgt::PI / 3.0) * cbArrowHeadSize);
        break;
    default:
        tgtAssert(false, "shouldn't get here");
    }
    return QPolygonF() << basePoint << arrowP1 << arrowP2 << basePoint;
}

QPolygonF createArrowHead(QPointF helpPoint, const QPointF& basePoint) {
    double angle = -tgt::PI / 2.0;

    QPointF arrowP1;
    QPointF arrowP2;
        arrowP1 = basePoint + QPointF(sin(angle - tgt::PI / 3.0) * cbArrowHeadSize          , cos(angle - tgt::PI / 3.f) * cbArrowHeadSize);
        arrowP2 = basePoint + QPointF(sin(angle - tgt::PI + tgt::PI / 3.0) * cbArrowHeadSize, cos(angle - tgt::PI + tgt::PI / 3.0) * cbArrowHeadSize);

    QPolygonF poly = QPolygonF() << basePoint << arrowP1 << arrowP2 << basePoint;

    QPointF p = helpPoint - basePoint;
    qreal sp = sqrt(p.x()*p.x()+p.y()*p.y());
    if(sp == 0)
        angle= 0;
    else if(helpPoint.x() < basePoint.x())
            angle = - acos(-p.y()/sp);
         else
            angle = acos(-p.y()/sp);

    QMatrix mat;
    mat.translate(basePoint.x(),basePoint.y());
    mat.rotate(angle * 180.f / tgt::PI);
    mat.translate(-basePoint.x(),-basePoint.y());

    return mat.map(poly);
}

    /*********************************************************************
     *                       Core Elements
     ********************************************************************/

//-------------------------------------
//      PortGraphicsItem
//-------------------------------------
QRectF NWEStyle_Classic::PortGI_boundingRect(const PortGraphicsItem* item) const {
    QRectF rect(0,0,10,10);
    if (item->getPort()->allowMultipleConnections() && !item->getPort()->isOutport()) {
        if (dynamic_cast<CoProcessorPort*>(item->getPort()))
            rect.setHeight(rect.height() * 1.75f);
        else
            rect.setWidth(rect.width() * 1.75f);
    }
    return rect;
}

QPainterPath NWEStyle_Classic::PortGI_shape(const PortGraphicsItem* item) const {
     QPainterPath path;
     path.addRect(PortGI_boundingRect(item));
     return path;
}

void NWEStyle_Classic::PortGI_initializePaintSettings(PortGraphicsItem* item) {}

void NWEStyle_Classic::PortGI_paint(PortGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    //get color and depth settings
    painter->setPen(QPen(setting.color1_, 2.0));
    painter->setBrush(setting.error_ != ES_NO_ERROR ? QBrush(Qt::red) : QBrush(Qt::lightGray));
    if (option->state & QStyle::State_Sunken)
        painter->setBrush(painter->brush().color().darker(150));
    if (option->state & QStyle::State_MouseOver)
        painter->setBrush(QBrush(setting.color1_));
    painter->setOpacity(setting.opacity_);
    //draw port
    QRectF boundRect = PortGI_boundingRect(item);
    //should be smaller to prevent trails in movement
    QRectF drawRect(boundRect.left()+1,boundRect.top()+1,boundRect.width()-2,boundRect.height()-2);
    if (item->getPort()->isLoopPort()) {
        if (item->getPort()->isOutport())
            painter->drawEllipse(drawRect);
        else {
            QPolygonF triangle;
            const QPointF& topLeftPoint = drawRect.topLeft();
            const QPointF bottomPoint = QPointF((drawRect.left() + drawRect.right()) / 2.f, drawRect.bottom());
            const QPointF& topRightPoint = drawRect.topRight();
            triangle << topLeftPoint << bottomPoint << topRightPoint;
            painter->drawConvexPolygon(triangle);
        }
    }
    else {
        painter->drawRect(drawRect);
    }
}

//-------------------------------------
//      ProcessorGraphicsItem
//-------------------------------------
QRectF NWEStyle_Classic::ProcessorGI_boundingRect(const ProcessorGraphicsItem* item) const{
    QRectF rect = item->getOrCreateNameLabel()->boundingRect();

    if (rect.width() < drawingRectMinimumWidth)
        rect.setWidth(drawingRectMinimumWidth);

    if (rect.width() < item->getMinimumWidthForPorts())
        rect.setWidth(item->getMinimumWidthForPorts());

    if (rect.height() < drawingRectMinimumHeight)
        rect.setHeight(drawingRectMinimumHeight);

    if (rect.width() < item->getMinimumHeightForPorts())
        rect.setWidth(item->getMinimumHeightForPorts());

    rect.setWidth(rect.width() + 20.f);

    return rect;
}

QPainterPath NWEStyle_Classic::ProcessorGI_shape(const ProcessorGraphicsItem* item) const{
     QPainterPath path;
     path.addRect(ProcessorGI_boundingRect(item));
     return path;
}

void NWEStyle_Classic::ProcessorGI_initializePaintSettings(ProcessorGraphicsItem* item) {
    //set text
    item->getOrCreateNameLabel()->setPlainText(item->getGuiName());
    item->getOrCreateNameLabel()->setFont(QFont("Helvetica", nwe_->getProcessorFontSizeScale()/10.f));
    item->getOrCreateNameLabel()->setDefaultTextColor(Qt::white);
}

void NWEStyle_Classic::ProcessorGI_paint(ProcessorGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    QRectF boundRect = ProcessorGI_boundingRect(item);
    //should be smaller to prevent trails in movement
    QRectF drawRect(boundRect.left()+1,boundRect.top()+1,boundRect.width()-2,boundRect.height()-2);

    //set color settings
    painter->setOpacity(setting.opacity_);
    item->setZValue(setting.zValue_);
        //gradient
    QLinearGradient gradient(0, 0, 0, drawRect.height());
    gradient.setSpread(QGradient::ReflectSpread);
    gradient.setColorAt(0.0, setting.color1_);
    gradient.setColorAt(0.4, Qt::black);
    gradient.setColorAt(0.6, Qt::black);
    gradient.setColorAt(1.0, setting.color1_);
    QBrush brush(gradient);
    painter->setBrush(brush);
    painter->setPen(QPen(QBrush(setting.color1_), 2.0));

    //draw processor
    painter->drawRect(drawRect);

    //glass highlight
    painter->setBrush(QBrush(Qt::white));
    painter->setPen(QPen(QBrush(Qt::white), 0.01));
    painter->setOpacity(0.30);
    drawRect.setHeight(drawRect.height()/2.0);
    painter->drawRect(drawRect);

    //draw error signe
    switch(setting.error_){
    case ES_ERROR_T1:{
        painter->setOpacity(1.0);
        drawRect.setHeight(drawRect.height()*2.0);
        qreal size = drawRect.height()/4.0;
        NWEStyle_Base::NWEStyle_Error1Renderer.render(painter,QRectF(drawRect.width()-size-buttonsOffsetX/2.0,drawRect.height()-size-buttonsOffsetY/2.0,size,size));
        break;}
    case ES_ERROR_T2:{
        painter->setOpacity(1.0);
        drawRect.setHeight(drawRect.height()*2.0);
        qreal size = drawRect.height()/4.0;
        NWEStyle_Base::NWEStyle_Error2Renderer.render(painter,QRectF(drawRect.width()-size-buttonsOffsetX/2.0,drawRect.height()-size-buttonsOffsetY/2.0,size,size));
        break;}
    default:
        break;
    }
}

//-------------------------------------
//      PropertyGraphicsItem
//-------------------------------------
QRectF NWEStyle_Classic::PropertyGI_boundingRect(const PropertyGraphicsItem* item) const {
    QRectF rect = item->getPropertyLabelItem()->boundingRect();
    rect.setWidth(item->getBoundingWidth());
    return rect;
}

QPainterPath NWEStyle_Classic::PropertyGI_shape(const PropertyGraphicsItem* item) const {
     QPainterPath path;
     path.addRect(PropertyGI_boundingRect(item));
     return path;
}

void NWEStyle_Classic::PropertyGI_initializePaintSettings(PropertyGraphicsItem* item) {
    //set text
    QString labelText = "";
    labelText.append(QString::fromStdString(item->getProperty()->getGuiName()));
    labelText.append("<br><\br>");
    QString type = QString::fromStdString(item->getProperty()->getTypeDescription());

    QString typeLine("<span style=\"font-size:7pt;color:#BDBDBD\"> %1 </span>");
    labelText.append(typeLine.arg(type));
    item->getPropertyLabelItem()->setHtml(labelText);
    item->resetSize(item->getPropertyLabelItem()->boundingRect().width());
}

void NWEStyle_Classic::PropertyGI_paint(PropertyGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    painter->setPen(QPen(Qt::black, 1));
    painter->setBrush(Qt::white);
    painter->drawRect(PropertyGI_boundingRect(item));
}

//-------------------------------------
//      PropertyListGraphicsItem
//-------------------------------------
QRectF NWEStyle_Classic::PropertyListGI_boundingRect(const PropertyListGraphicsItem* item) const {
    qreal height = 0.0;
    if(item->getHeaderLabelItem()->isVisible())
        height += item->getHeaderMode();
    foreach (PropertyGraphicsItem* propertyItem, item->getVisiblePropertyItems()) {
        height += propertyItem->boundingRect().height();
    }
    foreach(PropertyListGraphicsItem* listItem, item->getSubPropertyListItems()){
        height += listItem->boundingRect().height();
    }
    return QRect(0 ,0 , static_cast<int>(item->getCurrentWidth()), static_cast<int>(height));
}

QPainterPath NWEStyle_Classic::PropertyListGI_shape(const PropertyListGraphicsItem* item) const {
     QPainterPath path;
     path.addEllipse(PropertyListGI_boundingRect(item));
     return path;
}

void NWEStyle_Classic::PropertyListGI_initializePaintSettings(PropertyListGraphicsItem* item) {
    item->getHeaderLabelItem()->setFont(QFont("Helvetica", 10));
    item->getHeaderLabelItem()->setDefaultTextColor(Qt::white);
}

void NWEStyle_Classic::PropertyListGI_paint(PropertyListGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    if((item->getHeaderMode() != PropertyListGraphicsItem::NO_HEADER) && !item->getVisiblePropertyItems().empty()) {
        if(item->getHeaderMode() == PropertyListGraphicsItem::SMALL_HEADER){
            painter->setBrush(QBrush(Qt::black));
            painter->drawRect(0,0,item->getCurrentWidth(), PropertyListGraphicsItem::SMALL_HEADER);
        } else {
            painter->setBrush(QBrush(Qt::black));
            painter->drawRect(0,0,item->getCurrentWidth(), PropertyListGraphicsItem::BIG_HEADER);
            //glass highlight
            painter->setBrush(QBrush(Qt::white));
            painter->setPen(QPen(QBrush(Qt::white), 0.01));
            painter->setOpacity(0.30);
            painter->drawRect(0,0,item->getCurrentWidth(), PropertyListGraphicsItem::BIG_HEADER/2.0);
        }
    }
}

    /*********************************************************************
     *                       Util Elements
     ********************************************************************/

//-------------------------------------
//      ProgressBarGraphicsItem
//-------------------------------------
QRectF NWEStyle_Classic::ProgressBarGI_boundingRect(const ProgressBarGraphicsItem* item) const {
    return QRectF(0, 0, item->getWidth(), item->getHeight());
}

QPainterPath NWEStyle_Classic::ProgressBarGI_shape(const ProgressBarGraphicsItem* item) const {
  /*
   General Path
   c1     1________________________________________________________4      c4
         /                                                          \
        |                                                            |
        |                                                            |
         \__________________________________________________________/
   c2     2                                                        3      c3

        c1-c4: control points for splines constructing the curved edges
    */
    QRectF boundingRect_ = ProgressBarGI_boundingRect(item);
    QPainterPath path;
    const qreal cutOffValue = 5.0;
    const qreal clipRatio = 5.0;
    const QPointF p1 = QPointF(boundingRect_.x() + cutOffValue, boundingRect_.y());
    const QPointF p2 = QPointF(boundingRect_.x() + cutOffValue, boundingRect_.bottom());
    const QPointF p3 = QPointF(boundingRect_.right() - cutOffValue, boundingRect_.bottom());
    const QPointF p4 = QPointF(boundingRect_.right() - cutOffValue, boundingRect_.y());

    const QPointF c1 = QPointF(p1.x() - clipRatio, p1.y());
    const QPointF c2 = QPointF(p2.x() - clipRatio, p2.y());
    const QPointF c3 = QPointF(p3.x() + clipRatio, p3.y());
    const QPointF c4 = QPointF(p4.x() + clipRatio, p4.y());

    path.moveTo(p1);
    path.cubicTo(c1, c2, p2);
    path.lineTo(p3);
    path.cubicTo(c3, c4, p4);
    path.lineTo(p1);
    return path;
}

void NWEStyle_Classic::ProgressBarGI_initializePaintSettings(ProgressBarGraphicsItem* item) {}

void NWEStyle_Classic::ProgressBarGI_paint(ProgressBarGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {

    QRectF boundingRect_ = ProgressBarGI_boundingRect(item);
    painter->setClipPath(ProgressBarGI_shape(item));

    // paint the background
    QLinearGradient upperBackgroundGradient(boundingRect_.x(), boundingRect_.y(), boundingRect_.x(), boundingRect_.y() + boundingRect_.height() / 2.0);
    upperBackgroundGradient.setColorAt(0.0, item->getBackgroundColor1());
    upperBackgroundGradient.setColorAt(1.0, item->getBackgroundColor2());

    QBrush upperBackgroundBrush(upperBackgroundGradient);
    painter->setBrush(upperBackgroundBrush);
    painter->setPen(Qt::NoPen);
    painter->drawRect(boundingRect_.x(), boundingRect_.y(), boundingRect_.width(), boundingRect_.height() / 2.0);

    //painter->setBrush(QBrush(Qt::white));
    painter->setBrush(QBrush(Qt::black));
    painter->drawRect(boundingRect_.x(), boundingRect_.y() + boundingRect_.height() / 2.0, boundingRect_.width(), boundingRect_.height() / 2.0);

    // now the progress overlay
    QRectF progressRect = boundingRect_;
    progressRect.setWidth(boundingRect_.width() * item->getProgress());

    QLinearGradient upperForegroundGradient(progressRect.x(), progressRect.y(), progressRect.x(), progressRect.y() + progressRect.height() / 2.0);
    upperForegroundGradient.setColorAt(0.0, item->getUpperForegroundColor1());
    upperForegroundGradient.setColorAt(1.0, item->getUpperForegroundColor2());

    QBrush upperForegroundBrush(upperForegroundGradient);
    painter->setBrush(upperForegroundBrush);
    painter->drawRect(progressRect.x(), progressRect.y(), progressRect.width(), progressRect.height() / 2.0 + 1.0);

    QBrush lowerForegroundBrush(item->getLowerForegroundColor());
    painter->setBrush(lowerForegroundBrush);
    painter->drawRect(progressRect.x(), progressRect.y() + progressRect.height() / 2.0, progressRect.width(), progressRect.height() / 2.0);

}

//-------------------------------------
//      PropertyListButtonGraphicsItem
//-------------------------------------
QRectF NWEStyle_Classic::PropertyListButtonGI_boundingRect(const PropertyListButtonGraphicsItem* item) const {
    return QRect(0, 0, 10, 10);
}

QPainterPath NWEStyle_Classic::PropertyListButtonGI_shape(const PropertyListButtonGraphicsItem* item) const {
    QPainterPath path;
    path.addEllipse(PropertyListButtonGI_boundingRect(item));
    return path;
}

void NWEStyle_Classic::PropertyListButtonGI_initializePaintSettings(PropertyListButtonGraphicsItem* item) {}

void NWEStyle_Classic::PropertyListButtonGI_paint(PropertyListButtonGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    // add gloss effect
    QRectF glossArea = PropertyListButtonGI_boundingRect(item);
    glossArea.setLeft(glossArea.left() + glossArea.width() * 0.05);
    glossArea.setHeight(glossArea.height() * 0.5);
    glossArea.setWidth(glossArea.width() * 0.9);
    painter->setPen(QPen(QColor(0,0,0,0), 0));
    painter->setBrush(QColor(255,255,255,128));
    painter->drawEllipse(glossArea);

    painter->setPen(QPen(Qt::black, 1));
    painter->setBrush(QColor(255,255,255,128));

    if (option->state & QStyle::State_Sunken)
        painter->setBrush(QColor(170,100,100,128));
    if (option->state & QStyle::State_MouseOver)
        painter->setBrush(QColor(150,100,100,128));

    painter->drawEllipse(PropertyListButtonGI_boundingRect(item));
}

//-------------------------------------
//      WidgetToggleButtonGraphicsItem
//-------------------------------------
QRectF NWEStyle_Classic::WidgetToggleButtonGI_boundingRect(const WidgetToggleButtonGraphicsItem* item) const {
    return QRect(0, 0, 12, 8);
}

QPainterPath NWEStyle_Classic::WidgetToggleButtonGI_shape(const WidgetToggleButtonGraphicsItem* item) const {
    QPainterPath path;
    path.addRect(WidgetToggleButtonGI_boundingRect(item));
    return path;
}

void NWEStyle_Classic::WidgetToggleButtonGI_initializePaintSettings(WidgetToggleButtonGraphicsItem* item) {}

void NWEStyle_Classic::WidgetToggleButtonGI_paint(WidgetToggleButtonGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    painter->setPen(QPen(Qt::black, 3.5));
    painter->drawRect(WidgetToggleButtonGI_boundingRect(item));

    if (option->state & QStyle::State_MouseOver)
        painter->setPen(QPen(QColor(Qt::darkGray).darker(), 2));
    else
        painter->setPen(QPen(Qt::darkGray, 2));

    int visibleWidgetNum = 0;
    foreach (ProcessorWidget* widget, item->getWidgets()) {
        if (widget->isVisible())
            visibleWidgetNum++;
    }

    if (visibleWidgetNum > 0) {
        if (visibleWidgetNum == item->getWidgets().size())
            painter->setBrush(QColor(160, 241, 164, 175));
        else
            painter->setBrush(QColor(254, 196, 79, 175));
    }

    painter->drawRect(WidgetToggleButtonGI_boundingRect(item));
}

    /*********************************************************************
     *                       Connection Elements
     ********************************************************************/

//-------------------------------------
//      port arrow
//-------------------------------------
QRectF NWEStyle_Classic::PortArrowGI_boundingRect(const PortArrowGraphicsItem* item) const {
    QRectF rect = PortArrowGI_shape(item).boundingRect();
    // add a few extra pixels for the arrow and the pen
    if(item->getSourceItem()->getPort()->isLoopPort())
        return rect;

    qreal extra = (1 + cbArrowHeadSize) / 2.0;

    if (dynamic_cast<CoProcessorPort*>(item->getSourceItem()->getPort()))
            rect.adjust(0, -extra, 0, extra);
        else
            rect.adjust(-extra, 0, extra, 0);

    rect.adjust(-extra, -extra, extra, extra);

    return rect;
}

QPainterPath NWEStyle_Classic::PortArrowGI_shape(const PortArrowGraphicsItem* item) const {
    QPointF s = item->getSourcePoint();
    QPointF d = item->getDestinationPoint();

    QPainterPath path(s);
    qreal defl = fabs(s.y() - d.y()) / 2.f + fabs(s.x() - d.x()) / 2.f;
    PortGraphicsItem* portItem = item->getSourceItem();

    //qreal p = portArrowHoverSize;

    if (dynamic_cast<CoProcessorPort*>(portItem->getPort())) {
        d.setX(d.x()-cbArrowHeadSize);
        if(d.x() > s.x()){
            if(d.y() < s.y()) {
                    path.moveTo(s.x()-portArrowHoverSize,s.y()-portArrowHoverSize);
                    path.cubicTo(QPointF(s.x()-portArrowHoverSize+defl/2.0, s.y()-portArrowHoverSize),
                             QPointF(d.x()-portArrowHoverSize-defl/2.0, d.y()-portArrowHoverSize),
                             QPointF(d.x()-portArrowHoverSize, d.y()-portArrowHoverSize));
                    path.lineTo(d.x()+portArrowHoverSize,d.y()-portArrowHoverSize);
                    path.lineTo(d.x()+portArrowHoverSize,d.y()+portArrowHoverSize);
                    path.cubicTo(QPointF(d.x()+portArrowHoverSize-defl/2.0, d.y()+portArrowHoverSize),
                                QPointF(s.x()+portArrowHoverSize+defl/2.0, s.y()+portArrowHoverSize),
                                QPointF(s.x()+portArrowHoverSize, s.y()+portArrowHoverSize));
                    path.lineTo(s.x()-portArrowHoverSize,s.y()+portArrowHoverSize);
                    path.lineTo(s.x()-portArrowHoverSize,s.y()-portArrowHoverSize);
                }else {
                    path.moveTo(s.x()-portArrowHoverSize,s.y()-portArrowHoverSize);
                    path.lineTo(s.x()+portArrowHoverSize,s.y()-portArrowHoverSize);
                    path.cubicTo(QPointF(s.x()+portArrowHoverSize+defl/2.0, s.y()-portArrowHoverSize),
                             QPointF(d.x()+portArrowHoverSize-defl/2.0, d.y()-portArrowHoverSize),
                             QPointF(d.x()+portArrowHoverSize, d.y()-portArrowHoverSize));
                    path.lineTo(d.x()+portArrowHoverSize,d.y()+portArrowHoverSize);
                    path.lineTo(d.x()-portArrowHoverSize,d.y()+portArrowHoverSize);
                    path.cubicTo(QPointF(d.x()-portArrowHoverSize-defl/2.0, d.y()+portArrowHoverSize),
                                QPointF(s.x()-portArrowHoverSize+defl/2.0, s.y()+portArrowHoverSize),
                                QPointF(s.x()-portArrowHoverSize, s.y()+portArrowHoverSize));
                    path.lineTo(s.x()-portArrowHoverSize,s.y()-portArrowHoverSize);
                }
            } else
               if(d.y() > s.y()) {
                    path.moveTo(s.x()-portArrowHoverSize,s.y()-portArrowHoverSize);
                    path.cubicTo(QPointF(s.x()-portArrowHoverSize+defl/2.0, s.y()-portArrowHoverSize),
                             QPointF(d.x()-portArrowHoverSize-defl/2.0, d.y()-portArrowHoverSize),
                             QPointF(d.x()-portArrowHoverSize, d.y()-portArrowHoverSize));
                    path.lineTo(d.x()+portArrowHoverSize,d.y()-portArrowHoverSize);
                    path.lineTo(d.x()+portArrowHoverSize,d.y()+portArrowHoverSize);
                    path.cubicTo(QPointF(d.x()+portArrowHoverSize-defl/2.0, d.y()+portArrowHoverSize),
                                QPointF(s.x()+portArrowHoverSize+defl/2.0, s.y()+portArrowHoverSize),
                                QPointF(s.x()+portArrowHoverSize, s.y()+portArrowHoverSize));
                    path.lineTo(s.x()-portArrowHoverSize,s.y()+portArrowHoverSize);
                    path.lineTo(s.x()-portArrowHoverSize,s.y()-portArrowHoverSize);
                }else {
                    path.moveTo(s.x()-portArrowHoverSize,s.y()-portArrowHoverSize);
                    path.lineTo(s.x()+portArrowHoverSize,s.y()-portArrowHoverSize);
                    path.cubicTo(QPointF(s.x()+portArrowHoverSize+defl/2.0, s.y()-portArrowHoverSize),
                             QPointF(d.x()+portArrowHoverSize-defl/2.0, d.y()-portArrowHoverSize),
                             QPointF(d.x()+portArrowHoverSize, d.y()-portArrowHoverSize));
                    path.lineTo(d.x()+portArrowHoverSize,d.y()+portArrowHoverSize);
                    path.lineTo(d.x()-portArrowHoverSize,d.y()+portArrowHoverSize);
                    path.cubicTo(QPointF(d.x()-portArrowHoverSize-defl/2.0, d.y()+portArrowHoverSize),
                                QPointF(s.x()-portArrowHoverSize+defl/2.0, s.y()+portArrowHoverSize),
                                QPointF(s.x()-portArrowHoverSize, s.y()+portArrowHoverSize));
                    path.lineTo(s.x()-portArrowHoverSize,s.y()-portArrowHoverSize);
                }
    }
    else {
        Port* sourcePort = item->getSourceItem()->getPort();

        if (sourcePort->isLoopPort()) {
            defl *= 0.9;
            path.moveTo(s.x()-portArrowHoverSize,s.y());
            path.cubicTo(QPointF(s.x()-portArrowHoverSize+defl, s.y()+defl+5.0),
                            QPointF(d.x()-portArrowHoverSize+defl, d.y()-defl),
                            QPointF(d.x()-portArrowHoverSize, d.y()));
            path.lineTo(d.x()+portArrowHoverSize,d.y());
            path.cubicTo(QPointF(d.x()+portArrowHoverSize+defl, d.y()-defl),
                            QPointF(s.x()+portArrowHoverSize+defl, s.y()+defl+5.0),
                            QPointF(s.x()+portArrowHoverSize, s.y()));
            path.lineTo(s.x()-portArrowHoverSize,s.y());
        }
        else if(false) { //bundeled_
            /*QList<QPointF> bundlePoints = bundle_->getBundlePoints();

            float xoff = 5.f * (bundle_->getArrowIndex(this) - bundle_->arrowList_.size() / 2);
            for(int i = 0; i < bundlePoints.size(); i++)
                bundlePoints[i].rx() += xoff;
            float distS = 0.5f * (fabs(bundlePoints.front().x() + xoff - s.x()) + fabs(bundlePoints.front().y() - s.y()));

            path.cubicTo(s + QPointF(0.f, distS / 2.f),
                         bundlePoints.front() - QPointF(0.f, distS / 2.f),
                         bundlePoints.front());

            for(int i = 0; i < bundlePoints.size() - 1; i++) {
                QPointF bS = bundlePoints.at(i);
                QPointF bD = bundlePoints.at(i + 1);
                float distB = 0.5f * (fabs(bD.x() - bS.x()) + fabs(bD.y() - bS.y()));
                path.cubicTo(
                             bS + QPointF(0, distB /2.f),
                             bD - QPointF(0, distB /2.f),
                             bD);
            }

            float distD = 0.5f * (fabs(bundlePoints.back().x() - d.x()) + fabs(bundlePoints.back().y() - d.y()));
            path.cubicTo(bundlePoints.back() + QPointF(0, distD /2.0),
                         d - QPointF(0, distD /2.0),
                         d - QPointF(0, arrowHeadSize_));
            path.lineTo(d);*/
        }
        else {
            d.setY(d.y()-cbArrowHeadSize);
            if(d.y() > s.y()){
                if(d.x() < s.x()) {
                    path.moveTo(s.x()-portArrowHoverSize,s.y()-portArrowHoverSize);
                    path.cubicTo(QPointF(s.x()-portArrowHoverSize, s.y()+defl/2.0-portArrowHoverSize),
                             QPointF(d.x()-portArrowHoverSize, d.y()-defl/2.0-portArrowHoverSize),
                             QPointF(d.x()-portArrowHoverSize, d.y()-portArrowHoverSize));
                    path.lineTo(d.x()-portArrowHoverSize,d.y()+portArrowHoverSize);
                    path.lineTo(d.x()+portArrowHoverSize,d.y()+portArrowHoverSize);
                    path.cubicTo(QPointF(d.x()+portArrowHoverSize, d.y()-defl/2.0+portArrowHoverSize),
                                QPointF(s.x()+portArrowHoverSize, s.y()+defl/2.0+portArrowHoverSize),
                                QPointF(s.x()+portArrowHoverSize, s.y()+portArrowHoverSize));
                    path.lineTo(s.x()+portArrowHoverSize,s.y()-portArrowHoverSize);
                    path.lineTo(s.x()-portArrowHoverSize,s.y()-portArrowHoverSize);
                }else {
                    path.moveTo(s.x()-portArrowHoverSize,s.y()-portArrowHoverSize);
                    path.lineTo(s.x()-portArrowHoverSize,s.y()+portArrowHoverSize);
                    path.cubicTo(QPointF(s.x()-portArrowHoverSize, s.y()+defl/2.0+portArrowHoverSize),
                             QPointF(d.x()-portArrowHoverSize, d.y()-defl/2.0+portArrowHoverSize),
                             QPointF(d.x()-portArrowHoverSize, d.y()+portArrowHoverSize));
                    path.lineTo(d.x()+portArrowHoverSize,d.y()+portArrowHoverSize);
                    path.lineTo(d.x()+portArrowHoverSize,d.y()-portArrowHoverSize);
                    path.cubicTo(QPointF(d.x()+portArrowHoverSize, d.y()-defl/2.0-portArrowHoverSize),
                                QPointF(s.x()+portArrowHoverSize, s.y()+defl/2.0-portArrowHoverSize),
                                QPointF(s.x()+portArrowHoverSize, s.y()-portArrowHoverSize));
                    path.lineTo(s.x()-portArrowHoverSize,s.y()-portArrowHoverSize);
                }
            } else
                if(d.x() < s.x()) {
                    path.moveTo(s.x()-portArrowHoverSize,s.y()-portArrowHoverSize);
                    path.lineTo(s.x()-portArrowHoverSize,s.y()+portArrowHoverSize);
                    path.cubicTo(QPointF(s.x()-portArrowHoverSize, s.y()+defl/2.0+portArrowHoverSize),
                             QPointF(d.x()-portArrowHoverSize, d.y()-defl/2.0+portArrowHoverSize),
                             QPointF(d.x()-portArrowHoverSize, d.y()+portArrowHoverSize));
                    path.lineTo(d.x()+portArrowHoverSize,d.y()+portArrowHoverSize);
                    path.lineTo(d.x()+portArrowHoverSize,d.y()-portArrowHoverSize);
                    path.cubicTo(QPointF(d.x()+portArrowHoverSize, d.y()-defl/2.0-portArrowHoverSize),
                                QPointF(s.x()+portArrowHoverSize, s.y()+defl/2.0-portArrowHoverSize),
                                QPointF(s.x()+portArrowHoverSize, s.y()-portArrowHoverSize));
                    path.lineTo(s.x()-portArrowHoverSize,s.y()-portArrowHoverSize);
                }else {
                    path.moveTo(s.x()-portArrowHoverSize,s.y()-portArrowHoverSize);
                    path.cubicTo(QPointF(s.x()-portArrowHoverSize, s.y()+defl/2.0-portArrowHoverSize),
                             QPointF(d.x()-portArrowHoverSize, d.y()-defl/2.0-portArrowHoverSize),
                             QPointF(d.x()-portArrowHoverSize, d.y()-portArrowHoverSize));
                    path.lineTo(d.x()-portArrowHoverSize,d.y()+portArrowHoverSize);
                    path.lineTo(d.x()+portArrowHoverSize,d.y()+portArrowHoverSize);
                    path.cubicTo(QPointF(d.x()+portArrowHoverSize, d.y()-defl/2.0+portArrowHoverSize),
                                QPointF(s.x()+portArrowHoverSize, s.y()+defl/2.0+portArrowHoverSize),
                                QPointF(s.x()+portArrowHoverSize, s.y()+portArrowHoverSize));
                    path.lineTo(s.x()+portArrowHoverSize,s.y()-portArrowHoverSize);
                    path.lineTo(s.x()-portArrowHoverSize,s.y()-portArrowHoverSize);
                }
        }
    }
    //adding arrowhead
    if (item->getDestinationArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone) {
        QPolygonF poly = createArrowHead(item->getDestinationArrowHeadDirection(), item->getDestinationPoint());
        path.moveTo(poly.at(0));
        path.addPolygon(poly);
    }

    return path;
}

void NWEStyle_Classic::PortArrowGI_initializePaintSettings(const PortArrowGraphicsItem* item) {}

void NWEStyle_Classic::PortArrowGI_paint(PortArrowGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    QPointF s = item->getSourcePoint();
    QPointF d = item->getDestinationPoint();
    qreal defl = fabs(s.y() - d.y()) / 2.f + fabs(s.x() - d.x()) / 2.f;
    QPainterPath path(s);

    //draw shape
    if (dynamic_cast<CoProcessorPort*>(item->getSourceItem()->getPort())) {
        path.cubicTo(QPointF(s.x()+defl/2.0, s.y()),
                     QPointF(d.x()-defl/2.0-cbArrowHeadSize, d.y()),
                     QPointF(d.x()-cbArrowHeadSize,d.y()));
    }
    else {
        Port* sourcePort = item->getSourceItem()->getPort();

        if (sourcePort->isLoopPort()) {
            defl *= 0.9;
            path.cubicTo(QPointF(s.x()+defl, s.y()+ defl + 5.0),
                         QPointF(d.x()+defl, d.y()-defl),
                         d );
        }
        else if(false) { //bundeled_
            /*QList<QPointF> bundlePoints = bundle_->getBundlePoints();

            float xoff = 5.f * (bundle_->getArrowIndex(this) - bundle_->arrowList_.size() / 2);
            for(int i = 0; i < bundlePoints.size(); i++)
                bundlePoints[i].rx() += xoff;
            float distS = 0.5f * (fabs(bundlePoints.front().x() + xoff - s.x()) + fabs(bundlePoints.front().y() - s.y()));

            path.cubicTo(s + QPointF(0.f, distS / 2.f),
                         bundlePoints.front() - QPointF(0.f, distS / 2.f),
                         bundlePoints.front());

            for(int i = 0; i < bundlePoints.size() - 1; i++) {
                QPointF bS = bundlePoints.at(i);
                QPointF bD = bundlePoints.at(i + 1);
                float distB = 0.5f * (fabs(bD.x() - bS.x()) + fabs(bD.y() - bS.y()));
                path.cubicTo(
                             bS + QPointF(0, distB /2.f),
                             bD - QPointF(0, distB /2.f),
                             bD);
            }

            float distD = 0.5f * (fabs(bundlePoints.back().x() - d.x()) + fabs(bundlePoints.back().y() - d.y()));
            path.cubicTo(bundlePoints.back() + QPointF(0, distD /2.0),
                         d - QPointF(0, distD /2.0),
                         d - QPointF(0, arrowHeadSize_));
            path.lineTo(d);*/
        }
        else {
            path.cubicTo(QPointF(s.x(), s.y()+defl/2.0),
                     QPointF(d.x(), d.y()-cbArrowHeadSize-defl/2.0),
                     QPointF(d.x(), d.y()-cbArrowHeadSize));
        }
    }
    //get color and depth settings
    item->setZValue(setting.zValue_);
    painter->setOpacity(setting.opacity_);
    //draw path
    painter->setPen(QPen(setting.color1_, 2, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    painter->drawPath(path);
    //draw head
    if (item->getDestinationArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone){
        painter->setBrush(QBrush(setting.color1_));
        painter->drawPolygon(createArrowHead(item->getDestinationArrowHeadDirection(), item->getDestinationPoint()));
    }
}

//-------------------------------------
//      PortOwnerLinkArrowGraphicsItem
//-------------------------------------
QRectF NWEStyle_Classic::PortOwnerLinkArrowGI_boundingRect(const PortOwnerLinkArrowGraphicsItem* item) const {
    return PortOwnerLinkArrowGI_shape(item).boundingRect();
}

QPainterPath NWEStyle_Classic::PortOwnerLinkArrowGI_shape(const PortOwnerLinkArrowGraphicsItem* item) const {
    QPointF sourcePoint = item->getSourcePoint();
    QPointF destinationPoint = item->getDestinationPoint();

    qreal defl = fabs(sourcePoint.y() - destinationPoint.y()) / 2.0 + fabs(sourcePoint.x() - destinationPoint.x()) / 2.0;

    QPainterPath path(sourcePoint);

    if (item->getSourceItem() == item->getDestinationItem()) {
        path.cubicTo(QPointF(sourcePoint.x() - 20.f, sourcePoint.y() - 50.f),
                     QPointF(destinationPoint.x() + 20.f, destinationPoint.y() - 50.f),
                     QPointF(destinationPoint.x(), destinationPoint.y()));
        path.lineTo(destinationPoint);
    } else {
        path.lineTo(QPointF(sourcePoint.x(), sourcePoint.y() - 10.f));
        path.cubicTo(QPointF(sourcePoint.x(), sourcePoint.y() - defl / 2.f),
                     QPointF(destinationPoint.x(), destinationPoint.y() - defl / 2.f),
                     QPointF(destinationPoint.x(), destinationPoint.y() - 10.f));
        path.lineTo(destinationPoint);
    }

    return path;
}

void NWEStyle_Classic::PortOwnerLinkArrowGI_initializePaintSettings(PortOwnerLinkArrowGraphicsItem* item) {}

void NWEStyle_Classic::PortOwnerLinkArrowGI_paint(PortOwnerLinkArrowGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    QPointF sourcePoint = item->getSourcePoint();
    QPointF destinationPoint = item->getDestinationPoint();

    qreal defl = fabs(sourcePoint.y() - destinationPoint.y()) / 2.0 + fabs(sourcePoint.x() - destinationPoint.x()) / 2.0;

    QPainterPath path(sourcePoint);

    if (item->getSourceItem() == item->getDestinationItem()) {
        path.cubicTo(QPointF(sourcePoint.x() - 20.f, sourcePoint.y() - 50.f),
                     QPointF(destinationPoint.x() + 20.f, destinationPoint.y() - 50.f),
                     QPointF(destinationPoint.x(), destinationPoint.y()));
        path.lineTo(destinationPoint);
    } else {
        path.lineTo(QPointF(sourcePoint.x(), sourcePoint.y() - 10.f));
        path.cubicTo(QPointF(sourcePoint.x(), sourcePoint.y() - defl / 2.f),
                     QPointF(destinationPoint.x(), destinationPoint.y() - defl / 2.f),
                     QPointF(destinationPoint.x(), destinationPoint.y() - 10.f));
        path.lineTo(destinationPoint);
    }
    //set color and depth settings
    painter->setPen(QPen(setting.color1_, 2, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    painter->drawPath(path);
}

//-------------------------------------
//      PropertyLinkArrowrGraphicsItem
//-------------------------------------
QRectF NWEStyle_Classic::PropertyLinkArrowGI_boundingRect(const PropertyLinkArrowGraphicsItem* item) const {
    return PropertyLinkArrowGI_shape(item).boundingRect();
}

QPainterPath NWEStyle_Classic::PropertyLinkArrowGI_shape(const PropertyLinkArrowGraphicsItem* item) const {
    PropertyLinkArrowGraphicsItem::DockPositions dp = item->getDockPositions();
    QPointF s = item->getSourcePoint();
    QPointF d = item->getDestinationPoint();
    QPainterPath path(s);
    qreal defl = fabs(s.y() - d.y()) / 2.f + fabs(s.x() - d.x()) / 2.f;
    qreal p = portArrowHoverSize;

    switch (dp) {
    case PropertyLinkArrowGraphicsItem::DockPositionsLeftLeft:
        //setPos
        s.setX(s.x() - cbArrowHeadSize);
        d.setX(d.x() - cbArrowHeadSize);
        //src
        if (item->getSourceArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone) {
            QPolygonF poly = createArrowHead(ConnectionBaseGraphicsItem::ArrowHeadDirectionRight, item->getSourcePoint());
            path.addPolygon(poly);
            path.moveTo(s);
        } else {
            path.addRect(s.x(), s.y() - p, cbArrowHeadSize, 2*p);
            path.lineTo(s);
        }
        //curve
        if(s.y() >= d.y()){
            s.setY(s.y()-p);
            d.setY(d.y()+p);
            path.moveTo(s);
            path.cubicTo(QPointF(s.x()+p - defl, s.y()),
                         QPointF(d.x()+p - defl, d.y()),
                         d);
            s.setY(s.y()+2*p);
            d.setY(d.y()-2*p);
            path.lineTo(d);
            path.cubicTo(QPointF(d.x()-p - defl, d.y()),
                         QPointF(s.x()-p - defl, s.y()),
                         s);
            path.lineTo(s.x(),s.y()-2*p);
            path.moveTo(d.x(),d.y()+p);
        } else {
            s.setY(s.y()+p);
            d.setY(d.y()-p);
            path.moveTo(s);
            path.cubicTo(QPointF(s.x()+p - defl, s.y()),
                         QPointF(d.x()+p - defl, d.y()),
                         d);
            s.setY(s.y()-2*p);
            d.setY(d.y()+2*p);
            path.lineTo(d);
            path.cubicTo(QPointF(d.x()-p - defl, d.y()),
                         QPointF(s.x()-p - defl, s.y()),
                         s);
            path.lineTo(s.x(),s.y()+2*p);
            path.moveTo(d.x(),d.y()-p);
        }
        //dst
        if (item->getDestinationArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone) {
            QPolygonF poly = createArrowHead(ConnectionBaseGraphicsItem::ArrowHeadDirectionRight, item->getDestinationPoint());
            path.moveTo(poly.at(0));
            path.addPolygon(poly);
        } else {
            path.lineTo(item->getDestinationPoint());
        }
        break;
    case PropertyLinkArrowGraphicsItem::DockPositionsLeftRight:
        //setPos
        s.setX(s.x() - cbArrowHeadSize);
        d.setX(d.x() + cbArrowHeadSize);
        //src
        if (item->getSourceArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone) {
            QPolygonF poly = createArrowHead(ConnectionBaseGraphicsItem::ArrowHeadDirectionRight, item->getSourcePoint());
            path.addPolygon(poly);
            path.moveTo(s);
        } else {
            path.addRect(s.x(),s.y()-p,cbArrowHeadSize,2*p);
            path.lineTo(s);
        }
        //curve
        if(s.y() >= d.y()){
            s.setY(s.y()-p);
            d.setY(d.y()-p);
            path.moveTo(s);
            path.cubicTo(QPointF(s.x()+2*p - defl, s.y()),
                         QPointF(d.x()+2*p + defl, d.y()),
                         d);
            s.setY(s.y()+2*p);
            d.setY(d.y()+2*p);
            path.lineTo(d);
            path.cubicTo(QPointF(d.x()-2*p + defl, d.y()),
                         QPointF(s.x()-2*p - defl, s.y()),
                         s);
            path.lineTo(s.x(),s.y()-2*p);
            path.moveTo(d.x(),d.y()-p);
        } else {
            s.setY(s.y()-p);
            d.setY(d.y()-p);
            path.moveTo(s);
            path.cubicTo(QPointF(s.x()-2*p - defl, s.y()),
                         QPointF(d.x()-2*p + defl, d.y()),
                         d);
            s.setY(s.y()+2*p);
            d.setY(d.y()+2*p);
            path.lineTo(d);
            path.cubicTo(QPointF(d.x()+2*p + defl, d.y()),
                         QPointF(s.x()+2*p - defl, s.y()),
                         s);
            path.lineTo(s.x(),s.y()-2*p);
            path.moveTo(d.x(),d.y()-p);
        }
        //dst
        if (item->getDestinationArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone) {
            QPolygonF poly = createArrowHead(ConnectionBaseGraphicsItem::ArrowHeadDirectionLeft, item->getDestinationPoint());
            path.moveTo(poly.at(0));
            path.addPolygon(poly);
        } else {
            path.lineTo(item->getDestinationPoint());
        }
        break;
    case PropertyLinkArrowGraphicsItem::DockPositionsRightLeft:
        //setPos
        s.setX(s.x() + cbArrowHeadSize);
        d.setX(d.x() - cbArrowHeadSize);
        //src
        if (item->getSourceArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone) {
            QPolygonF poly = createArrowHead(ConnectionBaseGraphicsItem::ArrowHeadDirectionLeft, item->getSourcePoint());
            path.addPolygon(poly);
            path.moveTo(s);
        } else {
            path.addRect(s.x()-cbArrowHeadSize,s.y()-p,cbArrowHeadSize,2*p);
            path.moveTo(s);
        }
        //curve
        if(s.y() >= d.y()){
            s.setY(s.y()-p);
            d.setY(d.y()-p);
            path.moveTo(s);
            path.cubicTo(QPointF(s.x()-2*p + defl, s.y()),
                         QPointF(d.x()-2*p - defl, d.y()),
                         d);
            s.setY(s.y()+2*p);
            d.setY(d.y()+2*p);
            path.lineTo(d);
            path.cubicTo(QPointF(d.x()+2*p - defl, d.y()),
                         QPointF(s.x()+2*p + defl, s.y()),
                         s);
            path.lineTo(s.x(),s.y()-2*p);
            path.moveTo(d.x(),d.y()-p);
        } else {
            s.setY(s.y()-p);
            d.setY(d.y()-p);
            path.moveTo(s);
            path.cubicTo(QPointF(s.x()+2*p + defl, s.y()),
                         QPointF(d.x()+2*p - defl, d.y()),
                         d);
            s.setY(s.y()+2*p);
            d.setY(d.y()+2*p);
            path.lineTo(d);
            path.cubicTo(QPointF(d.x()-2*p - defl, d.y()),
                         QPointF(s.x()-2*p + defl, s.y()),
                         s);
            path.lineTo(s.x(),s.y()-2*p);
            path.moveTo(d.x(),d.y()-p);
        }
        //dst
        if (item->getDestinationArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone) {
            QPolygonF poly = createArrowHead(ConnectionBaseGraphicsItem::ArrowHeadDirectionRight, item->getDestinationPoint());
            path.moveTo(poly.at(0));
            path.addPolygon(poly);
        } else {
            path.lineTo(item->getDestinationPoint());
        }
        break;
    case PropertyLinkArrowGraphicsItem::DockPositionsRightRight:
        //setPos
        s.setX(s.x() + cbArrowHeadSize);
        d.setX(d.x() + cbArrowHeadSize);
        //src
        if (item->getSourceArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone) {
            QPolygonF poly = createArrowHead(ConnectionBaseGraphicsItem::ArrowHeadDirectionLeft, item->getSourcePoint());
            path.addPolygon(poly);
            path.moveTo(s);
        } else {
            path.addRect(s.x()-cbArrowHeadSize,s.y()-p,cbArrowHeadSize,2*p);
            path.moveTo(s);
        }
        //curve
       if(s.y() <= d.y()){
            s.setY(s.y()-p);
            d.setY(d.y()+p);
            path.moveTo(s);
            path.cubicTo(QPointF(s.x()+p + defl, s.y()),
                         QPointF(d.x()+p + defl, d.y()),
                         d);
            s.setY(s.y()+2*p);
            d.setY(d.y()-2*p);
            path.lineTo(d);
            path.cubicTo(QPointF(d.x()-p + defl, d.y()),
                         QPointF(s.x()-p + defl, s.y()),
                         s);
            path.lineTo(s.x(),s.y()-2*p);
            path.moveTo(d.x(),d.y()+p);
        } else {
            s.setY(s.y()+p);
            d.setY(d.y()-p);
            path.moveTo(s);
            path.cubicTo(QPointF(s.x()+p + defl, s.y()),
                         QPointF(d.x()+p + defl, d.y()),
                         d);
            s.setY(s.y()-2*p);
            d.setY(d.y()+2*p);
            path.lineTo(d);
            path.cubicTo(QPointF(d.x()-p + defl, d.y()),
                         QPointF(s.x()-p + defl, s.y()),
                         s);
            path.lineTo(s.x(),s.y()+2*p);
            path.moveTo(d.x(),d.y()-p);
        }
        //dst
        if (item->getDestinationArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone) {
            QPolygonF poly = createArrowHead(ConnectionBaseGraphicsItem::ArrowHeadDirectionLeft, item->getDestinationPoint());
            path.moveTo(poly.at(0));
            path.addPolygon(poly);
        } else {
            path.lineTo(item->getDestinationPoint());
        }
        break;
    default:
        tgtAssert(false, "shouldn't get here");
    }

    return path;
}

void NWEStyle_Classic::PropertyLinkArrowGI_initializePaintSettings(const PropertyLinkArrowGraphicsItem* item) {}

void NWEStyle_Classic::PropertyLinkArrowGI_paint(PropertyLinkArrowGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    //get color and depth settings
    item->setZValue(setting.zValue_);
    painter->setPen(QPen(setting.color1_,2));
    painter->setBrush(QBrush(setting.color1_));

    PropertyLinkArrowGraphicsItem::DockPositions dp = item->getDockPositions();
    QPointF s = item->getSourcePoint();
    QPointF d = item->getDestinationPoint();
    QPainterPath path(s);
    qreal defl = fabs(s.y() - d.y()) / 2.f + fabs(s.x() - d.x()) / 2.f;

    switch (dp) {
    case PropertyLinkArrowGraphicsItem::DockPositionsLeftLeft:
        //setPos
        s.setX(s.x() - cbArrowHeadSize);
        d.setX(d.x() - cbArrowHeadSize);
        //src
        if (item->getSourceArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone) {
            painter->drawPolygon(createArrowHead(ConnectionBaseGraphicsItem::ArrowHeadDirectionRight, item->getSourcePoint()));
            path.moveTo(s);
        } else {
            path.lineTo(s);
        }
        //curve
        path.cubicTo(QPointF(s.x() - defl, s.y()),
                     QPointF(d.x() - defl, d.y()),
                     d);
        //dst
        if (item->getDestinationArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone) {
            painter->drawPolygon(createArrowHead(ConnectionBaseGraphicsItem::ArrowHeadDirectionRight, item->getDestinationPoint()));
        } else {
            path.lineTo(item->getDestinationPoint());
        }
        break;
    case PropertyLinkArrowGraphicsItem::DockPositionsLeftRight:
        //setPos
        s.setX(s.x() - cbArrowHeadSize);
        d.setX(d.x() + cbArrowHeadSize);
        //src
        if (item->getSourceArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone) {
            painter->drawPolygon(createArrowHead(ConnectionBaseGraphicsItem::ArrowHeadDirectionRight, item->getSourcePoint()));
            path.moveTo(s);
        } else {
            path.lineTo(s);
        }
        //curve
        path.cubicTo(QPointF(s.x() - defl, s.y()),
                     QPointF(d.x() + defl, d.y()),
                     d);
        //dst
        if (item->getDestinationArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone) {
            painter->drawPolygon(createArrowHead(ConnectionBaseGraphicsItem::ArrowHeadDirectionLeft, item->getDestinationPoint()));
        } else {
            path.lineTo(item->getDestinationPoint());
        }
        break;
    case PropertyLinkArrowGraphicsItem::DockPositionsRightLeft:
        //setPos
        s.setX(s.x() + cbArrowHeadSize);
        d.setX(d.x() - cbArrowHeadSize);
        //src
        if (item->getSourceArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone) {
            painter->drawPolygon(createArrowHead(ConnectionBaseGraphicsItem::ArrowHeadDirectionLeft, item->getSourcePoint()));
            path.moveTo(s);
        } else {
            path.lineTo(s);
        }
        //curve
        path.cubicTo(QPointF(s.x() + defl, s.y()),
                     QPointF(d.x() - defl, d.y()),
                     d);
        //dst
        if (item->getDestinationArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone) {
            painter->drawPolygon(createArrowHead(ConnectionBaseGraphicsItem::ArrowHeadDirectionRight, item->getDestinationPoint()));
        } else {
            path.lineTo(item->getDestinationPoint());
        }
        break;
    case PropertyLinkArrowGraphicsItem::DockPositionsRightRight:
        //setPos
        s.setX(s.x() + cbArrowHeadSize);
        d.setX(d.x() + cbArrowHeadSize);
        //src
        if (item->getSourceArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone) {
            painter->drawPolygon(createArrowHead(ConnectionBaseGraphicsItem::ArrowHeadDirectionLeft, item->getSourcePoint()));
            path.moveTo(s);
        } else {
            path.lineTo(s);
        }
        //curve
        path.cubicTo(QPointF(s.x() + defl, s.y()),
                     QPointF(d.x() + defl, d.y()),
                     d);
        //dst
        if (item->getDestinationArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone) {
            painter->drawPolygon(createArrowHead(ConnectionBaseGraphicsItem::ArrowHeadDirectionLeft, item->getDestinationPoint()));
        } else {
            path.lineTo(item->getDestinationPoint());
        }
        break;
    default:
        tgtAssert(false, "shouldn't get here");
    }
    //draw path
    painter->setBrush(Qt::NoBrush);
    painter->drawPath(path);
}

//-------------------------------------
//      PortSizeLinkArrowGraphicsItem
//-------------------------------------
QRectF NWEStyle_Classic::PortSizeLinkArrowGI_boundingRect(const PortSizeLinkArrowGraphicsItem* item) const {
    return PortSizeLinkArrowGI_shape(item).boundingRect();
}

QPainterPath NWEStyle_Classic::PortSizeLinkArrowGI_shape(const PortSizeLinkArrowGraphicsItem* item) const {
    QPointF sourcePoint = item->getSourcePoint();
    QPointF destinationPoint = item->getDestinationPoint();

    qreal defl = fabs(sourcePoint.y() - destinationPoint.y()) / 2.1 + fabs(sourcePoint.x() - destinationPoint.x()) / 2.1;
    qreal p = portArrowHoverSize;

    QPainterPath path(sourcePoint);
    path.moveTo(sourcePoint.x()-p,sourcePoint.y());
    path.quadTo(QPointF((sourcePoint.x()+destinationPoint.x()) / 2.f-p, (sourcePoint.y()+destinationPoint.y()) / 2.f - defl),
                QPointF(destinationPoint.x()-p,destinationPoint.y()));
    path.lineTo(destinationPoint.x()+p,destinationPoint.y());
    path.quadTo(QPointF((sourcePoint.x()+destinationPoint.x()) / 2.f+p, (sourcePoint.y()+destinationPoint.y()) / 2.f - defl),
                QPointF(sourcePoint.x()+p,sourcePoint.y()));
    path.lineTo(sourcePoint.x()-p,sourcePoint.y());

    path.addPolygon(createArrowHead(QPointF((sourcePoint.x()+destinationPoint.x()) / 2.f, (sourcePoint.y()+destinationPoint.y()) / 2.f - defl)
                                         , item->getDestinationPoint()));

    return path;
}

void NWEStyle_Classic::PortSizeLinkArrowGI_initializePaintSettings(PortSizeLinkArrowGraphicsItem* item) {}

void NWEStyle_Classic::PortSizeLinkArrowGI_paint(PortSizeLinkArrowGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    QPointF sourcePoint = item->getSourcePoint();
    QPointF destinationPoint = item->getDestinationPoint();

    qreal defl = fabs(sourcePoint.y() - destinationPoint.y()) / 2.1 + fabs(sourcePoint.x() - destinationPoint.x()) / 2.1;

    QPainterPath path(sourcePoint);

    path.quadTo(QPointF((sourcePoint.x()+destinationPoint.x()) / 2.f, (sourcePoint.y()+destinationPoint.y()) / 2.f - defl),
                 destinationPoint);

    //get color and depth settings
    painter->setPen(QPen(setting.color1_, 2, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    painter->drawPath(path);
    painter->setBrush(QBrush(setting.color1_));
    painter->drawPolygon(createArrowHead(QPointF((sourcePoint.x()+destinationPoint.x()) / 2.f, (sourcePoint.y()+destinationPoint.y()) / 2.f - defl)
                                         , item->getDestinationPoint()));
}

    /*********************************************************************
     *                       ToolTip Elements
     ********************************************************************/

//-------------------------------------
//      ToolTipPortGraphicsItem
//-------------------------------------
QRectF NWEStyle_Classic::ToolTipPortGI_boundingRect(const ToolTipPortGraphicsItem* item) const {
    QRectF bRect = item->getToolTipTextItem()->boundingRect();
    if(item->getImage()) {
        float imageHeight = bRect.width() - 6.f;
        if(item->getImageAspect() > 1)
            imageHeight = imageHeight / item->getImageAspect();
        // create render tooltip
        bRect.setHeight(bRect.height() + imageHeight + 6.f);
    }
    return bRect;
}

QPainterPath NWEStyle_Classic::ToolTipPortGI_shape(const ToolTipPortGraphicsItem* item) const {
    QPainterPath path;
    path.addRect(ToolTipPortGI_boundingRect(item));
    return path;
}

void NWEStyle_Classic::ToolTipPortGI_initializePaintSettings(ToolTipPortGraphicsItem* item) {}

void NWEStyle_Classic::ToolTipPortGI_paint(ToolTipPortGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    QRectF bRect = ToolTipPortGI_boundingRect(item);
    painter->setBrush(setting.color1_);
    painter->drawRect(bRect);
    if (item->getImage()) {

        float imageHeight = bRect.width() - 6.f;
        float imageWidth = bRect.width() - 6.f;
        if(item->getImageAspect() > 1)
            imageHeight = imageHeight / item->getImageAspect();
        else
            imageWidth = imageWidth * item->getImageAspect();

        // checkers background
        QPixmap pm(20, 20);
        QPainter pmp(&pm);
        pmp.fillRect(0, 0, 10, 10, Qt::lightGray);
        pmp.fillRect(10, 10, 10, 10, Qt::lightGray);
        pmp.fillRect(0, 10, 10, 10, Qt::darkGray);
        pmp.fillRect(10, 0, 10, 10, Qt::darkGray);
        pmp.end();
        painter->setBrush(QBrush(pm));
        painter->drawRect((bRect.width() - imageWidth)/2.f,item->getToolTipTextItem()->boundingRect().height()+3.f,imageWidth, imageHeight);
        painter->drawImage(QRectF((bRect.width() - imageWidth)/2.f,item->getToolTipTextItem()->boundingRect().height()+3.f,imageWidth, imageHeight)
                            , *item->getImage());
    }
}

//-------------------------------------
//      ToolTipProcessorGraphicsItem
//-------------------------------------
QRectF NWEStyle_Classic::ToolTipProcessorGI_boundingRect(const ToolTipProcessorGraphicsItem* item) const {
    QRectF bRect = item->getToolTipTextItem()->boundingRect();
    return bRect;
}

QPainterPath NWEStyle_Classic::ToolTipProcessorGI_shape(const ToolTipProcessorGraphicsItem* item) const {
    QPainterPath path;
    path.addRect(ToolTipProcessorGI_boundingRect(item));
    return path;
}

void NWEStyle_Classic::ToolTipProcessorGI_initializePaintSettings(ToolTipProcessorGraphicsItem* item) {}

void NWEStyle_Classic::ToolTipProcessorGI_paint(ToolTipProcessorGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    QRectF bRect = ToolTipProcessorGI_boundingRect(item);
    painter->setBrush(NWEStyle_ToolTipBackgroundColor);
    painter->drawRect(bRect);
}

//-------------------------------------
//      TextBoxBaseGraphicsItem
//-------------------------------------
QRectF NWEStyle_Classic::TextBoxBaseGI_boundingRect(const TextBoxBaseGraphicsItem* item) const {
    QRectF bRect(item->mapFromScene(item->pos()),item->mapFromScene(item->pos()+QPointF(item->getCurrentSize().width(),item->getCurrentSize().height())));
    return bRect;
}

QPainterPath NWEStyle_Classic::TextBoxBaseGI_shape(const TextBoxBaseGraphicsItem* item) const {
    QPainterPath path;
    path.addRect(TextBoxBaseGI_boundingRect(item));
    return path;
}

void NWEStyle_Classic::TextBoxGI_initializePaintSettings(TextBoxGraphicsItem* item) {
    //same font settings
    QFont font(item->getContentEditor()->font());
    //set background color
    std::stringstream sstream;
    sstream << "background-color: rgb(" << NWEStyle_TextBoxBaseMainColor.red() << "," << NWEStyle_TextBoxBaseMainColor.green() <<
               "," << NWEStyle_TextBoxBaseMainColor.blue() << ")";
    item->getContentEditor()->setStyleSheet(sstream.str().c_str());
    //restore font size
    item->getContentEditor()->setFont(font);
}

void NWEStyle_Classic::TextBoxGI_paint(TextBoxGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    qreal border = item->getResizeBorder();

    item->setZValue(setting.zValue_);

    QRectF br = TextBoxBaseGI_boundingRect(item);

    //main
    painter->setBrush(NWEStyle_TextBoxBaseMainColor);
    painter->drawRect(br.x(),br.y(),br.width(),br.height());
        //highlight
    if(setting.color1_ != NWEStyle_TextBoxBaseMainColor) {
        painter->setPen(Qt::NoPen); //remove lines
        QLinearGradient gradientTop(0, 0, 0, border);
        gradientTop.setColorAt(0.0, setting.color1_);
        gradientTop.setColorAt(1.0, QColor(setting.color1_.red(),setting.color1_.green(),setting.color1_.blue(),0));
        QBrush brushTop(gradientTop);
        painter->setBrush(brushTop);
        painter->drawRect(br.x(),br.y(),br.width(),border);
        QLinearGradient gradientBottom(0, br.height()-border, 0, br.height());
        gradientBottom.setColorAt(0.0, QColor(setting.color1_.red(),setting.color1_.green(),setting.color1_.blue(),0));
        gradientBottom.setColorAt(1.0, setting.color1_);
        QBrush brushBottom(gradientBottom);
        painter->setBrush(brushBottom);
        painter->drawRect(br.x(),br.height()-border,br.width(),border);
        QLinearGradient gradientLeft(0, 0, border, 0);
        gradientLeft.setColorAt(0.0, setting.color1_);
        gradientLeft.setColorAt(1.0, QColor(setting.color1_.red(),setting.color1_.green(),setting.color1_.blue(),0));
        QBrush brushLeft(gradientLeft);
        painter->setBrush(brushLeft);
        painter->drawRect(br.x(),br.y(),border,br.height());
        QLinearGradient gradientRight(br.width()-border, 0, br.width(), 0);
        gradientRight.setColorAt(0.0, QColor(setting.color1_.red(),setting.color1_.green(),setting.color1_.blue(),0));
        gradientRight.setColorAt(1.0, setting.color1_);
        QBrush brushRight(gradientRight);
        painter->setBrush(brushRight);
        painter->drawRect(br.width()-border,br.y(),border,br.height());
    }
}

void NWEStyle_Classic::FrameBoxGI_initializePaintSettings(FrameBoxGraphicsItem* item) {}

void NWEStyle_Classic::FrameBoxGI_paint(FrameBoxGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    qreal border = item->getResizeBorder();
    qreal header = (item->getCaptionItem()->isVisible() ? item->getCaptionItem()->boundingRect().height() + 2*border : 0.f);
    QColor col = item->getBaseColor();
    item->setZValue(setting.zValue_);

    QRectF br = TextBoxBaseGI_boundingRect(item);

    //header (only frame)
    painter->setBrush(col);
    painter->drawRect(br.x(),br.y(),br.width(),header);


    //main
    painter->setBrush(col);
    painter->drawRect(br.x(),header,br.width(),br.height()-header);

    //highlight
    if(setting.color1_ != item->getBaseColor()) {
        painter->setPen(Qt::NoPen); //remove lines
        QLinearGradient gradientTop(0, 0, 0, border);
        gradientTop.setColorAt(0.0, setting.color1_);
        gradientTop.setColorAt(1.0, QColor(setting.color1_.red(),setting.color1_.green(),setting.color1_.blue(),0));
        QBrush brushTop(gradientTop);
        painter->setBrush(brushTop);
        painter->drawRect(br.x(),br.y(),br.width(),border);
        QLinearGradient gradientBottom(0, br.height()-border, 0, br.height());
        gradientBottom.setColorAt(0.0, QColor(setting.color1_.red(),setting.color1_.green(),setting.color1_.blue(),0));
        gradientBottom.setColorAt(1.0, setting.color1_);
        QBrush brushBottom(gradientBottom);
        painter->setBrush(brushBottom);
        painter->drawRect(br.x(),br.height()-border,br.width(),border);
        QLinearGradient gradientLeft(0, 0, border, 0);
        gradientLeft.setColorAt(0.0, setting.color1_);
        gradientLeft.setColorAt(1.0, QColor(setting.color1_.red(),setting.color1_.green(),setting.color1_.blue(),0));
        QBrush brushLeft(gradientLeft);
        painter->setBrush(brushLeft);
        painter->drawRect(br.x(),br.y(),border,br.height());
        QLinearGradient gradientRight(br.width()-border, 0, br.width(), 0);
        gradientRight.setColorAt(0.0, QColor(setting.color1_.red(),setting.color1_.green(),setting.color1_.blue(),0));
        gradientRight.setColorAt(1.0, setting.color1_);
        QBrush brushRight(gradientRight);
        painter->setBrush(brushRight);
        painter->drawRect(br.width()-border,br.y(),border,br.height());
    }

    //background (only frame)
    painter->setCompositionMode(QPainter::CompositionMode_Clear);
    painter->drawRect(br.x()+border,header+border,br.width()-2*border,br.height()-header-2*border);
    painter->setCompositionMode(QPainter::CompositionMode_SourceOver);
    painter->setBrush(col);
    painter->setOpacity(0.5f);
    painter->drawRect(br.x()+border,header+border,br.width()-2*border,br.height()-header-2*border);

}


} // namespace voreen
