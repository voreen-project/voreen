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

#include "voreen/qt/networkeditor/styles/nwestyle_material.h"
#include "voreen/qt/networkeditor/editor_settings.h"

#include "voreen/qt/voreenapplicationqt.h"

//gi
#include "voreen/qt/networkeditor/graphicitems/core/portgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/processorgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/utils/progressbargraphicsitem.h"


//qt
#include <QPainter>
#include <QPainterPath>
#include <QStyle>
#include <QStyleOption>

namespace voreen{

/*********************************************************************
 *                       General Color Defines
 ********************************************************************/
//editor
QBrush NWEStyle_Material::getBackgroundBrush() const {
    return QBrush(QColor(255, 255, 255, 255));
}

//general graphicsitem
QColor NWEStyle_Material::getSelectionColor() const {
    return QColor(160,160,160,255);
}
QColor NWEStyle_Material::getHoverColor() const {
    return QColor(200,200,200,255);
}

//port
//processor
QColor NWEStyle_Material::getProcessorColor1() const {
    return QColor(245,245,245,255);
}
//portarrow
QColor NWEStyle_Material::getPortArrowColor() const {
    return QColor(80, 80, 80, 255);
}

QColor NWEStyle_Material::getPropertyLinkArrowColor() const {
    return QColor(40, 40, 40, 255);
}
QColor NWEStyle_Material::getPortSizeLinkArrowColor() const {
    return QColor(20, 20, 20, 255);
}

// shadows
bool NWEStyle_Material::getShadowsEnabled() const {
    return true;
}

/*********************************************************************
 *                       General Functions
 ********************************************************************/
NWEStyle_Material::NWEStyle_Material(NetworkEditor* networkeditor)
    : NWEStyle_Classic(networkeditor)
{}

NWEStyle_Material::~NWEStyle_Material()
{}

/*********************************************************************
 *                       Core Elements
 ********************************************************************/

//-------------------------------------
//      PortGraphicsItem
//-------------------------------------
void NWEStyle_Material::PortGI_paint(PortGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    //get color and depth settings
    painter->setPen(Qt::NoPen);
    painter->setBrush(setting.error_ != ES_NO_ERROR ? QBrush(Qt::red) : QBrush(setting.color1_));
    if (option->state & QStyle::State_Sunken)
        painter->setBrush(painter->brush().color().darker(200));
    if (option->state & QStyle::State_MouseOver) {
        painter->setPen(QPen(getPortArrowColor(), 2.0f));
        painter->setBrush(QBrush(setting.color1_.lighter()));
    }
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
void NWEStyle_Material::ProcessorGI_initializePaintSettings(ProcessorGraphicsItem* item) {
    //set text
    item->getOrCreateNameLabel()->setPlainText(item->getGuiName());
    item->getOrCreateNameLabel()->setFont(QFont("Helvetica", nwe_->getProcessorFontSizeScale()/10.f));
    item->getOrCreateNameLabel()->setDefaultTextColor(Qt::black);
}

/*
// Helper function drawing more efficient & simple shadows than using QDropShadowEffect.
void drawShadow(QPainter* &painter, float margin, float radius,
                const QColor& start, const QColor& end, float startPosition, float endPosition0, float endPosition1, QRectF rect) {
    painter->setPen(Qt::NoPen);

    QLinearGradient gradient;
    gradient.setColorAt(startPosition, start);
    gradient.setColorAt(endPosition0, end);
    // Right
    QPointF right0(rect.width() - margin, rect.y() + rect.height() / 2);
    QPointF right1(rect.width(), rect.y() + rect.height() / 2);
    gradient.setStart(right0);
    gradient.setFinalStop(right1);
    painter->setBrush(QBrush(gradient));
    painter->drawRoundRect(QRectF(QPointF(rect.width() - margin * radius, margin), QPointF(rect.width(), rect.height() - margin)), 0, 0);
    // Left
    QPointF left0(margin, rect.y() + rect.height() / 2);
    QPointF left1(0, rect.y() + rect.height() / 2);
    gradient.setStart(left0);
    gradient.setFinalStop(left1);
    painter->setBrush(QBrush(gradient));
    painter->drawRoundRect(QRectF(QPointF(margin * radius, margin), QPointF(0, rect.height() - margin)), 0, 0);
    // Top
    QPointF top0(rect.width() / 2, rect.y() + margin);
    QPointF top1(rect.width() / 2, rect.y());
    gradient.setStart(top0);
    gradient.setFinalStop(top1);
    painter->setBrush(QBrush(gradient));
    painter->drawRoundRect(QRectF(QPointF(rect.width() - margin, 0), QPointF(margin, margin)), 0, 0);
    // Bottom
    QPointF bottom0(rect.x() + rect.width() / 2, rect.y() + rect.height() - margin);
    QPointF bottom1(rect.x() + rect.width() / 2, rect.y() + rect.height());
    gradient.setStart(bottom0);
    gradient.setFinalStop(bottom1);
    painter->setBrush(QBrush(gradient));
    painter->drawRoundRect(QRectF(QPointF(margin, rect.height() - margin), QPointF(rect.width() - margin, rect.height())), 0, 0);
    // BottomRight
    QPointF bottomright0(rect.x() + rect.width() - margin, rect.y() + rect.height() - margin);
    QPointF bottomright1(rect.x() + rect.width(), rect.y() + rect.height());
    gradient.setStart(bottomright0);
    gradient.setFinalStop(bottomright1);
    gradient.setColorAt(endPosition1, end);
    painter->setBrush(QBrush(gradient));
    painter->drawRoundRect(QRectF(bottomright0, bottomright1), 0.0, 0.0);
    // BottomLeft
    QPointF bottomleft0(rect.x() + margin, rect.y() + rect.height() - margin);
    QPointF bottomleft1(rect.x(), rect.y() + rect.height());
    gradient.setStart(bottomleft0);
    gradient.setFinalStop(bottomleft1);
    gradient.setColorAt(endPosition1, end);
    painter->setBrush(QBrush(gradient));
    painter->drawRoundRect(QRectF(bottomleft0, bottomleft1), 0.0, 0.0);
    // TopLeft
    QPointF topleft0(rect.x() + margin, rect.y() + margin);
    QPointF topleft1(rect.x(), rect.y());
    gradient.setStart(topleft0);
    gradient.setFinalStop(topleft1);
    gradient.setColorAt(endPosition1, end);
    painter->setBrush(QBrush(gradient));
    painter->drawRoundRect(QRectF(topleft0, topleft1), 0.0, 0.0);
    // TopRight
    QPointF topright0(rect.x() + rect.width() - margin, rect.y() + margin);
    QPointF topright1(rect.x() + rect.width(), rect.y());
    gradient.setStart(topright0);
    gradient.setFinalStop(topright1);
    gradient.setColorAt(endPosition1, end);
    painter->setBrush(QBrush(gradient));
    painter->drawRoundRect(QRectF(topright0, topright1), 0.0, 0.0);
}
*/

void NWEStyle_Material::ProcessorGI_paint(ProcessorGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    QRectF boundRect = ProcessorGI_boundingRect(item);
    //should be smaller to prevent trails in movement
    QRectF drawRect(boundRect.left()+1,boundRect.top()+1,boundRect.width()-2,boundRect.height()-2);

    //set color settings
    painter->setOpacity(setting.opacity_);
    item->setZValue(setting.zValue_);
    painter->setBrush(QBrush(setting.color1_));
    painter->setPen(Qt::NoPen);

    //draw processor
    painter->drawRoundedRect(drawRect, 3, 3);

    // TODO: shadows are currently drawn using QDropShadowEffect, which is slow. Try to mimic using gradients.
    //draw shadow
    //drawShadow(painter, 10, 1, QColor(120, 120, 120, 32), QColor(255, 255, 255, 0), 0.0f, 1.0f, 0.8f, boundRect);

    //draw error signe
    switch(setting.error_){
        case ES_ERROR_T1:{
            painter->setOpacity(1.0);
            qreal size = drawRect.height()/4.0;
            NWEStyle_Base::NWEStyle_Error1Renderer.render(painter,QRectF(drawRect.width()-size-buttonsOffsetX/2.0,drawRect.height()-size-buttonsOffsetY/2.0,size,size));
            break;}
        case ES_ERROR_T2:{
            painter->setOpacity(1.0);
            qreal size = drawRect.height()/4.0;
            NWEStyle_Base::NWEStyle_Error2Renderer.render(painter,QRectF(drawRect.width()-size-buttonsOffsetX/2.0,drawRect.height()-size-buttonsOffsetY/2.0,size,size));
            break;}
        default:
            break;
    }
}

/*********************************************************************
 *                       Util Elements
 ********************************************************************/

//-------------------------------------
//      ProgressBarGraphicsItem
//-------------------------------------
QPainterPath NWEStyle_Material::ProgressBarGI_shape(const ProgressBarGraphicsItem* item) const {
    QRectF boundingRect = ProgressBarGI_boundingRect(item);
    QPainterPath path;
    path.addRect(boundingRect);
    return path;
}

void NWEStyle_Material::ProgressBarGI_paint(ProgressBarGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {

    QRectF boundingRect_ = ProgressBarGI_boundingRect(item);
    painter->setClipPath(ProgressBarGI_shape(item));

    QBrush backgroundBrush(item->getBackgroundColor1());
    painter->setBrush(backgroundBrush);
    painter->setPen(Qt::NoPen);
    painter->drawRect(boundingRect_.x(), boundingRect_.y(), boundingRect_.width(), boundingRect_.height());

    // now the progress overlay
    QRectF progressRect = boundingRect_;
    progressRect.setWidth(boundingRect_.width() * item->getProgress());

    QBrush foregroundBrush(item->getLowerForegroundColor());
    painter->setBrush(foregroundBrush);
    painter->drawRect(progressRect.x(), progressRect.y(), progressRect.width(), progressRect.height());

}

//-------------------------------------
//      TextBoxBaseGraphicsItem
//-------------------------------------
void NWEStyle_Material::TextBoxGI_paint(TextBoxGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    item->setZValue(setting.zValue_);

    QRectF br = TextBoxBaseGI_boundingRect(item);

    //main
    painter->setPen(Qt::NoPen);
    painter->setBrush(getTextBoxBaseMainColor());
    painter->drawRect(br.x(),br.y(),br.width(),br.height());
}

void NWEStyle_Material::FrameBoxGI_paint(FrameBoxGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
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

    //background (only frame)
    painter->setCompositionMode(QPainter::CompositionMode_Clear);
    painter->drawRect(br.x()+border,header+border,br.width()-2*border,br.height()-header-2*border);
    painter->setCompositionMode(QPainter::CompositionMode_SourceOver);
    painter->setBrush(col);
    painter->setOpacity(0.5f);
    painter->drawRect(br.x()+border,header+border,br.width()-2*border,br.height()-header-2*border);
}


} // namespace voreen
