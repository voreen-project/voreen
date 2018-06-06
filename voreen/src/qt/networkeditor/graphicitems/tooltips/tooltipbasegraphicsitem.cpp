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

#include "voreen/qt/networkeditor/graphicitems/tooltips/tooltipbasegraphicsitem.h"

#include <QGraphicsTextItem>
#include <QTimer>

namespace voreen {

ToolTipBaseGraphicsItem::ToolTipBaseGraphicsItem(NWEBaseGraphicsItem* item)
    : NWEBaseGraphicsItem(item->getNetworkEditor())
    , toolTipText_(0)
    , timer_(new QTimer(this))
    , toolTipTimerTriggertMousePosition_()
{
    //text
    toolTipText_ = new QGraphicsTextItem(this);
    toolTipText_->setTextWidth(130);
    setParent(item);
    //timer
    timer_->setInterval(NWEToolTipTime);
    connect(timer_,SIGNAL(timeout()),this,SLOT(timerTimeout()));
    //item itself
    setFlag(QGraphicsItem::ItemIgnoresTransformations,true);
    setZValue(ZValuesToolTipGraphicsItem);
    item->scene()->addItem(this);
    hide();
}

//---------------------------------------------------------------------------------------------------------------
//                  getter and setter
//---------------------------------------------------------------------------------------------------------------
const QGraphicsTextItem* ToolTipBaseGraphicsItem::getToolTipTextItem() const{
    return toolTipText_;
}

void ToolTipBaseGraphicsItem::setToolTipTimerTriggertMousePosition(QPointF& point) {
    toolTipTimerTriggertMousePosition_ = point;
}

//---------------------------------------------------------------------------------------------------------------
//                  tooltip functions
//---------------------------------------------------------------------------------------------------------------
void ToolTipBaseGraphicsItem::startTimer() {
    timer_->start();
}

void ToolTipBaseGraphicsItem::stopTimer() {
    timer_->stop();
}

void ToolTipBaseGraphicsItem::timerTimeout() {
    timer_->stop();
    bool visible = updateToolTip(toolTipTimerTriggertMousePosition_.x(),toolTipTimerTriggertMousePosition_.y());
    setVisible(visible);
}

bool ToolTipBaseGraphicsItem::updateToolTip(qreal x, qreal y) {
    if(scene()->views().empty()) {
        tgtAssert(false,"Scene has no View!!!");
        setPos(x+1,y+1); //should not get here
        return false;
    }
    x = scene()->views()[0]->mapFromScene(x,y).x();
    y = scene()->views()[0]->mapFromScene(x,y).y();
    //bounding rect
    QRectF bRect = boundingRect();
    //view rect
    QRect vRect = scene()->views()[0]->rect();

    //set new x/y
    if((x+ bRect.right() + 1) <= vRect.right())
        x += 1;
    else
        x -= (bRect.width()) + 1;

    if((y + bRect.bottom() + 1) <= vRect.bottom())
        y += 1;
    else
        y -= (bRect.height()) + 1;

    setPos(scene()->views()[0]->mapToScene(x,y));
    return true;
}

} // namespace voreen
