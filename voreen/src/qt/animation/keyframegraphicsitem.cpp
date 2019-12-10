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

#include "voreen/qt/animation/animationexportwidget.h"
#include "voreen/qt/animation/keyframegraphicsitem.h"

#include <QPainter>
#include <QPen>
#include <QGraphicsSceneMouseEvent>
#include <iostream>

namespace voreen {

KeyframeGraphicsItem::KeyframeGraphicsItem(QRectF boundingRect, float duration)
    : boundingRect_(boundingRect)
    , duration_(duration)
    , moving_(false)
    , oldX_(0.f)
{
    setFlag(ItemIsMovable);
    setFlag(ItemIsSelectable);
    setCursor(Qt::ArrowCursor);
    setAcceptHoverEvents(true);
}

void KeyframeGraphicsItem::paint(QPainter *painter, const QStyleOptionGraphicsItem*  /*option*/, QWidget* /*widget*/) {
    QRectF itemRect(boundingRect().x()+1, boundingRect().y()+1, boundingRect().width()-1, boundingRect().height()-1);
    QLinearGradient gradient(0,0,0, 10);
    gradient.setSpread(QGradient::PadSpread);
    if(!isSelected()) {
        gradient.setColorAt(0.0, QColor(255, 255, 255, 255));
        gradient.setColorAt(1.0, QColor(0, 0, 150, 255));
    }
    else {
        gradient.setColorAt(0.0, QColor(255, 255, 255, 255));
        gradient.setColorAt(0.7, QColor(150, 0, 0, 255));
    }
    QBrush brush(gradient);
    painter->setBrush(gradient);
    painter->setOpacity(1.0);
    painter->drawRect(itemRect);
    setVisible(true);
}

QRectF KeyframeGraphicsItem::boundingRect() const {
    return boundingRect_;
    //return QRectF(-5, -10, 10, 18 );
}

void KeyframeGraphicsItem::mousePressEvent(QGraphicsSceneMouseEvent * event ) {

    if(event->modifiers() & Qt::CTRL) {
        event->ignore();
        return;
    }

    emit itemClicked(this);
    QGraphicsItem::mousePressEvent(event);
}

void KeyframeGraphicsItem::mouseMoveEvent (QGraphicsSceneMouseEvent* event ) {
    if (!moving_) {
        moving_ = true;
        oldX_ = pos().x();
    }

    if(event->scenePos().x() < 0) {
        setPos(0 , pos().y() /*10*/);
    }
    else if (event->scenePos().x() >= duration_) {
        setPos(duration_, /*10*/ pos().y());
    }
    else {
        setPos(event->scenePos().x(), pos().y() /*10*/);
    }
    emit itemMoving(this);
}

void KeyframeGraphicsItem::mouseReleaseEvent(QGraphicsSceneMouseEvent* e) {

    setPos(scenePos().x(), /*10*/ pos().y());
    if(e->modifiers() & Qt::ShiftModifier)
        emit itemReleased(this, true);
    else
        emit itemReleased(this, false);

    moving_ = false;

}

bool KeyframeGraphicsItem::moving() const {
    return moving_;
}

void KeyframeGraphicsItem::restoreOldPosition() {
    setPos(oldX_, pos().y());
}

void KeyframeGraphicsItem::setDuration(float duration) {
    duration_ = duration;
}

float KeyframeGraphicsItem::getDuration() const {
    return duration_;
}

} // namespace voreen
