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

#ifndef VRN_KEYFRAMEGRAPHICSITEM_H
#define VRN_KEYFRAMEGRAPHICSITEM_H

#include <QGraphicsItem>
#include <QGraphicsView>
#include <QLabel>
#include <QObject>

namespace voreen {

/**
* QGraphicsItem representing a KeyFrame
*/
#if (QT_VERSION >= 0x040600)
class KeyframeGraphicsItem : public QGraphicsObject {
#else
class KeyframeGraphicsItem : public QObject, public QGraphicsItem {
#endif
Q_OBJECT

public:
    KeyframeGraphicsItem(QRectF boundingRect = QRectF(-5, -10, 10, 18 ), float duration = 0.f);
    QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem* option, QWidget* widget);

    void setDuration(float duration);

    float getDuration() const;

    void restoreOldPosition();

    bool moving() const;

protected:

    virtual void mousePressEvent(QGraphicsSceneMouseEvent* event);
    virtual void mouseMoveEvent(QGraphicsSceneMouseEvent* event);
    virtual void mouseReleaseEvent(QGraphicsSceneMouseEvent* event);
    void hoverEnterEvent(QGraphicsSceneHoverEvent*){emit itemHoverEnter(this);}
    void hoverLeaveEvent(QGraphicsSceneHoverEvent*){emit itemHoverLeave(this);}
    void hoverMoveEvent(QGraphicsSceneHoverEvent*){emit itemHoverMove(this);}

private:
    /// actual value of the Keyframe
    int frame_;

    /// size of the item
    QRectF boundingRect_;

    float duration_;
    bool moving_;
    float oldX_;

signals:
    void itemClicked(KeyframeGraphicsItem*);
    void itemReleased(KeyframeGraphicsItem*, bool shift);
    void itemMoving(KeyframeGraphicsItem* kfgi);
    void itemHoverEnter(KeyframeGraphicsItem* kfgi);
    void itemHoverMove(KeyframeGraphicsItem* kfgi);
    void itemHoverLeave(KeyframeGraphicsItem* kfgi);
};

} // namespace voreen

#endif
