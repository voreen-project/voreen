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

#ifndef VRN_PROPERTYLINKARROWGRAPHICSITEM_H
#define VRN_PROPERTYLINKARROWGRAPHICSITEM_H

#include "connectionbasegraphicsitem.h"
#include "voreen/qt/networkeditor/editor_settings.h"

namespace voreen {

class PropertyGraphicsItem;

/**
 * Class to represent the link between two properties.
 * It is a connection between two propertygraphicsitems.
 * Each arrow is stored in a list in the PropertyGraphicsItem. The destination item
 * also stores an corresponding arrow. Each arrow has a pointer to its corresponding one.
 * with this pointer hover events are passed, which is needed in the PropertyLinkDialog.
 */
class PropertyLinkArrowGraphicsItem : public ConnectionBaseGraphicsItem {
Q_OBJECT
public:
    /**
     * Enum to determine from which side of the source propertygraphicsitem is connected
     * with which side of the destination propertygraphicsitem.
     * It is used in the paint functions to draw the connections shape.
     */
    enum DockPositions {
        DockPositionsLeftLeft,  ///< connection from left side to left side
        DockPositionsLeftRight, ///< connection from left side to right side
        DockPositionsRightLeft, ///< connection from right side to left side
        DockPositionsRightRight ///< connection from right side to right side
    };

    //constructor and destructor
    PropertyLinkArrowGraphicsItem(PropertyGraphicsItem* sourceItem = 0, PropertyGraphicsItem* destinationItem = 0);
    virtual ~PropertyLinkArrowGraphicsItem();

    //---------------------------------------------------------------------------------
    //      nwebasegraphicsitem functions
    //---------------------------------------------------------------------------------
    int type() const {return UserTypesPropertyLinkArrowGraphicsItem;}

    //style
    QRectF boundingRect() const;
    QPainterPath shape() const;
protected:
    void initializePaintSettings();
    void mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);

    //---------------------------------------------------------------------------------
    //      propertylinkarrowgraphicsitem functions
    //---------------------------------------------------------------------------------
public:
    //------------------------
    //  getter and setter
    //------------------------
    PropertyLinkArrowGraphicsItem* getCorrespondingPropertyLinkArrowItem() const;
    void setCorrespondingPropertyLinkArrowItem(PropertyLinkArrowGraphicsItem* arrow);
    PropertyGraphicsItem* getSourceItem() const;
    PropertyGraphicsItem* getDestinationItem() const;
    QPointF getSourcePoint() const;
    QPointF getDestinationPoint() const;
    DockPositions getDockPositions() const;

protected:
    //------------------------
    //  members
    //------------------------
    PropertyLinkArrowGraphicsItem* correspondingPropertyLinkArrowItem_;     ///< pointer to the corresponding arrow in the destination property
    QPointF clickedPoint_;      ///< ???
    bool initCurrentArrow_;     ///< ???
    bool movedAwayInEvent_;     ///< ???

    //---------------------------------------------------------------------------------
    //      event functions
    //---------------------------------------------------------------------------------
public:
    void mousePressEvent(QGraphicsSceneMouseEvent* event);
    void mouseMoveEvent(QGraphicsSceneMouseEvent* event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent* event);
protected:
    void hoverEnterEvent(QGraphicsSceneHoverEvent* event);
    void hoverLeaveEvent(QGraphicsSceneHoverEvent* event);
signals:
    /** Emmited in hoverEnterEvent. Not used yet. */
    void hoverEnter(PropertyLinkArrowGraphicsItem* item);
    /** Emmited in hoverLeaveEvent. Not used yet. */
    void hoverLeave(PropertyLinkArrowGraphicsItem* item);
};

} // namespace

#endif // VRN_PROPERTYLINKARROWGRAPHICSITEM_H
