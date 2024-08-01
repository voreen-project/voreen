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

#ifndef VRN_PORTGRAPHICSITEM_H
#define VRN_PORTGRAPHICSITEM_H

#include "propertyownergraphicsitem.h"
#include "voreen/core/ports/portwidget.h"

namespace voreen {

class Port;
class PortOwnerGraphicsItem;
class PortArrowGraphicsItem;
class PortSizeLinkArrowGraphicsItem;

class PortGraphicsItem : public PropertyOwnerGraphicsItem, protected PortWidget {
Q_OBJECT
public:
    //constructor + destructor
    PortGraphicsItem(Port* port, PortOwnerGraphicsItem* parent);
    virtual ~PortGraphicsItem();

    //---------------------------------------------------------------------------------
    //      nwebasegraphicsitem functions
    //---------------------------------------------------------------------------------
    int type() const {return UserTypesPortGraphicsItem; }
    void updateNWELayerAndCursor();

    //style functions
    QRectF boundingRect() const;
    QPainterPath shape() const;
protected:
    void initializePaintSettings();
    void mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);

    // Inherited from PortWidget
    virtual void updateFromPort();

    //handling child items (arrows)
public:
    void layoutChildItems();
protected:
    void createChildItems();
    void deleteChildItems();

    //menu
    void setContextMenuActions();
    //---------------------------------------------------------------------------------
    //      portgraphicsitem functions
    //---------------------------------------------------------------------------------
public:
    //getter and setter
    Port* getPort() const;
    PortOwnerGraphicsItem* getPortOwner() const;
    QList<PortGraphicsItem*> getConnectedOutportGraphicsItems() const;
    const QList<PortArrowGraphicsItem*>& getArrowList() const;
    QColor getColor() const;
    bool isOutport() const;
    void setVisible(bool visible);

    //port functions
    void update(const QRectF& rect = QRectF()); //needed?
        //port size arrows
    PortSizeLinkArrowGraphicsItem* setCurrentPortSizeLinkArrow(PortSizeLinkArrowGraphicsItem* arrow);
    void addPortSizeLinkArrow(PropertyGraphicsItem* src, PropertyGraphicsItem* dst);
    void removePortSizeLinkArrow(PropertyGraphicsItem* src, PropertyGraphicsItem* dst);
        //port arrows
    PortArrowGraphicsItem* setCurrentPortArrow(PortArrowGraphicsItem* arrow);
    void addGraphicalConnection(PortGraphicsItem* port);
    void removeGraphicalConnection(PortGraphicsItem* port);
protected:
    void addConnectedOutportGraphicsItem(PortGraphicsItem* port);
    void removeConnectedOutportGraphicsItem(PortGraphicsItem* port);

public:
    //events
    void mousePressEvent(QGraphicsSceneMouseEvent* event);
    void mouseMoveEvent(QGraphicsSceneMouseEvent* event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent* event);

signals:
    //signals to toggle tooltips
    void startedArrow();
    void endedArrow();

protected slots:
    //used in context menu
    void savePortData();

    //member
private:
    Port* port_;
    PortOwnerGraphicsItem* portOwner_;

    QList<PortArrowGraphicsItem*> arrowList_;
    QList<PortGraphicsItem*> connectedOutportGraphicsItems_;
    QList<PortSizeLinkArrowGraphicsItem*> linkArrowList_;

    PortArrowGraphicsItem* currentPortArrow_;
    PortSizeLinkArrowGraphicsItem* currentPortSizeLinkArrow_;
};

} // namespace

#endif // VRN_PORTGRAPHICSITEM_H




