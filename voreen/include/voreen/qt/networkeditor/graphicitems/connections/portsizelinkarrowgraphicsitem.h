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

#ifndef VRN_PORTSIZELINKARROWGRAPHICSITEM_H
#define VRN_PORTSIZELINKARROWGRAPHICSITEM_H

#include "connectionbasegraphicsitem.h"

#include "../../networkeditor.h"
#include "voreen/qt/networkeditor/editor_settings.h"

namespace voreen {

class PortGraphicsItem;
class PropertyGraphicsItem;

class PortSizeLinkArrowGraphicsItem : public ConnectionBaseGraphicsItem {
public:
    //constuctor and destructor
    PortSizeLinkArrowGraphicsItem(PropertyGraphicsItem* sourceItem, PortGraphicsItem* sourcePort);
    virtual ~PortSizeLinkArrowGraphicsItem();

    //---------------------------------------------------------------------------------
    //      nwebasegraphicsitem functions
    //---------------------------------------------------------------------------------
    int type() const {return UserTypesPortSizeLinkArrowGraphicsItem;}
    void updateNWELayerAndCursor();

    //style
    QRectF boundingRect() const;
    QPainterPath shape() const;
protected:
    void initializePaintSettings();
    void mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);

    //---------------------------------------------------------------------------------
    //      portsizelinkarrowgraphicsitem functions
    //---------------------------------------------------------------------------------
public:
    //getter and setter
    PropertyGraphicsItem* getSourceItem() const;
    PropertyGraphicsItem* getDestinationItem() const;
    PortGraphicsItem* getSourcePort() const;
    PortGraphicsItem* getDestinationPort() const;

    void setDestinationPort(PortGraphicsItem* dst);

    QPointF getSourcePoint() const;
    QPointF getDestinationPoint() const;

protected:
    //events
    void mousePressEvent(QGraphicsSceneMouseEvent* event);
    void mouseMoveEvent(QGraphicsSceneMouseEvent* event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent* event);

    //member
    PortGraphicsItem* sourcePort_;
    PortGraphicsItem* destinationPort_;

    QPointF clickedPoint_;
    bool arrowDraged_;
    bool dragedArrowMoved_;
    bool currentArrowInitialized_;

};

} // namespace

#endif // VRN_PORTSIZELINKARROWGRAPHICSITEM_H
