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

#ifndef VRN_PORTARROWGRAPHICSITEM_H
#define VRN_PORTARROWGRAPHICSITEM_H

#include "connectionbasegraphicsitem.h"
#include "voreen/qt/networkeditor/editor_settings.h"

#include "../../networkeditor.h"

namespace voreen {

class PortGraphicsItem;

class PortArrowGraphicsItem : public ConnectionBaseGraphicsItem {

public:
    //constructor and destructor
    PortArrowGraphicsItem(PortGraphicsItem* sourceItem);
    virtual ~PortArrowGraphicsItem();

    //---------------------------------------------------------------------------------
    //      nwebasegraphicsitem functions
    //---------------------------------------------------------------------------------
    void updateNWELayerAndCursor();
    int type() const {return UserTypesPortArrowGraphicsItem;}

    QRectF boundingRect() const;
    QPainterPath shape() const;
    void mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);
protected:
    void initializePaintSettings();

    //---------------------------------------------------------------------------------
    //      portarrowgraphicsitem functions
    //---------------------------------------------------------------------------------
public:
    //getter and setter
    PortGraphicsItem* getSourceItem() const;
    PortGraphicsItem* getDestinationItem() const;
    QPointF getSourcePoint() const;
    QPointF getDestinationPoint() const;

protected:
    //events
    void mousePressEvent(QGraphicsSceneMouseEvent* event);
    void mouseMoveEvent(QGraphicsSceneMouseEvent* event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent* event);

    QPointF clickedPoint_;      //position of left click
    bool arrowDraged_;
    bool dragedArrowMoved_;
    bool currentArrowInitialized_;     //indicates if the current arrow has been set

};

    //void setBundleInfo(bool enable, ConnectionBundle* bundle = 0);
    //bool isBundled() const;
    //ConnectionBundle* bundle_;
    //bool bundled_;

} // namespace

#endif // VRN_PORTARROWGRAPHICSITEM_H
