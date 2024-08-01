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

#ifndef VRN_PORTOWNERLINKARROWGRAPHICSITEM_H
#define VRN_PORTOWNERLINKARROWGRAPHICSITEM_H

#include "connectionbasegraphicsitem.h"
#include "voreen/qt/networkeditor/editor_settings.h"

namespace voreen {

class PropertyListButtonGraphicsItem;

class PortOwnerLinkArrowGraphicsItem : public ConnectionBaseGraphicsItem {
public:
    //constructor and destructor
    PortOwnerLinkArrowGraphicsItem(PropertyListButtonGraphicsItem* sourceItem, PropertyListButtonGraphicsItem* destinationItem = 0);
    virtual ~PortOwnerLinkArrowGraphicsItem();

    //---------------------------------------------------------------------------------
    //      nwebasegraphicsitem functions
    //---------------------------------------------------------------------------------
    int type() const {return UserTypesPortOwnerLinkArrowGraphicsItem;}
    void updateNWELayerAndCursor();

    //style functions
    QRectF boundingRect() const;
    QPainterPath shape() const;
protected:
    void initializePaintSettings();
    void mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);


    //---------------------------------------------------------------------------------
    //      portownerlinkarrowgraphicsitem functions
    //---------------------------------------------------------------------------------
public:
    //getter and setter
    PropertyListButtonGraphicsItem* getSourceItem() const;
    PropertyListButtonGraphicsItem* getDestinationItem() const;

    QPointF getSourcePoint() const;
    QPointF getDestinationPoint() const;

    void setVisible(bool visible);
};

} // namespace

#endif // VRN_PORTOWNERLINKARROWGRAPHICSITEM_H
