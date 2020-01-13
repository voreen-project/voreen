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

#ifndef VRN_CONNECTIONBASEGRAPHICSITEM_H
#define VRN_CONNECTIONBASEGRAPHICSITEM_H

#include "../nwebasegraphicsitem.h"

namespace voreen {

/**
 * Basic class of all connection graphicsitems.
 */
class ConnectionBaseGraphicsItem : public NWEBaseGraphicsItem {
public:
    /**
     * Enum used to determine the color of a connection based on the actual state.
     */
    enum ColorConnectableMode {
        CCM_DEFAULT,    ///< connection is in default state (default color)
        CCM_NO,         ///< error state (normally red)
        CCM_MAYBE,      ///< some modifications are needed (normally yellow)
        CCM_YES         ///< supported state (normally green)
    };

    /**
     * Enum determine the arrow head direction if it exists.
     */
    enum ArrowHeadDirection {
        ArrowHeadDirectionNone, ///< no arrow head
        ArrowHeadDirectionAny,  ///< not specified direction. E.g. used in PropertyLinkArrow
        ArrowHeadDirectionUp,   ///< up pointing arrow
        ArrowHeadDirectionDown, ///< down pointing arrow
        ArrowHeadDirectionLeft, ///< left pointing arrow
        ArrowHeadDirectionRight ///< right pointing arrow
    };

    /**
     * A connection is always between two NWEBaseGraphicsItems.
     * The second item can be zero in case of a not yet established connection.
     */
    ConnectionBaseGraphicsItem(NWEBaseGraphicsItem* sourceItem, NWEBaseGraphicsItem* destinationItem = 0);
    virtual ~ConnectionBaseGraphicsItem();

    //---------------------------------------------------------------------------------
    //      nwebasegraphicsitem functions
    //---------------------------------------------------------------------------------
    /** Connections are not layer sensitive */
    virtual void updateNWELayerAndCursor() {};
    /** Connections have no children */
    virtual void layoutChildItems() {};
protected:
    /** Connections have no children */
    virtual void createChildItems() {};
    /** Connections have no children */
    virtual void deleteChildItems() {};

    //---------------------------------------------------------------------------------
    //      connectionbaseebasegraphicsitem functions
    //---------------------------------------------------------------------------------
public:
    //------------------------
    //  getter and setter
    //------------------------
        //the connectable mode
    ColorConnectableMode getColorConnectableMode() const { return colorConnectableMode_;}
    void setColorConnectableMode(ColorConnectableMode ccm) {colorConnectableMode_ = ccm;}
        //arrow head directions
    ArrowHeadDirection getSourceArrowHeadDirection() const;
    void setSourceArrowHeadDirection(ArrowHeadDirection ahd);
    ArrowHeadDirection getDestinationArrowHeadDirection() const;
    void setDestinationArrowHeadDirection(ArrowHeadDirection ahd);
        //the source/destination items
    NWEBaseGraphicsItem* getSourceItem() const;
    void setSourceItem(NWEBaseGraphicsItem* item);
    NWEBaseGraphicsItem* getDestinationItem() const;
    void setDestinationItem(NWEBaseGraphicsItem* item);
        //source/destination points
    virtual QPointF getSourcePoint() const;
    void setSourcePoint(const QPointF& sourcePoint);
    virtual QPointF getDestinationPoint() const;
    void setDestinationPoint(const QPointF& destinationPoint);

    //------------------------
    //  hovered functions
    //------------------------
public:
    bool getIsHovered() const;
    void setIsHovered(bool b);
protected:
    /**
     * isHovered_ is used to simulate this connection would be hoverd.
     * It is used in the paint function to change the color of the connection
     * as if the mouse would be over it.
     */
    bool isHovered_; ///< member

    //------------------------
    //  members
    //------------------------
protected:
    ColorConnectableMode colorConnectableMode_;         ///< actual connectable mode
    ArrowHeadDirection sourceArrowHeadDirection_;       ///< source arrow head
    ArrowHeadDirection destinationArrowHeadDirection_;  ///< destination arrow head
    NWEBaseGraphicsItem* sourceItem_;                   ///< pointer to the source item
    NWEBaseGraphicsItem* destinationItem_;              ///< pointer to the destination item
    QPointF sourcePoint_;                               ///< the source point
    QPointF destinationPoint_;                          ///< the destination point
};

} // namespace voreen

#endif // VRN_CONNECTIONBASEGRAPHICSITEM_H
