/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_PROPERTYGRAPHICSITEM_H
#define VRN_PROPERTYGRAPHICSITEM_H

#include "../nwebasegraphicsitem.h"
#include "voreen/qt/networkeditor/editor_settings.h"

namespace voreen {

class Property;
class PropertyListGraphicsItem;
class PropertyOwnerGraphicsItem;
class PropertyLinkArrowGraphicsItem;

/**
 * This class represents a property of a PropertyOwner in the NetworkEditor.
 * It normally grouped with the other properties of the PropertyOwner in a PropertyListGraphicsItem.
 * This list is his parent and scene parent, i.e it will be drawn id the list is been painted.
 * The class is responsible for creating the graphical representation of its property links.
 */
class PropertyGraphicsItem : public NWEBaseGraphicsItem {
Q_OBJECT

public:
    //----------------------------
    //  constructor + destructor
    //----------------------------
    /**
     * Constructor
     * @param prop the core property represented by this.
     * @param parent the propertylist storing this.
     */
    PropertyGraphicsItem(const Property* prop, PropertyListGraphicsItem* parent);
    /** Destructor */
    virtual ~PropertyGraphicsItem();

//---------------------------------------------------------------------------------
//      nwebasegraphicsitem functions
//---------------------------------------------------------------------------------
    int type() const {return UserTypesPropertyGraphicsItem;}
    void updateNWELayerAndCursor() {};
    /** prepare change of property links. */
    void prepareGeometryChange();

    //style functions
    QRectF boundingRect() const;
    QPainterPath shape() const;
protected:
    void initializePaintSettings();
    void mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);

public:
    //----------------------------
    //  handling child items
    //----------------------------
    /** Property links lay themselves out. So nothing to do here. */
    void layoutChildItems(){};
protected:
    /**
     * Since properties are created without property links and the links are
     * added later on, we simulate this behavior. Nothing to do here.
     */
    void createChildItems(){};
    /** Deletes all property link items. */
    void deleteChildItems();

    //---------------------------------------------------------------------------------
    //      propertygraphicsitem functions
    //---------------------------------------------------------------------------------
public:
    //----------------------------
    //  getter and setter
    //----------------------------
        //property
    const Property* getProperty() const;
        //property list item
    PropertyListGraphicsItem* getPropertyListItem() const;
        //label
    QGraphicsTextItem* getPropertyLabelItem();
    const QGraphicsTextItem* getPropertyLabelItem() const;
        //property owner item
    PropertyOwnerGraphicsItem* getPropertyOwnerItem() const;
        //size management
    void resetSize(qreal width);
    qreal getBoundingWidth() const;
        //lists
    QList<PropertyLinkArrowGraphicsItem*> getSourceLinkList();
    QList<PropertyLinkArrowGraphicsItem*> getDestinationLinkList();
        //docking points
    QPointF getLeftDockingPoint() const;
    QPointF getRightDockingPoint() const;
        //arrow visibility
    void setShowLinkArrows(bool b);
        //linked to list
    PropertyListGraphicsItem* getPropertyListLinkedTo() const;
    void setPropertyListLinkedTo(PropertyListGraphicsItem* list);
        // set visibility of child items
    void setVisible(bool visible);

    //----------------------------
    //  linking functions
    //----------------------------
    /**
     * Adds a graphical link from this to propItem.
     * @note the 'core' link will not be created and is independent of the
     *       graphical representation.
     */
    PropertyLinkArrowGraphicsItem* addGraphicalLink(PropertyGraphicsItem* propItem);
    /**
     * Removes a graphical link from this to propItem.
     * @note the 'core' link will not be removed and is independent of the
     *       graphical representation.
     */
    void removeGraphicalLink(PropertyGraphicsItem* propItem);
protected:
    /**
     * After creating a graphical link, the destination item will be informed
     * to add a correspondent graphical link to its destination list
     */
    PropertyLinkArrowGraphicsItem* addDestinationLink(PropertyGraphicsItem* propItem);
    /**
     * After removing a graphical link, the destination item will be informed
     * to remove the correspondent graphical link from its destination list
     */
    void removeDestinationLink(PropertyGraphicsItem* propItem);

signals:
    /** Used in propertyLinkDialog */
    void linkProperties(PropertyGraphicsItem* src, PropertyGraphicsItem* dst);

    //----------------------------
    //  members
    //----------------------------
private:
    const Property* property_;                                  ///< the property represented by this
    PropertyListGraphicsItem* propertyListItem_;                ///< the PropertyList this is belonging to
    QGraphicsTextItem* propertyLabelItem_;                      ///< the property label item

    QList<PropertyLinkArrowGraphicsItem*> sourceLinkList_;      ///< list of property links beginning at this property
    QList<PropertyLinkArrowGraphicsItem*> destinationLinkList_; ///< list of property links ending at this property

    PropertyListGraphicsItem* propertyListLinkedTo_;            ///< only links to or from this list will be drawn
    bool showLinkArrows_;                                       ///< if false, no property links will be drawn

    qreal boundingWidth_;                                       ///< the width of this. Will be set from the PropertyList
};

} //namespace voreen

#endif // VRN_PROPERTYGRAPHICSITEM_H
