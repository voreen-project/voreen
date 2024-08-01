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

#ifndef VRN_PROPERTYLISTGRAPHICSITEM_H
#define VRN_PROPERTYLISTGRAPHICSITEM_H

#include "../nwebasegraphicsitem.h"
#include "voreen/qt/networkeditor/editor_settings.h"

namespace voreen {

class Property;
class PropertyGraphicsItem;
class PropertyOwnerGraphicsItem;

/**
 * This class is a list of PropertyGraphicsItems belonging to a PropertyOwner.
 * It is used in the PropertyLinkDialog.
 * It is a Qt child of its PropertyOwnerGraphicsItem i.e. it will be deleted automaticaly if
 * the parent is been deleted.
 * But it is not an GI child i.e. it is not drawen together with his parent. This is, because
 * the list should have another zValue than its owner.
 * It has a list of sub PropertyListGraphicsItems for recursion like a processor owning ports.
 */
class PropertyListGraphicsItem : public NWEBaseGraphicsItem {
Q_OBJECT
public:
    /**
     * This enum determines the representation of the PropertyOwner label.
     * The value equals the size of the header.
     */
    enum HeaderMode{
        NO_HEADER = 0,      ///< no header is used in showing the PropertyList in the NetworkEditor.
        SMALL_HEADER = 20,  ///< a small header is used, if this list is a sub-list (e.g. Port of PropertyOwner).
        BIG_HEADER = 40     ///< a big header is used in the PropertyLinkDialog.
    };

    /**
     * This enum filters which properties should be visible in the NetworkEditor/PropertyLinkDialog
     */
    enum PropertyVisibleMode{
        SHOW_ALL_PROPERTIES,        ///< show all properties
        SHOW_LINKED_PROPERTIES,     ///< show just linked properties
        SHOW_CAMERA_PROPERTIES      ///< show just camera properties (linked and unkinked)
    };

    /**
     * Hides same special properties. Used to keep the NetworkEditor cleen.
     */
    enum PropertyVisibleModeModifier{
        HIDE_NO_PROPERTIES,                 ///< hide nothing
        HIDE_CAMERA_PROPERTIES,             ///< hide camera properties
        HIDE_SIZE_PROPERTIES,               ///< hide size properties
        HIDE_CAMERA_AND_SIZE_PROPERTIES     ///< hide camera and size properties
    };

    /**
     * This enum filters which PropertyLinkArrows belonging to the Properties should be visible.
     */
    enum LinkArrowVisibleMode{
        SHOW_NO_ARROWS,             ///< hides all arrows.
        SHOW_ALL_ARROWS,            ///< show all arrows (used in the NetworkEditor).
        SHOW_FRIEND_LIST_ARROWS     ///< show only arrows to Properties in the friend list (used in PropertyLinkDialog).
    };

    //----------------------------
    //  constructor + destructor
    //----------------------------
    /**
     * Constructor
     * @param parent the parent GraphicsItem of this list.
     */
    PropertyListGraphicsItem(PropertyOwnerGraphicsItem* parent);
    /** Destructor */
    virtual ~PropertyListGraphicsItem();

//---------------------------------------------------------------------------------
//      nwebasegraphicsitem functions
//---------------------------------------------------------------------------------
    /** Type */
    int type() const {return UserTypesPropertyListGraphicsItem;}
    /** The representation is atm not layer sensitive. */
    void updateNWELayerAndCursor() {}; //TODO: set enums depending on layer?

    //----------------------------
    //  style functions
    //----------------------------
    /** Bounding rect of this list including all sub lists. */
    QRectF boundingRect() const;
    /** Shape of this list including all sub lists. */
    QPainterPath shape() const;
protected:
    /** Sets the header label font/color/. */
    void initializePaintSettings();
    /** Draws the header of the list */
    void mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);

    //----------------------------
    //  handling child elements
    //----------------------------
public:
    /** Lays all PropertyGraphicsItems and all sub lists out. */
    void layoutChildItems();
protected:
    /** Creates all PropertyGraphicsItems for each property of the parent. */
    void createChildItems();
    /** Deletes all PropertyGraphicsItems and sets the 'parent list' of all sub lists to null. */
    void deleteChildItems();

    //---------------------------------------------------------------------------------
    //      propertylistgraphicsitem functions
    //---------------------------------------------------------------------------------
public:
    //----------------------------
    //  getter and setter
    //----------------------------
        //parent item
    PropertyOwnerGraphicsItem* getParentPropertyOwnerItem() const;
        //super item
    void setSuperPropertyListItem(PropertyListGraphicsItem* plgi);
    PropertyListGraphicsItem* getSuperPropertyListItem() const;
        //sub lists
    QList<PropertyListGraphicsItem*> getSubPropertyListItems();
    const QList<PropertyListGraphicsItem*> getSubPropertyListItems() const;
        //header label
    const QGraphicsTextItem* getHeaderLabelItem() const;
    QGraphicsTextItem* getHeaderLabelItem();
    void updateHeaderLabelItem();
        //visible
    bool getIsVisibleInEditor() const;
    void setIsVisibleInEditor(bool b);
        //properties and lists
    PropertyGraphicsItem* getPropertyItem(const Property* prop);
    QList<PropertyGraphicsItem*> getPropertyItems();
    const QList<PropertyGraphicsItem*> getPropertyItems() const;
    QList<PropertyGraphicsItem*> getVisiblePropertyItems();
    const QList<PropertyGraphicsItem*> getVisiblePropertyItems() const;
    QList<PropertyGraphicsItem*> getAllPropertyItems();
    const QList<PropertyGraphicsItem*> getAllPropertyItems() const;
    QList<PropertyGraphicsItem*> getAllVisiblePropertyItems();
    const QList<PropertyGraphicsItem*> getAllVisiblePropertyItems() const;
        //modes
    void setPropertyVisibleMode(PropertyVisibleMode pvm);
    PropertyVisibleMode getPropertyVisibleMode() const;
    void setPropertyVisibleModifier(PropertyVisibleModeModifier pvmm);
    PropertyVisibleModeModifier getPropertyVisibleModeModifier() const;
    void setLinkArrowMode(LinkArrowVisibleMode lavm, PropertyListGraphicsItem* friendList = 0);
    LinkArrowVisibleMode getLinkArrowVisibleMode() const;
    void setHeaderMode(HeaderMode hm);
    HeaderMode getHeaderMode() const;
        //position
    void setPosition(const QPointF& pos);
    void setPosition(qreal x, qreal y);
        //width
    qreal getCurrentWidth() const;
    /**
     * Changes the scene of all PropertyItems and sub lists to the new scene.
     * It is used while opening/closing the PropertyLinkDialog.
     */
    void setScene(QGraphicsScene* scene);

    //----------------------------
    //  utils
    //----------------------------
    /** Passes change to properties and sub lists. */
    void prepareGeometryChange();
    /** Adds all properties and creates new items. */
    void addProperties(const QList<Property*>& items);
    /** Adds a sub list. */
    void addSubPropertyListItem(PropertyListGraphicsItem* list);
    /** Removes a sub list. */
    void removeSubPropertyListItem(PropertyListGraphicsItem* list);
    /** Used to return the current height for recursive layout. */
    float setChildListPosition(QPointF point);
    /** Sets the visiblility of all properties and sub lists. */
    void setVisible(bool visible);

    //----------------------------
    //  width management
    //----------------------------
    /** Find max width off all PropertyItems. */
    void findNewListWidth();
    /** Resize this and all sub lists. */
    void resizeListWidth(qreal width);

signals:
    /** Used for linking in the PropertyLinkDialog */
    void linkProperties(PropertyGraphicsItem* src, PropertyGraphicsItem* dst);

    //----------------------------
    //  members
    //----------------------------
private:
    PropertyOwnerGraphicsItem* parentPropertyOwnerItem_;    ///< PropertyOwnerGraphicsItem owning this list
    QGraphicsTextItem* headerLabelItem_;                    ///< label for the header
    QList<PropertyGraphicsItem*> propertyItems_;            ///< list of all PropertyGraphicsItems
    QList<PropertyListGraphicsItem*> subPropertyListItems_; ///< list of all sub-PropertyListGraphicsItems
    PropertyListGraphicsItem* superPropertyListItem_;       ///< pointer to the super list, if this is a sub list
    PropertyListGraphicsItem* friendPropertyListItem_;      ///< pointer to friend list used in PropertyLinkDialog

    HeaderMode currentHeaderMode_;                          ///< current mode, see enum definition
    PropertyVisibleMode currentPropertyMode_;               ///< current mode, see enum definition
    PropertyVisibleModeModifier currentPropertyModifier_;   ///< current mode, see enum definition
    LinkArrowVisibleMode currentLinkArrowMode_;             ///< current mode, see enum definition

    qreal currentHeight_;                                   ///< local value used in layoutChildren
    qreal currentWidth_;                                    ///< width of the this. Must be calculated

    bool isVisibleInEditor_;                                ///< stores, if this was visible in the NetworkEditor, if a PLDialog is being opened.
};

} // namespace

#endif // VRN_PROPERTYLISTGRAPHICSITEM_H
