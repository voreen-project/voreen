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

#ifndef VRN_PROPERTYOWNERWIDGET_H
#define VRN_PROPERTYOWNERWIDGET_H

#include "voreen/core/properties/property.h"
#include "voreen/qt/voreenqtapi.h"
#include <QWidget>
#include <QVBoxLayout>

namespace voreen {

class ExpandableHeaderButton;
class QPropertyWidget;
class GroupPropertyWidget;
class PropertyOwner;

/**
 * Widget containing a PropertyOwner's property widgets.
 * The title bar contains the name of the PropertyOwner and an +/- icon allowing to expand the properties.
 * The widget is used in the "PropertyListWidget" and settings dialogs of the mainwindow.
 */
class VRN_QT_API PropertyOwnerWidget : public QWidget, public PropertyOwnerObserver {
Q_OBJECT
public:
    /**
     * Constructor
     * @param title if title is empty, the name of the propertyower is taken
     * @param expanded the initial state of the widget
     * @param userExpandable enables the +/- in the header
     * @param addResetButton adds a reset button at the and of the widget
     */
    PropertyOwnerWidget(PropertyOwner* propertyOwner, QWidget* parent = 0, std::string title = "",
                        bool expanded = false, bool userExpandable = true, bool addResetButton = false);
    /** PropertyWidgets are deleted by hand cause of context problems. */
    virtual ~PropertyOwnerWidget();

    /** Returns the associated property owner */
    PropertyOwner* getPropertyOwner() const;
    /** Returns the current expanded state */
    bool isExpanded() const;
    /** Returns, if the header contains a +/- */
    bool isUserExpandable() const;

    /** PropertyOwner function. Updates the properties */
    virtual void propertiesChanged(const PropertyOwner*);

public slots:
    /** Changes, if the propertyWidget is expanded or not */
    void setExpanded(bool expanded);
    /** Changes, if the header should have a +/- */
    void setUserExpandable(bool expandable);
    /** Changes the expandation state */
    void toggleExpansionState();
    /** Updates the title (title is taken from propertyowner) */
    void updateHeaderTitle();
    /** Creates the widget and adds all propertywidgets to it */
    void instantiateWidgets();
    /** Updates the geoemtry. Should be called, if expadation changed or the property number */
    void updateState();
    /** Sets the level of detail */
    void setLevelOfDetail(Property::LevelOfDetail lod);
protected slots:
    /** Sets all property values of the property owner to their default state */
    void resetAllProperties();
    /** Updates the visibility of property groups. */
    void updateGroupVisibility();
    /** Sets the correct proprety visibility */
    void updatePropertyVisibility();
protected:
    /** Calls initialzeWidgets on the first call */
    virtual void showEvent(QShowEvent*);
    /** Get all properties of the property owner */
    virtual std::vector<Property*>* createPropertyList();
    //------------------
    //     member
    //------------------
    PropertyOwner* propertyOwner_;      ///< the property owner associated with this widget

    QVBoxLayout* mainLayout_;           ///< the main layout containing the header and propertywidget
    ExpandableHeaderButton* header_;    ///< customized header with +/- for expandation
    QWidget* propertyWidget_;           ///< widget containing all propertywidgets. visible, if the header is expanded

    std::map<std::string, GroupPropertyWidget*> propertyGroupsMap_; ///< map containing all property groups
    std::vector<QPropertyWidget*> widgets_;     ///< containing all created propertywidgets

    bool hasBeenInitialized_;           ///< determines, if the propertywidgets have been initialized
    bool addResetButton_;               ///< set in the constructor. If true, a reset button is added in "widgetInstantiation"
    Property::LevelOfDetail currentLOD_;///< the current LOD. Used in update Properties.
};

} // namespace

#endif
