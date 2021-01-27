/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_PROPERTYWIDGET_H
#define VRN_PROPERTYWIDGET_H

#include "voreen/core/properties/property.h"
#include "voreen/core/voreencoreapi.h"

namespace voreen {

/**
 * A PropertyWidget is a graphical Representation for a Property. It
 * allows manipulation of the Property's value.
 */
class VRN_CORE_API PropertyWidget {
public:
    PropertyWidget(Property* prop) : prop_(prop) {}
    virtual ~PropertyWidget() {}

    /**
     * The Widget should update itself from the Property's value when this is called.
     */
    virtual void updateFromProperty() = 0;

    virtual Property* getProperty() { return prop_; }

    /**
     * Sets the widgets enabled state. Disabled widgets can not be modified
     * by mouse or keyboard input.
     */
    virtual void setEnabled(bool enabled) = 0;

    /**
     * Toggles the widget's visibility state.
     */
    virtual void setVisible(bool state) = 0;

    /**
     * Updates the Gui Name
     */
    virtual void setPropertyGuiName(std::string) = 0;

    /**
     * Returns the Gui Name
     */
    virtual std::string getPropertyGuiName() = 0;

    /**
     * Updates the widget's visibility and read only state.
     */
    virtual void updateViewFlags(Property::ViewFlags flags) = 0;

    /**
     * When the Property is destroyed it calls disconnect so the Widget knows
     * the property is no longer around
     */
    virtual void disconnect() = 0;

    /**
     * Allows property widgets to serialize their state.
     *
     * This function is called by the Property's serialization routine.
     * If not overwritten in a subclass, it returns 0.
     */
    virtual MetaDataBase* getWidgetMetaData() const {
        return 0;
    }

    /**
     * Updates meta data containing property widget state for serialization.
     *
     * @note This function is called during property serialization. If not overwritten
     *       in a subclass it just does nothing, assuming that there is no meta data.
     */
    virtual void updateMetaData() const {}

protected:
    friend class PropertyWidgetFactory;
    /**
     * This function is called by PropertyWidgetFatory.
     * It can be used to initialize values which can not be initialized in the
     * constructor. An example is the call of virtual functions.
     */
    virtual void initialize() {}

    Property* prop_;
};


} // namespace voreen

#endif // VRN_PROPERTYWIDGET_H
