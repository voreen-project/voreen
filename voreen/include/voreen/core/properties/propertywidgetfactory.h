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

#ifndef VRN_PROPERTYWIDGETFACTORY_H
#define VRN_PROPERTYWIDGETFACTORY_H

#include "voreen/core/voreencoreapi.h"

namespace voreen {

class Property;
class PropertyWidget;

/**
 * Interface for factories that create GUI representations for properties.
 *
 * If a module contains custom properties, it is also expected to provide
 * a suitable PropertyWidgetFactory for these properties.
 *
 * @see VoreenModule::registerPropertyWidgetFactory
 */
class VRN_CORE_API PropertyWidgetFactory {
public:
    /** Constructor */
    PropertyWidgetFactory() {}
    /** Destructor */
    virtual ~PropertyWidgetFactory() {}

    /**
     * Creates the associated PropertyWidget to the passed Property.
     * Calls createAssociatedWidget internaly, which should be implemented in
     * derived classes. If createAssociatedWidget returns a non null value,
     * - initialize
     * - ...
     * will be called.
     */
    PropertyWidget* createWidget(Property* prop) const;

    //----------------------------
    // Functions to Override
    //----------------------------
protected:
    /**
     * Has to be implemented in sub-class.
     */
    virtual PropertyWidget* createAssociatedWidget(Property*) const = 0;
};

} // namespace

#endif // VRN_PROPERTYWIDGETFACTORY_H
