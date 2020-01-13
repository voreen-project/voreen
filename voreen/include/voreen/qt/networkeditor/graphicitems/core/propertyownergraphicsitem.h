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

#ifndef VRN_PROPERTYOWNERGRAPHICSITEM_H
#define VRN_PROPERTYOWNERGRAPHICSITEM_H

#include "../nwebasegraphicsitem.h"

#include "propertylistgraphicsitem.h"

namespace voreen {

class PropertyOwner;
class Property;
class PropertyGraphicsItem;
class PropertyListGraphicsItem;

/**
 * Qt graphic item class corresponding to voreen/core/propertyowner.
 */
class PropertyOwnerGraphicsItem : public NWEBaseGraphicsItem {
Q_OBJECT
public:
    //constructor + destructor
    PropertyOwnerGraphicsItem(PropertyOwner* propertyOwner, NetworkEditor* networkEditor);
    virtual ~PropertyOwnerGraphicsItem();

    //getter and setter
    PropertyOwner* getPropertyOwner() const;

    void setGuiName(QString str);
    QString getGuiName() const;

    PropertyListGraphicsItem* getPropertyList();

    //property functions
    virtual bool hasProperty(const Property* prop) const;

protected:
    //members
    PropertyOwner* propertyOwner_;  ///can be null
    QString guiName_;                   ///to get the guiname without casting
    PropertyListGraphicsItem propertyList_;      ///list to store the properties
};

} // namespace

#endif // VRN_PROPERTYOWNERGraphicsItem_H
