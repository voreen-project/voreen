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

#include "voreen/qt/networkeditor/graphicitems/core/propertyownergraphicsitem.h"

//core
#include "voreen/core/properties/propertyowner.h"


namespace voreen {

PropertyOwnerGraphicsItem::PropertyOwnerGraphicsItem(PropertyOwner* propertyOwner, NetworkEditor* networkEditor)
    : NWEBaseGraphicsItem(networkEditor)
    , propertyOwner_(propertyOwner)
    , guiName_(QString::fromStdString("NONAME"))
    , propertyList_(this)
{
    propertyList_.hide();
}

PropertyOwnerGraphicsItem::~PropertyOwnerGraphicsItem() {
}

//---------------------------------------------------------------------------------------------------------------
//                  getter and setter
//---------------------------------------------------------------------------------------------------------------
PropertyOwner* PropertyOwnerGraphicsItem::getPropertyOwner() const {
    return propertyOwner_;
}

PropertyListGraphicsItem* PropertyOwnerGraphicsItem::getPropertyList() {
    return &propertyList_;
}

void PropertyOwnerGraphicsItem::setGuiName(QString str) {
    guiName_ = str;
}

QString PropertyOwnerGraphicsItem::getGuiName() const {
    return guiName_;
}

//---------------------------------------------------------------------------------------------------------------
//                  property functions
//---------------------------------------------------------------------------------------------------------------
bool PropertyOwnerGraphicsItem::hasProperty(const Property* prop) const {
    if(!propertyOwner_) return false;
    const std::vector<Property*> properties = propertyOwner_->getProperties();
    for (size_t i = 0; i < properties.size(); ++i)
        if (properties[i] == prop)
            return true;
    return false;
}

} // namespace
