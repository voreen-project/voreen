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

#ifndef VRN_STRINGPROPERTY_H
#define VRN_STRINGPROPERTY_H

#include "voreen/core/properties/templateproperty.h"

namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API TemplateProperty<std::string>;
#endif

class VRN_CORE_API StringProperty : public TemplateProperty<std::string> {
public:
    StringProperty(const std::string& id, const std::string& guiText, const std::string& value = "",
        int invalidationLevel=Processor::INVALID_RESULT, Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    StringProperty();
    virtual ~StringProperty() {}

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "StringProperty"; }
    virtual std::string getTypeDescription() const { return "String"; }

    /**
     * @see Property::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Property::deserialize
     */
    virtual void deserialize(Deserializer& s);

    //! @see Property::setReadOnlyFlag
    void setReadOnly(bool readOnly);

    //! @see Property::isReadyOnlyFlagSet
    bool isReadOnly() const;

    void setInstantUpdate(bool instantUpdate);
    bool getInstantUpdate() const;

protected:
    bool instantUpdate_;

};

}

#endif
