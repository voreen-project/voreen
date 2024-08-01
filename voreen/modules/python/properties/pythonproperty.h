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

#ifndef VRN_PYTHONPROPERTY_H
#define VRN_PYTHONPROPERTY_H

#include "voreen/core/properties/templateproperty.h"
#include "modules/python/core/pythonscript.h"

namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API TemplateProperty<PythonScript>;
#endif

class VRN_CORE_API PythonProperty : public TemplateProperty<PythonScript> {
public:
    PythonProperty(const std::string& id, const std::string& guiText,
                   Processor::InvalidationLevel invalidationLevel=Processor::INVALID_PROGRAM,
                   Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    PythonProperty();
    ~PythonProperty();

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "PythonProperty"; }
    virtual std::string getTypeDescription() const { return "Python Script"; }

    /**
     * @see Property::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Property::deserialize
     */
    virtual void deserialize(Deserializer& s);

private:

    static const std::string loggerCat_;
};

}   // namespace

#endif
