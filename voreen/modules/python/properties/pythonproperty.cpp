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

#include "pythonproperty.h"

#include "voreen/core/utils/voreenfilepathhelper.h"
#include "tgt/filesystem.h"

namespace voreen {

const std::string PythonProperty::loggerCat_("voreen.python.PythonProperty");

//---------------------------------------------------------------------------------------------------------------

PythonProperty::PythonProperty(const std::string& id, const std::string& guiText, Processor::InvalidationLevel invalidationLevel, Property::LevelOfDetail lod)
    : TemplateProperty<PythonScript>(id, guiText, PythonScript(), invalidationLevel, lod)
{
}

PythonProperty::PythonProperty()
{
}

PythonProperty::~PythonProperty()
{
}

Property* PythonProperty::create() const {
    return new PythonProperty();
}

void PythonProperty::serialize(Serializer& s) const {
    Property::serialize(s);
    get().serialize(s);
}

void PythonProperty::deserialize(Deserializer& s) {
    Property::deserialize(s);

    PythonScript n = get();
    n.deserialize(s);
    set(n);

    invalidate();
    updateWidgets();
}


}   // namespace
