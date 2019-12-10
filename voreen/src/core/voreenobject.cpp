/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "voreen/core/voreenobject.h"

#include "tgt/logmanager.h"

namespace voreen {

const std::string VoreenObject::loggerCat_("voreen.VoreenObject");

VoreenObject::VoreenObject() {
    //tgtAssert(false, "VoreenObject constr.");
}

VoreenObject::VoreenObject(const std::string& id, const std::string& guiName)
    : id_(id)
    , guiName_(guiName)
{
    /*if (id_.empty() && tgt::Singleton<tgt::LogManager>::isInited())
        LWARNING("Empty id string passed to constructor"); */
}

std::string VoreenObject::getID() const {
    return id_;
}

std::string VoreenObject::getGuiName() const {
    return guiName_;
}

void VoreenObject::setGuiName(const std::string& guiName) {
    guiName_ = guiName;
}

//-----------------------------------------------------------------------------

VoreenSerializableObject::VoreenSerializableObject()
    : VoreenObject()
{}

VoreenSerializableObject::VoreenSerializableObject(const std::string& id, const std::string& guiName)
    : VoreenObject(id, guiName)
{}

} // namespace
