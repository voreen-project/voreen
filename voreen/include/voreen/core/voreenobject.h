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

#ifndef VRN_VOREENOBJECT_H
#define VRN_VOREENOBJECT_H

#include "voreen/core/voreencoreapi.h"
#include "voreen/core/io/serialization/abstractserializable.h"

#include <string>

namespace voreen {

/**
 * Base class for all Voreen types.
 *
 * Provides identifier handling.
 */
class VRN_CORE_API VoreenObject {

public:
    VoreenObject();

    /**
     * Standard constructor expecting the instance's id and GUI name.
     *
     * @param id identifier of the instance. Must not be empty and should not contain spaces.
     * @param guiName instance name that should be presentable to the user. May be empty.
     */
    VoreenObject(const std::string& id, const std::string& guiName);

    virtual ~VoreenObject() {}

    /**
     * Returns the instance's identifier.
     *
     * @note Overriding this function is usually not necessary.
     */
    virtual std::string getID() const;

    /**
     * Returns the instance's GUI name, which should be presentable to the user.
     *
     * @note Overriding this function is usually not necessary.
     */
    virtual std::string getGuiName() const;

    /**
     * Sets the instance's GUI name, which should be presentable to the user.
     *
     * @note Overriding this function is usually not necessary.
     */
    virtual void setGuiName(const std::string& guiName);

protected:
    std::string id_;         //< the instance's id
    std::string guiName_;    //< the instance's GUI name

    static const std::string loggerCat_;
};

//-----------------------------------------------------------------------------

/**
 * Abstract base class for serializable VoreenObjects, i.e., VoreenObject types
 * that can be handled by the serialization framework.
 *
 * @see Serializable
 */
class VRN_CORE_API VoreenSerializableObject : public AbstractSerializable, public VoreenObject {

public:
    VoreenSerializableObject();

    /**
     * Standard constructor expecting the instance's id and GUI name.
     *
     * @param id identifier of the instance. Must not be empty and should not contain spaces.
     * @param guiName instance name that should be presentable to the user. May be empty.
     */
    VoreenSerializableObject(const std::string& id, const std::string& guiName);

    /**
     * Virtual constructor: supposed to return an instance of the concrete subclass.
     */
    virtual VoreenSerializableObject* create() const = 0;

    /**
     * Returns the type name of the concrete subclass as string.
     * Necessary due to the lack of code reflection in C++.
     *
     * Is expected to be re-implemented by each concrete subclass.
     */
    virtual std::string getClassName() const = 0;
};

} // namespace

#endif // VRN_VOREENOBJECT_H
