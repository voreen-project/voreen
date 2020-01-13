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

#ifndef VRN_SERIALIZATIONEXCEPTIONS_H
#define VRN_SERIALIZATIONEXCEPTIONS_H

#include "voreen/core/utils/exception.h"

namespace voreen {

/**
 * Base class of all serialization exceptions.
 *
 * @note You should rather raise a derived exception,
 *       which is more specific to the error that occured,
 *       than a @c SerializationException.
 *
 * @see VoreenException
 */
class SerializationException : public VoreenException {
public:
    /**
     * @see VoreenException::VoreenException
     */
    SerializationException(const std::string& what = "") : VoreenException(what) {}
};

//----------------------------------------------------------------------------

/**
 * A @c SerializationInvalidOperationException is raised in case of a method-call in an
 * invalid, respectively unsupported way.
 *
 * @see SerializationException
 */
class SerializationInvalidOperationException : public SerializationException {
public:
    /**
     * @see SerializationException::SerializationException
     */
    SerializationInvalidOperationException(const std::string& what = "") : SerializationException(what) {}
};

//----------------------------------------------------------------------------

/**
 * A @c SerializationFormatException is raised when XML nodes does not fulfill expected format.
 *
 * @see SerializationException
 */
class SerializationFormatException : public SerializationException {
public:
    /**
     * @see SerializationException::SerializationException
     */
    SerializationFormatException(const std::string& what = "") : SerializationException(what) {}
};

//----------------------------------------------------------------------------

/**
 * A @c SerializationVersionMismatchException is raised in case of a version mismatch among
 * the XML document and the used @c XmlSerializer or @c XmlDeserializer.
 *
 * @see SerializationException
 */
class SerializationVersionMismatchException : public SerializationException {
public:
    /**
     * @see SerializationException::SerializationException
     */
    SerializationVersionMismatchException(const std::string& what = "") : SerializationException(what) {}
};

//----------------------------------------------------------------------------

/**
 * A @c SerializationNoSuchDataException is raised in case of searching for a XML node
 * by key that does not exists.
 *
 * @see SerializationException
 */
class SerializationNoSuchDataException : public SerializationException {
public:
    /**
     * @see SerializationException::SerializationException
     */
    SerializationNoSuchDataException(const std::string& what = "") : SerializationException(what) {}
};

//----------------------------------------------------------------------------

/**
 * A @c SerializationDuplicateIdException is raised in case of multiple XML nodes
 * share the same id attribute.
 *
 * @see SerializationException
 */
class SerializationDuplicateIdException : public SerializationException {
public:
    /**
     * @see SerializationException::SerializationException
     */
    SerializationDuplicateIdException(const std::string& what = "") : SerializationException(what) {}
};

//----------------------------------------------------------------------------

/**
 * A @c SerializationAttributeNamingException is raised in case of multiple XML attributes
 * sharing the same name or using a reserved attribute name.
 *
 * @see SerializationException
 */
class SerializationAttributeNamingException : public SerializationException {
public:
    /**
     * @see SerializationException::SerializationException
     */
    SerializationAttributeNamingException(const std::string& what = "") : SerializationException(what) {}
};

//----------------------------------------------------------------------------

/**
 * A @c SerializationReferenceResolvingException is raised in case of problems
 * concerning the reference resolving process.
 *
 * @see SerializationException
 */
class SerializationReferenceResolvingException : public SerializationException {
public:
    /**
     * @see SerializationException::SerializationException
     */
    SerializationReferenceResolvingException(const std::string& what = "")
        : SerializationException(what) {}
};

//----------------------------------------------------------------------------

/**
 * A @c SerializationMemoryAllocationException is raised in case of trying to allocate memory
 * for an @c AbstractSerializable.
 *
 * @see SerializationException
 */
class SerializationMemoryAllocationException : public SerializationException {
public:
    SerializationMemoryAllocationException(const std::string& what = "") : SerializationException(what) {}
};

} // namespace

#endif // VRN_SERIALIZATIONEXCEPTIONS_H
