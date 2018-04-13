/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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
 * Base class of all serialization exceptions concerning reading and writing XML.
 *
 * @note You should rather raise a derived exception,
 *       which is more specific to the error that occured,
 *       than a @c XmlSerializationException.
 *
 * @see SerializationException
 */
class XmlSerializationException : public SerializationException {
public:
    /**
     * @see SerializationException::SerializationException
     */
    XmlSerializationException(const std::string& what = "") : SerializationException(what) {}
};

//----------------------------------------------------------------------------

/**
 * A @c XmlSerializationInvalidOperationException is raised in case of a method-call in an
 * invalid, respectively unsupported way.
 *
 * @see XmlSerializationException
 */
class XmlSerializationInvalidOperationException : public XmlSerializationException {
public:
    /**
     * @see XmlSerializationException::XmlSerializationException
     */
    XmlSerializationInvalidOperationException(const std::string& what = "") : XmlSerializationException(what) {}
};

//----------------------------------------------------------------------------

/**
 * A @c XmlSerializationFormatException is raised when XML nodes does not fulfill expected format.
 *
 * @see XmlSerializationException
 */
class XmlSerializationFormatException : public XmlSerializationException {
public:
    /**
     * @see XmlSerializationException::XmlSerializationException
     */
    XmlSerializationFormatException(const std::string& what = "") : XmlSerializationException(what) {}
};

//----------------------------------------------------------------------------

/**
 * A @c XmlSerializationVersionMismatchException is raised in case of a version mismatch among
 * the XML document and the used @c XmlSerializer or @c XmlDeserializer.
 *
 * @see XmlSerializationException
 */
class XmlSerializationVersionMismatchException : public XmlSerializationException {
public:
    /**
     * @see XmlSerializationException::XmlSerializationException
     */
    XmlSerializationVersionMismatchException(const std::string& what = "") : XmlSerializationException(what) {}
};

//----------------------------------------------------------------------------

/**
 * A @c XmlSerializationNoSuchDataException is raised in case of searching for a XML node
 * by key that does not exists.
 *
 * @see XmlSerializationException
 */
class XmlSerializationNoSuchDataException : public XmlSerializationException {
public:
    /**
     * @see XmlSerializationException::XmlSerializationException
     */
    XmlSerializationNoSuchDataException(const std::string& what = "") : XmlSerializationException(what) {}
};

//----------------------------------------------------------------------------

/**
 * A @c XmlSerializationDuplicateIdException is raised in case of multiple XML nodes
 * share the same id attribute.
 *
 * @see XmlSerializationException
 */
class XmlSerializationDuplicateIdException : public XmlSerializationException {
public:
    /**
     * @see XmlSerializationException::XmlSerializationException
     */
    XmlSerializationDuplicateIdException(const std::string& what = "") : XmlSerializationException(what) {}
};

//----------------------------------------------------------------------------

/**
 * A @c XmlSerializationAttributeNamingException is raised in case of multiple XML attributes
 * sharing the same name or using a reserved attribute name.
 *
 * @see XmlSerializationException
 */
class XmlSerializationAttributeNamingException : public XmlSerializationException {
public:
    /**
     * @see XmlSerializationException::XmlSerializationException
     */
    XmlSerializationAttributeNamingException(const std::string& what = "") : XmlSerializationException(what) {}
};

//----------------------------------------------------------------------------

/**
 * A @c XmlSerializationReferenceResolvingException is raised in case of problems
 * concerning the reference resolving process.
 *
 * @see XmlSerializationException
 */
class XmlSerializationReferenceResolvingException : public XmlSerializationException {
public:
    /**
     * @see XmlSerializationException::XmlSerializationException
     */
    XmlSerializationReferenceResolvingException(const std::string& what = "")
        : XmlSerializationException(what) {}
};

//----------------------------------------------------------------------------

/**
 * A @c XmlSerializationMemoryAllocationException is raised in case of trying to allocate memory
 * for an @c AbstractSerializable.
 *
 * @see XmlSerializationException
 */
class XmlSerializationMemoryAllocationException : public XmlSerializationException {
public:
    XmlSerializationMemoryAllocationException(const std::string& what = "") : XmlSerializationException(what) {}
};

} // namespace

#endif // VRN_SERIALIZATIONEXCEPTIONS_H
