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

#ifndef VRN_SERIALIZABLE_H
#define VRN_SERIALIZABLE_H

#include "voreen/core/voreencoreapi.h"

namespace voreen {

class Serializer;
class Deserializer;

/**
 * @c Serializable declares the interface for user defined data types
 * which should allow (de-)serialization.
 *
 * @par
 * Just realize the @c Serializable interface in your user defined
 * class to make it (de-)serializable using @c XmlSerializer and @c XmlDeserializer.
 *
 * @par
 * Here is a short example of realizing the @c Serializable interface
 * to make a class (de-)serializable:
 * @code
 * class Data : public Serializable {
 * private:
 *     int i_;
 *
 * public:
 *     virtual void serialize(Serializer& s) const {
 *         s.serialize("i", i_);
 *     }
 *
 *     virtual void deserialize(Deserializer& s) {
 *         s.deserialize("i", i_);
 *     }
 * };
 * @endcode
 * For more complex examples and interaction with other classes of the serialization framework
 * see the &quot;serializertest&quot; application in &quot;apps/&quot; directory.
 *
 * @note You have to implement a default constructor for each class which is used as a pointer
 *       for serialization. Moreover, you have to make it available to the @c XmlDeserializer class
 *       and to all factories which needs a default constructor using a @c friend statement.
 *
 * @attention You have to ensure that @c Serializable respectively @c AbstractSerializable
 *            is the first class from which you derive. Otherwise there will be pointer
 *            mismatches due to type conversions.
 *
 * @attention Do not use global respectively static variables in the @c serialize
 *            and @c deserialize method, since this can lead to difficult bugs.
 *            Furthermore, you may expect that @c serialize and @c deserialize
 *            methods are executed sequentially, but some of them need a delayed execution.
 *
 * @attention If you want to make an abstract class serializable, please realize
 *            @c AbstractSerializable instead of @c Serializable as it prevents trying to
 *            make instances using the default constructor.
 *
 * @see AbstractSerializable
 * @see XmlSerializer
 * @see XmlDeserializer
 */
class VRN_CORE_API Serializable {
public:
    virtual ~Serializable() {}

    /**
     * Serializes all member variables using the given @c XmlSerializer.
     *
     * @note You have to implement the serialization in your @c Serializable realization.
     *
     * @param s the @c XmlSerializer
     */
    virtual void serialize(Serializer& s) const = 0;

    /**
     * Deserializes all member varaibles using the given @c XmlDeserializer.
     *
     * @note You have to implement the deserialization in your @c Serializable realization.
     *
     * @param s the @c XmlDeserializer
     */
    virtual void deserialize(Deserializer& s) = 0;
};

} // namespace

#endif // VRN_SERIALIZABLE_H
