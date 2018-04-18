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

#ifndef VRN_SERIALIZERBASE_H
#define VRN_SERIALIZERBASE_H

#include <vector>
#include <string>

#include "voreen/core/voreencoreapi.h"

namespace voreen {

// Implements common functionality for Voreens deserializers and serializers.
class VRN_CORE_API SerializerBase {
public:
    SerializerBase();

    /**
     * Sets whether to serialize pointer always as content instead of references.
     *
     * @attention This is not a cascading setting, which means that contained pointers
     *            are not serialized as content.
     *
     * @attention Serialization of all pointers as content can lead to redundant data.
     *
     * @param usePointerContentSerialization if @c true pointers are always serialized as content,
     *                                       otherwise as references when possible.
     */
    void setUsePointerContentSerialization(const bool& usePointerContentSerialization);

    /**
     * Returns whether pointers are always serialized as content instead of references.
     *
     * @return @c true if pointers are always serialized as content and @c false otherwise
     */
    bool getUsePointerContentSerialization() const;

    /**
     * Adds the given error @c message to the error list.
     *
     * @param message the error message
     */
    void addError(const std::string& message);

    /**
     * Adds the error message of the given @c exception to the error list.
     *
     * @param exception the exception
     */
    void addError(const std::exception& exception);

    /**
     * Removes the last error message from the error list.
     */
    void removeLastError();

    /**
     * Returns the error list.
     *
     * @return the error list
     */
    const std::vector<std::string>& getErrors() const;

    /**
     * Adds the error message from the given @c exception to the error list
     * and raise the exception afterwards.
     *
     * @tparam T exception type
     *
     * @param exception the exception
     *
     * @throws SerializationException the exception is always thrown
     */
    template<class T>
    void raise(const T& exception);

    /**
     * This is a helper class to ensure correct use pointer content serialization state.
     *
     * @note As C++ does not support a finally block statement, we need this
     *       class to ensure that cleanup code concerning the XML node
     *       for inserting or reading data is executed.
     */
    class VRN_CORE_API TemporaryUsePointerContentSerializationChanger {
    public:
        /**
         * Creates a @c TemporaryUsePointerContentSerializationChanger,
         * which changes the actual use pointer content serialization setting.
         *
         * @param serializer
         *     serializer whose use pointer content serialization setting should be changed
         * @param node the new setting
         */
        TemporaryUsePointerContentSerializationChanger(SerializerBase& serializer, const bool& usePointerContentSerialization);

        /**
         * Destructor ensures restoring the use pointer content serialization setting
         * which was set before this instance was created.
         */
        ~TemporaryUsePointerContentSerializationChanger();

    private:
        /**
         * Serializer whose use pointer content serialization setting is changed.
         */
        SerializerBase& serializer_;

        /**
         * Use pointer content setting which was set before
         * this @c TemporaryUsePointerContentSettingChanger was created.
         */
        const bool storedUsePointerContentSerialization_;
    };

protected:
    std::vector<std::string> errors_;
    bool usePointerContentSerialization_;
};


template<class T>
void SerializerBase::raise(const T& exception)
{
    addError(exception);
    throw exception;
}


} // namespace voreen

#endif // VRN_SERIALIZERBASE_H
