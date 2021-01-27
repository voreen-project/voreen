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

#ifndef VRN_VOREENFILEPATHHELPER_H
#define VRN_VOREENFILEPATHHELPER_H

#include "voreen/core/io/serialization/serializable.h"
#include "voreen/core/utils/exception.h"
#include "voreen/core/datastructures/volume/volumeurl.h"

#include <string>

namespace voreen {

    /**
     * This class should be used to serialize file paths.
     * The path will be serialized in different ways to make the deserialization flexible.
     * The path will be serialized:
     * - relative to the document
     * - relative to the voreen file
     * - absolut
     * During the deserialization the serialized path containing a valid file will be used.
     */
class VoreenFilePathHelper : public Serializable {
public:
    /**
     * Constructor
     * @param absFilePath The passed string should be the absolut path of the file.
     */
    VoreenFilePathHelper(const std::string& absFilePath = "");
    ~VoreenFilePathHelper();

    const std::string& getPath() const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

private:
    std::string path_;
};



} // namespace

#endif // VRN_VOREENFILEPATHHELPER_H
