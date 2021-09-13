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

#ifndef VRN_SAMPLEPOINTCONFIGLOADER_H
#define VRN_SAMPLEPOINTCONFIGLOADER_H

#include <string.h>

#include "voreen/core/io/serialization/serializable.h"
#include "tgt/vector.h"
#include <vector>

namespace voreen {


/**
 * Struct which represents a default samplepoint-category containing a name
 * and a default position.
 */
struct VRN_CORE_API SamplePointConfig : public Serializable {
    std::string name_;      // Name of the samplepoint category.
    tgt::vec3 defaultPos_;  // Default position of the samplepoints in this category.
    std::string textureFilename_; // Name of the file containing the texture for the region

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);
};

class VRN_CORE_API SamplePointConfigLoader {

public:
    std::vector<SamplePointConfig> loadSamplepointConfigFile(const std::string& filepath);

private:
    void trim(std::string& str);
};


    
} // END namespace vorren

#endif // VRN_SAMPLEPOINTCONFIGLOADER_H
