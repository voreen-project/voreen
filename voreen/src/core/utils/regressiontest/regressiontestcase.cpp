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

#include "voreen/core/utils/regressiontest/regressiontestcase.h"

#include "voreen/core/utils/stringutils.h"
#include "voreen/core/io/serialization/serializationexceptions.h"
#include "tgt/logmanager.h"
#include "tgt/vector.h"

namespace voreen {

std::ostream& operator<<(std::ostream &str, const RegressionTestFile& file) {
    str << file.filename_;
    return str;
}

bool RegressionTestCase::operator==(const RegressionTestCase& t2) const {
    return (moduleDir_ == t2.moduleDir_ &&
            group_ == t2.group_ &&
            name_ == t2.name_);
}

bool RegressionTestCase::operator<(const RegressionTestCase& t2) const {
    if (moduleDir_ != t2.moduleDir_)
        return moduleDir_ < t2.moduleDir_;
    else if (group_ != t2.group_)
        return group_ < t2.group_;
    else
        return name_ < t2.name_;
}

RegressionTestCaseConfiguration::RegressionTestCaseConfiguration()
    : Serializable()
    , enabled_(true)
    , ignored_(false)
    , timeout_(0)
    , pixelDiffTolerance_(-1.f)
    , maxErrorPixels_(-1)
    , pixelSearchNeighborhood_(-1)
    , voxelDiffTolerance_(-1.f)
    , maxErrorVoxels_(-1)
    , geometryDiffTolerance_(-1.f)
{}

void RegressionTestCaseConfiguration::serialize(Serializer& s) const  {
    s.serialize("enabled", enabled_);
    s.serialize("ignored", ignored_);
    s.serialize("timeout", timeout_);

    s.serialize("pixelDiffTolerance", pixelDiffTolerance_);
    s.serialize("maxErrorPixels", maxErrorPixels_);
    s.serialize("pixelSearchNeighborhood", pixelSearchNeighborhood_);

    s.serialize("voxelDiffTolerance", voxelDiffTolerance_);
    s.serialize("maxErrorVoxels", maxErrorVoxels_);

    s.serialize("geometryDiffTolerance", geometryDiffTolerance_);
}

void RegressionTestCaseConfiguration::deserialize(Deserializer& s) {
    try {
        try {
            s.deserialize("enabled", enabled_);
        }
        catch (SerializationNoSuchDataException& /*e*/) {}

        try {
            s.deserialize("ignored", ignored_);
        }
        catch (SerializationNoSuchDataException& /*e*/) {}

        try {
            s.deserialize("timeout", timeout_);
            if (timeout_ < 0) {
                LWARNINGC("regressiontest.RegressionTestCaseConfiguration", "Invalid timeout: " << timeout_);
                timeout_ = 0;
            }
        }
        catch (SerializationNoSuchDataException& /*e*/) {}

        try {
            s.deserialize("pixelDiffTolerance", pixelDiffTolerance_);
            if (pixelDiffTolerance_ < 0.f || pixelDiffTolerance_ > 1.f) {
                LWARNINGC("regressiontest.RegressionTestCaseConfiguration", "Invalid pixelDiffTolerance: " << pixelDiffTolerance_);
                pixelDiffTolerance_ = tgt::clamp(pixelDiffTolerance_, 0.f, 1.f);
            }
        }
        catch (SerializationNoSuchDataException& /*e*/) {}

        try {
            s.deserialize("maxErrorPixels", maxErrorPixels_);
            if (maxErrorPixels_ < 0) {
                LWARNINGC("regressiontest.RegressionTestCaseConfiguration", "Invalid maxErrorPixels: " << maxErrorPixels_);
                maxErrorPixels_ = 0;
            }
        }
        catch (SerializationNoSuchDataException& /*e*/) {}

        try {
            s.deserialize("pixelSearchNeighborhood", pixelSearchNeighborhood_);
            if (pixelSearchNeighborhood_ < 0) {
                LWARNINGC("regressiontest.RegressionTestCaseConfiguration", "Invalid pixelSearchNeighborhood: " << pixelSearchNeighborhood_);
                pixelSearchNeighborhood_ = 0;
            }
        }
        catch (SerializationNoSuchDataException& /*e*/) {}

        try {
            s.deserialize("voxelDiffTolerance", voxelDiffTolerance_);
            if (voxelDiffTolerance_ < 0.f || voxelDiffTolerance_ > 1.f) {
                LWARNINGC("regressiontest.RegressionTestCaseConfiguration", "Invalid voxelDiffTolerance: " << voxelDiffTolerance_);
                voxelDiffTolerance_ = tgt::clamp(voxelDiffTolerance_, 0.f, 1.f);
            }
        }
        catch (SerializationNoSuchDataException& /*e*/) {}

        try {
            s.deserialize("maxErrorVoxels", maxErrorVoxels_);
            if (maxErrorVoxels_ < 0) {
                LWARNINGC("regressiontest.RegressionTestCaseConfiguration", "Invalid maxErrorVoxels: " << maxErrorVoxels_);
                maxErrorVoxels_ = 0;
            }
        }
        catch (SerializationNoSuchDataException& /*e*/) {}

        try {
            s.deserialize("geometryDiffTolerance", geometryDiffTolerance_);
            if (geometryDiffTolerance_ < 0.f) {
                LWARNINGC("regressiontest.RegressionTestCaseConfiguration", "Invalid geometryDiffTolerance: " << geometryDiffTolerance_);
                geometryDiffTolerance_ = 0.f;
            }
        }
        catch (SerializationNoSuchDataException& /*e*/) {}

    }
    catch (SerializationException& e) {
        LERRORC("regressiontest.RegressionTestCaseConfiguration", "Deserialization failed: " << e.what());
    }
}

}
