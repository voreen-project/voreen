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

#include "voreen/core/datastructures/volume/volumerepresentation.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/volumegl.h"
#include "voreen/core/datastructures/octree/volumeoctreebase.h"

using tgt::vec3;
using tgt::svec3;

namespace voreen {

const std::string VolumeRepresentation::loggerCat_("voreen.VolumeRepresentation");

VolumeRepresentation::VolumeRepresentation(const svec3& dimensions)
    : dimensions_(dimensions)
{
    numVoxels_ = (size_t) dimensions_.x * (size_t) dimensions_.y * (size_t) dimensions_.z;
}

svec3 VolumeRepresentation::getDimensions() const {
    return dimensions_;
}

size_t VolumeRepresentation::getNumVoxels() const {
    return numVoxels_;
}

//---------------------------------------------------------------------------------

ConverterFactory::ConverterFactory() {

    // Considering the memory hierarchy: Disk / Octree -> RAM -> GL
    // we add the following converters:

    // Disk -> RAM
    addConverter(new RepresentationConverterLoadFromDisk());

    // Octree -> RAM
    addConverter(new RepresentationConverterOctreeToRAM());

    // RAM -> GL
    addConverter(new RepresentationConverterUploadGL());

    // We explicitly don't (/no longer) allow a conversion in the other direction, since Disk representations
    // should never be deleted anyway, since they require virtually no volatile memory and GL representations
    // require an active GL context which is not guaranteed to be available.

    // The 'direct' conversion Disk / Octree -> GL is performed by the memory manager, which will use the RAM
    // representation, if available, or will create one from Disk or Octree, otherwise which is then uploaded.
}

ConverterFactory::~ConverterFactory() {
    while(!converters_.empty()) {
        delete converters_.back();
        converters_.pop_back();
    }
}
} // namespace voreen
