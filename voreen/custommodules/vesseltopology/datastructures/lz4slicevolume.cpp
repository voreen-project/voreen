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

#include "lz4slicevolume.h"

#include <fstream>

namespace voreen {

const static std::string METADATA_ROOT_NODE_STRING = "metadata";

LZ4SliceVolumeMetadata::LZ4SliceVolumeMetadata(tgt::svec3 dimensions)
    : dimensions_(dimensions)
{
}

LZ4SliceVolumeMetadata::LZ4SliceVolumeMetadata(const std::string& xmlfile)
    : dimensions_(tgt::svec3::one)
{
    XmlDeserializer ds;
    std::ifstream filestream(xmlfile);
    ds.read(filestream);
    ds.deserialize(METADATA_ROOT_NODE_STRING, *this);
}

void LZ4SliceVolumeMetadata::serialize(Serializer& s) const {
    s.serialize("dimensions", dimensions_);
}
void LZ4SliceVolumeMetadata::deserialize(Deserializer& s) {
    s.deserialize("dimensions", dimensions_);
}

void LZ4SliceVolumeMetadata::save(const std::string& xmlfile) const {
    XmlSerializer ser;
    std::ofstream filestream(xmlfile);
    ser.serialize(METADATA_ROOT_NODE_STRING, *this);
}

}
