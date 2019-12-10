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

#include "vvdvolumereader.h"
#include "vvdformat.h"

#include <fstream>
#include <iostream>

#include "tgt/exception.h"
#include "tgt/vector.h"
#include "tgt/filesystem.h"

#include "voreen/core/utils/stringutils.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumedisk.h"

using tgt::vec3;
using tgt::ivec3;
using tgt::hor;
using tgt::lessThanEqual;

namespace voreen {

const std::string VvdVolumeReader::loggerCat_ = "voreen.io.VolumeReader.vvd";

VvdVolumeReader::VvdVolumeReader(ProgressBar* progress)
    : VolumeReader(progress)
{
    extensions_.push_back("vvd");
}

VolumeList* VvdVolumeReader::read(const std::string &url) {
    VolumeURL origin(url);
    std::string fileName = origin.getPath();

    // open file for reading
    std::fstream fileStream(fileName.c_str(), std::ios_base::in);
    if (fileStream.fail()) {
        throw tgt::FileException("Failed to open file '" + tgt::FileSystem::absolutePath(fileName) + "' for reading.");
    }

    // read data stream into deserializer
    XmlDeserializer d(fileName);
    d.setUseAttributes(true);
    try {
        d.read(fileStream);
    }
    catch (SerializationException& e) {
        throw tgt::FileException("SerializationException: Failed to read serialization data stream from file '"
                                     + fileName + "': " + e.what());
    }
    catch (...) {
        throw tgt::FileException("Failed to read serialization data stream from file '"
                                     + fileName + "' (unknown exception).");
    }

    std::vector<VvdObject> vec;
    Deserializer deserializer(d);
    // deserialize from data stream
    try {
        deserializer.deserialize("Volumes", vec, "Volume");
    }
    catch (std::exception& e) {
        throw tgt::FileException("Deserialization from file '" + fileName + "' failed: " + e.what());
    }
    catch (...) {
        throw tgt::FileException("Deserialization from file '" + fileName + "' failed (unknown exception).");
    }
    if (vec.empty())
        throw tgt::FileException("Deserialization from file '" + fileName + "' failed: no VvdObject deserialized");

    std::vector<std::unique_ptr<VolumeBase>> vc;
    for(size_t i=0; i<vec.size(); i++) {
        std::unique_ptr<Volume> vh(vec[i].createVolume());

        // If the volume does not have a RealWorldMapping, we add a denormalizing mapping for its data type.
        // As previously the default RealWorldMapping had a scale of one, and all DerivedData created at that
        // time used the default RealWorldMapping, we clear all DerivedData and thus force them to be recreated.
        if(!vh->hasMetaData(VolumeBase::META_DATA_NAME_REAL_WORLD_MAPPING)) {
            RealWorldMapping rwm;
            try {
                rwm = RealWorldMapping::createDenormalizingMapping(vh->getBaseType());
            } catch(tgt::UnsupportedFormatException&) {
                throw tgt::FileException("Deserialization from file '" + fileName + "' failed: Unsupported base type.");
            }

            // Only clear derived data, if real world mapping differs.
            if(rwm != vh->getRealWorldMapping()) {
                vh->clearDerivedData();
            }

            // Always set real world mapping to explicitly tell the volume has the default real world mapping.
            vh->setRealWorldMapping(rwm);
        }

        vh->setOrigin(origin);
        vc.push_back(std::move(vh));
    }

    VolumeList* list = new VolumeList();
    for(auto& volume : vc) {
        list->add(volume.release());
    }

    return list;
}

VolumeReader* VvdVolumeReader::create(ProgressBar* progress) const {
    return new VvdVolumeReader(progress);
}

} // namespace voreen
