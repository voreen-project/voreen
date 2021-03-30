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

#include "vvdvolumewriter.h"
#include "vvdformat.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/volume.h"

#include "tgt/filesystem.h"
#include "tgt/matrix.h"

namespace voreen {

const std::string VvdVolumeWriter::loggerCat_("voreen.io.VvdVolumeWriter");

VvdVolumeWriter::VvdVolumeWriter() {
    extensions_.push_back("vvd");
}

void VvdVolumeWriter::write(const std::string& filename, const VolumeBase* volumeHandle) {
    tgtAssert(volumeHandle, "No volume");

    if(!volumeHandle->hasRepresentation<VolumeDisk>() && !volumeHandle->hasRepresentation<VolumeRAM>()) {
        LWARNING("Neither volume disk nor RAM representation available");
        return;
    }

    std::string vvdname = tgt::FileSystem::cleanupPath(filename);
    std::string rawname = tgt::FileSystem::fullBaseName(vvdname) + ".raw";
    LINFO("saving " << vvdname << " and " << rawname);

    // VVD: ---------------------------

    XmlSerializer s(vvdname);
    s.setUseAttributes(true);

    VvdObject o = VvdObject(volumeHandle, rawname);
    std::vector<VvdObject> vec;
    vec.push_back(o);

    Serializer serializer(s);
    serializer.serialize("Volumes", vec, "Volume");

    //errorList_ = s.getErrors();

    // write serialization data to temporary string stream
    std::ostringstream textStream;
    try {
        s.write(textStream);
        if (textStream.fail())
            throw tgt::IOException("Failed to write serialization data to string stream.");
    }
    catch (std::exception& e) {
        throw tgt::IOException("Failed to write serialization data to string stream: " + std::string(e.what()));
    }
    catch (...) {
        throw tgt::IOException("Failed to write serialization data to string stream (unknown exception).");
    }

    // Now we have a valid StringStream containing the serialization data.
    // => Open output file and write it to the file.
    std::fstream fileStream(vvdname.c_str(), std::ios_base::out);
    if (fileStream.fail())
        throw tgt::IOException("Failed to open file '" + vvdname + "' for writing.");

    try {
        fileStream << textStream.str();
    }
    catch (std::exception& e) {
        throw tgt::IOException("Failed to write serialization data stream to file '"
                                     + vvdname + "': " + std::string(e.what()));
    }
    catch (...) {
        throw tgt::IOException("Failed to write serialization data stream to file '"
                                     + vvdname + "' (unknown exception).");
    }
    fileStream.close();

    std::fstream rawout(rawname.c_str(), std::ios::out | std::ios::binary);

    if (!rawout.is_open() || rawout.bad())
        throw tgt::IOException();

    // RAW: ---------------------------

    auto writeToDisk = [&rawout] (const VolumeRAM* volume) {
        const char* data = static_cast<const char*>(volume->getData());
        size_t numbytes = volume->getNumVoxels() * volume->getBytesPerVoxel();
        rawout.write(data, numbytes);
    };

    if (volumeHandle->hasRepresentation<VolumeRAM>()) {
        VolumeRAMRepresentationLock v(volumeHandle);
        tgtAssert(*v, "no volume");
        writeToDisk(*v);
    }
    else if (volumeHandle->hasRepresentation<VolumeDisk>()) {
        const VolumeDisk* volumeDisk = volumeHandle->getRepresentation<VolumeDisk>();
        tgtAssert(volumeDisk, "no disk volume");

        size_t numSlices = volumeHandle->getDimensions().z;
        tgtAssert(numSlices > 0, "empty volume");
        for (size_t slice=0; slice<numSlices; slice++) {
            try {
                std::unique_ptr<VolumeRAM> sliceVolume(volumeDisk->loadSlices(slice, slice));
                writeToDisk(sliceVolume.get());
            }
            catch (tgt::Exception& e) {
                LWARNING("Error writing slice to disk: failed to load slice from disk volume: " << e.what());
                break;
            }
        }
    }

    if (rawout.bad())
        throw tgt::IOException("Failed to write volume data to file (bad stream)");

    rawout.close();
}

VolumeWriter* VvdVolumeWriter::create(ProgressBar* /*progress*/) const {
    return new VvdVolumeWriter();
}

} // namespace voreen
