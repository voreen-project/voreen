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

#ifdef VRN_MODULE_FLOWREEN

#include "modules/flowreen/datastructures/deprecated/flow3d.h"
#include "modules/flowreen/utils/flowmath.h"
#include "flowreader.h"
#include "modules/flowreen/datastructures/deprecated/volumeflow3d.h"

#include <limits>

namespace voreen {

const std::string FlowReader::loggerCat_("voreen.io.FlowReader");

FlowReader::FlowReader(ProgressBar* const progress)
    : VolumeReader(progress)
{
    extensions_.push_back("flow");
}

VolumeList* FlowReader::read(const std::string& url) {
    VolumeURL origin(url);
    std::string fileName = origin.getPath();

    LINFO("reading flow file '" << fileName << "'...");

    // try to open the file
    //
    std::fstream ifs(fileName.c_str(), std::ios_base::in | std::ios_base::binary);
    if (ifs.good() == false)
        throw tgt::IOException("Unable to open flow file for reading", fileName);

    // read the magic number (string "VOREENFLOW")
    //
    char magicNumber[11] = {0};
    ifs.read(magicNumber, 11);
    std::string temp(magicNumber);
    if (temp != "VOREENFLOW")
        throw tgt::IOException("Missing magic number in flow file", fileName);

    // read file version (must currently be 1 or 2)
    //
    unsigned int fileVersion = 0;
    ifs.read(reinterpret_cast<char*>(&fileVersion), sizeof(unsigned int));
    LINFO("file version: " << fileVersion);
    if ((fileVersion < 1) || (fileVersion > 2))
        throw tgt::IOException("Unsupported file version of flow file", fileName);

    // read flow dimension (usually 3 for 3D flows)
    //
    unsigned int flowDimension = 0;
    ifs.read(reinterpret_cast<char*>(&flowDimension), sizeof(unsigned int));
    LINFO("flow dimension: " << flowDimension << "D");
    if (flowDimension != 3)
        throw tgt::IOException("Unsupported flow dimension in flow file", fileName);

    unsigned char dataOrientation = 0;
    unsigned char reverseSlicesMask = 0;

    ifs.read(reinterpret_cast<char*>(&dataOrientation), sizeof(unsigned char));
    if (fileVersion > 1) {
        ifs.read(reinterpret_cast<char*>(&reverseSlicesMask), sizeof(unsigned char));
    }

    // read the dimension of the volume data containing the flow
    //
    tgt::ivec3 dimensions;
    ifs.read(reinterpret_cast<char*>(&dimensions), sizeof(tgt::ivec3));
    LINFO("volume dimensions: " << dimensions);

    unsigned int byteSize = 0;
    ifs.read(reinterpret_cast<char*>(&byteSize), sizeof(unsigned int));
    LINFO("expected size of vector field: " << byteSize << " byte");

    VolumeFlow3D* volume = readConvert(dimensions, dataOrientation, ifs);
    ifs.close();

    if (volume == 0) {
        LERROR("an error occured during reading flow data! Proceeding impossible.");
        return 0;
    }

    if (reverseSlicesMask != 0)
        reverseSlices(volume, reverseSlicesMask);

    // TODO: volume container merge
    /*VolumeSet* volumeSet = new VolumeSet(fileName);
    VolumeSeries* volumeSeries = new VolumeSeries(Modality::MODALITY_FLOW.getName(),
                                                  Modality::MODALITY_FLOW);
    volumeSet->addSeries(volumeSeries);
    Volume* volumeHandle = new Volume(volume, 0.0f);

    volumeSeries->addVolumeHandle(volumeHandle); */

    VolumeList* collection = new VolumeList();
    Volume* volumeHandle = new Volume(volume, tgt::vec3(1.0f), tgt::vec3(0.0f));//FIXME: spacing?
    oldVolumePosition(volumeHandle);
    volumeHandle->setModality(Modality::MODALITY_FLOW);
    collection->add(volumeHandle);

    // TODO: origin does not save series and timestamp anymore
    //volumeHandle->setOrigin(fileName, Modality::MODALITY_FLOW.getName(), 0.0f);
    volumeHandle->setOrigin(fileName);

    return collection;
}

VolumeFlow3D* FlowReader::readConvert(const tgt::ivec3& dimensions, const BYTE orientation,
                                      std::fstream& ifs)
{
    tgt::ivec3 permutation = Flow3D::getAxisPermutation(orientation);
    switch (orientation) {
        case Flow3D::XZY:
            LINFO("voxel order: XZY (code = " << static_cast<int>(orientation) << ")");
            LINFO("changing voxel order from XZY to XYZ...");
            break;
        case Flow3D::YXZ:
            LINFO("voxel order: YXZ (code = " << static_cast<int>(orientation) << ")");
            LINFO("changing voxel order from YXZ to XYZ...");
            break;
        case Flow3D::YZX:
            LINFO("voxel order: YZX (code = " << static_cast<int>(orientation) << ")");
            LINFO("changing voxel order from YZX to XYZ...");
            break;
        case Flow3D::ZXY:
            LINFO("voxel order: ZXY (code = " << static_cast<int>(orientation) << ")");
            LINFO("changing voxel order from ZXY to XYZ...");
            break;
        case Flow3D::ZYX:
            LINFO("voxel order: ZYX (code = " << static_cast<int>(orientation) << ")");
            LINFO("changing voxel order from ZYX to XYZ...");
            break;
        case Flow3D::XYZ:
            LINFO("voxel order: XYZ (code = " << static_cast<int>(orientation) << ")");
            LINFO("no conversion of voxel order required.");
            break;
        default:
            LINFO("voxel order: unknown (code = " << static_cast<int>(orientation) << ")");
            LINFO("further processing not possible! canceling");
            return 0;
    }

    const size_t numVoxels = dimensions.x * dimensions.y * dimensions.z;
    const size_t byteSize = numVoxels * sizeof(tgt::vec3);

    tgt::vec3* voxels = 0;
    try {
        voxels = new tgt::vec3[numVoxels];
    } catch (std::bad_alloc) {
        throw;
    }

    const tgt::ivec3 xyz(0, 1, 2);  // target voxel order is XYZ
    const int max = tgt::max(dimensions);

    // no remainder possible because getNumBytes is a multiple of max
    //
    const size_t bufferByteSize = byteSize / static_cast<size_t>(max);
    const size_t numBufferElements = bufferByteSize / sizeof(tgt::vec3);
    tgt::vec3* buffer = new tgt::vec3[numBufferElements];
    memset(buffer, 0, bufferByteSize);

    ProgressBar* progress = getProgressBar();

    tgt::ivec3 pos(0, 0, 0);
    float maxValue = std::numeric_limits<float>::min();
    float minValue = std::numeric_limits<float>::max();
    float maxMagnitude = 0.0f;

    for (size_t i = 0; i < size_t(max); ++i) {
        ifs.read(reinterpret_cast<char*>(buffer), bufferByteSize);
        if (progress)
            progress->setProgress(static_cast<float>(i) / static_cast<float>(max));

        for (size_t j = 0; j < numBufferElements; ++j) {
            // get the number of the voxel in current data orientation
            //
            size_t voxelNumber = j + (i * numBufferElements);

            if (orientation != Flow3D::XYZ) {
                // convert the number into the position within the voxel
                // ((x, y, z)-coordinates in volume's bounding box)
                //
                pos = Flow3D::voxelNumberToPos(voxelNumber, permutation, dimensions);

                // re-calculate the number of the voxel for the desired
                // data orientation
                //
                size_t newNumber = Flow3D::posToVoxelNumber(pos, xyz, dimensions);
                voxels[newNumber] = buffer[j];
            } else
                voxels[voxelNumber] = buffer[j];

            // calculate max and min of current voxel and dataset's min and max.
            //
            float voxelMax = tgt::max(buffer[j]);
            float voxelMin = tgt::min(buffer[j]);
            float magnitude = tgt::length(buffer[j]);

            if (voxelMax > maxValue)
                maxValue = voxelMax;
            if (voxelMin < minValue)
                minValue = voxelMin;
            if (magnitude > maxMagnitude)
                maxMagnitude = magnitude;
        }
    }

    delete [] buffer;

    LINFO("min. among all components of all voxels: " << minValue);
    LINFO("max. among all components of all voxels: " << maxValue);
    LINFO("max. magnitude among all  voxels: " << maxMagnitude);

    return new VolumeFlow3D(voxels, dimensions, minValue, maxValue, maxMagnitude);
}

void FlowReader::reverseSlices(VolumeFlow3D* const volume, const BYTE sliceOrder) const {
    switch (sliceOrder) {
        case 'x':
        case 'X':
            LINFO("reversing slice order in x-direction...\n");
            reverseXSliceOrder(volume);
            break;
        case 'y':
        case 'Y':
            LINFO("reversing slice order in y-direction...\n");
            reverseYSliceOrder(volume);
            break;
        case 'z':
        case 'Z':
            LINFO("reversing slice order in z-direction...\n");
            reverseZSliceOrder(volume);
            break;
        default:
            break;
    }
}

VolumeReader* FlowReader::create(ProgressBar* progress) const {
    return new FlowReader(progress);
}

}   // namespace

#endif  // VRN_MODULE_FLOWREEN
