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

#include "voreen/core/datastructures/volume/volumepreview.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/octree/volumeoctreebase.h"

namespace voreen {

const std::string VolumePreview::loggerCat_("voreen.VolumePreview");

VolumePreview::VolumePreview()
    : VolumeDerivedData()
    , height_(0)
{}

VolumePreview::VolumePreview(int height, const std::vector<unsigned char>& data)
    : VolumeDerivedData()
    , height_(height)
    , prevData_(data)
{}

VolumeDerivedData* VolumePreview::create() const {
    return new VolumePreview();
}

VolumeDerivedData* VolumePreview::createFrom(const VolumeBase* handle) const {
    tgtAssert(handle, "no volume");

    int internHeight = 64;

    // gamma correction factor (amplifies low gray values)
    const float GAMMA = 1.8f;

    float xSpacing = handle->getSpacing()[0];
    float ySpacing = handle->getSpacing()[1];
    float xDimension = static_cast<float>(handle->getDimensions()[0]);
    float yDimension = static_cast<float>(handle->getDimensions()[1]);

    // determine offsets and scale factors for non-uniform aspect ratios
    float aspectRatio = (yDimension * ySpacing) / (xDimension * xSpacing);
    float xOffset, yOffset, xScale, yScale;
    if (aspectRatio <= 1.f) {
        xOffset = 0.f;
        yOffset = ((internHeight - 1.f) * (1.f - aspectRatio)) / 2.f;
        xScale = 1.f;
        yScale = 1.f / aspectRatio;
    }
    else {
        xOffset = ((internHeight -1.f) * (1.f - 1.f/aspectRatio)) / 2.f;
        yOffset = 0.f;
        xScale = aspectRatio;
        yScale = 1.f;
    }

    float maxVal, minVal;
    std::vector<float> prevData = std::vector<float>(internHeight * internHeight);

    const VolumeRAM* volumeRam = 0;
    const VolumeDisk* volumeDisk = 0;
    const VolumeOctreeBase* volumeOctree = 0;
    bool volumeRamCreated = false;
    if (handle->hasRepresentation<VolumeRAM>()) {
        volumeRam = handle->getRepresentation<VolumeRAM>();
        tgtAssert(volumeRam, "no volume");
    }
    else if (handle->hasRepresentation<VolumeDisk>()) {
        volumeDisk = handle->getRepresentation<VolumeDisk>();
        tgtAssert(volumeDisk, "no volume");
        volumeRamCreated = true;
    }
    else if (handle->hasRepresentation<VolumeOctreeBase>()) {
        volumeOctree = handle->getRepresentation<VolumeOctreeBase>();
        tgtAssert(volumeOctree, "no volume");
        volumeRamCreated = true;
    }

    tgt::vec3 position;
    int offset = static_cast<int>(handle->getDimensions().z - 1) / 2;

    if (volumeRam)
        position.z = static_cast<float>(offset);
    else if (volumeDisk) {
        try {
            volumeRam = volumeDisk->loadSlices(offset, offset);
        }
        catch (tgt::Exception& e) {
            LERROR(e.what());
        }
        if (!volumeRam) {
            LERROR("VolumeDisk::loadSlices failed to create a RAM volume");
            return 0;
        }
        position.z = 0;
    }
    else if (volumeOctree) {
        try {
            volumeRam = volumeOctree->createSlice(XY_PLANE, offset);
        }
        catch (tgt::Exception& e) {
            LERROR(e.what());
        }
        if (!volumeRam) {
            LERROR("VolumeOctree::loadSlices failed to create a RAM volume");
            return 0;
        }
        position.z = 0;
    }
    else {
        LERROR("Neither VolumeRAM nor VolumeDisk nor VolumeOctree available");
        return 0;
    }

    // generate preview in float buffer
    minVal = volumeRam->elementRange().y;
    maxVal = volumeRam->elementRange().x;
    for (int y = 0; y < internHeight; y++){
        for (int x = 0; x < internHeight; x++){
            position.x = ((x-xOffset) / (internHeight -1)) * xScale * (xDimension-1.f);
            position.y = ((y-yOffset) / (internHeight -1)) * yScale * (yDimension-1.f);
            int previewIndex = y*internHeight + x;
            float val = 0.f;
            if (position.x >= 0 && position.y >= 0 && position.x < xDimension && position.y < yDimension)
                val = volumeRam->getVoxelNormalizedLinear(position);
            prevData[previewIndex] = val;
            minVal = std::min(minVal, val);
            maxVal = std::max(maxVal, val);
        }
    }

    if (volumeRamCreated) {
        delete volumeRam;
        volumeRam = 0;
    }

    float valOffset = minVal;
    float valScale = maxVal - minVal;
    if (valScale > 0.f) {
        for (int y=0; y<internHeight; y++) {
            for (int x=0; x<internHeight; x++) {
                int previewIndex = y*internHeight+ x;
                prevData[previewIndex] = pow((prevData[previewIndex] - valOffset)  / valScale, 1.f / GAMMA);
            }
        }
    }

    std::vector<unsigned char> prevDataUInt8 = std::vector<unsigned char>(internHeight * internHeight);
    for (size_t i = 0; i < prevDataUInt8.size(); i++)
        prevDataUInt8[i] = tgt::iround(prevData[i] * 255.f);

    return new VolumePreview(internHeight, prevDataUInt8);
}

void VolumePreview::serialize(Serializer& s) const  {
    s.serialize("height", height_);
    s.serializeBinaryBlob("prevData", prevData_);
}

void VolumePreview::deserialize(Deserializer& s) {
    s.deserialize("height", height_);
    s.deserializeBinaryBlob("prevData", prevData_);
}

} // namespace voreen

