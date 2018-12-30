/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "volumecrop.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorsubset.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

namespace voreen {

using tgt::vec3;
using tgt::ivec3;
using tgt::mat4;

const std::string VolumeCrop::loggerCat_("voreen.base.VolumeCrop");

VolumeCrop::VolumeCrop()
    : CachingVolumeProcessor()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , boundingBoxPort_(Port::INPORT, "boundingbox.input", "Boundingbox Input (optional)")
    , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false)
    , clipRegion_("clipRegion", "Clip Region", tgt::IntBounds(tgt::ivec3::zero, tgt::ivec3::one), tgt::ivec3::zero, tgt::ivec3::one, tgt::ivec3::one)
    , continuousCropping_("continuousCropping", "Continuous Cropping", false)
    , button_("button", "Crop")
    , croppedDimensions_("croppedDimensions", "Cropped Dimensions",
        tgt::ivec3(0), tgt::ivec3(0), tgt::ivec3(100000))
    , clipRegionMatchesCroppedRegion_("clipRegionMatchesCroppedRegion", "Is crop current", false)
    , croppedSize_("croppedSize", "Cropped Size (MB)", 0, 0, 100000)
    , isCropped_(false)
{
    addPort(inport_);
    addPort(boundingBoxPort_);
    addPort(outport_);

    button_.onChange(MemberFunctionCallback<VolumeCrop>(this, &VolumeCrop::crop));
    clipRegion_.onChange(MemberFunctionCallback<VolumeCrop>(this, &VolumeCrop::updateInfoPropertys));
    inport_.onChange(MemberFunctionCallback<VolumeCrop>(this, &VolumeCrop::inportChanged));

    addProperty(clipRegion_);
    addProperty(continuousCropping_);
    addProperty(button_);
    addProperty(clipRegionMatchesCroppedRegion_);
    clipRegionMatchesCroppedRegion_.setReadOnlyFlag(true);

    addProperty(croppedDimensions_);
    addProperty(croppedSize_);

    croppedDimensions_.setReadOnlyFlag(true);
    croppedSize_.setReadOnlyFlag(true);

}

Processor* VolumeCrop::create() const {
    return new VolumeCrop();
}

void VolumeCrop::process() {
    tgtAssert(inport_.getData(), "no input volume");

    // adapt clipping plane properties on volume change

    tgt::svec3 inputDimensions = inport_.getData()->getDimensions();

    // compute size of cropped volume in each rendering pass
    tgt::IntBounds bounds = clipRegion_.get();
    tgt::ivec3 llf = bounds.getLLF();
    tgt::ivec3 urb = bounds.getURB();
    tgt::ivec3 dimensions = urb - llf + tgt::ivec3::one;
    croppedDimensions_.setMaxValue(tgt::ivec3(inputDimensions));
    croppedDimensions_.set(dimensions);

    size_t bpp = static_cast<size_t>(inport_.getData()->getBytesPerVoxel());
    size_t inputRamSize = hmul(inputDimensions) * bpp / (1024 * 1024);
    size_t croppedRamSize = hmul(dimensions) * bpp / (1024 * 1024);
    croppedSize_.setMaxValue(static_cast<int>(inputRamSize));
    croppedSize_.set(static_cast<int>(croppedRamSize));

    // actually crop volume only, if continuous cropping is enabled (potentially slow)
    if (continuousCropping_.get())
        crop();
}

bool VolumeCrop::isReady() const {
    return isInitialized() && inport_.isReady() && outport_.isReady();
}


void VolumeCrop::crop() {
    if (!inport_.hasData())
        return;

    tgt::IntBounds bounds = clipRegion_.get();
    tgt::ivec3 llf = bounds.getLLF();
    tgt::ivec3 urb = bounds.getURB();
    tgt::ivec3 start = llf;
    tgt::ivec3 dimensions = urb - llf + tgt::ivec3::one;

    croppedRegion_ = bounds;
    isCropped_ = true;

    Volume* outputVolume = VolumeOperatorSubset::APPLY_OP(inport_.getData(), start, dimensions, this);

    // assign computed volume to outport
    outport_.setData(outputVolume);
    updateInfoPropertys();
}

void VolumeCrop::adjustPropertiesToInput() {
    if(inport_.hasData()) {
        // adapt clipping plane properties to volume dimensions
        tgt::ivec3 dim = tgt::ivec3(inport_.getData()->getDimensions());

        clipRegion_.setMaxValue(dim-tgt::ivec3::one);
        clipRegion_.setMinRange(tgt::ivec3::zero);
        clipRegion_.setMaxRange(dim);

        if(boundingBoxPort_.hasData()) {
            std::unique_ptr<Geometry> geom = boundingBoxPort_.getData()->clone();
            geom->transform(inport_.getData()->getWorldToVoxelMatrix());

            // Geometry is now in voxel space of the input volume
            tgt::Bounds bounds = geom->getBoundingBox();
            clipRegion_.set(tgt::IntBounds(
                        tgt::max(tgt::ivec3::zero, tgt::ivec3(tgt::round(bounds.getLLF()))),
                        tgt::min(dim-tgt::ivec3::one, tgt::ivec3(tgt::round(bounds.getURB())))
                        ));
            clipRegion_.setReadOnlyFlag(true);
        } else {
            clipRegion_.setReadOnlyFlag(false);
        }
    }
}

void VolumeCrop::updateInfoPropertys()
{
    clipRegionMatchesCroppedRegion_.set(clipRegion_.get()==croppedRegion_ && isCropped_);
}

void VolumeCrop::inportChanged()
{
    isCropped_ = false;
    outport_.setData(0);
}

}   // namespace
