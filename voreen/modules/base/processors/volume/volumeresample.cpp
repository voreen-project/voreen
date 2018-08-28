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

#include "volumeresample.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorresample.h"

#include <climits>

namespace voreen {

const std::string VolumeResample::loggerCat_("voreen.base.VolumeResample");

VolumeResample::VolumeResample()
    : CachingVolumeProcessor()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false)
    , enableProcessingProp_("enableProcessing", "Enable")
    , allowUpsampling_("allowUpsampling", "Allow Upsampling", false)
    , keepVoxelRatio_("keepVoxelSizeRatio", "Keep Voxel Ratio", true)
    , resampleDimensionX_("resampleDimensionX", "Dimension X", 1, 1, 4096)
    , resampleDimensionY_("resampleDimensionY", "Dimension Y", 1, 1, 4096)
    , resampleDimensionZ_("resampleDimensionZ", "Dimension Z", 1, 1, 4096)
    , filteringMode_("filteringMode", "Filtering")
    , outputSizeVoxels_("outputSizeVoxels", "Output size (voxels)", 0, 0, INT_MAX)
    , outputSizeMB_("outputSizeMB", "Output size (MB)", 0, 0, INT_MAX)
    , forceUpdate_(true)
    , blockSignals_(false)
{
    addPort(inport_);
    addPort(outport_);

    filteringMode_.addOption("nearest", "Nearest");
    filteringMode_.addOption("linear",  "Linear ");
    filteringMode_.addOption("cubic",   "Cubic");
    filteringMode_.set("linear");

    enableProcessingProp_.onChange(
        MemberFunctionCallback<VolumeResample>(this, &VolumeResample::forceUpdate));
    resampleDimensionX_.onChange(
        MemberFunctionCallback<VolumeResample>(this, &VolumeResample::forceUpdate));
    resampleDimensionY_.onChange(
        MemberFunctionCallback<VolumeResample>(this, &VolumeResample::forceUpdate));
    resampleDimensionZ_.onChange(
        MemberFunctionCallback<VolumeResample>(this, &VolumeResample::forceUpdate));
    keepVoxelRatio_.onChange(
        MemberFunctionCallback<VolumeResample>(this, &VolumeResample::forceUpdate));
    filteringMode_.onChange(
        MemberFunctionCallback<VolumeResample>(this, &VolumeResample::forceUpdate));

    resampleDimensionX_.onChange(
        MemberFunctionOneParameterCallback<VolumeResample, int>(this, &VolumeResample::dimensionsChanged, 0));
    resampleDimensionY_.onChange(
        MemberFunctionOneParameterCallback<VolumeResample, int>(this, &VolumeResample::dimensionsChanged, 1));
    resampleDimensionZ_.onChange(
        MemberFunctionOneParameterCallback<VolumeResample, int>(this, &VolumeResample::dimensionsChanged, 2));

    allowUpsampling_.onChange(
        MemberFunctionCallback<VolumeResample>(this, &VolumeResample::allowUpsamplingChanged));
    keepVoxelRatio_.onChange(
        MemberFunctionCallback<VolumeResample>(this, &VolumeResample::keepVoxelRatioChanged));

    resampleDimensionX_.setTracking(false);
    resampleDimensionY_.setTracking(false);
    resampleDimensionZ_.setTracking(false);

    addProperty(enableProcessingProp_);
    addProperty(allowUpsampling_);
    addProperty(keepVoxelRatio_);
    addProperty(resampleDimensionX_);
    addProperty(resampleDimensionY_);
    addProperty(resampleDimensionZ_);
    addProperty(filteringMode_);

    addProperty(outputSizeVoxels_);
    addProperty(outputSizeMB_);
    outputSizeVoxels_.setReadOnlyFlag(true);
    outputSizeMB_.setReadOnlyFlag(true);
}

Processor* VolumeResample::create() const {
    return new VolumeResample();
}

void VolumeResample::process() {

    if (inport_.hasChanged())
        adjustDimensionProperties();

    // update output size properties
    const VolumeRAM* inputVolume = inport_.getData()->getRepresentation<VolumeRAM>();
    tgtAssert(inputVolume, "No input volume");
    tgt::svec3 inputDim = inputVolume->getDimensions();
    tgt::svec3 outputDim(resampleDimensionX_.get(), resampleDimensionY_.get(), resampleDimensionZ_.get());
    outputSizeVoxels_.setMaxValue(static_cast<int>(tgt::hmul(tgt::max(inputDim, outputDim))));
    outputSizeVoxels_.set(static_cast<int>(tgt::hmul(outputDim)));
    float mbPerVoxel = inputVolume->getBytesPerVoxel() / (1024.f * 1024.f);
    outputSizeMB_.setMaxValue(tgt::iround(tgt::hmul(tgt::max(inputDim, outputDim)) * mbPerVoxel));
    outputSizeMB_.set(tgt::iround(tgt::hmul(outputDim) * mbPerVoxel));

    // process input volume
    if (!enableProcessingProp_.get()) {
        outport_.setData(const_cast<VolumeBase*>(inport_.getData()), false);
    }
    else if (forceUpdate_ || inport_.hasChanged()) {
        resampleVolume();
    }
}

// private methods
//

void VolumeResample::resampleVolume() {
    tgtAssert(inport_.hasData(), "Inport has not data");
    forceUpdate_ = false;

    if (inport_.getData()->getRepresentation<VolumeRAM>()) {

        VolumeRAM::Filter filter;
        if (filteringMode_.isSelected("nearest"))
            filter = VolumeRAM::NEAREST;
        else if (filteringMode_.isSelected("linear"))
            filter = VolumeRAM::LINEAR;
        else if (filteringMode_.isSelected("cubic"))
            filter = VolumeRAM::CUBIC;
        else {
            LERROR("Unknown filter mode: " << filteringMode_.get());
            return;
        }

        tgt::ivec3 dimensions(resampleDimensionX_.get(), resampleDimensionY_.get(), resampleDimensionZ_.get());
        try {
            Volume* v = VolumeOperatorResample::APPLY_OP(inport_.getData(), dimensions, filter, this);
            outport_.setData(v);
        }
        catch (const std::bad_alloc&) {
            LERROR("resampleVolume(): bad allocation");
            outport_.setData(0);
        }
    }
    else {
        outport_.setData(0);
    }
}

void VolumeResample::adjustDimensionProperties() {
    if (!inport_.hasData() || !inport_.getData()->getRepresentation<VolumeRAM>())
        return;

    tgt::ivec3 volDim = inport_.getData()->getRepresentation<VolumeRAM>()->getDimensions();

    if (!allowUpsampling_.get()) {
        resampleDimensionX_.setMaxValue(volDim.x);
        resampleDimensionY_.setMaxValue(volDim.y);
        resampleDimensionZ_.setMaxValue(volDim.z);
    }
    else {
        resampleDimensionX_.setMaxValue(4096);
        resampleDimensionY_.setMaxValue(4096);
        resampleDimensionZ_.setMaxValue(4096);
    }
    resampleDimensionX_.updateWidgets();
    resampleDimensionY_.updateWidgets();
    resampleDimensionZ_.updateWidgets();

    if (resampleDimensionX_.get() == 1)
        resampleDimensionX_.set(volDim.x);
    if (resampleDimensionY_.get() == 1)
        resampleDimensionY_.set(volDim.y);
    if (resampleDimensionZ_.get() == 1)
        resampleDimensionZ_.set(volDim.z);

    resampleDimensionX_.set(tgt::clamp(resampleDimensionX_.get(), 1, resampleDimensionX_.getMaxValue()));
    resampleDimensionY_.set(tgt::clamp(resampleDimensionY_.get(), 1, resampleDimensionY_.getMaxValue()));
    resampleDimensionZ_.set(tgt::clamp(resampleDimensionZ_.get(), 1, resampleDimensionZ_.getMaxValue()));
}

void VolumeResample::forceUpdate() {
    if (blockSignals_)
        return;

    forceUpdate_ = true;
}

void VolumeResample::dimensionsChanged(int dim) {

    if (blockSignals_ || !keepVoxelRatio_.get() || !inport_.hasData() || !inport_.getData()->getRepresentation<VolumeRAM>())
        return;

    tgt::ivec3 volDim = inport_.getData()->getRepresentation<VolumeRAM>()->getDimensions();

    blockSignals_ = true;
    if (dim == 0) {
        int xDim = resampleDimensionX_.get();
        float yRatio = static_cast<float>(volDim.y) / volDim.x;
        float zRatio = static_cast<float>(volDim.z) / volDim.x;
        int yDim = tgt::iround(xDim * yRatio);
        int zDim = tgt::iround(xDim * zRatio);
        resampleDimensionY_.set(tgt::clamp(yDim, 1, resampleDimensionY_.getMaxValue()));
        resampleDimensionZ_.set(tgt::clamp(zDim, 1, resampleDimensionZ_.getMaxValue()));
    }
    else if (dim == 1) {
        int yDim = resampleDimensionY_.get();
        float xRatio = static_cast<float>(volDim.x) / volDim.y;
        float zRatio = static_cast<float>(volDim.z) / volDim.y;
        int xDim = tgt::iround(yDim * xRatio);
        int zDim = tgt::iround(yDim * zRatio);
        resampleDimensionX_.set(tgt::clamp(xDim, 1, resampleDimensionX_.getMaxValue()));
        resampleDimensionZ_.set(tgt::clamp(zDim, 1, resampleDimensionZ_.getMaxValue()));
    }
    else if (dim == 2) {
        int zDim = resampleDimensionZ_.get();
        float xRatio = static_cast<float>(volDim.x) / volDim.z;
        float yRatio = static_cast<float>(volDim.y) / volDim.z;
        int xDim = tgt::iround(zDim * xRatio);
        int yDim = tgt::iround(zDim * yRatio);
        resampleDimensionX_.set(tgt::clamp(xDim, 1, resampleDimensionX_.getMaxValue()));
        resampleDimensionY_.set(tgt::clamp(yDim, 1, resampleDimensionY_.getMaxValue()));
    }
    blockSignals_ = false;
}

void VolumeResample::allowUpsamplingChanged() {
    if (blockSignals_)
        return;

    adjustDimensionProperties();
}

void VolumeResample::keepVoxelRatioChanged() {
    if (blockSignals_)
        return;

    if (keepVoxelRatio_.get())
        dimensionsChanged(0);
}



}   // namespace
