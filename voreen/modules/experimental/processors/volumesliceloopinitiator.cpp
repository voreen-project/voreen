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

#include "volumesliceloopinitiator.h"

namespace voreen {

using tgt::Texture;
using tgt::ivec2;

const std::string VolumeSliceLoopInitiator::loggerCat_("voreen.experimental.VolumeSliceLoopInitiator");

VolumeSliceLoopInitiator::VolumeSliceLoopInitiator()
    : RenderProcessor(),
      sliceAlignment_("sliceAlignmentProp", "Slice Alignment"),
      minSliceIndex_("minSliceIndex", "Starting Slice Number ", 0, 0, 10000),
      maxSliceIndex_("maxSliceIndex", "Last Slice Number ", 0, 0, 10000),
      inport_(Port::INPORT, "volume.in"),
      outport_(Port::OUTPORT, "image.out"),
      loopInport_(Port::INPORT, "loop.in")
{
    loopInport_.setLoopPort(true);

    sliceAlignment_.addOption("xy-plane", "XY-Plane (axial)", XY_PLANE);
    sliceAlignment_.addOption("xz-plane", "XZ-Plane (coronal)", XZ_PLANE);
    sliceAlignment_.addOption("yz-plane", "YZ-Plane (sagittal)", YZ_PLANE);
    sliceAlignment_.onChange( CallMemberAction<VolumeSliceLoopInitiator>(this, &VolumeSliceLoopInitiator::updateSliceProperties) );
    addProperty(sliceAlignment_);

    addProperty(minSliceIndex_);
    addProperty(maxSliceIndex_);

    addPort(inport_);
    addPort(outport_);
    addPort(loopInport_);
}

Processor* VolumeSliceLoopInitiator::create() const {
    return new VolumeSliceLoopInitiator();
}

bool VolumeSliceLoopInitiator::isReady() const {
    return (inport_.isReady() && outport_.isReady());
}

void VolumeSliceLoopInitiator::process() {
    LGL_ERROR;

    const VolumeRAM* vol = inport_.getData()->getRepresentation<VolumeRAM>();
    tgt::ivec3 dims = vol->getDimensions();

    if (inport_.hasChanged())
        updateSliceProperties();  // validate the currently set values and adjust them if necessary

    switch(sliceAlignment_.getValue()) {
        case YZ_PLANE: {
                           int x = minSliceIndex_.get() + loopInport_.getLoopIteration();
                           if(x >= dims.x-1) {
                                LERROR("Out of bounds!");
                                return;
                           }

                           outport_.resize(dims.yz());
                           Texture* tex = outport_.getColorTexture();
                           tex->destroy();
                           tex->alloc();

                           for(int y=0; y<dims.y; y++) {
                               for(int z=0; z<dims.z; z++) {
                                   uint16_t value = tgt::iround(vol->getVoxelNormalized(x, y, z) * 65535.0f);
                                   tex->texel< tgt::Vector4<uint16_t> >(y, z) = tgt::Vector4<uint16_t>(value, value, value, 65535);
                               }
                           }

                           tex->uploadTexture();
                       }
                       break;
        case XZ_PLANE: {
                           int y = minSliceIndex_.get() + loopInport_.getLoopIteration();
                           if(y >= dims.y-1) {
                                LERROR("Out of bounds!");
                                return;
                           }

                           outport_.resize(ivec2(dims.x, dims.z));
                           Texture* tex = outport_.getColorTexture();
                           tex->destroy();
                           tex->alloc();

                           for(int x=0; x<dims.x; x++) {
                               for(int z=0; z<dims.z; z++) {
                                   uint16_t value = tgt::iround(vol->getVoxelNormalized(x, y, z) * 65535.0f);
                                   tex->texel< tgt::Vector4<uint16_t> >(x, z) = tgt::Vector4<uint16_t>(value, value, value, 65535);
                               }
                           }

                           tex->uploadTexture();
                       }
                       break;
        case XY_PLANE: {
                           int z = minSliceIndex_.get() + loopInport_.getLoopIteration();
                           if(z >= dims.z-1) {
                                LERROR("Out of bounds!");
                                return;
                           }

                           outport_.resize(dims.xy());
                           Texture* tex = outport_.getColorTexture();
                           tex->destroy();
                           tex->alloc();

                           for(int y=0; y<dims.y; y++) {
                               for(int x=0; x<dims.x; x++) {
                                   uint16_t value = tgt::iround(vol->getVoxelNormalized(x, y, z) * 65535.0f);
                                   tex->texel< tgt::Vector4<uint16_t> >(x, y) = tgt::Vector4<uint16_t>(value, value, value, 65535);
                               }
                           }

                           tex->uploadTexture();
                       }
                       break;
        default: tgtAssert(false, "should not get here!");
    }

    outport_.validateResult();
    LGL_ERROR;
}

void VolumeSliceLoopInitiator::invalidate(int inv) {
    RenderProcessor::invalidate(inv);

    if (!inport_.hasData())
        loopInport_.setNumLoopIterations(0);
    else if ((maxSliceIndex_.get() - minSliceIndex_.get() + 1) != loopInport_.getNumLoopIterations())
        loopInport_.setNumLoopIterations(maxSliceIndex_.get() - minSliceIndex_.get() + 1);
}

void VolumeSliceLoopInitiator::updateSliceProperties() {
    tgt::ivec3 volumeDim(0);
    if (inport_.getData() && inport_.getData()->getRepresentation<VolumeRAM>())
        volumeDim = inport_.getData()->getRepresentation<VolumeRAM>()->getDimensions();

    tgtAssert(sliceAlignment_.getValue() >= 0 && sliceAlignment_.getValue() <= 2, "Invalid alignment value");
    int numSlices = volumeDim[sliceAlignment_.getValue()];
    if (numSlices == 0)
        return;

    minSliceIndex_.setMaxValue(numSlices-1);
    maxSliceIndex_.setMaxValue(numSlices-1);
    //minSliceIndex_.set(0);
    if(maxSliceIndex_.get() > maxSliceIndex_.getMaxValue())
        maxSliceIndex_.set(numSlices-1);
}

} // voreen namespace
