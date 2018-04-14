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

#include "slicegenerator.h"

#include "tgt/glmath.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

using tgt::ivec2;
using tgt::vec3;
using tgt::mat4;
using tgt::Texture;

namespace voreen {

const std::string SliceGenerator::loggerCat_("voreen.SliceGenerator");

SliceGenerator::SliceGenerator()
    : RenderProcessor()
    , sliceAlignment_("sliceAlignmentProp", "Slice Alignment")
    , sliceIndex_("sliceIndex", "Slice Number ", 0, 0, 10000)
    , applyDatasetTransformationMatrix_("useDatasetTrafoMatrix", "Apply data set trafo matrix", true)
    , mwheelCycleHandler_("mouseWheelHandler", "Slice Cycling", &sliceIndex_)
    , outport_(Port::OUTPORT, "slice")
    , inport_(Port::INPORT, "volume")
    , geomPort_(Port::OUTPORT, "geometry")
{
    addInteractionHandler(mwheelCycleHandler_);

    sliceAlignment_.addOption("xy-plane", "XY-Plane (axial)", XY_PLANE);
    sliceAlignment_.addOption("xz-plane", "XZ-Plane (coronal)", XZ_PLANE);
    sliceAlignment_.addOption("yz-plane", "YZ-Plane (sagittal)", YZ_PLANE);
    sliceAlignment_.onChange( MemberFunctionCallback<SliceGenerator>(this, &SliceGenerator::updateSliceProperties) );
    addProperty(sliceAlignment_);

    addProperty(sliceIndex_);

    addPort(inport_);
    addPort(geomPort_);
    addPort(outport_);
    addProperty(applyDatasetTransformationMatrix_);
}

SliceGenerator::~SliceGenerator() {
}

std::string SliceGenerator::getProcessorInfo() const {
    return "Creates slices directly from volumes without using OpenGl.";
}

bool SliceGenerator::isReady() const {
    return (inport_.isReady() && outport_.isReady());
}

void SliceGenerator::process() {
    LGL_ERROR;

    const VolumeBase* volh = inport_.getData();
    const VolumeRAM* vol = inport_.getData()->getRepresentation<VolumeRAM>();
    tgt::ivec3 dims = vol->getDimensions();

    if (inport_.hasChanged())
        updateSliceProperties();  // validate the currently set values and adjust them if necessary

    FaceGeometry slice;
    vec3 urb = volh->getURB();
    vec3 llf = volh->getLLF();

    switch(sliceAlignment_.getValue()) {
        case YZ_PLANE: {
                           outport_.resize(dims.yz());
                           Texture* tex = outport_.getColorTexture();
                           tex->alloc(true);

                           int x = sliceIndex_.get();
                           for(int y=0; y<dims.y; y++) {
                               for(int z=0; z<dims.z; z++) {
                                   uint16_t value = tgt::iround(vol->getVoxelNormalized(x, y, z) * 65535.0f);
                                   tex->texel< tgt::Vector4<uint16_t> >(y, z) = tgt::Vector4<uint16_t>(value, value, value, 65535);
                               }
                           }

                           tex->uploadTexture();

                           float mix = x / (float)sliceIndex_.getMaxValue();
                           float xcoord = (mix*urb.x) + ((1.0f-mix) * llf.x);

                           slice.addVertex(VertexGeometry(tgt::vec3(xcoord, urb.y, urb.z), tgt::vec3(1.0f, 1.0f, 0.0f), tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f)));
                           slice.addVertex(VertexGeometry(tgt::vec3(xcoord, urb.y, llf.z), tgt::vec3(1.0f, 0.0f, 0.0f), tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f)));
                           slice.addVertex(VertexGeometry(tgt::vec3(xcoord, llf.y, llf.z), tgt::vec3(0.0f, 0.0f, 0.0f), tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f)));
                           slice.addVertex(VertexGeometry(tgt::vec3(xcoord, llf.y, urb.z), tgt::vec3(0.0f, 1.0f, 0.0f), tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f)));
                       }
                       break;
        case XZ_PLANE: {
                           outport_.resize(ivec2(dims.x, dims.z));
                           Texture* tex = outport_.getColorTexture();
                           tex->alloc(true);

                           int y = sliceIndex_.get();
                           for(int x=0; x<dims.x; x++) {
                               for(int z=0; z<dims.z; z++) {
                                   uint16_t value = tgt::iround(vol->getVoxelNormalized(x, y, z) * 65535.0f);
                                   tex->texel< tgt::Vector4<uint16_t> >(x, z) = tgt::Vector4<uint16_t>(value, value, value, 65535);
                               }
                           }

                           tex->uploadTexture();

                           float mix = y / (float)sliceIndex_.getMaxValue();
                           float ycoord = (mix*urb.y) + ((1.0f-mix) * llf.y);

                           slice.addVertex(VertexGeometry(tgt::vec3(urb.x, ycoord, urb.z), tgt::vec3(1.0f, 1.0f, 0.0f), tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f)));
                           slice.addVertex(VertexGeometry(tgt::vec3(urb.x, ycoord, llf.z), tgt::vec3(1.0f, 0.0f, 0.0f), tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f)));
                           slice.addVertex(VertexGeometry(tgt::vec3(llf.x, ycoord, llf.z), tgt::vec3(0.0f, 0.0f, 0.0f), tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f)));
                           slice.addVertex(VertexGeometry(tgt::vec3(llf.x, ycoord, urb.z), tgt::vec3(0.0f, 1.0f, 0.0f), tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f)));
                       }
                       break;
        case XY_PLANE: {
                           outport_.resize(dims.xy());
                           Texture* tex = outport_.getColorTexture();
                           tex->alloc(true);

                           int z = sliceIndex_.get();
                           for(int y=0; y<dims.y; y++) {
                               for(int x=0; x<dims.x; x++) {
                                   uint16_t value = tgt::iround(vol->getVoxelNormalized(x, y, z) * 65535.0f);
                                   tex->texel< tgt::Vector4<uint16_t> >(x, y) = tgt::Vector4<uint16_t>(value, value, value, 65535);
                               }
                           }

                           tex->uploadTexture();

                           float mix = z / (float)sliceIndex_.getMaxValue();
                           float zcoord = ((1.0f-mix)*llf.z) + (mix * urb.z);

                           slice.addVertex(VertexGeometry(tgt::vec3(urb.x, urb.y, zcoord), tgt::vec3(1.0f, 1.0f, 0.0f), tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f)));
                           slice.addVertex(VertexGeometry(tgt::vec3(urb.x, llf.y, zcoord), tgt::vec3(1.0f, 0.0f, 0.0f), tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f)));
                           slice.addVertex(VertexGeometry(tgt::vec3(llf.x, llf.y, zcoord), tgt::vec3(0.0f, 0.0f, 0.0f), tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f)));
                           slice.addVertex(VertexGeometry(tgt::vec3(llf.x, urb.y, zcoord), tgt::vec3(0.0f, 1.0f, 0.0f), tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f)));
                       }
                       break;
        default: tgtAssert(false, "should not get here!");
    }

    outport_.validateResult();
    LGL_ERROR;

    MeshGeometry mesh;
    mesh.addFace(slice);

    geometry_.clear();
    geometry_.addMesh(mesh);
    if (applyDatasetTransformationMatrix_.get())
        geometry_.transform(inport_.getData()->getPhysicalToWorldMatrix());
    geomPort_.setData(&geometry_, false);
}

void SliceGenerator::updateSliceProperties() {
    tgt::ivec3 volumeDim(0);
    if (inport_.getData() && inport_.getData()->getRepresentation<VolumeRAM>())
        volumeDim = inport_.getData()->getRepresentation<VolumeRAM>()->getDimensions();

    tgtAssert(sliceAlignment_.getValue() >= 0 && sliceAlignment_.getValue() <= 2, "Invalid alignment value");
    int numSlices = volumeDim[sliceAlignment_.getValue()];
    if (numSlices == 0)
        return;

    sliceIndex_.setMaxValue(numSlices-1);
    if (sliceIndex_.get() >= static_cast<int>(numSlices))
        sliceIndex_.set(static_cast<int>(numSlices / 2));
}

} // namespace voreen
