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

#include "cubeproxygeometry.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"
namespace voreen {

const std::string CubeProxyGeometry::loggerCat_("voreen.base.CubeProxyGeometry");

CubeProxyGeometry::CubeProxyGeometry()
    : Processor()
    , inport_(Port::INPORT, "volumehandle.volumehandle", "Volume Input")
    , outport_(Port::OUTPORT, "proxygeometry.geometry", "Proxy Geometry Output")
    , enableClipping_("useClipping", "Enable Clipping", true)
    , clipRegion_("clipRegion", "Clip Region", tgt::IntBounds(tgt::ivec3::zero, tgt::ivec3(10000)), tgt::ivec3::zero, tgt::ivec3(10000))
    , resetClipPlanes_("resetClipPlanes", "Reset Planes")
{

    addPort(inport_);
    addPort(outport_);

    enableClipping_.onChange(MemberFunctionCallback<CubeProxyGeometry>(this, &CubeProxyGeometry::adjustClipPropertiesVisibility));
    resetClipPlanes_.onChange(MemberFunctionCallback<CubeProxyGeometry>(this, &CubeProxyGeometry::resetClipPlanes));

    addProperty(enableClipping_);
    addProperty(clipRegion_);
    addProperty(resetClipPlanes_);

    clipRegion_.setGroupID("clipping");
    resetClipPlanes_.setGroupID("clipping");
    setPropertyGroupGuiName("clipping", "Clipping Planes");
    adjustClipPropertiesVisibility();

    oldVolumeDimensions_ = tgt::ivec3(0,0,0);
    inport_.onNewData(MemberFunctionCallback<CubeProxyGeometry>(this, &CubeProxyGeometry::adjustClipPropertiesRanges));
    inport_.onNewData(MemberFunctionCallback<CubeProxyGeometry>(this, &CubeProxyGeometry::adjustClippingToVolumeROI));
}

CubeProxyGeometry::~CubeProxyGeometry() {
}

Processor* CubeProxyGeometry::create() const {
    return new CubeProxyGeometry();
}

void CubeProxyGeometry::process() {
    tgtAssert(inport_.getData(), "no input volume");

    // get voxel space dimensions
    const VolumeBase* inputVolume = inport_.getData();
    tgt::vec3 voxelLLF = tgt::vec3(0.f);
    tgt::vec3 voxelURB = tgt::vec3(inputVolume->getDimensions() - tgt::svec3::one);

    // vertex coords of bounding box in voxel space without clipping (add a border of 0.5 voxels since a voxel is sampled at the center)
    tgt::vec3 coordLLF = voxelLLF - tgt::vec3(0.5f);
    tgt::vec3 coordURB = voxelURB + tgt::vec3(0.5f);

    // take into account the clipping
    if (enableClipping_.get()) {
        coordLLF = tgt::max(coordLLF, tgt::vec3(clipRegion_.get().getLLF()) - tgt::vec3(0.5f));
        coordURB = tgt::min(coordURB, tgt::vec3(clipRegion_.get().getURB()) + tgt::vec3(0.5f));
    }

    // compute the texture coordinates for the voxels
    tgt::vec3 texLLF = inputVolume->getVoxelToTextureMatrix() * coordLLF;
    tgt::vec3 texURB = inputVolume->getVoxelToTextureMatrix() * coordURB;

    // create output mesh
    TriangleMeshGeometryColorNormal* mesh = TriangleMeshGeometryColorNormal::createCube(texLLF, texURB, texLLF, texURB, 1.0f);
    mesh->transform(inputVolume->getTextureToWorldMatrix());

    outport_.setData(mesh, true);
}

void CubeProxyGeometry::resetClipPlanes() {
    clipRegion_.set(tgt::IntBounds(clipRegion_.getMinValue(), clipRegion_.getMaxValue()));
}

void CubeProxyGeometry::adjustClipPropertiesRanges() {
    tgtAssert(inport_.getData() && inport_.getData(), "No input volume");

    tgt::ivec3 oldRange = clipRegion_.getMaxValue();
    tgt::ivec3 newRange = tgt::ivec3(inport_.getData()->getDimensions())-tgt::ivec3::one;

    if (oldRange != newRange){
        clipRegion_.setMaxValue(newRange);
        clipRegion_.set(tgt::IntBounds(tgt::ivec3::zero, newRange));
    }
}

void CubeProxyGeometry::adjustClipPropertiesVisibility() {
    bool clipEnabled = enableClipping_.get();
    setPropertyGroupVisible("clipping", clipEnabled);
}

void CubeProxyGeometry::adjustClippingToVolumeROI() {
    // adjust clipping sliders to volume ROI, if set
    if (inport_.getData()->hasMetaData("RoiPixelOffset") && inport_.getData()->hasMetaData("RoiPixelLength")) {
        const tgt::ivec3 volumeDim = static_cast<tgt::ivec3>(inport_.getData()->getDimensions());
        tgt::ivec3 roiLlf = inport_.getData()->getMetaDataValue<IVec3MetaData>("RoiPixelOffset", tgt::ivec3(0));
        tgt::ivec3 roiLength = inport_.getData()->getMetaDataValue<IVec3MetaData>("RoiPixelLength", volumeDim);
        if (tgt::hor(tgt::lessThan(roiLlf, tgt::ivec3::zero)) || tgt::hor(tgt::greaterThanEqual(roiLlf, volumeDim))) {
            LWARNING("Invalid ROI offset: " << roiLlf);
        }
        else if (tgt::hor(tgt::lessThanEqual(roiLength, tgt::ivec3::zero))) {
            LWARNING("Invalid ROI length: " << roiLength);
        }
        else {
            tgt::ivec3 roiUrb = tgt::min(roiLlf + roiLength - 1, volumeDim - tgt::ivec3::one);

            LINFO("Applying volume ROI: llf=" << roiLlf << ", urb=" << roiUrb);

            clipRegion_.set(tgt::IntBounds(roiLlf, roiUrb));
        }
    }
}

} // namespace
