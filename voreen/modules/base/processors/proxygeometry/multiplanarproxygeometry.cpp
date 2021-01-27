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

#include "multiplanarproxygeometry.h"

namespace voreen {

using tgt::vec3;

MultiPlanarProxyGeometry::MultiPlanarProxyGeometry()
  : Processor()
  , inport_(Port::INPORT, "volume.in", "Volume Input")
  , outport_(Port::OUTPORT, "proxygeometry.out", "Proxy Geometry Output")
  , slicePosX_("slicePosX", "Slice Pos (X)", 0.5f, 0.f, 1.f)
  , slicePosY_("slicePosY", "Slice Pos (Y)", 0.5f, 0.f, 1.f)
  , slicePosZ_("slicePosZ", "Slice Pos (Z)", 0.5f, 0.f, 1.f)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(slicePosX_);
    addProperty(slicePosY_);
    addProperty(slicePosZ_);
}

Processor* MultiPlanarProxyGeometry::create() const {
    return new MultiPlanarProxyGeometry();
}

void MultiPlanarProxyGeometry::process() {
    // input volume
    const VolumeBase* inputVolume = inport_.getData();
    tgtAssert(inputVolume, "No input volume");
    tgt::vec3 llf = inputVolume->getLLF();
    tgt::vec3 urb = inputVolume->getURB();
    tgt::vec3 texLlf = tgt::vec3(0.f);
    tgt::vec3 texUrb = tgt::vec3(1.f);

    //
    // x-face
    //
    const float slicePosX = llf.x + slicePosX_.get()*(urb.x - llf.x);
    const float sliceTexX = slicePosX_.get();
    VertexGeometry xLL(tgt::vec3(slicePosX, llf.y, llf.z), tgt::vec3(sliceTexX, texLlf.y, texLlf.z));
    VertexGeometry xLR(tgt::vec3(slicePosX, urb.y, llf.z), tgt::vec3(sliceTexX, texUrb.y, texLlf.z));
    VertexGeometry xUR(tgt::vec3(slicePosX, urb.y, urb.z), tgt::vec3(sliceTexX, texUrb.y, texUrb.z));
    VertexGeometry xUL(tgt::vec3(slicePosX, llf.y, urb.z), tgt::vec3(sliceTexX, texLlf.y, texUrb.z));

    FaceGeometry xFace;
    xFace.addVertex(xLL);
    xFace.addVertex(xLR);
    xFace.addVertex(xUR);
    xFace.addVertex(xUL);

    // We need to double each of the three faces with antiparallel normal vectors,
    // in order to make sure that we have front and back faces. This is necessary,
    // because the MeshEntryExitPoints processor derives the exit points from back faces.

    // these offsets are added to the back faces' slice and tex coords for preventing
    // the ray direction (exitPoint - entryPoint) becoming the null vector.
    const float texCoordOffset = 0.001f;
    const float sliceCoordOffset = 0.00001f;

    const float slicePosXBack = slicePosX-sliceCoordOffset;
    const float sliceTexXBack = sliceTexX-texCoordOffset;
    VertexGeometry xLLBack(tgt::vec3(slicePosXBack, llf.y, llf.z), tgt::vec3(sliceTexXBack, texLlf.y, texLlf.z));
    VertexGeometry xLRBack(tgt::vec3(slicePosXBack, urb.y, llf.z), tgt::vec3(sliceTexXBack, texUrb.y, texLlf.z));
    VertexGeometry xURBack(tgt::vec3(slicePosXBack, urb.y, urb.z), tgt::vec3(sliceTexXBack, texUrb.y, texUrb.z));
    VertexGeometry xULBack(tgt::vec3(slicePosXBack, llf.y, urb.z), tgt::vec3(sliceTexXBack, texLlf.y, texUrb.z));
    FaceGeometry xFaceBack;   //< reverse order of vertices for back face
    xFaceBack.addVertex(xLLBack);
    xFaceBack.addVertex(xULBack);
    xFaceBack.addVertex(xURBack);
    xFaceBack.addVertex(xLRBack);


    //
    // y-face
    //
    const float slicePosY = llf.y + slicePosY_.get()*(urb.y - llf.y);
    const float sliceTexY = slicePosY_.get();
    VertexGeometry yLL(tgt::vec3(llf.x, slicePosY, llf.z), tgt::vec3(texLlf.x, sliceTexY, texLlf.z));
    VertexGeometry yLR(tgt::vec3(urb.x, slicePosY, llf.z), tgt::vec3(texUrb.x, sliceTexY, texLlf.z));
    VertexGeometry yUR(tgt::vec3(urb.x, slicePosY, urb.z), tgt::vec3(texUrb.x, sliceTexY, texUrb.z));
    VertexGeometry yUL(tgt::vec3(llf.x, slicePosY, urb.z), tgt::vec3(texLlf.x, sliceTexY, texUrb.z));
    FaceGeometry yFace;
    yFace.addVertex(yLL);
    yFace.addVertex(yUL);
    yFace.addVertex(yUR);
    yFace.addVertex(yLR);

    // y back face (see above)
    const float slicePosYBack = slicePosY-sliceCoordOffset;
    const float sliceTexYBack = sliceTexY-texCoordOffset;
    VertexGeometry yLLBack(tgt::vec3(llf.x, slicePosYBack, llf.z), tgt::vec3(texLlf.x, sliceTexYBack, texLlf.z));
    VertexGeometry yLRBack(tgt::vec3(urb.x, slicePosYBack, llf.z), tgt::vec3(texUrb.x, sliceTexYBack, texLlf.z));
    VertexGeometry yURBack(tgt::vec3(urb.x, slicePosYBack, urb.z), tgt::vec3(texUrb.x, sliceTexYBack, texUrb.z));
    VertexGeometry yULBack(tgt::vec3(llf.x, slicePosYBack, urb.z), tgt::vec3(texLlf.x, sliceTexYBack, texUrb.z));
    FaceGeometry yFaceBack;
    yFaceBack.addVertex(yLLBack);
    yFaceBack.addVertex(yLRBack);
    yFaceBack.addVertex(yURBack);
    yFaceBack.addVertex(yULBack);


    //
    // z-face
    //
    const float slicePosZ = llf.z + slicePosZ_.get()*(urb.z - llf.z);
    const float sliceTexZ = slicePosZ_.get();
    VertexGeometry zLL(tgt::vec3(llf.x, llf.y, slicePosZ), tgt::vec3(texLlf.x, texLlf.y, sliceTexZ));
    VertexGeometry zLR(tgt::vec3(urb.x, llf.y, slicePosZ), tgt::vec3(texUrb.x, texLlf.y, sliceTexZ));
    VertexGeometry zUR(tgt::vec3(urb.x, urb.y, slicePosZ), tgt::vec3(texUrb.x, texUrb.y, sliceTexZ));
    VertexGeometry zUL(tgt::vec3(llf.x, urb.y, slicePosZ), tgt::vec3(texLlf.x, texUrb.y, sliceTexZ));
    FaceGeometry zFace;
    zFace.addVertex(zLL);
    zFace.addVertex(zLR);
    zFace.addVertex(zUR);
    zFace.addVertex(zUL);

    // z back face (see above)
    const float slicePosZBack = slicePosZ-sliceCoordOffset;
    const float sliceTexZBack = sliceTexZ-texCoordOffset;
    VertexGeometry zLLBack(tgt::vec3(llf.x, llf.y, slicePosZBack), tgt::vec3(texLlf.x, texLlf.y, sliceTexZBack));
    VertexGeometry zLRBack(tgt::vec3(urb.x, llf.y, slicePosZBack), tgt::vec3(texUrb.x, texLlf.y, sliceTexZBack));
    VertexGeometry zURBack(tgt::vec3(urb.x, urb.y, slicePosZBack), tgt::vec3(texUrb.x, texUrb.y, sliceTexZBack));
    VertexGeometry zULBack(tgt::vec3(llf.x, urb.y, slicePosZBack), tgt::vec3(texLlf.x, texUrb.y, sliceTexZBack));
    FaceGeometry zFaceBack;
    zFaceBack.addVertex(zLLBack);
    zFaceBack.addVertex(zULBack);
    zFaceBack.addVertex(zURBack);
    zFaceBack.addVertex(zLRBack);

    // construct output mesh from faces
    MeshGeometry* geometry = new MeshGeometry();
    geometry->addFace(xFace);
    geometry->addFace(xFaceBack);
    geometry->addFace(yFace);
    geometry->addFace(yFaceBack);
    geometry->addFace(zFace);
    geometry->addFace(zFaceBack);

    // assign vertex tex coords as vertex colors so the mesh can be rendered directly (debugging, ...)
    for (size_t faceID=0; faceID<geometry->getFaceCount(); faceID++) {
        for (size_t vertexID=0; vertexID<geometry->getFace(faceID).getVertexCount(); vertexID++) {
            tgt::vec3 texCoords = geometry->getFace(faceID).getVertex(vertexID).getTexCoords();
            geometry->getFace(faceID).getVertex(vertexID).setColor(texCoords);
        }
    }

    geometry->setTransformationMatrix(inport_.getData()->getPhysicalToWorldMatrix());

    outport_.setData(geometry);
}

} // namespace voreen
