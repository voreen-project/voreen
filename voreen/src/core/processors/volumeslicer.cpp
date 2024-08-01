/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "voreen/core/processors/volumeslicer.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"

namespace voreen {

const std::string VolumeSlicer::loggerCat_("voreen.VolumeSlicer");

VolumeSlicer::VolumeSlicer()
    : VolumeRenderer()
    , samplingRate_("samplingRate", "Sampling rate", 4.f, 0.01f, 20.f)
    , transferFunc_("transferFunction", "Transfer function", Processor::INVALID_RESULT)//, TransFuncProperty::Editors::INTENSITY)
    , camera_("camera", "Camera", tgt::Camera(tgt::vec3(0.f, 0.f, 3.5f), tgt::vec3(0.f, 0.f, 0.f), tgt::vec3(0.f, 1.f, 0.f)))
    , cubeVertices_(0)
    , nSeq_(0)
    , v1_(0)
    , v2_(0)
    , geometryInport_(Port::INPORT, "proxgeometry.geometry", "Proxy Geometry Input")
    , volumeInport_(Port::INPORT, "volumehandle.volumehandle", "Volume Input", false, Processor::INVALID_PROGRAM)
    , maxLength_(0.f)
    , frontIdx_(0)
    , backIdx_(0)
    , sliceDistance_(0.1f)
{
    volumeInport_.addCondition(new PortConditionVolumeTypeGL());
    addPort(volumeInport_);
    addPort(geometryInport_);

    addProperty(samplingRate_);
    addProperty(transferFunc_);
    addProperty(camera_);
}

VolumeSlicer::~VolumeSlicer() {
    delete[] cubeVertices_;
    delete[] nSeq_;
    delete[] v1_;
    delete[] v2_;
}

void VolumeSlicer::initialize() {
    VolumeRenderer::initialize();

    cubeVertices_ = new tgt::vec3[8];

    nSeq_ = new int[64]; // permutation indices
    nSeq_[0]=0; nSeq_[1]=1; nSeq_[2]=2; nSeq_[3]=3; nSeq_[4]=4; nSeq_[5]=5; nSeq_[6]=6; nSeq_[7]=7;
    nSeq_[8]=1; nSeq_[9]=4; nSeq_[10]=5; nSeq_[11]=0; nSeq_[12]=3; nSeq_[13]=7; nSeq_[14]=2; nSeq_[15]=6;
    nSeq_[16]=2; nSeq_[17]=0; nSeq_[18]=5; nSeq_[19]=6; nSeq_[20]=3; nSeq_[21]=1; nSeq_[22]=7; nSeq_[23]=4;
    nSeq_[24]=3; nSeq_[25]=0; nSeq_[26]=6; nSeq_[27]=4; nSeq_[28]=1; nSeq_[29]=2; nSeq_[30]=7; nSeq_[31]=5;
    nSeq_[32]=4; nSeq_[33]=3; nSeq_[34]=7; nSeq_[35]=1; nSeq_[36]=0; nSeq_[37]=6; nSeq_[38]=5; nSeq_[39]=2;
    nSeq_[40]=5; nSeq_[41]=1; nSeq_[42]=7; nSeq_[43]=2; nSeq_[44]=0; nSeq_[45]=4; nSeq_[46]=6; nSeq_[47]=3;
    nSeq_[48]=6; nSeq_[49]=7; nSeq_[50]=3; nSeq_[51]=2; nSeq_[52]=5; nSeq_[53]=4; nSeq_[54]=0; nSeq_[55]=1;
    nSeq_[56]=7; nSeq_[57]=5; nSeq_[58]=4; nSeq_[59]=6; nSeq_[60]=2; nSeq_[61]=1; nSeq_[62]=3; nSeq_[63]=0;

    v1_ = new int[24]; // edge start indices
    v1_[0]=0; v1_[1]=1; v1_[2]=4; v1_[3]=4; // clip vertex 1
    v1_[4]=1; v1_[5]=0; v1_[6]=1; v1_[7]=4; // clip vertex 2
    v1_[8]=0; v1_[9]=2; v1_[10]=5; v1_[11]=5; // clip vertex 3
    v1_[12]=2; v1_[13]=0; v1_[14]=2; v1_[15]=5; // clip vertex 4
    v1_[16]=0; v1_[17]=3; v1_[18]=6; v1_[19]=6; // clip vertex 5
    v1_[20]=3; v1_[21]=0; v1_[22]=3; v1_[23]=6; // clip vertex 6

    v2_ = new int[24]; // edge end indices
    v2_[0]=1; v2_[1]=4; v2_[2]=7; v2_[3]=7; // clip vertex 1
    v2_[4]=5; v2_[5]=1; v2_[6]=4; v2_[7]=7; // clip vertex 2
    v2_[8]=2; v2_[9]=5; v2_[10]=7; v2_[11]=7; // clip vertex 3
    v2_[12]=6; v2_[13]=2; v2_[14]=5; v2_[15]=7; // clip vertex 4
    v2_[16]=3; v2_[17]=6; v2_[18]=7; v2_[19]=7; // clip vertex 5
    v2_[20]=4; v2_[21]=3; v2_[22]=6; v2_[23]=7; // clip vertex 6
}

void VolumeSlicer::setupUniforms(tgt::Shader* slicingPrg) {

    if(!volumeInport_.isReady())
        return;

    tgt::vec3 camPos  = camera_.get().getPosition();
    tgt::vec3 camLook = camera_.get().getLook();

    // Handle Transformation matrix of dataset.
    tgt::mat4 invTransMat = volumeInport_.getData()->getWorldToPhysicalMatrix();
    camPos = (invTransMat * tgt::vec4(camPos, 1.0)).xyz();
    camLook = normalize((invTransMat * tgt::vec4(camLook, 0.0)).xyz());

    tgt::vec3 urb;
    tgt::vec3 llf;

    // disable other than cube-formed proxygeometries for now
    if(!geometryInport_.isReady()) {
        MeshGeometry meshGeometry = volumeInport_.getData()->getBoundingBox(false);
        FaceGeometry f1 = meshGeometry[0];
        FaceGeometry f2 = meshGeometry[1];
        urb = meshGeometry.getTransformationMatrix() * f1[0].getCoords();
        llf = meshGeometry.getTransformationMatrix() * f2[0].getCoords();

    } else {
        //tgt::Bounds b = geometryInport_.getData()->getBoundingBox(true);
        tgt::Bounds b = geometryInport_.getData()->getBoundingBox(false);
        tgt::mat4 m = geometryInport_.getData()->getTransformationMatrix();
        urb = m * b.getURB();
        llf = m * b.getLLF();

        urb = (invTransMat * tgt::vec4(urb, 1.0)).xyz();
        llf = (invTransMat * tgt::vec4(llf, 1.0)).xyz();
    }

    cubeVertices_[0] = tgt::vec3(llf.x, urb.y, urb.z);
    cubeVertices_[1] = tgt::vec3(urb.x, urb.y, urb.z);
    cubeVertices_[2] = tgt::vec3(llf.x, urb.y, llf.z);
    cubeVertices_[3] = tgt::vec3(llf.x, llf.y, urb.z);
    cubeVertices_[4] = tgt::vec3(urb.x, llf.y, urb.z);
    cubeVertices_[5] = tgt::vec3(urb.x, urb.y, llf.z);
    cubeVertices_[6] = tgt::vec3(llf.x, llf.y, llf.z);
    cubeVertices_[7] = tgt::vec3(urb.x, llf.y, llf.z);

    // compute front index needed for clipping
    float minCameraDistance = std::numeric_limits<float>::max(); // distance between camera and closest vertex
    float maxCameraDistance = 0.0;

    for (unsigned int i=0;i<8;i++) {
        float cameraDistance = dot(camLook, cubeVertices_[i] - camPos);
        if (cameraDistance < minCameraDistance) {
            minCameraDistance = cameraDistance;
            frontIdx_ = i;
        }
        if (cameraDistance > maxCameraDistance) {
            maxCameraDistance = cameraDistance;
            backIdx_ = i;
        }
    }

    maxLength_ = dot(camLook, cubeVertices_[backIdx_] - cubeVertices_[frontIdx_]);

    // acquire volume size for sampling
    tgt::ivec3 dim = volumeInport_.getData()->getRepresentation<VolumeGL>()->getTexture()->getDimensions();

    // compute the distance between two adjacent slices in texture coordinates
    sliceDistance_ = maxLength_ / (tgt::max(dim) * samplingRate_.get());

    //if (interactionMode())
    //    sliceDistance /= interactionQuality_.get();

    //clipping uniforms
    slicingPrg->setUniform("frontIdx_", frontIdx_);
    slicingPrg->setUniform("vecView_", camLook);
    slicingPrg->setUniform("dPlaneStart_", dot(camLook, cubeVertices_[frontIdx_]));
    slicingPrg->setUniform("dPlaneIncr_", sliceDistance_);
    slicingPrg->setUniform("nSequence_", nSeq_, 64);
    slicingPrg->setUniform("vecVertices_", cubeVertices_, 8);
    slicingPrg->setUniform("v1_", v1_, 24);
    slicingPrg->setUniform("v2_", v2_, 24);
}


} // namespace voreen
