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

#include "volumelandmarkregistration.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"
#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"

// Eigen library
#include "Eigen/Dense"

namespace voreen {

using tgt::vec3;
using tgt::vec4;
using tgt::mat4;

const std::string VolumeLandmarkRegistration::loggerCat_("voreen.base.VolumeLandmarkRegistration");

VolumeLandmarkRegistration::VolumeLandmarkRegistration()
    : VolumeProcessor()
    , inportVolume_(Port::INPORT, "volume.in", "Volume Input")
    , inportReference_(Port::INPORT, "volume.reference", "Reference Volume Input")
    , inportLandmarks_(Port::INPORT, "landmarks.in", "Landmarks Input")
    , inportReferenceLandmarks_(Port::INPORT, "landmarks.reference", "Reference Landmarks Input")
    , outport_(Port::OUTPORT, "volume.registered", "Registered Volume Output", false)
    , enableProcessing_("enableProcessing", "Enable")
    , registrationMode_("registrationMode", "Registration Mode")
    , landmarkSpecification_("landmarkSpecification", "Landmark Specification")
    , p1_("p1", "Point 1", vec3(0.f), vec3(0.f), vec3(4096.f))
    , ref_p1_("ref_p1", "Point 1 (Reference)", vec3(0.f), vec3(0.f), vec3(4096.f))
    , p2_("p2", "Point 2", vec3(1.f, 0.f, 0.f), vec3(0.f), vec3(4096.f))
    , ref_p2_("ref_p2", "Point 2 (Reference)", vec3(1.f, 0.f, 0.f), vec3(0.f), vec3(4096.f))
    , p3_("p3", "Point 3", vec3(0.f, 1.f, 0.f), vec3(0.f), vec3(4096.f))
    , ref_p3_("ref_p3", "Point 3 (Reference)", vec3(0.f, 1.f, 0.f), vec3(0.f), vec3(4096.f))
    , p4_("p4", "Point 4", vec3(0.f, 0.f, 1.f), vec3(0.f), vec3(4096.f))
    , ref_p4_("ref_p4", "Point 4 (Reference)", vec3(0.f, 0.f, 1.f), vec3(0.f), vec3(4096.f))
    , registrationMatrix_("registrationMatrix", "Registration Matrix",
                          tgt::mat4::zero, tgt::mat4(-10000.f), tgt::mat4(10000.f), VALID)
{
    addPort(inportVolume_);
    addPort(inportReference_);
    addPort(inportLandmarks_);
    addPort(inportReferenceLandmarks_);
    addPort(outport_);

    registrationMode_.addOption("rigidBody", "Rigid Body");
    registrationMode_.addOption("affine",    "Affine");
    registrationMode_.selectByKey("affine");
    registrationMode_.onChange(MemberFunctionCallback<VolumeLandmarkRegistration>(this, &VolumeLandmarkRegistration::updateGuiState));

    landmarkSpecification_.addOption("properties", "Property Input");
    landmarkSpecification_.addOption("geometry",   "Geometry Ports");
    landmarkSpecification_.onChange(MemberFunctionCallback<VolumeLandmarkRegistration>(this, &VolumeLandmarkRegistration::updateGuiState));

    //registrationMatrix_.setReadOnlyFlag(true);

    addProperty(enableProcessing_);
    addProperty(registrationMode_);
    addProperty(landmarkSpecification_);

    addProperty(p1_);
    addProperty(ref_p1_);
    addProperty(p2_);
    addProperty(ref_p2_);
    addProperty(p3_);
    addProperty(ref_p3_);
    addProperty(p4_);
    addProperty(ref_p4_);
    addProperty(registrationMatrix_);
}

Processor* VolumeLandmarkRegistration::create() const {
    return new VolumeLandmarkRegistration();
}

bool VolumeLandmarkRegistration::isReady() const {
    bool mainPortsReady = inportVolume_.isReady() && inportReference_.isReady() && outport_.isConnected();
    bool geometryPortsReady = !landmarkSpecification_.isSelected("geometry") || (inportLandmarks_.isReady() && inportReferenceLandmarks_.isReady());
    return mainPortsReady && geometryPortsReady;
}

void VolumeLandmarkRegistration::process() {

    if (!enableProcessing_.get()) {
        outport_.setData(const_cast<VolumeBase*>(inportVolume_.getData()), false);
        registrationMatrix_.set(tgt::mat4::identity);
    }
    else {

        // perform registration
        tgt::mat4 registrationMatrix = tgt::mat4::identity;

        if (registrationMode_.isSelected("rigidBody")) {
            if (landmarkSpecification_.isSelected("properties")) {
                // landmarks from properties
                registrationMatrix = computeRigidBodyTransformation(inportVolume_.getData(), inportReference_.getData(),
                    p1_.get(), ref_p1_.get(), p2_.get(), ref_p2_.get(), p3_.get(), ref_p3_.get());
            }
            else if (landmarkSpecification_.isSelected("geometry")) {
                // landmarks from input geometries
                const PointListGeometryVec3* landmarksIn = dynamic_cast<const PointListGeometryVec3*>(inportLandmarks_.getData());
                const PointListGeometryVec3* landmarksReference = dynamic_cast<const PointListGeometryVec3*>(inportReferenceLandmarks_.getData());
                if (!landmarksIn && !landmarksReference) {
                    LWARNING("PointListGeometryVec3 expected.");
                }
                else if (landmarksIn->getNumPoints() < 3 || landmarksReference->getNumPoints() < 3) {
                    LWARNING("Three pairs of landmarks required for rigid-body registration.");
                }
                else {
                    tgt::vec3 p1 = landmarksIn->getData().at(0);
                    tgt::vec3 p2 = landmarksIn->getData().at(1);
                    tgt::vec3 p3 = landmarksIn->getData().at(2);
                    tgt::vec3 p1Ref = landmarksReference->getData().at(0);
                    tgt::vec3 p2Ref = landmarksReference->getData().at(1);
                    tgt::vec3 p3Ref = landmarksReference->getData().at(2);
                    registrationMatrix = computeRigidBodyTransformation(inportVolume_.getData(), inportReference_.getData(),
                        p1, p1Ref, p2, p2Ref, p3, p3Ref);
                }
            }
            else {
                LWARNING("Unknown landmark specification mode: " << landmarkSpecification_.get());
            }
        }

        else if (registrationMode_.isSelected("affine")) {
            if (landmarkSpecification_.isSelected("properties")) {
                // landmarks from properties
                std::vector<tgt::vec3> landmarksIn;
                std::vector<tgt::vec3> landmarksReference;
                landmarksIn.push_back(p1_.get());
                landmarksIn.push_back(p2_.get());
                landmarksIn.push_back(p3_.get());
                landmarksIn.push_back(p4_.get());
                landmarksReference.push_back(ref_p1_.get());
                landmarksReference.push_back(ref_p2_.get());
                landmarksReference.push_back(ref_p3_.get());
                landmarksReference.push_back(ref_p4_.get());
                registrationMatrix = computeAffineTransformation(inportVolume_.getData(), inportReference_.getData(),
                    landmarksIn, landmarksReference);
            }
            else if (landmarkSpecification_.isSelected("geometry")) {
                // landmarks from input geometries
                const PointListGeometryVec3* landmarksIn = dynamic_cast<const PointListGeometryVec3*>(inportLandmarks_.getData());
                const PointListGeometryVec3* landmarksReference = dynamic_cast<const PointListGeometryVec3*>(inportReferenceLandmarks_.getData());
                if (landmarksIn && landmarksReference) {
                    registrationMatrix = computeAffineTransformation(inportVolume_.getData(), inportReference_.getData(),
                        landmarksIn->getData(), landmarksReference->getData());
                }
                else {
                    LWARNING("PointListGeometryVec3 expected.");
                }
            }

        }
        else {
            LWARNING("Unknown registration mode: " << registrationMode_.get());
        }

        // assign computed registration matrix to output volume and respective property
        registrationMatrix_.set(registrationMatrix);

        VolumeBase* outputVolume =
            new VolumeDecoratorReplaceTransformation(inportVolume_.getData(), registrationMatrix);

        outport_.setData(outputVolume);
    }
}

void VolumeLandmarkRegistration::initialize() {
    VolumeProcessor::initialize();

    updateGuiState();
}

// private methods
//

tgt::mat4 VolumeLandmarkRegistration::computeRigidBodyTransformation(const VolumeBase* inputHandle, const VolumeBase* refHandle,
        tgt::vec3 p0, tgt::vec3 p0Ref, tgt::vec3 p1, tgt::vec3 p1Ref, tgt::vec3 p2, tgt::vec3 p2Ref) {

    tgtAssert(inputHandle && inputHandle->getRepresentation<VolumeRAM>(), "no volume");
    tgtAssert(refHandle && refHandle->getRepresentation<VolumeRAM>(), "no volume");

    //
    // transform reference points from voxel coordinates to volume coordinates
    //
    p0 = transformFromVoxelToVolumeCoords(p0, inputHandle);
    p0Ref = transformFromVoxelToVolumeCoords(p0Ref, refHandle);

    p1 = transformFromVoxelToVolumeCoords(p1, inputHandle);
    p1Ref = transformFromVoxelToVolumeCoords(p1Ref, refHandle);

    p2 = transformFromVoxelToVolumeCoords(p2, inputHandle);
    p2Ref = transformFromVoxelToVolumeCoords(p2Ref, refHandle);

    //
    // calculate uniform scaling matrix
    //
    float scale01 = length(p0Ref - p1Ref) / length(p0 - p1);
    float scale02 = length(p0Ref - p2Ref) / length(p0 - p2);
    float scale12 = length(p1Ref - p2Ref) / length(p1 - p2);
    LINFO("Scaling factors: " << scale01 << ", " << scale02 << ", " << scale12);

    float avgScale = (scale01 + scale02 + scale12) / 3.f; // values may differ due to inaccuracies of reference points
    if (!tgt::isFiniteNumber(avgScale) || avgScale == 0.f) {
        LWARNING("Rigid-body registration with uniform scaling impossible: invalid scaling");
        return tgt::mat4::identity;
    }
    mat4 inputScale = mat4::createScale(vec3(avgScale));

    //
    // calculate corresponding tripods in input and ref coordinate systems from reference points
    //
    vec3 p01 = p0 - p1;
    vec3 p01Ref = p0Ref - p1Ref;
    vec3 a = normalize(p01);
    vec3 aRef = normalize(p01Ref);

    vec3 p02 = p0 - p2;
    vec3 p02Ref = p0Ref - p2Ref;
    vec3 b = normalize(cross(p01, p02));
    vec3 bRef = normalize(cross(p01Ref, p02Ref));

    vec3 c = cross(b, a);
    vec3 cRef = cross(bRef, aRef);

    LINFO("Input tripod: " << a << ", " << b << ", " << c);
    LINFO("Reference tripod: " << aRef << ", " << bRef << ", " << cRef);
    if (!tgt::isFiniteNumber(length(a)) || length(a) == 0.f ||
        !tgt::isFiniteNumber(length(b)) || length(b) == 0.f ||
        !tgt::isFiniteNumber(length(c)) || length(c) == 0.f  ) {
        LWARNING("Rigid-body registration with uniform scaling impossible: invalid input tripod (points collinear?)");
        return tgt::mat4::identity;
    }
    if (!tgt::isFiniteNumber(length(aRef)) || length(aRef) == 0.f ||
        !tgt::isFiniteNumber(length(bRef)) || length(bRef) == 0.f ||
        !tgt::isFiniteNumber(length(cRef)) || length(cRef) == 0.f  ) {
        LWARNING("Rigid-body registration with uniform scaling impossible: invalid reference tripod (points collinear?)");
        return tgt::mat4::identity;
    }

    // rotation matrices mapping the tripod in source coordinates to the tripod in ref coordinates
    mat4 inputRot = mat4( vec4(a, 0.f),
                          vec4(b, 0.f),
                          vec4(c, 0.f),
                          vec4(vec3(0.f), 1.f)
                         );
    mat4 refRot = tgt::transpose( mat4( vec4(aRef, 0.f),
                                        vec4(bRef, 0.f),
                                        vec4(cRef, 0.f),
                                        vec4(vec3(0.f), 1.f)
                                       )
                                 );

    // translation matrices
    mat4 inputTranslate = mat4::createTranslation(-p0);
    mat4 refTranslate =  mat4::createTranslation(p0Ref);

    /* compose final registration matrix:
       1. translate p0 to origin
       2. rotate source tripod to be axis-aligned
       3. adjust (uniform!) scaling
       4. rotate axis-aligned tripod to fit ref tripod
       5. translate to p0ref
    */
    mat4 M = refTranslate * refRot * inputScale * inputRot * inputTranslate;

    LINFO("Registration matrix M mapping from input to reference volume coordinates: " << M);
    vec3 p0Trans = (M*vec4(p0, 1.f)).xyz();
    vec3 p1Trans = (M*vec4(p1, 1.f)).xyz();
    vec3 p2Trans = (M*vec4(p2, 1.f)).xyz();
    LINFO("Length of deviation vectors: ");
    LINFO("|M*P0-P0ref|: " << tgt::length(p0Trans - p0Ref));
    LINFO("|M*P1-P1ref|: " << tgt::length(p1Trans - p1Ref));
    LINFO("|M*P2-P2ref|: " << tgt::length(p2Trans - p2Ref));

    return (refHandle->getPhysicalToWorldMatrix() * M);
}

tgt::mat4 VolumeLandmarkRegistration::computeAffineTransformation(const VolumeBase* inputHandle, const VolumeBase* refHandle,
        std::vector<vec3> inputPoints, std::vector<vec3> refPoints) {

    tgtAssert(inputHandle, "no volume");
    tgtAssert(refHandle, "no volume");

    // check reference points
    int numPoints = static_cast<int>(std::min(inputPoints.size(), refPoints.size()));
    if (inputPoints.size() != refPoints.size()) {
        LWARNING("Unequal number of landmarks (input: " << inputPoints.size() << ", reference: " << refPoints.size() << "). "
            << "Using the first " << numPoints << " pairs.");
    }
    if (numPoints < 4) {
        LWARNING("Affine registration requires at least four pairs of landmarks (passed: " << numPoints << ")");
        return mat4::identity;
    }

    // transform reference points from voxel coordinates to volume coordinates
    for (int i=0; i<numPoints; i++) {
        inputPoints.at(i) = transformFromVoxelToVolumeCoords(inputPoints.at(i), inputHandle);
        refPoints.at(i) = transformFromVoxelToVolumeCoords(refPoints.at(i), refHandle);
    }

    // calculate least-squares solution
    Eigen::VectorXf y(numPoints*3);
    Eigen::MatrixXf X = Eigen::MatrixXf::Zero(numPoints*3,12);
    //result vector
    Eigen::VectorXf c(12);

    // Create matrix X and the corresponding QR-matrix from input reference points
    for (int i=0; i<numPoints; ++i) {
        vec3 inputPoint = inputPoints.at(i);
        // x-component
        X(3*i,0) = inputPoint.x;
        X(3*i,1) = inputPoint.y;
        X(3*i,2) = inputPoint.z;
        X(3*i,3) = 1;
        // y-component
        X(3*i+1,4) = inputPoint.x;
        X(3*i+1,5) = inputPoint.y;
        X(3*i+1,6) = inputPoint.z;
        X(3*i+1,7) = 1;
        // z-component
        X(3*i+2,8) = inputPoint.x;
        X(3*i+2,9) = inputPoint.y;
        X(3*i+2,10) = inputPoint.z;
        X(3*i+2,11) = 1;
    }

    Eigen::JacobiSVD<Eigen::MatrixXf> SVD = X.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);

    //test rank
    if (SVD.matrixU().cols() != 12) {
        LWARNING("The matrix has not full rank: reference points collinear?");
        return mat4::identity;
    }

    // create y-vector from dest reference points
    for (int i=0; i<numPoints; ++i) {
        vec3 destPoint = refPoints.at(i);
        y[3*i] = destPoint.x;
        y[3*i+1] = destPoint.y;
        y[3*i+2] = destPoint.z;
    }

    // least-squares solution
    c = SVD.solve(y);

    // copy solution to mat4
    mat4 M = mat4::identity;
    M[0][0] = c[0];
    M[0][1] = c[1];
    M[0][2] = c[2];
    M[0][3] = c[3];
    M[1][0] = c[4];
    M[1][1] = c[5];
    M[1][2] = c[6];
    M[1][3] = c[7];
    M[2][0] = c[8];
    M[2][1] = c[9];
    M[2][2] = c[10];
    M[2][3] = c[11];

    LINFO("Registration matrix M mapping from input to reference volume coordinates: " << M);
    LINFO("Length of deviation vectors: ");
    for (int i=0; i<numPoints; ++i) {
        vec3 pInput = inputPoints.at(i);
        vec3 pRef = refPoints.at(i);
        vec3 pTrans = (M*vec4(pInput, 1.f)).xyz();
        float dist = tgt::length(pTrans - pRef);
        if (dist < 0.05f)
            LINFO("|M*P" << i << "-P" << i << "ref|: " << dist);
        else
            LWARNING("|M*P" << i << "-P" << i << "ref|: " << dist);
    }

    return (refHandle->getPhysicalToWorldMatrix() * M);
}

vec3 VolumeLandmarkRegistration::transformFromVoxelToVolumeCoords(vec3 point, const VolumeBase* vol) {
    return vol->getVoxelToPhysicalMatrix() * point;
}

void VolumeLandmarkRegistration::updateGuiState() {
    bool affine = registrationMode_.isSelected("affine");
    bool propertyInput = landmarkSpecification_.isSelected("properties");

    p1_.setVisibleFlag(propertyInput);
    ref_p1_.setVisibleFlag(propertyInput);
    p2_.setVisibleFlag(propertyInput);
    ref_p2_.setVisibleFlag(propertyInput);
    p3_.setVisibleFlag(propertyInput);
    ref_p3_.setVisibleFlag(propertyInput);
    p4_.setVisibleFlag(affine && propertyInput);
    ref_p4_.setVisibleFlag(affine && propertyInput);
}


}   // namespace
