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

#include "commands_registration.h"
#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/datastructures/volume/volumecollection.h"
#include "tgt/tgt_math.h"

#include <vector>

#ifdef VRN_MODULE_JAMA
// for least-squares-fitting
#include "jama/jama_qr.h"
#endif

namespace voreen {

using tgt::vec3;
using tgt::vec4;
using tgt::ivec3;
using tgt::mat4;
using tgt::normalize;
using tgt::length;
using tgt::cross;

CommandRegistration::CommandRegistration(const std::string& name, const std::string& shortName, const std::string& info, const std::string& parameterList, const int argumentNum) /*:
    Command(name, shortName, info, parameterList, argumentNum) */
{}

vec3 CommandRegistration::transformFromVoxelToWorldCoords(vec3 point, const VolumeBase* vol) {
    return vol->getVoxelToPhysicalMatrix() * point;
}

//------------------------------------------------------------------------

CommandRegistrationUniformScaling::CommandRegistrationUniformScaling() /*:
    CommandRegistration("--register", "", "Registers two data sets by calculating an affine transformation with uniform scaling. \n \
                              The algorithm assumes that the registration can be done by a rigid-body transformation \
                              in combination with a uniform scaling. Therefore, it cannot handle the general case of an affine transformation (see register-affine).\n\n \
                              Input: The two Voreen dat volumes to register and three pairs (Pi, Pi') of reference points \
                              in *voxel* coordinates, where each Pi' in the source volume is to be mapped to its counterpart Pi in the dest volume. \
                              Each of the points has to be passed as a sequence of three integers or floats.\n\n \
                              Note: Since the registration matrix is intended to be used by Voreen, \
                              both volumes are assumed to be initially axis-aligned and centered around the origin \
                              with a largest side length of 2.0. Furthermore, the volumes' spacings are considered. \
                              The volumes' transformation matrices are ignored.",
                              "<DEST_VOL SRC_VOL P0 P0' P1 P1' P2 P2>", 20) */
{
//    loggerCat_ += "." + name_;
}

bool CommandRegistrationUniformScaling::checkParameters(const std::vector<std::string>& parameters) {
    bool result = (parameters.size() == 20);

/*    for (int i = 2; i < 20; ++i)
        result &= is<float>(parameters[i]); */

    return result;
}

bool CommandRegistrationUniformScaling::execute(const std::vector<std::string>& parameters) {
    //std::string operation = parameters[0];

    // read source, dest volume and reference points from parameters
    std::string destVolFilename = parameters[0];
    std::string srcVolFilename = parameters[1];

    vec3 p0dest; /* = vec3(cast<float>(parameters[2]), cast<float>(parameters[3]), cast<float>(parameters[4])); */
    vec3 p0src; /* = vec3(cast<float>(parameters[5]), cast<float>(parameters[6]), cast<float>(parameters[7])); */

    vec3 p1dest; /* = vec3(cast<float>(parameters[8]), cast<float>(parameters[9]), cast<float>(parameters[10])); */
    vec3 p1src; /* = vec3(cast<float>(parameters[11]), cast<float>(parameters[12]), cast<float>(parameters[13])); */

    vec3 p2dest; /* = vec3(cast<float>(parameters[14]), cast<float>(parameters[15]), cast<float>(parameters[16])); */
    vec3 p2src; /* = vec3(cast<float>(parameters[17]), cast<float>(parameters[18]), cast<float>(parameters[19])); */

    // load volumes
    VolumeSerializerPopulator volLoadPop;
    const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();
    const VolumeBase* destVol = serializer->read(destVolFilename)->first();
    VolumeBase* srcVol = serializer->read(srcVolFilename)->first();

    // transform reference points from voxel coordinates to world coordinates
    p0dest = transformFromVoxelToWorldCoords(p0dest, destVol);
    p0src = transformFromVoxelToWorldCoords(p0src, srcVol);

    p1dest = transformFromVoxelToWorldCoords(p1dest, destVol);
    p1src = transformFromVoxelToWorldCoords(p1src, srcVol);

    p2dest = transformFromVoxelToWorldCoords(p2dest, destVol);
    p2src = transformFromVoxelToWorldCoords(p2src, srcVol);

    // calculate scaling factor (note: the scaling is assumed to be isotrope)
    float scale01 = length(p0dest - p1dest) / length(p0src - p1src);
    float scale02 = length(p0dest - p2dest) / length(p0src - p2src);
    float scale12 = length(p1dest - p2dest) / length(p1src - p2src);
    float avgScale = (scale01 + scale02 + scale12) / 3.f; // values may differ due to inaccuracies of reference points

    // uniform scaling matrix
    mat4 srcScale = mat4::createScale(vec3(avgScale, avgScale, avgScale));

    // calculate corresponding tripods in source and dest coordinate systems from reference points
    vec3 p01dest = p0dest - p1dest;
    vec3 p01src = p0src - p1src;

    vec3 p02dest = p0dest - p2dest;
    vec3 p02src = p0src - p2src;

    //vec3 p12dest = p1dest - p2dest; //set but not used
    //vec3 p12src = p1src - p2src;

    vec3 aDest = normalize(p01dest);
    vec3 aSrc = normalize(p01src);

    vec3 bDest = normalize(cross(p01dest, p02dest));
    vec3 bSrc = normalize(cross(p01src, p02src));

    vec3 cDest = cross(bDest, aDest);
    vec3 cSrc = cross(bSrc, aSrc);

    // rotation matrices mapping the tripod in source coordinates to the tripod in dest coordinates
    mat4 srcRot = mat4( vec4(aSrc.x, aSrc.y, aSrc.z, 0.f),
                        vec4(bSrc.x, bSrc.y, bSrc.z, 0.f),
                        vec4(cSrc.x, cSrc.y, cSrc.z, 0.f),
                        vec4(0.f, 0.f, 0.f, 1.f) );
    mat4 destRot = tgt::transpose( mat4( vec4(aDest.x, aDest.y, aDest.z, 0.f),
                                         vec4(bDest.x, bDest.y, bDest.z, 0.f),
                                         vec4(cDest.x, cDest.y, cDest.z, 0.f),
                                         vec4(0.f, 0.f, 0.f, 1.f)) );

    // translation matrices
    mat4 srcTranslate = mat4::createTranslation(-p0src);
    mat4 destTranslate = mat4::createTranslation(p0dest);

    /* compose final registration matrix:
       1. translate p0source to origin
       2. rotate source tripod to be axis-aligned
       3. adjust (uniform!) scaling
       4. rotate axis-aligned tripod to fit dest tripod
       5. translate to p0dest
    */
    mat4 M = destTranslate * destRot * srcScale * srcRot * srcTranslate;

    std::cout << "\nRegistration matrix M mapping from source to dest volume coordinates: \n";
    std::cout << M << "\n";

    vec3 p0trans = (M*vec4(p0src.x, p0src.y, p0src.z, 1.f)).xyz();
    vec3 p1trans = (M*vec4(p1src.x, p1src.y, p1src.z, 1.f)).xyz();
    vec3 p2trans = (M*vec4(p2src.x, p2src.y, p2src.z, 1.f)).xyz();

    std::cout << "Deviation vectors (should be zero): \n";
    std::cout << "M*P0'-P0: " << (p0trans - p0dest) << "\n";
    std::cout << "M*P1'-P1: " << (p1trans - p1dest) << "\n";
    std::cout << "M*P2'-P2: " << (p2trans - p2dest) << "\n\n";

    std::string input = "";
    while (input == "") {
        std::cout << "Write registration matrix to source volume '" << srcVolFilename << "' (y or n)? ";
        std::cin >> input;
        if (input != "y" && input != "n")
            input = "";
        std::cout << "\n";
    }

    if (input == "y") {
        Volume* vh = static_cast<Volume*>(srcVol);
        vh->setPhysicalToWorldMatrix(M);
        serializer->write(srcVolFilename, vh);
    }

    return true;

}

//------------------------------------------------------------------------

CommandRegistrationAffine::CommandRegistrationAffine() /*:
    CommandRegistration("--register-affine", "", "Registers two data sets by calculating a general affine transformation. \n \
                                                 Input: The two Voreen dat volumes to register and at least four pairs (Pi, Pi') of reference points \
                                                 in *voxel* coordinates, where each Pi' in the source volume is to be mapped to its counterpart Pi in the dest volume. \
                                                 Each of the points has to be passed as a sequence of three integers or floats. The affine transformation matrix mapping \
                                                 from source volume coordinates to dest volume coordinates is calculated by a least-squares regression.\n\n \
                                                 Note: Since the registration matrix is intended to be used by Voreen, \
                                                 both volumes are assumed to be initially axis-aligned and centered around the origin \
                                                 with a largest side length of 2.0. Furthermore, the volumes' spacings are considered. \
                                                 The volumes' transformation matrices are ignored.") */
{
//    parameterList_ = "DEST_VOL SRC_VOL [Pi Pi']\n";
//    loggerCat_ += "." + name_;
}

bool CommandRegistrationAffine::checkParameters(const std::vector<std::string>& parameters) {
    // only accept 2 strings + a list of pairs of 3D vectors. At least four pairs of reference points are needed.
    return false;

    bool result = (parameters.size() >= 26) && ((parameters.size()-2) % 6 == 0);

    /*for (size_t i = 2; i < parameters.size(); ++i)
        result &= is<float>(parameters[i]); */

    return result;
}

#ifdef VRN_MODULE_JAMA
bool CommandRegistrationAffine::execute(const std::vector<std::string>& parameters) {
    // This doesn't work with the new commandline parser concept, because it requires doesn't support
    // a fixed number of parameters (ab)

    return false;

    // read source, dest volume and reference points from parameters
    std::string destVolFilename = parameters[0];
    std::string srcVolFilename = parameters[1];

    // read referemce points from parameters
    int numReferencePoints = (parameters.size()-2) / 6;
    std::vector<vec3> refPointsSrc;
    std::vector<vec3> refPointsDest;
/*    for (int i=0; i<numReferencePoints; ++i) {
        refPointsDest.push_back(vec3(cast<float>(parameters[2+i*6]), cast<float>(parameters[2+i*6+1]), cast<float>(parameters[2+i*6+2])));
        refPointsSrc.push_back(vec3(cast<float>(parameters[2+i*6+3]), cast<float>(parameters[2+i*6+4]), cast<float>(parameters[2+i*6+5])));
    } */

    // load volumes
    VolumeSerializerPopulator volLoadPop;
    const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();
    const VolumeBase* destVol = serializer->read(destVolFilename)->first();
    Volume* srcVol = dynamic_cast<Volume*>(serializer->read(srcVolFilename)->first());

    // transform reference points from voxel coordinates to world coordinates
    for (int i=0; i<numReferencePoints; ++i) {
        refPointsDest.at(i) = transformFromVoxelToWorldCoords(refPointsDest.at(i), destVol);
        refPointsSrc.at(i) = transformFromVoxelToWorldCoords(refPointsSrc.at(i), srcVol);
    }

    // calculate least-squares solution
    TNT::Array1D<float> y(numReferencePoints*3);
    TNT::Array2D<float> X(numReferencePoints*3, 12, 0.f);
    // result vector
    TNT::Array1D<float> c(12);

    // Create matrix X and the corresponding QR-matrix from source reference points
    for (int i=0; i<numReferencePoints; ++i) {
        vec3 srcPoint = refPointsSrc.at(i);
        // x-component
        X[3*i][0] = srcPoint.x;
        X[3*i][1] = srcPoint.y;
        X[3*i][2] = srcPoint.z;
        X[3*i][3] = 1;
        // y-component
        X[3*i+1][4] = srcPoint.x;
        X[3*i+1][5] = srcPoint.y;
        X[3*i+1][6] = srcPoint.z;
        X[3*i+1][7] = 1;
        // z-component
        X[3*i+2][8] = srcPoint.x;
        X[3*i+2][9] = srcPoint.y;
        X[3*i+2][10] = srcPoint.z;
        X[3*i+2][11] = 1;
    }
    JAMA::QR<float> QR_Matrix(X);

    if (!QR_Matrix.isFullRank()) {
        std::cout << "Transformation matrix could not be calculated. The reference points seem to be invalid.";
        return false;
    }

    // create y-vector from dest reference points
    for (int i=0; i<numReferencePoints; ++i) {
        vec3 destPoint = refPointsDest.at(i);
        y[3*i] = destPoint.x;
        y[3*i+1] = destPoint.y;
        y[3*i+2] = destPoint.z;
    }

    // least-squares solution
    c = QR_Matrix.solve(y);

    // copy solution to mat4
    mat4 M = mat4::createIdentity();
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

    std::cout << "\nRegistration matrix M mapping from source to dest volume coordinates: \n";
    std::cout << M << "\n";

    std::cout << "Deviation vectors (should be zero): \n";
    for (int i=0; i<numReferencePoints; ++i) {
        vec3 pSource = refPointsSrc.at(i);
        vec3 pDest = refPointsDest.at(i);
        vec3 pTrans = (M*vec4(pSource.x, pSource.y, pSource.z, 1.f)).xyz();
        std::cout << "M*Pi'-Pi: " << (pTrans - pDest) << "\n";
    }

    std::string input = "";
    while (input == "") {
        std::cout << "\nWrite registration matrix to source volume '" << srcVolFilename << "' (y or n)? ";
        std::cin >> input;
        if (input != "y" && input != "n")
            input = "";
        std::cout << "\n";
    }

    if (input == "y") {
        srcVol->setPhysicalToWorldMatrix(M);
        serializer->write(srcVolFilename, srcVol);
    }

    return true;

}
#else
bool CommandRegistrationAffine::execute(const std::vector<std::string>& /*parameters*/) {
    std::cout << "ERROR: CommandRegistrationAffine requires the JAMA module";
    return false;
}
#endif

} // namespace voreen
