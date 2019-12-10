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

#include "voreen/core/datastructures/geometry/meshgeometry.h"

#include "tgt/glmath.h"
#include "tgt/matrixstack.h"

#include "voreen/core/io/serialization/serialization.h"

using tgt::vec3;
using tgt::vec4;

namespace voreen {

MeshGeometry::MeshGeometry()
    : Geometry()
{}

MeshGeometry::~MeshGeometry() {
}

Geometry* MeshGeometry::create() const {
    return new MeshGeometry();
}

MeshGeometry MeshGeometry::createCube(vec3 coordLlf,
                                      vec3 coordUrb,
                                      vec3 texLlf,
                                      vec3 texUrb,
                                      vec3 colorLlf,
                                      vec3 colorUrb,
                                      float alpha)
{
    FaceGeometry top, front, left, back, right, bottom;
    createCubeFaces(top, front, left, back, right, bottom, coordLlf, coordUrb, texLlf, texUrb, colorLlf, colorUrb, alpha);

    MeshGeometry mesh;
    mesh.addFace(top);
    mesh.addFace(front);
    mesh.addFace(left);
    mesh.addFace(back);
    mesh.addFace(right);
    mesh.addFace(bottom);

    return mesh;
}

voreen::MeshGeometry MeshGeometry::createCube(tgt::vec3 coordLlf,
                                              tgt::vec3 coordUrb,
                                              tgt::vec3 texLlf,
                                              tgt::vec3 texUrb,
                                              tgt::vec3 colorLlf,
                                              tgt::vec3 colorUrb,
                                              tgt::vec3 normalTop,
                                              tgt::vec3 normalFront,
                                              tgt::vec3 normalLeft,
                                              tgt::vec3 normalBack,
                                              tgt::vec3 normalRight,
                                              tgt::vec3 normalBottom,
                                              float alpha)
{
    FaceGeometry top, front, left, back, right, bottom;
    createCubeFaces(top, front, left, back, right, bottom, coordLlf, coordUrb, texLlf, texUrb, colorLlf, colorUrb, alpha);
    top.setFaceNormal(normalTop);
    front.setFaceNormal(normalFront);
    left.setFaceNormal(normalLeft);
    back.setFaceNormal(normalBack);
    right.setFaceNormal(normalRight);
    bottom.setFaceNormal(normalBottom);

    MeshGeometry mesh;
    mesh.addFace(top);
    mesh.addFace(front);
    mesh.addFace(left);
    mesh.addFace(back);
    mesh.addFace(right);
    mesh.addFace(bottom);

    return mesh;
}

size_t MeshGeometry::getFaceCount() const {
    return faces_.size();
}

size_t MeshGeometry::getVertexCount() const {
    size_t result = 0;
    for (const_iterator it = begin(); it != end(); ++it)
        result += it->getVertexCount();
    return result;
}

bool MeshGeometry::empty() const {
    return (getVertexCount() == 0);
}

void MeshGeometry::addFace(const FaceGeometry& face) {
    faces_.push_back(face);
}

const FaceGeometry& MeshGeometry::getFace(size_t index) const {
    tgtAssert(index < faces_.size(), "Invalid face index");
    return faces_.at(index);
}

FaceGeometry& MeshGeometry::getFace(size_t index) {
    tgtAssert(index < faces_.size(), "Invalid face index");
    return faces_.at(index);
}

void MeshGeometry::clear() {
    faces_.clear();
}

MeshGeometry::iterator MeshGeometry::begin() {
    return faces_.begin();
}

MeshGeometry::const_iterator MeshGeometry::begin() const {
    return faces_.begin();
}

MeshGeometry::iterator MeshGeometry::end() {
    return faces_.end();
}

MeshGeometry::const_iterator MeshGeometry::end() const {
    return faces_.end();
}

FaceGeometry& MeshGeometry::operator[] (size_t index) {
    tgtAssert(index < faces_.size(), "Invalid face index");
    return faces_[index];
}

void MeshGeometry::render() const {
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.multMatrix(getTransformationMatrix());

    for (const_iterator it = begin(); it != end(); ++it)
        it->render();

    MatStack.popMatrix();
}

void MeshGeometry::clip(const tgt::plane& clipPlane, MeshGeometry& closingMesh, double epsilon) {
    tgtAssert(epsilon >= 0.0, "negative epsilon");

    tgt::plane pl = clipPlane.transform(getInvertedTransformationMatrix());

    // Clip all faces...
    for (iterator it = begin(); it != end(); ++it)
        it->clip(pl, epsilon);

    // Remove empty faces...
    for (size_t i = 0; i < faces_.size(); ++i) {
        // Is face empty?
        if (faces_.at(i).getVertexCount() < 3)
            faces_.erase(faces_.begin() + i--);
    }

    // Close convex polyhedron if necessary...
    typedef std::pair<VertexGeometry, VertexGeometry> EdgeType;
    typedef std::vector<EdgeType> EdgeListType;
    typedef std::vector<VertexGeometry> VertexListType;

    EdgeListType edgeList;
    FaceGeometry closingFace;

    // Search all face edges on the clipping plane...
    for (size_t i = 0; i < faces_.size(); ++i) {
        FaceGeometry face = faces_.at(i);

        VertexListType verticesOnClipplane;

        for (size_t j = 0; j < face.getVertexCount(); ++j) {
            if (face.getVertex(j).getDistanceToPlane(pl, epsilon) == 0)
                verticesOnClipplane.push_back(face.getVertex(j));

            // Is face in the same plane as the clipping plane?
            if (verticesOnClipplane.size() > 2)
                break;
        }

        // Does one face edge corresponds with clipping plane?
        if (verticesOnClipplane.size() == 2)
            edgeList.push_back(std::make_pair(verticesOnClipplane[0], verticesOnClipplane[1]));
    }

    // Is closing necessary?
    if (edgeList.size() > 1) {
        // Sort edges to produce contiguous vertex order...
        bool reverseLastEdge = false;
        for (size_t i = 0; i < edgeList.size() - 1; ++i) {
            for (size_t j = i + 1; j < edgeList.size(); ++j) {
                VertexGeometry connectionVertex;
                if (reverseLastEdge)
                    connectionVertex = edgeList.at(i).first;
                else
                    connectionVertex = edgeList.at(i).second;

                if (edgeList.at(j).first.equals(connectionVertex, epsilon)) {
                    std::swap(edgeList.at(i + 1), edgeList.at(j));
                    reverseLastEdge = false;
                    break;
                }
                else if (edgeList.at(j).second.equals(connectionVertex, epsilon)) {
                    std::swap(edgeList.at(i + 1), edgeList.at(j));
                    reverseLastEdge = true;
                    break;
                }
            }
        }

        // Convert sorted edge list to sorted vertex list...
        VertexListType closingFaceVertices;
        for (size_t i = 0; i < edgeList.size(); ++i) {
            bool reverseEdge = i != 0 && !closingFaceVertices.at(closingFaceVertices.size() - 1).equals(edgeList.at(i).first, epsilon);

            VertexGeometry first = (reverseEdge ? edgeList.at(i).second : edgeList.at(i).first);
            VertexGeometry second = (reverseEdge ? edgeList.at(i).first : edgeList.at(i).second);

            if (i == 0)
                closingFaceVertices.push_back(first);
            else
                closingFaceVertices.at(closingFaceVertices.size() - 1).combine(first);

            if (i < (edgeList.size() - 1))
                closingFaceVertices.push_back(second);
            else
                closingFaceVertices[0].combine(second);
        }

        // Convert vertex order to counter clockwise if necessary...
        vec3 closingFaceNormal(0, 0, 0);
        for (size_t i = 0; i < closingFaceVertices.size(); ++i)
            closingFaceNormal += tgt::cross(closingFaceVertices[i].getCoords(), closingFaceVertices[(i + 1) % closingFaceVertices.size()].getCoords());
        closingFaceNormal = tgt::normalize(closingFaceNormal);

        if (tgt::dot(pl.n, closingFaceNormal) < 0)
            std::reverse(closingFaceVertices.begin(), closingFaceVertices.end());

        // Close convex polyhedron...
        for (VertexListType::iterator it = closingFaceVertices.begin(); it != closingFaceVertices.end(); ++it) {
            // TODO(b_bolt01): Remove debug message...
            //std::cout << " cfv " << it->getCoords() << std::endl;
            closingFace.addVertex(*it);
        }
        addFace(closingFace);
    }

    // If there is only the clipplane left, erase it also...
    if (faces_.size() == 1)
        faces_.clear();

    if (closingFace.getVertexCount() > 0)
        closingMesh.addFace(closingFace);

    closingMesh.setTransformationMatrix(getTransformationMatrix());
}

void MeshGeometry::clip(const tgt::plane& clipPlane, double epsilon) {
    MeshGeometry dummy;
    clip(clipPlane, dummy, epsilon);
}

bool MeshGeometry::equals(const MeshGeometry& mesh, double epsilon /*= 1e-6*/) const {
    if (getFaceCount() != mesh.getFaceCount())
        return false;

    if(getTransformationMatrix() != mesh.getTransformationMatrix()) {
        LWARNINGC("", "Matrices differ" << std::endl << getTransformationMatrix() << std::endl << mesh.getTransformationMatrix());
        return false;
    }

    for (size_t i=0; i<getFaceCount(); i++) {
        if (!faces_[i].equals(mesh.faces_[i], epsilon))
            return false;
    }
    return true;
}

bool MeshGeometry::equals(const Geometry* geometry, double epsilon /*= 1e-6*/) const {
    const MeshGeometry* meshGeometry = dynamic_cast<const MeshGeometry*>(geometry);
    if (!meshGeometry)
        return false;
    else
        return equals(*meshGeometry, epsilon);
}

tgt::Bounds MeshGeometry::getBoundingBox(bool transformed) const {
    tgt::Bounds bounds;
    for (size_t i=0; i<faces_.size(); i++)
        bounds.addVolume(faces_[i].getBoundingBox());

    if(transformed)
        return bounds.transform(getTransformationMatrix());
    else
        return bounds;
}

void MeshGeometry::serialize(Serializer& s) const {
    s.serialize("faces", faces_);
    Geometry::serialize(s);
}

void MeshGeometry::deserialize(Deserializer& s) {
    s.deserialize("faces", faces_);
    Geometry::deserialize(s);
}

void MeshGeometry::createCubeFaces(FaceGeometry& topFace, FaceGeometry& frontFace, FaceGeometry& leftFace,
                                   FaceGeometry& backFace, FaceGeometry& rightFace, FaceGeometry& bottomFace,
                                   tgt::vec3 coordLlf, tgt::vec3 coordUrb, tgt::vec3 texLlf,
                                   tgt::vec3 texUrb, tgt::vec3 colorLlf, tgt::vec3 colorUrb, float alpha)
{
    // expecting coordLlf < coordUrb
    if (coordLlf.x > coordUrb.x) {
        std::swap(coordLlf.x, coordUrb.x);
        std::swap(texLlf.x, texUrb.x);
        std::swap(colorLlf.x, colorUrb.x);
    }
    if (coordLlf.y > coordUrb.y) {
        std::swap(coordLlf.y, coordUrb.y);
        std::swap(texLlf.y, texUrb.y);
        std::swap(colorLlf.y, colorUrb.y);
    }
    if (coordLlf.z > coordUrb.z) {
        std::swap(coordLlf.z, coordUrb.z);
        std::swap(texLlf.z, texUrb.z);
        std::swap(colorLlf.z, colorUrb.z);
    }

    VertexGeometry llf(vec3(coordLlf.x, coordLlf.y, coordLlf.z), vec3(texLlf.x, texLlf.y, texLlf.z), vec4(colorLlf.x, colorLlf.y, colorLlf.z, alpha));
    VertexGeometry lrf(vec3(coordUrb.x, coordLlf.y, coordLlf.z), vec3(texUrb.x, texLlf.y, texLlf.z), vec4(colorUrb.x, colorLlf.y, colorLlf.z, alpha));
    VertexGeometry lrb(vec3(coordUrb.x, coordLlf.y, coordUrb.z), vec3(texUrb.x, texLlf.y, texUrb.z), vec4(colorUrb.x, colorLlf.y, colorUrb.z, alpha));
    VertexGeometry llb(vec3(coordLlf.x, coordLlf.y, coordUrb.z), vec3(texLlf.x, texLlf.y, texUrb.z), vec4(colorLlf.x, colorLlf.y, colorUrb.z, alpha));

    VertexGeometry ulb(vec3(coordLlf.x, coordUrb.y, coordUrb.z), vec3(texLlf.x, texUrb.y, texUrb.z), vec4(colorLlf.x, colorUrb.y, colorUrb.z, alpha));
    VertexGeometry ulf(vec3(coordLlf.x, coordUrb.y, coordLlf.z), vec3(texLlf.x, texUrb.y, texLlf.z), vec4(colorLlf.x, colorUrb.y, colorLlf.z, alpha));
    VertexGeometry urf(vec3(coordUrb.x, coordUrb.y, coordLlf.z), vec3(texUrb.x, texUrb.y, texLlf.z), vec4(colorUrb.x, colorUrb.y, colorLlf.z, alpha));
    VertexGeometry urb(vec3(coordUrb.x, coordUrb.y, coordUrb.z), vec3(texUrb.x, texUrb.y, texUrb.z), vec4(colorUrb.x, colorUrb.y, colorUrb.z, alpha));

    topFace.addVertex(urb);
    topFace.addVertex(urf);
    topFace.addVertex(ulf);
    topFace.addVertex(ulb);

    frontFace.addVertex(llf);
    frontFace.addVertex(ulf);
    frontFace.addVertex(urf);
    frontFace.addVertex(lrf);

    leftFace.addVertex(llf);
    leftFace.addVertex(llb);
    leftFace.addVertex(ulb);
    leftFace.addVertex(ulf);

    backFace.addVertex(urb);
    backFace.addVertex(ulb);
    backFace.addVertex(llb);
    backFace.addVertex(lrb);

    rightFace.addVertex(urb);
    rightFace.addVertex(lrb);
    rightFace.addVertex(lrf);
    rightFace.addVertex(urf);

    bottomFace.addVertex(llf);
    bottomFace.addVertex(lrf);
    bottomFace.addVertex(lrb);
    bottomFace.addVertex(llb);
}

} // namespace
