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

#include "voreen/core/datastructures/geometry/facegeometry.h"

#include <map>
#include "tgt/glmath.h"
#include "tgt/immediatemode/immediatemode.h"

#include "voreen/core/io/serialization/serialization.h"

namespace voreen {

FaceGeometry::FaceGeometry()
    : normalIsSet_(false)
{}

FaceGeometry* FaceGeometry::create() const {
    return new FaceGeometry();
}

size_t FaceGeometry::getVertexCount() const {
    return vertices_.size();
}

bool FaceGeometry::empty() const {
    return (getVertexCount() == 0);
}

void FaceGeometry::addVertex(const VertexGeometry& vertex) {
    vertices_.push_back(vertex);
}

const VertexGeometry& FaceGeometry::getVertex(size_t index) const {
    return vertices_.at(index);
}

VertexGeometry& FaceGeometry::getVertex(size_t index) {
    return vertices_.at(index);
}


void FaceGeometry::setFaceNormal(const tgt::vec3& normal) {
    normalIsSet_ = true;
    normal_ = normal;
}

void FaceGeometry::clear() {
    vertices_.clear();
    normalIsSet_ = false;
}

FaceGeometry::iterator FaceGeometry::begin() {
    return vertices_.begin();
}

FaceGeometry::const_iterator FaceGeometry::begin() const {
    return vertices_.begin();
}

FaceGeometry::iterator FaceGeometry::end() {
    return vertices_.end();
}

FaceGeometry::const_iterator FaceGeometry::end() const {
    return vertices_.end();
}

VertexGeometry& FaceGeometry::operator[](size_t index) {
    return vertices_[index];
}

void FaceGeometry::render() const {
    //MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    //MatStack.pushMatrix();
    //MatStack.multMatrix(getTransformationMatrix());
    if(getVertexCount() >= 3) {
        GLint polygonMode[2];
        glGetIntegerv(GL_POLYGON_MODE, &polygonMode[0]);

        // values should be the same since only GL_FRONT_AND_BACK is supported for glPolygonMode since OpenGL 3.x
        if (polygonMode[0] == GL_LINE || polygonMode[1] == GL_LINE)
            IMode.begin(tgt::ImmediateMode::LINE_LOOP);
        else //use polygon in case of GL_POINTS or GL_FILL
            IMode.begin(tgt::ImmediateMode::POLYGON);

        for (const_iterator it = begin(); it != end(); it++) {
            IMode.color(it->getColor());
            IMode.texcoord(it->getTexCoords());

            if (normalIsSet_)
                IMode.normal(normal_);
            else if (it->isNormalDefined())
                IMode.normal(it->getNormal());

            IMode.vertex(it->getCoords());
        }
        IMode.end();
    }
    else {
        IMode.begin(tgt::ImmediateMode::LINES);
        for (const_iterator it = begin(); it != end(); it++) {
            IMode.color(it->getColor());
            IMode.texcoord(it->getTexCoords());

            if (normalIsSet_)
                IMode.normal(normal_);

            IMode.vertex(it->getCoords());
        }
        IMode.end();
    }

    IMode.color(1.f, 1.f, 1.f, 1.f);

    //MatStack.popMatrix();
}

void FaceGeometry::clip(const tgt::plane& clipPlane, double epsilon) {
    tgtAssert(epsilon >= 0.0, "negative epsilon");
    if (vertices_.size() < 2)
        return;

    // Since the clipped face vertices are appended to the vertex list, remember the current vertex list size.
    size_t vertexCount = vertices_.size();
    double lastDistance = vertices_.at(0).getDistanceToPlane(clipPlane, epsilon);

    // Process face edges...
    for (size_t i = 0; i < vertexCount; ++i) {
        double distance = vertices_.at((i + 1) % vertexCount).getDistanceToPlane(clipPlane, epsilon);

        // Keep both vertices?
        if (lastDistance <= 0 && distance <= 0) {
            // If processing the first edge, insert first vertex...
            if (i == 0)
                vertices_.push_back(vertices_.at(i));

            // If NOT processing the last edge, insert second vertex...
            if (i < (vertexCount - 1))
                vertices_.push_back(vertices_.at(i + 1));
        }
        // Discard first vertex, but keep second vertex?
        else if (lastDistance > 0 && distance <= 0) {
            // If NOT clipplane intersection vertex and second vertex are equal, insert clipplane intersection vertex...
            VertexGeometry intersectionVertex = VertexGeometry::interpolate(vertices_.at(i), vertices_.at((i + 1) % vertexCount), lastDistance / (lastDistance - distance));
            if (!vertices_.at((i + 1) % vertexCount).equals(intersectionVertex, epsilon))
                vertices_.push_back(intersectionVertex);

            // If NOT processing the last edge, insert second vertex...
            if (i < (vertexCount - 1))
                vertices_.push_back(vertices_.at(i + 1));
        }
        // Keep first vertex, but discard second vertex?
        else if (lastDistance <= 0 && distance > 0) {
            // If processing the first edge, insert first vertex...
            if (i == 0)
                vertices_.push_back(vertices_.at(i));

            // If NOT clipplane intersection vertex and first vertex are equal, insert clipplane intersection vertex...
            VertexGeometry intersectionVertex = VertexGeometry::interpolate(vertices_.at(i), vertices_.at((i + 1) % vertexCount), lastDistance / (lastDistance - distance));
            if (!vertices_.at(i).equals(intersectionVertex, epsilon))
                vertices_.push_back(intersectionVertex);
        }

        lastDistance = distance;
    }

    // Remove unclipped face vertices from the vertex list...
    vertices_.erase(vertices_.begin(), vertices_.begin() + vertexCount);
}

bool FaceGeometry::equals(const FaceGeometry& face, double epsilon /*= 1e-6*/) const {
    if (getVertexCount() != face.getVertexCount())
        return false;
    for (size_t i=0; i<getVertexCount(); i++) {
        if (!vertices_[i].equals(face.vertices_[i], epsilon))
            return false;
    }
    return true;
}

tgt::Bounds FaceGeometry::getBoundingBox() const {
    tgt::Bounds bounds;
    for (size_t i = 0; i < vertices_.size(); i++)
        bounds.addPoint(vertices_[i].getCoords());
    return bounds;
}

void FaceGeometry::serialize(Serializer& s) const {
    s.serialize("vertices", vertices_);
    if (normalIsSet_)
        s.serialize("normal", normal_);
}

void FaceGeometry::deserialize(Deserializer& s) {
    s.deserialize("vertices", vertices_);
    try {
        s.deserialize("normal", normal_);
        normalIsSet_ = true;
    }
    catch (...) {
        s.removeLastError();
        normalIsSet_ = false;
    }
}

bool FaceGeometry::operator==(const FaceGeometry& rhs) const {
    return (this->getVertexCount() == rhs.getVertexCount());
}

bool FaceGeometry::operator!=(const FaceGeometry& rhs) const {
    return !(*this == rhs);
}

} // namespace
