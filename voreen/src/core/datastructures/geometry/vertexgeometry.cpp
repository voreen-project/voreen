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

#include "voreen/core/datastructures/geometry/vertexgeometry.h"

#include "tgt/glmath.h"
#include "voreen/core/io/serialization/serialization.h"

#include <limits>

namespace voreen {

VertexGeometry::VertexGeometry(const tgt::vec3& coords, const tgt::vec3& texcoords, const tgt::vec4& color)
    : coords_(coords)
    , texcoords_(texcoords)
    , color_(color)
    , normalIsSet_(false)
    , normal_(tgt::vec3(0.f))
{}

VertexGeometry::VertexGeometry(const tgt::vec3& coords, const tgt::vec3& texcoords, const tgt::vec4& color, const tgt::vec3& normal)
    : coords_(coords)
    , texcoords_(texcoords)
    , color_(color)
    , normalIsSet_(true)
    , normal_(normal)
{}

VertexGeometry* VertexGeometry::create() const {
    return new VertexGeometry();
}

tgt::vec3 VertexGeometry::getCoords() const {
    return coords_;
}

void VertexGeometry::setCoords(const tgt::vec3& coords) {
    coords_ = coords;
}

tgt::vec3 VertexGeometry::getTexCoords() const {
    return texcoords_;
}

void VertexGeometry::setTexCoords(const tgt::vec3& texcoords) {
    texcoords_ = texcoords;
}

tgt::vec4 VertexGeometry::getColor() const {
    return color_;
}

void VertexGeometry::setColor(const tgt::vec4& color) {
    color_ = color;
}

void VertexGeometry::setColor(const tgt::vec3& color) {
    color_ = tgt::vec4(color, 1);
}

tgt::vec3 VertexGeometry::getNormal() const {
    tgtAssert(normalIsSet_, "Tried to access the normal when it was not set");
    return normal_;
}

void VertexGeometry::setNormal(const tgt::vec3& normal) {
    normal_ = normal;
    normalIsSet_ = true;
}

void VertexGeometry::removeNormal() {
    normal_ = tgt::vec3(0.f);
    normalIsSet_ = false;
}

bool VertexGeometry::isNormalDefined() const {
    return normalIsSet_;
}

double VertexGeometry::getLength() const {
    return tgt::length(coords_);
}

double VertexGeometry::getDistanceToPlane(const tgt::plane& plane, double epsilon) const {
    double distance = tgt::dot(plane.n, coords_) - plane.d;
    if (std::abs(distance) <= epsilon)
        return 0;
    else
        return distance;
}

double VertexGeometry::getDistance(const VertexGeometry& vertex) const {
    return tgt::length(vertex.coords_ - coords_);
}

void VertexGeometry::render() const {
    //MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    //MatStack.pushMatrix();
    //MatStack.multMatrix(getTransformationMatrix());
    glBegin(GL_POINTS);

    glColor4fv(color_.elem);
    if (normalIsSet_)
        tgt::normal(normal_);
    tgt::texCoord(texcoords_);
    tgt::vertex(coords_);

    glEnd();

    glColor4f(1.f, 1.f, 1.f, 1.f);
    //MatStack.popMatrix();
}

void VertexGeometry::combine(const VertexGeometry& vertex) {
    interpolate(vertex, 0.5);
}

VertexGeometry VertexGeometry::combine(const VertexGeometry& vertex1, const VertexGeometry& vertex2) {
    return interpolate(vertex1, vertex2, 0.5);
}

void VertexGeometry::interpolate(const VertexGeometry& vertex, double t) {
    coords_ += (vertex.coords_ - coords_) * static_cast<tgt::vec3::ElemType>(t);
    texcoords_ += (vertex.texcoords_ - texcoords_) * static_cast<tgt::vec3::ElemType>(t);
    color_ += (vertex.color_ - color_) * static_cast<tgt::vec4::ElemType>(t);
    normal_ += (vertex.normal_ - normal_) * static_cast<tgt::vec4::ElemType>(t);
}

VertexGeometry VertexGeometry::interpolate(const VertexGeometry& vertex1, const VertexGeometry& vertex2, double t) {
    VertexGeometry result = vertex1;
    result.interpolate(vertex2, t);
    return result;
}

void VertexGeometry::clip(const tgt::plane& clipPlane, double epsilon /*= 1e-6*/) {
    if (getDistanceToPlane(clipPlane, epsilon) > 0.0) {
        // invalidate vertex
        coords_ = tgt::vec3(std::numeric_limits<float>::quiet_NaN());
        texcoords_ = tgt::vec3(std::numeric_limits<float>::quiet_NaN());
        color_ = tgt::vec4(std::numeric_limits<float>::quiet_NaN());
        normal_ = tgt::vec3(std::numeric_limits<float>::quiet_NaN());
        normalIsSet_ = false;
    }
}

bool VertexGeometry::equals(const VertexGeometry& vertex, double epsilon) const {
    tgtAssert(epsilon >= 0.0, "negative epsilon");
    if (std::abs(getDistance(vertex)) > epsilon)
        return false;

    // also compare texcoords, color, and normal
    float epsilonSq = static_cast<float>(epsilon*epsilon);
    if (tgt::lengthSq(texcoords_ - vertex.texcoords_) > epsilonSq)
        return false;
    if (tgt::lengthSq(color_ - vertex.color_) > epsilonSq)
        return false;
    if (isNormalDefined() != vertex.isNormalDefined() ||
        (isNormalDefined() && tgt::lengthSq(normal_ - vertex.normal_) > epsilonSq))
        return false;

    return true;
}

tgt::Bounds VertexGeometry::getBoundingBox() const {
    return tgt::Bounds(coords_);
}

bool VertexGeometry::operator==(const VertexGeometry& vertex) const {
    return this->equals(vertex);
}

bool VertexGeometry::operator!=(const VertexGeometry& vertex) const {
    return !(*this == vertex);
}

void VertexGeometry::serialize(Serializer& s) const {
    s.serialize("coords", coords_);
    s.serialize("texcoords", texcoords_);
    s.serialize("color", color_);
    if (normalIsSet_)
        s.serialize("normal", normal_);
}

void VertexGeometry::deserialize(Deserializer& s) {
    s.deserialize("coords", coords_);
    s.deserialize("texcoords", texcoords_);
    s.deserialize("color", color_);
    try {
        s.deserialize("normal", normal_);
        normalIsSet_ = true;
    }
    catch (...) {
        normalIsSet_ = false;
        s.removeLastError();
    }
}

} // namespace
