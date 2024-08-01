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

#include "voreen/core/datastructures/geometry/meshlistgeometry.h"

#include "tgt/glmath.h"
#include "tgt/matrixstack.h"

#include "voreen/core/io/serialization/serialization.h"

namespace voreen {

MeshListGeometry::MeshListGeometry()
    : Geometry()
{}

MeshListGeometry::~MeshListGeometry() {
}

Geometry* MeshListGeometry::create() const {
    return new MeshListGeometry();
}

size_t MeshListGeometry::getMeshCount() const {
    return meshes_.size();
}

size_t MeshListGeometry::getFaceCount() const {
    size_t result = 0;
    for (const_iterator it = begin(); it != end(); ++it)
        result += it->getFaceCount();
    return result;
}

size_t MeshListGeometry::getVertexCount() const {
    size_t result = 0;
    for (const_iterator it = begin(); it != end(); ++it)
        result += it->getVertexCount();
    return result;
}

bool MeshListGeometry::empty() const {
    return (getVertexCount() == 0);
}

void MeshListGeometry::addMesh(const MeshGeometry& mesh) {
    meshes_.push_back(mesh);
}

const MeshGeometry& MeshListGeometry::getMesh(size_t index) const {
    tgtAssert(index < meshes_.size(), "Invalid index");
    return meshes_.at(index);
}

MeshGeometry& MeshListGeometry::getMesh(size_t index) {
    tgtAssert(index < meshes_.size(), "Invalid index");
    return meshes_.at(index);
}

void MeshListGeometry::clear() {
    meshes_.clear();
}

MeshListGeometry::iterator MeshListGeometry::begin() {
    return meshes_.begin();
}

MeshListGeometry::const_iterator MeshListGeometry::begin() const {
    return meshes_.begin();
}

MeshListGeometry::iterator MeshListGeometry::end() {
    return meshes_.end();
}

MeshListGeometry::const_iterator MeshListGeometry::end() const {
    return meshes_.end();
}

MeshGeometry& MeshListGeometry::operator[](size_t index) {
    tgtAssert(index < meshes_.size(), "Invalid index");
    return meshes_[index];
}

const MeshGeometry& MeshListGeometry::operator[](size_t index) const {
    tgtAssert(index < meshes_.size(), "Invalid index");
    return meshes_[index];
}

tgt::Bounds MeshListGeometry::getBoundingBox(bool transformed) const {
    tgt::Bounds bounds;
    for (size_t i=0; i<meshes_.size(); i++)
        bounds.addVolume(meshes_[i].getBoundingBox(transformed));

    if(transformed)
        return bounds.transform(getTransformationMatrix());
    else
        return bounds;
}

void MeshListGeometry::render() const {
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.multMatrix(getTransformationMatrix());

    for (const_iterator it = begin(); it != end(); ++it)
        it->render();

    MatStack.popMatrix();
}

void MeshListGeometry::clip(const tgt::plane& clipPlane, MeshListGeometry& closingFaces, double epsilon) {
    tgt::plane pl = clipPlane.transform(getInvertedTransformationMatrix());

    tgtAssert(epsilon >= 0.0, "negative epsilon");
    for (iterator it = begin(); it != end(); ++it) {
        MeshGeometry closingFace;
        it->clip(pl, closingFace, epsilon);
        if (!closingFace.empty())
            closingFaces.addMesh(closingFace);
    }

    closingFaces.setTransformationMatrix(getTransformationMatrix());
}

void MeshListGeometry::clip(const tgt::plane& clipPlane, double epsilon) {
    MeshListGeometry dummy;
    clip(clipPlane, dummy, epsilon);
}

bool MeshListGeometry::equals(const MeshListGeometry& meshList, double epsilon /*= 1e-6*/) const {
    if (getMeshCount() != meshList.getMeshCount())
        return false;
    for (size_t i=0; i<getMeshCount(); i++) {
        if (!meshes_[i].equals(meshList.meshes_[i], epsilon))
            return false;
    }
    return true;
}

bool MeshListGeometry::equals(const Geometry* geometry, double epsilon /*= 1e-6*/) const {
    const MeshListGeometry* meshListGeometry = dynamic_cast<const MeshListGeometry*>(geometry);
    if (!meshListGeometry)
        return false;
    else
        return equals(*meshListGeometry, epsilon);
}

void MeshListGeometry::serialize(Serializer& s) const {
    s.serialize("meshes", meshes_);
    Geometry::serialize(s);
}

void MeshListGeometry::deserialize(Deserializer& s) {
    s.deserialize("meshes", meshes_);
    Geometry::deserialize(s);
}

} // namespace
