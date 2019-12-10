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

#include "voreen/core/datastructures/geometry/geometrysequence.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"

#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/utils/hashing.h"
#include "tgt/glmath.h"
#include "tgt/matrixstack.h"
#include "tgt/immediatemode/immediatemode.h"

#include <sstream>

namespace voreen {

const std::string GeometrySequence::loggerCat_("voreen.GeometrySequence");

GeometrySequence::GeometrySequence(bool takeOwnership)
    : Geometry()
    , takeOwnership_(takeOwnership)
{}

GeometrySequence::~GeometrySequence() {
    if(takeOwnership_) {
        for(size_t i=0; i<geometries_.size(); i++)
            delete geometries_[i];
    }

    geometries_.clear();
}

void GeometrySequence::addGeometry(Geometry* g) {
    geometries_.push_back(g);
}

size_t GeometrySequence::getNumGeometries() const {
    return geometries_.size();
}

const Geometry* GeometrySequence::getGeometry(size_t i) const {
    return geometries_[i];
}

bool GeometrySequence::equals(const Geometry* geometry, double epsilon) const {
    tgtAssert(geometry, "null pointer passed");
    const GeometrySequence* gs = dynamic_cast<const GeometrySequence*>(geometry);
    if(gs) {
        if(geometries_.size() != gs->geometries_.size())
            return false;

        for(size_t i=0; i<geometries_.size(); i++) {
            if(!geometries_[i]->equals(gs->geometries_[i], epsilon))
                return false;
        }
        return true;
    }
    else
        return false;
}

void GeometrySequence::clip(const tgt::plane& clipPlane, double epsilon) {
    for(size_t i=0; i<geometries_.size(); i++)
        geometries_[i]->clip(clipPlane, epsilon);
}

void GeometrySequence::render() const {
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.multMatrix(getTransformationMatrix());

    for(size_t i=0; i<geometries_.size(); i++) {
        if(dynamic_cast<TriangleMeshGeometryBase*>(geometries_[i])) {
            MatStack.pushMatrix();
            MatStack.multMatrix(geometries_[i]->getTransformationMatrix());
            IMode.setMatstackUniforms();
            geometries_[i]->render();
            MatStack.popMatrix();
        }
        else
            geometries_[i]->render();
    }

    MatStack.popMatrix();
}

tgt::Bounds GeometrySequence::getBoundingBox(bool transformed) const {
    tgt::Bounds b;
    for(size_t i=0; i<geometries_.size(); i++)
        b.addVolume(geometries_[i]->getBoundingBox(true));//old version transformed

    if(transformed)
        return b.transform(getTransformationMatrix());
    else
        return b;
}

void GeometrySequence::serialize(Serializer& s) const {
    Geometry::serialize(s);
    s.serialize("geometries", geometries_);
}

void GeometrySequence::deserialize(Deserializer& s) {
    Geometry::deserialize(s);
    s.deserialize("geometries", geometries_);
}

} // namespace
