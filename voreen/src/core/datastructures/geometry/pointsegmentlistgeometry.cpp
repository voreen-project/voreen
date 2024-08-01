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

#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "voreen/core/datastructures/geometry/geometry.h"
#include "voreen/core/datastructures/geometry/geometrysequence.h"

namespace voreen {

static bool getSeedListsFromGeom(const Geometry* geom, const tgt::mat4& transform, PointSegmentListGeometryVec3& seeds);

static void getSeedListFromPointSegmentList(const PointSegmentListGeometryVec3& seedList, const tgt::mat4& transform, PointSegmentListGeometryVec3& seeds) {
    auto transformMat = transform * seedList.getTransformationMatrix();
    for (int j=0; j<seedList.getNumSegments(); j++) {
        std::vector<tgt::vec3> points;
        for(auto& vox : seedList.getSegment(j)) {
            points.push_back(transformMat.transform(vox));
        }
        seeds.addSegment(points);
    }
}
static bool getSeedListFromSequence(const GeometrySequence& geoms, const tgt::mat4& transform, PointSegmentListGeometryVec3& seeds) {
    bool allFine = true;
    for (int j=0; j<geoms.getNumGeometries(); j++) {
        allFine &= getSeedListsFromGeom(geoms.getGeometry(j), geoms.getTransformationMatrix(), seeds);
    }
    return allFine;
}
static bool getSeedListsFromGeom(const Geometry* geom, const tgt::mat4& transform, PointSegmentListGeometryVec3& seeds) {
    if(auto* seedList = dynamic_cast<const PointSegmentListGeometryVec3* >(geom)) {
        getSeedListFromPointSegmentList(*seedList, transform, seeds);
        return true;
    } else if (auto* geoms = dynamic_cast<const GeometrySequence* >(geom)) {
        return getSeedListFromSequence(*geoms, transform, seeds);
    } else {
        return false;
    }
}
bool PointSegmentListGeometryVec3::collectSegmentsFromGeometry(const Geometry& geom) {
    auto mat = tgt::mat4::createIdentity();
    return getSeedListsFromGeom(&geom, mat, *this);
}

std::unique_ptr<Geometry> PointSegmentListGeometryVec3::clone() const {
    PointSegmentListGeometryVec3* clone = new PointSegmentListGeometryVec3();
    clone->segmentList_ = segmentList_;
    clone->pointList_ = pointList_;
    clone->numPoints_ = numPoints_;
    clone->setTransformationMatrix(getTransformationMatrix());
    return std::unique_ptr<Geometry>(clone);
}
}
