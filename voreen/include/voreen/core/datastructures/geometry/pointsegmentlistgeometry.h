/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_SEGMENTLISTGEOMETRY_H
#define VRN_SEGMENTLISTGEOMETRY_H

#include "voreen/core/datastructures/geometry/geometry.h"
#include "tgt/glmath.h"
#include <vector>

#include "voreen/core/io/serialization/xmlserializer.h"
#include "voreen/core/io/serialization/xmldeserializer.h"

namespace voreen {

template<class T>
class PointSegmentListGeometry : public Geometry {

public:

    PointSegmentListGeometry() :
      Geometry(),
      numPoints_(0)
    { }

    virtual ~PointSegmentListGeometry() {
    }

    void setData(std::vector< std::vector<T> > segmentList) {

        segmentList_ = segmentList;
        pointList_.clear();
        numPoints_ = 0;
        for (size_t i=0; i<segmentList_.size(); ++i) {
            numPoints_ += static_cast<int>(segmentList_[i].size());
        }
    }

    void addSegment(std::vector<T> segment) {
        segmentList_.push_back(segment);
        pointList_.clear();
        numPoints_ += static_cast<int>(segment.size());
    }

    const std::vector<T>& getSegment(int i) const {
        tgtAssert(i >= 0 && i < getNumSegments(), "invalid index");
        return segmentList_.at(i);
    }

    const std::vector< std::vector<T> >& getData() const { return segmentList_; }

    const std::vector<T>& getPoints() const {

        // generate point list, if not present
        if (pointList_.empty() && numPoints_ > 0) {
            for (size_t i=0; i<segmentList_.size(); ++i) {
                pointList_.insert(pointList_.end(), segmentList_[i].begin(), segmentList_[i].end());
            }
        }

        return pointList_;
    }

    int getNumSegments() const { return static_cast<int>(segmentList_.size()); }

    int getNumPoints() const { return numPoints_; }

    virtual void render() const {
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.pushMatrix();
        MatStack.multMatrix(getTransformationMatrix());
        glBegin(GL_POINTS);
        for (size_t i=0; i < segmentList_.size(); ++i){
            for (size_t j=0; j < segmentList_[i].size(); ++j){
                // assuming data type stored in the segment list is compatible to tgt::vertex
                // if not: template instantiation will fail (compile error)
                tgt::vertex(segmentList_[i][j]);
            }
        }
        glEnd();

        MatStack.popMatrix();
    }

    /**
     * Returns true, if the passed Geometry is a PointSegmentListGeometry<T>
     * and all its segments/vertices are equal to this one's.
     *
     * @see Geometry::equals
     */
    virtual bool equals(const Geometry* geometry, double epsilon = 1e-6) const {
        const PointSegmentListGeometry<T>* pointSegListGeom = dynamic_cast<const PointSegmentListGeometry<T>* >(geometry);
        if (!pointSegListGeom)
            return false;
        if (getNumSegments() != pointSegListGeom->getNumSegments())
            return false;
        if (getNumPoints() != pointSegListGeom->getNumPoints())
            return false;

        float epsilonSq = static_cast<float>(epsilon*epsilon);
        for (size_t i=0; i < segmentList_.size(); ++i){
            if (segmentList_[i].size() != pointSegListGeom->segmentList_[i].size())
                return false;
            for (size_t j=0; j < segmentList_[i].size(); ++j){
                // assuming data type stored in the point list is a tgt vector type
                // if not: template instantiation will fail (compile error)
                if (tgt::lengthSq(segmentList_[i][j] - pointSegListGeom->segmentList_[i][j]) > epsilonSq)
                    return false;
            }
        }

        return true;
    }

    virtual void clip(const tgt::plane& clipPlane, double epsilon = 1e-6) {
        tgt::plane plane = clipPlane.transform(getInvertedTransformationMatrix());

        tgtAssert(epsilon, "negative epsilon");
        std::vector< std::vector<T> > clippedSegmentList;
        for (size_t i=0; i<segmentList_.size(); i++) {
            std::vector<T> clippedSegment;
            for (size_t j=0; j<segmentList_[i].size(); j++) {
                double distance = plane.distance(segmentList_[i][j]);
                if (distance <= epsilon)
                    clippedSegment.push_back(segmentList_[i][j]);
            }
            clippedSegmentList.push_back(clippedSegment);
        }
        segmentList_ = clippedSegmentList;
    }

    virtual tgt::Bounds getBoundingBox(bool transformed = true) const {
        tgt::Bounds bounds;
        for (size_t i=0; i<segmentList_.size(); i++)
            for (size_t j=0; j<segmentList_[i].size(); j++)
                bounds.addPoint(segmentList_[i][j]);

        if(transformed)
            return bounds.transform(getTransformationMatrix());
        else
            return bounds;
    }

    virtual void serialize(Serializer& s) const {
        s.serialize("segments", segmentList_);
        Geometry::serialize(s);
    }

    virtual void deserialize(Deserializer& s) {
        s.deserialize("segments", segmentList_);
        Geometry::deserialize(s);
    }

protected:
    // contains a list of segments, each segment consists of points
    std::vector< std::vector<T> > segmentList_;

    // contains the flattened segment list: a list of all points
    // is generated on first access
    mutable std::vector<T> pointList_;

    int numPoints_;

};

class PointSegmentListGeometryVec3 : public PointSegmentListGeometry<tgt::vec3> {
    virtual Geometry* create() const { return new PointSegmentListGeometryVec3(); }
    virtual std::string getClassName() const { return "PointSegmentListGeometryVec3"; }
    virtual std::unique_ptr<Geometry> clone() const;

public:
    // Collect and copy segments from all PointSegmentListGeometryVec3s in the
    // provided geometry. In particular, segments from all
    // PointSegmentListGeometryVec3s stored in GeometrySequences are collected.
    //
    // Returns true if all (potentially recursively) provided geometry were
    // compatible and could be processed, and false otherwise.
    bool collectSegmentsFromGeometry(const Geometry& other);
};

} // namespace

#endif // VRN_SEGMENTLISTGEOMETRY_H
