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

#ifndef VRN_POINTLISTGEOMETRY_H
#define VRN_POINTLISTGEOMETRY_H

#include "voreen/core/datastructures/geometry/geometry.h"
#include "tgt/vector.h"
#include "tgt/glmath.h"
#include "tgt/matrixstack.h"
#include <vector>

#include "voreen/core/io/serialization/xmlserializer.h"
#include "voreen/core/io/serialization/xmldeserializer.h"

namespace voreen {

template<class T>
class PointListGeometry : public Geometry {
public:
    PointListGeometry() : Geometry() {}
    virtual ~PointListGeometry() {}
    void addPoint(T point) {points_.push_back(point);}
    void removeLast() {points_.pop_back();}
    void clear() {points_.clear();}
    const std::vector<T>& getData() const { return points_; }
    void setData(std::vector<T> points) { points_ = points; }
    size_t getNumPoints() const { return points_.size(); }
    T getPoint(size_t i) const { return points_[i]; }

    typename std::vector<T>::iterator begin() { return points_.begin(); }
    typename std::vector<T>::const_iterator begin() const { return points_.begin(); }
    typename std::vector<T>::iterator end() { return points_.end(); }
    typename std::vector<T>::const_iterator end() const { return points_.end(); }

    T& operator[](size_t index) { return points_[index]; }
    const T& operator[](size_t index) const { return points_[index]; }

    virtual void render() const {
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.pushMatrix();
        MatStack.multMatrix(getTransformationMatrix());
        glBegin(GL_POINTS);
        for (size_t i=0; i < points_.size(); ++i){
            // assuming data type stored in the point list is compatible to tgt::vertex
            // if not: template instantiation will fail (compile error)
            tgt::vertex(points_[i]);
        }
        glEnd();

        MatStack.popMatrix();
    }

    /**
     * Returns true, if the passed Geometry is a PointListGeometry<T>
     * and all its vertices are equal to this one's.
     *
     * @see Geometry::equals
     */
    virtual bool equals(const Geometry* geometry, double epsilon = 1e-6) const {
        const PointListGeometry<T>* pointListGeom = dynamic_cast<const PointListGeometry<T>* >(geometry);
        if (!pointListGeom)
            return false;
        if (getNumPoints() != pointListGeom->getNumPoints())
            return false;

        if(getTransformationMatrix() != geometry->getTransformationMatrix())
            return false;

        float epsilonSq = static_cast<float>(epsilon*epsilon);
        for (size_t i=0; i<getNumPoints(); i++) {
            // assuming data type stored in the point list is a tgt vector type
            // if not: template instantiation will fail (compile error)
            if (tgt::lengthSq(points_[i] - pointListGeom->points_[i]) > epsilonSq)
                return false;
        }

        return true;
    }

    virtual void clip(const tgt::plane& clipPlane, double epsilon = 1e-6) {
        tgt::plane pl = clipPlane.transform(getInvertedTransformationMatrix());

        tgtAssert(epsilon, "negative epsilon");
        std::vector<T> clipped;
        for (size_t i=0; i<points_.size(); i++) {
            double distance = pl.distance(points_[i]);
            if (distance <= epsilon)
                clipped.push_back(points_[i]);
        }
        points_ = clipped;
    }

    virtual tgt::Bounds getBoundingBox(bool transformed = true) const {
        tgt::Bounds bounds;
        for (size_t i=0; i<points_.size(); i++)
            bounds.addPoint(points_[i]);

        if(transformed)
            return bounds.transform(getTransformationMatrix());
        else
            return bounds;
    }

    virtual void serialize(Serializer& s) const {
        s.serialize("points", points_);
        Geometry::serialize(s);
    }

    virtual void deserialize(Deserializer& s) {
        s.deserialize("points", points_);
        Geometry::deserialize(s);
    }

protected:
    std::vector<T> points_;

};

// -------------------------------------------------------------------------

class PointListGeometryVec3 : public PointListGeometry<tgt::vec3> {
public:
    PointListGeometryVec3() : PointListGeometry<tgt::vec3>() {}
    virtual Geometry* create() const { return new PointListGeometryVec3(); }
    virtual std::string getClassName() const { return "PointListGeometryVec3"; }
    virtual std::unique_ptr<Geometry> clone() const {
        PointListGeometryVec3* clone = new PointListGeometryVec3();
        clone->points_ = points_;
        clone->setTransformationMatrix(getTransformationMatrix());
        return std::unique_ptr<Geometry>(clone);
    }
};

//-----------------------------------------------------------------------------
//Interpolation helper functions:

template<typename T>
float getLength(const std::vector<T>& v) {
    float l = 0.0f;
    for(size_t i=1; i<v.size(); i++)
        l += distance(v[i], v[i-1]);
    return l;
}

template<typename T>
T h1(T t) {
    return (2.0f*t*t*t - 3.0f*t*t + 1.0f);
}

template<typename T>
T h2(T t) {
    return (-2.0f*t*t*t + 3.0f*t*t);
}

template<typename T>
T h3(T t) {
    return (t*t*t - 2.0f*t*t + t);
}

template<typename T>
T h4(T t) {
    return (t*t*t - t*t);
}

template<typename T>
T getChord1(const std::vector<T>& v, size_t i) {
    if(i == 0)
        return T(0.0f);
    else
        return v[i] - v[i-1];
}

template<typename T>
T getChord2(const std::vector<T>& v, size_t i) {
    if(i == (v.size()-1))
        return T(0.0f);
    else
        return v[i+1] - v[i];
}

template<typename T>
T getTangent(const std::vector<T>& v, size_t i) {
    return (getChord1(v, i) + getChord2(v, i)) * 0.5f;
}

template<typename T>
T hermite(T p0, T p1, T m0, T m1, float t) {
    return h1(t)*p0 + h3(t)*m0 + h4(t)*m1 + h2(t)*p1;
}

/// Use with PointListGeometry::getData()
template<typename T>
T catmullRom(const std::vector<T>& v, float t, float length) {
    if(v.size() == 0)
        return T(0.0f);
    if(v.size() == 1)
        return v[0];

    float dist = t*length;
    float d = 0.0f;
    for(size_t i=0; i<(v.size()-1); i++) {
        float dTemp = distance(v[i], v[i+1]);
        if((d+dTemp) >= dist) {
            T p0 = v[i];
            T p1 = v[i+1];
            T m0 = getTangent(v, i);
            T m1 = getTangent(v, i+1);

            float t2 = (dist - d) / dTemp;
            return hermite(p0, p1, m0, m1, t2);
        }

        d += dTemp;
    }
    return T(0.0f);
}

} // namespace

#endif // VRN_POINTLISTGEOMETRY_H

