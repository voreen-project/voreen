/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2018 University of Muenster, Germany,           *
 * Department of Computer Science.                                    *
 *                                                                    *
 * This file is part of the tgt library. This library is free         *
 * software; you can redistribute it and/or modify it under the terms *
 * of the GNU Lesser General Public License version 2.1 as published  *
 * by the Free Software Foundation.                                   *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       *
 * GNU Lesser General Public License for more details.                *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License in the file "LICENSE.txt" along with this library.         *
 * If not, see <http://www.gnu.org/licenses/>.                        *
 *                                                                    *
 **********************************************************************/

#ifndef TGT_BOUNDS_H
#define TGT_BOUNDS_H

#include "tgt/types.h"
#include "tgt/vector.h"
#include "tgt/matrix.h"
#include "tgt/plane.h"

namespace tgt {

/**
*   Axis-aligned bounding box
*/
template<typename T>
class TemplateBounds {
    typename tgt::Vector3<T> llf_; //lower left front
    typename tgt::Vector3<T> urb_; //upper right back

public:
    /**
     *   Constructs an undefined boundingbox
     */
    TemplateBounds() : llf_(tgt::Vector3<T>(std::numeric_limits<T>::max())), urb_(tgt::Vector3<T>(std::numeric_limits<T>::lowest())) {}

    /**
     *   Constructs an undefined boundingbox containing v
     */
    TemplateBounds(const typename tgt::Vector3<T>& v)
      : llf_(v),
        urb_(v)
    {}

    /**
     *   Constructs a bounding box containing v1 and v2
     */
    TemplateBounds(const typename tgt::Vector3<T>& v1, const typename tgt::Vector3<T>& v2)
      : llf_(v1),
        urb_(v1) {

        addPoint(v2);
    }

    /**
     *   Enlarges the box (if necessary) to contain v
     */
    void addPoint(const typename tgt::Vector3<T>& v);
    void addPoint(const typename tgt::Vector4<T>& v);

    /**
     *   Enlarges the box (if necessary) to contain b
     */
    void addVolume(const TemplateBounds& b) {
        addPoint(b.llf_);
        addPoint(b.urb_);
    }

    /**
     *   Sets the box to the intersection with b.
     *   This might result in an undefined box,
     *   if there is no intersection.
     */
    void intersectVolume(const TemplateBounds& b) {
        llf_ = max(llf_, b.llf_);
        urb_ = min(urb_, b.urb_);
    }

    /**
     *   Returns the Lower Left Front point
     */
    tgt::Vector3<T> getLLF() const { return llf_; }

    /**
     *   Returns the Upper Right Back point
     */
    tgt::Vector3<T> getURB() const { return urb_; }

    /**
     *   Returns the center of this box
     */
    tgt::Vector3<T> center() const {
        return (diagonal() / T(2) + llf_);
    }

    /**
     *   Returns the diagonal from llf to urb
    */
     tgt::Vector3<T> diagonal() const {
        return (urb_ - llf_);
    }

    /**
     *   Returns the volume of this box
     */
    T volume() const {
        T vol = (urb_.x - llf_.x) * (urb_.y - llf_.y) * (urb_.z - llf_.z);
        tgtAssert(vol >= 0, "Invalid bounds volume");
        return vol;
    }

    /**
     *   Returns true if box is defined:
     *      * not only a point
     *      * positive or zero volume
     */
    bool isDefined() const {
        if (llf_.x > urb_.x) return false;
        if (llf_.y > urb_.y) return false;
        if (llf_.z > urb_.z) return false;
        return !onlyPoint();
    }

    /**
     *   Returns true if box is only a point
     */
    bool onlyPoint() const { return urb_ == llf_; }

    /**
     *   Returns true if point is contained in this box
     */
    bool containsPoint(const typename tgt::Vector3<T>& p) const {
        return ( (p.x >= llf_.x) && (p.y >= llf_.y) && (p.z >= llf_.z)
                    && (p.x <= urb_.x) && (p.y <= urb_.y) && (p.z <= urb_.z) );
    }

    /**
     *   Returns true if b is contained in this box
     *   Box has to be defined!
     */
    bool containsVolume(const TemplateBounds& b) const {
        return ( containsPoint(b.llf_) && containsPoint(b.urb_) );
    }

    /**
     *   Returns true if the boxes intersect
     *   Box has to be defined!
     */
    bool intersects(const TemplateBounds& b) const {
        // Look for a separating axis on each box for each axis
        if ((llf_.x > b.urb_.x) || (b.llf_.x > urb_.x)) return false;
        if ((llf_.y > b.urb_.y) || (b.llf_.y > urb_.y)) return false;
        if ((llf_.z > b.urb_.z) || (b.llf_.z > urb_.z)) return false;

        // No separating axis ... they must intersect
        return true;
   }

    bool intersects(const plane& p) const;

    bool insideXZ(const TemplateBounds& bounds) const;
    bool insideXZ(const typename tgt::Vector3<T>& v) const;

    ///Returns true if bounds is inside this.
    bool inside(const TemplateBounds& bounds) const;
    ///Returns true if v is inside this.
    bool inside(const typename tgt::Vector3<T>& v) const;

    /**
     * Transform this BB using m.
     *
     * @return A BoundingBox containing the 8 vertices defined by this BB, transformed using m.
     */
    TemplateBounds transform(const mat4& m) const;

    bool operator==(const tgt::TemplateBounds<T>& other) const{
        if (llf_ != other.llf_)
            return false;
        if (urb_ != other.urb_)
            return false;
        return true;
    }

    bool operator!=(const tgt::TemplateBounds<T>& other) const{
        return !(*this == other);
    }
};

template<typename T>
void TemplateBounds<T>::addPoint(const typename tgt::Vector3<T>& v) {
    llf_ = min(llf_, v);
    urb_ = max(urb_, v);
}

template<typename T>
void TemplateBounds<T>::addPoint(const typename tgt::Vector4<T>& v) {
    addPoint(v.xyz());
}

template<typename T>
bool TemplateBounds<T>::insideXZ(const TemplateBounds& bounds) const {
    tgtAssert(       isDefined(), "This Box ist not defined.");
    tgtAssert(bounds.isDefined(), "Box b ist not defined.");

    typename tgt::Vector3<T> llfb = bounds.getLLF();
    typename tgt::Vector3<T> urbb = bounds.getURB();

    float r0x = min(llf_[0], urb_[0]);
    float r1x = max(llf_[0], urb_[0]);
    float r0y = min(llf_[2], urb_[2]);
    float r1y = max(llf_[2], urb_[2]);
    float r2x = min(llfb[0], urbb[0]);
    float r3x = max(llfb[0], urbb[0]);
    float r2y = min(llfb[2], urbb[2]);
    float r3y = max(llfb[2], urbb[2]);

    return (r0x >= r2x) && (r0y >= r2y)
        && (r1x <= r3x) && (r1y <= r3y);
}

template<typename T>
bool TemplateBounds<T>::insideXZ(const typename tgt::Vector3<T>& v) const {
    tgtAssert(  isDefined(), "This Box ist not defined.");

    return (llf_[0] <= v[0]) && (v[0] <= urb_[0])
        && (llf_[2] <= v[2]) && (v[2] <= urb_[2]);
}

template<typename T>
bool TemplateBounds<T>::inside(const TemplateBounds& bounds) const {
    tgtAssert(       isDefined(), "This Box ist not defined.");
    tgtAssert(bounds.isDefined(), "Box b ist not defined.");

    typename tgt::Vector3<T> llfb = bounds.getLLF();
    typename tgt::Vector3<T> urbb = bounds.getURB();

    float r0x = min(llf_[0], urb_[0]);
    float r1x = max(llf_[0], urb_[0]);
    float r0y = min(llf_[1], urb_[1]);
    float r1y = max(llf_[1], urb_[1]);
    float r0z = min(llf_[2], urb_[2]);
    float r1z = max(llf_[2], urb_[2]);

    float r2x = min(llfb[0], urbb[0]);
    float r3x = max(llfb[0], urbb[0]);
    float r2y = min(llfb[1], urbb[1]);
    float r3y = max(llfb[1], urbb[1]);
    float r2z = min(llfb[2], urbb[2]);
    float r3z = max(llfb[2], urbb[2]);

    return (r0x >= r2x) && (r1x <= r3x)
        && (r0y >= r2y) && (r1y <= r3y)
        && (r0z >= r2z) && (r1z <= r3z);
}

template<typename T>
bool TemplateBounds<T>::inside(const typename tgt::Vector3<T>& v) const {
    tgtAssert(  isDefined(), "This Box ist not defined.");

    return (llf_[0] <= v[0]) && (v[0] <= urb_[0])
        && (llf_[1] <= v[1]) && (v[1] <= urb_[1])
        && (llf_[2] <= v[2]) && (v[2] <= urb_[2]);
}

template<typename T>
TemplateBounds<T> TemplateBounds<T>::transform(const mat4& m) const {
    TemplateBounds b;
    b.addPoint(m * typename tgt::Vector3<T>(llf_.x, llf_.y, llf_.z));
    b.addPoint(m * typename tgt::Vector3<T>(urb_.x, llf_.y, llf_.z));
    b.addPoint(m * typename tgt::Vector3<T>(llf_.x, urb_.y, llf_.z));
    b.addPoint(m * typename tgt::Vector3<T>(llf_.x, llf_.y, urb_.z));
    b.addPoint(m * typename tgt::Vector3<T>(urb_.x, urb_.y, llf_.z));
    b.addPoint(m * typename tgt::Vector3<T>(llf_.x, urb_.y, urb_.z));
    b.addPoint(m * typename tgt::Vector3<T>(urb_.x, llf_.y, urb_.z));
    b.addPoint(m * typename tgt::Vector3<T>(urb_.x, urb_.y, urb_.z));
    return b;
}

template<typename T>
bool TemplateBounds<T>::intersects(const plane& p) const {
    bool pointsNeg = false;
    bool pointsPos = false;

    float d = p.distance(typename tgt::Vector3<T>(llf_.x, llf_.y, llf_.z));
    if(d < 0.0f) pointsNeg = true; else if(d > 0.0f) pointsPos = true;
    d = p.distance(typename tgt::Vector3<T>(urb_.x, llf_.y, llf_.z));
    if(d < 0.0f) pointsNeg = true; else if(d > 0.0f) pointsPos = true;
    d = p.distance(typename tgt::Vector3<T>(llf_.x, urb_.y, llf_.z));
    if(d < 0.0f) pointsNeg = true; else if(d > 0.0f) pointsPos = true;
    d = p.distance(typename tgt::Vector3<T>(llf_.x, llf_.y, urb_.z));
    if(d < 0.0f) pointsNeg = true; else if(d > 0.0f) pointsPos = true;
    d = p.distance(typename tgt::Vector3<T>(urb_.x, urb_.y, llf_.z));
    if(d < 0.0f) pointsNeg = true; else if(d > 0.0f) pointsPos = true;
    d = p.distance(typename tgt::Vector3<T>(llf_.x, urb_.y, urb_.z));
    if(d < 0.0f) pointsNeg = true; else if(d > 0.0f) pointsPos = true;
    d = p.distance(typename tgt::Vector3<T>(urb_.x, llf_.y, urb_.z));
    if(d < 0.0f) pointsNeg = true; else if(d > 0.0f) pointsPos = true;
    d = p.distance(typename tgt::Vector3<T>(urb_.x, urb_.y, urb_.z));
    if(d < 0.0f) pointsNeg = true; else if(d > 0.0f) pointsPos = true;

    return (pointsNeg && pointsPos);
}

/// ostream-operator
template<typename T>
std::ostream& operator<< (std::ostream& o, const TemplateBounds<T>& b) {
    return (o << "(llf: " << b.getLLF() << " urb: " << b.getURB() << ")");
}

typedef TemplateBounds<float> Bounds;
typedef TemplateBounds<int> IntBounds;
typedef TemplateBounds<size_t> SBounds;

template<typename T>
class TGT_API HasTemplateBounds {
public:
    HasTemplateBounds(const TemplateBounds<T>& bounds)
      : boundingBox_(bounds)
    {}

    HasTemplateBounds()
      : boundingBox_(TemplateBounds<T>())
    {}

    /**
    *   Returns the boundingbox.
    */
    const TemplateBounds<T>& getTemplateBounds() const {
        return boundingBox_;
    }

protected:
    TemplateBounds<T> boundingBox_;
};

typedef HasTemplateBounds<float> HasBounds;


} // namespace

#endif //TGT_BOUNDS_H
