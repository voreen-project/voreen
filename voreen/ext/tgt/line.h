/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2019 University of Muenster, Germany,           *
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

#ifndef TGT_LINE_H
#define TGT_LINE_H

#include "tgt/vector.h"

namespace tgt {

///This implements a 2D line(segment).
template <class T>
class Line2 {
    Vector2<T> p_; /// position
    Vector2<T> d_; /// direction

public:
    /// trivial standard constructor
    Line2() {}

    /// construct a line from starting and end point
    Line2(const Vector2<T>& start, const Vector2<T>& end)
        : p_(start)
        , d_(end - start)
    {}

    /**
    * Interpolates between start and end
    * @param t Should be betweeen 0 and 1 for line segments
    * @return (start*(1-t) + end*t)
    */
    Vector2<T> getFromParam(T t) {
        return p_ + t * d_;
    }

    Vector2<T> getStart() const { return p_; }
    Vector2<T> getEnd() const { return p_ + d_; }
    ///Set startingpoint to start, leaves end unchanged.
    void setStart(const Vector2<T>& start) { d_ = (p_ + d_) - start; p_ = start; }
    ///Set endpoint to end, leaves start unchanged.
    void setEnd(const Vector2<T>& end) { d_ = (end - p_); }

    ///@return Vector from start to end
    Vector2<T> getDirection() const { return d_; }
    ///Set direction, changing endpoint, leaves start unchanged.
    void setDirection(const Vector2<T>& dir) { d_ = dir; }

    /**
    * intersect a linesegment with another linesegment
    * See RTR p617.
    *
    * @param l the line to intersect this object with
    * @param v holds the result if an intersection took place,
    * otherwise it remains unchanged
    * @return true if these lines intersect, false otherwise
    */
    bool intersectSegments(const Line2<T>& l, Vector2<T>& v) const;

    /**
    * intersect a line with another line
    *
    * @param l the line to intersect this object with
    * @param v holds the result if an intersection took place,
    * otherwise it remains unchanged
    * @return true if these lines intersect, false otherwise
    */
    bool intersect(const Line2<T>& l, Vector2<T>& v) const;
};

// non inline implementation

template <class T>
bool Line2<T>::intersectSegments(const Line2<T>& l, Vector2<T>& v) const {
    T f = dot(l.d_, Vector2<T>(-d_.y, d_.x));

    T e = dot((p_-l.p_), Vector2<T>(-d_.y, d_.x));
    T d = dot((p_-l.p_), Vector2<T>(-l.d_.y, l.d_.x));

//     std::cout << "d: " << d << " e: "<< e << " f: " << f << "\n";

    //test if s and t are in range [0, 1]:
    if(f > 0.0) {
        if((d < 0.0) || (d > f))
            return false;
    }
    else {
        if((d > 0.0) || (d < f))
            return false;
    }

    if(f > 0.0) {
        if((e < 0.0) || (e > f))
            return false;
    }
    else {
        if((e > 0.0) || (e < f))
            return false;
    }

//     T t = e/f;
    T s = d/f;
//     std::cout << "s: " << s << " t: "<< t << "\n";

//     v = l.p_ + t * l.d_;
//     std::cout << "v: " << v << "\n";
    v = p_ + (s * d_);
//     std::cout << "v: " << v << "\n";
    return true;
}


template <class T>
bool Line2<T>::intersect(const Line2<T>& l, Vector2<T>& v) const {
    // check whether these lines interesect
    T q = d_.y * l.d_.x - d_.x * l.d_.y;

    if (q == 0.f)
        return false; // these lines a collinear or stack

    // calculate paramter for this line where the lines intersect
    T t = -( l.d_.x * (p_.y - l.p_.y) + l.d_.y * l.p_.x - l.d_.y * p_.x ) / q;
    v = p_ + t * d_;

    return true;
}
//-------------------------------------------------------------------------------------------------

///This implements a 3D line(segment).
template <class T>
class Line3 {
    Vector3<T> p_; /// position
    Vector3<T> d_; /// direction

public:
    /// trivial standard constructor
    Line3() {}

    /// construct a line from two points on a line
    Line3(const Vector3<T>& start, const Vector3<T>& end)
        : p_(start)
        , d_(end - start)
    {}

    /**
    * Interpolates between start and end
    * @param t Should be betweeen 0 and 1.
    * @return (start*(1-t) + end*t)
    */
    Vector3<T> getFromParam(T t) const {
        return p_ + t * d_;
    }

    Vector3<T> getStart() const { return p_; }
    Vector3<T> getEnd() const { return p_ + d_; }
    ///Set startingpoint to start, leaves end unchanged.
    void setStart(const Vector3<T>& start) { d_ = (p_ + d_) - start; p_ = start; }
    ///Set endpoint to end, leaves start unchanged.
    void setEnd(const Vector3<T>& end) { d_ = (end - p_); }

    ///@return Vector from start to end
    Vector3<T> getDirection() const { return d_; }
    ///Set direction, changing endpoint, leaves start unchanged.
    void setDirection(const Vector3<T>& dir) { d_ = dir; }

    /**
     * intersect a linesegment with another linesegment
     * See RTR p618
     *
     * @param l the line to intersect this object with
     * @param v holds the result if an intersection took place,
     * otherwise it remains unchanged
     * @return true if these lines intersect
    */
    bool intersectSegments(const Line3<T>& l, Vector3<T>& v) const;

    /**
     * intersect a line with another line
     * See RTR p618
     *
     * @param l the line to intersect this object with
     * @param v holds the result if an intersection took place,
     * otherwise it remains unchanged
     * @return true if these lines intersect
    */
    bool intersect(const Line3<T>& l, Vector3<T>& v) const;
};

//------------------------------------------------------------------------------
//  non inline implementation
//------------------------------------------------------------------------------

template <class T>
bool Line3<T>::intersectSegments(const Line3<T>& l, Vector3<T>& v) const {
    Vector3<T> c = cross(d_, l.d_);
    T q = lengthSq(c);

    if (std::abs(q) < 0.00001f)
        return false;
    else {
        Vector3<T> o = l.p_ - p_;
        Vector3<T> d = l.d_;

        T t = o.x*d_.y*c.z + d_.x*c.y*o.z + c.x*o.y*d_.z
            - c.x*d_.y*o.z - d_.x*o.y*c.z - o.x*c.y*d_.z;

        T s = o.x*l.d_.y*c.z + l.d_.x*c.y*o.z + c.x*o.y*l.d_.z
            - c.x*l.d_.y*o.z - l.d_.x*o.y*c.z - o.x*c.y*l.d_.z;

        t /= q;
        s /= q;

        v = p_ + t*d_;
        if( (s >= 0.0) && (s <= 1.0)
          &&(t >= 0.0) && (t <= 1.0) )
            return true;
        else
            return false;
    }
}

template <class T>
bool Line3<T>::intersect(const Line3<T>& l, Vector3<T>& v) const {
    Vector3<T> c = cross(d_, l.d_);
    T q = lengthSq(c);

    if (std::abs(q) < 0.00001f)
        return false;
    else {
        Vector3<T> o = l.p_ - p_;
        Vector3<T> d = l.d_;

        T t = o.x*d_.y*c.z + d_.x*c.y*o.z + c.x*o.y*d_.z
            - c.x*d_.y*o.z - d_.x*o.y*c.z - o.x*c.y*d_.z;

        t /= q;
        v = p_ + t*d_;

        return true;
    }
}

/// ostream-operator
template <class T>
std::ostream& operator << (std::ostream& s, const Line2<T>& l) {
    return (s << "x = " << l.getStart() << " + t * " << l.getDirection());
}

template <class T>
std::ostream& operator << (std::ostream& s, const Line3<T>& l) {
    return (s << "x = " << l.getStart() << " + t * " << l.getDirection());
}

// typedefs for easy usage
// defined types with a small 'l' like other absolute math classes in tgt
typedef Line2<float>   line2;
typedef Line3<float>   line3;
typedef Line2<double> dline2;
typedef Line3<double> dline3;

} // namespace tgt

#endif //TGT_LINE_H
