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

#ifndef VRN_INTERVAL_H
#define VRN_INTERVAL_H

#include "voreen/core/io/serialization/serialization.h"
#include <limits>
#include "tgt/tgt_math.h"

namespace voreen {

//namespace {
//
//template<typename T>
//bool isNan(const T& value) {
//    return (value != value);
//};
//
//}

class VoreenException;
/**
 * \brief templated class for handling intervals of numeric types
 *
 * Offers support for common interval arithmetic.
 **/
template<typename T>
class Interval : public Serializable  {
public:
    /// default constructor (creates the empty interval (0,0))
    Interval();

    /// copy constructor
    Interval(const Interval<T>& src);

    /// constructor
    Interval(T left, T right, bool leftOpen=false, bool rightOpen=false);

    /// destructor
    virtual ~Interval();

    /// returns left interval border
    inline const T& getLeft() const;
    /// returns right interval border
    inline const T& getRight() const;

    /// returns whether left interval border ist open
    inline bool getLeftOpen() const;
    /// returns whether right interval border ist open
    inline bool getRightOpen() const;

    /// returns whether interval is empty
    bool empty() const;

    /// returns the size of the interval (right - left)
    T size() const;

    /// returns whether the interval contains the value \a value
    bool contains(const T& value) const;
    /// returns whether the interval contains the interval \a interval
    bool contains(const Interval<T>& interval) const;


    /// intersects this interval with \a rhs
    void intersectWith(const Interval<T>& rhs);

    /// builds the smallest interval containing both intervals *this and \a rhs
    /// if one of these intervals is empty, the other is returned
    void unionWith(const Interval<T>& rhs);

    /// sets all infinite, NaN or std::numeric_limits<T>::max() values to according values in rhs
    void clampInfinitesTo(const Interval<T>& rhs);


    /// enlarges the current interval by factor \a factor
    void enlarge(const T& factor);

    /// enlarges current interval as much as it contains \a value
    void nibble(const T& value);

    /// returns \a value clamped to current interval
    T clampValue(const T& value) const;

    std::string toString() const;


    /**
    * @see Property::serialize
    */
    virtual void serialize(Serializer& s) const;

    /**
    * @see Property::deserialize
    */
    virtual void deserialize(Deserializer& s);


    /// returns an interval whose borders are incremented by \a rhs
    Interval<T> operator+(const T& rhs) const;
    /// returns an interval whose borders are decremented by \a rhs
    Interval<T> operator-(const T& rhs) const;
    /// returns an interval whose borders are multiplied by \a rhs
    Interval<T> operator*(const T& rhs) const;
    /// returns an interval whose borders are divided by \a rhs
    Interval<T> operator/(const T& rhs) const;


    bool operator==(const Interval<T>& rhs) const;
    bool operator!=(const Interval<T>& rhs) const;

private:
    T left_;                ///< left interval border
    T right_;               ///< right interval border

    bool leftOpen_;         ///< flag whether left interval border is open
    bool rightOpen_;        ///< flag whether right interval border is open

};



// ---------------------------------------- implementation ----------------------------------------

template<typename T>
Interval<T>::Interval()
    : Serializable()
    , left_(0)
    , right_(0)
    , leftOpen_(true)
    , rightOpen_(true)
{}

template<typename T>
Interval<T>::Interval(const Interval<T>& src)
    : Serializable()
    , left_(src.left_)
    , right_(src.right_)
    , leftOpen_(src.leftOpen_)
    , rightOpen_(src.rightOpen_)
{}


template<typename T>
Interval<T>::Interval(T left, T right, bool leftOpen, bool rightOpen)
    : left_(left)
    , right_(right)
    , leftOpen_(leftOpen)
    , rightOpen_(rightOpen)
{}

template<typename T>
Interval<T>::~Interval() {
}

template<typename T>
const T& Interval<T>::getLeft() const {
    return left_;
}

template<typename T>
const T& Interval<T>::getRight() const {
    return right_;
}

template<typename T>
bool Interval<T>::getLeftOpen() const {
    return leftOpen_;
}

template<typename T>
bool Interval<T>::getRightOpen() const {
    return rightOpen_;
}

template<typename T>
T Interval<T>::size() const {
    return right_ - left_;
}

template<typename T>
bool Interval<T>::empty() const {
    if (left_ > right_)
        return true;
    if (left_ == right_)
        return (leftOpen_ || rightOpen_);
    return false;
}

template<typename T>
bool Interval<T>::contains(const T& value) const {
    if (tgt::isNaN(value))
        return false;
    if ((value < left_) || (leftOpen_ && value == left_))
        return false;
    if ((value > right_) || (rightOpen_ && value == right_))
        return false;
    return true;
}

template<typename T>
bool Interval<T>::contains(const Interval<T>& rhs) const {
    if ((left_ > rhs.left_) || ((leftOpen_ || rhs.leftOpen_) && left_ == rhs.left_))
        return false;
    if ((right_ < rhs.right_) || ((rightOpen_ || rhs.rightOpen_) && right_ == rhs.right_))
        return false;
    return true;
}

template<typename T>
void Interval<T>::intersectWith(const Interval<T>& rhs) {
    if (rhs.empty())
        return;
    if (empty()) {
        left_ = rhs.left_;
        leftOpen_ = rhs.leftOpen_;
        right_ = rhs.right_;
        rightOpen_ = rhs.rightOpen_;
        return;
    }

    if (left_ < rhs.left_) {
        left_ = rhs.left_;
        leftOpen_ = rhs.leftOpen_;
    }
    else if (left_ == rhs.left_) {
        leftOpen_ = leftOpen_ || rhs.leftOpen_;
    }
    else {
        // do nothing, keep values
    }

    if (right_ > rhs.right_) {
        right_ = rhs.right_;
        rightOpen_ = rhs.rightOpen_;
    }
    else if (right_ == rhs.right_) {
        rightOpen_ = rightOpen_ || rhs.rightOpen_;
    }
    else {
        // do nothing, keep values
    }
}

template<typename T>
void Interval<T>::unionWith(const Interval<T>& rhs) {
    // check if one of the intervals is empty
    if (empty()) {
        left_ = rhs.left_;
        leftOpen_ = rhs.leftOpen_;
        right_ = rhs.right_;
        rightOpen_ = rhs.rightOpen_;
        return;
    }
    if (rhs.empty())
        return;

    if (left_ < rhs.left_) {
        // do nothing, keep values
    }
    else if (left_ == rhs.left_) {
        leftOpen_ = leftOpen_ && rhs.leftOpen_;
    }
    else {
        left_ = rhs.left_;
        leftOpen_ = rhs.leftOpen_;
    }

    if (right_ > rhs.right_) {
        // do nothing, keep values
    }
    else if (right_ == rhs.right_) {
        rightOpen_ = rightOpen_ && rhs.rightOpen_;
    }
    else {
        right_ = rhs.right_;
        rightOpen_ = rhs.rightOpen_;
    }
}

template<typename T>
void Interval<T>::clampInfinitesTo(const Interval<T>& rhs) {
    if (rhs.empty())
        return;
    if (empty()) {
        left_ = rhs.left_;
        leftOpen_ = rhs.leftOpen_;
        right_ = rhs.right_;
        rightOpen_ = rhs.rightOpen_;
        return;
    }

    if (tgt::isNaN(left_) || left_ == std::numeric_limits<T>::infinity() || left_ == -std::numeric_limits<T>::max()) {
        left_ = rhs.left_;
        leftOpen_ = rhs.leftOpen_;
    }
    if (tgt::isNaN(right_) || right_ == std::numeric_limits<T>::infinity() || right_ == std::numeric_limits<T>::max()) {
        right_ = rhs.right_;
        rightOpen_ = rhs.rightOpen_;
    }
}

template<typename T>
void Interval<T>::serialize(Serializer& s) const {
    s.serialize("left", left_);
    s.serialize("right", right_);
    s.serialize("leftOpen", leftOpen_);
    s.serialize("rightOpen", rightOpen_);
}

template<typename T>
void Interval<T>::deserialize(Deserializer& s) {
    s.deserialize("left", left_);
    s.deserialize("right", right_);
    s.deserialize("leftOpen", leftOpen_);
    s.deserialize("rightOpen", rightOpen_);
}

template<typename T>
Interval<T> Interval<T>::operator+(const T& rhs) const {
    return Interval<T>(left_+rhs, right_+rhs, leftOpen_, rightOpen_);
}

template<typename T>
Interval<T> Interval<T>::operator-(const T& rhs) const {
    return Interval<T>(left_-rhs, right_-rhs, leftOpen_, rightOpen_);
}

template<typename T>
Interval<T> Interval<T>::operator*(const T& rhs) const {
    return Interval<T>(left_*rhs, right_*rhs, leftOpen_, rightOpen_);
}

template<typename T>
Interval<T> Interval<T>::operator/(const T& rhs) const {
    if (rhs == 0)
        throw VoreenException("division with zero");

    return Interval<T>(left_/rhs, right_/rhs, leftOpen_, rightOpen_);
}

template<typename T>
void Interval<T>::enlarge(const T& factor) {
    T length = right_-left_;
    T offset = length*(factor-1.0)/2.0;
    left_ -= offset;
    right_ += offset;
}

template<typename T>
void Interval<T>::nibble(const T& value) {
    if (empty()) {
        left_ = value;
        right_ = value;
        leftOpen_ = false;
        rightOpen_ = false;
    }
    else if (leftOpen_ && left_ == value)
        leftOpen_ = false;
    else if (rightOpen_ && right_ == value)
        rightOpen_ = false;
    else if (left_ > value) {
        left_ = value;
        leftOpen_ = false;
    }
    else if (right_ < value) {
        right_ = value;
        rightOpen_ = false;
    }
}

template<typename T>
T Interval<T>::clampValue(const T& value) const {
    if (value < left_)
        return left_;
    else if (value > right_)
        return right_;
    else // contains(value)
        return value;
}

template<typename T>
std::string Interval<T>::toString() const {
    std::stringstream result;
    result << "Left: " << left_ << ", Right: " << right_ << ", LeftOpen: " << leftOpen_<< ", RightOpen: " << rightOpen_;
    return result.str();
}

template<typename T>
bool Interval<T>::operator==(const Interval<T>& rhs) const {
    return (left_ == rhs.left_
         && right_ == rhs.right_
         && leftOpen_ == rhs.leftOpen_
         && rightOpen_ == rhs.rightOpen_);
}

template<typename T>
bool Interval<T>::operator!=(const Interval<T>& rhs) const {
    return (left_ != rhs.left_
         || right_ != rhs.right_
         || leftOpen_ != rhs.leftOpen_
         || rightOpen_ != rhs.rightOpen_);
}

}// namespace voreen


#endif // VRN_INTERVAL_H
