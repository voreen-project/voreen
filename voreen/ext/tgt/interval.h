/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2020 University of Muenster, Germany,           *
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

#ifndef TGT_INTERVAL_H
#define TGT_INTERVAL_H

#include <limits>
namespace tgt{
    template<typename T>
    class Interval{
    public:
        Interval(tgt::Vector2<T> val = tgt::Vector2<T>(T(0), T(1)), T minValue = 0, T maxValue = T(1),
            T minRange=T(0), T maxRange=std::numeric_limits<T>::max());

        bool operator==(const tgt::Interval<T>& other);
        bool operator!=(const tgt::Interval<T>& other);

        void setMinValue(T val);
        void setMaxValue(T val);
        void setMinRange(T val);
        void setMaxRange(T val);

        T getMinValue() const;
        T getMaxValue() const;
        T getMinRange() const;
        T getMaxRange() const;

        tgt::Vector2<T> get() const;
        void set(tgt::Vector2<T> val);
        bool isValid(tgt::Vector2<T> val);

        tgt::Vector2<T> fixupMin(tgt::Vector2<T> val) const;
        tgt::Vector2<T> fixupMax(tgt::Vector2<T> val) const;

        void fixup();
    private:
        tgt::Vector2<T> value_;
        T minValue_;
        T maxValue_;
        T minRange_;
        T maxRange_;
    };

}

namespace tgt{

template<typename T>
Interval<T>::Interval(tgt::Vector2<T> val, T minValue, T maxValue,
            T minRange, T maxRange)
            : value_(val)
            , minValue_(minValue)
            , maxValue_(maxValue)
            , minRange_(minRange)
            , maxRange_(maxRange){
    fixup();
}

template<typename T>
bool Interval<T>::operator==(const tgt::Interval<T>& other){
    if (value_ != other.value_)
        return false;
    if (minValue_ != other.minValue_ || maxValue_ != other.maxValue_)
        return false;
    if (minRange_ != other.minRange_|| maxRange_!= other.maxRange_)
        return false;
    return true;
}

template<typename T>
bool Interval<T>::operator!=(const tgt::Interval<T>& other){
    return !(*this == other);
}

template<typename T>
void Interval<T>::setMinValue(T val){
    minValue_ = val;
    fixup();
}

template<typename T>
void Interval<T>::setMaxValue(T val){
    maxValue_ = val;
    fixup();
}

template<typename T>
void Interval<T>::setMinRange(T val){
    minRange_ = val;
    fixup();
}

template<typename T>
void Interval<T>::setMaxRange(T val){
    maxRange_ = val;
    fixup();
}

template<typename T>
T Interval<T>::getMinValue() const{
    return minValue_;
}

template<typename T>
T Interval<T>::getMaxValue() const{
    return maxValue_;
}

template<typename T>
T Interval<T>::getMinRange() const{
    return minRange_;
}

template<typename T>
T Interval<T>::getMaxRange() const{
    return maxRange_;
}

template<typename T>
tgt::Vector2<T> Interval<T>::get() const{
    return value_;
}

template<typename T>
void Interval<T>::set(tgt::Vector2<T> val){
    value_ = val;
    fixup();
}

template<typename T>
bool Interval<T>::isValid(tgt::Vector2<T> val){
    T min = val.x;
    T max = val.y;
    T range = max-min;

    if (range < minRange_)
        return false;
    if (range > maxRange_)
        return false;
    if (min < minValue_)
        return false;
    if (max > maxValue_)
        return false;
    return true;
}

template<typename T>
tgt::Vector2<T> Interval<T>::fixupMin(tgt::Vector2<T> minmax) const{
    T min = minmax.x;
    T max = minmax.y;

    if (min+minRange_ > max){
        min = max-minRange_;
    }else if (maxRange_ != std::numeric_limits<T>::max() &&  min+maxRange_ < max){
        // HACK(apv):
        // if maxRange_ is the maximum value, we assume
        // the maxRange is infinite to prevent integer overflows
        min = max-maxRange_;
    }

    if (min > maxValue_){
        min = maxValue_-minRange_;
    }
    if (min < minValue_){
        min  = minValue_;
    }
    return typename tgt::Vector2<T> (min, max);
}

template<typename T>
tgt::Vector2<T> Interval<T>::fixupMax(tgt::Vector2<T> minmax) const{
    T min = minmax.x;
    T max = minmax.y;

    if (min+minRange_ > max){
        max = min+minRange_;
    }else if (maxRange_ != std::numeric_limits<T>::max() && min+maxRange_ < max){
        // HACK(apv):
        // if maxRange_ is the maximum value, we assume
        // the maxRange is infinite to prevent integer overflows
        max = min+maxRange_;
    }

    if (max > maxValue_){
        max = maxValue_;
    }
    if (max < minValue_){
        max  = minValue_+minRange_;
    }
    return typename tgt::Vector2<T> (min, max);
}

template<typename T>
void Interval<T>::fixup(){
    tgt::Vector2<T> v = value_;
    v = fixupMin(v);
    v = fixupMax(v);
    value_ = v;
}

}

#endif
