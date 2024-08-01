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

#ifndef VRN_VOLUMEELEMENT
#define VRN_VOLUMEELEMENT

#include <limits>

#include "voreen/core/datastructures/tensor.h"
#include "voreen/core/voreencoreapi.h"

// Note: please ensure that no OpenGL dependencies are added into this file

namespace voreen {


namespace {
/**
 * Used internally by VolumeElement.
 */
template<class T>
T getTypeUpperLimit() {
    return std::numeric_limits<T>::max();
}

template<class T>
T getTypeLowerLimit() {
    return std::numeric_limits<T>::min();
}

template<class T>
bool isTypeSigned() {
    return std::numeric_limits<T>::is_signed;
}

template<class T>
bool isTypeInteger() {
    return std::numeric_limits<T>::is_integer;
}

template<>
float getTypeUpperLimit() {
    return 1.f;
}

template<>
float getTypeLowerLimit() {
    return 0.f;
}

template<>
double getTypeUpperLimit() {
    return 1.0;
}

template<>
double getTypeLowerLimit() {
    return 0.0;
}

uint8_t swapEndian(uint8_t value) {
    return value;
}
int8_t swapEndian(int8_t value) {
    return value;
}

uint16_t swapEndian(uint16_t value) {
    return (value>>8) | (value<<8);
}
int16_t swapEndian(int16_t value) {
    return (value>>8) | (value<<8);
}

uint32_t swapEndian(uint32_t value) {
    return (value>>24) |
        ((value<<8) & 0x00FF0000) |
        ((value>>8) & 0x0000FF00) |
        (value<<24);
}

int32_t swapEndian(int32_t value) {
    return (value>>24) |
        ((value<<8) & 0x00FF0000) |
        ((value>>8) & 0x0000FF00) |
        (value<<24);
}

uint64_t swapEndian(uint64_t value) {
    return (value>>56) |
        ((value<<40) & 0x00FF000000000000) |
        ((value<<24) & 0x0000FF0000000000) |
        ((value<<8)  & 0x000000FF00000000) |
        ((value>>8)  & 0x00000000FF000000) |
        ((value>>24) & 0x0000000000FF0000) |
        ((value>>40) & 0x000000000000FF00) |
        (value<<56);
}

int64_t swapEndian(int64_t value) {
    return (value>>56) |
        ((value<<40) & 0x00FF000000000000) |
        ((value<<24) & 0x0000FF0000000000) |
        ((value<<8)  & 0x000000FF00000000) |
        ((value>>8)  & 0x00000000FF000000) |
        ((value>>24) & 0x0000000000FF0000) |
        ((value>>40) & 0x000000000000FF00) |
        (value<<56);
}

float swapEndian(float value) {
    float result;
    char* src = (char*)&value;
    char* dst = (char*)&result;
    dst[0] = src[3];
    dst[1] = src[2];
    dst[2] = src[1];
    dst[3] = src[0];
    return result;
}

double swapEndian(double value) {
    double result;
    char* src = (char*)&value;
    char* dst = (char*)&result;
    dst[0] = src[7];
    dst[1] = src[6];
    dst[2] = src[5];
    dst[3] = src[4];
    dst[4] = src[3];
    dst[5] = src[2];
    dst[6] = src[1];
    dst[7] = src[0];
    return result;
}

}

/**
 * Helper class for getting sth like the maximum value
 */
template<class T>
class VolumeElement {
public:
    inline static T getZero();

    inline static T getChannel(const T& t, size_t channel);
    inline static void setChannel(T value, T& t, size_t channel);

    /**
     * Returns whether the data type is signed.
     */
    inline bool static isSigned();

    /**
     * Returns whether the data type is bounded,
     * which is true for integer types and false for
     * floating point types.
     */
    inline bool static isInteger();

    /**
     * Returns the maximal possible value of T
     */
    inline static T rangeMax();

    /**
     * Returns the minimal possible value of T
     */
    inline static T rangeMin();

    /**
    * Returns the maximum value an element can have. In contrast to max() no
    * vector is returned but the maximum value an element in the vector can have.
    * Since we cannot make any assumption about the element's data type,
    * we have to use float as fall back.
    */
    inline static float rangeMaxElement();

    /**
     * Returns the minimum value an element can have.
     */
    inline static float rangeMinElement();

    /**
     * Extracts the minimal element from the passed value
     * and converts it to float.
     */
    inline static float getMinElement(const T& value);

    /**
     * Extracts the maximal element from the passed value
     * and converts it to float.
     */
    inline static float getMaxElement(const T& value);

    inline static int getNumChannels();

    inline static T swapEndianness(T v);

    /**
    * Calculates the squared difference between val1 and val2.
    */
    static double calcSquaredDifference(T val1, T val2);

    typedef double DoubleType;
    typedef T BaseType;
};

template<class T>
class VolumeElement<tgt::Vector2<T> > {
public:
    inline static tgt::Vector2<T> getZero();

    inline static T getChannel(const tgt::Vector2<T>& t, size_t channel);
    inline static void setChannel(T value, tgt::Vector2<T>& t, size_t channel);

    inline static bool isSigned();

    inline bool static isInteger();

    inline static tgt::Vector2<T> rangeMax();
    inline static tgt::Vector2<T> rangeMin();

    /**
    * Returns the maximum value an element can have. In contrast to max() no
    * vector is returned but the maximum value an element in the vector can have.
    */
    inline static float rangeMaxElement();
    inline static float rangeMinElement();

    inline static float getMaxElement(const tgt::Vector2<T>& value);
    inline static float getMinElement(const tgt::Vector2<T>& value);

    inline static int getNumChannels();

    inline static tgt::Vector2<T> swapEndianness(tgt::Vector2<T> v);

    /**
    * Calculates the squared difference between val1 and val2. Since val1 and val2 are 2-dimensional
    * vectors, each squared difference is multiplied by 1/2.
    */
    static double calcSquaredDifference(tgt::Vector2<T> val1, tgt::Vector2<T> val2);

    typedef tgt::Vector2<double> DoubleType;
    typedef T BaseType;
};

template<class T>
class VolumeElement<tgt::Vector3<T> > {
public:
    inline static tgt::Vector3<T> getZero();

    inline static T getChannel(const tgt::Vector3<T>& t, size_t channel);
    inline static void setChannel(T value, tgt::Vector3<T>& t, size_t channel);

    inline static bool isSigned();

    inline bool static isInteger();

    inline static tgt::Vector3<T> rangeMax();
    inline static tgt::Vector3<T> rangeMin();

    /**
    * Returns the maximum value an element can have. In contrast to max() no
    * vector is returned but the maximum value an element in the vector can have.
    */
    inline static float rangeMaxElement();
    inline static float rangeMinElement();

    inline static float getMaxElement(const tgt::Vector3<T>& value);
    inline static float getMinElement(const tgt::Vector3<T>& value);

    inline static int getNumChannels();

    inline static tgt::Vector3<T> swapEndianness(tgt::Vector3<T> v);

    /**
    * Calculates the squared difference between val1 and val2. Since val1 and val2 are 3-dimensional
    * vectors, each squared difference is multiplied by 1/3.
    */
    static double calcSquaredDifference(tgt::Vector3<T> val1, tgt::Vector3<T> val2);

    typedef tgt::Vector3<double> DoubleType;
    typedef T BaseType;
};

template<class T>
class VolumeElement<tgt::Vector4<T> > {
public:
    inline static tgt::Vector4<T> getZero();

    inline static T getChannel(const tgt::Vector4<T>& t, size_t channel);

    inline static void setChannel(T value, tgt::Vector4<T>& t, size_t channel);

    inline static bool isSigned();

    inline bool static isInteger();

    inline static tgt::Vector4<T> rangeMax();
    inline static tgt::Vector4<T> rangeMin();

    /**
    * Returns the maximum value an element can have. In contrast to max() no
    * vector is returned but the maximum value an element in the vector can have.
    */
    inline static float rangeMaxElement();
    inline static float rangeMinElement();

    inline static float getMaxElement(const tgt::Vector4<T>& value);
    inline static float getMinElement(const tgt::Vector4<T>& value);

    inline static int getNumChannels();

    inline static tgt::Vector4<T> swapEndianness(tgt::Vector4<T> v);

    /**
    * Calculates the squared difference between val1 and val2. Since val1 and val2 are 4-dimensional
    * vectors, each squared difference is multiplied by 1/4.
    */
    static double calcSquaredDifference(tgt::Vector4<T> val1, tgt::Vector4<T> val2);

    typedef tgt::Vector4<double> DoubleType;
    typedef T BaseType;
};

template<class T>
class VolumeElement<tgt::Matrix3<T> > {
public:
    inline static tgt::Matrix3<T> getZero();

    inline static T getChannel(const tgt::Matrix3<T>& t, size_t channel);
    inline static void setChannel(T value, tgt::Matrix3<T>& t, size_t channel);

    inline static bool isSigned();

    inline bool static isInteger();

    inline static tgt::Matrix3<T> rangeMax();
    inline static tgt::Matrix3<T> rangeMin();

    /**
    * Returns the maximum value an element can have. In contrast to max() no
    * vector is returned but the maximum value an element in the vector can have.
    */
    inline static float rangeMaxElement();
    inline static float rangeMinElement();

    inline static float getMaxElement(const tgt::Matrix3<T>& value);
    inline static float getMinElement(const tgt::Matrix3<T>& value);

    inline static int getNumChannels();

    inline static tgt::Matrix3<T> swapEndianness(tgt::Matrix3<T> v);

    /**
    * Calculates the squared difference between val1 and val2. Since val1 and val2 are 3-dimensional
    * vectors, each squared difference is multiplied by 1/3.
    */
    static double calcSquaredDifference(tgt::Matrix3<T> val1, tgt::Matrix3<T> val2);

    typedef tgt::Matrix3<double> DoubleType;
    typedef T BaseType;
};

template<class T>
class VolumeElement<tgt::Matrix4<T> > {
public:
    inline static tgt::Matrix4<T> getZero();

    inline static T getChannel(const tgt::Matrix4<T>& t, size_t channel);
    inline static void setChannel(T value, tgt::Matrix4<T>& t, size_t channel);

    inline static bool isSigned();

    inline bool static isInteger();

    inline static tgt::Matrix4<T> rangeMax();
    inline static tgt::Matrix4<T> rangeMin();

    /**
    * Returns the maximum value an element can have. In contrast to max() no
    * vector is returned but the maximum value an element in the vector can have.
    */
    inline static float rangeMaxElement();
    inline static float rangeMinElement();

    inline static float getMaxElement(const tgt::Matrix4<T>& value);
    inline static float getMinElement(const tgt::Matrix4<T>& value);

    inline static int getNumChannels();

    inline static tgt::Matrix4<T> swapEndianness(tgt::Matrix4<T> v);

    /**
    * Calculates the squared difference between val1 and val2. Since val1 and val2 are 3-dimensional
    * vectors, each squared difference is multiplied by 1/3.
    */
    static double calcSquaredDifference(tgt::Matrix4<T> val1, tgt::Matrix4<T> val2);

    typedef tgt::Matrix4<double> DoubleType;
    typedef T BaseType;
};

template<class T>
class VolumeElement<Tensor2<T> > {
public:
    inline static Tensor2<T> getZero();

    inline static T getChannel(const Tensor2<T>& t, size_t channel);
    inline static void setChannel(T value, Tensor2<T>& t, size_t channel);

    inline static bool isSigned();

    inline bool static isInteger();

    inline static Tensor2<T> rangeMax();
    inline static Tensor2<T> rangeMin();

    /**
    * Returns the maximum value an element can have. In contrast to max() no
    * vector is returned but the maximum value an element in the vector can have.
    */
    inline static float rangeMaxElement();
    inline static float rangeMinElement();

    inline static float getMaxElement(const Tensor2<T>& value);
    inline static float getMinElement(const Tensor2<T>& value);

    inline static int getNumChannels();

    inline static Tensor2<T> swapEndianness(Tensor2<T> v);

    /**
    * Calculates the squared difference between val1 and val2. Since val1 and val2 are 3-dimensional
    * vectors, each squared difference is multiplied by 1/3.
    */
    static double calcSquaredDifference(Tensor2<T> val1, Tensor2<T> val2);

    typedef Tensor2<double> DoubleType;
    typedef T BaseType;
};

template<class T>
float getTypeAsFloat(T /*value*/) {
    tgtAssert(false, "Missing template specialization! Should not get here!");
    //If you got this assertion add a template specialization below.
    return 0.0f;
}

template<>
VRN_CORE_API inline float getTypeAsFloat(float value);
template<>
VRN_CORE_API inline float getTypeAsFloat(double value);

template<>
VRN_CORE_API inline float getTypeAsFloat(uint8_t value);
template<>
VRN_CORE_API inline float getTypeAsFloat(int8_t value);

template<>
VRN_CORE_API inline float getTypeAsFloat(uint16_t value);
template<>
VRN_CORE_API inline float getTypeAsFloat(int16_t value);

template<>
VRN_CORE_API inline float getTypeAsFloat(uint32_t value);
template<>
VRN_CORE_API inline float getTypeAsFloat(int32_t value);

template<>
VRN_CORE_API inline float getTypeAsFloat(uint64_t value);
template<>
VRN_CORE_API inline float getTypeAsFloat(int64_t value);

template<class T>
T getFloatAsType(float /*value*/) {
    tgtAssert(false, "Missing template specialization! Should not get here!");
    //If you got this assertion add a template specialization below.
    return T(0);
}


template<>
VRN_CORE_API inline float getFloatAsType(float value);
template<>
VRN_CORE_API inline double getFloatAsType(float value);

template<>
VRN_CORE_API inline uint8_t getFloatAsType(float value);
template<>
VRN_CORE_API inline int8_t getFloatAsType(float value);

template<>
VRN_CORE_API inline uint16_t getFloatAsType(float value);
template<>
VRN_CORE_API inline int16_t getFloatAsType(float value);

template<>
VRN_CORE_API inline uint32_t getFloatAsType(float value);
template<>
VRN_CORE_API inline int32_t getFloatAsType(float value);

template<>
VRN_CORE_API inline uint64_t getFloatAsType(float value);
template<>
VRN_CORE_API inline int64_t getFloatAsType(float value);



// ==============
// implementation
// ==============

template<class T>
inline T VolumeElement<T>::getZero() {
    return T(0);
}

template<class T>
#ifdef VRN_DEBUG
inline T VolumeElement<T>::getChannel(const T& t, size_t channel) {
    tgtAssert(channel == 0, "channel must be zero");
#else
inline  T VolumeElement<T>::getChannel(const T& t, size_t /*channel*/) {
#endif
    return t;
}

template<class T>
#ifdef VRN_DEBUG
inline void VolumeElement<T>::setChannel(T value, T& t, size_t channel) {
    tgtAssert(channel == 0, "channel must be zero");
#else
inline void VolumeElement<T>::setChannel(T value, T& t, size_t /*channel*/) {
#endif
    t = value;
}

template<class T>
inline bool VolumeElement<T>::isSigned() {
    return isTypeSigned<T>();
}

template<class T>
inline bool VolumeElement<T>::isInteger() {
    return isTypeInteger<T>();
}

template<class T>
inline T VolumeElement<T>::rangeMax() {
    return static_cast<T>( getTypeUpperLimit<T>() );
}

template<class T>
inline T VolumeElement<T>::rangeMin() {
    return static_cast<T>( getTypeLowerLimit<T>() );
}

/**
* Returns the maximum value an element can have. In contrast to max() no
* vector is returned but the maximum value an element in the vector can have.
*/
template<class T>
inline float VolumeElement<T>::rangeMaxElement() {
    return static_cast<float>( getTypeUpperLimit<T>() );
}

template<class T>
float voreen::VolumeElement<T>::getMaxElement(const T& value) {
    return static_cast<float>(value);
}

template<class T>
float voreen::VolumeElement<T>::getMinElement(const T& value) {
    return static_cast<float>(value);
}

template<class T>
inline float VolumeElement<T>::rangeMinElement() {
    return static_cast<float>( getTypeLowerLimit<T>() );
}

template<class T>
inline int VolumeElement<T>::getNumChannels() {
    return 1;
}

template<class T>
inline T VolumeElement<T>::swapEndianness(T v) {
    return swapEndian(v);
}

/**
* Calculates the squared difference between val1 and val2.
*/
template<class T>
double VolumeElement<T>::calcSquaredDifference(T val1, T val2) {
    double result;

    result = static_cast<double>(val1) - static_cast<double>(val2);
    result = result * result;

    return result;
}

// --- Vector2 ----

template<class T>
inline tgt::Vector2<T> VolumeElement<tgt::Vector2<T> >::getZero() {
    return tgt::Vector2<T>(T(0), T(0));
}

template<class T>
inline T VolumeElement<tgt::Vector2<T> >::getChannel(const tgt::Vector2<T>& t, size_t channel) {
    tgtAssert(channel < 2, "channel must be less than 2");
    if (channel < 2)
        return t[channel];
    else
        return 0;
}

template<class T>
inline void VolumeElement<tgt::Vector2<T> >::setChannel(T value, tgt::Vector2<T>& t, size_t channel) {
    tgtAssert(channel < 2, "channel must be less than 2");
    t[channel] = value;
}

template<class T>
inline bool VolumeElement<tgt::Vector2<T> >::isSigned() {
    return isTypeSigned<T>();
}

template<class T>
inline bool VolumeElement<tgt::Vector2<T> >::isInteger() {
    return isTypeInteger<T>();
}

template<class T>
inline tgt::Vector2<T> VolumeElement<tgt::Vector2<T> >::rangeMax() {
    return tgt::Vector2<T>(getTypeUpperLimit<T>());
}

template<class T>
inline tgt::Vector2<T> VolumeElement<tgt::Vector2<T> >::rangeMin() {
    return tgt::Vector2<T>(getTypeLowerLimit<T>());
}

/**
* Returns the maximum value an element can have. In contrast to max() no
* vector is returned but the maximum value an element in the vector can have.
*/
template<class T>
inline float VolumeElement<tgt::Vector2<T> >::rangeMaxElement() {
    return static_cast<float>(getTypeUpperLimit<T>());
}

template<class T>
inline float VolumeElement<tgt::Vector2<T> >::rangeMinElement() {
    return static_cast<float>(getTypeLowerLimit<T>());
}

template<class T>
float voreen::VolumeElement<tgt::Vector2<T> >::getMaxElement(const tgt::Vector2<T>& value) {
    return static_cast<float>(value.x);
}

template<class T>
float voreen::VolumeElement<tgt::Vector2<T> >::getMinElement(const tgt::Vector2<T>& value) {
    return static_cast<float>(value.x);
}

template<class T>
inline int VolumeElement<tgt::Vector2<T> >::getNumChannels() {
    return 2;
}

template<class T>
inline tgt::Vector2<T> VolumeElement<tgt::Vector2<T> >::swapEndianness(tgt::Vector2<T> v) {
    v[0] = VolumeElement<T>::swapEndianness(v[0]);
    v[1] = VolumeElement<T>::swapEndianness(v[1]);
    return v;
}

/**
* Calculates the squared difference between val1 and val2. Since val1 and val2 are 3-dimensional
* vectors, each squared difference is multiplied by 1/3.
*/
template<class T>
double VolumeElement<tgt::Vector2<T> >::calcSquaredDifference(tgt::Vector2<T> val1, tgt::Vector2<T> val2) {
    double result;
    tgt::Vector2d difference;

    difference.x = (double)val1.x - (double)val2.x;
    difference.y = (double)val1.y - (double)val2.y;

    difference.x = difference.x * difference.x;
    difference.y = difference.y * difference.y;

    result = difference.x + difference.y;

    return result;
}

// --- Vector3 ----

template<class T>
inline tgt::Vector3<T> VolumeElement<tgt::Vector3<T> >::getZero() {
    return tgt::Vector3<T>(T(0), T(0), T(0));
}

template<class T>
inline T VolumeElement<tgt::Vector3<T> >::getChannel(const tgt::Vector3<T>& t, size_t channel) {
    tgtAssert(channel < 3, "channel must be less than 3");
    if (channel < 3)
        return t[channel];
    else
        return 0;
}

template<class T>
inline void VolumeElement<tgt::Vector3<T> >::setChannel(T value, tgt::Vector3<T>& t, size_t channel) {
    tgtAssert(channel < 3, "channel must be less than 3");
    t[channel] = value;
}

template<class T>
inline bool VolumeElement<tgt::Vector3<T> >::isSigned() {
    return isTypeSigned<T>();
}

template<class T>
inline bool VolumeElement<tgt::Vector3<T> >::isInteger() {
    return isTypeInteger<T>();
}

template<class T>
inline tgt::Vector3<T> VolumeElement<tgt::Vector3<T> >::rangeMax() {
    return tgt::Vector3<T>(getTypeUpperLimit<T>());
}

template<class T>
inline tgt::Vector3<T> VolumeElement<tgt::Vector3<T> >::rangeMin() {
    return tgt::Vector3<T>(getTypeLowerLimit<T>());
}

/**
* Returns the maximum value an element can have. In contrast to max() no
* vector is returned but the maximum value an element in the vector can have.
*/
template<class T>
inline float VolumeElement<tgt::Vector3<T> >::rangeMaxElement() {
    return static_cast<float>(getTypeUpperLimit<T>());
}

template<class T>
inline float VolumeElement<tgt::Vector3<T> >::rangeMinElement() {
    return static_cast<float>(getTypeLowerLimit<T>());
}

template<class T>
float voreen::VolumeElement<tgt::Vector3<T> >::getMaxElement(const tgt::Vector3<T>& value) {
    return static_cast<float>(value.x);
}

template<class T>
float voreen::VolumeElement<tgt::Vector3<T> >::getMinElement(const tgt::Vector3<T>& value) {
    return static_cast<float>(value.y);
}

template<class T>
inline int VolumeElement<tgt::Vector3<T> >::getNumChannels() {
    return 3;
}

template<class T>
inline tgt::Vector3<T> VolumeElement<tgt::Vector3<T> >::swapEndianness(tgt::Vector3<T> v) {
    v[0] = VolumeElement<T>::swapEndianness(v[0]);
    v[1] = VolumeElement<T>::swapEndianness(v[1]);
    v[2] = VolumeElement<T>::swapEndianness(v[2]);
    return v;
}

/**
* Calculates the squared difference between val1 and val2. Since val1 and val2 are 3-dimensional
* vectors, each squared difference is multiplied by 1/3.
*/
template<class T>
double VolumeElement<tgt::Vector3<T> >::calcSquaredDifference(tgt::Vector3<T> val1, tgt::Vector3<T> val2) {
    double result;
    tgt::Vector3d difference;

    difference.x = (double)val1.x - (double)val2.x;
    difference.y = (double)val1.y - (double)val2.y;
    difference.z = (double)val1.z - (double)val2.z;

    difference.x = difference.x * difference.x;
    difference.y = difference.y * difference.y;
    difference.z = difference.z * difference.z;

    result = difference.x + difference.y + difference.z;

    return result;
}

// --- Vector4 ----

template<class T>
inline tgt::Vector4<T> VolumeElement<tgt::Vector4<T> >::getZero() {
    return tgt::Vector4<T>(T(0), T(0), T(0), T(0));
}

template<class T>
inline T VolumeElement<tgt::Vector4<T> >::getChannel(const tgt::Vector4<T>& t, size_t channel) {
    tgtAssert(channel < 4, "channel must be less than 4");
    if (channel < 4)
        return t[channel];
    else
        return 0;
}

template<class T>
inline void VolumeElement<tgt::Vector4<T> >::setChannel(T value, tgt::Vector4<T>& t, size_t channel) {
    tgtAssert(channel < 4, "channel must be less than 4");
    t[channel] = value;
}

template<class T>
inline bool VolumeElement<tgt::Vector4<T> >::isSigned() {
    return isTypeSigned<T>();
}

template<class T>
inline bool VolumeElement<tgt::Vector4<T> >::isInteger() {
    return isTypeInteger<T>();
}

template<class T>
inline tgt::Vector4<T> VolumeElement<tgt::Vector4<T> >::rangeMax() {
    return tgt::Vector4<T>(getTypeUpperLimit<T>());
}

template<class T>
inline tgt::Vector4<T> VolumeElement<tgt::Vector4<T> >::rangeMin() {
    return tgt::Vector4<T>(getTypeLowerLimit<T>());
}

/**
* Returns the maximum value an element can have. In contrast to max() no
* vector is returned but the maximum value an element in the vector can have.
*/
template<class T>
inline float VolumeElement<tgt::Vector4<T> >::rangeMaxElement() {
    return static_cast<float>(getTypeUpperLimit<T>());
}

template<class T>
inline float VolumeElement<tgt::Vector4<T> >::rangeMinElement() {
    return static_cast<float>(getTypeLowerLimit<T>());
}

template<class T>
float voreen::VolumeElement<tgt::Vector4<T> >::getMaxElement(const tgt::Vector4<T>& value) {
    return static_cast<float>(value.x);
}

template<class T>
float voreen::VolumeElement<tgt::Vector4<T> >::getMinElement(const tgt::Vector4<T>& value) {
    return static_cast<float>(value.x);
}

template<class T>
inline int VolumeElement<tgt::Vector4<T> >::getNumChannels() {
    return 4;
}

template<class T>
inline tgt::Vector4<T> VolumeElement<tgt::Vector4<T> >::swapEndianness(tgt::Vector4<T> v) {
    v[0] = VolumeElement<T>::swapEndianness(v[0]);
    v[1] = VolumeElement<T>::swapEndianness(v[1]);
    v[2] = VolumeElement<T>::swapEndianness(v[2]);
    v[3] = VolumeElement<T>::swapEndianness(v[3]);
    return v;
}

/**
* Calculates the squared difference between val1 and val2. Since val1 and val2 are 4-dimensional
* vectors, each squared difference is multiplied by 1/4.
*/
template<class T>
double VolumeElement<tgt::Vector4<T> >::calcSquaredDifference(tgt::Vector4<T> val1, tgt::Vector4<T> val2) {
    double result;
    tgt::Vector4d difference;

    difference.x = static_cast<double>(val1.x) - static_cast<double>(val2.x);
    difference.y = static_cast<double>(val1.y) - static_cast<double>(val2.y);
    difference.z = static_cast<double>(val1.z) - static_cast<double>(val2.z);
    difference.w = static_cast<double>(val1.w) - static_cast<double>(val2.w);

    difference.x = difference.x*difference.x;
    difference.y = difference.y*difference.y;
    difference.z = difference.z*difference.z;
    difference.w = difference.w*difference.w;

    result = difference.x + difference.y + difference.z + difference.w;

    return result;
}

// --- Matrix3 ----

template<class T>
inline tgt::Matrix3<T> VolumeElement<tgt::Matrix3<T> >::getZero() {
    return tgt::Matrix3<T>(T(0));
}

template<class T>
inline T VolumeElement<tgt::Matrix3<T> >::getChannel(const tgt::Matrix3<T>& t, size_t channel) {
    tgtAssert(channel < 9, "channel must be less than 9");
    if (channel < 9)
        return t.elem[channel];
    else
        return 0;
}

template<class T>
inline void VolumeElement<tgt::Matrix3<T> >::setChannel(T value, tgt::Matrix3<T>& t, size_t channel) {
    tgtAssert(channel < 9, "channel must be less than 9");
    t.elem[channel] = value;
}

template<class T>
inline bool VolumeElement<tgt::Matrix3<T> >::isSigned() {
    return isTypeSigned<T>();
}

template<class T>
inline bool VolumeElement<tgt::Matrix3<T> >::isInteger() {
    return isTypeInteger<T>();
}

template<class T>
inline tgt::Matrix3<T> VolumeElement<tgt::Matrix3<T> >::rangeMax() {
    return tgt::Matrix3<T>(getTypeUpperLimit<T>());
}

template<class T>
inline tgt::Matrix3<T> VolumeElement<tgt::Matrix3<T> >::rangeMin() {
    return tgt::Matrix3<T>(getTypeLowerLimit<T>());
}

/**
* Returns the maximum value an element can have. In contrast to max() no
* vector is returned but the maximum value an element in the vector can have.
*/
template<class T>
inline float VolumeElement<tgt::Matrix3<T> >::rangeMaxElement() {
    return static_cast<float>(getTypeUpperLimit<T>());
}

template<class T>
inline float VolumeElement<tgt::Matrix3<T> >::rangeMinElement() {
    return static_cast<float>(getTypeLowerLimit<T>());
}

template<class T>
float voreen::VolumeElement<tgt::Matrix3<T> >::getMaxElement(const tgt::Matrix3<T>& value) {
    return static_cast<float>(value[0]);//TODO: correct?
}

template<class T>
float voreen::VolumeElement<tgt::Matrix3<T> >::getMinElement(const tgt::Matrix3<T>& value) {
    return static_cast<float>(value[0]);//TODO: correct?
}

template<class T>
inline int VolumeElement<tgt::Matrix3<T> >::getNumChannels() {
    return 9;
}

template<class T>
inline tgt::Matrix3<T> VolumeElement<tgt::Matrix3<T> >::swapEndianness(tgt::Matrix3<T> v) {
    for(int i=0; i<VolumeElement<tgt::Matrix3<T> >::getNumChannels(); i++)
        v[i] = VolumeElement<T>::swapEndianness(v[i]);
    return v;
}

/**
* Calculates the squared difference between val1 and val2. Since val1 and val2 are 3-dimensional
* vectors, each squared difference is multiplied by 1/3.
*/
template<class T>
double VolumeElement<tgt::Matrix3<T> >::calcSquaredDifference(tgt::Matrix3<T> val1, tgt::Matrix3<T> val2) {
    double result = 0.0;
    //FIXME
    //tgt::Vector3d difference;

    //difference.x = (double)val1.x - (double)val2.x;
    //difference.y = (double)val1.y - (double)val2.y;
    //difference.z = (double)val1.z - (double)val2.z;

    //difference.x = difference.x * difference.x;
    //difference.y = difference.y * difference.y;
    //difference.z = difference.z * difference.z;

    //result = difference.x + difference.y + difference.z;

    return result;
}

// --- Matrix4 ----

template<class T>
inline tgt::Matrix4<T> VolumeElement<tgt::Matrix4<T> >::getZero() {
    return tgt::Matrix4<T>(T(0));
}

template<class T>
inline T VolumeElement<tgt::Matrix4<T> >::getChannel(const tgt::Matrix4<T>& t, size_t channel) {
    tgtAssert(channel < 16, "channel must be less than 16");
    if (channel < 16)
        return t.elem[channel];
    else
        return 0;
}

template<class T>
inline void VolumeElement<tgt::Matrix4<T> >::setChannel(T value, tgt::Matrix4<T>& t, size_t channel) {
    tgtAssert(channel < 16, "channel must be less than 16");
    t.elem[channel] = value;
}

template<class T>
inline bool VolumeElement<tgt::Matrix4<T> >::isSigned() {
    return isTypeSigned<T>();
}

template<class T>
inline bool VolumeElement<tgt::Matrix4<T> >::isInteger() {
    return isTypeInteger<T>();
}

template<class T>
inline tgt::Matrix4<T> VolumeElement<tgt::Matrix4<T> >::rangeMax() {
    return tgt::Matrix4<T>(getTypeUpperLimit<T>());
}

template<class T>
inline tgt::Matrix4<T> VolumeElement<tgt::Matrix4<T> >::rangeMin() {
    return tgt::Matrix4<T>(getTypeLowerLimit<T>());
}

/**
* Returns the maximum value an element can have. In contrast to max() no
* vector is returned but the maximum value an element in the vector can have.
*/
template<class T>
inline float VolumeElement<tgt::Matrix4<T> >::rangeMaxElement() {
    return static_cast<float>(getTypeUpperLimit<T>());
}

template<class T>
inline float VolumeElement<tgt::Matrix4<T> >::rangeMinElement() {
    return static_cast<float>(getTypeLowerLimit<T>());
}

template<class T>
float voreen::VolumeElement<tgt::Matrix4<T> >::getMaxElement(const tgt::Matrix4<T>& value) {
    return static_cast<float>(value[0]);//TODO: correct?
}

template<class T>
float voreen::VolumeElement<tgt::Matrix4<T> >::getMinElement(const tgt::Matrix4<T>& value) {
    return static_cast<float>(value[0]);//TODO: correct?
}

template<class T>
inline int VolumeElement<tgt::Matrix4<T> >::getNumChannels() {
    return 16;
}

template<class T>
inline tgt::Matrix4<T> VolumeElement<tgt::Matrix4<T> >::swapEndianness(tgt::Matrix4<T> v) {
    for(int i=0; i<VolumeElement<tgt::Matrix4<T> >::getNumChannels(); i++)
        v[i] = VolumeElement<T>::swapEndianness(v[i]);
    return v;
}

/**
* Calculates the squared difference between val1 and val2. Since val1 and val2 are 3-dimensional
* vectors, each squared difference is multiplied by 1/3.
*/
template<class T>
double VolumeElement<tgt::Matrix4<T> >::calcSquaredDifference(tgt::Matrix4<T> val1, tgt::Matrix4<T> val2) {
    //FIXME
    double result = 0.0;
    //tgt::Vector3d difference;

    //difference.x = (double)val1.x - (double)val2.x;
    //difference.y = (double)val1.y - (double)val2.y;
    //difference.z = (double)val1.z - (double)val2.z;

    //difference.x = difference.x * difference.x;
    //difference.y = difference.y * difference.y;
    //difference.z = difference.z * difference.z;

    //result = difference.x + difference.y + difference.z;

    return result;
}


// --- Tensor2 ----

template<class T>
inline Tensor2<T> VolumeElement<Tensor2<T> >::getZero() {
    return Tensor2<T>(T(0));
}

template<class T>
inline T VolumeElement<Tensor2<T> >::getChannel(const Tensor2<T>& t, size_t channel) {
    tgtAssert(channel < 6, "channel must be less than 6");
    if (channel < 6)
        return t.elem[channel];
    else
        return 0;
}

template<class T>
inline void VolumeElement<Tensor2<T> >::setChannel(T value, Tensor2<T>& t, size_t channel) {
    tgtAssert(channel < 6, "channel must be less than 6");
    t.elem[channel] = value;
}

template<class T>
inline bool VolumeElement<Tensor2<T> >::isSigned() {
    return isTypeSigned<T>();
}

template<class T>
inline bool VolumeElement<Tensor2<T> >::isInteger() {
    return isTypeInteger<T>();
}

template<class T>
inline Tensor2<T> VolumeElement<Tensor2<T> >::rangeMax() {
    return Tensor2<T>(getTypeUpperLimit<T>());
}

template<class T>
inline Tensor2<T> VolumeElement<Tensor2<T> >::rangeMin() {
    return Tensor2<T>(getTypeLowerLimit<T>());
}

/**
* Returns the maximum value an element can have. In contrast to max() no
* vector is returned but the maximum value an element in the vector can have.
*/
template<class T>
inline float VolumeElement<Tensor2<T> >::rangeMaxElement() {
    return static_cast<float>(getTypeUpperLimit<T>());
}

template<class T>
inline float VolumeElement<Tensor2<T> >::rangeMinElement() {
    return static_cast<float>(getTypeLowerLimit<T>());
}

template<class T>
float voreen::VolumeElement<Tensor2<T> >::getMaxElement(const Tensor2<T>& value) {
    return static_cast<float>(value[0]);//TODO: correct?
}

template<class T>
float voreen::VolumeElement<Tensor2<T> >::getMinElement(const Tensor2<T>& value) {
    return static_cast<float>(value[0]);//TODO: correct?
}

template<class T>
inline int VolumeElement<Tensor2<T> >::getNumChannels() {
    return 6;
}

template<class T>
inline Tensor2<T> VolumeElement<Tensor2<T> >::swapEndianness(Tensor2<T> v) {
    for(int i=0; i<VolumeElement<Tensor2<T> >::getNumChannels(); i++)
        v[i] = VolumeElement<T>::swapEndianness(v[i]);
    return v;
}

/**
* Calculates the squared difference between val1 and val2. Since val1 and val2 are 6-dimensional
* vectors, each squared difference is multiplied by 1/6.
*/
template<class T>
double VolumeElement<Tensor2<T> >::calcSquaredDifference(Tensor2<T> val1, Tensor2<T> val2) {
    double result;
    Tensor2<double> difference;

    difference.Dxx = static_cast<double>(val1.Dxx) - static_cast<double>(val2.Dxx);
    difference.Dxy = static_cast<double>(val1.Dxy) - static_cast<double>(val2.Dxy);
    difference.Dxz = static_cast<double>(val1.Dxz) - static_cast<double>(val2.Dxz);
    difference.Dyy = static_cast<double>(val1.Dyy) - static_cast<double>(val2.Dyy);
    difference.Dyz = static_cast<double>(val1.Dyz) - static_cast<double>(val2.Dyz);
    difference.Dzz = static_cast<double>(val1.Dzz) - static_cast<double>(val2.Dzz);

    difference.Dxx = difference.Dxx*difference.Dxx;
    difference.Dxy = difference.Dxy*difference.Dxy;
    difference.Dxz = difference.Dxz*difference.Dxz;
    difference.Dyy = difference.Dyy*difference.Dyy;
    difference.Dyz = difference.Dyz*difference.Dyz;
    difference.Dzz = difference.Dzz*difference.Dzz;

    result = difference.Dxx + difference.Dxy + difference.Dxz + difference.Dyy + difference.Dyz + difference.Dzz;

    return result;
}

template<>
float getTypeAsFloat(float value) {
    return value;
}

template<>
float getTypeAsFloat(double value) {
    return static_cast<float>(value);
}

template<>
float getTypeAsFloat(int8_t value) {
    if (value >= 0)
        return static_cast<float>(value) / VolumeElement<int8_t>::rangeMaxElement();
    else
        return static_cast<float>(value) / -VolumeElement<int8_t>::rangeMinElement();
}

template<>
float getTypeAsFloat(uint8_t value) {
    return static_cast<float>(value) / VolumeElement<uint8_t>::rangeMaxElement();
}

template<>
float getTypeAsFloat(int16_t value) {
    if (value >= 0)
        return static_cast<float>(value) / VolumeElement<int16_t>::rangeMaxElement();
    else
        return static_cast<float>(value) / -VolumeElement<int16_t>::rangeMinElement();
}

template<>
float getTypeAsFloat(uint16_t value) {
    return static_cast<float>(value) / VolumeElement<uint16_t>::rangeMaxElement();
}

template<>
float getTypeAsFloat(int32_t value) {
    if (value >= 0)
        return static_cast<float>(value) / VolumeElement<int32_t>::rangeMaxElement();
    else
        return static_cast<float>(value) / -VolumeElement<int32_t>::rangeMinElement();
}

template<>
float getTypeAsFloat(uint32_t value) {
    return static_cast<float>(value) / VolumeElement<uint32_t>::rangeMaxElement();
}

template<>
float getTypeAsFloat(int64_t value) {
    if (value >= 0)
        return static_cast<float>(value) / VolumeElement<int64_t>::rangeMaxElement();
    else
        return static_cast<float>(value) / -VolumeElement<int64_t>::rangeMinElement();
}

template<>
float getTypeAsFloat(uint64_t value) {
    return static_cast<float>(value) / VolumeElement<uint64_t>::rangeMaxElement();
}

template<>
float getFloatAsType(float value) {
    return value;
}
template<>
double getFloatAsType(float value) {
    return static_cast<double>(value);
}

template<>
int8_t getFloatAsType(float value) {
    value = tgt::clamp(value, -1.0f, 1.0f);
    if(value >= 0.0f)
        return static_cast<int8_t>(value*VolumeElement<int8_t>::rangeMaxElement());
    else
        return static_cast<int8_t>(value*-VolumeElement<int8_t>::rangeMinElement());
}

template<>
uint8_t getFloatAsType(float value) {
    value = tgt::clamp(value, 0.0f, 1.0f);
    return static_cast<uint8_t>(value*VolumeElement<uint8_t>::rangeMaxElement());
}

template<>
int16_t getFloatAsType(float value) {
    value = tgt::clamp(value, -1.0f, 1.0f);
    if(value >= 0.0f)
        return static_cast<int16_t>(value*VolumeElement<int16_t>::rangeMaxElement());
    else
        return static_cast<int16_t>(value*-VolumeElement<int16_t>::rangeMinElement());
}
template<>
uint16_t getFloatAsType(float value) {
    value = tgt::clamp(value, 0.0f, 1.0f);
    return static_cast<uint16_t>(value*VolumeElement<uint16_t>::rangeMaxElement());
}

template<>
int32_t getFloatAsType(float value) {
    // uint32_MAX does not fit into a float => We lose precision and an input value of 1.0f is mapped onto uint32_MAX+1 == 0
    // Therefore we cast value to the larger int64_t and clamp then.
    if(value >= 0.0f) {
        int64_t vali = value*VolumeElement<int32_t>::rangeMaxElement();
        return std::min(vali, static_cast<int64_t>(std::numeric_limits<int32_t>::max()));
    } else {
        int64_t vali = value*-VolumeElement<int32_t>::rangeMinElement();
        return std::max(vali, static_cast<int64_t>(std::numeric_limits<int32_t>::min()));
    }
}
template<>
uint32_t getFloatAsType(float value) {
    // uint32_MAX does not fit into a float => We lose precision and an input value of 1.0f is mapped onto uint32_MAX+1 == 0
    // Therefore we cast value to the larger int64_t and clamp then.
    int64_t vali = value*VolumeElement<uint32_t>::rangeMaxElement();
    return tgt::clamp(vali, int64_t(0), static_cast<int64_t>(std::numeric_limits<uint32_t>::max()));
}

template<>
int64_t getFloatAsType(float value) {
    value = tgt::clamp(value, -1.0f, 1.0f);
    if(value >= 0.0f)
        return static_cast<int64_t>(value*VolumeElement<int64_t>::rangeMaxElement());
    else
        return static_cast<int64_t>(-value*VolumeElement<int64_t>::rangeMinElement());
}
template<>
uint64_t getFloatAsType(float value) {
    value = tgt::clamp(value, 0.0f, 1.0f);
    return static_cast<uint64_t>(value*VolumeElement<uint64_t>::rangeMaxElement());
}

} // namespace voreen

#endif // VRN_VOLUMEELEMENT
