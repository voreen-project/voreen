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

#ifndef VRN_HDF5_UTILS_H
#define VRN_HDF5_UTILS_H

#if defined(__clang__)
// Get rid of stupid "hidden overloaded virtual function" warnings in hdf5 lib with clang
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Woverloaded-virtual"
#include "H5Cpp.h"
#pragma clang diagnostic pop
#else
#include "H5Cpp.h"
#endif

#include "voreen/core/datastructures/volume/volume.h"

#include <tgt/vector.h>
#include <tgt/exception.h>
#include <boost/thread.hpp>

// Apparently semver is not a thing...
#if H5_VERSION_GE(1, 10, 1) || H5_VERSION_GE(1, 8, 21) && !H5_VERSION_GE(1, 10, 0)
#define H5_STUPID_LOCATION_API_CHANGES
#endif

namespace voreen {

/**
 * HDF5-cxx-libs are NOT threadsafe!
 *
 * Therefore use this mutex to avoid concurrency when using anything
 * from H5Cpp.h or similar.
 */
extern boost::recursive_mutex hdf5libMutex;


/**
 * Get the respective H5::PredType for the given template Parameter.
 * The general case is not implemented and thus generates a linker
 * error when using this function for a type we did not specify a HDF5
 * type for (yet).
 * @return PredType for given native type.
 * @note: This implementation Assumes little endianess.
 */
template<typename C> H5::PredType getPredType();
template<> inline H5::PredType getPredType<float>();
template<> inline H5::PredType getPredType<double>();
template<> inline H5::PredType getPredType<uint8_t>();
template<> inline H5::PredType getPredType<uint16_t>();
template<> inline H5::PredType getPredType<uint32_t>();
template<> inline H5::PredType getPredType<uint64_t>();
template<> inline H5::PredType getPredType<int8_t>();
template<> inline H5::PredType getPredType<int16_t>();
template<> inline H5::PredType getPredType<int32_t>();
template<> inline H5::PredType getPredType<int64_t>();

/**
 * Get the appropriate HDF5 PredType by the Voreen style base type string.
 * @baseType Voreen style base type string.
 * @return PredType for the base type string.
 */
H5::PredType getPredType(const std::string& baseType);

/**
 * Get the appropriate Voreen style base type string for a given HDF5 PredType.
 * @format HDF5 PredType to be converted to Voreen style base type string.
 * @return Voreen style base type string according to the given Predtype.
 */
std::string getBaseTypeFromDataType(const H5::DataType& type);

/**
 * Create a HDF5 DataSpace suitable to represent and store the data of the given volume.
 * @return A HDF5-DataSpace according to the given voreen style parameters.
 */
H5::DataSpace createDataSpace(tgt::svec3 volumeDimensions, size_t numberOfChannels = 1);

/**
 * Create an array attribute at the given location and write values to it.
 */
template<typename T>
void writeArrayAttribute(
#ifdef H5_STUPID_LOCATION_API_CHANGES
        const H5::H5Object& loc,
#else
        const H5::H5Location& loc,
#endif
        const H5std_string& name, const T* values, hsize_t numberOfValues);


/**
 * Create a vec3 attribute at the given location and write vec to it.
 */
template<typename T>
void writeVec3Attribute(
#ifdef H5_STUPID_LOCATION_API_CHANGES
        const H5::H5Object& loc,
#else
        const H5::H5Location& loc,
#endif
        const H5std_string& name, const tgt::Vector3<T>& vec);


/**
 * Create a scalar attribute at the given location and write scalar to it.
 */
template<typename T>
void writeScalarAttribute(
#ifdef H5_STUPID_LOCATION_API_CHANGES
        const H5::H5Object& loc,
#else
        const H5::H5Location& loc,
#endif
        const H5std_string& name, const T scalar);


/**
 * Create a string attribute at the given location and write str to it.
 */
void writeStringAttribute(
#ifdef H5_STUPID_LOCATION_API_CHANGES
        const H5::H5Object& loc,
#else
        const H5::H5Location& loc,
#endif
        const H5std_string& name, const H5std_string& str);


/**
 * Read an array attribute with the given name from the given location.
 */
template<typename T>
void readArrayAttribute(
#ifdef H5_STUPID_LOCATION_API_CHANGES
        const H5::H5Object& loc,
#else
        const H5::H5Location& loc,
#endif
        const H5std_string& name, T* values, hsize_t numValues);


/**
 * Read a vec3 attribute with the given name from the given location.
 */
template<typename T>
tgt::Vector3<T> readVec3Attribute(
#ifdef H5_STUPID_LOCATION_API_CHANGES
        const H5::H5Object& loc,
#else
        const H5::H5Location& loc,
#endif
        const H5std_string name);


/**
 * Read a scalar attribute with the given name from the given location.
 */
template<typename T>
T readScalarAttribute(
#ifdef H5_STUPID_LOCATION_API_CHANGES
        const H5::H5Object& loc,
#else
        const H5::H5Location& loc,
#endif
        const H5std_string name);


/**
 * Read a string attribute with the given name from the given location.
 */
H5std_string readStringAttribute(
#ifdef H5_STUPID_LOCATION_API_CHANGES
        const H5::H5Object& loc,
#else
        const H5::H5Location& loc,
#endif
        const H5std_string name);


/**
 * Convert HDF5 style vector to tgt vector.
 * @param vec HDF5 style array vector representation.
 * @return a tgt size vector equivalent to the HDF5 representation.
 */
template<typename t>
inline tgt::Vector2<t> vec2HDF5ToTgt(const t vec[2]);
template<typename t>
inline tgt::Vector3<t> vec3HDF5ToTgt(const t vec[3]);
template<typename t>
inline tgt::Vector4<t> vec4HDF5ToTgt(const t vec[4]);

/**
 * Convert tgt vector to HDF5 style vector.
 * @param v a tgt size vector.
 * @param vec HDF5 style array vector representation the data from
 *          v will be written to.
 */
template<typename t>
inline void vec2TgtToHDF5(const tgt::Vector2<t>& v, t vec[2]);
template<typename t>
inline void vec3TgtToHDF5(const tgt::Vector3<t>& v, t vec[3]);
template<typename t>
inline void vec4TgtToHDF5(const tgt::Vector4<t>& v, t vec[4]);

/**
 * Convert HDF5 style dimensions to tgt dimensions.
 * @param dimensions HDF5 style array dimension representation.
 * @return a tgt size vector equivalent to the HDF5 representation.
 */
tgt::svec2 vec2HDF5ToTgt(const hsize_t dimensions[2]);
tgt::svec3 vec3HDF5ToTgt(const hsize_t dimensions[3]);
tgt::svec4 vec4HDF5ToTgt(const hsize_t dimensions[4]);

/**
 * Convert tgt dimensions to HDF5 style dimensions.
 * @param a tgt size vector.
 * @param dimensions HDF5 style array dimension representation the data from
 *          vec will be written to.
 */
void vec2TgtToHDF5(const tgt::svec2& vec, hsize_t dimensions[2]);
void vec3TgtToHDF5(const tgt::svec3& vec, hsize_t dimensions[3]);
void vec4TgtToHDF5(const tgt::svec4& vec, hsize_t dimensions[4]);


// ====== Implementation of template functions =================================================================
template<typename T>
void writeArrayAttribute(
#ifdef H5_STUPID_LOCATION_API_CHANGES
        const H5::H5Object& loc,
#else
        const H5::H5Location& loc,
#endif
        const H5std_string& name, const T* values, hsize_t numberOfValues) {

    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);
    // Create Attribute
    const H5::DataSpace dataspace(1, &numberOfValues);
    H5::Attribute attribute = loc.createAttribute(name, getPredType<T>(), dataspace);

    attribute.write(getPredType<T>(), values);
}

template<typename T>
void writeVec3Attribute(
#ifdef H5_STUPID_LOCATION_API_CHANGES
        const H5::H5Object& loc,
#else
        const H5::H5Location& loc,
#endif
        const H5std_string& name, const tgt::Vector3<T>& vec) {

    // Convert tgt vector to hdf style
    T h5vec[3];
    vec3TgtToHDF5(vec, h5vec);

    // Write the attribute
    writeArrayAttribute(loc, name, h5vec, 3);
}

template<typename T>
void writeScalarAttribute(
#ifdef H5_STUPID_LOCATION_API_CHANGES
        const H5::H5Object& loc,
#else
        const H5::H5Location& loc,
#endif
        const H5std_string& name, const T scalar) {

    //Write the attribute
    writeArrayAttribute(loc, name, &scalar, 1);
}

template<typename T>
void readArrayAttribute(
#ifdef H5_STUPID_LOCATION_API_CHANGES
        const H5::H5Object& loc,
#else
        const H5::H5Location& loc,
#endif
        const H5std_string& name, T* values, hsize_t numValues) {

    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);
    // Try to open element size attribute.
    H5::Attribute attr = loc.openAttribute(name);
    H5::DataSpace dataSpace = attr.getSpace();

    // Check dimensions
    if(static_cast<hsize_t>(dataSpace.getSimpleExtentNpoints()) != numValues) {
        throw H5::AttributeIException("readArrayAttribute", "Attribute to be read from HDF5 location contains " + std::to_string(dataSpace.getSimpleExtentNpoints()) + "Element(s). Expected: " + std::to_string(numValues) + ".");
    }

    // Check type
    H5::PredType expectedType = getPredType<T>();
    if(!(attr.getDataType() == expectedType)) {
        throw H5::AttributeIException("readArrayAttribute", "Attribute's type is not float.");
    }
    // Read attribute
    attr.read(expectedType, values);
}

template<typename T>
tgt::Vector3<T> readVec3Attribute(
#ifdef H5_STUPID_LOCATION_API_CHANGES
        const H5::H5Object& loc,
#else
        const H5::H5Location& loc,
#endif
        const H5std_string name) {

    T values[3];
    readArrayAttribute(loc, name, values, 3);
    return vec3HDF5ToTgt(values);
}

template<typename T>
T readScalarAttribute(
#ifdef H5_STUPID_LOCATION_API_CHANGES
        const H5::H5Object& loc,
#else
        const H5::H5Location& loc,
#endif
        const H5std_string name) {

    T value;
    readArrayAttribute(loc, name, &value, 1);
    return value;
}

template<typename t>
inline tgt::Vector2<t> vec2HDF5ToTgt(const t vec[2]) {
    // HDF5 Style indexing (slowest changing first): yx!
    return tgt::Vector2<t>(vec[1], vec[0]);
}

template<typename t>
inline tgt::Vector3<t> vec3HDF5ToTgt(const t vec[3]) {
    // HDF5 Style indexing (slowest changing first): zyx!
    return tgt::Vector3<t>(vec[2], vec[1], vec[0]);
}

template<typename t>
inline tgt::Vector4<t> vec4HDF5ToTgt(const t vec[4]) {
    // HDF5 Style indexing (slowest changing first): wzyx!
    return tgt::Vector4<t>(vec[3], vec[2], vec[1], vec[0]);
}

template<typename t>
inline void vec2TgtToHDF5(const tgt::Vector2<t>& v, t vec[2]) {
    // HDF5 Style indexing (slowest changing first): yx!
    vec[0] = v.y;
    vec[1] = v.x;
}

template<typename t>
inline void vec3TgtToHDF5(const tgt::Vector3<t>& v, t vec[3]) {
    // HDF5 Style indexing (slowest changing first): zyx!
    vec[0] = v.z;
    vec[1] = v.y;
    vec[2] = v.x;
}

template<typename t>
inline void vec4TgtToHDF5(const tgt::Vector4<t>& v, t vec[4]) {
    // HDF5 Style indexing (slowest changing first): wzyx!
    vec[0] = v.w;
    vec[1] = v.z;
    vec[2] = v.y;
    vec[3] = v.x;
}

template<>
H5::PredType getPredType<float>() {
    return H5::PredType::IEEE_F32LE;
}
template<>
H5::PredType getPredType<double>() {
    return H5::PredType::IEEE_F64LE;
}
template<>
H5::PredType getPredType<uint8_t>() {
    return H5::PredType::STD_U8LE;
}
template<>
H5::PredType getPredType<uint16_t>() {
    return H5::PredType::STD_U16LE;
}
template<>
H5::PredType getPredType<uint32_t>() {
    return H5::PredType::STD_U32LE;
}
template<>
H5::PredType getPredType<uint64_t>() {
    return H5::PredType::STD_U64LE;
}
template<>
H5::PredType getPredType<int8_t>() {
    return H5::PredType::STD_I8LE;
}
template<>
H5::PredType getPredType<int16_t>() {
    return H5::PredType::STD_I16LE;
}
template<>
H5::PredType getPredType<int32_t>() {
    return H5::PredType::STD_I32LE;
}
template<>
H5::PredType getPredType<int64_t>() {
    return H5::PredType::STD_I64LE;
}
} // namespace voreen
#endif //VRN_HDF5_UTILS_H
