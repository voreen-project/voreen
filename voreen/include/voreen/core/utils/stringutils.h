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

#ifndef VRN_STRINGUTILS_H
#define VRN_STRINGUTILS_H

#include "ext/tgt/vector.h"
#include "voreen/core/voreencoreapi.h"
#include "voreen/core/utils/exception.h"

#include "tgt/vector.h"

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

namespace voreen {

/// Converts an int to a string.
VRN_CORE_API std::string itos(int i, int stringLength = -1, char fillChar = '0');

#ifdef __APPLE__
/// Long unsigned int / size_t is an own type on MacOS, provide special implementation to avoid ambiguousness.
VRN_CORE_API std::string itos(long unsigned int i, int stringLength = -1, char fillChar = '0');
#endif

/// Converts an uint32_t to a string.
VRN_CORE_API std::string itos(uint32_t i, int stringLength = -1, char fillChar = '0');

/// Converts an uint64_t to a string.
VRN_CORE_API std::string itos(uint64_t i, int stringLength = -1, char fillChar = '0');

/**
 * Converts a float to a string.
 *
 * @param precision number of decimals to print.
 *  For precision=-1, sprintf's standard floating point formatting is used (%f).
 */
VRN_CORE_API std::string ftos(float f, int precision = -1);

/**
 * Converts a double to a string.
 *
 * @param precision number of decimals to print.
 *  For precision=-1, sprintf's standard floating point formatting is used (%f).
 */
VRN_CORE_API std::string dtos(double d, int precision = -1);

/**
 * Converts the string to a null-terminated char array with length s.size()+1.
 * Deleting the allocated memory is up to the caller.
 */
VRN_CORE_API char* strToChr(const std::string& s);

/// Generic to-string conversion using a stringstream.
template<class T>
std::string genericToString(const T& value);

/**
 * Generic from-string conversion using a stringstream.
 *
 * @throw VoreenException if the conversion failed
 */
template<class T>
T genericFromString(const std::string& str);

/// MSVC++ 2010 and later already provides the sto* functions in the <string> header
#if (!defined(_MSC_VER) || (_MSC_VER < 1600)) && !defined(__GXX_EXPERIMENTAL_CXX0X__) && !defined(__APPLE__)
/// Converts a string to an int using stringstreams
VRN_CORE_API int stoi(const std::string& s);

/// Converts a string to a float using stringstreams
VRN_CORE_API float stof(const std::string& s);

/// Converts a string to a double using stringstreams
VRN_CORE_API double stod(const std::string& s);
#endif // !defined(_MSC_VER) || (_MSC_VER < 1600)

/**
 * Returns a copy of \p str where all occurrences of \p from have been replaced by \p to.
 */
VRN_CORE_API std::string strReplaceAll(const std::string& str, const std::string& from, const std::string& to);

/// Returns a copy of \p str where the first occurrence of \p from has been replaced by \p to.
VRN_CORE_API std::string strReplaceFirst(const std::string& str, const std::string& from, const std::string& to);

/// Returns a copy of \p str where the last occurrence of \p from has been replaced by \p to.
VRN_CORE_API std::string strReplaceLast(const std::string& str, const std::string& from, const std::string& to);

/**
 * Removes whitespaces from beginning and end of a string.
 *
 * @param str Input string.
 * @param charlist Characters to remove, defaults to space, tab, newline, carriage return, 0, vertical tab.
 */
VRN_CORE_API std::string trim(std::string str, const std::string& charlist = " \t\n\r\0\x0B");

/// Converts a string to lower case.
VRN_CORE_API std::string toLower(const std::string& str);

/// Converts a string to upper case.
VRN_CORE_API std::string toUpper(const std::string& str);

/**
 * Splits a string by the specified delimiter and returns the items in a vector.
 */
VRN_CORE_API std::vector<std::string> strSplit(const std::string& str, char delim);
/// @overload
VRN_CORE_API std::vector<std::string> strSplit(const std::string& str, const std::string& delim);

/**
 * Joins a sequence of tokens to a string. The converted tokens
 * are separated by the specified delimiter.
 */
template<typename T>
std::string strJoin(const std::vector<T>& tokens, const std::string& delim);

/// Returns true if \p input ends with \p ending
VRN_CORE_API bool endsWith(const std::string& input, const std::string& ending);

/// Returns true if \p input starts with \p ending
VRN_CORE_API bool startsWith(const std::string& input, const std::string& start);

/// Formats the passed byte size as string using the appropriate unit (bytes/kB/MB/GB).
VRN_CORE_API std::string formatMemorySize(uint64_t bytes);

/// Formats the passed length in mm as string using the appropriate unit (km/m/cm/mm/micron/nm).
VRN_CORE_API std::string formatSpatialLength(float base_mm);

/// Formats the passed spatial position in mm as vector string using the appropriate unit (km/m/cm/mm/micron/nm).
VRN_CORE_API std::string formatSpatialLength(const tgt::vec3& position);

/// Formats the passed time in ms as string using the appropriate unit (min/sec/ms).
VRN_CORE_API std::string formatTime(const size_t time);

/**
 * Returns a copy of \p str where all Windows style newlines (\r\n) have been replaced by Unix style newlines (\n)
 */
VRN_CORE_API std::string convertNewlinesWindowsToUnix(const std::string& str);

/**
 * Returns a copy of \p str where all Unix style newlines (\n) have been replaced by Windows style newlines (\r\n)
 */
VRN_CORE_API std::string convertNewlinesUnixToWindows(const std::string& str);

// ----------------------------------------------------------------------------
// template definitions

template<class T>
std::string genericToString(const T& value) {
    std::ostringstream stream;
    stream << value;
    return stream.str();
}

template<class T>
T genericFromString(const std::string& str) {
    T result;
    std::istringstream stream;
    stream.str(str);
    stream >> result;
    if (stream.fail())
        throw VoreenException("failed to convert string '" + str + "'");
    return result;
}

template<typename T>
std::string strJoin(const std::vector<T>& tokens, const std::string& delim) {
    if (tokens.empty())
        return "";
    std::stringstream s;
    s << tokens.at(0);
    for (size_t i=1; i<tokens.size(); i++)
        s << delim << tokens.at(i);
    return s.str();
}

// ----------------------------------------------------------------------------
// Functions related to conversion from BaseType/Format to compile time types
// (i.e. uint*_t, int*_t, float, double and the tgt::Vector?<T> variants)

// Get the BaseType string (i.e., uint*, int*, float or double) from the given Type
template<typename T>
std::string getBaseTypeFromType();

// Get the Format string (i.e., uint*, int*, float, double or Vector*(...)) from the given Type
template<typename T>
std::string getFormatFromType();

// Get the Format string from the given base type and number of channels
std::string getFormatFromBaseTypeAndChannels(const std::string baseType, size_t numChannels);

// Macro: DISPATCH_FOR_FORMAT(format, FUNC, ARGUMENTS...)

// Dispatch to a template function matching the given format with the above macro.
//
// Example usage:
//
// template<typename Voxel>
// static void createVolumeRAM(tgt::svec3 dimensions, std::unique_ptr<VolumeRAM>& res) {
//     res.reset(new VolumeAtomic<Voxel>(dimensions));
// }
//
// std::unique_ptr<VolumeRAM> createVolumeRAMForFormat(std::string format, tgt::svec3 dimensions) {
//     std::unique_ptr<VolumeRAM> res(nullptr);
//     DISPATCH_FOR_FORMAT(format, createVolumeRAM, dimensions, res);
//     return res;
// }

/// Implementation -------------------------------------------------------------

#define IMPL_FORMAT_CONVERSION_FOR_TYPE(TYPE) \
template<> \
inline std::string getBaseTypeFromType<TYPE ## _t>() { \
    return #TYPE; \
}\
template<> \
inline std::string getBaseTypeFromType<tgt::Vector2<TYPE ## _t>>() { \
    return #TYPE; \
}\
template<> \
inline std::string getBaseTypeFromType<tgt::Vector3<TYPE ## _t>>() { \
    return #TYPE; \
}\
template<> \
inline std::string getBaseTypeFromType<tgt::Vector4<TYPE ## _t>>() { \
    return #TYPE; \
}\
\
template<> \
inline std::string getFormatFromType<TYPE ## _t>() { \
    return #TYPE; \
}\
template<> \
inline std::string getFormatFromType<tgt::Vector2<TYPE ## _t>>() { \
    return "Vector2(" #TYPE ")"; \
}\
template<> \
inline std::string getFormatFromType<tgt::Vector3<TYPE ## _t>>() { \
    return "Vector3(" #TYPE ")"; \
}\
template<> \
inline std::string getFormatFromType<tgt::Vector4<TYPE ## _t>>() { \
    return "Vector4(" #TYPE ")"; \
}\

IMPL_FORMAT_CONVERSION_FOR_TYPE(uint8)
IMPL_FORMAT_CONVERSION_FOR_TYPE(uint16)
IMPL_FORMAT_CONVERSION_FOR_TYPE(uint32)
IMPL_FORMAT_CONVERSION_FOR_TYPE(uint64)
IMPL_FORMAT_CONVERSION_FOR_TYPE(int8)
IMPL_FORMAT_CONVERSION_FOR_TYPE(int16)
IMPL_FORMAT_CONVERSION_FOR_TYPE(int32)
IMPL_FORMAT_CONVERSION_FOR_TYPE(int64)
IMPL_FORMAT_CONVERSION_FOR_TYPE(float)
IMPL_FORMAT_CONVERSION_FOR_TYPE(double)

#define IMPL_DISPATCH_FOR_FORMAT_OPTION(TYPE, format, FUNC, ...) \
    if ((format) == getFormatFromType<TYPE ## _t>()) { \
        FUNC<TYPE ## _t>(__VA_ARGS__); \
    } else \
    if ((format) == getFormatFromType<tgt::Vector2<TYPE ## _t>>() ) { \
        FUNC<tgt::Vector2<TYPE ## _t>>(__VA_ARGS__); \
    } else \
    if ((format) == getFormatFromType<tgt::Vector3<TYPE ## _t>>()) { \
        FUNC<tgt::Vector3<TYPE ## _t>>(__VA_ARGS__); \
    } else \
    if ((format) == getFormatFromType<tgt::Vector4<TYPE ## _t>>()) { \
        FUNC<tgt::Vector4<TYPE ## _t>>(__VA_ARGS__); \
    } else \

#define DISPATCH_FOR_FORMAT(format, FUNC, ...) \
    IMPL_DISPATCH_FOR_FORMAT_OPTION(uint8, format, FUNC, __VA_ARGS__) \
    IMPL_DISPATCH_FOR_FORMAT_OPTION(uint16, format, FUNC, __VA_ARGS__) \
    IMPL_DISPATCH_FOR_FORMAT_OPTION(uint32, format, FUNC, __VA_ARGS__) \
    IMPL_DISPATCH_FOR_FORMAT_OPTION(uint64, format, FUNC, __VA_ARGS__) \
    IMPL_DISPATCH_FOR_FORMAT_OPTION(int8, format, FUNC, __VA_ARGS__) \
    IMPL_DISPATCH_FOR_FORMAT_OPTION(int16, format, FUNC, __VA_ARGS__) \
    IMPL_DISPATCH_FOR_FORMAT_OPTION(int32, format, FUNC, __VA_ARGS__) \
    IMPL_DISPATCH_FOR_FORMAT_OPTION(int64, format, FUNC, __VA_ARGS__) \
    IMPL_DISPATCH_FOR_FORMAT_OPTION(float, format, FUNC, __VA_ARGS__) \
    IMPL_DISPATCH_FOR_FORMAT_OPTION(double, format, FUNC, __VA_ARGS__) \
    { \
        tgtAssert(false, "Invalid format"); \
    } \

} // namespace

#endif // VRN_STRINGUTILS_H
