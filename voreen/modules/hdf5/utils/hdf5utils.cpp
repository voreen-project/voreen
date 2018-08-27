/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "hdf5utils.h"

#include "tgt/logmanager.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace {

// Loggercat for hdf5 error stack printing
const std::string loggerCat_("voreen.hdf5.errorstack");

#define LOG LDEBUG

// Print an entry of the hdf5 error stack - just like the default error handling routine, but to LDEBUG
// Based on:: https://www.hdfgroup.org/HDF5/doc/UG/13_ErrorHandling.html
herr_t printErrorDescription(unsigned n, const H5E_error2_t *err_desc, void*)
{
#define MSG_SIZE       64
    char                maj[MSG_SIZE];
    char                min[MSG_SIZE];
    char                cls[MSG_SIZE];
    const std::string indent = "  ";

    LOG(indent + "#" + std::to_string(n) + ": " + err_desc->file_name
                + " line " + std::to_string(err_desc->line)
                + " in " + err_desc->func_name + "(): "
                + err_desc->desc);

    if(H5Eget_class_name(err_desc->cls_id, cls, MSG_SIZE)<0)
        return -1;

    if(H5Eget_msg(err_desc->maj_num, NULL, maj, MSG_SIZE)<0)
        return -1;

    if(H5Eget_msg(err_desc->min_num, NULL, min, MSG_SIZE)<0)
        return -1;

    LOG(indent + indent + "class: " + cls);
    LOG(indent + indent + "major: " + maj);
    LOG(indent + indent + "minor: " + min);

    return 0;
}

// Handler for HDF5 errors
herr_t hdf5ErrorHandler(hid_t, void *)
{
    // Lock library for HDF5 operations
    boost::lock_guard<boost::recursive_mutex> lock(voreen::hdf5libMutex);

    LDEBUG("A HDF5 error was detected.");

    // Walk the error stack and print each entry using our defined function
    // First make a copy of the current stack, because it will be cleared when
    // using other hdf5lib-functions, which is problematic when walking the stack
    // in the next step.
    hid_t stack = H5Eget_current_stack();
    H5Ewalk2(stack, H5E_WALK_UPWARD, printErrorDescription, NULL);
    // We are done: Close the stack copy and release its resources.
    H5Eclose_stack(stack);

    //asm("int3");
    return 0;
}

bool initHDF5Lib() {
    // Lock library for HDF5 operations
    boost::lock_guard<boost::recursive_mutex> lock(voreen::hdf5libMutex);

    //Setup error handling: Only print error stack using LDEBUG
    H5Eset_auto(H5E_DEFAULT, hdf5ErrorHandler, NULL);
    return true;
}

// Initialize HDF5 library before starting voreen
bool init = initHDF5Lib();

} //anonymous namespace


namespace voreen {

boost::recursive_mutex hdf5libMutex;

H5::StrType getUTF8StrType() {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);
    // According to https://www.hdfgroup.org/HDF5/doc/Advanced/UsingUnicode/index.html
    // and the example at https://support.hdfgroup.org/ftp/HDF5/examples/misc-examples/vlstratt.c
    H5::StrType utf8StrType(H5::PredType::C_S1);
    utf8StrType.setCset(H5T_CSET_UTF8);
    utf8StrType.setSize(H5T_VARIABLE);
    return utf8StrType;
}

H5::PredType getPredType(const std::string& baseType) {
    if(baseType == "uint8") {               // Unsigned integer types
        return H5::PredType::STD_U8LE;
    } else if(baseType == "uint16") {
        return H5::PredType::STD_U16LE;
    } else if(baseType == "uint32") {
        return H5::PredType::STD_U32LE;
    } else if(baseType == "uint64") {
        return H5::PredType::STD_U64LE;
    } else if(baseType == "int8") {         // Signed integer types
        return H5::PredType::STD_I8LE;
    } else if(baseType == "int16") {
        return H5::PredType::STD_I16LE;
    } else if(baseType == "int32") {
        return H5::PredType::STD_I32LE;
    } else if(baseType == "int64") {
        return H5::PredType::STD_I64LE;
    } else if(baseType == "float") {        // Floating point types
        return H5::PredType::IEEE_F32LE;
    } else if(baseType == "double") {
        return H5::PredType::IEEE_F64LE;
    } else {
        throw tgt::IOException("Cannot generate HDF5 baseType for volume base type \"" + baseType + "\".");
    }
}

std::string getBaseTypeFromDataType(const H5::DataType& type) {
    if(type == H5::PredType::STD_U8LE) {          // Unsigned integer types
        return "uint8";
    } else if(type == H5::PredType::STD_U16LE) {
        return "uint16";
    } else if(type == H5::PredType::STD_U32LE) {
        return "uint32";
    } else if(type == H5::PredType::STD_U64LE) {
        return "uint64";
    } else if(type == H5::PredType::STD_I8LE) {   // Signed integer types
        return "int8";
    } else if(type == H5::PredType::STD_I16LE) {
        return "int16";
    } else if(type == H5::PredType::STD_I32LE) {
        return "int32";
    } else if(type == H5::PredType::STD_I64LE) {
        return "int64";
    } else if(type == H5::PredType::IEEE_F32LE) { // Floating point types
        return "float";
    } else if(type == H5::PredType::IEEE_F64LE) {
        return "double";
    } else if(type == getUTF8StrType()) {
        return "string";
    } else {
        throw tgt::IOException("Unsupported PredType: ID=" + std::to_string(type.getId()));
    }
}

H5::DataSpace createDataSpace(tgt::svec3 volumeDimensions, size_t numberOfChannels /*= 1*/) {
    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);
    // If there is only one channel, no extra DataSpace dimension is needed
    if(numberOfChannels == 1) {
        hsize_t dimensions[3];
        vec3TgtToHDF5(volumeDimensions, dimensions);
        return H5::DataSpace(3, dimensions);
    } else {
        // Voreen only supports 1 to 4 channels.
        tgtAssert(numberOfChannels > 1 && numberOfChannels <= 4, "Invalid number of volume channels");

        hsize_t dimensions[4];
        // Channel is the fastest changing dimension in the DataSpace
        vec4TgtToHDF5(tgt::svec4(numberOfChannels, volumeDimensions), dimensions);
        return H5::DataSpace(4, dimensions);
    }
}

void writeStringAttribute(
#ifdef H5_STUPID_LOCATION_API_CHANGES
        const H5::H5Object& loc,
#else
        const H5::H5Location& loc,
#endif
        const H5std_string& name, const H5std_string& str) {

    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);
    // Create Attribute
    static const H5::DataSpace dataspace(H5S_SCALAR);
    H5::StrType strDataType = getUTF8StrType();

    H5::Attribute attribute = loc.createAttribute(name, strDataType, dataspace);

    // Write to the attribute
    attribute.write(strDataType, str);
}

H5std_string readStringAttribute(
#ifdef H5_STUPID_LOCATION_API_CHANGES
        const H5::H5Object& loc,
#else
        const H5::H5Location& loc,
#endif
        const H5std_string name) {

    boost::lock_guard<boost::recursive_mutex> lock(hdf5libMutex);
    // Try to open element size attribute.
    H5::Attribute attr = loc.openAttribute(name);
    H5::DataSpace dataSpace = attr.getSpace();

    // Check dimensions
    if(static_cast<hsize_t>(dataSpace.getSimpleExtentNpoints()) != 1) {
        throw H5::AttributeIException("readStringAttribute", "Attribute does not contain exactly one string.");
    }

    // Read attribute
    H5std_string str;
    attr.read(attr.getDataType(), str);
    return str;
}

tgt::svec2 vec2HDF5ToTgt(const hsize_t dimensions[2]) {
    // HDF5 Style indexing (slowest changing first): yx!
    return tgt::svec2(dimensions[1], dimensions[0]);
}
tgt::svec3 vec3HDF5ToTgt(const hsize_t dimensions[3]) {
    // HDF5 Style indexing (slowest changing first): zyx!
    return tgt::svec3(dimensions[2], dimensions[1], dimensions[0]);
}
tgt::svec4 vec4HDF5ToTgt(const hsize_t dimensions[4]) {
    // HDF5 Style indexing (slowest changing first): wzyx!
    return tgt::svec4(dimensions[3], dimensions[2], dimensions[1], dimensions[0]);
}

void vec2TgtToHDF5(const tgt::svec2& vec, hsize_t dimensions[2]) {
    // HDF5 Style indexing (slowest changing first): yx!
    dimensions[0] = vec.y;
    dimensions[1] = vec.x;
}
void vec3TgtToHDF5(const tgt::svec3& vec, hsize_t dimensions[3]) {
    // HDF5 Style indexing (slowest changing first): zyx!
    dimensions[0] = vec.z;
    dimensions[1] = vec.y;
    dimensions[2] = vec.x;
}
void vec4TgtToHDF5(const tgt::svec4& vec, hsize_t dimensions[4]) {
    // HDF5 Style indexing (slowest changing first): wzyx!
    dimensions[0] = vec.w;
    dimensions[1] = vec.z;
    dimensions[2] = vec.y;
    dimensions[3] = vec.x;
}

} // namespace voreen
