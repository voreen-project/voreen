/*
 *  segyvolumereader.h
 *  voreen_core
 *
 *  Created by Aqeel Al-Naser on 23/05/2010.
 *  Copyright 2010 The University of Manchester. All rights reserved.
 *
 */

#ifndef VRN_SEGYVOLUMEREADER_H
#define VRN_SEGYVOLUMEREADER_H

// ::: Size of Headers in SEG-Y file :::
#define SEGY_TEXTUAL_HEADER_SIZE 3200
#define SEGY_BINARY_HEADER_SIZE 400
#define SEGY_TRACE_HEADER_SIZE 240
// -------------------------------------
// ::: Byte no. of some interesting values from SEGY Binary Header :::
// NOTE: address = byte no. - 1
#define SEGY_SAMPLE_INTERVAL_BYTE_NUM 3217            // 2 bytes
#define SEGY_SAMPLES_PER_DATA_TRACE_BYTE_NUM 3221     // 2 bytes
#define SEGY_DATA_SAMPLE_FORMAT_BYTE_NUM 3225         // 2 bytes
#define SEGY_FIXED_LENGTH_TRACE_FLAG_BYTE_NUM 3225    // 2 bytes
#define SEGY_NUM_OF_EXTENDED_TEXTUAL_FILE_HEADER_BYTE_NUM 3505 // 2 bytes
#define SEGY_DATA_SAMPLE_FORMAT_BYTE_NUM 3225         // 2 bytes
// -------------------------------------
// ::: Byte no. of some interesting values from SEGY Trace Header :::
// NOTE: address = address of current trace header + byte no. - 1
#define SEGY_TRACE_SEQUENCE_NUM_WITHIN_LINE_BYTE_NUM 1    // 4 bytes
// -------------------------------------

#include <string>

#include "tgt/vector.h"
#include "tgt/matrix.h"

#include "voreen/core/io/volumereader.h"
#include "voreen/core/datastructures/volume/modality.h"

namespace voreen {
/**
    * Reader for <tt>.segy</tt> file containing a seismic volume dataset.
    */
class SEGYVolumeReader : public VolumeReader {

public:
    SEGYVolumeReader(ProgressBar* progress = 0);
    virtual VolumeReader* create(ProgressBar* progress = 0) const;

    virtual std::string getClassName() const   { return "SEGYVolumeReader"; }
    virtual std::string getFormatDescription() const { return "SEGY seismic data"; }

    virtual VolumeList* read(const std::string& fileName);

    virtual VolumeList* readSlices(const std::string& fileName, size_t firstSlice=0, size_t lastSlice=0);

    // SEGY DATA FORMAT
    enum {
        SEGY_IBM_FLOAT=1,    // 4-byte IBM floating-point
        SEGY_INT32=2,        // 4-byte two'scomplement integer
        SEGY_INT16=3,        // 2-byte two'scomplement integer
        SEGY_IEEE_FLOAT=5,   // 4-byte IEEE floating-point
        SEGY_INT8=8,         // 1-byte two'scomplement integer
    };

private:

    static const std::string loggerCat_;

    tgt::ivec3 dimensions_;
    tgt::vec3 spacing_;
    tgt::mat4 transformation_;
    //float timeStep_; //Currently not used
    Modality modality_;
    std::string metaString_;
    std::string unit_;

    // these are the header info required to be retrieved before any reading attempt:
    short samplesPerDataTrace_;
    short dataSampleFormat_;
    short extendedTextualFileHeaderRecords_;
    size_t sizeOfSample_;

    // retrieves header info for a given SEGY file:
    virtual void readHeaderInfo(const std::string& fileName);


    // ::::::: Converters To/From Big-Endian/Little-Endian :::::::
    // -----------------------------------------------------------
    // convert between big-endian and little-endian values (16-bit unsigned short):
    void endian_swap (unsigned short& x);

    // convert between big-endian and little-endian values (32-bit unsigned int):
    void endian_swap (unsigned int& x);

    // convert between big-endian and little-endian values (32-bit float):
    void endian_swap (float& x);

    // convert between big-endian and little-endian of any 32-bit memory value:
    void endian_swap_32 (void* x);

    // convert between big-endian and little-endian of any 16-bit memory value:
    void endian_swap_16 (void* x);

}; // class SEGYVolumeReader


} // namespace voreen

#endif    // VRN_SEGYVOLUMEREADER_H
