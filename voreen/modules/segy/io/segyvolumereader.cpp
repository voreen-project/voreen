/*
 *  segyvolumereader.cpp
 *  voreen_core
 *
 *  Created by Aqeel Al-Naser on 24/05/2010.
 *  Copyright 2010 The University of Manchester. All rights reserved.
 *
 */

#include "segyvolumereader.h"

#include "tgt/exception.h"
#include "tgt/filesystem.h"
#include <sys/types.h>

#include "voreen/core/io/progressbar.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumeoperator.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"

#if WIN32
#ifndef __MINGW32__
    #define fseeko _fseeki64
#else
    #define fseeko fseeko64
#endif
#endif

using tgt::ivec3;   // for int
using tgt::vec3;    // for float

namespace voreen {

    const std::string SEGYVolumeReader::loggerCat_ = "voreen.segy.SEGYVolumeReader";

    // ::::::::::::::: Default Values :::::::::::::::
    // >>>>>>>>>> change later if needed >>>>>>>>>>>>
    tgt::vec3 spacing_ = ivec3(1, 1, 1);
    float timeStep_ = -1.0f;
    std::string metaString_ = "";
    std::string unit_ = "";
    // ----------------------------------------------

    // constructor
    SEGYVolumeReader::SEGYVolumeReader(ProgressBar* progress)
    : VolumeReader(progress)
    {
        extensions_.push_back("sgy");
        extensions_.push_back("segy");
    }

    /**********************************************************
     * ::::::: Overriding Functions from VolumeReader ::::::: *
     **********************************************************/

    VolumeList* SEGYVolumeReader::read(const std::string& fileName) {
        return readSlices(fileName, 0, 0);
    } // read


    // >>>>>>> TO DO: update progress per upload step if required >>>>>
    // >>>>>>> TO DO: change spacing if required >>>>>>>>>>>>>>>>>>>>>>
    // >>>>>>> TO DO: find out why converting gives better results? >>>
    // >>>>>>> TO DO: check if slice order need change >>>>>>>>>>>>>>>>
    // >>>>>>> assumming identity matrix for transformation >>>>>>>>>>>
    // --------------------------------------------------------------
    VolumeList* SEGYVolumeReader::readSlices(const std::string& fileName, size_t firstSlice, size_t lastSlice) {
        // retrieve header info:
        readHeaderInfo(fileName);

        // check if we have to read only some slices instead of the whole volume.
        if ( ! (firstSlice==0 && lastSlice==0)) {
            if (lastSlice > firstSlice) {
                dimensions_.z = static_cast<int>(lastSlice - firstSlice);
            }
        }//if

        std::string info = "Loading SEG-Y file " + fileName + " ";

        if (dimensions_ == tgt::ivec3::zero) {
            throw tgt::CorruptedFileException("No readHints set.", fileName);
        }

        FILE* fin;
        fin = fopen(fileName.c_str(),"rb");

        if (fin==NULL) {
            throw tgt::IOException("Unable to open raw file for reading", fileName);
        }


        // --------------------------------------------------------------------

        // Now create proper volume type:

        VolumeRAM* volume;

        // >>>>>>>>>> assuming default spacing, i.e. 1.0 >>>>>>>>>>>>

        if (dataSampleFormat_ == SEGY_IBM_FLOAT ||
            dataSampleFormat_ == SEGY_IEEE_FLOAT)
        {
            LINFO(info << "(4-byte float)");
            volume = new VolumeRAM_Float(dimensions_);
        }
        else if (dataSampleFormat_ == SEGY_INT32)
        {
            LINFO(info << "(4-byte int)");
            volume = new VolumeRAM_Int32(dimensions_);
        }
        else if (dataSampleFormat_ == SEGY_INT16)
        {
            LINFO(info << "(2-byte int)");
            volume = new VolumeRAM_Int16(dimensions_);
        }
        else if (dataSampleFormat_ == SEGY_INT8)
        {
            LINFO(info << "(1-byte int)");
            volume = new VolumeRAM_Int8(dimensions_);
        }
        else {
            fclose(fin);
            throw tgt::CorruptedFileException("Unsupported format code # " + std::to_string(dataSampleFormat_), fileName);
        }


        // --------------------------------------------------------------------

        // Now upload data from file into volume


        // clear volume to be ready for uploading data from file
        volume->clear();

        if (getProgressBar()) {
            getProgressBar()->setTitle("Loading volume");
            getProgressBar()->setProgressMessage("Loading volume: " + tgt::FileSystem::fileName(fileName));
        }

        // >>>>>>>>>>>>> TO DO: update progress per upload step >>>>>>>>>>>

        size_t offset;

        // Reposition stream position indicator to first slice:
        offset = SEGY_TEXTUAL_HEADER_SIZE + SEGY_BINARY_HEADER_SIZE
                    + extendedTextualFileHeaderRecords_ * SEGY_TEXTUAL_HEADER_SIZE
                    + (SEGY_TRACE_HEADER_SIZE + samplesPerDataTrace_ * sizeOfSample_) // i.e. dimension_.x
                        * dimensions_.y // i.e. number of cross-lines
                        * firstSlice;

        fseeko (fin, offset, SEEK_SET);

        //std::cout << "First slide .. position of file = " << ftell(fin) << "\n"; // >>>>>>>>>>>>> TESTING >>>>>>>>>>>

        // read each trace and upload it to volume
        for (size_t i = 0; i < (static_cast<size_t>(dimensions_.y) * static_cast<size_t>(dimensions_.z)); i++) {

            //std::cout << "trace num = " << i << "\n"; // >>>>>>>>>>>>> TESTING >>>>>>>>>>>

            // Reposition stream position indicator to trace data - skipping the trace header:
            fseeko ( fin, SEGY_TRACE_HEADER_SIZE, SEEK_CUR);

            //std::cout << "\t\tposition of file = " << ftell(fin) << "\n"; // >>>>>>>>>>>>> TESTING >>>>>>>>>>>

            // read all samples of this trace into volume data
            if (fread(reinterpret_cast<char*>(volume->getData())
                      + samplesPerDataTrace_ * sizeOfSample_ * i,
                      1, samplesPerDataTrace_ * sizeOfSample_, fin) == 0)
                LWARNING("fread() failed");

            // convert each sample value to little-endian based on its type:
            if (dataSampleFormat_ == SEGY_IEEE_FLOAT)
            {
                float* samples = reinterpret_cast<float*>(volume->getData())
                                + samplesPerDataTrace_ * i;

                for (int y=0; y < samplesPerDataTrace_; y++) {
                    endian_swap_32(&samples[y]);
                }//for
            }//if
            else if (dataSampleFormat_ == SEGY_INT32)
            {
                int* samples = reinterpret_cast<int*>(volume->getData())
                                + samplesPerDataTrace_ * i;

                for (int y=0; y < samplesPerDataTrace_; y++) {
                    endian_swap_32(&samples[y]);
                }//for
            }//else if
            else if (dataSampleFormat_ == SEGY_INT16)
            {
                short* samples = reinterpret_cast<short*>(volume->getData())
                                + samplesPerDataTrace_ * i;

                for (int y=0; y < samplesPerDataTrace_; y++) {
                    endian_swap_16(&samples[y]);
                }//for
            }//else if

        }//for

        if (lastSlice == 0) {
            if (feof(fin) ) {
                delete volume;
                // throw exception
                throw tgt::CorruptedFileException("unexpected EOF", fileName);
            }
        }

        fclose(fin);

        // --------------------------------------------------------------------

        // >>>>>>>>>>> slice order ???!!! >>>>>>>>>>>>

        // --------------------------------------------------------------------

        VolumeList* volumeList = new VolumeList();
        Volume* volumeHandle = new Volume(volume, spacing_, vec3(0.0f));

        // encode raw parameters into search string
        std::ostringstream searchStream;
        searchStream << "objectModel=" << "I" << "&";
        searchStream << "format=" << dataSampleFormat_ << "&";
        searchStream << "dim_x=" << dimensions_.x << "&";
        searchStream << "dim_y=" << dimensions_.y << "&";
        searchStream << "dim_z=" << dimensions_.z << "&";
        searchStream << "spacing_x=" << spacing_.x << "&";
        searchStream << "spacing_y=" << spacing_.y << "&";
        searchStream << "spacing_z=" << spacing_.z << "&";
        //searchStream << "headerskip=" << headerskip_;
        volumeHandle->setRealWorldMapping(RealWorldMapping::createDenormalizingMapping(volumeHandle->getBaseType()));
        volumeHandle->setOrigin(VolumeURL("segy", fileName, searchStream.str()));
        oldVolumePosition(volumeHandle);

        volumeList->add(volumeHandle);

        return volumeList;


    } // readSlices

    /**********************************************************
     * ::::::::::::::::::: Helper Methods ::::::::::::::::::: *
     **********************************************************/

    // >>>>>> TO DO: handle exceptions correctly
    // current version assumes the followings:
    // 1. no extended textual header records
    // 2. all traces have same number of samples
    // 3. all in-lines have same number of x-lines
    void SEGYVolumeReader::readHeaderInfo(const std::string& fileName)
    {
        // create and open the file:
        FILE* fin;
        fin = fopen(fileName.c_str(),"rb");

        size_t result; // to be used for fread() return value

        // --------------------------------------------------------------------

        // check for Extended Textual File Header Records
        // >>>>>>>>> for now, assume no extra headers >>>>>>>>>>
        extendedTextualFileHeaderRecords_ = 0;

        // --------------------------------------------------------------------

        // Now, get number of samples per data trace from Binary Header,
        // >>>>>>>>>>>> assuming all traces have same number of samples >>>>>>>>>>

        // Reposition stream position indicator to required byte no within Binary File Header:
        fseeko ( fin, SEGY_SAMPLES_PER_DATA_TRACE_BYTE_NUM - 1, SEEK_SET);

        // copy the data, i.e. number of samples per trace:
        result = fread (&samplesPerDataTrace_,sizeof(samplesPerDataTrace_),1,fin);
        if (result != 1) {
            throw tgt::CorruptedFileException("SEGY: Reading error");
        }

        // convert to little-endian:
        endian_swap_16(&samplesPerDataTrace_);

        // --------------------------------------------------------------------

        // Now get data sample format and update size of sample

        // Reposition stream position indicator to required byte no within Binary File Header:
        fseeko ( fin, SEGY_DATA_SAMPLE_FORMAT_BYTE_NUM - 1, SEEK_SET);

        // copy the data, i.e. number of samples per trace:
        result = fread (&dataSampleFormat_,sizeof(dataSampleFormat_),1,fin);
        if (result != 1) {
            throw tgt::CorruptedFileException("SEGY: Reading error");
        }

        // convert to little-endian:
        endian_swap_16(&dataSampleFormat_);

        // calculate size of sample in bytes:
        if (dataSampleFormat_ == SEGY_IBM_FLOAT ||
            dataSampleFormat_ == SEGY_INT32 ||
            dataSampleFormat_ == SEGY_IEEE_FLOAT)
        {
            sizeOfSample_ = 4; //bytes
        }//if

        else if (dataSampleFormat_ == SEGY_INT16)
        {
            sizeOfSample_ = 2; //bytes
        }//else if

        else if (dataSampleFormat_ == SEGY_INT8)
        {
            sizeOfSample_ = 1; //bytes
        }//else if

        else {
            throw tgt::CorruptedFileException("SEGY: Unsupported format code # " + std::to_string(dataSampleFormat_), fileName);
        }//else


        // --------------------------------------------------------------------

        // Now, get number of in-lines and x-lines by accessing each trace header
        // >>>>>>>>> assuming all in-lines have same number of cross-lines >>>>>>>>>>>

        // trace sequence number within line to be read from file
        unsigned int traceSeqNum = 0; //initial value

        // in-line number
        unsigned int inLine = 0;

        // cross-line number
        unsigned int xLine = 0;

        // Reposition stream position indicator to required byte num. of first trace header:
        fseeko ( fin,
               SEGY_TEXTUAL_HEADER_SIZE + SEGY_BINARY_HEADER_SIZE
               + extendedTextualFileHeaderRecords_ * SEGY_TEXTUAL_HEADER_SIZE
               + SEGY_TRACE_SEQUENCE_NUM_WITHIN_LINE_BYTE_NUM - 1,
               SEEK_SET);

        // go through each trace till end of file
        while (1) {

            // copy the data, i.e. trace sequence num. within line
            result = fread (&traceSeqNum,sizeof(traceSeqNum),1,fin);
            if (result != 1) {
                //fputs ("Reading error",stderr);
                break;
            }

            // convert to little-endian:
            endian_swap(traceSeqNum);

            // calculate number of in-lines; once trace sequence number is 1 then we start a new in-line
            if (traceSeqNum == 1) {
                inLine++;
            }//if

            // calculate (max) number of cross-lines
            if (traceSeqNum > xLine) {
                xLine = traceSeqNum;
            }//if


            // Reposition stream position indicator for next trace:
            fseeko ( fin,
                   SEGY_TRACE_HEADER_SIZE
                   + samplesPerDataTrace_ * sizeOfSample_
                   - sizeof(traceSeqNum),
                   SEEK_CUR);

        }//while


        // --------------------------------------------------------------------

        // Now, update dimension of this volume:
        dimensions_ = ivec3(samplesPerDataTrace_, xLine, inLine);

        // --------------------------------------------------------------------

        // terminate
        fclose (fin);


        // print result:
        // >>>>>>>>>>>>>>> Remove Later >>>>>>>>>>>>>
        //std::cout << "No. of samples per trace: " << samplesPerDataTrace_ << "\n";
        //std::cout << "No. of in-lines: " << inLine << "\n";
        //std::cout << "No. of cross-lines: " << xLine << "\n";

    } //readHeaderInfo



    // ::::::: Converters To/From Big-Endian/Little-Endian :::::::
    // -----------------------------------------------------------

    // convert between big-endian and little-endian values (16-bit unsigned short):
    void SEGYVolumeReader::endian_swap (unsigned short& x)
    {
        x = (x>>8) |
        (x<<8);
    }//endian_swap

    // convert between big-endian and little-endian values (32-bit unsigned int):
    void SEGYVolumeReader::endian_swap (unsigned int& x)
    {
        x = (x>>24) |
        ((x<<8) & 0x00FF0000) |
        ((x>>8) & 0x0000FF00) |
        (x<<24);
    }//endian_swap

    // convert between big-endian and little-endian values (32-bit float):
    void SEGYVolumeReader::endian_swap (float& x)
    {
        // temp char to help is swap
        char temp;

        // char pointer pointing to the same float to allow swapping between bytes
        char* charX = (char*) &x;

        // swap bytes 0 and 3:
        temp = charX[3];
        charX[3] = charX[0];
        charX[0] = temp;

        // swap bytes 1 and 2
        temp = charX[2];
        charX[2] = charX[1];
        charX[1] = temp;

    }//endian_swap

    // convert between big-endian and little-endian of 32-bit memory value:
    void SEGYVolumeReader::endian_swap_32 (void* x)
    {
        // char pointer pointing to same location as x, to allow swapping between bytes
        char* charX = (char*) x;

        // temp char to help is swap
        char temp;

        // swap bytes 0 and 3:
        temp = charX[3];
        charX[3] = charX[0];
        charX[0] = temp;

        // swap bytes 1 and 2
        temp = charX[2];
        charX[2] = charX[1];
        charX[1] = temp;

    }//endian_swap_32

    // convert between big-endian and little-endian of 16-bit memory value:
    void SEGYVolumeReader::endian_swap_16 (void* x)
    {
        // char pointer pointing to same location as x, to allow swapping between bytes
        char* charX = (char*) x;

        // temp char to help is swap
        char temp;

        // swap:
        temp = charX[1];
        charX[1] = charX[0];
        charX[0] = temp;

    }//endian_swap_16


    VolumeReader* SEGYVolumeReader::create(ProgressBar* progress) const {
        return new SEGYVolumeReader(progress);
    }

} // namespace voreen
