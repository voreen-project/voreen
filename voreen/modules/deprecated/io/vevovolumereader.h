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

#ifndef VRN_VEVOVOLUMEREADER_H
#define VRN_VEVOVOLUMEREADER_H

#include "voreen/core/io/volumereader.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "tgt/logmanager.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

#include <string>
#include <vector>
#include <iomanip>

namespace voreen {

enum ECGPeak {
    PEAK1_START,
    PEAK1,
    NEG_PEAK,
    PEAK2,
    PEAK2_END
};

struct ECGframe {
    int16_t val;
    double tstamp;
    ECGPeak peakType;
    size_t fileid;

    ECGframe();
    ECGframe(int v, double t, ECGPeak p, size_t f)
        : val(v), tstamp(t), peakType(p), fileid(f) {}
};

struct VevoInfoFile {
    size_t num_frames;              // num frames in raw data .rdb file
    size_t dataXRes;                // x resolution
    float  dataXSize;               // x frame size
    size_t dataYRes;                // y resolution
    float  dataYSize;               // y frame size
    size_t ECGDataSize;             // data size of ECG region
    size_t ECGDataOffset;           // offset of ECG raw frames in .rdb file
    double ECGDataTstamp;           // timestamp of the ECG samples
    size_t ECGFreq;
    float  HPlen;                   // length of heart period
};

class Volume;

class VevoFrame {
public:
    int num_;
    size_t size_;
    unsigned int offset_;
    float tstamp_;          // frame's timestamp in msecs
    std::string file_;      // .rdb frame container
    size_t segID_;          // segmented id of the frame depending on the heart period segmentation granularity n;
                            // its range is in [0:n-1].
    size_t fileid_;         // in which consecutive file is the frame located
    float hpBegin_;         // begin of the heart period R-Wave as a timestamp value
    float hpEnd_;           // end of the heart period R-Wave as a timestamp value

    unsigned short int dir_;//  the direction the transducer arm within the RMV
                            //  scanhead was moving when this frame of B-Mode
                            //  image data was acquired (1 = "forward",
                            //  2 = "reverse")
    int16_t ecgval_;

    // init an invalid frame per default
    VevoFrame();
    VevoFrame(int n, size_t s, int o, double t, std::string file, size_t fileid, unsigned short int dir, int16_t ecg = 0);
    VevoFrame(const VevoFrame& v);
    ~VevoFrame() {}

    VevoFrame& operator=(const VevoFrame& v);

    friend std::ostream& operator<<(std::ostream& os, const VevoFrame& v) {
        return os
        << "++++++++++++++++++++\n"
        << "Frame #: " << v.num_ << "\n"
        << "Offset: " << v.offset_ << "\n"
        << "ECG Value: " << v.ecgval_ << "\n"
        << "Size: " << v.size_ << "\n"
        << "Begin of HP R-Wave: " << std::setprecision(10) << v.hpBegin_ << "\n"
        << "Timestamp: "          << std::setprecision(10) << v.tstamp_ << "\n"
        << "End of HP R-Wave: "   << std::setprecision(10) << v.hpEnd_ << "\n"
        << "Contained in: " << v.file_ << "\n"
        << "with File ID: " << v.fileid_ << "\n"
        << "RMV Scanhead Direction: " << v.dir_ << "\n"
        << "====================\n";
    }
};

/**
 * An image data reader class that handles loading ultrasound image data acquired
 * from the VisualSonics Vevo Ultrasound Imaging System. For more info see
 * http://www.visualsonics.com or check the local manuals.
 */
class VevoVolumeReader : public VolumeReader {
    bool loadedValidFile_;
    VolumeRAM* volDS_;
    float zspacing_;
    bool syncDataset_;

    std::string renameSingle(std::string fname);

    // splits a line into its parts at the char * split points
    std::vector<std::string> splitLine(std::string, const char *);

    // delete all occurences of a string c in the string s
    std::string& deleteAll(std::string& s, const std::string& c);

public:
    VevoVolumeReader();

    virtual VolumeList* read(const std::string &fileName);

    struct VevoInfoFile* readCSVFile(const std::string& fname);

    void readVevoFramesCollection(std::string fname,
            std::vector<std::string>& filenames);

    void loadFramesDescrFromFile(const std::string& fname, size_t fileid,
            std::vector<VevoFrame>& frames);

    void loadECGValues(std::string f, size_t ecgds, size_t ecgoffset,
        std::vector<int16_t>& ecgvals);

    size_t loadRawFrame(VevoFrame& frame, char *to);

    size_t loadInvertedRawFrame(VevoFrame& frame, char *to, size_t xres);

    size_t appendToVolume(VevoFrame& v, VolumeRAM_UInt16* ds, size_t offset, size_t xres);

    bool isUSVolLoaded() const { return loadedValidFile_; }

    void setDataset(VolumeRAM* ds) { volDS_ = ds; }

    float getZSpacing() { return zspacing_; }

    bool getSyncDataset() { return syncDataset_; }

    ~VevoVolumeReader() {
    }

private:
    static const std::string loggerCat_;
};

} // namespace voreen

#endif //VRN_VEVOVOLUMEREADER_H
