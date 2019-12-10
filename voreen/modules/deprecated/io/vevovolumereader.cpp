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

#include "vevovolumereader.h"

#include <fstream>
#include <iostream>
#include <string>

#include "tgt/assert.h"
#include "tgt/exception.h"
#include "tgt/logmanager.h"
#include "tgt/vector.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"

using std::string;
using tgt::vec3;
using tgt::ivec3;

namespace voreen {

namespace {

/*
 * A class for signaling bad conversions of types. For more info see:
 * http://www.parashift.com/c++-faq-lite/misc-technical-issues.html#faq-39.2
 */
class BadConversion : public std::runtime_error {
public:
    BadConversion(const std::string& s)
        : std::runtime_error(s)
    { }
};

inline double toDouble(const std::string& s) {
    std::istringstream i(s);
    double x;
    if (!(i >> x))
        throw BadConversion("toDouble(\"" + s + "\")");
    return x;
}

inline int toInt(const std::string& s) {
    std::istringstream i(s);
    int x;
    if (!(i >> x))
        throw BadConversion("toInt(\"" + s + "\")");
    return x;
}

} // namespace

//-------------------------------------------------------------------------------------------------

VevoFrame::VevoFrame() {
    num_ = -1;
}

VevoFrame::VevoFrame(int n, size_t s, int o, double t, std::string file, size_t fileid, unsigned short int dir, int16_t ecg)
    : num_(n)
    , size_(s)
    , offset_(o)
    , tstamp_(static_cast<float>(t))
    , file_(file)
    , fileid_(fileid)
    , dir_(dir)
    , ecgval_(ecg)
{}

VevoFrame::VevoFrame(const VevoFrame& v)
    : num_(v.num_)
    , size_(v.size_)
    , offset_(v.offset_)
    , tstamp_(v.tstamp_)
    , file_(v.file_)
    , segID_(v.segID_)
    , fileid_(v.fileid_)
    , hpBegin_(v.hpBegin_)
    , hpEnd_(v.hpEnd_)
    , dir_(v.dir_)
    , ecgval_(v.ecgval_)
{}

VevoFrame& VevoFrame::operator=(const VevoFrame& v) {
    this->~VevoFrame();
    new(this) VevoFrame(v);
    return *this;
}

//-------------------------------------------------------------------------------------------------

const std::string VevoVolumeReader::loggerCat_ = "voreen.io.VolumeReader.Vevo";

VevoVolumeReader::VevoVolumeReader()
  : loadedValidFile_(false),
    zspacing_(1.0),
    syncDataset_(false)
{
    extensions_.push_back("rdi");
    extensions_.push_back("rdm");
}

VolumeList* VevoVolumeReader::read(const string &fname) {
    tgt::ivec3 dimensions(0);
    tgt::vec3 spacing;
    string fileName;
    size_t posRDI;
    size_t posRDM;

    posRDI = fname.rfind("rdi");
    posRDM = fname.rfind("rdm");
    if (posRDI != string::npos) {
        LERROR("Opening single 2D frames stack, no volume dataset can be"
               "reconstructed from one file.");
        exit(EXIT_FAILURE);
    }
    else if (posRDM != string::npos) {
        std::vector<string> filenames;

        readVevoFramesCollection(fname, filenames);
        volDS_ = new VolumeRAM_UInt8(tgt::vec3(0, 0, 0));
    }
    else {
        LERROR("unknown filename " << fname);
    }

    if (dimensions == tgt::ivec3::zero) {
        //FIXME: this is alway true, something missing?
        LERROR("Wrong dimensions acquired from Vevo info file: " << fname);
        volDS_ = new VolumeRAM_UInt8(tgt::vec3(8));
    }

    VolumeList* volumeList = new VolumeList();
    Volume* volumeHandle = new Volume(volDS_, spacing, tgt::vec3(0.0f));
    volumeHandle->setRealWorldMapping(RealWorldMapping::createDenormalizingMapping(volumeHandle->getBaseType()));
    oldVolumePosition(volumeHandle);
    volumeList->add(volumeHandle);

    return volumeList;
}

std::vector<string> VevoVolumeReader::splitLine(string s, const char *ch) {
    std::vector<string> res;
    size_t offset = 0;
    std::size_t idx = 0;

    while ((idx = s.find(ch, offset)) != string::npos) {
        string s_tmp = s.substr(offset, idx - offset);
        // remove quotation marks
        string s2_tmp = deleteAll(s_tmp, "\"");
        res.push_back(s2_tmp);
        offset = idx + 1;
    }

    string last = s.substr(offset);
    res.push_back(deleteAll(last, "\""));

    return res;
}

// delete all occurences of a string c in the string s.
string& VevoVolumeReader::deleteAll(string& s, const string& c) {
    size_t idx = 0;

    while ((idx = s.find(c, idx)) != string::npos) {
        s.erase(idx, c.size());
    }

    return s;
}

void VevoVolumeReader::loadFramesDescrFromFile(const string& fname,
        size_t f_id, std::vector<VevoFrame>& frames) {

    std::ifstream in(fname.c_str());

    if (!in) {
        LERROR("Cannot open .rdi input file :" << fname);
        return;
    }

    string s;

    while (getline(in, s)) {
        // remove \n
        s.erase(s.end() - 1);
        std::vector<string> v = splitLine(s, ",");

        if (!v.empty()) {
            size_t frameNum = 0, frameOffset = 0, frameSize = 0;
            double frameTstamp = 0.0;

            // check the order of the different attributes in the .rdi file
            // before adding a frame to the framesPerFile_ vector
            if (v[0].substr(0, 17) == "Image Data Offset") {
                frameOffset = atoi(v[1].c_str());
                frameNum = atoi((v[0].substr(26)).c_str());
            }

            if (v[0].substr(0, 15) == "Image Data Size") {
                frameSize = atoi(v[1].c_str());
            }

            if (v[0].substr(0, 20) == "Image Data Timestamp") {
                frameTstamp = atof(v[1].c_str());
            }

            if (v[0].substr(0, 26) == "Image Data Probe Direction") {
                unsigned short f_direction = atoi(v[1].c_str());
                // we are now ready to add a frame:
                string tmp(fname);
                tmp.replace(fname.length() - 3, 3, "rdb");
                VevoFrame v(static_cast<int>(frameNum), frameSize, static_cast<int>(frameOffset), frameTstamp, tmp, f_id, f_direction);
                frames.push_back(v);
            }
        }
    }

    in.close();
}

/*
 * Caller has to cleanup the *VevoInfoFile allocated memory.
 */
struct VevoInfoFile* VevoVolumeReader::readCSVFile(const string& fname) {

    std::ifstream in(fname.c_str());

    if (!in) {
        LERROR("Cannot open rdi input file :" << fname);
        return NULL;
    }

    struct VevoInfoFile *vif = new VevoInfoFile();

    string s;

    while (getline(in, s)) {
        // remove \n
        s.erase(s.end() - 1);
        std::vector<string> v = splitLine(s, ",");

        if (!v.empty()) {
            if (v[0].substr(0, 12) == "Image Length") {
                vif->num_frames = atoi(v[1].c_str());
            }

            if (v[0].substr(0, 23) == "B-Mode/RX/AD-Gate-Width") {
                vif->dataXRes = atoi(v[1].c_str());
            }

            if (v[0].substr(0, 24) == "B-Mode/TX/Trig-Tbl-Trigs") {
                vif->dataYRes = atoi(v[1].c_str());
            }

            if (v[0].substr(0, 13) == "ECG Data Size") {
                vif->ECGDataSize = atoi(v[1].c_str());
            }

            if (v[0].substr(0, 15) == "ECG Data Offset") {
                vif->ECGDataOffset = atoi(v[1].c_str());
            }

            if (v[0].substr(0, 18) == "ECG Data Timestamp") {
                vif->ECGDataTstamp = atof(v[1].c_str());
            }

            if (v[0].substr(0, 20) == "B-Mode/ECG/Frequency") {
                vif->ECGFreq = atoi(v[1].c_str());
            }

            if (v[0].substr(0, 23) == "B-Mode/ECG/Heart-Period") {
                vif->HPlen = static_cast<float>(atof(v[1].c_str()));
            }

            if (v[0].substr(0, 29) == "B-Mode/RX/Sector-Width-Target") {
                vif->dataXSize = static_cast<float>(atof(v[1].c_str()));
            }

            if (v[0].substr(0, 30) == "B-Mode/RX/Sector-Height-Target") {
                vif->dataYSize = static_cast<float>(atof(v[1].c_str()));
            }
        }
    }

    in.close();

    return vif;
}

void VevoVolumeReader::loadECGValues(string f, size_t ecgds, size_t ecgoffset,
                                     std::vector<int16_t>& ecgvals)
{
    int16_t buf;

    string fname = renameSingle(f);

    std::fstream fin(fname.c_str(), std::ios::in | std::ios::binary);

    if (!fin.good()) {
        LERROR("loadECGValues(): bad file I/O stream\n");
        return ;
    }

    fin.seekg(ecgoffset, std::ios::beg);

    if (ecgds > 0) {
        for (size_t j = 0; j < ecgds / 2; j++) {
            fin.read(reinterpret_cast<char*>(&buf), 2);
            ecgvals.push_back(buf);
        }

        LINFO("Loaded " << ecgvals.size() << " ECG samples.");
    }

    fin.close();
}

/**
 * Read all filenames from
 * @param fname and save them into
 * @param filenames
 */
void VevoVolumeReader::readVevoFramesCollection(string fname,
        std::vector<string>& filenames)
{
    std::ifstream in(fname.c_str());

    if (!in) {
        //LERROR("Cannot open input file: " << fname);
        return ;
    }

#ifdef WIN32
    static const string FILE_SEP("\\");

#else
    static const string FILE_SEP("/");

#endif

    // get filename's absolute path
    string tmp_fname(fname);

    size_t pos = tmp_fname.rfind(FILE_SEP);

    tmp_fname.erase(pos + 1);

    string s;

    while (getline(in, s)) {

        if (s == "LATERAL POSITIONS") {
            break;
        }
        else if (s.substr(0, 9) == "z-spacing") {
            std::vector<string> v = splitLine(s, " ");
            zspacing_ = static_cast<float>(atof(v[1].c_str()));
        }
        else if (s.substr(0, 4) == "sync") {
            syncDataset_ = true;
        }
    }

    while (getline(in, s))
        filenames.push_back(tmp_fname + s);

    LINFO(filenames.size() << " files loaded from RDM collection.");

    in.close();
}

string VevoVolumeReader::renameSingle(string fname) {
    string fileName;
    size_t pos;

    if ((pos = fname.rfind("rdi")) != string::npos) {
        fileName = fname;

        if (pos != string::npos)
            fileName.replace(pos, 3, "rdb");

        return fileName;
    }
    else {
        LERROR("unknown filename " << fname);
        return NULL;
    }
}

/*
 * returns the size of the appended frame
 */
size_t VevoVolumeReader::appendToVolume(VevoFrame& v, VolumeRAM_UInt16* ds,
        size_t ds_offset, size_t xres)
{
    size_t ret;
    char *to = reinterpret_cast<char*>(ds->voxel() + ds_offset);

    if (v.dir_ == 1)
        ret = loadRawFrame(v, to);
    else
        ret = loadInvertedRawFrame(v, to, xres);

    return ret;
}

size_t VevoVolumeReader::loadInvertedRawFrame(VevoFrame& frame, char *to, size_t xres) {
    std::fstream fin((frame.file_).c_str(), std::ios::in | std::ios::binary);

    if (!fin) {
        LERROR("Cannot open input file " << frame.file_);
        return 0;
    }
    else {
        LDEBUG("Will read frame from : " << frame.file_);
    }

    if (frame.offset_) {
        fin.seekg(frame.offset_);
        LDEBUG("Skipped to offset: " << fin.tellg());
        LDEBUG(", frame's offset: " << frame.offset_);

        if (!fin.good())
            LERROR("Error positioning stream at offset!");
    }

    // jump to the last line of the frame
    size_t f_begin = frame.offset_;

    size_t f_end = frame.offset_ + frame.size_;

    int line_incr = static_cast<int>(xres) * 2;

    size_t bytes_read = 0;

    fin.seekg(f_end);

    if (!fin.good())
        LERROR("Error positioning stream at end of frame!");

    // read frame backwards linewise
    while (f_end > f_begin) {
        fin.seekg( -line_incr, std::ios::cur);

        fin.read(to, line_incr);
        to += line_incr;
        bytes_read += static_cast<size_t>(fin.gcount());

        fin.seekg( -line_incr, std::ios::cur);
        f_end -= line_incr;

    }

    if (bytes_read < frame.size_)
        LERROR("Truncated input while reading inverted frame: " << bytes_read);

    return bytes_read;
}

size_t VevoVolumeReader::loadRawFrame(VevoFrame& frame, char *to) {
    std::fstream fin((frame.file_).c_str(), std::ios::in | std::ios::binary);

    if (!fin) {
        LERROR("Cannot open input file " << frame.file_);
        return 0;
    }
    else {
        LDEBUG("Will read frame from : " << frame.file_);
    }

    if (frame.offset_) {
        fin.seekg(frame.offset_);
        LDEBUG("Skipped to offset: " << fin.tellg());
        LDEBUG(", frame's offset: " << frame.offset_);

        if (!fin.good())
            LERROR("Error positioning stream!");
    }

    fin.read(to, frame.size_);

    if (!fin.good()) {
        LERROR("Error reading frame: ");

        if (fin.rdstate() & std::ios_base::failbit)
            LERROR("\tFormat error during reading (failbit)");
        else
            if (fin.rdstate() & std::ios_base::badbit)
                LERROR("\tStream corrupted, i.e. invalid offset (badbit)");
            else
                LERROR("\tStream totally b0rked!");

            fin.clear();
    }

    size_t read = (size_t)fin.gcount();

    if (read < frame.size_)
        LERROR("Truncated input while reading a frame!");

    fin.close();

    return read;
}

} // namespace voreen
