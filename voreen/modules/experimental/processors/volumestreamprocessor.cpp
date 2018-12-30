/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "volumestreamprocessor.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/processors/processorwidgetfactory.h"
#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/io/textfilereader.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumegl.h"

#ifdef VRN_WITH_LZO
#include "lzo/lzo1x.h"
#endif

#ifdef __unix__
#include <time.h>
#endif

namespace voreen {

namespace {

uint64_t getTicks() {
    uint64_t ticks;

#ifdef __unix__
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    ticks = now.tv_sec * 1000000 + now.tv_nsec / 1000;
#else
    ticks = 0;
#endif
    return ticks;
}

} // namespace

const std::string VolumeStreamProcessor::loggerCat_("voreen.experimental.VolumeStreamProcessor");

VolumeStreamProcessor::VolumeStreamProcessor()
    : Processor(),
      volume_(0),
      filename_("streamfile", "Stream File", "Select Stream File",
                VoreenApplication::app()->getUserDataPath("volumes"), "Stream Files (*.dat)",
                FileDialogProperty::OPEN_FILE, Processor::INVALID_RESULT),
      step_("step", "Time Step", 0, 0, 1000),
      stream_("stream"),
      outport_(Port::OUTPORT, "volumehandle.volumehandle", "volumehandle.volumehandle", false),
      blockSize_(0),
      lastStep_(-1),
      streamFile_(0),
      rawSize_(0)
{
    addProperty(filename_);
    filename_.onChange(MemberFunctionCallback<VolumeStreamProcessor>(this, &VolumeStreamProcessor::loadStream));

    addProperty(step_);
    step_.onChange(MemberFunctionCallback<VolumeStreamProcessor>(this, &VolumeStreamProcessor::loadStep));

    addProperty(stream_);

    addPort(outport_);
}

VolumeStreamProcessor::~VolumeStreamProcessor() {}

Processor* VolumeStreamProcessor::create() const {
    return new VolumeStreamProcessor();
}

void VolumeStreamProcessor::process() {
}

void VolumeStreamProcessor::initialize() {
    Processor::initialize();

    volume_ = new Volume(0, tgt::vec3(0.f), tgt::vec3(0.f));

    outport_.setData(volume_, false);
    loadStep();
}

void VolumeStreamProcessor::deinitialize() {
    delete streamFile_;

    delete volume_;

    Processor::deinitialize();
}

void VolumeStreamProcessor::loadStep() {
    const VolumeRAM* v = volume_->getRepresentation<VolumeRAM>();
    if (!v)
        return;

    if (streamFileName_.empty())
        loadStepFromFile();
    else
        loadStepFromStream();
}

void VolumeStreamProcessor::loadStepFromStream() {
    int step = step_.get();
    if (step >= static_cast<int>(steps_.size()))
        return;


    if (blockSize_ == 0) {
        const std::vector<BlockInfo>& blocks = steps_.at(step).blocks_;
        streamFile_->seekg(static_cast<std::streamoff>(blocks.at(0).offset_));
        streamFile_->read(reinterpret_cast<char*>(volume_->getWritableRepresentation<VolumeRAM>()->getData()),
                          blocks.at(0).size_);

        invalidate(INVALID_RESULT);
        lastStep_ = step;
        return;
    }

    int isDelta = 0;
    int readstep = step;
    if (step > 0 && step == lastStep_ + 1) {
        isDelta = 1;
    } else if (step > 0 && step == lastStep_ - 1) {
        readstep = step + 1;
        isDelta = -1;
    }


    if (!steps_[step].ready_) {
        std::cout << "step not ready " << step << std::endl;
        uint64_t tick = getTicks();
        readBlocks(readstep);
        uint64_t tick_disk = getTicks()- tick;
        tick = getTicks();
        uncompressBlocks(readstep);
        uint64_t tick_uncompress = getTicks()- tick;

        tick = getTicks();
        includeBlocks(readstep, isDelta);
        uint64_t tick_include = getTicks()- tick;

        std::cout << "timings: disk " << tick_disk << ", decompress " << tick_uncompress << ", include "
        << tick_include << std::endl;
    }

    showStep(readstep);
}

void VolumeStreamProcessor::showStep(int step) {
    std::vector<BlockData>& datas = steps_[step].datas_;
    for (size_t j=0; j < datas.size(); j++) {
        delete[] datas[j].rawData_;
        datas[j].rawData_ = 0;
        delete[] datas[j].data_;
        datas[j].data_ = 0;
    }
    steps_[step].ready_ = false;

    invalidate(INVALID_RESULT);   //FIXME: parameter should not be necessary!?!?!?!?!?!?!!!!!!!!!!!!!!!!!!
    lastStep_ = step;
}

void VolumeStreamProcessor::readBlocks(int step) {
    const std::vector<BlockInfo>& blocks = steps_[step].blocks_;
    std::vector<BlockData>& datas = steps_[step].datas_;

    // load raw block data
    streamFile_->seekg(static_cast<std::streamoff>(blocks[0].offset_));

    int i = 0;
    tgt::ivec3 p;
    for (p.x = 0; p.x < blockCount_.x; p.x++) {
        for (p.y = 0; p.y < blockCount_.y; p.y++) {
            for (p.z = 0; p.z < blockCount_.z; p.z++) {
                if (!streamFile_->good()) {
                    LERROR("Seek to step " << step << " in stream file failed: " << streamFileName_);
                    return;
                }

                char* load = 0;
                datas[i].data_ = new uint16_t[blockSize_ * blockSize_ * blockSize_];
                if (compression_ == UNCOMPRESSED) {
                    datas[i].rawData_ = 0;
                    load = reinterpret_cast<char*>(datas[i].data_);
                } else {
                    datas[i].rawData_ = new char[blocks[i].size_];
                    load = reinterpret_cast<char*>(datas[i].rawData_);
                }

                streamFile_->read(load, blocks[i].size_);
                if (!streamFile_->good()) {
                    LERROR("Read from stream file failed: " << streamFileName_);
                    return;
                }
                datas[i].loaded_ = true;
                if (compression_ == UNCOMPRESSED)
                    datas[i].ready_ = true;

                i++;
            }
        }
    }
}

void VolumeStreamProcessor::uncompressBlocks(int step) {
    const std::vector<BlockInfo>& blocks = steps_[step].blocks_;
    std::vector<BlockData>& datas = steps_[step].datas_;

    // uncompress data if necessary
    if (compression_ == LZO) {
        for (size_t j=0; j < blocks.size(); j++) {
            while (!datas[j].loaded_) {
            }
#ifdef VRN_WITH_LZO
            lzo_bytep in = reinterpret_cast<lzo_bytep>(datas[j].rawData_);
            lzo_bytep in_stop = in + blocks[j].size_;
            lzo_bytep out = reinterpret_cast<lzo_bytep>(datas[j].data_);
            datas[j].finalSize_ = 0;
            while (in < in_stop) {
                uint32_t blocksize_uncompressed = *reinterpret_cast<uint32_t*>(in);
                in += sizeof(blocksize_uncompressed);
                uint32_t blocksize_compressed = *reinterpret_cast<uint32_t*>(in);
                in += sizeof(blocksize_compressed);

                lzo_uint new_len = blocksize_uncompressed;

                if (blocksize_compressed == 0 && blocksize_uncompressed > 0) {
                    // uncompressible block, stored as-is
                    for (size_t i=0; i < blocksize_uncompressed; i++)
                        out[i] = in[i];
                    blocksize_compressed = blocksize_uncompressed;
                } else {
                    lzo1x_decompress(in, blocksize_compressed, out, &new_len, 0);
                }

                out += new_len;
                datas[j].finalSize_ += new_len;
                in += blocksize_compressed;
            }
            datas[j].ready_ = true;
#endif
        }
    }
}

void VolumeStreamProcessor::includeBlocks(int step, int isDelta) {
    const std::vector<BlockInfo>& blocks = steps_[step].blocks_;
    std::vector<BlockData>& datas = steps_[step].datas_;
    VolumeRAM* v = volume_->getWritableRepresentation<VolumeRAM>();

    if (prediction_ != DELTA)
        isDelta = 0;

    // copy blocks into final, working on bit-packing and prediction
    uint16_t* v16 = reinterpret_cast<uint16_t*>(v->getData());
    const int bs = blockSize_;
    const tgt::ivec3 s = blockCount_ * blockSize_;

    int i = 0;
    tgt::ivec3 p;
    for (p.x = 0; p.x < blockCount_.x; p.x++) {
        for (p.y = 0; p.y < blockCount_.y; p.y++) {
            for (p.z = 0; p.z < blockCount_.z; p.z++) {

                while (!datas[i].ready_) {
                }

                uint16_t* tmp = datas[i].data_;
                if (blocks[i].bits_ == 8) {
                    uint16_t ref = blocks[i].referenceValue_;
                    uint8_t* tmp8 = reinterpret_cast<uint8_t*>(tmp);

                    if (isDelta != 0) {
                        for (int z = 0; z < bs; z++) {
                            for (int y = 0; y < bs; y++) {
                                for (int x = 0; x < bs; x++) {
                                    v16[(p.x * bs + x) + (p.y * bs + y) * s.x + (p.z * bs + z) * s.x * s.y]
                                        += ref + (tmp8[x + y * bs + z * bs * bs] * isDelta);
                                }
                            }
                        }
                    } else {
                        for (int z = 0; z < bs; z++) {
                            for (int y = 0; y < bs; y++) {
                                for (int x = 0; x < bs; x++) {
                                    v16[(p.x * bs + x) + (p.y * bs + y) * s.x + (p.z * bs + z) * s.x * s.y]
                                        = ref + tmp8[x + y * bs + z * bs * bs];
                                }
                            }
                        }
                    }
                } else if (blocks[i].bits_ == 16) {
                    if (isDelta != 0) {
                        for (int z = 0; z < bs; z++) {
                            for (int y = 0; y < bs; y++) {
                                for (int x = 0; x < bs; x++) {
                                    v16[(p.x * bs + x) + (p.y * bs + y) * s.x + (p.z * bs + z) * s.x * s.y]
                                        += tmp[x + y * bs + z * bs * bs] * isDelta;
                                }
                            }
                        }
                    } else {
                        for (int z = 0; z < bs; z++) {
                            for (int y = 0; y < bs; y++) {
                                for (int x = 0; x < bs; x++) {
                                    v16[(p.x * bs + x) + (p.y * bs + y) * s.x + (p.z * bs + z) * s.x * s.y]
                                        = tmp[x + y * bs + z * bs * bs];
                                }
                            }
                        }
                    }
                } else if (blocks[i].bits_ == 4) {
                    uint16_t ref = blocks[i].referenceValue_;
                    uint8_t* tmp8 = reinterpret_cast<uint8_t*>(tmp);
                    if (isDelta != 0) {
                        for (int z = 0; z < bs; z++) {
                            for (int y = 0; y < bs; y++) {
                                for (int x = 0; x < bs / 2; x++) {
                                    v16[(p.x * bs + x*2) + (p.y * bs + y) * s.x + (p.z * bs + z) * s.x * s.y]
                                        += ref + ((tmp8[x + y * bs/2 + z * bs/2 * bs] >> 4) * isDelta);
                                    v16[(p.x * bs + x*2+1) + (p.y * bs + y) * s.x + (p.z * bs + z) * s.x * s.y]
                                        += ref + ((tmp8[x + y * bs/2 + z * bs/2 * bs] & 0x0F) * isDelta);
                                }
                            }
                        }
                    } else {
                        for (int z = 0; z < bs; z++) {
                            for (int y = 0; y < bs; y++) {
                                for (int x = 0; x < bs/2; x++) {
                                    v16[(p.x * bs + x*2) + (p.y * bs + y) * s.x + (p.z * bs + z) * s.x * s.y]
                                        = ref + ((tmp8[x + y * bs/2 + z * bs/2 * bs] >> 4) * isDelta);
                                    v16[(p.x * bs + x*2+1) + (p.y * bs + y) * s.x + (p.z * bs + z) * s.x * s.y]
                                        = ref + ((tmp8[x + y * bs/2 + z * bs/2 * bs] & 0x0F) * isDelta);
                                }
                            }
                        }
                    }
                } else if (blocks[i].bits_ == 0 && blocks[i].referenceValue_ != 0) {
                    uint16_t ref = blocks[i].referenceValue_;
                    if (isDelta != 0) {
                        for (int z = 0; z < bs; z++) {
                            for (int y = 0; y < bs; y++) {
                                for (int x = 0; x < bs; x++) {
                                    v16[(p.x * bs + x) + (p.y * bs + y) * s.x + (p.z * bs + z) * s.x * s.y]
                                        += ref * isDelta;
                                }
                            }
                        }
                    } else {
                        for (int z = 0; z < bs; z++) {
                            for (int y = 0; y < bs; y++) {
                                for (int x = 0; x < bs; x++) {
                                    v16[(p.x * bs + x) + (p.y * bs + y) * s.x + (p.z * bs + z) * s.x * s.y]
                                        = ref;
                                }
                            }
                        }
                    }
                }

                i++;
            }
        }
    }
    steps_[step].ready_ = true;
}

void VolumeStreamProcessor::loadStepFromFile() {
    VolumeRAM* v = volume_->getWritableRepresentation<VolumeRAM>();
    int step = step_.get();

    std::string filename;

    if (files_.size() > 0 && step < static_cast<int>(files_.size()))
        filename = files_[step];
    else
        return;

    if (step > 0 && step == lastStep_ + 1 && static_cast<int>(deltas_.size()) >= step) {
        LINFO("using delta " << lastStep_ << " -> " << step);
    }

    std::ifstream f(filename.c_str(), std::ios_base::binary);
    if (!f) {
        LERROR("Could not load from file " << filename);
        return;
    }

    std::cout << "reading" << std::endl;
    if (blockSize_ == 0) {
        f.read(reinterpret_cast<char*>(v->getData()), v->getNumBytes());
    } else {
        // read block
    }

    if (!f.good()) {
        LERROR("Read from file failed: " << filename);
        return;
    }

    VolumeRAM_Float* vf = dynamic_cast<VolumeRAM_Float*>(v);
    if (vf && spreadMin_ != spreadMax_) {
        LINFO("Normalizing float volume with min/max: " << tgt::vec2(vf->min(), vf->max()));
        static float gmin = 1000.f;
        static float gmax = -1000.f;
        if (vf->min() < gmin)
            gmin = vf->min();
        if (vf->max() > gmax)
            gmax = vf->max();
        LINFO("Global min/max so far: " << tgt::vec2(gmin, gmax));

        const size_t n = vf->getNumVoxels();

        // use spread values  if available
        if (spreadMin_ != spreadMax_) {
            LINFO("Using spread " << tgt::vec2(spreadMin_, spreadMax_));

            const float d = spreadMax_ - spreadMin_;
            float* voxel = vf->voxel();
            for (size_t i = 0; i < n; ++i)
                voxel[i] = (voxel[i] - spreadMin_) / d;
        } else {
            const float d = vf->max() - vf->min();
            const float p = vf->min();
            float* voxel = vf->voxel();
            for (size_t i = 0; i < n; ++i)
                voxel[i] = (voxel[i] - p) / d;
        }
        vf->invalidate();
    }
    std::cout << "finished loading" << std::endl;
    invalidate(INVALID_RESULT);   //FIXME: parameter should not be necessary!?!?!?!?!?!?!!!!!!!!!!!!!!!!!!

    lastStep_ = step;
}

void VolumeStreamProcessor::loadStream() {
#ifdef VRN_WITH_LZO
    lzo_init();
#endif

    try {
        lastStep_ = -1;

        TextFileReader reader(filename_.get());
        if (!reader)
            throw tgt::FileNotFoundException("reading dat file", filename_.get());

        std::string type;
        std::istringstream args;
        std::string mode, format, compression;
        int steps = 0;
        tgt::ivec3 resolution(0);
        tgt::vec3 sliceThickness(1.f);
        spreadMin_ = 0.f;
        spreadMax_ = 0.f;
        blockSize_ = 0;
        compression_ = UNCOMPRESSED;
        prediction_ = DIRECT;

        steps_.clear();
        files_.clear();
        deltas_.clear();
        streamFileName_ = "";
        std::string blockfile = "";

        while (reader.getNextLine(type, args, false)) {
            if (type == "StreamFileName:") {
                args >> streamFileName_;
                LDEBUG(type << " " << streamFileName_);
                streamFileName_ = tgt::FileSystem::dirName(filename_.get()) + "/" + streamFileName_;
            } else if (type == "BlockFileName:") {
                args >> blockfile;
                LDEBUG(type << " " << blockfile);
                blockfile = tgt::FileSystem::dirName(filename_.get()) + "/" + blockfile;
            } else if (type == "ObjectFileName:") {
                std::string s;
                args >> s;
                LDEBUG(type << " " << s);
                files_.push_back(tgt::FileSystem::dirName(filename_.get()) + "/" + s);
            } else if (type == "Delta:") {
                std::string s;
                args >> s;
                LDEBUG(type << " " << s);
                deltas_.push_back(tgt::FileSystem::dirName(filename_.get()) + "/" + s);
            } else if (type == "Resolution:") {
                args >> resolution[0];
                args >> resolution[1];
                args >> resolution[2];
                LDEBUG(type << " " << resolution[0] << " x " <<
                       resolution[1] << " x " << resolution[2]);
            } else if (type == "SliceThickness:") {
                args >> sliceThickness[0] >> sliceThickness[1] >> sliceThickness[2];
                LDEBUG(type << " " << sliceThickness[0] << " "
                       << sliceThickness[1] << " " << sliceThickness[2]);
            } else if (type == "Blocksize:") {
                args >> blockSize_;
                LDEBUG(type << " " << blockSize_);
            } else if (type == "Mode:") {
                args >> mode;
                LDEBUG(type << " " << mode);
                if (mode == "delta" || mode == "deltac")
                    prediction_ = DELTA;
            } else if (type == "Format:") {
                args >> format;
                LDEBUG(type << " " << format);
            } else if (type == "Compression:") {
                args >> compression;
                LDEBUG(type << " " << compression);

                if (compression == "lzo" || compression == "lzohi")
                    compression_ = LZO;
            } else if (type == "Steps:") {
                args >> steps;
                LDEBUG(type << " " << steps);
            } else if (type == "Blockcount:") {
                args >> blockCount_.x >> blockCount_.y >> blockCount_.z;
                LDEBUG(type << " " << blockCount_);
            } else if (type == "Spread:") {
                args >> spreadMin_ >> spreadMax_;
                LDEBUG(type << " " << spreadMin_ << " " << spreadMax_);
            } else {
                LERROR("Unknown type: " << type);
            }

            if (args.fail()) {
                LERROR("Format error");
            }
        }

        if (!streamFileName_.empty() && blockfile.empty()) {
            LERROR("missing blockfile!");
        } else {
            std::ifstream f(blockfile.c_str(), std::ios_base::binary);
            if (!f) {
                LERROR("Could not load from block file '" << blockfile << "'");
                return;
            }
            for (int step=0; step < steps; step++) {
                StepInfo stepinfo;
                for (int block=0; block < tgt::hmul(blockCount_); block++) {
                    BlockInfo blockinfo;
                    f.read(reinterpret_cast<char*>(&blockinfo), sizeof(blockinfo));
                    stepinfo.blocks_.push_back(blockinfo);
                    BlockData data;
                    data.info_ = &stepinfo.blocks_.back();
                    stepinfo.datas_.push_back(data);
                }
                steps_.push_back(stepinfo);
            }
        }

        volume_->removeRepresentation<VolumeGL>();

        if (format == "UCHAR") {
            VolumeRAM_UInt8* vol8 = new VolumeRAM_UInt8(resolution);
            volume_->setVolumeData(vol8);
            volume_->setSpacing(sliceThickness);
        }
        else if (format == "USHORT") {
            VolumeRAM_UInt16* vol16 = new VolumeRAM_UInt16(resolution);
            volume_->setVolumeData(vol16);
            volume_->setSpacing(sliceThickness);
        }
        else if (format == "FLOAT") {
            VolumeRAM_Float* vf = new VolumeRAM_Float(resolution);
            volume_->setVolumeData(vf);
            volume_->setSpacing(sliceThickness);
        }

        delete streamFile_;
        streamFile_ = 0;
        if (!streamFileName_.empty()) {
            streamFile_ = new std::ifstream(streamFileName_.c_str(), std::ios_base::binary);
            if (!streamFile_->good()) {
                LERROR("Could not open stream file " << streamFileName_);
                delete streamFile_;
                streamFile_ = 0;
            }
            streamFile_->seekg(0, std::ios::end);
            rawSize_ = streamFile_->tellg();
            streamFile_->seekg(0);
        }

        step_.set(0);
        step_.setMaxValue(steps - 1);
    }
    catch (tgt::FileException) {
        LERROR("loading failed");
    }
    catch (std::bad_alloc) {
        LERROR("out of memory");
    }

    loadStep();
}

void VolumeStreamProcessor::setStep(int step) {
    step_.set(step);
}

int VolumeStreamProcessor::getStepCount()  const {
    return step_.getMaxValue();
}

uint64_t VolumeStreamProcessor::getRawSize() const {
    return rawSize_;
}

void VolumeStreamProcessor::flush() {
    for (size_t i=0; i < steps_.size(); i++) {
        std::vector<BlockData>& datas = steps_[i].datas_;
        for (size_t j=0; j < datas.size(); j++) {
            delete[] datas[j].rawData_;
            datas[j].rawData_ = 0;
            delete[] datas[j].data_;
            datas[j].data_ = 0;
            datas[j].loaded_ = false;
            datas[j].ready_ = false;
        }
        steps_[i].ready_ = false;
    }
}

int VolumeStreamProcessor::getBlockSize() const {
    return blockSize_;
}



} // namespace
