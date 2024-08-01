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

#include "commands_streaming.h"
#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/io/volumeserializerpopulator.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorsubset.h"
#include "voreen/core/datastructures/volume/volumecollection.h"
#include "voreen/core/io/textfilereader.h"
#include "modules/experimental/processors/volumestreamprocessor.h"

#include "tgt/exception.h"

#include <iostream>
#include <limits>
#include <sstream>

#ifdef VRN_WITH_LZO
#include "lzo/lzo1x.h"
#endif

#ifdef __unix__
#include <time.h>
#endif

using std::string;
using tgt::vec3;
using tgt::ivec3;
using tgt::vec2;

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


CommandStreaming::CommandStreaming() /*:
    Command("--stream", "", "Create data stream out of multiple data sets.\n"
            "\tblocksize:\tset to 0 to disable blocking\n"
            "\tprediction:\tdirect: save data directly\n"
            "\t           \tdelta: save difference data\n"
            "\t           \tdeltac: save difference data and bit-compress\n"
            "\tcompression:\tnone, lzo, lzohi\n"
            "\toutformat:\toutput format: float, ushort",
            "<block-size> <prediction> <compression> <outformat> <in> <out>", -1) */
{
//    loggerCat_ += "stream";
}

bool CommandStreaming::checkParameters(const std::vector<std::string>& /*parameters*/) {
    /*if (parameters.size() != 6)
        return false;

    std::istringstream s(parameters[0]);
    int blocksize;
    if (!((s >> blocksize) && (blocksize >= 0)))
        return false;

    std::set<std::string> modes;
    modes.insert("direct");
    modes.insert("delta");
    modes.insert("deltac");
    if (!isValueInSet(parameters[1], modes))
        return false;

    std::set<std::string> compression;
    compression.insert("none");
    compression.insert("lzo");
    compression.insert("lzohi");
    if (!isValueInSet(parameters[2], compression))
        return false;

    std::set<std::string> outformats;
    outformats.insert("float");
    outformats.insert("ushort");
    if (!isValueInSet(parameters[3], outformats))
        return false; */

    return true;
}

namespace {

void analyze(VolumeRAM_Float* vf) {
    const size_t n = vf->getNumVoxels();

    int buckets[65536];
    for (int i=0; i < 65536; i++)
        buckets[i] = 0;

    float* voxel = vf->voxel();
    float sum = 0.f;
    for (size_t i = 0; i < n; ++i) {
        float f = voxel[i];
        sum += f;
        buckets[static_cast<int>(f * 65535.f)]++;
    }

    int used = 0;
    for (int i=0; i < 65536; i++)
        if (buckets[i] > 0)
            used++;

    std::cout << vf->min() << " / " << vf->max() << " : " << (sum / n) << ", "
              << used << " => " << ((used * 100) / 65536) << "%"  << std::endl;
}

int usedBits16(uint32_t value) {
    for (uint32_t i=0; i < 16; i++) {
        if (static_cast<uint32_t>(1 << i) >= value)
            return i;
    }
    return 16;
}

inline void endian_swap(unsigned int& x) {
    x = (x>>24) |
        ((x<<8) & 0x00FF0000) |
        ((x>>8) & 0x0000FF00) |
        (x<<24);
}

} // namespace



bool CommandStreaming::execute(const std::vector<std::string>& /*parameters*/) {
/*    VolumeSerializerPopulator volLoadPop;
    //const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();

    int blocksize;
    std::istringstream s(parameters[0]);
    s >> blocksize;
    const string mode = parameters[1];
    const string compression = parameters[2];
    const string outformat = parameters[3];
    const string in = parameters[4];
    const string out = parameters[5];
    LINFO("got input " << in << ", block size: " << blocksize << ", "
          << "mode: " << mode << ", writing to '" << out << "'");

#ifdef VRN_WITH_LZO
    const size_t LZO_BLOCKSIZE = 1048576;
    lzo_byte lzo_buf[LZO_BLOCKSIZE + LZO_BLOCKSIZE / 16 + 64 + 3];
    lzo_bytep lzo_wrkmem[LZO1X_1_MEM_COMPRESS];
    lzo_bytep lzo_wrkmem_hi[LZO1X_999_MEM_COMPRESS];

    if (compression == "lzo" || compression == "lzohi")
        lzo_init();
#else
    if (compression == "lzo" || compression == "lzohi") {
        LERROR("Voreen was compiled without VRN_WITH_LZO, lzo compression not available.");
        return false;
    }
#endif

    string filename_dat = out + ".dat";
    string filename_stream = out + ".str";
    string filename_blocks = out + ".blk";
    string filename_log = out + ".log";

    std::ofstream fstr(filename_stream.c_str(), std::ios_base::binary);

    std::ofstream flog(filename_log.c_str());
    flog << "params: voltool --stream ";
    for (size_t i=0; i < parameters.size(); i++)
        flog << parameters[i] << " ";
    flog << std::endl;

    uint64_t startTime = getTicks();
    int log_totalblocks = 0;
    int log_incompressible = 0;
    int log_bitblocks[17];
    uint64_t log_bitcompressedsize[17];
    uint64_t bitcompressedsize = 0;
    int log_actualbitblocks[17];
    int log_nonZeroUniformBlocks = 0;

    for (int i=0; i <= 16; i++) {
        log_bitblocks[i] = 0;
        log_bitcompressedsize[i] = 0;
        log_actualbitblocks[i] = 0;
    }

    try {
        TextFileReader reader(in);
        if (!reader)
            throw tgt::FileNotFoundException("reading dat file", in);

        std::string type;
        std::istringstream args;
        std::string format;
        tgt::ivec3 resolution(0);
        tgt::vec3 sliceThickness(1.f);
        tgt::ivec3 resize(0);
        float spreadMin = 0.f;
        float spreadMax = 0.f;
        std::vector<std::string> files;
        std::string byteOrder;

        files.clear();
        std::string streamfile;

        // parse input dat files
        while (reader.getNextLine(type, args, false)) {
            if (type == "ObjectFileName:") {
                std::string s;
                args >> s;
                LDEBUG(type << " " << s);
                files.push_back( tgt::FileSystem::dirName(in) + "/" +  s);
            } else if (type == "Resolution:") {
                args >> resolution[0];
                args >> resolution[1];
                args >> resolution[2];
                LDEBUG(type << " " << resolution[0] << " x " <<
                       resolution[1] << " x " << resolution[2]);
            } else if (type == "SliceThickness:") {
                args >> sliceThickness[0] >> sliceThickness[1] >> sliceThickness[2];
                LDEBUG(type << " " << sliceThickness[0] << " " <<
                       sliceThickness[1] << " " << sliceThickness[2]);
            } else if (type == "Mode:") {
                // args >> mode;
                // LDEBUG(type << " " << mode);
            } else if (type == "Format:") {
                args >> format;
                LDEBUG(type << " " << format);
            } else if (type == "Spread:") {
                args >> spreadMin >> spreadMax;
                LDEBUG(type << " " << spreadMin << " " << spreadMax);
            } else if (type == "Resize:") {
                args >> resize.x >> resize.y >> resize.z;
            } else if (type == "ByteOrder:") {
                args >> byteOrder;
            } else {
                LERROR("Unknown type: " << type);
            }

            if (args.fail()) {
                LERROR("Format error");
            }
        }

        if (tgt::hmul(resize) == 0)
            resize = resolution;

        // allocate input volume
        VolumeRAM_UInt8* vol8 = 0;
        VolumeRAM_Float* vf = 0;
        VolumeRAM* v = 0;
        if (format == "UCHAR") {
            vol8 = new VolumeRAM_UInt8(resolution);
            v = vol8;
        }
        else if (format == "FLOAT") {
            vf = new VolumeRAM_Float(resolution);
            v = vf;
        }

        VolumeRAM_UInt16* vol16 = new VolumeRAM_UInt16(resolution);
        Volume vol16h(vol16, sliceThickness, vec3(0.0f));

        ivec3 size(blocksize);
        ivec3 blockcount = tgt::iceil(vec3(resize) / vec3((float)blocksize));
        if (blocksize == 0) {
            size = resize;
            blockcount = ivec3(1);
        }

        // write output dat file
        std::ofstream fdat(filename_dat.c_str());
        fdat << "Steps:\t\t" << files.size() << std::endl;
        fdat << "Mode:\t\t" << mode << std::endl;
        fdat << "Blocksize:\t" << blocksize << std::endl;
        fdat << "Blockcount:\t" << blockcount.x << " " << blockcount.y << " "
             << blockcount.z << std::endl;
        fdat << "Resolution:\t" << resize.x << " " << resize.y << " " << resize.z << std::endl;
        fdat << "SliceThickness:\t" << sliceThickness.x << " " << sliceThickness.y << " "
             << sliceThickness.z << std::endl;

        if (outformat == "float")
            fdat << "Format:\t\t" << "FLOAT" << std::endl;
        else if (outformat == "ushort")
            fdat << "Format:\t\t" << "USHORT" << std::endl;
        fdat << "Compression:\t" << compression << std::endl;
        fdat << "StreamFileName:\t" << tgt::FileSystem::fileName(filename_stream) << std::endl;
        fdat << "BlockFileName:\t" << tgt::FileSystem::fileName(filename_blocks) << std::endl;


        // open block info file for output
        std::ofstream fblocks(filename_blocks.c_str(), std::ios_base::binary);

        LINFO("creating " << tgt::hmul(blockcount) << " blocks");

        float gmin = 10000.f;
        float gmax = -10000.f;
        long rawsize = 0;

        uint16_t* tmp = new uint16_t[tgt::hmul(size)];

        // work on input files
        std::vector<VolumeRAM_UInt16*> blocks_last;
        for (size_t t=0; t < files.size(); t++) {
            LINFO("file " << files[t]);
            std::ifstream f(files[t].c_str(), std::ios_base::binary);
            if (!f) {
                LERROR("Could not open file " << files[t]);
                return false;
            }

            f.read(reinterpret_cast<char*>(v->getData()), v->getNumBytes());
            if (!f.good()) {
                LERROR("Read from file failed: " << files[t]);
                return false;
            }

            if (byteOrder == "big-endian") {
                LINFO("swapping endianess");
                const size_t n = vf->getNumVoxels();
                float* fvoxel = vf->voxel();
                for (size_t i = 0; i < n; i++) {
                    unsigned int& v = *reinterpret_cast<unsigned int*>(&fvoxel[i]);
                    endian_swap(v);
                }

                // special handling for hurricane dataset
                for (size_t i = 0; i < n; i++) {
                    if (fvoxel[i] >= 1e35)
                        fvoxel[i] = 0.f;
                }
                vf->invalidate();
            }

            // normalize floats to [0.0; 1.0]
            if (vf && spreadMin != spreadMax) {
                LINFO("Normalizing float volume with min/max: " << tgt::vec2(vf->min(), vf->max()));
                if (vf->min() < gmin)
                    gmin = vf->min();
                if (vf->max() > gmax)
                    gmax = vf->max();

                const size_t n = vf->getNumVoxels();

                // use spread values if available
                if (spreadMin != spreadMax) {
                    LINFO("Using spread " << tgt::vec2(spreadMin, spreadMax));

                    const float d = spreadMax - spreadMin;
                    float* voxel = vf->voxel();
                    for (size_t i = 0; i < n; ++i)
                        voxel[i] = (voxel[i] - spreadMin) / d;
                } else {
                    const float d = vf->max() - vf->min();
                    const float p = vf->min();
                    float* voxel = vf->voxel();
                    for (size_t i = 0; i < n; ++i)
                        voxel[i] = (voxel[i] - p) / d;
                }
                vf->invalidate();
            }

            if (vf) {
                const size_t n = vf->getNumVoxels();
                float* fvoxel = vf->voxel();
                uint16_t* ivoxel = vol16->voxel();
                for (size_t i = 0; i < n; ++i)
                    ivoxel[i] = static_cast<uint16_t>(fvoxel[i] * std::numeric_limits<uint16_t>::max());
                vol16->invalidate();
            }


            std::vector<VolumeRAM_UInt16*> blocks;
            VolumeStreamProcessor::BlockInfo blockinfo;

            int blocknum = 0;
            ivec3 p;

            if (blockcount == ivec3(1)) {
                blockinfo.offset_ = fstr.tellp();

                Volume* block = VolumeOperatorSubset::APPLY_OP(&vol16h, tgt::ivec3(0), resize);
                VolumeRAM* v = block->getWritableRepresentation<VolumeRAM>();
                fstr.write(reinterpret_cast<char*>(v->getData()), v->getNumBytes());
                blockinfo.size_ = v->getNumBytes();
                delete block;

                blockinfo.bits_ = 16;
                fblocks.write(reinterpret_cast<char*>(&blockinfo), sizeof(blockinfo));
            }
            else
            for (p.x = 0; p.x < blockcount.x; p.x++) {
                for (p.y = 0; p.y < blockcount.y; p.y++) {
                    for (p.z = 0; p.z < blockcount.z; p.z++) {
                        log_totalblocks++;
                        Volume* h= VolumeOperatorSubset::APPLY_OP(&vol16h, p * size, size);
                        VolumeRAM_UInt16* block = dynamic_cast<VolumeRAM_UInt16*>(h->getWritableRepresentation<VolumeRAM>());
                        blocks.push_back(block);

                        blockinfo.offset_ = fstr.tellp();
                        blockinfo.minValue_ = block->min();
                        blockinfo.maxValue_ = block->max();

                        // std::cout << "min/max: " << tgt::ivec2(block->min(), block->max())
                        //           << " usage: " << block->max() - block->min() << " => "
                        //           << static_cast<float>(block->max() - block->min()) / std::numeric_limits<uint16_t>::max()
                        //           << std::endl;

                        rawsize += block->getNumBytes();
                        size_t blockBytes = block->getNumBytes();

                        uint16_t* data = 0;
                        if (t == 0 || mode == "direct") {
                            blockinfo.prediction_ = VolumeStreamProcessor::DIRECT;

                            const size_t n = block->getNumVoxels();
                            uint16_t* voxel = block->voxel();

                            uint16_t lmax = std::numeric_limits<uint16_t>::min();
                            uint16_t lmin = std::numeric_limits<uint16_t>::max();


                            for (size_t i = 0; i < n; ++i) {
                                uint16_t stmp = voxel[i];
                                if (stmp < lmin)
                                    lmin = stmp;
                                if (stmp > lmax)
                                    lmax = stmp;

                            }

                            // blockinfo.referenceValue_ = lmin;
                            // for (size_t i = 0; i < n; ++i) {
                            //     voxel[i] = voxel[i] - lmin;
                            // }

                            int bits = usedBits16(lmax - lmin + 1);
                            int storebits = 16;
                            blockinfo.bits_ = 16;
                            blockinfo.referenceValue_ = 0;

                            std::cout << "lmin/max:" << tgt::ivec2(lmin, lmax)
                                      << " usage: " << lmax - lmin << " => "
                                      << bits << " bits (store " << storebits << ")" << std::endl;

                            data = block->voxel();
                            bitcompressedsize += block->getNumBytes();
                            for (int i=0; i <= 16; i++)
                                log_bitcompressedsize[i] += block->getNumBytes();
                            log_bitblocks[bits]++;
                            log_actualbitblocks[16]++;
                        } else {
                            blockinfo.prediction_ = VolumeStreamProcessor::DELTA;
                            const size_t n = block->getNumVoxels();
                            uint16_t* voxel = block->voxel();
                            uint16_t* lastvoxel = blocks_last[blocknum]->voxel();

                            int16_t lmax = std::numeric_limits<int16_t>::min();
                            int16_t lmin = std::numeric_limits<int16_t>::max();

                            for (size_t i = 0; i < n; ++i) {
                                tmp[i] = voxel[i] - lastvoxel[i];
                                int16_t stmp = tmp[i];
                                if (stmp < lmin)
                                    lmin = stmp;
                                if (stmp > lmax)
                                    lmax = stmp;
                            }

                            int bits = usedBits16(lmax - lmin + 1);
                            int storebits = 16;

                            log_bitblocks[bits]++;

                            if (bits == 0 && mode == "deltac") {
                                storebits = 0;
                                blockinfo.bits_ = 0;
                                blockinfo.referenceValue_ = lmin;
                                blockBytes = 0;
                                if (lmin != 0)
                                    log_nonZeroUniformBlocks++;
                            }
                            else if (bits <= 4 && mode == "deltac") {
                                storebits = 4;
                                blockinfo.bits_ = 4;
                                blockinfo.referenceValue_ = lmin;

                                uint8_t* data8 = reinterpret_cast<uint8_t*>(tmp);
                                for (size_t i = 0; i < (n / 2); ++i) {
                                    int16_t o1 = tmp[2*i];
                                    int16_t o2 = tmp[2*i+1];
                                    data8[i] =   (static_cast<uint8_t>(o1 - lmin) << 4) |
                                                 (static_cast<uint8_t>(o2 - lmin) & 0x0F);
                                }
                                blockBytes /= 4;
                            } else if (bits <= 8 && mode == "deltac") {
                                storebits = 8;
                                blockinfo.bits_ = 8;
                                blockinfo.referenceValue_ = lmin;

                                uint8_t* data8 = reinterpret_cast<uint8_t*>(tmp);
                                for (size_t i = 0; i < n; ++i) {
                                    data8[i] = tmp[i] - lmin;
                                }
                                blockBytes /= 2;
                            } else {
                                blockinfo.bits_ = 16;
                                blockinfo.referenceValue_ = 0;
                            }

                            log_actualbitblocks[storebits]++;
                            bitcompressedsize += (block->getNumVoxels() * storebits) / 8;

                            for (int i=0; i <= 16; i++)
                                if (bits <= i)
                                    log_bitcompressedsize[i] += (block->getNumVoxels() * i) / 8;
                                else
                                    log_bitcompressedsize[i] += (block->getNumVoxels() * bits) / 8;


                            std::cout << "lmin/max:" << tgt::ivec2(lmin, lmax)
                                      << " usage: " << lmax - lmin << " => "
                                      << bits << " bits (store " << storebits << ")" << std::endl;

                            data = tmp;
                        }

#ifdef VRN_WITH_LZO
                        if (compression == "lzo" || compression == "lzohi") {
                            blockinfo.size_ = 0;

                            lzo_bytep in = (lzo_bytep)data;
                            size_t pos = 0;
                            while (pos < blockBytes) {
                                size_t read = LZO_BLOCKSIZE;
                                if (pos + LZO_BLOCKSIZE > blockBytes)
                                    read = blockBytes - pos;

                                lzo_uint out_len;
                                if (compression == "lzohi")
                                    lzo1x_999_compress(in, read, lzo_buf, &out_len, lzo_wrkmem_hi);
                                else
                                    lzo1x_1_compress(in, read, lzo_buf, &out_len, lzo_wrkmem);

                                if (out_len > read) {
                                    log_incompressible++;
                                    // store uncompressible blocks unchanged
                                    std::cout << "compressed larger than uncompressed, storing as-is! "
                                              << out_len << " > " << read << std::endl;

                                    uint32_t blocksize_compressed = 0;
                                    uint32_t blocksize_uncompressed = read;
                                    fstr.write(reinterpret_cast<char*>(&blocksize_uncompressed), sizeof(blocksize_uncompressed));
                                    fstr.write(reinterpret_cast<char*>(&blocksize_compressed), sizeof(blocksize_compressed));

                                    fstr.write(reinterpret_cast<char*>(in), read);
                                    blockinfo.size_ += read + sizeof(blocksize_compressed) + sizeof(blocksize_uncompressed);

                                } else {
                                    uint32_t blocksize_compressed = out_len;
                                    uint32_t blocksize_uncompressed = read;
                                    fstr.write(reinterpret_cast<char*>(&blocksize_uncompressed), sizeof(blocksize_uncompressed));
                                    fstr.write(reinterpret_cast<char*>(&blocksize_compressed), sizeof(blocksize_compressed));

                                    fstr.write(reinterpret_cast<char*>(lzo_buf), out_len);
                                    blockinfo.size_ += out_len + sizeof(blocksize_compressed) + sizeof(blocksize_uncompressed);
                                }
                                in += read;
                                pos += read;
                            }
                        }
                        else
#endif // VRN_WITH_LZO
                        {
                            fstr.write(reinterpret_cast<char*>(data), blockBytes);
                            blockinfo.size_ = blockBytes;
                        }

                        fblocks.write(reinterpret_cast<char*>(&blockinfo), sizeof(blockinfo));
                        blocknum++;
                    }
                }
            }
            for (size_t i=0; i < blocks_last.size(); i++)
                delete blocks_last[i];
            blocks_last = blocks;
        }
        for (size_t i=0; i < blocks_last.size(); i++)
            delete blocks_last[i];
        delete[] tmp;

        LINFO("Global min/max: " << tgt::vec2(gmin, gmax));

        uint16_t imin = static_cast<uint16_t>(((gmin - spreadMin) / (spreadMax - spreadMin)) * std::numeric_limits<uint16_t>::max());
        uint16_t imax = static_cast<uint16_t>(((gmax - spreadMin) / (spreadMax - spreadMin)) * std::numeric_limits<uint16_t>::max());
        LINFO("    => min/max: " << tgt::ivec2(imin, imax)
              << " usage: " << imax - imin << " => " << static_cast<float>(imax - imin) / std::numeric_limits<uint16_t>::max());

        if (mode == "deltac")
            LINFO("variable word length compression: " << bitcompressedsize
                  << " (" << (1.f - static_cast<float>(bitcompressedsize) / rawsize) << "%)");

        delete v;

        int secs = static_cast<int>((getTicks() - startTime) / 1000000);
        flog << "runtime: " << secs << " sec (" << secs / 60 << " min)" << std::endl;

        uint64_t totalsize = static_cast<uint64_t>(log_totalblocks) * blocksize * blocksize * blocksize * sizeof(uint16_t);
        uint64_t compressed = fstr.tellp();
        flog << "raw size:   " << totalsize << std::endl
             << "compressed: " << compressed << std::endl
             << "ratio:      x" << totalsize / (float)compressed << ", "
             << 100.f - (compressed / (float)totalsize) * 100.f << "%"<< std::endl;

        flog << "spread: " << vec2(spreadMin, spreadMax) << std::endl
             << "min/max: " << tgt::vec2(gmin, gmax) << std::endl
             << "usage: " << imax - imin << " ("
             << 100.f * static_cast<float>(imax - imin) / std::numeric_limits<uint16_t>::max() << "%, "
             << usedBits16(imax - imin + 1) << " bits)" << std::endl;


        flog << "bitcompressed: " << bitcompressedsize << "\t("
             << 100.f - (bitcompressedsize / (float)totalsize) * 100.f << "%)" << std::endl;

        for (int i=16; i >= 0; i--)
            flog << "bitcompressed " << i << ": " << log_bitcompressedsize[i] << "\t("
             << 100.f - (log_bitcompressedsize[i] / (float)totalsize) * 100.f << "%)" << std::endl;


        flog << "total blocks: " << log_totalblocks << std::endl;

        for (int i=16; i >= 0; i--) {
            if (log_actualbitblocks[i] > 0)
                flog << "stored blocks with " << i << " bits: "
                     << log_actualbitblocks[i] << " ("
                     << 100.f * (log_actualbitblocks[i] / (float)log_totalblocks) << "%)" << std::endl;
        }

        for (int i=16; i >= 0; i--) {
            flog << "blocks " << i << " bit: " << log_bitblocks[i];
            if (log_bitblocks[i] > 0)
                flog << "\t(" << (log_bitblocks[i] / (float)log_totalblocks) * 100.f << "%)";
            flog << std::endl;
        }

        flog << "non-zero uniform blocks: " << log_nonZeroUniformBlocks << std::endl;
        flog << "incompressible blocks: " << log_incompressible << std::endl;

    }
    catch (tgt::FileException) {
        LERROR("loading failed");
        return false;
    }
    catch (std::bad_alloc) {
        LERROR("out of memory");
        return false;
    } */

    return true;
}

} // namespace
