/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_VOLUMESTREAMPROCESSOR_H
#define VRN_VOLUMESTREAMPROCESSOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/filedialogproperty.h"

#include "../properties/volumestreamproperty.h"

#include <deque>

namespace voreen {

class Volume;

/**
 * Volume data set supplier in the network.
 *
 */
class VRN_CORE_API VolumeStreamProcessor : public Processor {
public:
    VolumeStreamProcessor();
    ~VolumeStreamProcessor();

    virtual std::string getCategory() const { return "Data Source"; }
    virtual std::string getClassName() const { return "VolumeStream"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }

    virtual Processor* create() const;

    virtual void process();

    virtual void initialize();
    virtual void deinitialize();

    void setStep(int step);
    int getStepCount() const;
    uint64_t getRawSize() const;

    struct StreamStep {
        StreamStep(int step = -1)
            : step_(step), loaded_(false), processed_(false), uploaded_(false), data_(0)
        {}

        int step_;
        bool loaded_;
        bool processed_;
        bool uploaded_;
        char* data_;
    };

    enum CompressionType { UNCOMPRESSED = 0, LZO };
    enum PredictionType { DIRECT = 0, DELTA };

#pragma pack(push)  // push current alignment to stack
#pragma pack(1)     // set alignment to 1 byte boundary

    struct BlockInfo {
        uint64_t offset_;
        uint32_t size_;
        uint16_t minValue_;
        uint16_t maxValue_;
        uint16_t referenceValue_;
        uint8_t prediction_;
        uint8_t bits_;
    };

#pragma pack(pop)   // restore original alignment from stack

    struct BlockData {
        BlockData()
            : info_(0), rawData_(0), data_(0), finalSize_(0), loaded_(false), ready_(false)
        {}

        BlockInfo* info_;
        char* rawData_;
        uint16_t* data_;
        int finalSize_;
        bool loaded_;
        bool ready_;
    };

    struct StepInfo {
        StepInfo()
            : ready_(false)
        {}
        std::vector<BlockInfo> blocks_;
        std::vector<BlockData> datas_;
        bool ready_;
    };

    std::deque<StreamStep>& getWorkList() { return workList_; }

    void readBlocks(int step);
    void uncompressBlocks(int step);
    void includeBlocks(int step, int isDelta);
    void showStep(int step);

    void flush();

    void loadStep();
    int getBlockSize() const;

protected:
    virtual void setDescriptions() {
        setDescription("");
    }
    void loadStream();

    void loadStepFromStream();
    void loadStepFromFile();

    Volume* volume_;

    FileDialogProperty filename_;
    IntProperty step_;
    VolumeStreamProperty stream_;

    VolumePort outport_;

    std::string streamFileName_;
    std::vector<std::string> files_;
    std::vector<std::string> deltas_;
    float spreadMin_, spreadMax_;
    int blockSize_;
    tgt::ivec3 blockCount_;
    int lastStep_;
    CompressionType compression_;
    PredictionType prediction_;

    std::vector<StepInfo> steps_;
    std::deque<StreamStep> workList_;

    std::ifstream* streamFile_;

    uint64_t rawSize_;

    static const std::string loggerCat_;
};

} // namespace

#endif
