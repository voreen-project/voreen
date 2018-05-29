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

#ifndef VRN_BINARYMEDIAN_H
#define VRN_BINARYMEDIAN_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/temppathproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/ports/volumeport.h"

#include "../volumefiltering/volumefilter.h"
#include "../volumefiltering/slicereader.h"

namespace voreen {

struct BinaryMedianInput {
    std::unique_ptr<SliceReader> sliceReader;
    std::unique_ptr<HDF5FileVolume> outputVolume;

     BinaryMedianInput(std::unique_ptr<SliceReader>&& pSliceReader, std::unique_ptr<HDF5FileVolume>&& pOutputVolume)
         : sliceReader(std::move(pSliceReader))
         , outputVolume(std::move(pOutputVolume))
     {
     }

     BinaryMedianInput(const BinaryMedianInput&) = delete;
     BinaryMedianInput(BinaryMedianInput&& old)
         : sliceReader(old.sliceReader.release())
         , outputVolume(old.outputVolume.release())
     {
     }
};
struct BinaryMedianOutput {
    std::string outputVolumeFilePath;
};

class BinaryMedian : public AsyncComputeProcessor<BinaryMedianInput, BinaryMedianOutput> {
public:
    BinaryMedian();
    virtual ~BinaryMedian();

    virtual std::string getClassName() const         { return "BinaryMedian";      }
    virtual std::string getCategory() const       { return "Volume Processing"; }
    virtual bool isEndProcessor() const       { return true; }
    virtual VoreenSerializableObject* create() const;

    virtual void setDescriptions() {
        setDescription("Processor that performs a binary median operation on a binary input volume");
    }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;   }

    virtual void initialize();
    virtual bool isReady() const;

    virtual BinaryMedianInput prepareComputeInput();
    virtual BinaryMedianOutput compute(BinaryMedianInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(BinaryMedianOutput output);

protected:
    virtual void adjustPropertiesToInput();
    void updateObjectVoxelThreshold();
    void updateKernelSizeDisplay();

private:

    // Ports
    VolumePort inport_;
    VolumePort outport_;

    // General properties
    TempPathProperty outputVolumeFilePath_;
    IntProperty outputVolumeDeflateLevel_;

    FloatProperty  binarizationThreshold_;
    OptionProperty<SamplingStrategyType> samplingStrategyType_;
    IntProperty outsideVolumeValue_;

    BoolProperty enabled_;
    BoolProperty useUniformExtent_;
    IntProperty extentX_;
    IntProperty extentY_;
    IntProperty extentZ_;
    StringProperty kernelSizeDisplay_;

    BoolProperty forceMedian_;
    IntProperty objectVoxelThreshold_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_BINARYMEDIAN_H
