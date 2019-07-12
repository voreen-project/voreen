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

#ifndef VRN_LARGEVOLUMEFORMATCONVERSION_H
#define VRN_LARGEVOLUMEFORMATCONVERSION_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include <string>
#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/temppathproperty.h"
#include "../datastructures/lz4slicevolume.h"

namespace voreen {

struct LargeVolumeFormatConversionInput {
    std::string outputPath_;
    std::string outputBaseType_;
    const VolumeBase* inputVolume_;
};
struct LargeVolumeFormatConversionOutput {
    std::unique_ptr<Volume> outputVolume_;
};
class VRN_CORE_API LargeVolumeFormatConversion : public AsyncComputeProcessor<LargeVolumeFormatConversionInput, LargeVolumeFormatConversionOutput> {
public:
    typedef LargeVolumeFormatConversionOutput ComputeOutput;
    typedef LargeVolumeFormatConversionInput ComputeInput;

    LargeVolumeFormatConversion();
    virtual ~LargeVolumeFormatConversion();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "LargeVolumeFormatConversion"; }
    virtual std::string getCategory() const  { return "Volume Processing";      }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL;        }

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

    virtual bool usesExpensiveComputation() const { return true; }
    virtual void adjustPropertiesToInput();

protected:
    virtual void setDescriptions() {
        setDescription(
                "Converts a multi-channel input volume to an arbitrary multi-channel target format. "
                "This processor operates on on-disk volumes and is thus capable of processing volumes that exceed the main memory in size."
                );
    }
private:
    VolumePort inport_;
    VolumePort outport_;

    BoolProperty enableProcessing_;
    StringOptionProperty targetBaseType_;
    IntProperty numChannels_;
    TempPathProperty outputVolumeFilePath_;

    static const std::string loggerCat_; ///< category used in logging
};

}   //namespace

#endif // VRN_LARGEVOLUMEFORMATCONVERSION_H
