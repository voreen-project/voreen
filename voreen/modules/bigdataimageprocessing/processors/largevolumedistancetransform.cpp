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

#include "largevolumedistancetransform.h"
#include "voreen/core/utils/stringutils.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"

#include "../algorithm/distancetransform.h"

#include "tgt/vector.h"

namespace voreen {

const std::string LargeVolumeDistanceTransform::loggerCat_("voreen.bigdataimageprocessing.LargeVolumeDistanceTransform");

LargeVolumeDistanceTransform::LargeVolumeDistanceTransform()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output",false)
    , binarizationThreshold_("binarizationThreshold", "Threshold", 0.5f, 0.0f, std::numeric_limits<float>::max(), Processor::INVALID_RESULT, FloatProperty::STATIC, Property::LOD_ADVANCED)
    , outputVolumeFilePath_("outputVolumeFilePath", "Output volume file path", "Output volume file path", "", LZ4SliceVolumeBase::FILE_EXTENSION)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(binarizationThreshold_);
        binarizationThreshold_.setTracking(false);
        binarizationThreshold_.setGroupID("binarization");

    addProperty(outputVolumeFilePath_);
}

LargeVolumeDistanceTransform::~LargeVolumeDistanceTransform() {}

Processor* LargeVolumeDistanceTransform::create() const {
    return new LargeVolumeDistanceTransform();
}

void LargeVolumeDistanceTransform::adjustPropertiesToInput() {
    const VolumeBase* input = inport_.getData();
    if(!input) {
        return;
    }
    if(!input->hasDerivedData<VolumeMinMax>()) {
        LINFO("Calculating VolumeMinMax. This may take a while...");
    }

    const VolumeMinMax* mm = input->getDerivedData<VolumeMinMax>();
    float min = mm->getMin();
    float max = mm->getMax();
    if(min != max) {
        binarizationThreshold_.setMinValue(min);
        binarizationThreshold_.setMaxValue(max);
        binarizationThreshold_.adaptDecimalsToRange(4);
    }
}

LargeVolumeDistanceTransform::ComputeInput LargeVolumeDistanceTransform::prepareComputeInput() {
    if(outputVolumeFilePath_.get().empty()) {
        throw InvalidInputException("No volume path specified", InvalidInputException::S_ERROR);
    }

    if(!inport_.hasData()) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    return LargeVolumeDistanceTransform::ComputeInput {
        inport_.getThreadSafeData(),
        outputVolumeFilePath_.get(),
        binarizationThreshold_.get(),
    };
}


LargeVolumeDistanceTransform::ComputeOutput LargeVolumeDistanceTransform::compute(LargeVolumeDistanceTransform::ComputeInput input, ProgressReporter& progressReporter) const {
    auto res = compute_distance_transform(*input.inputVolume_, input.binarizationThreshold_, input.outputPath_, progressReporter);

    return {
        std::move(res).toVolume(),
    };
}
void LargeVolumeDistanceTransform::processComputeOutput(LargeVolumeDistanceTransform::ComputeOutput output) {
    outport_.setData(output.outputVolume_.release());
}

}   // namespace
