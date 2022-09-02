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

#include "../algorithm/distancetransform.h"

#include "tgt/vector.h"

namespace voreen {

const std::string LargeVolumeDistanceTransform::loggerCat_("voreen.bigdataimageprocessing.LargeVolumeDistanceTransform");

LargeVolumeDistanceTransform::LargeVolumeDistanceTransform()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output",false)
    , outputVolumeFilePath_("outputVolumeFilePath", "Output volume file path", "Output volume file path", "", LZ4SliceVolumeBase::FILE_EXTENSION)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(outputVolumeFilePath_);
}

LargeVolumeDistanceTransform::~LargeVolumeDistanceTransform() {}

Processor* LargeVolumeDistanceTransform::create() const {
    return new LargeVolumeDistanceTransform();
}

void LargeVolumeDistanceTransform::adjustPropertiesToInput() {
    //
}

LargeVolumeDistanceTransform::ComputeInput LargeVolumeDistanceTransform::prepareComputeInput() {
    if(outputVolumeFilePath_.get().empty()) {
        throw InvalidInputException("No volume path specified", InvalidInputException::S_ERROR);
    }

    if(!inport_.hasData()) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    return LargeVolumeDistanceTransform::ComputeInput {
        outputVolumeFilePath_.get(),
        inport_.getThreadSafeData()
    };
}


LargeVolumeDistanceTransform::ComputeOutput LargeVolumeDistanceTransform::compute(LargeVolumeDistanceTransform::ComputeInput input, ProgressReporter& progressReporter) const {
    const float binarizationThreshold = 0.5f;
    auto res = compute_distance_transform(*input.inputVolume_, binarizationThreshold, input.outputPath_, progressReporter);

    return {
        std::move(res).toVolume(),
    };
}
void LargeVolumeDistanceTransform::processComputeOutput(LargeVolumeDistanceTransform::ComputeOutput output) {
    outport_.setData(output.outputVolume_.release());
}

}   // namespace
