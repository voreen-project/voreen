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

#include "volumeinformation.h"

#include "voreen/core/datastructures/volume/histogram.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatornumsignificant.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

namespace voreen {

VolumeInformation::VolumeInformation()
    : VolumeProcessor()
    , volume_(Port::INPORT, "volume.inport", "Volume Input", false, Processor::INVALID_RESULT)
    , computeButton_("compute", "Compute information")
    , computeContinuously_("computeContinuously", "Compute Continuously")
    , numVoxels_("numVoxels", "Num voxels", 0, 0, 1<<30)
    , numSignificant_("numSignificant", "Num significant voxels", 0, 0, 1<<30)
    , minValue_("min", "Min intensity", 0.f, -1e16, 1e16)
    , maxValue_("max", "Max intensity", 0.f, -1e16, 1e16)
    , meanValue_("mean", "Mean intensity", 0.f, -1e16, 1e16)
    , standardDeviation_("dev", "Standard deviation", 0.f, 0.f, 1e16)
    , entropy_("entropy", "Entropy", 0.f, 0.f, 1e16)
{
    addPort(volume_);

    numVoxels_.setReadOnlyFlag(true);
    numSignificant_.setReadOnlyFlag(true);
    minValue_.setReadOnlyFlag(true);
    maxValue_.setReadOnlyFlag(true);
    entropy_.setReadOnlyFlag(true);
    meanValue_.setReadOnlyFlag(true);
    standardDeviation_.setReadOnlyFlag(true);
    computeButton_.onChange(MemberFunctionCallback<VolumeInformation>(this, &VolumeInformation::computeInformation));

    addProperty(computeButton_);
    addProperty(computeContinuously_);
    addProperty(numVoxels_);
    addProperty(numSignificant_);
    addProperty(minValue_);
    addProperty(maxValue_);
    addProperty(meanValue_);
    addProperty(standardDeviation_);
    addProperty(entropy_);
}

Processor* VolumeInformation::create() const {
    return new VolumeInformation;
}

void VolumeInformation::process() {
    if (volume_.hasChanged() && computeContinuously_.get())
        computeInformation();
}

void VolumeInformation::computeInformation() {

    if (!volume_.hasData()) {
        numVoxels_.set(0);
        numSignificant_.set(0);
        minValue_.set(0.f);
        maxValue_.set(0.f);
        meanValue_.set(0.f);
        standardDeviation_.set(0.f);
        entropy_.set(0.f);
        return;
    }

    const VolumeRAM* volume = volume_.getData()->getRepresentation<VolumeRAM>();
    tgtAssert(volume, "No input volume");
    size_t numVoxels = volume->getNumVoxels();
    tgt::vec2 intensityRange = volume->elementRange();

    numVoxels_.setMaxValue(static_cast<int>(numVoxels));
    numVoxels_.set(static_cast<int>(numVoxels));

    numSignificant_.setMaxValue(static_cast<int>(numVoxels));
    numSignificant_.set(static_cast<int>(VolumeOperatorNumSignificant::APPLY_OP(volume_.getData())));

    // compute the entropy of the volume
    // E(H) = -\sum_{i=0}^{K-1}{H_p(i) \log_2(H_p(i)}
    // with H_p(i): normalized histrogram with intensity value i
    int histMax = tgt::iround(intensityRange.y);
    // handle float data as if it was 16 bit to prevent overflow
    histMax = tgt::clamp(histMax, 0, 1<<16);//TODO

    Histogram1D histogram = createHistogram1DFromVolume(volume_.getData(), histMax);
    double entropy = 0.0;
    for (size_t i = 0; i < static_cast<size_t>(histMax); ++i) {
        uint64_t value = histogram.getBucket(static_cast<int>(i));

        if (value == 0)
            continue;

        double hi = static_cast<double>(value) / static_cast<double>(numVoxels);
        double loghi = log(hi);

        entropy += hi * loghi;
    }
    entropy *= -1;
    entropy_.setMaxValue(intensityRange.y);
    entropy_.set(static_cast<float>(entropy));

    // min value
    minValue_.setMinValue(intensityRange.x);
    minValue_.setMaxValue(intensityRange.y);
    //minValue_.set(static_cast<float>(histogram.getSignificantRange().x));
    float min = volume_.getData()->getDerivedData<VolumeMinMax>()->getMinNormalized();
    minValue_.set(min);

    //max value
    maxValue_.setMinValue(intensityRange.x);
    maxValue_.setMaxValue(intensityRange.y);
    //maxValue_.set(static_cast<float>(histogram.getSignificantRange().y));
    float max = volume_.getData()->getDerivedData<VolumeMinMax>()->getMaxNormalized();
    maxValue_.set(max);

    // mean value
    double mean = 0.0;
    for (size_t i = 0; i < numVoxels; ++i) {
        float value = volume->getVoxelNormalized(i);
        mean += value;
    }
    mean /= numVoxels;
    meanValue_.setMinValue(intensityRange.x);
    meanValue_.setMaxValue(intensityRange.y);
    meanValue_.set(static_cast<float>(mean));

    // variance
    double variance = 0.0;
    for (size_t i = 0; i < numVoxels; ++i) {
        float value = volume->getVoxelNormalized(i);
        variance += pow(static_cast<double>(value) - mean, 2.0);
    }
    variance /= numVoxels;
    standardDeviation_.setMaxValue(intensityRange.y);
    standardDeviation_.set(static_cast<float>(sqrt(variance)));
}

} // namespace
