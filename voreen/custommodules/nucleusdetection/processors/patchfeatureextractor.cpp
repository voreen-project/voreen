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

#include "patchfeatureextractor.h"
//#include "voreen/core/io/progressbar.h"

#include "voreen/core/datastructures/volume/volumebase.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/histogram.h"

#include "tgt/stopwatch.h"

#include "../util/featureextractors.h"
#include "../util/cosinetransform.h"

namespace voreen {

const std::string PatchFeatureExtractor::loggerCat_("voreen.nucleusdetection.PatchFeatureExtractor");
const std::vector<int> PatchFeatureExtractor::scales_{static_cast<int>(SCALE_ORIGINAL), static_cast<int>(SCALE_SECOND), static_cast<int>(SCALE_THIRD), static_cast<int>(SCALE_FOURTH)};

PatchFeatureExtractor::PatchFeatureExtractor()
    : VolumeProcessor()
    //, dataInport_(Port::INPORT, "volumehandle.input", "Volume Input", false)
    , filterInport1_(Port::INPORT, "filterlist.input1", "Filter Kernel List Input (Original resolution)", false)
    , filterInport2_(Port::INPORT, "filterlist.input2", "Filter Kernel List Input (Downsampling factor 2)", false)
    , filterInport3_(Port::INPORT, "filterlist.input3", "Filter Kernel List Input (Downsampling factor 3)", false)
    , filterInport4_(Port::INPORT, "filterlist.input4", "Filter Kernel List Input (Downsampling factor 4)", false)
    // patch settings  
    , patchSize_("patchSize", "3D Patch Size")
    , normalizePatches_("normalize", "Normalize patch intensities locally")
    // feature selection
    , rawPatch2D_("rawPatch2D", "Raw 2D Patches")
    , rawPatch3D_("rawPatch3D", "Raw 3D Patches")
    , filteredPatch3D_("filteredPatch3D", "Filtered 3D Patches")
    , cosineTransform1D_("dct1D", "1D DCT")
    , cosineTransform2D_("dct2D", "2D DCT")
    , intensityMean_("intensityMean", "Mean of Intensities\n(computed before local normalization)")
    , intensityStandardDeviation_("intensityStdDev", "Standard Deviation of Intensities\n(computed before local normalization)")
    // additional progress information
    , progressProperty_("progressProperty", "Progress")
{
    //addPort(dataInport_);
    
    // inports for filter kernels
    addPort(filterInport1_);
    addPort(filterInport2_);
    addPort(filterInport3_);
    addPort(filterInport4_);

    // patch settings (size and normalization of the scales)
    patchSize_.addOption("3",  "3x3x3",    3);
    patchSize_.addOption("5",  "5x5x5",    5);
    patchSize_.addOption("7",  "7x7x7",    7);
    patchSize_.addOption("9",  "9x9x9",    9);
    patchSize_.addOption("11", "11x11x11", 11);
    patchSize_.addOption("13", "13x13x13", 13);
    patchSize_.addOption("15", "15x15x15", 15);
    patchSize_.selectByKey("7");
    addProperty(patchSize_);

    fillSamplingScaleOptions(normalizePatches_);
    normalizePatches_.selectByKey("all");
    addProperty(normalizePatches_);

    patchSize_.setGroupID("patch");
    normalizePatches_.setGroupID("patch");
    setPropertyGroupGuiName("patch", "Patch Settings");

    // feature selection
    fillSamplingScaleOptions(rawPatch2D_);
    addProperty(rawPatch2D_);
    rawPatch2D_.setGroupID("features");

    fillSamplingScaleOptions(rawPatch3D_);
    addProperty(rawPatch3D_);
    rawPatch3D_.setGroupID("features");

    fillSamplingScaleOptions(filteredPatch3D_);
    addProperty(filteredPatch3D_);
    filteredPatch3D_.setGroupID("features");

    fillSamplingScaleOptions(cosineTransform1D_);
    addProperty(cosineTransform1D_);
    cosineTransform1D_.setGroupID("features");

    fillSamplingScaleOptions(cosineTransform2D_);
    addProperty(cosineTransform2D_);
    cosineTransform2D_.setGroupID("features");

    fillSamplingScaleOptions(intensityMean_);
    addProperty(intensityMean_);
    intensityMean_.setGroupID("features");

    fillSamplingScaleOptions(intensityStandardDeviation_);
    addProperty(intensityStandardDeviation_);
    intensityStandardDeviation_.setGroupID("features");

    setPropertyGroupGuiName("features", "Feature Selection");

    // progress information
    addProperty(progressProperty_);
    addProgressBar(&progressProperty_);
}

bool PatchFeatureExtractor::isReady() const {
    if (!isInitialized())
        return false;

    return true;
}

void PatchFeatureExtractor::fillSamplingScaleOptions(OptionProperty<int>& prop) {
    // only one scale
    prop.addOption("None",  "None",                         static_cast<int>(SCALE_NONE));
    prop.addOption("1",     "1 (Original Resolution)",      static_cast<int>(SCALE_ORIGINAL));
    prop.addOption("2",     "2 (Downsampling Factor 2)",    static_cast<int>(SCALE_SECOND));
    prop.addOption("3",     "3 (Downsampling Factor 3)",    static_cast<int>(SCALE_THIRD));
    prop.addOption("4",     "4 (Downsampling Factor 4)",    static_cast<int>(SCALE_FOURTH));
    // two scales
    prop.addOption("1&2",   "1 & 2",                        static_cast<int>(SCALE_ORIGINAL) | static_cast<int>(SCALE_SECOND));
    prop.addOption("1&3",   "1 & 3",                        static_cast<int>(SCALE_ORIGINAL) | static_cast<int>(SCALE_THIRD));
    prop.addOption("1&4",   "1 & 4",                        static_cast<int>(SCALE_ORIGINAL) | static_cast<int>(SCALE_FOURTH));
    prop.addOption("2&3",   "2 & 3",                        static_cast<int>(SCALE_SECOND)  | static_cast<int>(SCALE_THIRD));
    prop.addOption("2&4",   "2 & 4",                        static_cast<int>(SCALE_SECOND)  | static_cast<int>(SCALE_FOURTH));
    prop.addOption("3&4",   "3 & 4",                        static_cast<int>(SCALE_THIRD)   | static_cast<int>(SCALE_FOURTH));
    // three scales
    prop.addOption("1&2&3",   "1 & 2 & 3",                  static_cast<int>(SCALE_ORIGINAL) | static_cast<int>(SCALE_SECOND) | static_cast<int>(SCALE_THIRD));
    prop.addOption("1&2&4",   "1 & 2 & 4",                  static_cast<int>(SCALE_ORIGINAL) | static_cast<int>(SCALE_SECOND) | static_cast<int>(SCALE_FOURTH));
    prop.addOption("1&3&4",   "1 & 3 & 4",                  static_cast<int>(SCALE_ORIGINAL) | static_cast<int>(SCALE_THIRD) | static_cast<int>(SCALE_FOURTH));
    prop.addOption("2&3&4",   "2 & 3 & 4",                  static_cast<int>(SCALE_SECOND)  | static_cast<int>(SCALE_THIRD) | static_cast<int>(SCALE_FOURTH));
    // four scales
    prop.addOption("all",     "All (1 & 2 & 3 & 4)",        static_cast<int>(SCALE_ORIGINAL) | static_cast<int>(SCALE_SECOND) | static_cast<int>(SCALE_THIRD) | static_cast<int>(SCALE_FOURTH));

    prop.selectByKey("None");
}

size_t PatchFeatureExtractor::computeNumFeatures(std::vector<std::pair<size_t,float*> > filterKernels) const {
    // compute sizes of 1d, 2d, and 3d patches
    size_t patchSize = static_cast<size_t>(patchSize_.getValue());
    size_t patchSize2d = patchSize * patchSize;
    size_t patchSize3d = patchSize2d * patchSize; 

    size_t features = 0;
    int currentFeature = 0;

    // raw patches (2D)
    currentFeature = rawPatch2D_.getValue();
    for (auto scale : scales_) {
        if (currentFeature & scale)  // 3 x 2D patch size
            features += 3 * patchSize2d;
    }

    // raw patches (3D)
    currentFeature = rawPatch3D_.getValue();
    for (auto scale : scales_) {
        if (currentFeature & scale)  // 3D patch size
            features += patchSize3d;
    }

    // filtered 3D patches
    for (auto filter : filterKernels) {
        features += filter.first;
    }

    // DCT (1D)
    currentFeature = cosineTransform1D_.getValue();
    for (auto scale : scales_) {
        if (currentFeature & scale)  // 3 x 1D patch size
            features += 3 * patchSize;
    }

    // DCT (2D)
    currentFeature = cosineTransform2D_.getValue();
    for (auto scale : scales_) {
        if (currentFeature & scale)  // 3 x 2D patch size
            features += 3 * patchSize2d;
    }

    // mean and std dev
    currentFeature = intensityMean_.getValue();
    for (auto scale : scales_) {
        if (currentFeature & scale)
            features++;
    }
    currentFeature = intensityStandardDeviation_.getValue();
    for (auto scale : scales_) {
        if (currentFeature & scale)
            features++;
    }

    return features;
}

tgt::bvec4 PatchFeatureExtractor::determineNecessaryScales() const {
    tgt::bvec4 result = tgt::bvec4(false);
    int scales = rawPatch2D_.getValue() | rawPatch3D_.getValue() | filteredPatch3D_.getValue() | cosineTransform1D_.getValue() | cosineTransform2D_.getValue() | intensityMean_.getValue() | intensityStandardDeviation_.getValue(); 
    for (size_t i = 0; i < 4; ++i) {
        if (scales & scales_.at(i))
            result[i] = true;
    }
    return result;
}

size_t PatchFeatureExtractor::computeNecessaryBorder(tgt::bvec4 necessaryScales) const {
    size_t radius = static_cast<size_t>(patchSize_.getValue()) / 2;
    if (necessaryScales[3])
        radius *= 4;
    else if (necessaryScales[2])
        radius *= 3;
    else if (necessaryScales[1])
        radius *= 2;

    return radius;
}

std::vector<float*> PatchFeatureExtractor::extractScalePatches(const VolumeRAM* vol, tgt::svec3 voxelPos, tgt::vec2 scalingMinMax, tgt::bvec4 necessaryScales) const {
    std::vector<float*> scales = {nullptr, nullptr, nullptr, nullptr};
    size_t radius = static_cast<size_t>(patchSize_.getValue()) / 2;

    for (size_t i = 0; i < 4; ++i) { 
        if (!necessaryScales[i])
            continue;
        else {
            scales.at(i) = new float[patchSize_.getValue() * patchSize_.getValue() * patchSize_.getValue()];
            float* patch = scales.at(i);
            // extract patch on the current scale
            size_t stride = i + 1;
            size_t patchIndex = 0;
            for (size_t cz = voxelPos.z - (radius * stride); cz <= voxelPos.z + (radius * stride); cz = cz + stride) {
                for (size_t cy = voxelPos.y - (radius * stride); cy <= voxelPos.y + (radius * stride); cy = cy + stride) {
                    for (size_t cx = voxelPos.x - (radius * stride); cx <= voxelPos.x + (radius * stride); cx = cx + stride) {
                        // get the raw value and normalize it with the min and max (99%) intensity
                        patch[patchIndex++] = scaleIntensity(vol->getVoxelNormalized(cx, cy, cz), scalingMinMax);
                    }
                }
            }
        }
    }

    return scales;
}

void PatchFeatureExtractor::freePatchMemory(std::vector<float*> scales) const {
    for (auto& p : scales) {
        delete[] p;
        p = nullptr;
    }
} 

void PatchFeatureExtractor::computeFeaturesForVoxel(Sample& sample, const VolumeRAM* volume, tgt::svec3 voxelPos, std::vector<std::pair<size_t, float*> > filterKernels, float* dctCoefficients, tgt::bvec4 necessaryScales, tgt::vec2 minMaxIntensity)  const {
    // get necessary patches
    std::vector<float*> patches = extractScalePatches(volume, voxelPos, minMaxIntensity, necessaryScales);

    size_t patchSize = patchSize_.getValue();
    size_t patchSize2D = patchSize * patchSize;
    size_t patchSize3D = patchSize2D * patchSize;

    size_t currentFeatureOffset = 0;
    for (size_t i = 0; i < 4; ++i) {
        // scale not necessary -> no features
        if (!necessaryScales[i])
            continue;

        float* patch = patches.at(i);

        // compute mean and sigma if necessary
        float mean = 0.f; float sigma = 1.f;
        if (scales_.at(i) & (intensityMean_.getValue() | intensityStandardDeviation_.getValue() | normalizePatches_.getValue())) {
            float sum = 0.f;
            // get mean and standard deviation to normalize the patch
            for (size_t index = 0; index < patchSize3D; ++index) {
                sum += patch[index];
            }
            mean = sum / static_cast<float>(patchSize3D);

            sigma = 0.f;
            for (size_t index = 0; index < patchSize3D; ++index) {
                sigma += std::pow(patch[index] - mean, 2.f);
            }
            sigma *= 1.f / static_cast<float>(patchSize3D - 1);
            sigma = std::sqrt(sigma);

            // if sigma is near 0: set to 1
            if (sigma < 1e-6f)
                sigma = 1.f;
        }

        if (scales_.at(i) & normalizePatches_.getValue()) {
            // use mean and standard deviation to normalize the patch
            for (size_t index = 0; index < patchSize3D; ++index) {
                float value = (patch[index] - mean) / sigma;
                patch[index] = value;
            }
        }

        // now extract all features and write them into the sample
        size_t offset2_5D = 0;

        // 2.5D patches
        if (scales_.at(i) & rawPatch2D_.getValue()) {
            offset2_5D = currentFeatureOffset;
            currentFeatureOffset += sample2_5raw(patchSize, patch, &sample.features_[currentFeatureOffset]); 
        }

        // 3D patches
        if (scales_.at(i) & rawPatch3D_.getValue()) {
            currentFeatureOffset += sample3raw(patchSize, patch, &sample.features_[currentFeatureOffset]);
        }

        if (scales_.at(i) & filteredPatch3D_.getValue()) {
            currentFeatureOffset += sample3filtered(patchSize, patch, filterKernels.at(i).first, filterKernels.at(i).second, //
                                                    &sample.features_[currentFeatureOffset]);
        }
        
        if (scales_.at(i) & cosineTransform1D_.getValue()) {
            currentFeatureOffset += sample1dct(patchSize, patch, dctCoefficients,//
                                               &sample.features_[currentFeatureOffset]);
        }

        if (scales_.at(i) & cosineTransform2D_.getValue()) {
            float* patches2_5D = 0;
            if (scales_.at(i) & rawPatch2D_.getValue())
                patches2_5D = &sample.features_[offset2_5D];
            else {
                patches2_5D = new float[3 * patchSize2D];
                sample2_5raw(patchSize, patch, patches2_5D);
            }

            currentFeatureOffset += sample2_5dct(patchSize, patches2_5D, dctCoefficients,//
                                                 &sample.features_[currentFeatureOffset]);

            if (!(scales_.at(i) & rawPatch2D_.getValue()))
                delete[] patches2_5D;
        }

        // append mean and standard deviation
        if (scales_.at(i) & intensityMean_.getValue()) 
            sample.features_[currentFeatureOffset++] = mean;

        if (scales_.at(i) & intensityStandardDeviation_.getValue()) 
            sample.features_[currentFeatureOffset++] = sigma;
    }

    freePatchMemory(patches);
}

std::vector<std::pair<size_t,float*> > PatchFeatureExtractor::getFilterKernels(size_t patchSize) {
    std::vector<std::pair<size_t, float*> > kernels; // = {0, 0, 0, 0};

    if (!checkForRequiredFilterKernels()) {
        LERROR("Cannot get required filter kernels!");
        return kernels;
    }

    LINFO("Retrieving filter kernels...");
    tgt::Stopwatch timer;
    timer.start();
    
    // kernels are all provided
    if (filteredPatch3D_.getValue() & static_cast<int>(SCALE_ORIGINAL)) {
        size_t numFilters = filterInport1_.getData()->size();
        float* filters = new float[patchSize * patchSize * patchSize * numFilters];

        // get all of the filters
        for (size_t i = 0; i < numFilters; ++i) {
            const VolumeBase* v = filterInport1_.getData()->at(i);
            if (v->getDimensions() != tgt::svec3(patchSize,patchSize,patchSize) || v->getFormat() != "float") {
                LERROR(" float filter kernels of size " << patchSize << "^3 required ");
                timer.stop();
                delete[] filters;
                return kernels;
            }
            const VolumeRAM* ramFilter = v->getRepresentation<VolumeRAM>();
            if (!ramFilter) {
                LERROR("Could not get RAM representation of filter kernel " << i);
                timer.stop();
                delete[] filters;
                return kernels;
            }
            std::memcpy((void*) &filters[patchSize * patchSize * patchSize * i], ramFilter->getData(), patchSize * patchSize * patchSize * sizeof(float));
        }

        kernels.push_back(std::make_pair(numFilters, filters));
    }
    else
        kernels.push_back(std::make_pair(0,nullptr));

    if (filteredPatch3D_.getValue() & static_cast<int>(SCALE_SECOND)) {
        size_t numFilters = filterInport2_.getData()->size();
        float* filters = new float[patchSize * patchSize * patchSize * numFilters];

        // get all of the filters
        for (size_t i = 0; i < numFilters; ++i) {
            const VolumeBase* v = filterInport2_.getData()->at(i);
            if (v->getDimensions() != tgt::svec3(patchSize,patchSize,patchSize) || v->getFormat() != "float") {
                LERROR(" float filter kernels of size " << patchSize << "^3 required ");
                timer.stop();
                delete[] filters;
                return kernels;
            }
            const VolumeRAM* ramFilter = v->getRepresentation<VolumeRAM>();
            if (!ramFilter) {
                LERROR("Could not get RAM representation of filter kernel " << i);
                timer.stop();
                delete[] filters;
                return kernels;
            }
            std::memcpy((void*) &filters[patchSize * patchSize * patchSize * i], ramFilter->getData(), patchSize * patchSize * patchSize * sizeof(float));
        }

        kernels.push_back(std::make_pair(numFilters, filters));
    }
    else
        kernels.push_back(std::make_pair(0, nullptr));

    if (filteredPatch3D_.getValue() & static_cast<int>(SCALE_THIRD)) {
        size_t numFilters = filterInport3_.getData()->size();
        float* filters = new float[patchSize * patchSize * patchSize * numFilters];

        // get all of the filters
        for (size_t i = 0; i < numFilters; ++i) {
            const VolumeBase* v = filterInport3_.getData()->at(i);
            if (v->getDimensions() != tgt::svec3(patchSize,patchSize,patchSize) || v->getFormat() != "float") {
                LERROR(" float filter kernels of size " << patchSize << "^3 required ");
                timer.stop();
                delete[] filters;
                return kernels;
            }
            const VolumeRAM* ramFilter = v->getRepresentation<VolumeRAM>();
            if (!ramFilter) {
                LERROR("Could not get RAM representation of filter kernel " << i);
                timer.stop();
                delete[] filters;
                return kernels;
            }
            std::memcpy((void*) &filters[patchSize * patchSize * patchSize * i], ramFilter->getData(), patchSize * patchSize * patchSize * sizeof(float));
        }

        kernels.push_back(std::make_pair(numFilters, filters));
    }
    else
        kernels.push_back(std::make_pair(0,nullptr));

    if (filteredPatch3D_.getValue() & static_cast<int>(SCALE_FOURTH)) {
        size_t numFilters = filterInport4_.getData()->size();
        float* filters = new float[patchSize * patchSize * patchSize * numFilters];

        // get all of the filters
        for (size_t i = 0; i < numFilters; ++i) {
            const VolumeBase* v = filterInport4_.getData()->at(i);
            if (v->getDimensions() != tgt::svec3(patchSize,patchSize,patchSize) || v->getFormat() != "float") {
                LERROR(" float filter kernels of size " << patchSize << "^3 required ");
                timer.stop();
                delete[] filters;
                return kernels;
            }
            const VolumeRAM* ramFilter = v->getRepresentation<VolumeRAM>();
            if (!ramFilter) {
                LERROR("Could not get RAM representation of filter kernel " << i);
                timer.stop();
                delete[] filters;
                return kernels;
            }
            std::memcpy((void*) &filters[patchSize * patchSize * patchSize * i], ramFilter->getData(), patchSize * patchSize * patchSize * sizeof(float));
        }

        kernels.push_back(std::make_pair(numFilters, filters));
    }
    else
        kernels.push_back(std::make_pair(0,nullptr));

    if (kernels.size() != 4) {
        LERROR("Could not load all filter kernels");
        timer.stop();
    }
    else {
        timer.stop();
        LINFO("Time for retrieving filter kernels: " << (float) timer.getRuntime() / 1000.f << " seconds");
    }
    
    return kernels;
}

void PatchFeatureExtractor::deleteFilterKernels(std::vector<std::pair<size_t,float*> > kernels) {
    for (auto& kernel : kernels) {
        delete[] kernel.second;
        kernel.first = 0;
        kernel.second = nullptr;
    }
} 

bool PatchFeatureExtractor::checkForRequiredFilterKernels() const {
    if (filteredPatch3D_.getValue() & static_cast<int>(SCALE_ORIGINAL)) {
        if (!filterInport1_.isReady() || !filterInport1_.getData() || filterInport1_.getData()->empty()) {
            LERROR("Found no filter kernels for original resolution - cannot determine number of features");
            return false;
        }
    }

    if (filteredPatch3D_.getValue() & static_cast<int>(SCALE_SECOND)) {
        if (!filterInport2_.isReady() || !filterInport2_.getData() || filterInport2_.getData()->empty()) {
            LERROR("Found no filter kernels for downsampling 2 - cannot determine number of features");
            return false;
        }
    }

    if (filteredPatch3D_.getValue() & static_cast<int>(SCALE_THIRD)) { 
        if (!filterInport3_.isReady() || !filterInport3_.getData() || filterInport3_.getData()->empty()) {
            LERROR("Found no filter kernels for downsampling 3 - cannot determine number of features");
            return false;
        }
    }

    if (filteredPatch3D_.getValue() & static_cast<int>(SCALE_FOURTH)) { 
        if (!filterInport4_.isReady() || !filterInport4_.getData() || filterInport4_.getData()->empty()) {
            LERROR("Found no filter kernels for downsampling 4 - cannot determine number of features");
            return false;
        }
    }

    return true;
}

tgt::vec2 PatchFeatureExtractor::getScalingMinMaxIntensity(const VolumeBase* volume) const {
    // TODO: channel selection would be nice...
    float minIntensity = volume->getDerivedData<VolumeMinMax>()->getMinNormalized();
    float maxTmp = volume->getDerivedData<VolumeMinMax>()->getMaxNormalized();
    // for the maximum, we only want to use 99% of the data to ignore outliers
    VolumeHistogramIntensity* histogram = volume->getDerivedData<VolumeHistogramIntensity>();
    float sumBuckets = 0;
    for (size_t i = 0; i < histogram->getBucketCount(); ++i) {
        sumBuckets += histogram->getNormalized(static_cast<int>(i));
    }
    int i = -1;
    float currentSum = 0.f;
    while (((currentSum / sumBuckets) < 0.99f) && i < (int) histogram->getBucketCount()) {
        i += 1;
        currentSum += histogram->getNormalized(i);
    }
    float maxIntensity = (static_cast<float>(i) / static_cast<float>(histogram->getBucketCount() - 1)) * (maxTmp - minIntensity) + minIntensity;
    tgtAssert(maxIntensity - minIntensity > 0, "max = min intensity");

    return tgt::vec2(minIntensity, maxIntensity);
}

float PatchFeatureExtractor::scaleIntensity(float value, tgt::vec2 minMaxValue) const {
    // scale, do not clamp to 1
    return (value - minMaxValue.x) / (minMaxValue.y - minMaxValue.x);
}

void PatchFeatureExtractor::deactivateProperties() {
    patchSize_.setReadOnlyFlag(true);
    normalizePatches_.setReadOnlyFlag(true);
    rawPatch2D_.setReadOnlyFlag(true);
    rawPatch3D_.setReadOnlyFlag(true);
    filteredPatch3D_.setReadOnlyFlag(true);
    cosineTransform1D_.setReadOnlyFlag(true);
    cosineTransform2D_.setReadOnlyFlag(true);
    intensityMean_.setReadOnlyFlag(true);
    intensityStandardDeviation_.setReadOnlyFlag(true);
}

void PatchFeatureExtractor::reactivateProperties() {
    patchSize_.setReadOnlyFlag(false);
    normalizePatches_.setReadOnlyFlag(false);
    rawPatch2D_.setReadOnlyFlag(false);
    rawPatch3D_.setReadOnlyFlag(false);
    filteredPatch3D_.setReadOnlyFlag(false);
    cosineTransform1D_.setReadOnlyFlag(false);
    cosineTransform2D_.setReadOnlyFlag(false);
    intensityMean_.setReadOnlyFlag(false);
    intensityStandardDeviation_.setReadOnlyFlag(false);
}

} // namespace
