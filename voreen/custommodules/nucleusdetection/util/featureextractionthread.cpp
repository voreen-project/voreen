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

#include "featureextractionthread.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeram.h"

#include "featureextractors.h"
#include "cosinetransform.h"

namespace voreen {

FeatureExtractionThread::FeatureExtractionThread(const PatchCaffeClassifier* processor, size_t localNumVoxels, tgt::svec3 outputVolumeDim, size_t globalOffset, size_t localOffset, size_t border,
                            size_t numFeatures, float* currentBatch, const VolumeRAM* volume, std::vector<std::pair<size_t, float*> > filterKernels, float* dctCoefficients, tgt::bvec4 necessaryScales, tgt::vec2 minMaxIntensity)
    : BackgroundThread()
    , processor_(processor) 
    , localNumVoxels_(localNumVoxels)
    , outputVolumeDim_(outputVolumeDim)
    , globalOffset_(globalOffset)
    , localOffset_(localOffset)
    , border_(border)
    , numFeatures_(numFeatures)
    , currentBatch_(currentBatch)
    , volume_(volume)
    , filterKernels_(filterKernels)
    , dctCoefficients_(dctCoefficients)
    , necessaryScales_(necessaryScales)
    , minMaxIntensity_(minMaxIntensity)
    , localOutputData_(0)
{ }

void FeatureExtractionThread::threadMain() {
    // allocate local result memory
    localOutputData_ = new float[localNumVoxels_ * numFeatures_];

    // compute total number of voxels in the output volume
    size_t totalNumVoxels = tgt::hmul(outputVolumeDim_);

    PatchFeatureExtractor::Sample sample(numFeatures_);

    size_t actualNumVoxels = 0;
    for (size_t i = 0; i < localNumVoxels_; ++i) {
        // compute position of voxel in output volume
        tgt::svec3 voxel = tgt::svec3::zero;
        size_t totalOffset = globalOffset_ + localOffset_ + i;
        
        if (totalOffset >= totalNumVoxels)
            break;
        
        voxel.z = totalOffset / (outputVolumeDim_.x * outputVolumeDim_.y);   // slice has x * y voxels        
        size_t myXYOffset = totalOffset - voxel.z * (outputVolumeDim_.x * outputVolumeDim_.y);
        voxel.y = myXYOffset / outputVolumeDim_.x; 
        size_t myXOffset = myXYOffset - voxel.y * outputVolumeDim_.x;
        voxel.x = myXOffset;
            
        tgtAssert(!tgt::hor(tgt::greaterThanEqual(voxel, outputVolumeDim_)), "Sampling outside allowed range in input volume!");  
        
        // compute position in input volume
        tgt::svec3 currentPos = voxel + tgt::svec3(border_);  
    
        processor_->computeFeaturesForVoxel(sample, volume_, currentPos, filterKernels_, dctCoefficients_, necessaryScales_, minMaxIntensity_);

        // copy the feature vector into the local results
        std::memcpy(&localOutputData_[numFeatures_ * i], &sample.features_[0], numFeatures_ * sizeof(float));

        actualNumVoxels++;
    }

    // copy the local results into the current batch
    if (actualNumVoxels > 0)
        std::memcpy(&currentBatch_[localOffset_ * numFeatures_], &localOutputData_[0], actualNumVoxels * numFeatures_ * sizeof(float));
    
    delete[] localOutputData_; 
    localOutputData_ = 0;
}

}   // namespace
