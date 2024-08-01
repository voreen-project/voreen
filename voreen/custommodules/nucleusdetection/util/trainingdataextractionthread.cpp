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

#include "trainingdataextractionthread.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeram.h"

#include "featureextractors.h"
#include "cosinetransform.h"

namespace voreen {

TrainingDataExtractionThread::TrainingDataExtractionThread(const PatchTrainingDataExtractor* processor, tgt::svec3 inputVolumeDim, size_t firstPosition, size_t lastPosition, const std::vector<tgt::svec3>& positions, 
        size_t border, size_t numFeatures, const VolumeRAM* volume, const VolumeRAM* labelVolume, std::vector<std::pair<size_t, float*> > filterKernels, float* dctCoefficients, tgt::bvec4 necessaryScales, tgt::vec2 minMaxIntensity)
    : BackgroundThread()
    , processor_(processor) 
    , inputVolumeDim_(inputVolumeDim)
    , firstPosition_(firstPosition)
    , lastPosition_(lastPosition)
    , positions_(positions)
    , border_(border)
    , numFeatures_(numFeatures)
    , volume_(volume)
    , labelVolume_(labelVolume)
    , filterKernels_(filterKernels)
    , dctCoefficients_(dctCoefficients)
    , necessaryScales_(necessaryScales)
    , minMaxIntensity_(minMaxIntensity)
{ }

std::vector<PatchFeatureExtractor::Sample>& TrainingDataExtractionThread::getForegroundVoxels() {
    boost::mutex::scoped_lock lock(resultsMutex_);
    return foregroundVoxels_;
}

std::vector<PatchFeatureExtractor::Sample>& TrainingDataExtractionThread::getBackgroundVoxels() {
    boost::mutex::scoped_lock lock(resultsMutex_);
    return backgroundVoxels_;
}

void TrainingDataExtractionThread::threadMain() {

    // gets its own mutex to prevent other threads from readin unfinished results
    boost::mutex::scoped_lock lock(resultsMutex_);

    PatchFeatureExtractor::Sample sample(numFeatures_);

    // only process local positions
    for (size_t i = firstPosition_; i <= lastPosition_; ++i) {
        tgt::svec3 currentPos = positions_.at(i);

        processor_->computeFeaturesForVoxel(sample, volume_, currentPos, filterKernels_, dctCoefficients_, necessaryScales_, minMaxIntensity_);

        // get the label and append to one of the lists
        float currentLabel = labelVolume_->getVoxelNormalized(currentPos);
        if (currentLabel == 0.f)
            backgroundVoxels_.push_back(sample);
        else
            foregroundVoxels_.push_back(sample);
    }
}

}   // namespace
