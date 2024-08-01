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

#ifndef VRN_TRAININGDATAEXTRACTIONTHREAD
#define VRN_TRAININGDATAEXTRACTIONTHREAD

#include "voreen/core/voreencoreapi.h"

#include "voreen/core/utils/backgroundthread.h"

#include "../processors/patchtrainingdataextractor.h"

namespace voreen {

class VRN_CORE_API TrainingDataExtractionThread : public BackgroundThread { 

public:

    TrainingDataExtractionThread(const PatchTrainingDataExtractor* processor, tgt::svec3 inputVolumeDim, size_t firstPosition, size_t lastPosition, const std::vector<tgt::svec3>& positions, size_t border, size_t numFeatures, 
            const VolumeRAM* volume, const VolumeRAM* labelVolume, std::vector<std::pair<size_t, float*> > filterKernels, float* dctCoefficients, tgt::bvec4 necessaryScales, tgt::vec2 minMaxIntensity); 

    /// getter for foreground voxels, only call after the thread is finished!
    std::vector<PatchFeatureExtractor::Sample>& getForegroundVoxels();

    /// getter for background voxels, only call after the thread is finished!
    std::vector<PatchFeatureExtractor::Sample>& getBackgroundVoxels();

protected:

    // do not use the default constructor
    TrainingDataExtractionThread();

    virtual void threadMain();

    const PatchTrainingDataExtractor* processor_;   ///< classifier processor this thread is started from
    tgt::svec3 inputVolumeDim_;                     ///< dimensions of the input volume
    size_t firstPosition_;                             ///< first slice to be processed by this thread
    size_t lastPosition_;                              ///< last slice to be processed by this thread
    const std::vector<tgt::svec3>& positions_; 
    size_t border_;                                 ///< border of the volume (= necessary radius for feature extraction)

    size_t numFeatures_;                            ///< size of a single feature vector

    const VolumeRAM* volume_;                       ///< input volume
    const VolumeRAM* labelVolume_;                  ///< foreground / background labels
    std::vector<std::pair<size_t, float*> > filterKernels_;
    float* dctCoefficients_;
    tgt::bvec4 necessaryScales_;
    tgt::vec2 minMaxIntensity_;

    boost::mutex resultsMutex_;  ///< mutex to secure that results can only be read after the main function of the thread is finished
    std::vector<PatchFeatureExtractor::Sample> foregroundVoxels_;   // foreground results
    std::vector<PatchFeatureExtractor::Sample> backgroundVoxels_;   // background results
};

} // namespace

#endif 
