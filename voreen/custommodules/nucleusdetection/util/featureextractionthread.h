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

#ifndef VRN_FEATUREEXTRACTIONTHREAD
#define VRN_FEATUREEXTRACTIONTHREAD

#include "voreen/core/voreencoreapi.h"

#include "voreen/core/utils/backgroundthread.h"

#include "../processors/patchcaffeclassifier.h"

namespace voreen {

class VRN_CORE_API FeatureExtractionThread : public BackgroundThread { 

public:

    FeatureExtractionThread(const PatchCaffeClassifier* processor, size_t localNumVoxels, tgt::svec3 outputVolumeDim, size_t globalOffset, size_t localOffset, size_t border,
                            size_t numFeatures, float* currentBatch, const VolumeRAM* volume, std::vector<std::pair<size_t, float*> > filterKernels, float* dctCoefficients, tgt::bvec4 necessaryScales, tgt::vec2 minMaxIntensity); 

protected:

    // do not use the default constructor
    FeatureExtractionThread();

    virtual void threadMain();

    // no interruption handling since the thread is not meant to be interrupted
    //virtual void handleInterruption();

    float* localOutputData_;

    const PatchCaffeClassifier* processor_;   ///< classifier processor this thread is started from
    size_t localNumVoxels_;             ///< number of voxels to process in this thread
    tgt::svec3 outputVolumeDim_;        ///< dimensions of the output volume
    size_t globalOffset_;               ///< current global offset depending on the batch number
    size_t localOffset_;                ///< additional local offset (for this thread)
    size_t border_;                     ///< border of the volume (= necessary radius for feature extraction)

    size_t numFeatures_;                ///< size of a single feature vector
    float* currentBatch_;               ///< current (global) batch where the results have to be copied to

    const VolumeRAM* volume_;           ///< input volume
    std::vector<std::pair<size_t, float*> > filterKernels_;
    float* dctCoefficients_;
    tgt::bvec4 necessaryScales_;
    tgt::vec2 minMaxIntensity_;
};

} // namespace

#endif 
