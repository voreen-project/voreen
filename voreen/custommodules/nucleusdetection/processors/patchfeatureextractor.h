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

#ifndef VRN_PATCHFEATUREEXTRACTOR_H
#define VRN_PATCHFEATUREEXTRACTOR_H

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/volumeinfoproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"

#include "voreen/core/properties/optionproperty.h"

namespace voreen {

/**
 * Abstract base class for processors that perform an extraction of local features.
 * Can be used for classification and / or segmentation
 */
class VRN_CORE_API PatchFeatureExtractor : public VolumeProcessor {
public:
    PatchFeatureExtractor();

    virtual std::string getClassName() const  { return "PatchFeatureExtractor";      }
    virtual std::string getCategory() const   { return "Classification";          }

    virtual bool isReady() const;

protected:

    friend class FeatureExtractionThread;
    friend class TrainingDataExtractionThread;

    /// enum for determining, on which scales the features are computed
    enum SamplingScales {
        SCALE_NONE =        0,      // do not compute this feature
        SCALE_ORIGINAL =    1,      // compute feature in original resolution
        SCALE_SECOND =      1 << 1, // compute feature using every second voxel
        SCALE_THIRD =       1 << 2, // compute feature using every third voxel
        SCALE_FOURTH =      1 << 3  // compute feature using every fourth voxel
    };

    /**
     * struct for specifying a single feature vector
     */
    struct Sample {
        Sample(size_t numFeatures) {
            features_.resize(numFeatures);
        }

        std::vector<float> features_;   // actual features of the sample
    };

    /// deactivates all settings (i.e., sets them read-only)
    virtual void deactivateProperties();

    /// activates all properties (i.e., removes read-only flag)
    virtual void reactivateProperties();

    virtual size_t computeNumFeatures(std::vector<std::pair<size_t,float*> > filterKernels) const;

    virtual tgt::bvec4 determineNecessaryScales() const;

    virtual size_t computeNecessaryBorder(tgt::bvec4 necessaryScales) const;

    virtual std::vector<float*> extractScalePatches(const VolumeRAM* vol, tgt::svec3 voxelPos, tgt::vec2 scalingMinMax, tgt::bvec4 necessaryScales) const;

    virtual void freePatchMemory(std::vector<float*> scales) const;

    virtual void computeFeaturesForVoxel(Sample& sample, const VolumeRAM* volume, tgt::svec3 voxelPos, std::vector<std::pair<size_t, float*> > filterKernels, float* dctCoefficients, tgt::bvec4 necessaryScales, tgt::vec2 minMaxIntensity) const;

    virtual std::vector<std::pair<size_t,float*> > getFilterKernels(size_t patchSize);

    virtual void deleteFilterKernels(std::vector<std::pair<size_t,float*> > kernels);

    virtual bool checkForRequiredFilterKernels() const;

    virtual tgt::vec2 getScalingMinMaxIntensity(const VolumeBase* volume) const;

    virtual float scaleIntensity(float value, tgt::vec2 minMaxValue) const;


    /// Fills the option property with the entries for all scales
    void fillSamplingScaleOptions(OptionProperty<int>& prop);

    // input: volume data, and k-means filter kernels
    //VolumeListPort dataInport_;
    //VolumeListPort labelInport_;
    VolumeListPort filterInport1_;  ///< filter kernel input for original resolution
    VolumeListPort filterInport2_;  ///< filter kernel input for half resolution
    VolumeListPort filterInport3_;  ///< filter kernel input for patches that use every third voxel 
    VolumeListPort filterInport4_;  ///< filter kernel 

    OptionProperty<int> patchSize_; ///< edge length of 3D patch for local feature extraction (the same on every scale)
    OptionProperty<int> normalizePatches_; ///< normalize patch intensities locally, i.e. set the patch mean to 0 and patch standard deviation to 1?
 
    OptionProperty<int> rawPatch2D_;
    OptionProperty<int> rawPatch3D_;
    OptionProperty<int> filteredPatch3D_;
    OptionProperty<int> cosineTransform1D_;
    OptionProperty<int> cosineTransform2D_;
    OptionProperty<int> intensityMean_;
    OptionProperty<int> intensityStandardDeviation_;

private:
    static const std::string loggerCat_;
    static const std::vector<int> scales_;

    ProgressProperty progressProperty_;
};

} // namespace
#endif 
