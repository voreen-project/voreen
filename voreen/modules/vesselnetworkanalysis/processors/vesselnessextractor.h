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

#ifndef VRN_VESSELNESS_EXTRACTOR_H
#define VRN_VESSELNESS_EXTRACTOR_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/temppathproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/ports/volumeport.h"

#include "modules/bigdataimageprocessing/volumefiltering/slicereader.h"
#include "modules/bigdataimageprocessing/volumefiltering/volumefilter.h"
#include "modules/bigdataimageprocessing/volumefiltering/parallelvolumefilter.h"

namespace voreen {

struct VesselnessExtractorInput {
    const VolumeBase& input;
    std::unique_ptr<HDF5FileVolume> output;
    int scaleSpaceSteps;
    tgt::vec3 minStandardDeviationVec;
    tgt::vec3 maxStandardDeviationVec;

    VesselnessExtractorInput(
              const VolumeBase& input
            , std::unique_ptr<HDF5FileVolume>&& output
            , int scaleSpaceSteps
            , tgt::vec3 minStandardDeviationVec
            , tgt::vec3 maxStandardDeviationVec
            )
        : input(input)
        , output(std::move(output))
        , scaleSpaceSteps(scaleSpaceSteps)
        , minStandardDeviationVec(minStandardDeviationVec)
        , maxStandardDeviationVec(maxStandardDeviationVec)
    {
    }

    VesselnessExtractorInput(const VesselnessExtractorInput&) = delete;
    VesselnessExtractorInput(VesselnessExtractorInput&& old)
        : input(old.input)
        , output(std::move(old.output))
        , scaleSpaceSteps(old.scaleSpaceSteps)
        , minStandardDeviationVec(old.minStandardDeviationVec)
        , maxStandardDeviationVec(old.maxStandardDeviationVec)
    {
    }

    tgt::vec3 getStandardDeviationForStep(int step) const {
        float alpha = static_cast<float>(step)/(scaleSpaceSteps-1);
        tgt::vec3 minLog(
                std::log(minStandardDeviationVec.x),
                std::log(minStandardDeviationVec.y),
                std::log(minStandardDeviationVec.z));
        tgt::vec3 maxLog(
                std::log(maxStandardDeviationVec.x),
                std::log(maxStandardDeviationVec.y),
                std::log(maxStandardDeviationVec.z));

        tgt::vec3 interLog = minLog*(1.0f - alpha) + maxLog*alpha;
        return tgt::vec3(
                std::exp(interLog.x),
                std::exp(interLog.y),
                std::exp(interLog.z));
    }
};

typedef std::string VesselnessExtractorOutput;

// Used to extract "vesselness"-features from input volumes.
// Currently live switching between vesselness measures is not supported,
// so that only directional uniformity aware Sato-vesselness can be extracted.
class VesselnessExtractor : public AsyncComputeProcessor<VesselnessExtractorInput, VesselnessExtractorOutput> {
public:
    VesselnessExtractor();
    virtual ~VesselnessExtractor();
    virtual Processor* create() const;

    virtual std::string getClassName() const      { return "VesselnessExtractor";  }
    virtual std::string getCategory() const       { return "Volume Processing"; }
    virtual CodeState getCodeState() const        { return CODE_STATE_TESTING;   }
    virtual bool usesExpensiveComputation() const { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Transforms in input volume into an output volume containing the voxel wise vesselness of the input. "
                "The implementation is based on \"Three-dimensional multi-scale line filter for segmentation andvisualization of curvilinear structures in medical images\" by Sato et al.");
        vesselRadiusRangeRW_.setDescription("The vessel radius range in mm. "
                "Correct specification of the vessel radius is essential for good results. "
                "Consider using a <b>SliceViewer</b> to visually gauge vessel radii if not known from other sources.");
        scaleSpaceSteps_.setDescription("The number of (logarithmic) scale space steps that will be accumulated in the specified vessel radius range. "
                "Higher values generally yield butter results, but require more computation time."
                );
    }

    virtual bool isReady() const;
    virtual bool isEndProcessor() const       { return true; }
    virtual void initialize();

    virtual void process();

    virtual VesselnessExtractorInput prepareComputeInput();
    virtual VesselnessExtractorOutput compute(VesselnessExtractorInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(VesselnessExtractorOutput output);

    virtual void adjustPropertiesToInput();
    void updateSmoothingProperties();

private:
    VolumePort inport_;
    VolumePort outport_;

    BoolProperty enabled_;
    TempPathProperty outputVolumeFilePath_;
    FloatIntervalProperty vesselRadiusRangeRW_;
    IntProperty scaleSpaceSteps_;
    FloatVec3Property minStandardDeviationVec_;
    FloatVec3Property maxStandardDeviationVec_;
    IntVec3Property minSmoothingKernelSize_;
    IntVec3Property maxSmoothingKernelSize_;

    // Disable properties when processor is not enabled:
    PropertyDisabler propertyDisabler_;

    static const std::string loggerCat_;
};

}   //namespace


#endif //VRN_VESSELNESS_EXTRACTOR_H
