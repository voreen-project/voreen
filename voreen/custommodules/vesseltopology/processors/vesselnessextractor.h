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

#include "custommodules/bigdataimageprocessing/volumefiltering/slicereader.h"
#include "custommodules/bigdataimageprocessing/volumefiltering/volumefilter.h"
#include "custommodules/bigdataimageprocessing/volumefiltering/parallelvolumefilter.h"

namespace voreen {

struct VesselnessExtractorInput {
    const VolumeBase& input;
    std::unique_ptr<HDF5FileVolume> output;
    int scaleSpaceSteps;
    tgt::vec3 minStandardDeviationVec;
    tgt::vec3 maxStandardDeviationVec;
    std::string baseType;

    VesselnessExtractorInput(
              const VolumeBase& input
            , std::unique_ptr<HDF5FileVolume>&& output
            , int scaleSpaceSteps
            , tgt::vec3 minStandardDeviationVec
            , tgt::vec3 maxStandardDeviationVec
            , std::string baseType
            )
        : input(input)
        , output(std::move(output))
        , scaleSpaceSteps(scaleSpaceSteps)
        , minStandardDeviationVec(minStandardDeviationVec)
        , maxStandardDeviationVec(maxStandardDeviationVec)
        , baseType(baseType)
    {
    }

    VesselnessExtractorInput(const VesselnessExtractorInput&) = delete;
    VesselnessExtractorInput(VesselnessExtractorInput&& old)
        : input(old.input)
        , output(std::move(old.output))
        , scaleSpaceSteps(old.scaleSpaceSteps)
        , minStandardDeviationVec(old.minStandardDeviationVec)
        , maxStandardDeviationVec(old.maxStandardDeviationVec)
        , baseType(old.baseType)
    {
    }

    tgt::vec3 getStandardDeviationForStep(int step) const {
        float alpha = static_cast<float>(step)/(scaleSpaceSteps-1);
        return minStandardDeviationVec*(1.0f - alpha) + maxStandardDeviationVec*alpha;
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
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;   }
    virtual bool usesExpensiveComputation() const { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Transforms in input volume into an output volume containing the voxel wise vesselness of the input.");
    }

    virtual bool isReady() const;
    virtual bool isEndProcessor() const       { return true; }

    virtual VesselnessExtractorInput prepareComputeInput();
    virtual VesselnessExtractorOutput compute(VesselnessExtractorInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(VesselnessExtractorOutput output);

    virtual void adjustPropertiesToInput();
    void updateSmoothingProperties();

private:
    VolumePort inport_;
    VolumePort outport_;

    TempPathProperty outputVolumeFilePath_;
    FloatIntervalProperty vesselRadiusRangeRW_;
    IntProperty scaleSpaceSteps_;
    FloatVec3Property minStandardDeviationVec_;
    FloatVec3Property maxStandardDeviationVec_;
    IntVec3Property minSmoothingKernelSize_;
    IntVec3Property maxSmoothingKernelSize_;
    // TODO: Optionproperty<?> vesselType (>, <, egal)

    static const std::string loggerCat_;
};

}   //namespace


#endif //VRN_VESSELNESS_EXTRACTOR_H
