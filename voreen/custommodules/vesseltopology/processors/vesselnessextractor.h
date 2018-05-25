#ifndef VRN_VESSELNESS_EXTRACTOR_H
#define VRN_VESSELNESS_EXTRACTOR_H

#include "custommodules/bigdataimageprocessing/processors/largedataprocessor.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/ports/volumeport.h"

#include "custommodules/bigdataimageprocessing/volumefiltering/slicereader.h"
#include "custommodules/bigdataimageprocessing/volumefiltering/volumefilter.h"
#include "custommodules/bigdataimageprocessing/volumefiltering/parallelvolumefilter.h"

namespace voreen {

// Used to extract "vesselness"-features from input volumes.
// Currently live switching between vesselness measures is not supported,
// so that only directional uniformity aware Sato-vesselness can be extracted.
class VesselnessExtractor : public LargeDataProcessor {
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

    virtual void process();
    virtual void adjustPropertiesToInput();
    void updateSmoothingProperties();

    std::unique_ptr<SliceReader> buildStack(const tgt::vec3& standardDeviationVec, const std::string& baseType);
    tgt::vec3 getStandardDeviationForStep(int step) const;

private:
    VolumePort inport_;
    VolumePort outport_;

    FileDialogProperty outputVolumeFilePath_;
    FloatIntervalProperty vesselRadiusRangeRW_;
    IntProperty scaleSpaceSteps_;
    FloatVec3Property minStandardDeviationVec_;
    FloatVec3Property maxStandardDeviationVec_;
    IntVec3Property minSmoothingKernelSize_;
    IntVec3Property maxSmoothingKernelSize_;
    FloatProperty blobRejectorWeight_;
    FloatProperty planeRejectorWeight_;
    FloatProperty intensityThreshold_;
    // TODO: Optionproperty<?> vesselType (>, <, egal)

    static const std::string loggerCat_;
};

}   //namespace


#endif //VRN_VESSELNESS_EXTRACTOR_H
