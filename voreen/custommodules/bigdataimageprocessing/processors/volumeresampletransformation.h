/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_VOLUMERESAMPLETRANSFORMATION_H
#define VRN_VOLUMERESAMPLETRANSFORMATION_H

#include "voreen/core/processors/asynccomputeprocessor.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/matrixproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/temppathproperty.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"
#include "voreen/core/datastructures/volume/volumefactory.h"

#include "modules/hdf5/io/hdf5filevolume.h"

#ifdef VRN_MODULE_OPENMP
    #include "omp.h"
#endif

namespace {
    //Filtering
    template<class O>
    struct LinearFiltering;
    template<class O>
    struct NearestFiltering;

    // Outside volume handling
    struct Clamping;
    struct ConstantValue;
} // Anonymous namespace


namespace voreen {

enum OutsideVolumeHandling {
    CLAMP,
    CONSTANT_VALUE
};

enum FilteringMode {
    NEAREST,
    LINEAR
};

struct VolumeResampleTransformationInput {
    const VolumeBase& input;
    std::unique_ptr<HDF5FileVolume> outputVolume;
    const tgt::mat4 voxelOutToVoxelIn;
    RealWorldMapping inputNormalized2outputNormalized;
    FilteringMode filteringMode;
    OutsideVolumeHandling outsideVolumeHandling;
    float outsideVolumeValue;

    VolumeResampleTransformationInput(
              const VolumeBase& input
            , std::unique_ptr<HDF5FileVolume>&& outputVolume
            , tgt::mat4 voxelOutToVoxelIn
            , RealWorldMapping inputNormalized2outputNormalized
            , FilteringMode filteringMode
            , OutsideVolumeHandling outsideVolumeHandling
            , float outsideVolumeValue
            )
        : input(input)
        , outputVolume(std::move(outputVolume))
        , voxelOutToVoxelIn(voxelOutToVoxelIn)
        , inputNormalized2outputNormalized(inputNormalized2outputNormalized)
        , filteringMode(filteringMode)
        , outsideVolumeHandling(outsideVolumeHandling)
        , outsideVolumeValue(outsideVolumeValue)
    {
    }

    VolumeResampleTransformationInput(const VolumeResampleTransformationInput&) = delete;
    VolumeResampleTransformationInput(VolumeResampleTransformationInput&& old)
        : input(old.input)
        , outputVolume(std::move(old.outputVolume))
        , voxelOutToVoxelIn(old.voxelOutToVoxelIn)
        , inputNormalized2outputNormalized(old.inputNormalized2outputNormalized)
        , filteringMode(old.filteringMode)
        , outsideVolumeHandling(old.outsideVolumeHandling)
        , outsideVolumeValue(old.outsideVolumeValue)
    {
    }
};
struct VolumeResampleTransformationOutput {
    std::string outputVolumeFilePath;
};

/**
 * Resizes the input volume to the specified dimensions
 * by using a selectable filtering mode.
 */
class VRN_CORE_API VolumeResampleTransformation : public AsyncComputeProcessor<VolumeResampleTransformationInput, VolumeResampleTransformationOutput> {
public:
    VolumeResampleTransformation();
    Processor* create() const;

    std::string getClassName() const      { return "VolumeResampleTransformation";    }
    std::string getCategory() const       { return "Volume Processing"; }
    CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;   }
    bool usesExpensiveComputation() const { return true; }
    virtual bool isReady() const;

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

protected:
    void setDescriptions() {
        setDescription(
                "Transforms the input volume and resamples it so that the resulting voxel grid is aligned with the voxel grid of the input volume. "
                "The transformation is applied in world space. "
                );
    }

    void adjustPropertiesToInput();

private:
    void recalculateOutputProperties();
    void updateSpacingProperties();

    // Calculate properties of a hypothetical output volume given the input volume
    // and the current processor configuration.
    tgt::Bounds getPhysicalOutputBounds(const VolumeBase& input) const;
    tgt::ivec3 getOutputDimensions(const VolumeBase& input) const;
    tgt::vec3 getOutputOffset(const VolumeBase& input) const;

    template<class S>
    void resample(const VolumeBase& input, std::unique_ptr<HDF5FileVolume> output, const tgt::mat4& voxelOutToVoxelIn, RealWorldMapping inputNormalized2outputNormalized, S samplingStrategy, ProgressReporter& progress) const;

    tgt::vec3 getCurrentSpacing() const;


    VolumePort inport_;
    VolumePort outport_;

    enum SpacingHandling {
        KEEP,
        SET,
        MULTIPLY,
    };

    OptionProperty<SpacingHandling> spacingHandling_;
    FloatVec3Property outputSpacing_;
    FloatMat4Property transformMatrix_;
    OptionProperty<FilteringMode> filteringMode_;
    OptionProperty<OutsideVolumeHandling> outsideVolumeHandling_;
    FloatProperty outsideVolumeValue_; //Only visible if CONSTANT_VALUE is selected

    TempPathProperty outputVolumeFilePath_;

    /// Read-only properties displaying the data size of the resampled volume in MB/voxels
    IntVec3Property outputDimensions_;
    StringProperty outputSizeMB_;

    static const std::string loggerCat_; ///< category used in logging
};

// ------------------------------------
// Implementation ---------------------
// ------------------------------------

template<class S>
void VolumeResampleTransformation::resample(const VolumeBase& inputVol, std::unique_ptr<HDF5FileVolume> output, const tgt::mat4& voxelOutToVoxelIn, RealWorldMapping inputNormalized2outputNormalized, S samplingStrategy, ProgressReporter& progress) const {
    tgt::svec3 outputDim = output->getDimensions();
    tgt::ivec3 inputMaxVoxelPos(tgt::ivec3(inputVol.getDimensions()) - tgt::ivec3(1));


    std::string baseType = output->getBaseType();

    const VolumeRAM* input = inputVol.getRepresentation<VolumeRAM>();
    tgtAssert(input, "Could not create volumeram");

    ThreadedTaskProgressReporter parallelProgress(progress, outputDim.z);

    bool aborted = false;
#ifdef VRN_MODULE_OPENMP
        #pragma omp parallel for
#endif
    for(size_t z=0; z<outputDim.z; ++z) {

#ifdef VRN_MODULE_OPENMP
        if(aborted) {
            continue;
        }
#endif

        std::unique_ptr<VolumeRAM> outputSlice(VolumeFactory().create(baseType, tgt::svec3(outputDim.xy(), 1)));

        for (long y = 0; y<static_cast<long>(outputDim.y); ++y) { // MSVC requires y to be signed.
            for(size_t x=0; x<outputDim.x; ++x) {
                tgt::vec4 outputVoxelPos(x,y,z,1);
                tgt::vec4 inputVoxelPos = voxelOutToVoxelIn*outputVoxelPos;
                float inputNormalized = samplingStrategy.sample(*input, inputVoxelPos.xyz());
                float outputNormalized = inputNormalized2outputNormalized.realWorldToNormalized(inputNormalized);

                outputSlice->setVoxelNormalized(outputNormalized, outputVoxelPos.x, outputVoxelPos.y, 0);
            }
        }
        output->writeSlices(outputSlice.get(), z);

        if(parallelProgress.reportStepDone()) {
#ifdef VRN_MODULE_OPENMP
            #pragma omp critical
            {
                aborted = true;
            }
#else
            aborted = true;
            break;
#endif
        }
    }
    if(aborted) {
        throw boost::thread_interrupted();
    }
}

} //namespace

namespace {
    //Filtering
    template<class O>
    struct LinearFiltering {
        LinearFiltering(O outsideVolumeHandler)
            : outsideVolumeHandler_(outsideVolumeHandler)
        {
        }

        float sample(const voreen::VolumeRAM& input, const tgt::vec3& pos) {
            tgt::vec3 wc = pos - tgt::vec3(tgt::ifloor(pos));
            tgt::vec3 wf = tgt::vec3(1.0f) - wc;
            tgt::ivec3 f(tgt::floor(pos));
            tgt::ivec3 c(tgt::iceil(pos));
            float val = 0;
            val += wf.x * wf.y * wf.z * outsideVolumeHandler_.get(input, tgt::ivec3(f.x, f.y, f.z));
            val += wc.x * wf.y * wf.z * outsideVolumeHandler_.get(input, tgt::ivec3(c.x, f.y, f.z));
            val += wf.x * wc.y * wf.z * outsideVolumeHandler_.get(input, tgt::ivec3(f.x, c.y, f.z));
            val += wc.x * wc.y * wf.z * outsideVolumeHandler_.get(input, tgt::ivec3(c.x, c.y, f.z));
            val += wf.x * wf.y * wc.z * outsideVolumeHandler_.get(input, tgt::ivec3(f.x, f.y, c.z));
            val += wc.x * wf.y * wc.z * outsideVolumeHandler_.get(input, tgt::ivec3(c.x, f.y, c.z));
            val += wf.x * wc.y * wc.z * outsideVolumeHandler_.get(input, tgt::ivec3(f.x, c.y, c.z));
            val += wc.x * wc.y * wc.z * outsideVolumeHandler_.get(input, tgt::ivec3(c.x, c.y, c.z));
            return val;
        }

        O outsideVolumeHandler_;
    };

    template<class O>
    struct NearestFiltering {
        NearestFiltering(O outsideVolumeHandler)
            : outsideVolumeHandler_(outsideVolumeHandler)
        {
        }

        float sample(const voreen::VolumeRAM& input, const tgt::vec3& pos) {
            return outsideVolumeHandler_.get(input, tgt::iround(pos));
        }

        O outsideVolumeHandler_;
    };

    // Outside volume handling
    struct Clamping {
        Clamping(tgt::ivec3& inputMaxVoxelPos)
            : inputMaxVoxelPos_(inputMaxVoxelPos)
        {
        }
        float get(const voreen::VolumeRAM& input, const tgt::ivec3& pos) {
            tgt::ivec3 inputSamplePos = tgt::clamp(pos, tgt::ivec3(0,0,0), inputMaxVoxelPos_);
            return input.getVoxelNormalized(inputSamplePos);
        }

        tgt::ivec3 inputMaxVoxelPos_;
    };
    struct ConstantValue {
        ConstantValue(tgt::ivec3& inputMaxVoxelPos, float outsideVolumeValue)
            : inputMaxVoxelPos_(inputMaxVoxelPos)
            , outsideVolumeValue_(outsideVolumeValue)
        {
        }
        float get(const voreen::VolumeRAM& input, const tgt::ivec3& pos) {
            if(   0 <= pos.x && pos.x <= inputMaxVoxelPos_.x
               && 0 <= pos.y && pos.y <= inputMaxVoxelPos_.y
               && 0 <= pos.z && pos.z <= inputMaxVoxelPos_.z) {
                return input.getVoxelNormalized(pos);
            } else {
                return outsideVolumeValue_;
            }
        }

        tgt::ivec3 inputMaxVoxelPos_;
        float outsideVolumeValue_;
    };
} // Anonymous namespace

#endif // VRN_VOLUMERESAMPLETRANSFORMATION_H
