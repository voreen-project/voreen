/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_VESSELGRAPHCREATOR_H
#define VRN_VESSELGRAPHCREATOR_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/datastructures/volume/volumelist.h"
#include "voreen/core/properties/filedialogproperty.h"

#include "voreen/core/ports/geometryport.h"
#include "voreen/core/ports/genericport.h"

#include "../ports/vesselgraphport.h"
#include "../ports/vesselgraphlistport.h"

namespace voreen {

struct VesselGraphCreatorInput {
    const VolumeBase& segmentation;
    const VolumeBase* sampleMask;
    std::vector<tgt::vec3> fixedForegroundPoints;
    float binarizationThresholdSegmentationNormalized;
    int numRefinementIterations;
    int minVoxelLength;
    float minElongation;
    float minBulgeSize;
    bool saveDebugData;

    VesselGraphCreatorInput(
            const VolumeBase& segmentation,
            const VolumeBase* sampleMask,
            std::vector<tgt::vec3> fixedForegroundPoints,
            float binarizationThresholdSegmentationNormalized,
            int numRefinementIterations,
            int minVoxelLength,
            float minElongation,
            float minBulgeSize,
            bool saveDebugData
            )
        : segmentation(segmentation)
        , sampleMask(sampleMask)
        , fixedForegroundPoints(fixedForegroundPoints)
        , binarizationThresholdSegmentationNormalized(binarizationThresholdSegmentationNormalized)
        , numRefinementIterations(numRefinementIterations)
        , minVoxelLength(minVoxelLength)
        , minElongation(minElongation)
        , minBulgeSize(minBulgeSize)
        , saveDebugData(saveDebugData)
    {
        tgtAssert(!sampleMask || segmentation.getDimensions() == sampleMask->getDimensions(), "Sample mask dimension mismatch");
    }

    VesselGraphCreatorInput(const VesselGraphCreatorInput&) = delete;
    VesselGraphCreatorInput(VesselGraphCreatorInput&& old)
        : segmentation(old.segmentation)
        , sampleMask(old.sampleMask)
        , fixedForegroundPoints(std::move(old.fixedForegroundPoints))
        , binarizationThresholdSegmentationNormalized(old.binarizationThresholdSegmentationNormalized)
        , numRefinementIterations(old.numRefinementIterations)
        , minVoxelLength(old.minVoxelLength)
        , minElongation(old.minElongation)
        , minBulgeSize(old.minBulgeSize)
        , saveDebugData(old.saveDebugData)
    {
    }
};
struct VesselGraphCreatorOutput {
    std::unique_ptr<VesselGraph> graph;
    std::vector<std::unique_ptr<VolumeBase>> generatedVolumes;
    std::unique_ptr<std::vector<VesselGraph>> generatedGraphs;
};

// A processor that extracts a vessel graph from a voxel skeleton volume and annotates
// edges and nodes with features derived from the supplied segmentation.
class VesselGraphCreator : public AsyncComputeProcessor<VesselGraphCreatorInput, VesselGraphCreatorOutput> {
public:
    VesselGraphCreator();
    virtual ~VesselGraphCreator();

    virtual std::string getClassName() const         { return "VesselGraphCreator";      }
    virtual std::string getCategory() const       { return "Volume Processing"; }
    virtual VoreenSerializableObject* create() const;
    virtual void setDescriptions() {
        setDescription(
                "Create a VesselGraph from a binary input volume (and optionally a mask that specifies the sample volume). "
                "This is the original implementation of the algorithm described in \"Scalable Robust Graph and Feature Extraction for Arbitrary Vessel Networks in Volumetric Datasets\" by Drees et al. "
                "The resulting VesselGraph can be rendered using <b>VesselGraphRenderer</b> or exported using <b>VesselGraphGlobalStats</b> or <b>VesselGraphSave</b>. "
                "The centerlines of the graph can be extracted using a <b>VesselGraphCenterlineConverter</b>."
                );
        binarizationThresholdSegmentation_.setDescription("Values above this threshold will be considered foreground, others background. If the input volume is not binary already, this property can therefore be used for thresholding.");
        numRefinementIterations_.setDescription("Maximum number of refinement iterations. Note that this value can generally be set to a very high value as the computation is interrupted automatically once a fixed point is reached, i.e., when the refinement does not make progress anymore.");

        minBulgeSize_.setDescription("Edges with a bulge size below this threshold will be considered for deletion during the refinement. A bulge size of 1.0 roughly corresponds to hemisphere-shaped bulge.");
    }
    virtual CodeState getCodeState() const        { return CODE_STATE_TESTING; }
    virtual bool isReady() const;

    virtual VesselGraphCreatorInput prepareComputeInput();
    virtual VesselGraphCreatorOutput compute(VesselGraphCreatorInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(VesselGraphCreatorOutput output);

protected:
    virtual void adjustPropertiesToInput();

private:

    // Ports
    VolumePort segmentedVolumeInport_;
    VolumePort sampleMaskInport_;
    GeometryPort fixedForegroundPointInport_;
    VesselGraphPort graphOutport_;
    VolumeListPort generatedVolumesOutport_; //For debug purposes
    std::vector<std::unique_ptr<VolumeBase>> lastGeneratedVolumes_; //owns volumes as generated volumes only holds references
    VesselGraphListPort generatedGraphsOutport_; //For debug purposes

    // Binarization
    FloatProperty binarizationThresholdSegmentation_;

    // Refinement
    IntProperty numRefinementIterations_;
    IntProperty minVoxelLength_;
    FloatProperty minElongation_;
    FloatProperty minBulgeSize_;
    BoolProperty saveDebugData_;

    // Info
    StringProperty tmpStorageSizeInfo_;

    static const std::string loggerCat_;
};


} // namespace voreen

#endif // VRN_VESSELGRAPHCREATOR_H
