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

#ifndef VRN_OCTREEWALKER_H
#define VRN_OCTREEWALKER_H

#include "voreen/core/processors/asynccomputeprocessor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/temppathproperty.h"
#include "voreen/core/datastructures/octree/volumeoctree.h"

#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "voreen/core/utils/voreenblas/voreenblascpu.h"

#ifdef VRN_MODULE_OPENMP
#include "modules/openmp/include/voreenblasmp.h"
#endif
#ifdef VRN_MODULE_OPENCL
#include "modules/opencl/utils/voreenblascl.h"
#endif

#include "../util/noisemodel.h"

#include <string>
#include <chrono>
#include <unordered_set>

namespace voreen {

class Volume;
class RandomWalkerSolver;
class RandomWalkerSeeds;
class RandomWalkerWeights;
class OctreeBrickPoolManagerMmap;

struct OctreeWalkerPreviousResult {
public:
    VolumeOctree& octree();
    VolumeBase& volume();

    OctreeWalkerPreviousResult(OctreeWalkerPreviousResult&& other);
    OctreeWalkerPreviousResult& operator=(OctreeWalkerPreviousResult&& other);
    OctreeWalkerPreviousResult(VolumeOctree& octree, std::unique_ptr<VolumeBase>&& volume, PointSegmentListGeometryVec3 foregroundSeeds, PointSegmentListGeometryVec3 backgroundSeeds);
    ~OctreeWalkerPreviousResult();

    void destroyButRetainNodes(std::unordered_set<const VolumeOctreeNode*>& nodesToSave) &&;
    bool isPresent() const;

    PointSegmentListGeometryVec3 foregroundSeeds_;
    PointSegmentListGeometryVec3 backgroundSeeds_;
private:
    VolumeOctree* octree_;                  // NEVER owns its own brickpool manager, ALWAYS a representation of previousVolume
    std::unique_ptr<VolumeBase> volume_;    // ALWAYS stores reference to representation in previousOctree_, never null
};

struct OctreeWalkerInput {
    VolumeOctree* previousResult_;
    PointSegmentListGeometryVec3* previousForegroundSeeds_;
    PointSegmentListGeometryVec3* previousBackgroundSeeds_;
    std::unique_ptr<OctreeBrickPoolManagerMmap>& brickPoolManager_;
    const VolumeBase& volume_;
    const VolumeOctree& octree_;
    std::vector<PortDataPointer<Geometry>> foregroundGeomSeeds_;
    std::vector<PortDataPointer<Geometry>> backgroundGeomSeeds_;
    int minWeight_;
    int betaBias_;
    const VoreenBlas* blas_;
    VoreenBlas::ConjGradPreconditioner precond_;
    float errorThreshold_;
    int maxIterations_;
    float homogeneityThreshold_;
    float incrementalSimilarityThreshold_;
    RWNoiseModel noiseModel_;
};

struct OctreeWalkerOutput {
    OctreeWalkerOutput(
        OctreeWalkerPreviousResult&& result,
        std::unordered_set<const VolumeOctreeNode*>&& sharedNodes,
        std::chrono::duration<float> duration
    );
    ~OctreeWalkerOutput();
    OctreeWalkerOutput(OctreeWalkerOutput&& other);

    OctreeWalkerPreviousResult result_;
    std::unordered_set<const VolumeOctreeNode*> sharedNodes_;
    std::chrono::duration<float> duration_;
    bool movedOut_;
};

class OctreeWalker : public AsyncComputeProcessor<OctreeWalkerInput, OctreeWalkerOutput> {
public:
    OctreeWalker();
    virtual ~OctreeWalker();
    virtual Processor* create() const;

    virtual std::string getCategory() const             { return "Volume Processing"; }
    virtual std::string getClassName() const            { return "OctreeWalker";      }
    virtual Processor::CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL;  }

    virtual bool isReady() const;

    static const std::string loggerCat_; ///< category used in logging

protected:
    virtual void setDescriptions() {
        setDescription("Performs a semi-automatic octree volume segmentation using a hierarchical 3D random walker algorithm.");
    }

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

    virtual void initialize();
    virtual void deinitialize();

    void clearPreviousResults();

private:
    const VoreenBlas* getVoreenBlasFromProperties() const;

    VolumePort inportVolume_;
    GeometryPort inportForegroundSeeds_;
    GeometryPort inportBackgroundSeeds_;
    VolumePort outportProbabilities_;

    OptionProperty<RWNoiseModel> noiseModel_;
    IntProperty minEdgeWeight_;
    IntProperty betaBias_;
    StringOptionProperty preconditioner_;
    IntProperty errorThreshold_;
    IntProperty maxIterations_;
    StringOptionProperty conjGradImplementation_;
    FloatProperty homogeneityThreshold_;
    FloatProperty incrementalSimilarityThreshold_;

    VoreenBlasCPU voreenBlasCPU_;
#ifdef VRN_MODULE_OPENMP
    VoreenBlasMP voreenBlasMP_;
#endif
#ifdef VRN_MODULE_OPENCL
    VoreenBlasCL voreenBlasCL_;
#endif

    ButtonProperty clearResult_;
    TempPathProperty resultPath_;
    std::string prevResultPath_;

    boost::optional<OctreeWalkerPreviousResult> previousResult_;
    std::unique_ptr<OctreeBrickPoolManagerMmap> brickPoolManager_;

    // Clock and duration used for time keeping
    typedef std::chrono::steady_clock clock;
};

} //namespace

#endif
