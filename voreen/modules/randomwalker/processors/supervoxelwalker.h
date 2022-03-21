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

#ifndef VRN_SUPERVOXELWALKER_H
#define VRN_SUPERVOXELWALKER_H

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

#include "../util/memprofiler.h"

#include <string>
#include <chrono>
#include <unordered_set>

namespace voreen {

class Volume;
class RandomWalkerSolver;
class RandomWalkerSeeds;
class RandomWalkerWeights;
class OctreeBrickPoolManagerMmap;

typedef uint32_t SuperVoxelID;

struct SuperVoxelWalkerPreprocessingResult {
public:
    SuperVoxelWalkerPreprocessingResult(SuperVoxelWalkerPreprocessingResult&& other) = default;
    SuperVoxelWalkerPreprocessingResult& operator=(SuperVoxelWalkerPreprocessingResult&& other) = default;
    SuperVoxelWalkerPreprocessingResult(std::vector<std::vector<SuperVoxelID>>&& edges, std::vector<float>&& regionMeans, VolumeAtomic<SuperVoxelID>&& labels);
    ~SuperVoxelWalkerPreprocessingResult();

    std::vector<std::vector<SuperVoxelID>> edges_;
    std::vector<float> regionMeans_;
    VolumeAtomic<SuperVoxelID> labels_;
    size_t maxConnectivity_;
};

struct SuperVoxelWalkerInput {
    const SuperVoxelWalkerPreprocessingResult* preprocessingResult_; //optional
    const std::vector<float>* previousSolution_; //optional
    const VolumeBase& volume_;
    const VolumeRAM& volram_;
    std::vector<PortDataPointer<Geometry>> foregroundGeomSeeds_;
    std::vector<PortDataPointer<Geometry>> backgroundGeomSeeds_;
    float minWeight_;
    float beta_;
    float errorThreshold_;
    int maxIterations_;
    int maxSuperVoxelVolume_;
    float maxSuperVoxelIntensityDifference_;
    bool generateDebugVolume_;
    bool useMagmaSolver_;
};

struct SuperVoxelWalkerOutput {
    ~SuperVoxelWalkerOutput();
    SuperVoxelWalkerOutput(SuperVoxelWalkerOutput&& other);
    SuperVoxelWalkerOutput(std::unique_ptr<VolumeBase>&& result, std::unique_ptr<VolumeBase>&& resultSuperVoxels, boost::optional<SuperVoxelWalkerPreprocessingResult>&& preprocessingResult, std::vector<float> previousSolution, std::chrono::duration<float> duration);

    std::unique_ptr<VolumeBase> result_;
    std::unique_ptr<VolumeBase> resultSuperVoxels_;
    boost::optional<SuperVoxelWalkerPreprocessingResult> preprocessingResult_;
    std::vector<float> previousSolution_;
    std::chrono::duration<float> duration_;
    bool movedOut_;
};

class SuperVoxelWalker : public AsyncComputeProcessor<SuperVoxelWalkerInput, SuperVoxelWalkerOutput> {
public:
    SuperVoxelWalker();
    virtual ~SuperVoxelWalker();
    virtual Processor* create() const;

    virtual std::string getCategory() const             { return "Volume Processing"; }
    virtual std::string getClassName() const            { return "SuperVoxelWalker";      }
    virtual Processor::CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL;  }

    virtual bool isReady() const;

    static const std::string loggerCat_; ///< category used in logging

protected:
    virtual void setDescriptions() {
        setDescription("Implements the super voxel based random walker method by Fabijańska and Gocławski.");
    }

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

    virtual void initialize();
    virtual void deinitialize();

    void clearPreviousResults();

private:
    VolumePort inportVolume_;
    GeometryPort inportForegroundSeeds_;
    GeometryPort inportBackgroundSeeds_;
    VolumePort outportProbabilities_;
    VolumePort outportSuperVoxels_;

    IntProperty minEdgeWeight_;
    IntProperty beta_;
    IntProperty errorThreshold_;
    IntProperty maxIterations_;
    StringOptionProperty conjGradImplementation_;
    IntProperty maxSuperVoxelVolume_;
    FloatProperty maxSuperVoxelIntensityDifference_;
    BoolProperty generateDebugVolume_;

    ButtonProperty clearResult_;

    ProfileDataCollector ramProfiler_;
    ProfileDataCollector vramProfiler_;

    boost::optional<SuperVoxelWalkerPreprocessingResult> preprocessingResult_;
    boost::optional<std::vector<float>> previousSolution_;

    // Clock and duration used for time keeping
    typedef std::chrono::steady_clock clock;
};

} //namespace

#endif
