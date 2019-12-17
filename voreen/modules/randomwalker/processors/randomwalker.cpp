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

#include "randomwalker.h"

#include "../solver/randomwalkersolver.h"
#include "../solver/randomwalkerseeds.h"
#include "../solver/randomwalkerweights.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatormorphology.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorresample.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatornumsignificant.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "tgt/vector.h"

#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"

#include <climits>

namespace voreen {

const std::string RandomWalker::loggerCat_("voreen.RandomWalker.RandomWalker");
using tgt::vec3;

RandomWalker::RandomWalker()
#ifdef VRN_MODULE_OPENCL
    : cl::OpenCLProcessor<AsyncComputeProcessor<RandomWalkerInput, RandomWalkerOutput>>(),
#else
    : AsyncComputeProcessor<RandomWalkerInput, RandomWalkerOutput>(),
#endif
    inportVolume_(Port::INPORT, "volume.input"),
    inportForegroundSeeds_(Port::INPORT, "geometry.seedsForeground", "geometry.seedsForeground", true),
    inportBackgroundSeeds_(Port::INPORT, "geometry.seedsBackground", "geometry.seedsBackground", true),
    inportForegroundSeedsVolume_(Port::INPORT, "volume.seedsForeground"),
    inportBackgroundSeedsVolume_(Port::INPORT, "volume.seedsBackground"),
    outportSegmentation_(Port::OUTPORT, "volume.segmentation", "volume.segmentation", false),
    outportProbabilities_(Port::OUTPORT, "volume.probabilities", "volume.probabilities", false),
    outportEdgeWeights_(Port::OUTPORT, "volume.edgeweights", "volume.edgeweights", false),
    usePrevProbAsInitialization_("usePrevProbAsInitialization", "Use Previous Probabilities as Initialization", false, Processor::VALID, Property::LOD_ADVANCED),
    beta_("beta", "Edge Weight Scale: 2^beta", 12, 0, 20),
    minEdgeWeight_("minEdgeWeight", "Min Edge Weight: 10^(-t)", 5, 0, 10),
    preconditioner_("preconditioner", "Preconditioner"),
    errorThreshold_("errorThreshold", "Error Threshold: 10^(-t)", 2, 0, 10),
    maxIterations_("conjGradIterations", "Max Iterations", 1000, 1, 5000),
    conjGradImplementation_("conjGradImplementation", "Implementation"),
    enableLevelOfDetail_("enableLevelOfDetail", "Enable", false),
    lodMinLevel_("lodMinLevel", "Min Level", 0, 0, 5),
    lodMaxLevel_("lodMaxLevel", "Max Level", 2, 0, 5),
    lodForegroundSeedThresh_("lodForegroundSeedThresh", "Foreground Seed Threshold", 0.99f, 0.5f, 1.f),
    lodBackgroundSeedThresh_("lodBackgroundSeedThresh", "Background Seed Threshold", 0.01f, 0.0f, 5.f),
    lodSeedErosionKernelSize_("lodSeedErosionKernelSize", "Seed Erosion Kernel Size"),
    lodMinResolution_("lodMinResolution", "Min Resolution", tgt::ivec3(0), tgt::ivec3(0), tgt::ivec3(INT_MAX)),
    lodMaxResolution_("lodMaxResolution", "Max Resolution", tgt::ivec3(0), tgt::ivec3(0), tgt::ivec3(INT_MAX)),
    enableClipping_("enableClipping", "Enable Clipping", false),
    clipRegion_("clipRegion", "Clip Region", tgt::IntBounds(tgt::ivec3(0), tgt::ivec3(1)), tgt::ivec3(0), tgt::ivec3(10000)),
    enableTransFunc_("enableTransFunc", "Enable TransFunc", false),
    edgeWeightTransFunc_("edgeWeightTransFunc", "Edge Weight TransFunc"),
    edgeWeightBalance_("edgeWeightBalance", "Edge Weight Balance", 0.3f, 0.f, 1.f),
    foregroundThreshold_("foregroundThreshold", "Foreground Threshold", 0.5f, 0.f, 1.f),
    resampleOutputVolumes_("resampleOutputVolumes", "Resample to Input Dimensions", true),
    currentInputVolume_(0)
{
    // ports
    addPort(inportVolume_);
    addPort(inportForegroundSeeds_);
    addPort(inportBackgroundSeeds_);
    addPort(inportForegroundSeedsVolume_);
    addPort(inportBackgroundSeedsVolume_);
    addPort(outportSegmentation_);
    addPort(outportProbabilities_);
    addPort(outportEdgeWeights_);

    addProperty(usePrevProbAsInitialization_);

    // random walker properties
    addProperty(beta_);
    addProperty(minEdgeWeight_);
    beta_.setGroupID("rwparam");
    minEdgeWeight_.setGroupID("rwparam");
    setPropertyGroupGuiName("rwparam", "Random Walker Parametrization");

    // level of detail
    addProperty(enableLevelOfDetail_);
    addProperty(lodMinLevel_);
    addProperty(lodMaxLevel_);
    addProperty(lodForegroundSeedThresh_);
    addProperty(lodBackgroundSeedThresh_);
    lodSeedErosionKernelSize_.addOption("3",  "3x3x3",    3);
    lodSeedErosionKernelSize_.addOption("5",  "5x5x5",    5);
    lodSeedErosionKernelSize_.addOption("7",  "7x7x7",    7);
    lodSeedErosionKernelSize_.addOption("9",  "9x9x9",    9);
    lodSeedErosionKernelSize_.addOption("15", "15x15x15", 15);
    lodSeedErosionKernelSize_.addOption("25", "25x25x25", 25);
    lodSeedErosionKernelSize_.addOption("35", "35x35x35", 35);
    lodSeedErosionKernelSize_.addOption("45", "45x45x45", 45);
    addProperty(lodSeedErosionKernelSize_);
    lodMinResolution_.setReadOnlyFlag(true);
    lodMaxResolution_.setReadOnlyFlag(true);
    addProperty(lodMinResolution_);
    addProperty(lodMaxResolution_);
    enableLevelOfDetail_.setGroupID("levelOfDetail");
    lodMinLevel_.setGroupID("levelOfDetail");
    lodMaxLevel_.setGroupID("levelOfDetail");
    lodForegroundSeedThresh_.setGroupID("levelOfDetail");
    lodBackgroundSeedThresh_.setGroupID("levelOfDetail");
    lodSeedErosionKernelSize_.setGroupID("levelOfDetail");
    lodMinResolution_.setGroupID("levelOfDetail");
    lodMaxResolution_.setGroupID("levelOfDetail");
    setPropertyGroupGuiName("levelOfDetail", "Level of Detail");
    enableLevelOfDetail_.onChange(MemberFunctionCallback<RandomWalker>(this, &RandomWalker::updateGuiState));
    lodMinLevel_.onChange(MemberFunctionCallback<RandomWalker>(this, &RandomWalker::lodMinLevelChanged));
    lodMaxLevel_.onChange(MemberFunctionCallback<RandomWalker>(this, &RandomWalker::lodMaxLevelChanged));

    // conjugate gradient solver
    preconditioner_.addOption("none", "None");
    preconditioner_.addOption("jacobi", "Jacobi");
    preconditioner_.select("jacobi");
    addProperty(preconditioner_);
    addProperty(errorThreshold_);
    addProperty(maxIterations_);
    conjGradImplementation_.addOption("blasCPU", "CPU");
#ifdef VRN_MODULE_OPENMP
    conjGradImplementation_.addOption("blasMP", "OpenMP");
    conjGradImplementation_.select("blasMP");
#endif
#ifdef VRN_MODULE_OPENCL
    conjGradImplementation_.addOption("blasCL", "OpenCL");
    conjGradImplementation_.select("blasCL");
#endif
    addProperty(conjGradImplementation_);
    preconditioner_.setGroupID("conjGrad");
    errorThreshold_.setGroupID("conjGrad");
    maxIterations_.setGroupID("conjGrad");
    conjGradImplementation_.setGroupID("conjGrad");
    setPropertyGroupGuiName("conjGrad", "Conjugate Gradient Solver");

    // clipping
    addProperty(enableClipping_);
    addProperty(clipRegion_);
    enableClipping_.setGroupID("clipping");
    clipRegion_.setGroupID("clipping");
    setPropertyGroupGuiName("clipping", "Clipping");
    enableClipping_.onChange(MemberFunctionCallback<RandomWalker>(this, &RandomWalker::updateGuiState));

    // transfer functions
    addProperty(enableTransFunc_);
    addProperty(edgeWeightTransFunc_);
    addProperty(edgeWeightBalance_);
    enableTransFunc_.setGroupID("classificationIncorporation");
    edgeWeightBalance_.setGroupID("classificationIncorporation");
    edgeWeightTransFunc_.setGroupID("classificationIncorporation");
    setPropertyGroupGuiName("classificationIncorporation", "Classification");
    enableTransFunc_.onChange(MemberFunctionCallback<RandomWalker>(this, &RandomWalker::updateGuiState));

    // output volumes
    addProperty(foregroundThreshold_);
    foregroundThreshold_.setGroupID("output");
    addProperty(resampleOutputVolumes_);
    resampleOutputVolumes_.setGroupID("output");
    setPropertyGroupGuiName("output", "Output");
}

RandomWalker::~RandomWalker() {
}

Processor* RandomWalker::create() const {
    return new RandomWalker();
}

void RandomWalker::initialize() {
#ifdef VRN_MODULE_OPENCL
    if(OpenCLModule::getInstance()->isInitialized()) {
        cl::OpenCLProcessor<AsyncComputeProcessor>::initialize();
    } else {
        conjGradImplementation_.removeOption("blasCL");
        AsyncComputeProcessor::initialize();
    }
#else
    AsyncComputeProcessor::initialize();
#endif

    updateGuiState();
}

void RandomWalker::deinitialize() {
    // clear lod volumes
    for (size_t i=0; i<lodVolumes_.size(); i++)
        delete lodVolumes_.at(i);
    lodVolumes_.clear();

#ifdef VRN_MODULE_OPENCL
    if(OpenCLModule::getInstance()->isInitialized()) {
        cl::OpenCLProcessor<AsyncComputeProcessor>::deinitialize();
    } else {
        AsyncComputeProcessor::deinitialize();
    }
#else
    AsyncComputeProcessor::deinitialize();
#endif
}

#ifdef VRN_MODULE_OPENCL
void RandomWalker::initializeCL() {
    voreenBlasCL_.initialize();
    invalidate();
}

void RandomWalker::deinitializeCL() {
    interruptComputation();
    voreenBlasCL_.deinitialize();
}

bool RandomWalker::isDeviceChangeSupported() const {
    return true;
}
#endif

bool RandomWalker::isReady() const {
    bool ready = false;
    ready |= outportSegmentation_.isConnected();
    ready |= outportProbabilities_.isConnected();
    ready |= outportEdgeWeights_.isConnected();
    ready &= inportVolume_.isReady();
    ready &= inportForegroundSeeds_.isReady();
    ready &= inportBackgroundSeeds_.isReady();
    return ready;
}

RandomWalker::ComputeInput RandomWalker::prepareComputeInput() {
    edgeWeightTransFunc_.setVolume(inportVolume_.getData());


    tgtAssert(inportVolume_.hasData(), "no input volume");

    // clear previous results and update property ranges, if input volume has changed
    if (inportVolume_.hasChanged()) {
        outportSegmentation_.setData(0);
        outportProbabilities_.setData(0);
        outportEdgeWeights_.setData(0);

        // clear lod volumes
        for (size_t i=0; i<lodVolumes_.size(); i++)
            delete lodVolumes_.at(i);
        lodVolumes_.clear();

        // adjust clip plane properties
        tgt::ivec3 volDim = inportVolume_.getData()->getRepresentation<VolumeRAM>()->getDimensions();
        clipRegion_.setMaxValue(volDim - tgt::ivec3::one);

        lodMinResolution_.setMaxValue(volDim);
        float scaleFactorMin = static_cast<float>(1 << lodMaxLevel_.get());
        tgtAssert(scaleFactorMin >= 1.f, "invalid scale factor");
        tgt::ivec3 minDim = tgt::iround(tgt::vec3(volDim) / scaleFactorMin);
        lodMinResolution_.set(minDim);

        lodMaxResolution_.setMaxValue(volDim);
        float scaleFactorMax = static_cast<float>(1 << lodMinLevel_.get());
        tgtAssert(scaleFactorMax >= 1.f, "invalid scale factor");
        tgt::ivec3 maxDim = tgt::iround(tgt::vec3(volDim) / scaleFactorMax);
        lodMaxResolution_.set(maxDim);

        prevProbabilities_ = std::vector<float>();
    }

    const int startLevel = enableLevelOfDetail_.get() ? lodMaxLevel_.get() : 0;
    const int endLevel = enableLevelOfDetail_.get() ? lodMinLevel_.get() : 0;
    tgtAssert(startLevel-endLevel >= 0, "invalid level range");

    boost::optional<tgt::IntBounds> bounds = boost::none;
    if (enableClipping_.get()) {
        bounds = clipRegion_.get();
    }

    // 2. Edge weight calculator (independent from scale level)
    std::unique_ptr<RandomWalkerWeights> weights(getEdgeWeightsFromProperties());

    // select BLAS implementation and preconditioner
    const VoreenBlas* voreenBlas = getVoreenBlasFromProperties();
    VoreenBlas::ConjGradPreconditioner precond = VoreenBlas::NoPreconditioner;
    if (preconditioner_.isSelected("jacobi"))
        precond = VoreenBlas::Jacobi;


    std::vector<float> prevProbs;
    if(usePrevProbAsInitialization_.get()) {
        prevProbs = std::vector<float>(prevProbabilities_);
    }

    float lodForegroundSeedThresh = lodForegroundSeedThresh_.get();
    float lodBackgroundSeedThresh = lodBackgroundSeedThresh_.get();

    float errorThresh = 1.f / pow(10.f, static_cast<float>(errorThreshold_.get()));
    int maxIterations = maxIterations_.get();
    int lodSeedErosionKernelSize = lodSeedErosionKernelSize_.getValue();

    return ComputeInput {
        inportVolume_.getThreadSafeData(),
        inportForegroundSeeds_.getThreadSafeAllData(),
        inportBackgroundSeeds_.getThreadSafeAllData(),
        inportForegroundSeedsVolume_.getThreadSafeData(),
        inportBackgroundSeedsVolume_.getThreadSafeData(),
        std::move(lodVolumes_),
        prevProbs,
        startLevel,
        endLevel,
        bounds,
        std::move(weights),
        voreenBlas,
        precond,
        lodForegroundSeedThresh,
        lodBackgroundSeedThresh,
        errorThresh,
        maxIterations,
        lodSeedErosionKernelSize,
    };
}

namespace {
    struct LoopRecord {
        int iteration;
        int level;
        float scaleFactor;
        tgt::ivec3 workDim;
        size_t numSeeds;
        size_t numForegroundSeeds;
        size_t numBackgroundSeeds;
        tgt::vec2 probabilityRange;
        int numIterations;
        std::chrono::duration<float> timeIteration;
        std::chrono::duration<float> timeSetup;
        std::chrono::duration<float> timeSolving;
        std::chrono::duration<float> timeSeedAnalysis;

        void print() {
            size_t numVoxels = tgt::hmul(workDim);
            std::string cat = "voreen.RandomWalker.RandomWalker";
            LINFOC(cat, iteration << ". Iteration: level=" << level << ", scaleFactor=" << scaleFactor
                << ", dim=" << workDim);
            LINFOC(cat, "* num voxels: " << numVoxels << ", num seeds:  " << numSeeds << " (ratio: " << (float)numSeeds/tgt::hmul(workDim) << ")");
            LINFOC(cat, "* num unseeded: " << numVoxels-numSeeds);
            LINFOC(cat, "* probability range: " << probabilityRange);
            LINFOC(cat, "* runtime: " << timeIteration.count() << " sec");
            LINFOC(cat, "  - system setup: " << timeSetup.count() << " sec");
            LINFOC(cat, "  - solving: " << timeSolving.count() << " sec"
                << " (iterations: " << numIterations << ")");
            LINFOC(cat, "  - seed analysis: " << timeSeedAnalysis.count() << " sec");
        }
    };
} // namespace anonymous

static void getSeedListsFromPorts(std::vector<PortDataPointer<Geometry>>& geom, PointSegmentListGeometry<tgt::vec3>& seeds) {

    for (size_t i=0; i<geom.size(); i++) {
        const PointSegmentListGeometry<tgt::vec3>* seedList = dynamic_cast<const PointSegmentListGeometry<tgt::vec3>* >(geom.at(i).get());
        if (!seedList)
            LWARNINGC("voreen.RandomWalker.RandomWalker", "Invalid geometry. PointSegmentListGeometry<vec3> expected.");
        else {
            auto transformMat = seedList->getTransformationMatrix();
            for (int j=0; j<seedList->getNumSegments(); j++) {
                std::vector<tgt::vec3> points;
                for(auto& vox : seedList->getSegment(j)) {
                    points.push_back(transformMat.transform(vox));
                }
                seeds.addSegment(points);
            }
        }
    }
}

RandomWalker::ComputeOutput RandomWalker::compute(ComputeInput input, ProgressReporter& progressReporter) const {
    RandomWalkerOutput invalidResult = RandomWalkerOutput {
        std::unique_ptr<RandomWalkerSolver>(nullptr),
        std::chrono::duration<float>(0),
        std::vector<const Volume*>(),
        std::vector<float>(),
    };

    progressReporter.setProgress(0.0);
    SubtaskProgressReporterCollection<3> globalProgressSteps(progressReporter, {0.01, 0.01, 0.98});

    auto start = clock::now();
    tgtAssert(input.inputHandle_, "No input Volume");

    const VolumeBase* inputHandle = input.inputHandle_;
    const VolumeRAM* inputVolume = inputHandle->getRepresentation<VolumeRAM>();

    // create lod volumes, if not present
    bool computeLODs = false;
    for (int level=std::max(input.endLevel_, 1); level<=input.startLevel_ && !computeLODs; level++)
        computeLODs |= ((level-1) >= (int)input.lodVolumes_.size()) || (input.lodVolumes_.at(level-1) == 0);
    if (computeLODs) {
        LINFO("Computing level of detail volumes...");
        auto start = clock::now();
        for (int level=std::max(input.endLevel_, 1); level<=input.startLevel_; level++) {
            while ((level-1) >= (int)input.lodVolumes_.size())
                input.lodVolumes_.push_back(0);
            if (input.lodVolumes_.at(level-1) == 0) {
                float scaleFactor = static_cast<float>(1 << level);
                tgtAssert(scaleFactor >= 1.f, "invalid scale factor");
                tgt::ivec3 levelDim = tgt::iround(tgt::vec3(inputVolume->getDimensions()) / scaleFactor);
                try {
                    Volume* levelVolume = VolumeOperatorResample::APPLY_OP(inputHandle, levelDim, VolumeRAM::LINEAR);
                    input.lodVolumes_.at(level-1) = levelVolume;
                }
                catch (std::bad_alloc& e) {
                    LERROR("Failed create level-" << level << " volume (dim=" << levelDim << ") : bad allocation");
                    return invalidResult;
                }
            }
        }
        auto end = clock::now();
        LINFO("...finished (" << (end-start).count() << " sec)");
    }
    globalProgressSteps.get<0>().setProgress(1.0);

    // work resources
    std::unique_ptr<RandomWalkerSolver> solver = 0;
    const VolumeBase* workVolume = 0;
    std::unique_ptr<VolumeRAM_UInt8> foregroundSeedVol = 0;
    std::unique_ptr<VolumeRAM_UInt8> backgroundSeedVol = 0;

    // input seed lists
    PointSegmentListGeometryVec3 foregroundSeedsPort;
    PointSegmentListGeometryVec3 backgroundSeedsPort;
    getSeedListsFromPorts(input.foregroundGeomSeeds_, foregroundSeedsPort);
    getSeedListsFromPorts(input.backgroundGeomSeeds_, backgroundSeedsPort);

    // input seed volumes
    try {
        if (input.foregroundVolSeeds_) {
            foregroundSeedVol.reset(dynamic_cast<const VolumeRAM_UInt8*>(input.foregroundVolSeeds_->getRepresentation<VolumeRAM>())->clone());

            if (!foregroundSeedVol) {
                VolumeOperatorConvert converter;
                std::unique_ptr<Volume> h(converter.apply<uint8_t>(input.foregroundVolSeeds_));
                foregroundSeedVol.reset(dynamic_cast<VolumeRAM_UInt8*>(h->getWritableRepresentation<VolumeRAM>()));
                h->releaseAllRepresentations();
            }
        }
        if (input.backgroundVolSeeds_) {
            backgroundSeedVol.reset(dynamic_cast<const VolumeRAM_UInt8*>(input.backgroundVolSeeds_->getRepresentation<VolumeRAM>())->clone());

            if (!backgroundSeedVol) {
                VolumeOperatorConvert converter;
                std::unique_ptr<Volume> h(converter.apply<uint8_t>(input.backgroundVolSeeds_));
                backgroundSeedVol.reset(dynamic_cast<VolumeRAM_UInt8*>(h->getWritableRepresentation<VolumeRAM>()));
                h->releaseAllRepresentations();
            }
        }
    }
    catch (std::bad_alloc& e) {
        LERROR("Failed to convert input seed volumes: bad allocation");
        return invalidResult;
    }
    globalProgressSteps.get<1>().setProgress(1.0);

    std::vector<float> newProbabilities;

    //
    // Multi Scale Loop
    //
    auto& loopProgress = globalProgressSteps.get<2>();
    int numSteps = input.startLevel_ - input.endLevel_ + 1;
    int step = 0;
    for (int level = input.startLevel_; level >= input.endLevel_; level--) {
        SubtaskProgressReporter iterationProgress(loopProgress, tgt::vec2(static_cast<float>(step)/numSteps, static_cast<float>(step+1)/numSteps));
        loopProgress.setProgress(0.0);
        ++step;


        LoopRecord loopRecord;
        loopRecord.iteration = input.startLevel_-level+1;
        loopRecord.level = level;
        auto iterationStart = clock::now();

        /*
         * 0. Current scale factor and work volume
         */
        float scaleFactor = static_cast<float>(1 << level);
        tgtAssert(scaleFactor >= 1.f, "invalid scale factor");
        tgt::ivec3 workDim = tgt::iround(tgt::vec3(inputVolume->getDimensions()) / scaleFactor);
        //LINFO("Scale Factor: " << scaleFactor << ", Work dim: " << workDim);
        loopRecord.scaleFactor = scaleFactor;
        loopRecord.workDim = workDim;

        // get current work volume
        if (level == 0)
            workVolume = inputHandle;
        else {
            tgtAssert((level-1) < (int)input.lodVolumes_.size() && input.lodVolumes_.at(level-1), "lod volume missing");
            workVolume = input.lodVolumes_.at(level-1);
        }

        /*
         * 1. Seed points
         */
        // Convert input seed point list according to current scale level
        PointSegmentListGeometryVec3 foregroundSeeds;
        PointSegmentListGeometryVec3 backgroundSeeds;
        for (int s=0; s<foregroundSeedsPort.getNumSegments(); s++) {
            std::vector<tgt::vec3> segment;
            for (size_t i=0; i<foregroundSeedsPort.getSegment(s).size(); i++)
                segment.push_back(foregroundSeedsPort.getSegment(s).at(i) / scaleFactor);
            foregroundSeeds.addSegment(segment);
        }
        for (int s=0; s<backgroundSeedsPort.getNumSegments(); s++) {
            std::vector<tgt::vec3> segment;
            for (size_t i=0; i<backgroundSeedsPort.getSegment(s).size(); i++)
                segment.push_back(backgroundSeedsPort.getSegment(s).at(i) / scaleFactor);
            backgroundSeeds.addSegment(segment);
        }

        // adapted clipping planes
        tgt::ivec3 clipLLF(-1);
        tgt::ivec3 clipURB(-1);
        if (input.clipRegion_) {
            clipLLF = tgt::iround(tgt::vec3(input.clipRegion_->getLLF()) / scaleFactor);
            clipURB = tgt::iround(tgt::vec3(input.clipRegion_->getURB()) / scaleFactor);
        }

        std::unique_ptr<RandomWalkerSeeds> seeds(new RandomWalkerTwoLabelSeeds(foregroundSeeds, backgroundSeeds,
            foregroundSeedVol.get(), backgroundSeedVol.get(), clipLLF, clipURB));

        /*
         * 2. Set up Random Walker system.
         */

        const RandomWalkerSeeds* seedsRef = seeds.get(); // only valid until solver is deleted.
        solver.reset(new RandomWalkerSolver(workVolume, seeds.release(), input.weights_.release()));
        try {
            auto start = clock::now();
            solver->setupEquationSystem();
            auto finish = clock::now();
            //LINFO("...finished: " << std::chrono::duration<float>(finish-start).count() << " sec");
            loopRecord.timeSetup = finish - start;
        }
        catch (tgt::Exception& e) {
            LERROR("Failed to setup Random Walker equation system: " << e.what());
            return invalidResult;
        }
        loopRecord.numSeeds = seedsRef->getNumSeeds();
        if (const RandomWalkerTwoLabelSeeds* twoLabelSeeds = dynamic_cast<const RandomWalkerTwoLabelSeeds*>(seedsRef)) {
            loopRecord.numForegroundSeeds = twoLabelSeeds->getNumForegroundSeeds();
            loopRecord.numBackgroundSeeds = twoLabelSeeds->getNumBackgroundSeeds();
        }

        /*
         * 3. Compute Random Walker solution.
         */

        // solve
        try {
            float* initialization = nullptr;
            std::vector<float> initializationStorage;
            initializationStorage.reserve(solver->getSystemSize());
            if(!input.prevProbabilities_.empty()) {
                tgtAssert(input.prevProbabilities_.size() == solver->getNumVoxels(), "Old and new probalities size missmatch");
                // This is highly depends on the implementation of RandomWalkerSolver!
                for(size_t i = 0; i < input.prevProbabilities_.size(); ++i) {
                    if(!solver->isSeedPoint(i)) {
                        initializationStorage.push_back(input.prevProbabilities_[i]);
                    }
                }
                LINFO("Using existing probalities as initialization");
                initialization = initializationStorage.data();
            }

            auto start = clock::now();
            int iterations = solver->solve(input.blas_, initialization, input.precond_, input.errorThreshold_, input.maxIterations_, loopProgress);
            auto finish = clock::now();
            loopRecord.timeSolving = finish - start;
            loopRecord.numIterations = iterations;

            tgtAssert(solver->getSystemState() == RandomWalkerSolver::Solved, "System not solved");
            newProbabilities = std::vector<float>();
            newProbabilities.reserve(solver->getNumVoxels());
            for(size_t i = 0; i < solver->getNumVoxels(); ++i) {
                newProbabilities.push_back(solver->getProbabilityValue(i));
            }
        }
        catch (VoreenException& e) {
            LERROR("Failed to compute Random Walker solution: " << e.what());
            return invalidResult;
        }

        loopRecord.probabilityRange = solver->getProbabilityRange();

        /*
         * 4. Derive seed volumes for next iteration
         */
        bool lastIteration = (level == input.endLevel_);
        if (!lastIteration) {

            auto start = clock::now();

            // generate probability map from current solution
            std::unique_ptr<VolumeRAM_UInt16> probabilityVolume = 0;
            try {
                probabilityVolume.reset(solver->generateProbabilityVolume<VolumeRAM_UInt16>());
            }
            catch (VoreenException& e) {
                LERROR("computeRandomWalkerSolution() Failed to generate probability volume: " << e.what());
                return invalidResult;
            }

            // allocate seed volumes for next iteration
            try {
                foregroundSeedVol.reset(new VolumeRAM_UInt8(workDim));
                backgroundSeedVol.reset(new VolumeRAM_UInt8(workDim));
            }
            catch (std::bad_alloc& e) {
                LERROR("computeRandomWalkerSolution() Failed to create seed volumes: bad allocation");
                return invalidResult;
            }

            // threshold probability map to derive seeds
            float maxProbValue = probabilityVolume->elementRange().y;
            int foregroundThresh = tgt::iround(input.lodForegroundSeedThresh_*maxProbValue);
            int backgroundThresh = tgt::iround(input.lodBackgroundSeedThresh_*maxProbValue);
            foregroundSeedVol->clear();
            backgroundSeedVol->clear();

            for (size_t i=0; i<probabilityVolume->getNumVoxels(); i++) {
                int probValue = probabilityVolume->voxel(i);
                if (probValue <= backgroundThresh)
                    backgroundSeedVol->voxel(i) = 255;
                else if (probValue >= foregroundThresh)
                    foregroundSeedVol->voxel(i) = 255;
            }

            // erode obtained fore- and background seed volumes
            std::unique_ptr<Volume> erh(VolumeOperatorCubeErosion::APPLY_OP(new Volume(foregroundSeedVol.get(), vec3(1.0f), vec3(0.0f)), input.lodSeedErosionKernelSize_)); //FIXME: small memory leak
            foregroundSeedVol.reset(static_cast<VolumeRAM_UInt8*>(erh->getWritableRepresentation<VolumeRAM>()));
            erh->releaseAllRepresentations();
            erh.reset(VolumeOperatorCubeErosion::APPLY_OP(new Volume(backgroundSeedVol.get(), vec3(1.0f), vec3(0.0f)), input.lodSeedErosionKernelSize_));
            backgroundSeedVol.reset(static_cast<VolumeRAM_UInt8*>(erh->getWritableRepresentation<VolumeRAM>()));
            erh->releaseAllRepresentations();
            for (size_t i=0; i<foregroundSeedVol->getNumVoxels(); i++) {
                    if (foregroundSeedVol->voxel(i) < 255)
                        foregroundSeedVol->voxel(i) = 0;
            }
            for (size_t i=0; i<backgroundSeedVol->getNumVoxels(); i++) {
                if (backgroundSeedVol->voxel(i) < 255)
                    backgroundSeedVol->voxel(i) = 0;
            }

            auto finish = clock::now();
            loopRecord.timeSeedAnalysis = finish - start;
        } // !lastIteration
        else {
            loopRecord.timeSeedAnalysis = std::chrono::duration<float>(0);
        }

        auto iterationEnd = clock::now();
        loopRecord.timeIteration = iterationEnd - iterationStart;

        loopRecord.print();

        loopProgress.setProgress(1.0);
    } // end loop

    tgtAssert(solver, "no random walker solver");

    auto finish = clock::now();
    return ComputeOutput {
        std::move(solver), finish - start, std::move(input.lodVolumes_), newProbabilities,
    };
}
void RandomWalker::processComputeOutput(ComputeOutput output) {
    if (output.solver_ && output.solver_->getSystemState() == RandomWalkerSolver::Solved) {
        prevProbabilities_ = output.newProbabilities_;
        lodVolumes_ = output.lodVolumes_;

        LINFO("Total runtime: " << output.duration_.count() << " sec");

        // put out results
        if (outportSegmentation_.isConnected())
            putOutSegmentation(output.solver_.get());
        else
            outportSegmentation_.setData(0);

        if (outportProbabilities_.isConnected())
            putOutProbabilities(output.solver_.get());
        else
            outportProbabilities_.setData(0);

        if (outportEdgeWeights_.isConnected())
            putOutEdgeWeights(output.solver_.get());
        else
            outportEdgeWeights_.setData(0);
    }
    else {
        LERROR("Failed to compute Random Walker solution");
    }
}

void RandomWalker::putOutSegmentation(const RandomWalkerSolver* solver) {
     tgtAssert(solver, "null pointer passed");
     tgtAssert(inportVolume_.hasData(), "no input data");
     tgtAssert(solver->getSystemState() == RandomWalkerSolver::Solved, "system not solved");
     const VolumeRAM* inputVolume = inportVolume_.getData()->getRepresentation<VolumeRAM>();

    outportSegmentation_.setData(0);

    VolumeRAM_UInt8* segVolume = 0;
    try {
        segVolume = solver->generateBinarySegmentation<VolumeRAM_UInt8>(foregroundThreshold_.get());
        Volume* segHandle = new Volume(segVolume, inportVolume_.getData());
        segHandle->setRealWorldMapping(RealWorldMapping());
        if (segVolume->getDimensions() != inportVolume_.getData()->getDimensions()) {
            tgt::vec3 spacingScale = tgt::vec3(inportVolume_.getData()->getDimensions()) / tgt::vec3(segVolume->getDimensions());
            segHandle->setSpacing(inportVolume_.getData()->getSpacing() * spacingScale);
        }
        size_t numForeground = VolumeOperatorNumSignificant::APPLY_OP(segHandle);
        LINFO("Foreground ratio: " << (float)numForeground / segVolume->getNumVoxels());

        if (resampleOutputVolumes_.get() && segVolume->getDimensions() != inputVolume->getDimensions()) {
                LINFO("Resampling segmentation volume to input volume's dimensions");
                Volume* t = VolumeOperatorResample::APPLY_OP(segHandle, inputVolume->getDimensions(), VolumeRAM::NEAREST);
                delete segHandle;
                segHandle = t;
        }
        outportSegmentation_.setData(segHandle);
    }
    catch (VoreenException& e) {
        LERROR("Failed to generate segmentation volume: " << e.what());
        delete segVolume;
    }
}


RandomWalkerWeights* RandomWalker::getEdgeWeightsFromProperties() const {
    float beta = static_cast<float>(1<<beta_.get());
    float minWeight = 1.f / pow(10.f, static_cast<float>(minEdgeWeight_.get()));
    float tfBlendFactor = edgeWeightBalance_.get();

    if (enableTransFunc_.get())
        return new RandomWalkerWeightsTransFunc(edgeWeightTransFunc_.get(), beta, tfBlendFactor, minWeight);
    else
        return new RandomWalkerWeightsIntensity(beta, minWeight);
}

const VoreenBlas* RandomWalker::getVoreenBlasFromProperties() const {

#ifdef VRN_MODULE_OPENMP
    if (conjGradImplementation_.isSelected("blasMP")) {
        return &voreenBlasMP_;
    }
#endif
#ifdef VRN_MODULE_OPENCL
    if (conjGradImplementation_.isSelected("blasCL")) {
        return &voreenBlasCL_;
    }
#endif

    return &voreenBlasCPU_;
}

void RandomWalker::putOutProbabilities(const RandomWalkerSolver* solver) {
    tgtAssert(solver, "null pointer passed");
    tgtAssert(solver->getSystemState() == RandomWalkerSolver::Solved, "system not solved");
    tgtAssert(inportVolume_.hasData(), "no input data");
    const VolumeRAM* inputVolume = inportVolume_.getData()->getRepresentation<VolumeRAM>();

    outportProbabilities_.setData(0);

    VolumeRAM_UInt16* probabilityVolume = 0;
    try {
        probabilityVolume = solver->generateProbabilityVolume<VolumeRAM_UInt16>();
        Volume* probabilityHandle = new Volume(probabilityVolume, inportVolume_.getData());
        probabilityHandle->setRealWorldMapping(RealWorldMapping());
        if (probabilityVolume->getDimensions() != inportVolume_.getData()->getDimensions()) {
            tgt::vec3 spacingScale = tgt::vec3(inportVolume_.getData()->getDimensions()) / tgt::vec3(probabilityVolume->getDimensions());
            probabilityHandle->setSpacing(inportVolume_.getData()->getSpacing() * spacingScale);
        }

        if (resampleOutputVolumes_.get() && probabilityVolume->getDimensions() != inputVolume->getDimensions()) {
            LINFO("Resampling probability volume to input volume's dimensions");
            Volume* t = VolumeOperatorResample::APPLY_OP(probabilityHandle, inputVolume->getDimensions(), VolumeRAM::NEAREST);
            delete probabilityHandle;
            probabilityHandle = t;
        }

        outportProbabilities_.setData(probabilityHandle);
    }
    catch (VoreenException& e) {
        LERROR("Failed to generate probability volume: " << e.what());
        delete probabilityVolume;
    }
}

void RandomWalker::putOutEdgeWeights(const RandomWalkerSolver* solver) {
    tgtAssert(solver, "null pointer passed");
    tgtAssert(solver->getSystemState() == RandomWalkerSolver::Solved, "system not solved");
    tgtAssert(inportVolume_.hasData(), "no input data");
    const VolumeRAM* inputVolume = inportVolume_.getData()->getRepresentation<VolumeRAM>();

    outportEdgeWeights_.setData(0);

    VolumeRAM_UInt8* edgeWeights = 0;
    try {
        edgeWeights = new VolumeRAM_UInt8(solver->getVolumeDimensions());
    }
    catch (VoreenException& e) {
        LERROR("Failed to generate edge weight volume: " << e.what());
        return;
    }

    const EllpackMatrix<float>& mat = solver->getMatrix();

    float maxWeightSum = 0.f;
    for (size_t row=0; row<solver->getSystemSize(); row++) {
        maxWeightSum = std::max(maxWeightSum, mat.getValue(row, row));
    }

    for (size_t i=0; i<solver->getNumVoxels(); i++) {
        float weight;
        if (solver->isSeedPoint(i)) {
            weight = solver->getSeedValue(i);
        }
        else {
            size_t row = solver->getRowIndex(i); // volIndexToRow[i];
            tgtAssert(row < solver->getSystemSize(), "Invalid row");
            weight = logf(1.f + (mat.getValue(row, row) / maxWeightSum)*1e3f) / logf(1e3f);
            weight = 1.f - weight;
        }

        edgeWeights->voxel(i) = tgt::clamp(tgt::iround(weight * 255.f), 0, 255);
    }

    Volume* ewHandle = new Volume(edgeWeights, inportVolume_.getData());
    ewHandle->setRealWorldMapping(RealWorldMapping());
    if (edgeWeights->getDimensions() != inportVolume_.getData()->getDimensions()) {
        tgt::vec3 spacingScale = tgt::vec3(inportVolume_.getData()->getDimensions()) / tgt::vec3(edgeWeights->getDimensions());
        ewHandle->setSpacing(inportVolume_.getData()->getSpacing() * spacingScale);
    }

    if (resampleOutputVolumes_.get() && edgeWeights->getDimensions() != inputVolume->getDimensions()) {
        try {
            LINFO("Resampling edge weight volume to input volume's dimensions");
            Volume* newHandle = VolumeOperatorResample::APPLY_OP(ewHandle, inputVolume->getDimensions(), VolumeRAM::NEAREST);
            delete ewHandle;
            ewHandle = newHandle;
            edgeWeights = dynamic_cast<VolumeRAM_UInt8*>(ewHandle->getWritableRepresentation<VolumeRAM>());
        }
        catch (std::bad_alloc&) {
            LERROR("Failed to resample edge weight volume: bad allocation");
            delete edgeWeights;
            return;
        }
    }

    outportEdgeWeights_.setData(ewHandle);
}

void RandomWalker::lodMinLevelChanged() {
    lodMaxLevel_.set(std::max(lodMaxLevel_.get(), lodMinLevel_.get()));

    if (inportVolume_.hasData()) {
        float scaleFactor = static_cast<float>(1 << lodMinLevel_.get());
        tgtAssert(scaleFactor >= 1.f, "invalid scale factor");
        tgt::ivec3 volDim = inportVolume_.getData()->getRepresentation<VolumeRAM>()->getDimensions();
        tgt::ivec3 maxDim = tgt::iround(tgt::vec3(volDim) / scaleFactor);
        lodMaxResolution_.set(maxDim);
    }
}

void RandomWalker::lodMaxLevelChanged() {
    lodMinLevel_.set(std::min(lodMinLevel_.get(), lodMaxLevel_.get()));

    if (inportVolume_.hasData()) {
        float scaleFactor = static_cast<float>(1 << lodMaxLevel_.get());
        tgtAssert(scaleFactor >= 1.f, "invalid scale factor");
        tgt::ivec3 volDim = inportVolume_.getData()->getRepresentation<VolumeRAM>()->getDimensions();
        tgt::ivec3 minDim = tgt::iround(tgt::vec3(volDim) / scaleFactor);
        lodMinResolution_.set(minDim);
    }
}

void RandomWalker::updateGuiState() {
    bool clipping = enableClipping_.get();
    clipRegion_.setVisibleFlag(clipping);

    bool useTransFunc = enableTransFunc_.get();
    edgeWeightTransFunc_.setVisibleFlag(useTransFunc);
    edgeWeightBalance_.setVisibleFlag(useTransFunc);

    bool lodEnabled = enableLevelOfDetail_.get();
    lodMaxLevel_.setVisibleFlag(lodEnabled);
    lodMinLevel_.setVisibleFlag(lodEnabled);
    lodForegroundSeedThresh_.setVisibleFlag(lodEnabled);
    lodBackgroundSeedThresh_.setVisibleFlag(lodEnabled);
    lodSeedErosionKernelSize_.setVisibleFlag(lodEnabled);
    lodMinResolution_.setVisibleFlag(lodEnabled);
    lodMaxResolution_.setVisibleFlag(lodEnabled);

    resampleOutputVolumes_.setVisibleFlag(lodEnabled);
}

}   // namespace
