/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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
    : VolumeProcessor(),
    inportVolume_(Port::INPORT, "volume.input"),
    inportForegroundSeeds_(Port::INPORT, "geometry.seedsForeground", "geometry.seedsForeground", true),
    inportBackgroundSeeds_(Port::INPORT, "geometry.seedsBackground", "geometry.seedsBackground", true),
    inportForegroundSeedsVolume_(Port::INPORT, "volume.seedsForeground"),
    inportBackgroundSeedsVolume_(Port::INPORT, "volume.seedsBackground"),
    outportSegmentation_(Port::OUTPORT, "volume.segmentation", "volume.segmentation", false),
    outportProbabilities_(Port::OUTPORT, "volume.probabilities", "volume.probabilities", false),
    outportEdgeWeights_(Port::OUTPORT, "volume.edgeweights", "volume.edgeweights", false),
    computeButton_("computeButton", "Compute"),
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
    useCaching_("useCaching", "Use Cache", false, VALID),
    clearCache_("clearCache", "Clear Cache", VALID),
    runOnInit_("runOnInit", "Run On Initialization", false),
    cache_(this),
    recomputeRandomWalker_(false),
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

    computeButton_.onClick(MemberFunctionCallback<RandomWalker>(this, &RandomWalker::computeButtonClicked));
    addProperty(computeButton_);
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
    addProperty(useCaching_);
    useCaching_.setGroupID("output");
    addProperty(clearCache_);
    clearCache_.setGroupID("output");
    setPropertyGroupGuiName("output", "Output");
    clearCache_.onClick(MemberFunctionCallback<RandomWalker>(this, &RandomWalker::clearCache));

    addProperty(runOnInit_);
}

RandomWalker::~RandomWalker() {
}

Processor* RandomWalker::create() const {
    return new RandomWalker();
}

void RandomWalker::initialize() {
    VolumeProcessor::initialize();

#ifdef VRN_MODULE_OPENCL
    voreenBlasCL_.initialize();
#endif

    cache_.addAllInports();
    cache_.addAllOutports();
    cache_.addAllProperties();
    cache_.initialize();

    updateGuiState();

    if (runOnInit_.get())
        computeButton_.clicked();
}

void RandomWalker::deinitialize() {
    // clear lod volumes
    for (size_t i=0; i<lodVolumes_.size(); i++)
        delete lodVolumes_.at(i);
    lodVolumes_.clear();

    VolumeProcessor::deinitialize();
}

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

void RandomWalker::beforeProcess() {
    VolumeProcessor::beforeProcess();

    if (useCaching_.get() && recomputeRandomWalker_) {
        if (cache_.restore()) {
            recomputeRandomWalker_ = false;
            setValid();
            currentInputVolume_ = inportVolume_.getData()->getRepresentation<VolumeRAM>();
            LINFO("Restored from cache");
        }
    }

    // assign volume to transfer function
    edgeWeightTransFunc_.setVolume(inportVolume_.getData());
    LGL_ERROR;
}

void RandomWalker::process() {

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
    currentInputVolume_ = inportVolume_.getData()->getRepresentation<VolumeRAM>();

    if (recomputeRandomWalker_) {
        auto start = clock::now();
        RandomWalkerSolver* solver = computeRandomWalkerSolution();
        if (solver && solver->getSystemState() == RandomWalkerSolver::Solved) {
            auto finish = clock::now();
            LINFO("Total runtime: " << std::chrono::duration<float>(finish - start).count() << " sec");

            // put out results
            if (outportSegmentation_.isConnected())
                putOutSegmentation(solver);
            else
                outportSegmentation_.setData(0);

            if (outportProbabilities_.isConnected())
                putOutProbabilities(solver);
            else
                outportProbabilities_.setData(0);

            if (outportEdgeWeights_.isConnected())
                putOutEdgeWeights(solver);
            else
                outportEdgeWeights_.setData(0);
        }
        else {
            LERROR("Failed to compute Random Walker solution");
        }
        delete solver;
        recomputeRandomWalker_ = false;

        if (useCaching_.get())
            cache_.store();
    }

}

void RandomWalker::invalidate(int inv) {

    if (!recomputeRandomWalker_ && !inportVolume_.hasChanged())
        return;

    VolumeProcessor::invalidate(inv);
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

RandomWalkerSolver* RandomWalker::computeRandomWalkerSolution() {

    if (!inportVolume_.hasData() || !inportVolume_.getData()->getRepresentation<VolumeRAM>()) {
        LWARNING("No volume");
        return 0;
    }
    const VolumeBase* inputHandle = inportVolume_.getData();
    const VolumeRAM* inputVolume = inputHandle->getRepresentation<VolumeRAM>();

    const int startLevel = enableLevelOfDetail_.get() ? lodMaxLevel_.get() : 0;
    const int endLevel = enableLevelOfDetail_.get() ? lodMinLevel_.get() : 0;
    tgtAssert(startLevel-endLevel >= 0, "invalid level range");

    // create lod volumes, if not present
    bool computeLODs = false;
    for (int level=std::max(endLevel, 1); level<=startLevel && !computeLODs; level++)
        computeLODs |= ((level-1) >= (int)lodVolumes_.size()) || (lodVolumes_.at(level-1) == 0);
    if (computeLODs) {
        LINFO("Computing level of detail volumes...");
        auto start = clock::now();
        for (int level=std::max(endLevel, 1); level<=startLevel; level++) {
            while ((level-1) >= (int)lodVolumes_.size())
                lodVolumes_.push_back(0);
            if (lodVolumes_.at(level-1) == 0) {
                float scaleFactor = static_cast<float>(1 << level);
                tgtAssert(scaleFactor >= 1.f, "invalid scale factor");
                tgt::ivec3 levelDim = tgt::iround(tgt::vec3(inputVolume->getDimensions()) / scaleFactor);
                try {
                    Volume* levelVolume = VolumeOperatorResample::APPLY_OP(inputHandle, levelDim, VolumeRAM::LINEAR);
                    lodVolumes_.at(level-1) = levelVolume;
                }
                catch (std::bad_alloc&) {
                    LERROR("Failed create level-" << level << " volume (dim=" << levelDim << ") : bad allocation");
                    return 0;
                }
            }
        }
        auto end = clock::now();
        LINFO("...finished (" << (end-start).count() << " sec)");
    }

    // work resources
    RandomWalkerSolver* solver = 0;
    const VolumeBase* workVolume = 0;
    VolumeRAM_UInt8* foregroundSeedVol = 0;
    VolumeRAM_UInt8* backgroundSeedVol = 0;

    // input seed lists
    PointSegmentListGeometryVec3 foregroundSeedsPort;
    PointSegmentListGeometryVec3 backgroundSeedsPort;
    getSeedListsFromPorts(foregroundSeedsPort, backgroundSeedsPort);

    // input seed volumes
    try {
        if (inportForegroundSeedsVolume_.hasData()) {
            const VolumeBase* inputSeedVol = inportForegroundSeedsVolume_.getData();
            foregroundSeedVol = dynamic_cast<const VolumeRAM_UInt8*>(inputSeedVol->getRepresentation<VolumeRAM>())->clone();

            if (!foregroundSeedVol) {
                VolumeOperatorConvert converter;
                Volume* h = converter.apply<uint8_t>(inportForegroundSeedsVolume_.getData());
                foregroundSeedVol = dynamic_cast<VolumeRAM_UInt8*>(h->getWritableRepresentation<VolumeRAM>());
                h->releaseAllRepresentations();
                delete h;
            }
        }
        if (inportBackgroundSeedsVolume_.hasData()) {
            const VolumeBase* inputSeedVol = inportBackgroundSeedsVolume_.getData();
            backgroundSeedVol = dynamic_cast<const VolumeRAM_UInt8*>(inputSeedVol->getRepresentation<VolumeRAM>())->clone();

            if (!backgroundSeedVol) {
                VolumeOperatorConvert converter;
                Volume* h = converter.apply<uint8_t>(inportBackgroundSeedsVolume_.getData());
                backgroundSeedVol = dynamic_cast<VolumeRAM_UInt8*>(h->getWritableRepresentation<VolumeRAM>());
                h->releaseAllRepresentations();
                delete h;
            }
        }
    }
    catch (std::bad_alloc&) {
        LERROR("Failed to convert input seed volumes: bad allocation");
        goto finalize;
    }

    //
    // Multi Scale Loop
    //
    for (int level = startLevel; level >= endLevel; level--) {

        LoopRecord loopRecord;
        loopRecord.iteration = startLevel-level+1;
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
            tgtAssert((level-1) < (int)lodVolumes_.size() && lodVolumes_.at(level-1), "lod volume missing");
            workVolume = lodVolumes_.at(level-1);
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
        if (enableClipping_.get()) {
            clipLLF = tgt::iround(tgt::vec3(clipRegion_.get().getLLF()) / scaleFactor);
            clipURB = tgt::iround(tgt::vec3(clipRegion_.get().getURB()) / scaleFactor);
        }

        RandomWalkerSeeds* seeds = new RandomWalkerTwoLabelSeeds(foregroundSeeds, backgroundSeeds,
            foregroundSeedVol, backgroundSeedVol, clipLLF, clipURB);

        /*
         * 2. Edge weight calculator (independent from scale level)
         */
        RandomWalkerWeights* weights = getEdgeWeightsFromProperties();

        /*
         * 3. Set up Random Walker system.
         */
        //LINFO("Constructing Random Walker equation system...");
        delete solver;
        solver = new RandomWalkerSolver(workVolume, seeds, weights);
        try {
            auto start = clock::now();
            solver->setupEquationSystem();
            auto finish = clock::now();
            //LINFO("...finished: " << std::chrono::duration<float>(finish-start).count() << " sec");
            loopRecord.timeSetup = finish - start;
            loopRecord.numSeeds = seeds->getNumSeeds();
            if (RandomWalkerTwoLabelSeeds* twoLabelSeeds = dynamic_cast<RandomWalkerTwoLabelSeeds*>(seeds)) {
                loopRecord.numForegroundSeeds = twoLabelSeeds->getNumForegroundSeeds();
                loopRecord.numBackgroundSeeds = twoLabelSeeds->getNumBackgroundSeeds();
            }
        }
        catch (tgt::Exception& e) {
            LERROR("Failed to setup Random Walker equation system: " << e.what());
            goto finalize;
        }

        /*
         * 4. Compute Random Walker solution.
         */
        // select BLAS implementation and preconditioner
        const VoreenBlas* voreenBlas = getVoreenBlasFromProperties();
        VoreenBlas::ConjGradPreconditioner precond = VoreenBlas::NoPreconditioner;
        if (preconditioner_.isSelected("jacobi"))
            precond = VoreenBlas::Jacobi;

        // solve
        float errorThresh = 1.f / pow(10.f, static_cast<float>(errorThreshold_.get()));
        try {
            float* initialization = nullptr;
            std::vector<float> initializationStorage;
            initializationStorage.reserve(solver->getSystemSize());
            if(usePrevProbAsInitialization_.get() && prevProbabilities_.size() > 0) {
                tgtAssert(prevProbabilities_.size() == solver->getNumVoxels(), "Old and new probalities size missmatch");
                // This is highly depends on the implementation of RandomWalkerSolver!
                for(size_t i = 0; i < prevProbabilities_.size(); ++i) {
                    if(!solver->isSeedPoint(i)) {
                        initializationStorage.push_back(prevProbabilities_[i]);
                    }
                }
                LINFO("Using existing probalities as initialization");
                initialization = initializationStorage.data();
            }

            auto start = clock::now();
            int iterations = solver->solve(voreenBlas, initialization, precond, errorThresh, maxIterations_.get());
            auto finish = clock::now();
            loopRecord.timeSolving = finish - start;
            loopRecord.numIterations = iterations;

            tgtAssert(solver->getSystemState() == RandomWalkerSolver::Solved, "System not solved");
            prevProbabilities_ = std::vector<float>();
            prevProbabilities_.reserve(solver->getNumVoxels());
            for(size_t i = 0; i < solver->getNumVoxels(); ++i) {
                prevProbabilities_.push_back(solver->getProbabilityValue(i));
            }
        }
        catch (VoreenException& e) {
            LERROR("Failed to compute Random Walker solution: " << e.what());
            goto finalize;
        }

        loopRecord.probabilityRange = solver->getProbabilityRange();

        /*
         * 5. Derive seed volumes for next iteration
         */
        bool lastIteration = (level == lodMinLevel_.get());
        if (!lastIteration) {

            auto start = clock::now();

            // generate probability map from current solution
            VolumeRAM_UInt16* probabilityVolume = 0;
            try {
                probabilityVolume = solver->generateProbabilityVolume<VolumeRAM_UInt16>();
            }
            catch (VoreenException& e) {
                LERROR("computeRandomWalkerSolution() Failed to generate probability volume: " << e.what());
                goto finalize;
            }

            // previous seed volumes not needed anymore
            delete foregroundSeedVol;
            delete backgroundSeedVol;
            foregroundSeedVol = 0;
            backgroundSeedVol = 0;

            // allocate seed volumes for next iteration
            try {
                foregroundSeedVol = new VolumeRAM_UInt8(workDim);
                backgroundSeedVol = new VolumeRAM_UInt8(workDim);
            }
            catch (std::bad_alloc&) {
                LERROR("computeRandomWalkerSolution() Failed to create seed volumes: bad allocation");
                delete probabilityVolume;
                goto finalize;
            }

            // threshold probability map to derive seeds
            float maxProbValue = probabilityVolume->elementRange().y;
            int foregroundThresh = tgt::iround(lodForegroundSeedThresh_.get()*maxProbValue);
            int backgroundThresh = tgt::iround(lodBackgroundSeedThresh_.get()*maxProbValue);
            foregroundSeedVol->clear();
            backgroundSeedVol->clear();

            for (size_t i=0; i<probabilityVolume->getNumVoxels(); i++) {
                int probValue = probabilityVolume->voxel(i);
                if (probValue <= backgroundThresh)
                    backgroundSeedVol->voxel(i) = 255;
                else if (probValue >= foregroundThresh)
                    foregroundSeedVol->voxel(i) = 255;
            }
            delete probabilityVolume;
            probabilityVolume = 0;

            // erode obtained fore- and background seed volumes
            Volume* erh = VolumeOperatorCubeErosion::APPLY_OP(new Volume(foregroundSeedVol, vec3(1.0f), vec3(0.0f)), lodSeedErosionKernelSize_.getValue()); //FIXME: small memory leak
            VolumeRAM_UInt8* er = static_cast<VolumeRAM_UInt8*>(erh->getWritableRepresentation<VolumeRAM>());
            delete foregroundSeedVol;
            foregroundSeedVol = er;
            erh->releaseAllRepresentations();
            delete erh;
            erh = VolumeOperatorCubeErosion::APPLY_OP(new Volume(backgroundSeedVol, vec3(1.0f), vec3(0.0f)), lodSeedErosionKernelSize_.getValue()); //FIXME: small memory leak
            er = static_cast<VolumeRAM_UInt8*>(erh->getWritableRepresentation<VolumeRAM>());
            erh->releaseAllRepresentations();
            delete erh;
            erh = 0;
            delete backgroundSeedVol;
            backgroundSeedVol = er;
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

    } // end loop

    tgtAssert(solver, "no random walker solver");

finalize:
    // we can delete these because they are cloned/converted
    delete foregroundSeedVol;
    delete backgroundSeedVol;

    // return solver carrying the computed random walker solution
    return solver;
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

void RandomWalker::getSeedListsFromPorts(PointSegmentListGeometry<tgt::vec3>& foregroundSeeds,
    PointSegmentListGeometry<tgt::vec3>& backgroundSeeds) const {

    std::vector<const Geometry*> foregroundGeom = inportForegroundSeeds_.getAllData();
    std::vector<const Geometry*> backgroundGeom = inportBackgroundSeeds_.getAllData();

    for (size_t i=0; i<foregroundGeom.size(); i++) {
        const PointSegmentListGeometry<tgt::vec3>* seedList = dynamic_cast<const PointSegmentListGeometry<tgt::vec3>* >(foregroundGeom.at(i));
        if (!seedList)
            LWARNING("Invalid geometry. PointSegmentListGeometry<vec3> expected.");
        else {
            auto transformMat = seedList->getTransformationMatrix();
            for (int j=0; j<seedList->getNumSegments(); j++) {
                std::vector<tgt::vec3> points;
                for(auto& vox : seedList->getSegment(j)) {
                    points.push_back(transformMat.transform(vox));
                }
                foregroundSeeds.addSegment(points);
            }
        }
    }

    for (size_t i=0; i<backgroundGeom.size(); i++) {
        const PointSegmentListGeometry<tgt::vec3>* seedList = dynamic_cast<const PointSegmentListGeometry<tgt::vec3>* >(backgroundGeom.at(i));
        if (!seedList)
            LWARNING("Invalid geometry. PointSegmentListGeometry<vec3> expected.");
        else {
            auto transformMat = seedList->getTransformationMatrix();
            for (int j=0; j<seedList->getNumSegments(); j++) {
                std::vector<tgt::vec3> points;
                for(auto& vox : seedList->getSegment(j)) {
                    points.push_back(transformMat.transform(vox));
                }
                backgroundSeeds.addSegment(points);
            }
        }
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

void RandomWalker::computeButtonClicked() {
    recomputeRandomWalker_ = true;
    invalidate(Processor::INVALID_RESULT);
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

void RandomWalker::clearCache()  {
    cache_.clearCache();
}

}   // namespace
