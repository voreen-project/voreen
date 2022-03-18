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

#include "supervoxelwalker.h"

#include "../solver/randomwalkerweights.h"
#include "../util/preprocessing.h"
#include "../util/seeds.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/octree/octreeutils.h"
#include "voreen/core/utils/hashing.h"
#include "tgt/vector.h"
#include "tgt/memory.h"

#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"

#include <climits>
#include <random>
#include <eigen3/Eigen/Eigen>

#ifdef VRN_MODULE_OPENMP
#include <omp.h>
#endif

namespace voreen {

namespace {

inline size_t volumeCoordsToIndex(int x, int y, int z, const tgt::ivec3& dim) {
    return z*dim.y*dim.x + y*dim.x + x;
}

inline size_t volumeCoordsToIndex(const tgt::ivec3& coords, const tgt::ivec3& dim) {
    return coords.z*dim.y*dim.x + coords.y*dim.x + coords.x;
}

}

SuperVoxelWalkerOutput::SuperVoxelWalkerOutput(
        std::unique_ptr<VolumeBase>&& result,
        std::unique_ptr<VolumeBase>&& resultSuperVoxels,
        boost::optional<SuperVoxelWalkerPreprocessingResult>&& preprocessingResult,
        std::vector<float> previousSolution,
        std::chrono::duration<float> duration
)   : result_(std::move(result))
    , resultSuperVoxels_(std::move(resultSuperVoxels))
    , preprocessingResult_(std::move(preprocessingResult))
    , previousSolution_(previousSolution)
    , duration_(duration)
{ }

SuperVoxelWalkerOutput::SuperVoxelWalkerOutput(SuperVoxelWalkerOutput&& other)
    : result_(std::move(other.result_))
    , resultSuperVoxels_(std::move(other.resultSuperVoxels_))
    , preprocessingResult_(std::move(other.preprocessingResult_))
    , duration_(other.duration_)
    , previousSolution_(other.previousSolution_)
    , movedOut_(other.movedOut_)
{
}

SuperVoxelWalkerOutput::~SuperVoxelWalkerOutput() {
}

const std::string SuperVoxelWalker::loggerCat_("voreen.RandomWalker.SuperVoxelWalker");

SuperVoxelWalker::SuperVoxelWalker()
    : AsyncComputeProcessor<SuperVoxelWalkerInput, SuperVoxelWalkerOutput>()
    , inportVolume_(Port::INPORT, "volume.input")
    , inportForegroundSeeds_(Port::INPORT, "geometry.seedsForeground", "geometry.seedsForeground", true)
    , inportBackgroundSeeds_(Port::INPORT, "geometry.seedsBackground", "geometry.seedsBackground", true)
    , outportProbabilities_(Port::OUTPORT, "volume.probabilities", "volume.probabilities", false)
    , outportSuperVoxels_(Port::OUTPORT, "volume.supervoxels", "volume.supervoxels", false)
    , minEdgeWeight_("minEdgeWeight", "Min Edge Weight: 10^(-t)", 5, 0, 10)
    , beta_("beta", "Beta: 2^v", 0, 0, 20, Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_DEBUG)
    , preconditioner_("preconditioner", "Preconditioner")
    , errorThreshold_("errorThreshold", "Error Threshold: 10^(-t)", 2, 0, 10)
    , maxIterations_("conjGradIterations", "Max Iterations", 1000, 1, 5000)
    , conjGradImplementation_("conjGradImplementation", "Implementation")
    , maxSuperVoxelVolume_("maxSuperVoxelVolume", "Max Super Voxel Volume", 10000, 1, INT_MAX)
    , maxSuperVoxelIntensityDifference_("maxSuperVoxelIntensityDifference", "Max Super Voxel Intensity Difference", 1.0, 0.0, 10000000000.0)
    , clearResult_("clearResult", "Clear Result", Processor::INVALID_RESULT, Property::LOD_DEBUG)
    , preprocessingResult_(boost::none)
    , previousSolution_(boost::none)
{
    // ports
    addPort(inportVolume_);
    addPort(inportForegroundSeeds_);
    addPort(inportBackgroundSeeds_);
    addPort(outportProbabilities_);
    addPort(outportSuperVoxels_);
    ON_CHANGE_LAMBDA(inportVolume_, [this] () {
            clearPreviousResults();
            if(inportVolume_.hasData()) {
                auto vol = inportVolume_.getData();
                auto rwm = vol->getRealWorldMapping();
                float maxRange = std::abs(rwm.normalizedToRealWorld(1.0f) - rwm.normalizedToRealWorld(0.0f));
                maxSuperVoxelIntensityDifference_.setMaxValue(maxRange);

                size_t maxVol = vol->getNumVoxels()/(5*5*5);
                maxSuperVoxelVolume_.setMaxValue(std::min(static_cast<size_t>(INT_MAX), maxVol));
            }
            });

    // random walker properties
    addProperty(minEdgeWeight_);
        minEdgeWeight_.setGroupID("rwparam");
        minEdgeWeight_.setTracking(false);
    addProperty(beta_);
        beta_.setGroupID("rwparam");
        beta_.setTracking(false);
    addProperty(maxSuperVoxelVolume_);
        maxSuperVoxelVolume_.setGroupID("rwparam");
        maxSuperVoxelVolume_.setTracking(false);
        ON_CHANGE(maxSuperVoxelVolume_, SuperVoxelWalker, clearPreviousResults);
    addProperty(maxSuperVoxelIntensityDifference_);
        maxSuperVoxelIntensityDifference_.setGroupID("rwparam");
        maxSuperVoxelIntensityDifference_.adaptDecimalsToRange(5);
        maxSuperVoxelIntensityDifference_.setTracking(false);
        ON_CHANGE(maxSuperVoxelIntensityDifference_, SuperVoxelWalker, clearPreviousResults);
    setPropertyGroupGuiName("rwparam", "Random Walker Parametrization");

    // conjugate gradient solver
    addProperty(preconditioner_);
        preconditioner_.addOption("none", "None");
        preconditioner_.addOption("jacobi", "Jacobi");
        preconditioner_.select("jacobi");
        preconditioner_.setGroupID("conjGrad");
    addProperty(errorThreshold_);
        errorThreshold_.setGroupID("conjGrad");
        errorThreshold_.setTracking(false);
    addProperty(maxIterations_);
        maxIterations_.setGroupID("conjGrad");
        maxIterations_.setTracking(false);
    addProperty(conjGradImplementation_);
        conjGradImplementation_.addOption("blasCPU", "CPU");
#ifdef VRN_MODULE_OPENMP
        conjGradImplementation_.addOption("blasMP", "OpenMP");
        conjGradImplementation_.select("blasMP");
#endif
#ifdef VRN_MODULE_OPENCL
        conjGradImplementation_.addOption("blasCL", "OpenCL");
        conjGradImplementation_.select("blasCL");
#endif
        conjGradImplementation_.setGroupID("conjGrad");
    setPropertyGroupGuiName("conjGrad", "Conjugate Gradient Solver");

    // conjugate gradient solver
    addProperty(clearResult_);
        ON_CHANGE_LAMBDA(clearResult_, [this] () {
                interruptComputation();
                clearPreviousResults();
                });
        clearResult_.setGroupID("resultcache");
    setPropertyGroupGuiName("resultcache", "Result Cache");
}

SuperVoxelWalker::~SuperVoxelWalker() {
}

Processor* SuperVoxelWalker::create() const {
    return new SuperVoxelWalker();
}

void SuperVoxelWalker::initialize() {
    AsyncComputeProcessor::initialize();

#ifdef VRN_MODULE_OPENCL
    voreenBlasCL_.initialize();
#endif
}

void SuperVoxelWalker::deinitialize() {
    interruptComputation();
    clearPreviousResults();

    AsyncComputeProcessor::deinitialize();
}

bool SuperVoxelWalker::isReady() const {
    bool ready = false;
    ready |= outportProbabilities_.isConnected();
    ready &= inportVolume_.isReady();
    ready &= inportForegroundSeeds_.isReady();
    ready &= inportBackgroundSeeds_.isReady();
    return ready;
}

SuperVoxelWalker::ComputeInput SuperVoxelWalker::prepareComputeInput() {
    tgtAssert(inportVolume_.hasData(), "no input volume");

    // clear previous results and update property ranges, if input volume has changed
    if (inportVolume_.hasChanged()) {
        outportProbabilities_.setData(0);
    }
    auto vol = inportVolume_.getThreadSafeData();

    if (vol->getNumChannels() != 1) {
        throw InvalidInputException("Only single channel volumes are supported.", InvalidInputException::S_ERROR);
    }

    const VolumeRAM* volram = vol->getRepresentation<VolumeRAM>();
    // clear previous results and update property ranges, if input volume has changed
    if (!volram) {
        throw InvalidInputException("Failed to create VolumeRAM representation", InvalidInputException::S_ERROR);
    }

    // select BLAS implementation and preconditioner
    const VoreenBlas* voreenBlas = getVoreenBlasFromProperties();
    VoreenBlas::ConjGradPreconditioner precond = VoreenBlas::NoPreconditioner;
    if (preconditioner_.isSelected("jacobi"))
        precond = VoreenBlas::Jacobi;

    float errorThresh = 1.f / pow(10.f, static_cast<float>(errorThreshold_.get()));
    int maxIterations = maxIterations_.get();

    float beta = std::pow(2.0f, beta_.get());
    float minWeight = std::pow(10.0f, -minEdgeWeight_.get());

    return ComputeInput {
        preprocessingResult_ ? &*preprocessingResult_ : nullptr,
        previousSolution_ ? &*previousSolution_ : nullptr,
        *vol,
        *volram,
        inportForegroundSeeds_.getThreadSafeAllData(),
        inportBackgroundSeeds_.getThreadSafeAllData(),
        minWeight,
        beta,
        voreenBlas,
        precond,
        errorThresh,
        maxIterations,
        maxSuperVoxelVolume_.get(),
        maxSuperVoxelIntensityDifference_.get(),
    };
}

const SuperVoxelID UNLABELED = 0;
const SuperVoxelID LABEL_START = 1;

static void growRegion(const VolumeAtomic<float>& vol, VolumeAtomic<SuperVoxelID>& labels, std::vector<float>& regionMeans, tgt::ivec3 start, SuperVoxelID label, size_t maxVol, float maxDiff, RealWorldMapping rwm) {
    std::vector<tgt::ivec3> voxelQueue;
    voxelQueue.push_back(start);
    std::vector<tgt::ivec3> neighborQueue;

    size_t numVoxels = 0;
    float intensitySum = 0.0f;
    tgt::ivec3 dim = vol.getDimensions();
    tgtAssert(tgt::ivec3(labels.getDimensions()) == dim, "Dimension mismatch");

    do {
        for(auto& p : voxelQueue) {
            if(labels.voxel(p) != UNLABELED) {
                continue;
            }
            float intensity = rwm.normalizedToRealWorld(vol.voxel(p));
            intensitySum += intensity;

            labels.voxel(p) = label;
            numVoxels += 1;


            for(int oz : {-1, 0, 1}) {
                for(int oy : {-1, 0, 1}) {
                    for(int ox : {-1, 0, 1}) {
                        tgt::ivec3 neighbor = p;
                        neighbor.x += ox;
                        neighbor.y += oy;
                        neighbor.z += oz;

                        if(tgt::hand(tgt::greaterThanEqual(neighbor, tgt::ivec3::zero)) && tgt::hand(tgt::lessThan(neighbor, dim))
                                && labels.voxel(neighbor) == UNLABELED
                                && tgt::abs(rwm.normalizedToRealWorld(vol.voxel(neighbor)) - intensity) < maxDiff) {

                            neighborQueue.push_back(neighbor);
                        }
                    }
                }
            }

            if(numVoxels == maxVol) {
                goto done;
            }
        }
        std::swap(voxelQueue, neighborQueue);
        neighborQueue.clear();
    } while(!voxelQueue.empty());
done:

    regionMeans.push_back(intensitySum/numVoxels);
}

static SuperVoxelWalkerPreprocessingResult preprocess(const SuperVoxelWalkerInput& input, ProgressReporter& progressReporter) {
    SubtaskProgressReporterCollection<2> globalProgressSteps(progressReporter, {0.5, 1.0});
    VolumeAtomic<float> vol = toVolumeAtomicFloat(input.volram_);
    tgt::svec3 dim = vol.getDimensions();
    VolumeAtomic<SuperVoxelID> labels(vol.getDimensions());
    labels.fill(UNLABELED);
    std::vector<float> regionMeans;
    regionMeans.push_back(0.0); //Dummy for UNLABELED
    auto rwm = input.volume_.getRealWorldMapping();

    size_t numVoxels = vol.getNumVoxels();

    std::vector<size_t> indices;
    indices.reserve(numVoxels);
    for(size_t i=0; i<numVoxels; ++i) {
        indices.push_back(i);
    }

    //TODO make random seed configurable?
    int randomSeed = 0;
    std::mt19937 g(randomSeed);
    std::shuffle(indices.begin(), indices.end(), g);

    SuperVoxelID label = LABEL_START;
    int i=0;
    for(auto index : indices) {
        if(labels.voxel(index) == UNLABELED) {
            tgt::ivec3 p = linearCoordToCubic(index, dim);
            growRegion(vol, labels, regionMeans, p, label, input.maxSuperVoxelVolume_, input.maxSuperVoxelIntensityDifference_, rwm);
            label += 1;
        }
        globalProgressSteps.get<0>().setProgress(static_cast<float>(i++)/numVoxels);
    }
    indices.clear();

    size_t numSuperVoxels = regionMeans.size() - 1;

    std::vector<std::vector<SuperVoxelID>> edges(numSuperVoxels+1);

    i=0;
    VRN_FOR_EACH_VOXEL(p, tgt::ivec3(0,0,0), tgt::ivec3(dim)) {
        SuperVoxelID cls = labels.voxel(p);
        std::vector<SuperVoxelID>& clsEdges = edges[cls];
        for(int d=0; d<3; ++d) {
            for(int o : {-1, 1}) {
                tgt::ivec3 neighbor = p;
                neighbor[d] += o;

                if(neighbor[d] < 0 || neighbor[d] == dim[d]) {
                    continue;
                }

                SuperVoxelID neighCls = labels.voxel(neighbor);
                if(cls == neighCls) {
                    continue;
                }

                bool found = false;
                for(auto i : clsEdges) {
                    if(i == neighCls) {
                        found = true;
                        break;
                    }
                }
                if(!found) {
                    clsEdges.push_back(neighCls);
                    edges[neighCls].push_back(cls);
                }
            }
        }
        globalProgressSteps.get<1>().setProgress(static_cast<float>(i++)/numVoxels);
    }

    return SuperVoxelWalkerPreprocessingResult {
        std::move(edges),
        std::move(regionMeans),
        std::move(labels),
    };
}
static std::unique_ptr<VolumeBase> createSuperVoxelVolume(const SuperVoxelWalkerInput& input, const SuperVoxelWalkerPreprocessingResult& preprocessingResult) {
    tgt::svec3 dim = input.volume_.getDimensions();
    auto meanvol = tgt::make_unique<VolumeAtomic<float>>(dim);
    VRN_FOR_EACH_VOXEL(p, tgt::ivec3(0,0,0), tgt::ivec3(dim)) {
        meanvol->voxel(p) = preprocessingResult.regionMeans_[preprocessingResult.labels_.voxel(p)];
    }
    auto output = tgt::make_unique<Volume>(meanvol.release(), input.volume_.getSpacing(), input.volume_.getOffset(), input.volume_.getPhysicalToWorldMatrix());

    //auto output = tgt::make_unique<Volume>(preprocessingResult.labels_.clone(), input.volume_.getSpacing(), input.volume_.getOffset(), input.volume_.getPhysicalToWorldMatrix());
    //output->setRealWorldMapping(RealWorldMapping::createDenormalizingMapping<SuperVoxelID>());
    return output;
}

static Eigen::Matrix<float, Eigen::Dynamic, 1> solveSystem(const SuperVoxelWalkerInput& input, const Eigen::SparseMatrix<float, Eigen::ColMajor>& mat, float* vec, float* init) {
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, 1>> btms(vec, mat.rows());
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, 1>> init_eigen(init, mat.rows());

    Eigen::ConjugateGradient<Eigen::SparseMatrix<float>, Eigen::Lower|Eigen::Upper> solver;
    solver.setTolerance(input.errorThreshold_);
    solver.setMaxIterations(input.maxIterations_);
    solver.compute(mat);

    Eigen::Matrix<float, Eigen::Dynamic, 1> solution = solver.solveWithGuess(btms, init_eigen);
    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error()      << std::endl;


    //Eigen::SimplicialLDLT<Eigen::SparseMatrix<float, Eigen::RowMajor>> solver;
    //solver.analyzePattern(mat);
    //solver.factorize(mat);

    //Eigen::Matrix<float, Eigen::Dynamic, 1> solution = solver.solve(btms);


    return solution;
}

static std::pair<std::unique_ptr<VolumeAtomic<float>>, std::vector<float>> solve(const SuperVoxelWalkerInput& input, const SuperVoxelWalkerPreprocessingResult& preprocessingResult, ProgressReporter& progressReporter) {
    tgt::svec3 dim = input.volume_.getDimensions();

    PointSegmentListGeometryVec3 foregroundSeeds;
    PointSegmentListGeometryVec3 backgroundSeeds;
    getSeedListsFromPorts(input.foregroundGeomSeeds_, foregroundSeeds);
    getSeedListsFromPorts(input.backgroundGeomSeeds_, backgroundSeeds);

    const int UNSEEDED = -1;
    const int SEED_FG = 1;
    const int SEED_BG = 0;

    size_t numSuperVoxels = preprocessingResult.edges_.size() - 1;

    std::vector<int> superVoxelSeedClasses(preprocessingResult.edges_.size(), UNSEEDED);
    {
        RandomWalkerTwoLabelSeeds seeds(dim,foregroundSeeds,backgroundSeeds);
        seeds.initialize();

        if(seeds.getNumSeeds() == 0) {
            LERRORC(SuperVoxelWalker::loggerCat_, "No seeds");

            auto outvol = tgt::make_unique<VolumeAtomic<float>>(dim);
            outvol->fill(0.5);
            std::vector<float> previousSolution(preprocessingResult.edges_.size(), 0.5);
            return std::make_pair(std::move(outvol), previousSolution);
        }

        VRN_FOR_EACH_VOXEL(p, tgt::ivec3(0,0,0), tgt::ivec3(dim)) {
            SuperVoxelID cls = preprocessingResult.labels_.voxel(p);
            tgtAssert(cls != UNLABELED, "All voxels should be labeled");

            if(seeds.isSeedPoint(p)) {
                int seedClass = seeds.getSeedValue(p) == 1.0 ? SEED_FG : SEED_BG;
                int& existingSeedClass = superVoxelSeedClasses[cls];

                if(existingSeedClass != seedClass && existingSeedClass != UNSEEDED) {
                    LWARNINGC(SuperVoxelWalker::loggerCat_, "Super voxel seed class conflict");
                }
                existingSeedClass = seedClass;
            }
        }
    }

    std::vector<int> superVoxelToMatTable;
    superVoxelToMatTable.reserve(superVoxelSeedClasses.size());

    size_t seededIndex = 0;
    size_t unseededIndex = 0;
    for(int i = 1; i < superVoxelSeedClasses.size(); ++i) {
        if(superVoxelSeedClasses[i] != UNSEEDED) {
            superVoxelToMatTable[i] = seededIndex;
            ++seededIndex;
        } else {
            superVoxelToMatTable[i] = unseededIndex;
            ++unseededIndex;
        }
    }
    size_t unseeded = unseededIndex;
    size_t seeded = seededIndex;
    tgtAssert(seeded + unseeded == numSuperVoxels, "Seeded/Unseeded/total count mismatch");


    float squareDiffStd;
    {
        // Find normalization factor:
        std::vector<float> squareDifferences;
        float sum = 0.0f;
        for(size_t i = 1; i < superVoxelSeedClasses.size(); ++i) {
            const std::vector<SuperVoxelID>& edges = preprocessingResult.edges_[i];

            float weightSum = 0;
            for(auto j : edges) {
                float diff = preprocessingResult.regionMeans_[i] - preprocessingResult.regionMeans_[j];
                float diffSquared = diff*diff;
                squareDifferences.push_back(diffSquared);
                sum += diffSquared;
            }
        }
        float mean = sum/squareDifferences.size();
        float varsum = 0.0f;
        for(float val : squareDifferences) {
            float diff = mean-val;
            varsum += diff*diff;
        }
        float variance = varsum/(squareDifferences.size() -1);
        squareDiffStd = std::sqrt(variance);
    }

    std::vector<Eigen::Triplet<float>> triplets_lu;
    std::vector<float> diagonal(unseeded, 0.0f);


    std::vector<float> vec(unseeded, 0.0f);
    for(size_t i = 1; i < superVoxelSeedClasses.size(); ++i) {
        const std::vector<SuperVoxelID>& edges = preprocessingResult.edges_[i];
        int iSeed = superVoxelSeedClasses[i];
        size_t iRow = superVoxelToMatTable[i];

        float weightSum = 0;
        for(auto j : edges) {
            if(j > i) {
                continue;
            }
            int jSeed = superVoxelSeedClasses[j];
            float diff = preprocessingResult.regionMeans_[i] - preprocessingResult.regionMeans_[j];
            float w = std::max(std::exp(-input.beta_*diff*diff / squareDiffStd), input.minWeight_);

            tgtAssert(!std::isnan(w) && std::isfinite(w), "Invalid weight");

            if(jSeed != UNSEEDED) {
                if(iSeed == UNSEEDED) {
                    vec[iRow] += w * jSeed;
                }
            } else {
                size_t jRow = superVoxelToMatTable[j];
                if(iSeed == UNSEEDED) {
                    triplets_lu.emplace_back(iRow, jRow, -w);
                    triplets_lu.emplace_back(jRow, iRow, -w);
                } else {
                    vec[jRow] += w * iSeed;
                }

                // Update weight sum of neighbor with smaller index.
                diagonal[jRow] += w;
            }
            weightSum += w;
        }

        if(iSeed == UNSEEDED) {
            // This is the first time writing to mat at this location, so overwriting is fine.
            diagonal[iRow] += weightSum;
        }
    }

    for(int i=0; i<unseeded; ++i) {
        triplets_lu.emplace_back(i,i,diagonal[i]);
    }
    diagonal.clear();

    Eigen::SparseMatrix<float, Eigen::ColMajor> lu(unseeded,unseeded);

    lu.setFromTriplets(triplets_lu.begin(), triplets_lu.end());
    triplets_lu.clear();

    lu.makeCompressed();

    std::vector<float> initialization(unseeded);
    if(input.previousSolution_) {
        for(size_t i = 1; i < superVoxelSeedClasses.size(); ++i) {
            size_t iRow = superVoxelToMatTable[i];
            int iSeed = superVoxelSeedClasses[i];
            if(iSeed == UNSEEDED) {
                initialization[iRow] = (*input.previousSolution_)[i];
            }
        }
    } else {
        std::fill(initialization.begin(), initialization.end(), 0.5f);
    }

    auto solution = solveSystem(input, lu, vec.data(), initialization.data());

    auto outvol = tgt::make_unique<VolumeAtomic<float>>(dim);

    std::vector<float> futurePrevSolution(preprocessingResult.edges_.size());
    for(size_t i = 1; i < superVoxelSeedClasses.size(); ++i) {
        size_t row = superVoxelToMatTable[i];
        int seed = superVoxelSeedClasses[i];

        float& out = futurePrevSolution[i];
        if(seed == UNSEEDED) {
            out = solution[row];
        } else {
            out = seed;
        }
    }

    VRN_FOR_EACH_VOXEL(p, tgt::ivec3(0,0,0), tgt::ivec3(dim)) {
        SuperVoxelID cls = preprocessingResult.labels_.voxel(p);
        tgtAssert(cls != UNLABELED, "All voxels should be labeled");

        outvol->voxel(p) = futurePrevSolution[cls];
    }

    return std::make_pair(std::move(outvol), futurePrevSolution);
}

static std::pair<std::unique_ptr<VolumeBase>, std::vector<float>> runRW(const SuperVoxelWalkerInput& input, const SuperVoxelWalkerPreprocessingResult& preprocessingResult, ProgressReporter& progressReporter) {
    auto o = solve(input, preprocessingResult, progressReporter);
    auto output = tgt::make_unique<Volume>(o.first.release(), input.volume_.getSpacing(), input.volume_.getOffset(), input.volume_.getPhysicalToWorldMatrix());
    return std::make_pair(std::move(output), o.second);
}

SuperVoxelWalker::ComputeOutput SuperVoxelWalker::compute(ComputeInput input, ProgressReporter& progressReporter) const {
    progressReporter.setProgress(0.0);

    auto start = clock::now();

    SubtaskProgressReporterCollection<2> globalProgressSteps(progressReporter, {0.1, 0.9});

    const tgt::svec3 volumeDim = input.volume_.getDimensions();
    const size_t numChannels = 1;

    boost::optional<SuperVoxelWalkerPreprocessingResult> ppNew = boost::none;


    const SuperVoxelWalkerPreprocessingResult* preprocessingResult;
    if(input.preprocessingResult_) {
        preprocessingResult = input.preprocessingResult_;
    } else {
        ppNew = preprocess(input, globalProgressSteps.get<0>());
        preprocessingResult = &*ppNew;
    }

    auto conn = preprocessingResult->maxConnectivity_;
    auto numSuperVoxels = preprocessingResult->edges_.size() - 1;
    auto totalNumEntries = 0;
    for(const auto& e : preprocessingResult->edges_) {
        totalNumEntries += e.size();
    }
    LINFO("Total number of supervoxels: " << numSuperVoxels);
    LINFO("Maximum supervoxel connectivity: " << conn);
    LINFO("Estimated memory for matrix: " << formatMemorySize(conn * numSuperVoxels * 4));
    LINFO("Total number of matrix entries: " << totalNumEntries);

    //std::unique_ptr<VolumeBase> output = nullptr;

    auto outputSuperVoxels = createSuperVoxelVolume(input, *preprocessingResult);
    auto output = runRW(input, *preprocessingResult, globalProgressSteps.get<1>());
    auto finish = clock::now();

    progressReporter.setProgress(1.0);
    return ComputeOutput (
        std::move(output.first),
        std::move(outputSuperVoxels),
        std::move(ppNew),
        std::move(output.second),
        finish - start
    );
}

void SuperVoxelWalker::processComputeOutput(ComputeOutput output) {
    LINFO("Total runtime: " << output.duration_.count() << " sec");

    if(output.preprocessingResult_) {
        preprocessingResult_ = std::move(output.preprocessingResult_);
    }

    outportProbabilities_.setData(output.result_.release());
    outportSuperVoxels_.setData(output.resultSuperVoxels_.release());
}
void SuperVoxelWalker::clearPreviousResults() {
    preprocessingResult_ = boost::none;
    previousSolution_ = boost::none;
}

const VoreenBlas* SuperVoxelWalker::getVoreenBlasFromProperties() const {

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

SuperVoxelWalkerPreprocessingResult::~SuperVoxelWalkerPreprocessingResult() {
}
SuperVoxelWalkerPreprocessingResult::SuperVoxelWalkerPreprocessingResult(std::vector<std::vector<SuperVoxelID>>&& edges, std::vector<float>&& regionMeans, VolumeAtomic<SuperVoxelID>&& labels)
    : edges_(std::move(edges))
    , regionMeans_(regionMeans)
    , labels_(std::move(labels))
    , maxConnectivity_(0)
{
    for(auto& e : edges_) {
        maxConnectivity_ = std::max(maxConnectivity_, e.size());
    }
}

}   // namespace
