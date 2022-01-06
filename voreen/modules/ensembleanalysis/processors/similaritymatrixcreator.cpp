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

#include "similaritymatrixcreator.h"

#include "../utils/ensemblehash.h"
#include "../utils/utils.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"

#include <random>

// Determine if we want to use memory mapped files for storing the flags.
// This in general is necessary since matrices will get too big for large ensembles.
#define USE_MEMORY_MAPPED_FILES

#ifdef USE_MEMORY_MAPPED_FILES
#include "voreen/core/datastructures/diskarraystorage.h"
#endif

namespace voreen {

const std::string SimilarityMatrixCreator::loggerCat_("voreen.ensembleanalysis.SimilarityMatrixCreator");

SimilarityMatrixCreator::SimilarityMatrixCreator()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , inport_(Port::INPORT, "inport", "Ensemble Datastructure Input", false)
    , seedMask_(Port::INPORT, "seedmask", "Seed Mask Input (optional)")
    , outport_(Port::OUTPORT, "outport", "Similarity Matrix Output", false)
    , singleChannelSimilarityMeasure_("singleChannelSimilarityMeasure", "Single Field Similarity Measure")
    , isoValue_("isoValue", "Iso-Value", 0.5f, 0.0f, 1.0f)
    , multiChannelSimilarityMeasure_("multiChannelSimilarityMeasure", "Multi Field Similarity Measure")
    , weight_("weight", "Weight", 0.5f, 0.0f, 1.0f)
    , numSeedPoints_("numSeedPoints", "Number of Seed Points", 8192, 1, 131072)
    , seedTime_("seedTime", "Current Random Seed", static_cast<int>(time(0)), std::numeric_limits<int>::min(), std::numeric_limits<int>::max())
    , clearCache_("clearCache", "Clear Cache", Processor::VALID, Property::LOD_ADVANCED)
{
    // Ports
    addPort(inport_);
    addPort(seedMask_);
    addPort(outport_);

    // Calculation
    addProperty(singleChannelSimilarityMeasure_);
    singleChannelSimilarityMeasure_.addOption("isovalue", "Iso-Contours", MEASURE_ISOCONTOURS);
    singleChannelSimilarityMeasure_.addOption("generalized", "Generalized", MEASURE_GENERALIZED);
    singleChannelSimilarityMeasure_.addOption("avgDifference", "Avg. Difference", MEASURE_AVG_DIFFERENCE);
    singleChannelSimilarityMeasure_.set("generalized");
    ON_CHANGE_LAMBDA(singleChannelSimilarityMeasure_, [this] {
       isoValue_.setVisibleFlag(singleChannelSimilarityMeasure_.getValue() == MEASURE_ISOCONTOURS);
    });

    addProperty(isoValue_);
    isoValue_.setVisibleFlag(false);

    addProperty(multiChannelSimilarityMeasure_);
    multiChannelSimilarityMeasure_.addOption("magnitude", "Magnitude", MEASURE_MAGNITUDE);
    multiChannelSimilarityMeasure_.addOption("angleDifference", "Angle", MEASURE_ANGLE);
    multiChannelSimilarityMeasure_.addOption("jiang", "Jiang", MEASURE_JIANG);
    multiChannelSimilarityMeasure_.addOption("crossproduct", "Crossproduct Magnitude", MEASURE_CROSSPRODUCT);
    multiChannelSimilarityMeasure_.addOption("split_channels", "Split Channels", MEASURE_SPLIT_CHANNELS);
    multiChannelSimilarityMeasure_.addOption("euclidean_norm", "Euclidean norm", MEASURE_EUCLIDEAN_NORM);
    multiChannelSimilarityMeasure_.set("euclidean_norm");
    ON_CHANGE_LAMBDA(multiChannelSimilarityMeasure_, [this] {
        weight_.setVisibleFlag(multiChannelSimilarityMeasure_.getValue() == MEASURE_JIANG);
    });

    addProperty(weight_);
    weight_.setVisibleFlag(false);

    addProperty(numSeedPoints_);
    addProperty(seedTime_);

    addProperty(clearCache_);
    ON_CHANGE_LAMBDA(clearCache_, [this] {
        tgt::FileSystem::deleteDirectoryRecursive(getCachePath());
    });
}

Processor* SimilarityMatrixCreator::create() const {
    return new SimilarityMatrixCreator();
}


bool SimilarityMatrixCreator::isReady() const {
    if (!isInitialized()) {
        setNotReadyErrorMessage("Not initialized.");
        return false;
    }
    if (!inport_.isReady()) {
        setNotReadyErrorMessage("Inport not ready.");
        return false;
    }

    // Note: Seed Mask is optional!

    if(!outport_.isReady()) {
        setNotReadyErrorMessage("Output not connected");
        return false;
    }

    return true;
}

std::vector<std::reference_wrapper<Port>> SimilarityMatrixCreator::getCriticalPorts() {
    auto criticalPorts = AsyncComputeProcessor<ComputeInput, ComputeOutput>::getCriticalPorts();
    criticalPorts.erase(std::remove_if(criticalPorts.begin(), criticalPorts.end(), [this] (const std::reference_wrapper<Port>& port){
       return port.get().getID() == seedMask_.getID();
    }), criticalPorts.end());
    return criticalPorts;
}

void SimilarityMatrixCreator::adjustPropertiesToInput() {
    //TODO: implement heuristic for auto-selecting number of seed points depending on dataset
    //numSeedPoints_.set(32768);
}

std::string SimilarityMatrixCreator::calculateHash() const {
    std::string hash;

    hash = EnsembleHash(*inport_.getData()).getHash();
    hash += seedMask_.getHash();
    hash += std::to_string(seedTime_.get());
    hash += std::to_string(numSeedPoints_.get());

    return VoreenHash::getHash(hash);
}

SimilarityMatrixCreatorInput SimilarityMatrixCreator::prepareComputeInput() {
    auto ensemble = inport_.getThreadSafeData();
    if (!ensemble)
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);

    for(const std::string& fieldName : ensemble->getCommonFieldNames()) {
        size_t numChannels = ensemble->getNumChannels(fieldName);
        if(numChannels != 1 && numChannels != 3) {
            throw InvalidInputException("Only volumes with 1 or 3 channels supported", InvalidInputException::S_ERROR);
        }
    }

    tgt::Bounds bounds = ensemble->getCommonBounds();
    if(!bounds.isDefined()) {
        throw InvalidInputException("Bounding box is empty", InvalidInputException::S_ERROR);
    }

    // Set up random generator.
    std::function<float()> rnd(
            std::bind(std::uniform_real_distribution<float>(0.0f, 1.0f), std::mt19937(seedTime_.get())));

    const VolumeBase* seedMask = seedMask_.getThreadSafeData();
    auto numSeedPoints = static_cast<size_t>(numSeedPoints_.get());
    std::vector<tgt::vec3> seedPoints;
    seedPoints.reserve(numSeedPoints);
    if (seedMask) {
        bounds.intersectVolume(seedMask->getBoundingBox().getBoundingBox());
        if(!bounds.isDefined()) {
            throw InvalidInputException("Seed Mask does not overlap with ensemble ROI", InvalidInputException::S_ERROR);
        }

        VolumeRAMRepresentationLock seedMaskLock(seedMask);
        tgt::mat4 seedMaskVoxelToWorldMatrix = seedMask->getVoxelToWorldMatrix();
        tgt::svec3 dim = seedMaskLock->getDimensions();
        for(size_t z=0; z < dim.z; z++) {
            for(size_t y=0; y < dim.y; y++) {
                for(size_t x=0; x < dim.x; x++) {
                    if(seedMaskLock->getVoxelNormalized(x, y, z) != 0.0f) {
                        tgt::vec3 pos = seedMaskVoxelToWorldMatrix * tgt::vec3(x, y, z);
                        if(bounds.containsPoint(pos)) {
                            seedPoints.emplace_back(pos);
                        }
                    }
                }
            }
        }

        if (seedPoints.empty()) {
            throw InvalidInputException("No seed points found in ROI", InvalidInputException::S_ERROR);
        }

        tgt::shuffle(seedPoints.begin(), seedPoints.end(), std::mt19937(seedTime_.get()));
        seedPoints.resize(std::min(seedPoints.size(), numSeedPoints));

        LINFO("Restricting seed points to volume mask using " << seedPoints.size() << " seeds");
    }
    else {
        // Without a seed mask, we uniformly sample the whole space enclosed by the roi.
        for (size_t k = 0; k<numSeedPoints; k++) {
            // Since argument evaluation order is unspecified in c++, we need to ensure the order manually.
            float x = rnd();
            float y = rnd();
            float z = rnd();

            tgt::vec3 seedPoint(x, y, z);
            seedPoint = bounds.getLLF() + seedPoint * bounds.diagonal();
            seedPoints.push_back(seedPoint);
        }
    }

    tgtAssert(!seedPoints.empty(), "no seed points found");
    if (seedPoints.empty()) {
        throw InvalidInputException("No seed points found", InvalidInputException::S_ERROR);
    }

    std::unique_ptr<SimilarityMatrixList> outputMatrices(new SimilarityMatrixList(*ensemble));

    return SimilarityMatrixCreatorInput{
            std::move(ensemble),
            std::move(outputMatrices),
            std::move(seedPoints),
            singleChannelSimilarityMeasure_.getValue(),
            isoValue_.get(),
            multiChannelSimilarityMeasure_.getValue(),
            weight_.get(),
            calculateHash()
    };
}

SimilarityMatrixCreatorOutput SimilarityMatrixCreator::compute(SimilarityMatrixCreatorInput input, ProgressReporter& progress) const {

    std::unique_ptr<SimilarityMatrixList> similarityMatrices = std::move(input.outputMatrices);

    // First, try to load from cache.
    std::string cachePath = getCachePath() + "/" + input.hash + ".vsm";
    if(tgt::FileSystem::fileExists(cachePath)) {
        try {
            std::ifstream stream(cachePath);
            if (!stream)
                throw tgt::CorruptedFileException("Could not read file: " + cachePath);

            JsonDeserializer json;
            json.read(stream, true);
            Deserializer s(json);
            s.deserialize("similarity", *similarityMatrices);
            LINFO("Cache loaded successfully!");

            return SimilarityMatrixCreatorOutput{
                    std::move(similarityMatrices)
            };
        } catch (std::exception& e) {
            LWARNING("Could not load cache: " << e.what());
        }
    }
    else {
        LINFO("No cache file found.");
    }

    // Create similarity matrix.
    std::vector<tgt::vec3> seedPoints = std::move(input.seedPoints);

    progress.setProgress(0.0f);

    const std::vector<std::string>& fieldNames = input.ensemble->getCommonFieldNames();
    for (size_t fi=0; fi<fieldNames.size(); fi++) {

        const std::string& fieldName = fieldNames[fi];
        size_t numChannels = input.ensemble->getNumChannels(fieldName);
        tgt::vec2 valueRange;
        if(numChannels == 1) {
             valueRange = input.ensemble->getValueRange(fieldName);
        }
        else {
            // If we use multi-channel volumes, we need to calculate the min. and max. magnitude in order
            // to properly scale values later on to generate a matrix whose values are within [0, 1].
            valueRange = input.ensemble->getMagnitudeRange(fieldName);
        }

        // Init empty flags.
        const size_t numElements = input.ensemble->getTotalNumTimeSteps() * seedPoints.size() * numChannels;
#ifdef USE_MEMORY_MAPPED_FILES
        std::unique_ptr<DiskArrayStorage<float>> flagStorage;
        DiskArray<float> Flags;
        auto storeValue = [&] (float value) { flagStorage->storeElement(value); };
#else
        std::vector<float> Flags;
        Flags.reserve(numElements);
        auto storeValue = [&] (float value) { Flags.push_back(value); };
#endif

        std::string tmpPath =  VoreenApplication::app()->getUniqueTmpFilePath();
        std::string flagFile = input.hash + "_" + fieldName  + ".flags";

        std::string cachePath = tgt::FileSystem::cleanupPath(getCachePath() + "/" + flagFile);
        bool cachedFileFound = VoreenApplication::app()->useCaching() && tgt::FileSystem::fileExists(cachePath);

        // Load cache file, if found one.
        if(cachedFileFound) {
            LINFO("Found cached flag file for field " << fieldName);
#ifdef USE_MEMORY_MAPPED_FILES
            // Reuse memory mapped file. We need to create a copy because the file
            // gets removed after DiskArrayStorage is destructed.
            tgt::FileSystem::copyFile(cachePath, tmpPath);
            flagStorage.reset(new DiskArrayStorage<float>(tmpPath, numElements));
            Flags = flagStorage->asArray();
#else
            Flags.resize(numElements); // Actually set size.
            std::fstream file(cachePath, std::ios_base::binary | std::ios_base::in);
            file.read(reinterpret_cast<char*>(Flags.data()), Flags.size() * sizeof(float) / sizeof(char));
            file.close();
#endif
            progress.setProgress(0.1f * (fi+1) / fieldNames.size());
        }
        else { // Or freshly create Flag array.

            LINFO("Creating flag file for " << fieldName);

#ifdef USE_MEMORY_MAPPED_FILES
            flagStorage.reset(new DiskArrayStorage<float>(tmpPath));
#endif

            SubtaskProgressReporter memberProgressReporter(progress, tgt::vec2(fi, fi+0.7f) / tgt::vec2(fieldNames.size()));
            float progressPerTimeStep = 1.0f / (input.ensemble->getTotalNumTimeSteps());
            size_t index = 0;
            for (const EnsembleMember& member : input.ensemble->getMembers()) {
                for (const TimeStep& timeStep : member.getTimeSteps()) {

                    const VolumeBase* volume = timeStep.getVolume(fieldName);
                    if(!volume) {
                        LERROR("Could not load volume " << timeStep.getURL(fieldName).getURL() << ". Filling with " << MissingValue);
                        size_t numMissingValues = seedPoints.size() * numChannels;
                        for(size_t i=0; i<numMissingValues; i++) {
                            storeValue(MissingValue);
                        }
                        continue;
                    }

                    tgt::Bounds bounds = volume->getBoundingBox().getBoundingBox();
                    tgt::mat4 worldToVoxelMatrix = volume->getWorldToVoxelMatrix();
                    RealWorldMapping rwm = volume->getRealWorldMapping();

                    VolumeRAMRepresentationLock lock(volume);
                    for (const tgt::vec3& seedPoint : seedPoints) {

                        // Skip points outside bounds.
                        if(!bounds.containsPoint(seedPoint)) {
                            for (size_t channel = 0; channel < numChannels; channel++) {
                                storeValue(MissingValue);
                            }
                            continue;
                        }

                        tgt::vec3 pos = worldToVoxelMatrix * seedPoint;
                        for (size_t channel = 0; channel < numChannels; channel++) {
                            float value = lock->getVoxelNormalizedLinear(pos, channel);
                            value = rwm.normalizedToRealWorld(value);
                            storeValue(value);
                        }
                    }

                    // Update progress.
                    index++;
                    memberProgressReporter.setProgress(index * progressPerTimeStep);
                }
            }

#ifdef USE_MEMORY_MAPPED_FILES
            Flags = flagStorage->asArray();
#endif

            // If caching is enabled, store the Flag file in the cache directory.
            if (VoreenApplication::app()->useCaching()) {
                tgt::FileSystem::createDirectoryRecursive(tgt::FileSystem::dirName(cachePath));
#ifdef USE_MEMORY_MAPPED_FILES
                // Once we are done, copy the tmp file to the cache folder.
                try {
                    tgt::FileSystem::copyFile(tmpPath, cachePath);
                }
                catch (tgt::FileException& e) {
                    LWARNING("Could not store Cache file of field " << fieldName << " - " << e.what());
                }
#else
                std::fstream file(cachePath, std::ios_base::binary | std::ios_base::out);
                file.write(reinterpret_cast<char*>(Flags.data()), Flags.size() * sizeof(float) / sizeof(char));
                file.close();
#endif
            }
        }

        ////////////////////////////////////////////////////////////
        //Calculate distances for up-right corner and reflect them//
        ////////////////////////////////////////////////////////////

        LINFO("Calculating Distance Matrix for " << fieldName);

        auto index = [&] (size_t timeStepIndex, size_t seedIndex, size_t channel = 0) {
            return timeStepIndex * seedPoints.size() * numChannels + seedIndex * numChannels + channel;
        };

        SimilarityMatrix& distanceMatrix = similarityMatrices->getSimilarityMatrix(fieldName);
        long size = static_cast<long>(distanceMatrix.getSize());

        SubtaskProgressReporter flagsProgress(progress,tgt::vec2((fi+0.7f), (fi+1.0f)) / tgt::vec2(fieldNames.size()));
        ThreadedTaskProgressReporter threadedProgress(flagsProgress, size);
        bool aborted = false;

#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (long i = 0; i < size; i++) {
#ifdef VRN_MODULE_OPENMP
            if(aborted) {
                continue;
            }
#else
            if(aborted) {
                break;
            }
#endif

            for (long j = 0; j <= i; j++) {
                if(numChannels == 1 || input.multiChannelSimilarityMeasure == MEASURE_MAGNITUDE || input.multiChannelSimilarityMeasure == MEASURE_SPLIT_CHANNELS) {

                    float intersectionSamples = 0.0f;
                    float unionSamples = 0.0f;

                    // If we decide to split the channels, we consider each channel as flag.
                    size_t numFlags = seedPoints.size();
                    if(input.multiChannelSimilarityMeasure == MEASURE_SPLIT_CHANNELS) {
                        numFlags *= numChannels;
                    }

                    for (size_t k = 0; k < numFlags; k++) {

                        float a = 0.0f;
                        float b = 0.0f;

                        if(numChannels > 1 && input.multiChannelSimilarityMeasure == MEASURE_MAGNITUDE) {
                            // Calculate length.
                            for (size_t channel = 0; channel < numChannels; channel++) {
                                float flagA = Flags[index(i, k, channel)];
                                a += flagA * flagA;

                                float flagB = Flags[index(j, k, channel)];
                                b += flagB * flagB;
                            }

                            a = std::sqrt(a);
                            b = std::sqrt(b);
                        }
                        else if(numChannels > 1 && input.multiChannelSimilarityMeasure == MEASURE_SPLIT_CHANNELS) {
                            size_t newK = k/numChannels;                //k of the current channel
                            size_t channel = k / (numFlags/numChannels);
                            a = Flags[index(i, newK, channel)];
                            b = Flags[index(j, newK, channel)];
                        }
                        else {
                            a = Flags[index(i, k)];
                            b = Flags[index(j, k)];
                        }

                        // Normalize range to interval [0, 1].
                        a = mapRange(a, valueRange.x, valueRange.y, 0.0f, 1.0f);
                        b = mapRange(b, valueRange.x, valueRange.y, 0.0f, 1.0f);

                        if(input.singleChannelSimilarityMeasure == MEASURE_AVG_DIFFERENCE) {
                            intersectionSamples = 1.0f - std::abs(a - b);
                            unionSamples += 1.0f;
                        }
                        else {
                            if (input.singleChannelSimilarityMeasure == MEASURE_ISOCONTOURS) {
                                a = a < input.isoValue ? 1.0f : 0.0f;
                                b = b < input.isoValue ? 1.0f : 0.0f;
                            }

                            intersectionSamples += (1.0f - std::max(a, b));
                            unionSamples += (1.0f - std::min(a, b));
                        }
                    }

                    if (unionSamples > 0.0f)
                        distanceMatrix(i, j) = (unionSamples - intersectionSamples) / unionSamples;
                    else
                        distanceMatrix(i, j) = 1.0f;
                }
                else {

                    Statistics statistics(false);

                    for (size_t k = 0; k < seedPoints.size(); k++) {

                        tgt::vec4 vector_i = tgt::vec4::zero;
                        tgt::vec4 vector_j = tgt::vec4::zero;

                        for (size_t channel = 0; channel < numChannels; channel++) {
                            vector_i[channel] = Flags[index(i, k, channel)];
                            vector_j[channel] = Flags[index(j, k, channel)];
                        }

                        if(input.multiChannelSimilarityMeasure == MEASURE_ANGLE) {
                            if (vector_i != tgt::vec4::zero && vector_j != tgt::vec4::zero) {
                                tgt::vec4 normVector_i = tgt::normalize(vector_i);
                                tgt::vec4 normVector_j = tgt::normalize(vector_j);

                                float dot = tgt::dot(normVector_i, normVector_j);
                                float angle = std::acos(tgt::clamp(dot, -1.0f, 1.0f)) / tgt::PIf;
                                if(!tgt::isNaN(angle)) {
                                    statistics.addSample(angle);
                                }
                                else {
                                    tgtAssert(false, "NaN value");
                                }

                            }
                            else if (vector_i == tgt::vec4::zero && vector_j == tgt::vec4::zero) {
                                statistics.addSample(0.0f);
                            }
                            else {
                                statistics.addSample(1.0f);
                            }
                        }
                        else if(input.multiChannelSimilarityMeasure == MEASURE_JIANG) {
                            float a = tgt::length(vector_i);
                            float b = tgt::length(vector_j);

                            if (a > 0.0f && b > 0.0f) {
                                tgt::vec4 normVector_i = vector_i / a;
                                tgt::vec4 normVector_j = vector_j / b;

                                float dot = tgt::dot(normVector_i, normVector_j);
                                float angle = std::asin(tgt::clamp(dot, -1.0f, 1.0f));
                                tgtAssert(!tgt::isNaN(angle), "NaN value");

                                // We don't use the lower bound of the value range on purpose here!
                                float magnitude = mapRange(std::abs(a - b), 0.0f, valueRange.y, 0.0f, 1.0f);
                                statistics.addSample(1.0f - ((1.0f - input.weight) * std::exp(-magnitude) + input.weight * std::exp(-2.0f*angle)));
                            }
                            else if (a == 0.0f && b == 0.0f) {
                                statistics.addSample(0.0f);
                            }
                            else { // Exactly one vector was zero.

                                // We add a 'maximally different' sample (which leads, however, to a discontinuity).
                                statistics.addSample(1.0f);
                                // Instead, we could also fallback to the magnitude difference
                                //statistics.addSample(std::abs(a-b)/valueRange.y);
                            }
                        }
                        else if(input.multiChannelSimilarityMeasure == MEASURE_CROSSPRODUCT) {
                            if (vector_i == tgt::vec4::zero && vector_j == tgt::vec4::zero) {
                                statistics.addSample(0.0f);
                            }
                            else if (vector_i != tgt::vec4::zero && vector_j != tgt::vec4::zero) {
                                // Normalize vectors according to max magnitude within data set.
                                tgt::vec3 a = vector_i.xyz() / valueRange.y;
                                tgt::vec3 b = vector_j.xyz() / valueRange.y;

                                float area = tgt::length(tgt::cross(a, b));
                                // In case area is 0, we have to account for colinear vectors.
                                if(area < std::numeric_limits<float>::epsilon()) {
                                    float length_a = tgt::length(a);
                                    float length_b = tgt::length(b);

                                    tgt::vec3 normA = a / length_a;
                                    tgt::vec3 normB = b / length_b;

                                    // Determine direction of collinearity.
                                    float dot = tgt::dot(normA, normB);
                                    float angle = std::acos(tgt::clamp(dot, -1.0f, 1.0f));
                                    if(angle > tgt::PIf*0.5f) {
                                        statistics.addSample(tgt::abs(length_a + length_b) * 0.5f);
                                    }
                                    else {
                                        statistics.addSample(tgt::abs(length_a - length_b) * 0.5f);
                                    }
                                }
                                else {
                                    statistics.addSample(area);
                                }
                            }
                            else {
                                statistics.addSample(0.0f);
                            }
                        }
                        else if(input.multiChannelSimilarityMeasure == MEASURE_EUCLIDEAN_NORM) {
                            statistics.addSample(tgt::length(vector_i - vector_j) / (2*valueRange.y));
                        }
                    }

                    distanceMatrix(i, j) = statistics.getMean();
                    //DistanceMatrix(i, j) = statistics.getMedian(); // Needs collecting samples enabled
                    //DistanceMatrix(i, j) = statistics.getRelStdDev();
                }
            }


            if (threadedProgress.reportStepDone()) {
#ifdef VRN_MODULE_OPENMP
                #pragma omp critical
                aborted = true;
#else
                aborted = true;
                break;
#endif
            }
        }

        if (aborted) {
            throw boost::thread_interrupted();
        }
    }

    // Save to cache.
    try {
        std::ofstream stream(cachePath);
        JsonSerializer json;
        Serializer s(json);
        s.serialize("similarity", *similarityMatrices);
        json.write(stream, false, true);
        LINFO("Cache saved successfully!");
    } catch(std::exception& e) {
        LWARNING("Could not store cache: " << e.what());
    }

    progress.setProgress(1.0f);

    return SimilarityMatrixCreatorOutput{
            std::move(similarityMatrices)
    };
}

void SimilarityMatrixCreator::processComputeOutput(ComputeOutput output) {
    outport_.setData(output.outputMatrices.release(), true);
}


} // namespace
