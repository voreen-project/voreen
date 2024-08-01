/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "patchtrainingdataextractor.h"

#include "voreen/core/datastructures/volume/volumebase.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "../util/cosinetransform.h"
#include "modules/hdf5/utils/hdf5utils.h"

#include "tgt/stopwatch.h"

#include <algorithm>

#include "../util/trainingdataextractionthread.h"

namespace voreen {

const std::string PatchTrainingDataExtractor::loggerCat_("voreen.nucleusdetection.PatchTrainingDataExtractor");

PatchTrainingDataExtractor::PatchTrainingDataExtractor()
    : PatchFeatureExtractor()
    , dataInport_(Port::INPORT, "volumedata.input", "Volume Data Input", false)
    , labelInport_(Port::INPORT, "volumelabels.input", "Label Volume Input", false)
    , intensityReferenceVolume_(Port::INPORT, "intensityreference.input", "Intensity Reference Volume Input", false)
    , filename_("outputFilename", "Training Data Output File", "Select file...", "",
            "HDF5 (*.h5)", FileDialogProperty::SAVE_FILE, Processor::INVALID_PATH)
    , foregroundModulus_("foregroundmodulus", "n (each n-th sample is a foreground sample)", 2, 2, 25)
    , capBackground_("capBackground", "Cap Background", false)
    , numSamplesPerFile_("numsamples", "Number of Samples per HDF5-File", 20000, 1000, 500000)
    , computeButton_("train", "Extract Training Data")
    , autoCompute_("autoCompute", "Auto Compute")
    , timeEstimation_("timeestimation", "Runtime Estimation", " - ", Processor::VALID)
    , startComputation_(false)
{
    addPort(dataInport_);
    addPort(labelInport_);
    addPort(intensityReferenceVolume_);

    computeButton_.onChange(MemberFunctionCallback<PatchTrainingDataExtractor>(this, &PatchTrainingDataExtractor::computePatches));

    timeEstimation_.setEditable(false);
    //timeEstimation_.setInstantUpdate(true);
    addProperty(timeEstimation_);

    addProperty(filename_);
    addProperty(foregroundModulus_);
    addProperty(capBackground_);
    addProperty(numSamplesPerFile_);

    addProperty(computeButton_);
    addProperty(autoCompute_);
}

Processor* PatchTrainingDataExtractor::create() const {
    return new PatchTrainingDataExtractor();
}

void PatchTrainingDataExtractor::process() {
    if (autoCompute_.get())
        computePatches();
}

void PatchTrainingDataExtractor::computePatches() {

    if (!isInitialized())
        return;

    timeEstimation_.set(" - ");

    if (!dataInport_.isReady() || !dataInport_.getData() || dataInport_.getData()->empty()) {
        LERROR(" no volume ");
        return;
    }

    if (!labelInport_.isReady() || !labelInport_.getData() || labelInport_.getData()->empty()) {
        LERROR(" no labels ");
        return;
    }

    if (!checkForRequiredFilterKernels()) {
        LERROR(" Missing filter kernels! ");
        return;
    }

    if (!(dataInport_.getData()->size() == labelInport_.getData()->size())) {
        LERROR(" Different number of data and label volumes! ");
        return;
    }

    // start
    deactivateProperties();
    setProgress(0.f);

    size_t patchSize = static_cast<size_t>(patchSize_.getValue());

    // get required filter kernels
    std::vector<std::pair<size_t,float*> > filterKernels = getFilterKernels(patchSize);
    if (filterKernels.size() != 4) {
        LERROR("Cannot load filter kernels!");
        deleteFilterKernels(filterKernels);
        reactivateProperties();
        return;
    }

    LINFO("Extracting training data...");
    tgt::Stopwatch timer;
    timer.start();

    // compute number of features and offset indices for features
    size_t numFeatures = computeNumFeatures(filterKernels);

    if (!numFeatures) {
        LERROR("No features selected!");
        deleteFilterKernels(filterKernels);
        reactivateProperties();
        return;
    }

    // if necessary, pre-compute coefficients for 2D cosine transform
    float* dctCoefficients = 0;
    if (cosineTransform1D_.getValue() || cosineTransform2D_.getValue()) {
        dctCoefficients = new float[patchSize * patchSize];                
        generateDCTCoefficients(patchSize, dctCoefficients);
    }

    tgt::bvec4 necessaryScales = determineNecessaryScales();
    size_t border = computeNecessaryBorder(necessaryScales);


    // compute the total number of voxels that have to be processed over all volumes for progress information
    double numVoxelsTotal = 0.0;
    double numVoxelsProcessed = 0.0;
    for (size_t volumeIndex = 0; volumeIndex < dataInport_.getData()->size(); ++volumeIndex) {
        // get dimensions of this volume
        tgt::svec3 volDim = dataInport_.getData()->at(volumeIndex)->getDimensions();
        // subtract border
        for (size_t dim = 0; dim < 3; ++dim) 
            volDim.elem[dim] = (volDim.elem[dim] > 2 * border) ? volDim.elem[dim] - 2*border : 0;
        size_t voxelsInVolume = tgt::hmul(volDim);
        numVoxelsTotal += static_cast<double>(voxelsInVolume);
    }


    // get number of CPU cores to start the right amount of threads for feature extraction
    size_t numCores = static_cast<size_t>(boost::thread::hardware_concurrency());
    if (numCores == 0) {
        numCores = 1;
        LWARNING("Could not detect number of CPU cores. Using (slow) single threaded feature extraction.");
    }
    else 
        LINFO("Using " << numCores << " threads for feature extraction");
        
    // timer for estimating the remaining time
    tgt::Stopwatch estimationTimer;        
    estimationTimer.start();

    for (size_t volumeIndex = 0; volumeIndex < dataInport_.getData()->size(); ++volumeIndex) {

        tgt::svec3 volDim = dataInport_.getData()->at(volumeIndex)->getDimensions();
        const VolumeRAM* dataVolume = dataInport_.getData()->at(volumeIndex)->getRepresentation<VolumeRAM>();
        const VolumeRAM* labelVolume = labelInport_.getData()->at(volumeIndex)->getRepresentation<VolumeRAM>();

        if (!(dataInport_.getData()->at(volumeIndex)->getDimensions() == labelInport_.getData()->at(volumeIndex)->getDimensions())) {
            LERROR(" Dimensions of label volume and input volume " << volumeIndex << " are different! ");
            deleteFilterKernels(filterKernels);
            reactivateProperties();
            delete[] dctCoefficients;
            return;
        }

        if (!dataVolume || !labelVolume) {
            LERROR("No VolumeRAM available!");
            deleteFilterKernels(filterKernels);
            reactivateProperties();
            delete[] dctCoefficients;
            return;
        }

        // usually data is normalized over uint16_t range, but we want to normalize over the volume's histogram to account for 
        // variation of fluorescence intensity across different data sets -> use histogram information to scale values
        tgt::vec2 scalingMinMax = tgt::vec2(0.f, 1.f);
        // if present, we use the intensity reference volume
        if (intensityReferenceVolume_.getData())
            scalingMinMax = getScalingMinMaxIntensity(intensityReferenceVolume_.getData());
        else
            scalingMinMax = getScalingMinMaxIntensity(dataInport_.getData()->at(volumeIndex));
 
        // 1. compute all positions of foreground and background voxels
        std::vector<tgt::svec3> foregroundPositions;
        std::vector<tgt::svec3> backgroundPositions;

        for (tgt::svec3 pos = tgt::svec3(border); pos.z < volDim.z - border; ++pos.z) {
            for (pos.y = border; pos.y < volDim.y - border; ++pos.y) {
                for (pos.x = border; pos.x < volDim.x - border; ++pos.x) {
                    float currentLabel = labelVolume->getVoxelNormalized(pos);
                    if (currentLabel ==  0.f)
                        backgroundPositions.push_back(pos);
                    else
                        foregroundPositions.push_back(pos); 
                }
            }
        } 

        // 2. randomly shuffle the positions of each class
        std::random_shuffle(foregroundPositions.begin(), foregroundPositions.end());
        std::random_shuffle(backgroundPositions.begin(), backgroundPositions.end());

        // 2b. If necessary, throw away background voxels to cap to the number of foreground voxels
        if (!foregroundPositions.empty() && capBackground_.get()) {
            size_t cap = foregroundPositions.size();
            backgroundPositions.resize(cap);
        }

        // 3. determine parameters to write out the data into several files
        size_t numForegroundSamples = foregroundPositions.size();
        size_t numBackgroundSamples = backgroundPositions.size();

        // number of samples per file
        size_t numSamples = static_cast<size_t>(numSamplesPerFile_.get());

        // determine number of files (each file should be balanced and we want to use all of the data)
        size_t numBlocks = static_cast<size_t>(std::ceil(static_cast<double>(numBackgroundSamples) / static_cast<double>(foregroundModulus_.get() - 1)));   
        size_t numFiles = numBlocks * static_cast<size_t>(foregroundModulus_.get()) / numSamples;
        if ((numBlocks * static_cast<size_t>(foregroundModulus_.get())) % numSamples != 0)
            numFiles++;
            
        if (numBackgroundSamples == 0) {
            numFiles = 0;
            LWARNING(" Cannot export training data for data set " << volumeIndex << ": found no background samples!");
        }

        // determine base file name
        std::string baseFileName = tgt::FileSystem::fullBaseName(tgt::FileSystem::cleanupPath(filename_.get()));

        size_t currentForegroundIndex = 0;
        size_t currentBackgroundIndex = 0;

        // create a data matrix and label vector which will be filled for each file   
        float* featureData = new float[numSamples * numFeatures];
        float* labelData = new float[numSamples];

        size_t foregroundModulus = static_cast<size_t>(foregroundModulus_.get()); 
        
        for (size_t fileNumber = 0; fileNumber < numFiles; ++fileNumber) {
            // positions to evaluate for the current file
            std::vector<tgt::svec3> filePositions;

            // fill list of positions for this file
            for (size_t position = 0; position < numSamples; ++position) {
                // get the current position
                if ((position % foregroundModulus == 0) && (numForegroundSamples > 0)) { // foreground
                    filePositions.push_back(foregroundPositions.at(currentForegroundIndex));
                    currentForegroundIndex = (currentForegroundIndex + 1) % numForegroundSamples;
                }
                else { // background
                    filePositions.push_back(backgroundPositions.at(currentBackgroundIndex));
                    currentBackgroundIndex = (currentBackgroundIndex + 1) % numBackgroundSamples;
                }
            }

            // now we want to compute the features for the given positions
            std::vector<Sample> foregroundVoxels;
            std::vector<Sample> backgroundVoxels;

            // compute number of positions that each thread has to process
            size_t voxelPackageSize = filePositions.size() / numCores;
            if (voxelPackageSize == 0)
                voxelPackageSize = 1;

            // start threads to extract the features
            std::vector<TrainingDataExtractionThread*> workerThreads;
            size_t startVoxel = 0;
            size_t lastVoxel = filePositions.size() - 1;
            for (size_t threadCount = 0; threadCount < numCores; ++threadCount) {
                // compute the actual slice range for this thread
                size_t localStartVoxel = startVoxel + voxelPackageSize * threadCount;
                if (localStartVoxel > lastVoxel)
                    break;
                size_t localEndVoxel = (threadCount == numCores - 1) ? lastVoxel : startVoxel + voxelPackageSize * (threadCount + 1) - 1;
                if (localEndVoxel > lastVoxel)
                    localEndVoxel = lastVoxel;
                TrainingDataExtractionThread* thread = new TrainingDataExtractionThread(this, volDim, localStartVoxel, localEndVoxel, filePositions, border, 
                                                                                numFeatures, dataVolume, labelVolume, filterKernels, dctCoefficients, necessaryScales, scalingMinMax);
                workerThreads.push_back(thread);
                thread->run();
            }

            // join each thread, retrieve its data and delete it
            for (auto& t : workerThreads) {
                t->join();

                // copy data (std::move is required to prevent mysterious memory leaks)
                std::vector<Sample>& localForegroundVoxels = t->getForegroundVoxels();
                for (auto& s : localForegroundVoxels)
                    foregroundVoxels.push_back(std::move(s));
                std::vector<Sample>& localBackgroundVoxels = t->getBackgroundVoxels();
                for (auto& s : localBackgroundVoxels)
                    backgroundVoxels.push_back(std::move(s));

                delete t;
            }

            // fill matrix for writing out the data
            Sample currentSample(numFeatures);
            float currentLabel;
            size_t localForegroundIndex = 0;
            size_t localBackgroundIndex = 0;

            tgtAssert(numSamples == foregroundVoxels.size() + backgroundVoxels.size(), "wrong number of samples");

            // we know that the data is balanced in the desired way
            for (size_t sample = 0; sample < numSamples; ++sample) {
                // get the current sample
                if ((sample % foregroundModulus == 0) && (numForegroundSamples > 0)) { // foreground
                    tgtAssert(localForegroundIndex < foregroundVoxels.size(), "not enough foreground samples");

                    currentSample = foregroundVoxels.at(localForegroundIndex);
                    currentLabel = 1.f;
                    localForegroundIndex++;
                }
                else { // background
                    tgtAssert(localBackgroundIndex < backgroundVoxels.size(), "not enough background voxels");

                    currentSample = backgroundVoxels.at(localBackgroundIndex);
                    currentLabel = 0.f;
                    localBackgroundIndex++;
                }
                // copy the feature vector and label
                std::memcpy(&featureData[sample * numFeatures], &currentSample.features_[0], numFeatures * sizeof(float));
                labelData[sample] = currentLabel;
            }

            std::stringstream fileNameStream;
            fileNameStream << baseFileName << "_" << volumeIndex << "_" << fileNumber << ".h5";

            using namespace H5;
            const H5std_string  FILE_NAME( fileNameStream.str() );
            const H5std_string  DATA_NAME( "data" );
            const H5std_string  OUTPUT_NAME( "label" );
            const int featureLength = numFeatures;
            const int sampleSize =  numSamples;
            const int outputLength = 1; // labels->cols

            H5File file(FILE_NAME, H5F_ACC_TRUNC);
            hsize_t dimsf[2];
            dimsf[0] = sampleSize;
            dimsf[1] = featureLength;
            DataSpace dataSpace(2,dimsf);
            dimsf[1] = outputLength;
            DataSpace outputSpace(2,dimsf);
            FloatType dataType(PredType::NATIVE_FLOAT);
            DataSet featureSet = file.createDataSet( DATA_NAME, dataType, dataSpace);
            DataSet labelSet = file.createDataSet(OUTPUT_NAME, dataType, outputSpace);
            featureSet.write(featureData, PredType::NATIVE_FLOAT);
            labelSet.write(labelData,PredType::NATIVE_FLOAT);
            file.close();
        }

        delete[] featureData;
        delete[] labelData;

        // we have processed another volume -> we add its voxels to the number of processed voxels and compute the new time estimation and progress information
        tgt::svec3 processedDimensions = volDim - tgt::svec3(2*border);
        numVoxelsProcessed += (double) tgt::hmul(processedDimensions);
        double processedPart = numVoxelsProcessed / numVoxelsTotal;
        double currentTime = estimationTimer.getRuntime() / (1000.0 * 60.0); // runtime in minutes
        double estimate = currentTime * ((numVoxelsTotal-numVoxelsProcessed) / numVoxelsProcessed);
        std::stringstream strstr;
        if (tgt::iround(estimate) < 1)
            strstr << "< 1 min";
        else
            strstr << tgt::iround(estimate) << " min";
        timeEstimation_.set(strstr.str());
        setProgress(std::min(0.99f, static_cast<float>(processedPart)));
    }

    // filter functions are not needed anymore
    deleteFilterKernels(filterKernels);
    delete[] dctCoefficients;

    estimationTimer.stop();
    timer.stop();
    LINFO("Time for generating and exporting training data: " << (float) timer.getRuntime() / 1000.f << " seconds");

    setProgress(1.f);
    timeEstimation_.set("Finished");

    reactivateProperties();
}

void PatchTrainingDataExtractor::deactivateProperties() {
    PatchFeatureExtractor::deactivateProperties();
    filename_.setReadOnlyFlag(true);
    foregroundModulus_.setReadOnlyFlag(true);
    numSamplesPerFile_.setReadOnlyFlag(true);
    computeButton_.setReadOnlyFlag(true);
    autoCompute_.setReadOnlyFlag(true);
}

void PatchTrainingDataExtractor::reactivateProperties() {
    PatchFeatureExtractor::reactivateProperties();
    filename_.setReadOnlyFlag(false);
    foregroundModulus_.setReadOnlyFlag(false);
    numSamplesPerFile_.setReadOnlyFlag(false);
    computeButton_.setReadOnlyFlag(false);
    autoCompute_.setReadOnlyFlag(false);
}

} // namespace
