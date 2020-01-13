/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "patchcaffeclassifier.h"

#include "voreen/core/datastructures/volume/volumebase.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/histogram.h"

#include "tgt/stopwatch.h"

#include <algorithm>

#include <caffe/caffe.hpp>

#include "../util/cosinetransform.h"
#include "../util/featureextractionthread.h"

namespace voreen {

const std::string PatchCaffeClassifier::loggerCat_("voreen.nucleusdetection.PatchCaffeClassifier");

PatchCaffeClassifier::PatchCaffeClassifier()
    : PatchFeatureExtractor()
    , dataInport_(Port::INPORT, "volumehandle.input", "Volume Input", false)
    , classifierOutput_(Port::OUTPORT, "classifier.output", "Classifier Output", false)
    , protoFilename_("protoFilename", "Prototxt File Name", "Select file...", "",
            "Prototxt (*.prototxt)", FileDialogProperty::OPEN_FILE, Processor::INVALID_PATH)
    , modelFilename_("modelFilename", "Model File Name", "Select file...", "",
            "Caffe Model (*.caffemodel)", FileDialogProperty::OPEN_FILE, Processor::INVALID_PATH)
    , computeButton_("classifyButton", "Classify Data")
    , timeEstimation_("timeestimation", "Runtime Estimation", " - ", Processor::VALID)
    , autoCompute_("autoCompute", "Auto Compute")
    , startComputation_(false)
{
    addPort(dataInport_);

    addPort(classifierOutput_);

    computeButton_.onChange(MemberFunctionCallback<PatchCaffeClassifier>(this, &PatchCaffeClassifier::computePatches));

    addProperty(protoFilename_);
    addProperty(modelFilename_);

    addProperty(computeButton_);
    
    timeEstimation_.setEditable(false);
    //timeEstimation_.setInstantUpdate(true);
    addProperty(timeEstimation_);

    addProperty(autoCompute_);
}

Processor* PatchCaffeClassifier::create() const {
    return new PatchCaffeClassifier();
}

void PatchCaffeClassifier::process() {
    if (autoCompute_.get())
        computePatches();
}

void PatchCaffeClassifier::computePatches() {

    if (!isInitialized())
        return;
    
    timeEstimation_.set(" - ");

    classifierOutput_.clear();

    if (!dataInport_.isReady()) {
        LERROR(" no volume ");
        return;
    }

    if (!checkForRequiredFilterKernels()) {
        LERROR(" Missing filter kernels! ");
        return;
    }

    const VolumeRAM* dataVolume = dataInport_.getData()->getRepresentation<VolumeRAM>();

    if (!dataVolume) {
        LERROR("No VolumeRAM available!");
        return;
    }

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

    tgt::svec3 volDim = dataInport_.getData()->getDimensions();

    LINFO("Classifying...");
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
    
    if (tgt::hor(tgt::lessThanEqual(tgt::ivec3(volDim) - tgt::ivec3(2 * border), tgt::ivec3::zero))) {
        LERROR("Input volume is too small for classification");
        deleteFilterKernels(filterKernels);
        delete[] dctCoefficients;
        reactivateProperties();
    }

    // create output volume
    VolumeRAM_2xFloat* volumeram = new VolumeRAM_2xFloat(tgt::svec3((volDim.x - 2*border) , (volDim.y - 2*border) , (volDim.z - 2*border)));
    Volume* volume = new Volume(volumeram, dataInport_.getData()->getSpacing(), dataInport_.getData()->getOffset());
    tgt::svec3 outputDim = volDim - tgt::svec3(2 * border);
    size_t numOutputVoxels = tgt::hmul(outputDim);  // = (volDim.x - 2*border) * (volDim.y - 2*border) * (volDim.z - 2*border);
    tgt::vec2* outputBuffer = reinterpret_cast<tgt::vec2*>(volumeram->getData());

    // set GPU mode for Caffe
    caffe::Caffe::set_mode(caffe::Caffe::GPU);
    int device_id = 0;
    caffe::Caffe::SetDevice(device_id);

    // load classifier
    std::shared_ptr<caffe::Net<float>> net;
    net.reset(new caffe::Net<float>(protoFilename_.get(), caffe::TEST));
    net->CopyTrainedLayersFrom(modelFilename_.get());
    
    caffe::Blob<float>* inputLayer = net->input_blobs()[0];
    size_t inputBatchSize = static_cast<size_t>(inputLayer->shape(0));

    // determine number of batches for classification
    size_t numBatches = numOutputVoxels / inputBatchSize;
    if (numOutputVoxels % inputBatchSize != 0)
        numBatches++;

    // usually data is normalized over uint16_t range, but we want to normalize over the volume's histogram to account for 
    // variation of fluorescence intensity across different data sets -> use histogram information to scale values
    tgt::vec2 scalingMinMax = getScalingMinMaxIntensity(dataInport_.getData());

    // memory for a batch that is given as input to the classifier and the corresponding results
    float* currentBatch = new float[numFeatures * inputBatchSize];

    // get number of CPU cores to start the right amount of threads for feature extraction
    size_t numCores = static_cast<size_t>(boost::thread::hardware_concurrency());
    if (numCores == 0) {
        numCores = 1;
        LWARNING("Could not detect number of CPU cores. Using (slow) single threaded feature extraction.");
    }
    else 
        LINFO("Using " << numCores << " threads for feature extraction");

    // determine local batch size for the threads (last thread might differ!)
    size_t localBatchSize = inputBatchSize / numCores;
    size_t lastThreadBatchSize = inputBatchSize - (localBatchSize * (numCores - 1));

    double avgFeatureExtractionTime = 0.0;
    double avgClassificationTime = 0.0;
    tgt::Stopwatch performanceTimer;
    tgt::Stopwatch estimationTimer;
    double timeEstimation = 0.0;
    
    // Classify data: iterate over batches
    for (size_t currentBatchNumber = 0; currentBatchNumber < numBatches; ++currentBatchNumber) {
        estimationTimer.start();
        performanceTimer.start();
        
        size_t currentIndexOffset = currentBatchNumber * inputBatchSize;

        // start multiple threads to extract the feature vectors of all voxels
        std::vector<FeatureExtractionThread*> workerThreads;
        for (size_t threadCount = 0; threadCount < numCores; ++threadCount) {
            size_t localNumVoxels = (threadCount == numCores - 1) ? lastThreadBatchSize : localBatchSize;
            FeatureExtractionThread* thread = new FeatureExtractionThread(this, localNumVoxels, outputDim, currentIndexOffset, threadCount * localBatchSize, border, 
                                                                            numFeatures, currentBatch, dataVolume, filterKernels, dctCoefficients, necessaryScales, scalingMinMax);
            workerThreads.push_back(thread);
            thread->run();
        }

        // join and delete all threads
        for (auto& t : workerThreads) {
            t->join();
            delete t;
        }
        
        performanceTimer.stop();
        avgFeatureExtractionTime += static_cast<double>(performanceTimer.getRuntime());
        performanceTimer.reset();
        performanceTimer.start();
        
        // classify the batch
        caffe::Blob<float>* input_layer = net->input_blobs()[0];
        float* input_data = input_layer->mutable_cpu_data();
        std::memcpy(&input_data[0], &currentBatch[0], inputBatchSize * numFeatures * sizeof(float));

        net->Forward(net->input_blobs());

        // retrieve result
        caffe::Blob<float>* output_layer = net->output_blobs()[0];
        const float* output_data = output_layer->cpu_data();
        
        performanceTimer.stop();
        avgClassificationTime += static_cast<double>(performanceTimer.getRuntime());
        performanceTimer.reset();
        
        // copy result into output volume and handle the case of a partially used last batch
        size_t resultsToCopy = inputBatchSize;
        if ((currentBatchNumber == numBatches - 1) && (numOutputVoxels % inputBatchSize != 0))
            resultsToCopy = numOutputVoxels % inputBatchSize;

        std::memcpy((void*) &outputBuffer[inputBatchSize * currentBatchNumber], &output_data[0], resultsToCopy * 2 * sizeof(float));

        // set the progress
        setProgress(std::min(0.99f, static_cast<float>(currentBatchNumber) / static_cast<float>(numBatches)));
        
        estimationTimer.stop();
        // update every fifth batch
        if (currentBatchNumber > 0 && currentBatchNumber % 5 == 0) {
            // time per batch in minutes
            double timePerBatch = static_cast<double>(estimationTimer.getRuntime()) / (5.0 * 1000 * 60);
            double estimate = static_cast<double>(numBatches - currentBatchNumber) * timePerBatch;
            std::stringstream strstr;
            if (tgt::iround(estimate) < 1)
                strstr << "< 1 min";
            else
                strstr << tgt::iround(estimate) << " min";
            timeEstimation_.set(strstr.str());
            estimationTimer.reset();
        }
    }

    deleteFilterKernels(filterKernels);
    delete[] dctCoefficients;

    delete[] currentBatch;

    classifierOutput_.setData(volume, true);

    setProgress(1.f);
    timeEstimation_.set("Finished");

    timer.stop();
    LINFO("Time for classification: " << (float) timer.getRuntime() / 1000.f << " seconds");
    
    LINFO("Average feature extraction time per batch (ms): " << avgFeatureExtractionTime / static_cast<double>(numBatches));
    LINFO("Average classification time per batch (ms): " << avgClassificationTime / static_cast<double>(numBatches));

    reactivateProperties();
}

void PatchCaffeClassifier::deactivateProperties() {
    PatchFeatureExtractor::deactivateProperties();
    protoFilename_.setReadOnlyFlag(true);
    modelFilename_.setReadOnlyFlag(true);
    computeButton_.setReadOnlyFlag(true);
    autoCompute_.setReadOnlyFlag(true);
}

void PatchCaffeClassifier::reactivateProperties() {
    PatchFeatureExtractor::reactivateProperties();
    protoFilename_.setReadOnlyFlag(false);
    modelFilename_.setReadOnlyFlag(false);
    computeButton_.setReadOnlyFlag(false);
    autoCompute_.setReadOnlyFlag(false);
}



} // namespace
