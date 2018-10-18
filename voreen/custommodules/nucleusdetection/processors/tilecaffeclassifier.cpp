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

#include "tilecaffeclassifier.h"

#include "voreen/core/datastructures/volume/volumebase.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
/*#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/histogram.h"*/

#include "tgt/stopwatch.h"

//#include <algorithm>

#include <caffe/caffe.hpp>

namespace voreen {

const std::string TileCaffeClassifier::loggerCat_("voreen.nucleusdetection.TileCaffeClassifier");

TileCaffeClassifier::TileCaffeClassifier()
    : VolumeProcessor()
    , dataInport_(Port::INPORT, "volumehandle.input", "Volume Input", false)
    , classifierOutput_(Port::OUTPORT, "classifier.output", "Classifier Output", false)
    , protoFilename_("protoFilename", "Prototxt File Name", "Select file...", "",
            "Prototxt (*.prototxt)", FileDialogProperty::OPEN_FILE, Processor::INVALID_PATH)
    , modelFilename_("modelFilename", "Model File Name", "Select file...", "",
            "Caffe Model (*.caffemodel)", FileDialogProperty::OPEN_FILE, Processor::INVALID_PATH)
    , tileSize_("tilesize", "2D Tile Size", 244, 4, 1000)
    , borderSize_("bordersize", "Border Size (on EACH side of tile)", 92, 0, 500) 
    , channelWidth_("channelwidth", "Channel Width (3rd dimension)", 5, 1, 25, 2)
    , useMirroring_("usemirroring", "Use Mirroring instead of zero padding", true)
    , computeButton_("classifyButton", "Classify Data")
    , timeEstimation_("timeestimation", "Runtime Estimation", " - ", Processor::VALID)
    , autoCompute_("autoCompute", "Auto Compute")
    , startComputation_(false)
{
    addPort(dataInport_);

    addPort(classifierOutput_);

    computeButton_.onChange(MemberFunctionCallback<TileCaffeClassifier>(this, &TileCaffeClassifier::classify));

    addProperty(protoFilename_);
    addProperty(modelFilename_);

    addProperty(tileSize_);
    addProperty(borderSize_);
    addProperty(channelWidth_);
    addProperty(useMirroring_);

    addProperty(computeButton_);
    
    timeEstimation_.setReadOnly(true);
    //timeEstimation_.setInstantUpdate(true);
    addProperty(timeEstimation_);

    addProperty(autoCompute_);
}

Processor* TileCaffeClassifier::create() const {
    return new TileCaffeClassifier();
}

void TileCaffeClassifier::process() {
    if (autoCompute_.get())
        classify();
}

void TileCaffeClassifier::classify() {

    if (!isInitialized())
        return;
    
    timeEstimation_.set(" - ");

    classifierOutput_.clear();

    if (!dataInport_.isReady()) {
        LERROR(" no volume ");
        return;
    }

    const VolumeRAM* dataVolume = dataInport_.getData()->getRepresentation<VolumeRAM>();

    if (!dataVolume) {
        LERROR("No VolumeRAM available!");
        return;
    }

    deactivateProperties();
    setProgress(0.f);

    // timer for estimating the remaining time
    tgt::Stopwatch estimationTimer;        
    estimationTimer.start();

    // check input volume
    if (dataInport_.getData()->getNumChannels() > 1) {
        LERROR("Cannot extract tiles because volume has multiple channels");
        reactivateProperties();
        estimationTimer.stop();
        return;
    }
    
    if (dataInport_.getData()->getDimensions().z < static_cast<size_t>(channelWidth_.get())) {
        LERROR("Input volume does not have enough z slices for selected channels");
        reactivateProperties();
        estimationTimer.stop();
        return;
    }

    // set GPU mode for Caffe (TODO: optionally allow for CPU mode)
    caffe::Caffe::set_mode(caffe::Caffe::GPU);
    int device_id = 0;
    caffe::Caffe::SetDevice(device_id);

    // load classifier
    std::shared_ptr<caffe::Net<float>> net;
    net.reset(new caffe::Net<float>(protoFilename_.get(), caffe::TEST));
    net->CopyTrainedLayersFrom(modelFilename_.get());
    
    caffe::Blob<float>* inputLayer = net->input_blobs()[0];
    size_t inputBatchSize = static_cast<size_t>(inputLayer->shape(0));
    size_t inputChannels = static_cast<size_t>(inputLayer->shape(1));
    size_t inputX = static_cast<size_t>(inputLayer->shape(2));
    size_t inputY = static_cast<size_t>(inputLayer->shape(3));

    // TODO: check if input sizes match the tile and channel settings
    // TODO: automatically set the tile settings from the loaded network!

    if (inputBatchSize != 1) {
        LERROR("Net does not have a batch size of 1 - aborting");
        reactivateProperties();
        estimationTimer.stop();
        return;
    }

    // get settings for extracting the tiles
    bool mirroring = useMirroring_.get();

    // get 2D size of a tile
    size_t tileSize = static_cast<size_t>(tileSize_.get());
    size_t borderSize = static_cast<size_t>(borderSize_.get());

    // get channel width (3rd dimension)
    size_t channelWidth = static_cast<size_t>(channelWidth_.get());
    if (channelWidth % 2 != 1) {
        LERROR("Unsupported channel width: only odd numbers of channels are supported!");
        reactivateProperties();
    }

    // compute the total number of tiles that have to be processed for progress information
    double numTilesTotal = 0.0;
    double numTilesProcessed = 0.0;

    tgt::svec3 volDim = dataVolume->getDimensions();
    // calculate num slices, subtract upper and lower border due to channel width
    size_t numSlices = volDim.z - (channelWidth - 1);
    // calculate numTiles for each slice
    size_t tilesX = (volDim.x / tileSize) + (volDim.x % tileSize > 0 ? 1 : 0);
    size_t tilesY = (volDim.y / tileSize) + (volDim.y % tileSize > 0 ? 1 : 0);

    numTilesTotal = static_cast<double>(numSlices * tilesX * tilesY);



    // since there may be a tile that does not entirely fit, we make it fit and distribute the overlap evenly among the data
    tgt::svec2 lastTileStart = tgt::svec2(volDim.x, volDim.y) - tgt::svec2(tileSize);
    tgt::svec2 offsets = (lastTileStart + tgt::svec2::one) / tgt::svec2(tilesX, tilesY) + tgt::svec2( ((lastTileStart.x+1) % tilesX > 0) ? 1 : 0,  ((lastTileStart.y+1) % tilesY > 0) ? 1 : 0);

    // create tile buffer
    float* tileBuffer = new float[channelWidth * (tileSize + 2* borderSize) * (tileSize + 2* borderSize)];

    // iterate over slices
    for (size_t z = channelWidth / 2; z < volDim.z - channelWidth / 2; ++z) {

        //iterate over tiles in x and y direction
        for (size_t y = 0; y < tilesY; ++y) {
            for (size_t x = 0; x < tilesX; ++x) {

                tgt::svec3 tileStartVoxel = tgt::svec3(std::min(x * offsets.x, lastTileStart.x), std::min(y * offsets.y, lastTileStart.y), z - channelWidth / 2);

                // we need to take into account the border (positions may be negative)
                tgt::ivec3 tileWithBorderStartVoxel = tgt::ivec3(static_cast<int>(tileStartVoxel.x) - static_cast<int>(borderSize),
                                                                 static_cast<int>(tileStartVoxel.y) - static_cast<int>(borderSize),
                                                                 static_cast<int>(tileStartVoxel.z));


                // now extract the tile
                //size_t currentPosition = 0;
                // TODO: include mirroring in Z-direction
                for (size_t currentZ = 0; currentZ < channelWidth; ++currentZ) {
                    for (size_t currentY = 0; currentY < tileSize + 2* borderSize; ++currentY) {
                        for (size_t currentX = 0; currentX < tileSize + 2* borderSize; ++currentX) {

                            tgt::ivec3 voxelCoordinate = tgt::ivec3(currentX, currentY, currentZ) + tileWithBorderStartVoxel;

                            float currentIntensity = 0.f;
                            // check if coordinate is inside volume
                            if (tgt::hand(tgt::greaterThanEqual(voxelCoordinate, tgt::ivec3::zero)) && tgt::hand(tgt::lessThan(voxelCoordinate, tgt::ivec3(volDim))))
                                currentIntensity = dataVolume->getVoxelNormalized(tgt::svec3(voxelCoordinate)); 
                            else if (mirroring) {
                                // use mirroring instead of padding with zeros
                                if (voxelCoordinate.x < 0)
                                    voxelCoordinate.x = std::abs(voxelCoordinate.x) % static_cast<int>(volDim.x);

                                if (voxelCoordinate.y < 0)
                                    voxelCoordinate.y = std::abs(voxelCoordinate.y) % static_cast<int>(volDim.y);

                                if (voxelCoordinate.x >= static_cast<int>(volDim.x)) {
                                    voxelCoordinate.x = (static_cast<int>(volDim.x) - 1) - (voxelCoordinate.x % (static_cast<int>(volDim.x) - 1));
                                }

                                if (voxelCoordinate.y >= static_cast<int>(volDim.y)) {
                                    voxelCoordinate.y = (static_cast<int>(volDim.y) - 1) - (voxelCoordinate.y % (static_cast<int>(volDim.y) - 1));
                                }

                                tgtAssert((tgt::hand(tgt::greaterThanEqual(voxelCoordinate, tgt::ivec3::zero)) && tgt::hand(tgt::lessThan(voxelCoordinate, tgt::ivec3(volDim)))), "voxel outside volume");
                                
                                currentIntensity = dataVolume->getVoxelNormalized(tgt::svec3(voxelCoordinate));
                            }

                            // set current intensity to a tile buffer
                            size_t index = currentZ * (tileSize + 2* borderSize) * (tileSize + 2* borderSize) + currentY * (tileSize + 2* borderSize) + currentX;
                            tileBuffer[index] = currentIntensity;
                        }
                    }
                }              

                // classify the tile buffer and write it to the output volume
                caffe::Blob<float>* input_layer = net->input_blobs()[0];
                float* input_data = input_layer->mutable_cpu_data();
                std::memcpy(&input_data[0], &tileBuffer[0], channelWidth * (tileSize + 2* borderSize) * (tileSize + 2* borderSize) * sizeof(float));

                net->Forward(net->input_blobs());

                // retrieve result
                caffe::Blob<float>* output_layer = net->output_blobs()[0];
                const float* output_data = output_layer->cpu_data();           
                
                // TODO: set classified tile to output volume


                      
                numTilesProcessed += 1.0;

            }
        }

        // we have processed another slice -> we compute the new time estimation and progress information
        double processedPart = numTilesProcessed / numTilesTotal;
        double currentTime = estimationTimer.getRuntime() / (1000.0 * 60.0); // runtime in minutes
        double estimate = currentTime * ((numTilesTotal-numTilesProcessed) / numTilesProcessed);
        std::stringstream strstr;
        if (tgt::iround(estimate) < 1)
            strstr << "< 1 min";
        else
            strstr << tgt::iround(estimate) << " min";
        timeEstimation_.set(strstr.str());
        setProgress(std::min(0.99f, static_cast<float>(processedPart)));
    }

    // delete tile buffer
    delete[] tileBuffer;

    // TODO: set output volume
    //classifierOutput_.setData(volume, true);

    estimationTimer.stop();
    LINFO("Time for classification: " << (float) estimationTimer.getRuntime() / 1000.f << " seconds");

    setProgress(1.f);
    timeEstimation_.set("Finished");

    reactivateProperties();
}

void TileCaffeClassifier::deactivateProperties() {
    protoFilename_.setReadOnlyFlag(true);
    modelFilename_.setReadOnlyFlag(true);
    tileSize_.setReadOnlyFlag(true);
    borderSize_.setReadOnlyFlag(true);
    channelWidth_.setReadOnlyFlag(true);
    useMirroring_.setReadOnlyFlag(true);
    computeButton_.setReadOnlyFlag(true);
    autoCompute_.setReadOnlyFlag(true);
}

void TileCaffeClassifier::reactivateProperties() {
    protoFilename_.setReadOnlyFlag(false);
    modelFilename_.setReadOnlyFlag(false);
    tileSize_.setReadOnlyFlag(false);
    borderSize_.setReadOnlyFlag(false);
    channelWidth_.setReadOnlyFlag(false);
    useMirroring_.setReadOnlyFlag(false);
    computeButton_.setReadOnlyFlag(false);
    autoCompute_.setReadOnlyFlag(false);
}



} // namespace
