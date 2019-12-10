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

#include "tileextractor.h"


#include "voreen/core/datastructures/volume/volumebase.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/histogram.h"
#include "modules/hdf5/utils/hdf5utils.h"

#include "tgt/stopwatch.h"

namespace voreen {

const std::string TileExtractor::loggerCat_("voreen.nucleusdetection.TileExtractor");

TileExtractor::TileExtractor()
    : VolumeProcessor()
    , dataInport_(Port::INPORT, "volumedata.input", "Volume Data Input", false)
    , labelInport_(Port::INPORT, "volumelabels.input", "Label Volume Input", false)
    , intensityReferenceVolume_(Port::INPORT, "intensityreference.input", "Intensity Reference Volume Input", false)
    , tileSize_("tilesize", "2D Tile Size", 244, 20, 1000)
    , borderSize_("bordersize", "Border Size (on EACH side of tile)", 92, 0, 500) 
    , channelWidth_("channelwidth", "Channel Width (3rd dimension)", 5, 1, 25, 2)
    , filename_("outputFilename", "Training Data Output File", "Select file...", "",
            "HDF5 (*.h5)", FileDialogProperty::SAVE_FILE, Processor::INVALID_PATH)
    , numTilesPerFile_("numtiles", "Number of Tiles per HDF5-File", 50, 1, 500000, Processor::INVALID_RESULT, NumericProperty<int>::DYNAMIC)
    , useMirroring_("usemirroring", "Use Mirroring instead of zero padding", true)
    , computeButton_("train", "Extract Training Data")
    , autoCompute_("autoCompute", "Auto Compute")
    , timeEstimation_("timeestimation", "Runtime Estimation", " - ", Processor::VALID)
    , startComputation_(false)
    , progressProperty_("progressProperty", "Progress")
{
    addPort(dataInport_);
    addPort(labelInport_);
    addPort(intensityReferenceVolume_);

    addProperty(tileSize_);
    addProperty(borderSize_);
    addProperty(channelWidth_);

    tileSize_.onChange(MemberFunctionCallback<TileExtractor>(this, &TileExtractor::adjustTilesPerFile));
    borderSize_.onChange(MemberFunctionCallback<TileExtractor>(this, &TileExtractor::adjustTilesPerFile));
    channelWidth_.onChange(MemberFunctionCallback<TileExtractor>(this, &TileExtractor::adjustTilesPerFile));

    computeButton_.onChange(MemberFunctionCallback<TileExtractor>(this, &TileExtractor::computeTiles));

    timeEstimation_.setEditable(false);
    //timeEstimation_.setInstantUpdate(true);
    addProperty(timeEstimation_);

    addProperty(filename_);
    addProperty(numTilesPerFile_);

    addProperty(useMirroring_);

    addProperty(computeButton_);
    addProperty(autoCompute_);

    // progress information
    addProperty(progressProperty_);
    addProgressBar(&progressProperty_);

    adjustTilesPerFile();
}

Processor* TileExtractor::create() const {
    return new TileExtractor();
}

void TileExtractor::process() {
    if (autoCompute_.get())
        computeTiles();
}

bool TileExtractor::isReady() const {
    if (!isInitialized())
        return false;

    return true;
}

void TileExtractor::adjustTilesPerFile() {
    // get 2D size of a tile
    size_t tileSize = static_cast<size_t>(tileSize_.get());
    size_t borderSize = static_cast<size_t>(borderSize_.get());

    // get channel width (3rd dimension)
    size_t channelWidth = static_cast<size_t>(channelWidth_.get());

    // compute memory for a single tile (including labels)
    size_t memPerTile = sizeof(float) * channelWidth * (tileSize + 2 * borderSize) * (tileSize + 2 * borderSize) + sizeof(float) * tileSize * tileSize;

    // divide 2 GB by memory of a tile to find the maximum number
    size_t caffeLimit = 2 * 1000 * 1000 * 1000;

    size_t maxTiles = caffeLimit / memPerTile;

    if (maxTiles == 0) {
        LERROR("Memory size of a single tile exceeds Caffe limit of 2GB");
        maxTiles = 1;
    }

    numTilesPerFile_.setMaxValue(static_cast<int>(std::min(maxTiles, static_cast<size_t>(500000))));
}

void TileExtractor::computeTiles() {

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

    if (!(dataInport_.getData()->size() == labelInport_.getData()->size())) {
        LERROR(" Different number of data and label volumes! ");
        return;
    }

    // start
    deactivateProperties();
    setProgress(0.f);

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

    // compute the total number of tiles that have to be processed over all volumes for progress information
    double numTilesTotal = 0.0;
    double numTilesProcessed = 0.0;
    for (size_t volumeIndex = 0; volumeIndex < dataInport_.getData()->size(); ++volumeIndex) {
        // get dimensions of this volume
        tgt::svec3 volDim = dataInport_.getData()->at(volumeIndex)->getDimensions();
        // calculate num slices, subtract upper and lower border due to channel width
        size_t numSlices = 0;
        if (volDim.z > channelWidth - 1)
            numSlices = volDim.z - (channelWidth - 1);
        // calculate numTiles for each slice
        size_t tilesX = (volDim.x / tileSize) + (volDim.x % tileSize > 0 ? 1 : 0);
        size_t tilesY = (volDim.y / tileSize) + (volDim.y % tileSize > 0 ? 1 : 0);

        numTilesTotal += static_cast<double>(numSlices * tilesX * tilesY);
    }

    // timer for estimating the remaining time
    tgt::Stopwatch estimationTimer;        
    estimationTimer.start();

    size_t numForegroundVoxels = 0;
    size_t numBackgroundVoxels = 0;

    for (size_t volumeIndex = 0; volumeIndex < dataInport_.getData()->size(); ++volumeIndex) {

        tgt::svec3 volDim = dataInport_.getData()->at(volumeIndex)->getDimensions();
        const VolumeRAM* dataVolume = dataInport_.getData()->at(volumeIndex)->getRepresentation<VolumeRAM>();
        const VolumeRAM* labelVolume = labelInport_.getData()->at(volumeIndex)->getRepresentation<VolumeRAM>();

        if (!(dataInport_.getData()->at(volumeIndex)->getDimensions() == labelInport_.getData()->at(volumeIndex)->getDimensions())) {
            LERROR(" Dimensions of label volume and input volume " << volumeIndex << " are different! ");
            reactivateProperties();
            estimationTimer.stop();
            return;
        }

        if (!dataVolume || !labelVolume) {
            LERROR("No VolumeRAM available!");
            reactivateProperties();
            estimationTimer.stop();
            return;
        }

        if (dataInport_.getData()->at(volumeIndex)->getNumChannels() > 1) {
            LERROR("Cannot extract tiles for volume " << volumeIndex << " because volume has multiple channels");
            continue;
        }

        // usually data is normalized over uint16_t range, but we want to normalize over the volume's histogram to account for 
        // variation of fluorescence intensity across different data sets -> use histogram information to scale values
        tgt::vec2 scalingMinMax = tgt::vec2(0.f, 1.f);
        // if present, we use the intensity reference volume
        if (intensityReferenceVolume_.getData())
            scalingMinMax = getScalingMinMaxIntensity(intensityReferenceVolume_.getData());
        else
            scalingMinMax = getScalingMinMaxIntensity(dataInport_.getData()->at(volumeIndex));


        // number of tiles per file
        size_t numTiles = static_cast<size_t>(numTilesPerFile_.get());

        // determine number of files for this volume
        if (volDim.z <= channelWidth) {
            LERROR("Cannot extract tiles for volume " << volumeIndex << " because not enough slices are present for this channel width");
            continue;
        }
        
        if (volDim.x < tileSize) {
            LERROR("Cannot extract tiles for volume " << volumeIndex << " because x dimension is smaller than tile size");
            continue;
        }
        
        if (volDim.y < tileSize) {
            LERROR("Cannot extract tiles for volume " << volumeIndex << " because y dimension is smaller than tile size");
            continue;
        }
        
        size_t localSlices = volDim.z - (channelWidth - 1);
        size_t localTilesX = (volDim.x / tileSize) + (volDim.x % tileSize > 0 ? 1 : 0);
        size_t localTilesY = (volDim.y / tileSize) + (volDim.y % tileSize > 0 ? 1 : 0);

        // since there may be a tile that does not entirely fit, we want to make it fit and distribute the overlap evenly among the data
        tgt::svec2 lastTileStart = tgt::svec2(volDim.x, volDim.y) - tgt::svec2(tileSize);
        tgt::svec2 offsets = (lastTileStart + tgt::svec2::one) / tgt::svec2(localTilesX, localTilesY) + tgt::svec2( ((lastTileStart.x+1) % localTilesX > 0) ? 1 : 0,  ((lastTileStart.y+1) % localTilesY > 0) ? 1 : 0);

        // determine base file name
        std::string baseFileName = tgt::FileSystem::fullBaseName(tgt::FileSystem::cleanupPath(filename_.get()));

        // create a data matrix and label vector which will be filled for each file   
        float* featureData = new float[numTiles * channelWidth * (tileSize + 2 * borderSize) * (tileSize + 2 * borderSize)];
        float* labelData = new float[numTiles * tileSize * tileSize];

        size_t tilePosition = 0;
        size_t localFileNumber = 0;

        // iterate over slices
        for (size_t z = channelWidth / 2; z < volDim.z - channelWidth / 2; ++z) {
        
            //iterate over tiles in x and y direction
            for (size_t y = 0; y < localTilesY; ++y) {
                for (size_t x = 0; x < localTilesX; ++x) {

                    tgt::svec3 tileStartVoxel = tgt::svec3(std::min(x * offsets.x, lastTileStart.x), std::min(y * offsets.y, lastTileStart.y), z - channelWidth / 2);

                    // we need to take into account the border (positions may be negative)
                    tgt::ivec3 tileWithBorderStartVoxel = tgt::ivec3(static_cast<int>(tileStartVoxel.x) - static_cast<int>(borderSize),
                                                                     static_cast<int>(tileStartVoxel.y) - static_cast<int>(borderSize),
                                                                     static_cast<int>(tileStartVoxel.z));

                    // now extract the tile
                    size_t currentPosition = 0;
                    for (size_t currentZ = 0; currentZ < channelWidth; ++currentZ) {
                        for (size_t currentY = 0; currentY < tileSize + 2* borderSize; ++currentY) {
                            for (size_t currentX = 0; currentX < tileSize + 2* borderSize; ++currentX) {
                            
                                tgt::ivec3 voxelCoordinate = tgt::ivec3(currentX, currentY, currentZ) + tileWithBorderStartVoxel;

                                float currentIntensity = 0.f;
                                // check if coordinate is inside volume
                                if (tgt::hand(tgt::greaterThanEqual(voxelCoordinate, tgt::ivec3::zero)) && tgt::hand(tgt::lessThan(voxelCoordinate, tgt::ivec3(volDim))))
                                    currentIntensity = scaleIntensity(dataVolume->getVoxelNormalized(tgt::svec3(voxelCoordinate)), scalingMinMax); 
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
                                    
                                    currentIntensity = scaleIntensity(dataVolume->getVoxelNormalized(tgt::svec3(voxelCoordinate)), scalingMinMax);
                                }

                                size_t index = tilePosition * channelWidth * (tileSize + 2 * borderSize) * (tileSize + 2 * borderSize) + currentPosition;
                                tgtAssert(index < numTiles * channelWidth * (tileSize + 2 * borderSize) * (tileSize + 2 * borderSize), "Index out of bounds");
                                featureData[index] = currentIntensity;

                                currentPosition++;
                            }
                        }
                    }

                    currentPosition = 0;
                    // now extract the labels (without the border and only the middle slice we want to classify)
                    for (size_t currentY = 0; currentY < tileSize; ++currentY) {
                        for (size_t currentX = 0; currentX < tileSize; ++currentX) {
                        
                            tgt::svec3 voxelCoordinate = tgt::svec3(currentX + tileStartVoxel.x, currentY + tileStartVoxel.y, z); 

                            float currentLabel = 0.f;
                            // check if coordinate is inside volume
                            if (tgt::hand(tgt::greaterThanEqual(voxelCoordinate, tgt::svec3::zero)) && tgt::hand(tgt::lessThan(voxelCoordinate, volDim))) {
                                // inside volume -> read label
                                if (labelVolume->getVoxelNormalized(voxelCoordinate) != 0) {
                                    currentLabel = 1.f;
                                    numForegroundVoxels++;
                                }
                                else
                                    numBackgroundVoxels++; 
                            }
                            else if (mirroring) {
                                // use mirroring instead of padding with zeros
                                // negative values cannot occur
                                //if (voxelCoordinate.x < 0)
                                //    voxelCoordinate.x = std::abs(voxelCoordinate.x) % static_cast<int>(volDim.x);

                                //if (voxelCoordinate.y < 0)
                                //    voxelCoordinate.y = std::abs(voxelCoordinate.y) % static_cast<int>(volDim.y);

                                if (voxelCoordinate.x >= volDim.x) {
                                    voxelCoordinate.x = (volDim.x - 1) - (voxelCoordinate.x % (volDim.x - 1));
                                }

                                if (voxelCoordinate.y >= volDim.y) {
                                    voxelCoordinate.y = (volDim.y - 1) - (voxelCoordinate.y % (volDim.y - 1));
                                }

                                tgtAssert((tgt::hand(tgt::greaterThanEqual(voxelCoordinate, tgt::svec3::zero)) && tgt::hand(tgt::lessThan(voxelCoordinate, volDim))), "voxel outside label volume");
                                
                                if (labelVolume->getVoxelNormalized(voxelCoordinate) != 0) {
                                    currentLabel = 1.f;
                                    numForegroundVoxels++;
                                }
                                else
                                    numBackgroundVoxels++;
                            }

                            size_t index = tilePosition * tileSize * tileSize + currentPosition;
                            tgtAssert(index < numTiles * tileSize * tileSize, "Index out of bounds");
                            labelData[index] = currentLabel;

                            currentPosition++;
                        }
                    }


                    tilePosition++;
                    numTilesProcessed += 1.0;

                    if (tilePosition == numTiles) {
                        // a file is complete -> write it out
                        std::stringstream fileNameStream;
                        fileNameStream << baseFileName << "_" << volumeIndex << "_" << localFileNumber << ".h5";

                        using namespace H5;
                        const H5std_string  FILE_NAME( fileNameStream.str() );
                        const H5std_string  DATA_NAME( "data" );
                        const H5std_string  OUTPUT_NAME( "label" );

                        H5File file(FILE_NAME, H5F_ACC_TRUNC);
                        hsize_t dataDimsf[4];
                        dataDimsf[0] = numTiles;
                        dataDimsf[1] = channelWidth;
                        dataDimsf[2] = (tileSize + 2 * borderSize);
                        dataDimsf[3] = (tileSize + 2 * borderSize);
                        DataSpace dataSpace(4,dataDimsf);
                        hsize_t labelDimsf[4];
                        labelDimsf[0] = tilePosition;
                        labelDimsf[1] = 1;
                        labelDimsf[2] = tileSize;
                        labelDimsf[3] = tileSize;
                        DataSpace outputSpace(4,labelDimsf);
                        FloatType dataType(PredType::NATIVE_FLOAT);
                        DataSet featureSet = file.createDataSet( DATA_NAME, dataType, dataSpace);
                        DataSet labelSet = file.createDataSet(OUTPUT_NAME, dataType, outputSpace);
                        featureSet.write(featureData, PredType::NATIVE_FLOAT);
                        labelSet.write(labelData,PredType::NATIVE_FLOAT);
                        file.close();

                        tilePosition = 0;
                        localFileNumber++;
                    }
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

        // write out the rest of the data in a last file
        if (tilePosition > 0) {
            std::stringstream fileNameStream;
            fileNameStream << baseFileName << "_" << volumeIndex << "_" << localFileNumber << ".h5";

            using namespace H5;
            const H5std_string  FILE_NAME( fileNameStream.str() );
            const H5std_string  DATA_NAME( "data" );
            const H5std_string  OUTPUT_NAME( "label" );

            H5File file(FILE_NAME, H5F_ACC_TRUNC);
            hsize_t dataDimsf[4];
            dataDimsf[0] = tilePosition;
            dataDimsf[1] = channelWidth;
            dataDimsf[2] = (tileSize + 2 * borderSize);
            dataDimsf[3] = (tileSize + 2 * borderSize);
            DataSpace dataSpace(4,dataDimsf);
            hsize_t labelDimsf[4];
            labelDimsf[0] = tilePosition;
            labelDimsf[1] = 1;
            labelDimsf[2] = tileSize;
            labelDimsf[3] = tileSize;
            DataSpace outputSpace(4,labelDimsf);
            FloatType dataType(PredType::NATIVE_FLOAT);
            DataSet featureSet = file.createDataSet( DATA_NAME, dataType, dataSpace);
            DataSet labelSet = file.createDataSet(OUTPUT_NAME, dataType, outputSpace);
            featureSet.write(featureData, PredType::NATIVE_FLOAT);
            labelSet.write(labelData,PredType::NATIVE_FLOAT);
            file.close();

            tilePosition = 0;
            localFileNumber++;
        }

        delete[] featureData;
        delete[] labelData;

    }

    estimationTimer.stop();
    LINFO("Time for generating and exporting training data: " << (float) estimationTimer.getRuntime() / 1000.f << " seconds");

    setProgress(1.f);
    timeEstimation_.set("Finished");

    LINFO("Wrote out training data containing " << numForegroundVoxels << " foreground voxels and " << numBackgroundVoxels << " background voxels");

    reactivateProperties();
}

tgt::vec2 TileExtractor::getScalingMinMaxIntensity(const VolumeBase* volume) const {
    // TODO: channel selection would be nice...
    float minIntensity = volume->getDerivedData<VolumeMinMax>()->getMinNormalized();
    float maxTmp = volume->getDerivedData<VolumeMinMax>()->getMaxNormalized();
    // for the maximum, we only want to use 99% of the data to ignore outliers
    VolumeHistogramIntensity* histogram = volume->getDerivedData<VolumeHistogramIntensity>();
    float sumBuckets = 0;
    for (size_t i = 0; i < histogram->getBucketCount(); ++i) {
        sumBuckets += histogram->getNormalized(static_cast<int>(i));
    }
    int i = -1;
    float currentSum = 0.f;
    while (((currentSum / sumBuckets) < 0.99f) && i < (int) histogram->getBucketCount()) {
        i += 1;
        currentSum += histogram->getNormalized(i);
    }
    float maxIntensity = (static_cast<float>(i) / static_cast<float>(histogram->getBucketCount() - 1)) * (maxTmp - minIntensity) + minIntensity;
    tgtAssert(maxIntensity - minIntensity > 0, "max = min intensity");

    return tgt::vec2(minIntensity, maxIntensity);
}

float TileExtractor::scaleIntensity(float value, tgt::vec2 minMaxValue) const {
    // scale, do not clamp to 1
    return (value - minMaxValue.x) / (minMaxValue.y - minMaxValue.x);
}

void TileExtractor::deactivateProperties() {
    tileSize_.setReadOnlyFlag(true);
    borderSize_.setReadOnlyFlag(true);
    channelWidth_.setReadOnlyFlag(true);
    //normalizeTiles_.setReadOnlyFlag(true);
    filename_.setReadOnlyFlag(true);
    numTilesPerFile_.setReadOnlyFlag(true);
    useMirroring_.setReadOnlyFlag(true);
    computeButton_.setReadOnlyFlag(true);
    autoCompute_.setReadOnlyFlag(true);
}

void TileExtractor::reactivateProperties() {
    tileSize_.setReadOnlyFlag(false);
    borderSize_.setReadOnlyFlag(false);
    channelWidth_.setReadOnlyFlag(false);    
    // normalizeTiles_.setReadOnlyFlag(false);
    filename_.setReadOnlyFlag(false);
    numTilesPerFile_.setReadOnlyFlag(false);
    useMirroring_.setReadOnlyFlag(false);
    computeButton_.setReadOnlyFlag(false);
    autoCompute_.setReadOnlyFlag(false);
}

} // namespace
