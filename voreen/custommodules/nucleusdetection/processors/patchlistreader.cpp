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

#include "patchlistreader.h"

#include "voreen/core/datastructures/volume/volumebase.h"
#include "voreen/core/datastructures/volume/volumeram.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "tgt/filesystem.h"

#include "modules/hdf5/utils/hdf5utils.h"

namespace voreen {

const std::string PatchListReader::loggerCat_("voreen.nucleusdetection.PatchListReader");

PatchListReader::PatchListReader()
    : VolumeProcessor()
    , outport_(Port::OUTPORT, "volumelist.output", "Volume List Outport")
    , filename_("inputFilename", "File", "Select file...", "",
            "HDF5 file (*.h5 *.hdf5)", FileDialogProperty::OPEN_FILE, Processor::INVALID_PATH)
    , normalizePatches_("normalize", "Normalize Patches", true)
    , allowAdditionalFeaturesAtEnd_("allowAtEnd", "Allow additional two trailing features", true) 
    , loadButton_("load", "Load")
    , autoCompute_("autoLoad", "Auto Load", false)
    , outputList_(0)
{
    addPort(outport_);

    loadButton_.onChange(MemberFunctionCallback<PatchListReader>(this, &PatchListReader::loadPatches));

    addProperty(filename_);
    addProperty(normalizePatches_);
    addProperty(allowAdditionalFeaturesAtEnd_);
    addProperty(loadButton_);
    
    addProperty(autoCompute_);
}

Processor* PatchListReader::create() const {
    return new PatchListReader();
}

void PatchListReader::deinitialize() {
    clearOutputList();
    VolumeProcessor::deinitialize();
}

void PatchListReader::clearOutputList() {
    if (!outputList_)
        return;

    if (outport_.getData())
        outport_.setData(0);

    // delete all volumes and delete list
    while (!outputList_->empty()) {
        VolumeBase* v = outputList_->first();
        outputList_->remove(v);
        delete v;
    }
    
    delete outputList_;
    outputList_ = 0;
}

void PatchListReader::process() {
    //clearOutputList();
    if (autoCompute_.get())
        loadPatches();
}

void PatchListReader::loadPatches() {
    if (!isInitialized())
        return;

    if (filename_.get().empty()) {
        LERROR("no filename");
        return;
    }

    if (outputList_)
        clearOutputList();

    std::string filename = tgt::FileSystem::cleanupPath(filename_.get());
    std::string datasetName = "filters";
    
    H5::H5File file;
    try {
        file = H5::H5File(filename, H5F_ACC_RDONLY );
    }
    catch (H5::FileIException&) {
        LERROR("Could not open HDF5 input file - not reading patch list");
        return;
    }
    
    size_t objsInFile = file.getNumObjs();
    std::vector<std::string> datasetNames;

    for (size_t item=0; item < objsInFile; item++){
        datasetNames.push_back(file.getObjnameByIdx(item));
    }

    std::vector<std::string>::iterator dataSetIterator = find(datasetNames.begin(),datasetNames.end(), datasetName);

    if (dataSetIterator == datasetNames.end()) {
        LERROR("Expected dataset \"filters\"... abort");
        file.close();
        return;
    }

    H5::DataSet dataSet = file.openDataSet(*dataSetIterator);
    H5::DataSpace dataSpace = dataSet.getSpace();
    hsize_t dims[2]; dims[0] = 0; dims[1] = 0;
    size_t rank = dataSpace.getSimpleExtentNdims();
    float* patches = 0;
    if (rank==2) {
        dataSpace.getSimpleExtentDims(dims);
        
        hsize_t start[2];
        start[0] = 0;
        start[1] = 0;

        H5::DataSpace memspace(2, dims);
        patches = new float[dims[0]*dims[1]];
        dataSpace.selectHyperslab(H5S_SELECT_SET, dims, start);
        dataSet.read(patches,H5::PredType::NATIVE_FLOAT, memspace, dataSpace);
    }
    else {
        LERROR("Wrong dimensions... expected two-dimensonal data");
        file.close();
        return;
    }

    file.close();
    
    size_t numFeatures = dims[1];
    size_t edgeLength = std::cbrt(numFeatures);

    size_t numFeaturesResultingFromEdgeLength = std::pow(edgeLength, 3);
    size_t featureDifference = numFeatures - numFeaturesResultingFromEdgeLength;
    
    if (!allowAdditionalFeaturesAtEnd_.get() && numFeaturesResultingFromEdgeLength != numFeatures) {
        LERROR("Expected uniform edge length in patch... feature size is " << dims[1] << "... aborting ");
        return;
    }
    
    if (featureDifference != 0 && featureDifference != 2) {
        LERROR("Detected additional number of trailing features - expected " << numFeaturesResultingFromEdgeLength << " or " << numFeaturesResultingFromEdgeLength + 2 << ", found " << numFeatures);
        return;
    }
    

    if (dims[0] < 1) {
        LERROR("Found no filters... aborting");
        return;
    }

    size_t numPatches = dims[0];

    outputList_ = new VolumeList();

    for (size_t i = 0; i < numPatches; ++i) {
        VolumeRAM_Float* volumeRAM = new VolumeRAM_Float(tgt::svec3(edgeLength), true);

        if (normalizePatches_.get()) {
            size_t baseIndex = i * numFeatures;
            // get mean and standard deviation to normalize the patch
            float mean = 0.f; float sigma = 0.f;
            float sum = 0.f;
            for (size_t index = 0; index < numFeatures - featureDifference; ++index) {
                size_t inputIndex = baseIndex + index; 
                sum += patches[inputIndex];
            }
            mean = sum / (static_cast<float>(numFeatures - featureDifference));

            for (size_t index = 0; index < numFeatures - featureDifference; ++index) {
                size_t inputIndex = baseIndex + index; 
                sigma += std::pow(patches[inputIndex] - mean, 2.f);
            }
            sigma *= 1.f / (static_cast<float>(numFeatures - featureDifference) - 1.f);
            sigma = std::sqrt(sigma);
            
            // if sigma is near 0: set to 1
            if (sigma < 1e-6f)
                sigma = 1.f;

            // use mean and standard deviation to normalize the patch
            for (size_t index = 0; index < numFeatures - featureDifference; ++index) {
                size_t inputIndex = baseIndex + index; 
                float value = (patches[inputIndex] - mean) / sigma;
                volumeRAM->voxel(index) = value;
            }
        }
        else {
            // no normalization -> just copy the data
            for (size_t index = 0; index < numFeatures - featureDifference; ++index) {
                size_t inputIndex = i * numFeatures + index; 
                volumeRAM->voxel(index) = patches[inputIndex];
            }
        }
        Volume* outputVolume = new Volume(volumeRAM, tgt::vec3(1.f), tgt::vec3(0.f));
        outputList_->add(outputVolume);
    }

    delete[] patches;

    outport_.setData(outputList_, false);
}

} // namespace
