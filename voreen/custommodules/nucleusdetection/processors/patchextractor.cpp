/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "patchextractor.h"
#include "voreen/core/io/progressbar.h"

#include "voreen/core/datastructures/volume/volumebase.h"
#include "voreen/core/datastructures/volume/volumeram.h"

#include "tgt/filesystem.h"
#include "modules/hdf5/utils/hdf5utils.h"
namespace voreen {

const std::string PatchExtractor::loggerCat_("voreen.nucleusdetection.PatchExtractor");

PatchExtractor::PatchExtractor()
    : VolumeProcessor()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input", false)
    , filename_("outputFilename", "File", "Select file...", "",
            "HDF5 (*.h5)", FileDialogProperty::SAVE_FILE, Processor::INVALID_PATH)
    , normalizePatches_("normalize", "Normalize Patches", true)
    , computeButton_("save", "Save")
    , startComputation_(false)
    , progressProperty_("progressProperty", "Progress")
{
    addPort(inport_);

    computeButton_.onChange(MemberFunctionCallback<PatchExtractor>(this, &PatchExtractor::computePatches));

    addProperty(filename_);
    addProperty(normalizePatches_);
    addProperty(computeButton_);
    addProperty(progressProperty_);
    addProgressBar(&progressProperty_);
}

Processor* PatchExtractor::create() const {
    return new PatchExtractor();
}

void PatchExtractor::process() {
    // do nothing
}

void PatchExtractor::computePatches() {
    if (!isInitialized())
        return;

    if (!inport_.isReady() || !inport_.getData() || inport_.getData()->empty()) {
        LERROR(" no volume data ");
        return;
    }

    if (filename_.get().empty()) {
        LERROR("no filename");
        return;
    }
    
    // check volumes and compute number of patches
    size_t numPatchesTotal = 0;
    for (size_t i = 0; i < inport_.getData()->size(); ++i) {
        const VolumeBase* v = inport_.getData()->at(i);
        tgt::svec3 volDim = v->getDimensions();
        
        if (volDim.x <= 6 || volDim.y <= 6 || volDim.z <= 6) {
            LERROR("Volume " << i << " is too small!");
            return;
        }
        
        numPatchesTotal += (volDim.x - 6) * (volDim.y - 6) * (volDim.z - 6);
        
        const VolumeRAM* volume = v->getRepresentation<VolumeRAM>();
        
        if (!volume) {
            LERROR("No VolumeRAM for volume " << i << " available!");
            return;
        }
    }

    //const VolumeRAM* volume = inport_.getData()->getRepresentation<VolumeRAM>();
    //tgt::svec3 volDim = inport_.getData()->getDimensions();
    
    //if (volDim.x <= 6 || volDim.y <= 6 || volDim.z <= 6) {
    //    LERROR("Volume is too small!");
    //    return;
    //}
    
    
    // deactivate settings
    filename_.setReadOnlyFlag(true);
    normalizePatches_.setReadOnlyFlag(true);
    computeButton_.setReadOnlyFlag(true);
    
    //size_t numPatchesTotal = (volDim.x - 6) * (volDim.y - 6) * (volDim.z - 6); // number of patches that will be written in total
    size_t numPatches = 0;      // number of patches that has been written
    float patch[7 * 7 * 7];     // a single patch

    float* featureData = new float[7 * 7 * 7 * numPatchesTotal];

    // now create an output file
    std::string filename = tgt::FileSystem::cleanupPath(filename_.get());
    
    using namespace H5;
    const H5std_string FILE_NAME(filename);
    const H5std_string DATA_NAME( "data" );
    const int featureLength = 7 * 7 * 7;
    const int sampleSize = numPatchesTotal;
        
    H5File file(FILE_NAME, H5F_ACC_TRUNC);
    hsize_t dimsf[2];
    dimsf[0] = sampleSize;
    dimsf[1] = featureLength;
    DataSpace dataSpace(2,dimsf);
    FloatType dataType(PredType::NATIVE_FLOAT);
    DataSet featureSet = file.createDataSet( DATA_NAME, dataType, dataSpace);

    for (size_t i = 0; i < inport_.getData()->size(); ++i) {
        const VolumeBase* v = inport_.getData()->at(i);
        tgt::svec3 volDim = v->getDimensions();
        const VolumeRAM* volume = v->getRepresentation<VolumeRAM>();
        
        if (!volume) {
            LERROR("No VolumeRAM for volume " << i << " available!");
            return;
        }
        
        // extract patches for this volume

        for (size_t z = 3; z + 3 < volDim.z; ++z) {
                    
            setProgress(std::min(static_cast<float>(numPatches) / static_cast<float>(numPatchesTotal), 0.99f));
            
            for (size_t y = 3; y + 3 < volDim.y; ++y) {
                for (size_t x = 3; x + 3 < volDim.x; ++x) {

                    size_t patchIndex = 0;
                    // (x,y,z) is our current voxel -> write out a patch of its local neighborhood
                    for (size_t cz = z-3; cz <= z + 3; ++cz) {
                        for (size_t cy = y - 3; cy <= y + 3; ++cy) {
                            for (size_t cx = x - 3; cx <= x + 3; ++cx) {
                                // write out as raw binary
                                patch[patchIndex++] = volume->getVoxelNormalized(cx, cy, cz);
                            }
                        }
                    }

                    // normalize patch to mean 0 and standard deviation 1
                    if (normalizePatches_.get()) {
                        float sum = 0.f;
                        // get mean and standard deviation to normalize the patch
                        for (size_t index = 0; index < 7 * 7 * 7; ++index) {
                            sum += patch[index];
                        }
                        float mean = sum / (7.f * 7.f * 7.f);

                        float sigma = 0.f;
                        for (size_t index = 0; index < 7 * 7 * 7; ++index) {
                            sigma += std::pow(patch[index] - mean, 2.f);
                        }
                        sigma *= 1.f / (7.f * 7.f * 7.f - 1.f);
                        sigma = std::sqrt(sigma);

                        // if sigma is near 0: set to 1
                        if (sigma < 1e-6f)
                            sigma = 1.f;

                        // use mean and standard deviation to normalize the patch
                        for (size_t index = 0; index < 7 * 7 * 7; ++index) {
                            float value = (patch[index] - mean) / sigma;
                            patch[index] = value;
                        }
                    }

                    //rawout.write((char*) &patch[0], (std::streamsize) sizeof(patch));
                    std::memcpy(&featureData[numPatches * 7 * 7 * 7], &patch[0], 7 * 7 * 7 * sizeof(float));
                    numPatches++;
                }
            }
        }
    
    }

    featureSet.write(featureData, PredType::NATIVE_FLOAT);
    //file.close();
    file.close();

    delete[] featureData;

    // re-activate settings
    filename_.setReadOnlyFlag(false);
    normalizePatches_.setReadOnlyFlag(false);
    computeButton_.setReadOnlyFlag(false);

    setProgress(1.f);
}

} // namespace
