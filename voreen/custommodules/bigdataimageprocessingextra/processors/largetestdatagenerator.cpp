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

#include "largetestdatagenerator.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/octree/volumeoctree.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/utils/hashing.h"

#include "modules/hdf5/io/hdf5volumewriter.h"
#include "modules/hdf5/io/hdf5volumereader.h"

#include "tgt/filesystem.h"

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/lexical_cast.hpp>

#include <memory>

namespace voreen {

const std::string LargeTestDataGenerator::loggerCat_("voreen.bigdataimageprocessingextra.largetestdatagenerator");

LargeTestDataGenerator::LargeTestDataGenerator()
    : AsyncComputeProcessor<LargeTestDataGeneratorInput, LargeTestDataGeneratorOutput>()
    , outportNoisy_(Port::OUTPORT, "largetestdatagenerator.outport", "Volume Output Noisy", true, Processor::VALID)
    , outportGT_(Port::OUTPORT, "largetestdatagenerator.outportgt", "Volume Output Ground Truth", true, Processor::VALID)
    , foregroundLabelsPort_(Port::OUTPORT, "largetestdatagenerator.foregroundLabelsPort", "Foreground Seeds Outport", true, Processor::VALID)
    , backgroundLabelsPort_(Port::OUTPORT, "largetestdatagenerator.backgroundLabelsPort", "Background Seeds Outport", true, Processor::VALID)
    , volumeDimensions_("volumeDimensions", "Volume Dimensions", tgt::ivec3(2), tgt::ivec3(2), tgt::ivec3(10000))
    , structureSizeRange_("structureSizeRange", "Structure Size", 1, 1, 10000)
    , noiseLevel_("noiseLevel", "Noise Level", 0.01, 0.0, 1.0)
    , seed_("seed", "RNG Seed", 0, 0, std::numeric_limits<int>::max())
    , outputVolumeNoisyFilePath_("outputVolumeFilePath", "Volume Noisy Output", "Path", "", "HDF5 (*.h5)", FileDialogProperty::SAVE_FILE, Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , outputVolumeGTFilePath_("outputVolumeFilePathgt", "GT Volume Output", "Path", "", "HDF5 (*.h5)", FileDialogProperty::SAVE_FILE, Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , scenario_("scenario", "Test Data Scenario")
{
    addPort(outportNoisy_);
    addPort(outportGT_);
    addPort(foregroundLabelsPort_);
    addPort(backgroundLabelsPort_);

    addProperty(volumeDimensions_);
        ON_CHANGE_LAMBDA(volumeDimensions_, [this] () {
            structureSizeRange_.setMaxValue(tgt::min(volumeDimensions_.get()));
        });
        volumeDimensions_.setTracking(false);
    addProperty(structureSizeRange_);
        structureSizeRange_.setMaxValue(tgt::min(volumeDimensions_.get()));
    addProperty(noiseLevel_);
        noiseLevel_.setTracking(false);
    addProperty(seed_);
        seed_.setTracking(false);
    addProperty(scenario_);
        scenario_.addOption("cells", "Cells", LargeTestDataGeneratorInput::CELLS);
        scenario_.addOption("vessels", "Vessels", LargeTestDataGeneratorInput::VESSELS);
        scenario_.selectByValue(LargeTestDataGeneratorInput::CELLS);
    addProperty(outputVolumeNoisyFilePath_);
    addProperty(outputVolumeGTFilePath_);
}

LargeTestDataGeneratorInput LargeTestDataGenerator::prepareComputeInput() {
    // Reset output volume to make sure it (and the hdf5filevolume) are not used any more
    outportGT_.setData(nullptr);
    outportNoisy_.setData(nullptr);

    const std::string volumeNoisyFilePath = outputVolumeNoisyFilePath_.get();
    const std::string volumeGTFilePath = outputVolumeGTFilePath_.get();
    const std::string volumeLocation = HDF5VolumeWriter::VOLUME_DATASET_NAME;

    tgt::svec3 dim = volumeDimensions_.get();

    float noiseRange = noiseLevel_.get();

    LargeTestDataGeneratorInput::random_engine_type randomEngine(randomDevice());
    randomEngine.seed(seed_.get());

    const std::string baseType = "uint8";
    const tgt::vec3 spacing = tgt::vec3::one;
    const tgt::vec3 offset = tgt::vec3::zero;
    const RealWorldMapping rwm = RealWorldMapping(1,0,"");
    const int deflateLevel = 1;

    if(volumeNoisyFilePath.empty() || volumeGTFilePath.empty()) {
        throw InvalidInputException("No volume file path specified!", InvalidInputException::S_ERROR);
    }

    std::unique_ptr<HDF5FileVolume> outputVolumeNoisy = nullptr;
    std::unique_ptr<HDF5FileVolume> outputVolumeGT = nullptr;
    try {
        outputVolumeNoisy = std::unique_ptr<HDF5FileVolume>(HDF5FileVolume::createVolume(volumeNoisyFilePath, volumeLocation, baseType, dim, 1, true, deflateLevel, tgt::svec3(dim.x, dim.y, 1), false));
        outputVolumeGT = std::unique_ptr<HDF5FileVolume>(HDF5FileVolume::createVolume(volumeGTFilePath, volumeLocation, baseType, dim, 1, true, deflateLevel, tgt::svec3(dim.x, dim.y, 1), false));
    } catch(tgt::IOException& e) {
        throw InvalidInputException("Could not create output volume.", InvalidInputException::S_ERROR);
    }


    outputVolumeNoisy->writeSpacing(spacing);
    outputVolumeNoisy->writeOffset(offset);
    outputVolumeNoisy->writeRealWorldMapping(rwm);

    outputVolumeGT->writeSpacing(spacing);
    outputVolumeGT->writeOffset(offset);
    outputVolumeGT->writeRealWorldMapping(rwm);

    LINFO("Using structure size range: " << structureSizeRange_.get());
    LINFO("Using voldim: " << dim);
    LINFO("Using noise: " << noiseRange);
    LINFO("Using seed: " << seed_.get());

    return LargeTestDataGeneratorInput(
        scenario_.getValue(),
        std::move(outputVolumeNoisy),
        std::move(outputVolumeGT),
        randomEngine,
        noiseRange,
        structureSizeRange_.get()
    );
}

struct Balls {
    Balls()
        : center_()
        , radius_()
    {
    }

    void add(tgt::ivec3 center, int radius) {
        center_.push_back(center);
        radius_.push_back(radius);
    }
    size_t size() const {
        tgtAssert(center_.size() == radius_.size(), "ball component sizes mismatch");
        return center_.size();
    }

    bool inside(tgt::ivec3 p) const {
        tgtAssert(center_.size() == radius_.size(), "ball component sizes mismatch");
        int size = center_.size();
        int ret = 0;

        auto rad = radius_.begin();
        auto cen = center_.begin();
#ifdef VRN_MODULE_OPENMP
#pragma omp simd
#endif
        for(int i=0; i < size; ++i) {
            int radius = rad[i];
            tgt::ivec3 center = cen[i];
            tgt::vec3 diff = p - center;
            if((diff.x * diff.x + diff.y * diff.y + diff.z * diff.z) < radius * radius) {
                ret += 1;
            }
        }
        return ret > 0;
    }

private:
    std::vector<tgt::ivec3> center_;
    std::vector<int> radius_;
};
LargeTestDataGeneratorOutput LargeTestDataGenerator::compute(LargeTestDataGeneratorInput input, ProgressReporter& progress) const {
    tgtAssert(input.outputVolumeNoisy, "No outputVolume");
    tgtAssert(input.outputVolumeGT, "No outputVolume");
    const tgt::svec3 dim = input.outputVolumeNoisy->getDimensions();
    tgtAssert(input.outputVolumeGT->getDimensions() == dim, "Dimension mismatch");
    Balls balls{};

    std::unique_ptr<PointSegmentListGeometryVec3> foregroundLabels(new PointSegmentListGeometryVec3());
    std::unique_ptr<PointSegmentListGeometryVec3> backgroundLabels(new PointSegmentListGeometryVec3());

    int minRadius = std::max(1, input.structureSizeRange.x/2);
    int maxRadius = std::max(1, input.structureSizeRange.y/2);

    int elementVolumeEstimate;
    switch(input.scenario) {
        case LargeTestDataGeneratorInput::CELLS: {
            float a = minRadius;
            float b = maxRadius;
            //if(a == b) {
            //    elementVolumeEstimate = a;
            //} else {
                elementVolumeEstimate = tgt::round(tgt::PIf*(b*b*b*b-a*a*a*a)/(3 * (b - a)));
            //}
            break;
        }
        case LargeTestDataGeneratorInput::VESSELS: {
            tgtAssert(false, "Unimplemented radius estimation for vessel scenario");
            elementVolumeEstimate = 1;
            break;
        }
        default: {
            tgtAssert(false, "Invalid scenario");
        }
    }

    int totalVolume = tgt::hmul(dim);
    int numElements = std::max(1, totalVolume/elementVolumeEstimate/10);

    LINFO("Placing " << numElements << " Elements");

    std::uniform_int_distribution<> xDistr(0, dim.x-1);
    std::uniform_int_distribution<> yDistr(0, dim.y-1);
    std::uniform_int_distribution<> zDistr(0, dim.z-1);
    std::uniform_int_distribution<> rDistr(minRadius, maxRadius);

    std::normal_distribution<float> noiseDistr(0.0, input.noiseRange);

    size_t numObjects = 0;

    for(int i=0; i<numElements; ++i) {
        tgt::ivec3 p(xDistr(input.randomEngine), yDistr(input.randomEngine), zDistr(input.randomEngine));
        balls.add(p, rDistr(input.randomEngine));

        std::vector<tgt::vec3> segment;
        segment.push_back(p);
        segment.push_back(tgt::vec3(p)+tgt::vec3(0.001));
        foregroundLabels->addSegment(segment);
    }

    numObjects += balls.size();

    int max_tries = 10000;
    int tries = max_tries;
    int i=0;
    for(; i<numObjects && tries > 0;) {
        tgt::ivec3 p(xDistr(input.randomEngine), yDistr(input.randomEngine), zDistr(input.randomEngine));

        if(!balls.inside(p)) {
            std::vector<tgt::vec3> segment;
            segment.push_back(p);
            segment.push_back(tgt::vec3(p)+tgt::vec3(0.001));
            backgroundLabels->addSegment(segment);
            ++i;

            tries = max_tries;
        } else {
            tries--;
        }
    }

    if(tries == 0) {
        LWARNING("Failed to position " << (numObjects-i) << "background seeds");
    }

    const float insideBase = 0.7;
    const float outsideBase = 0.3;
    for(int z=0; z<dim.z; ++z) {
        progress.setProgress(static_cast<float>(z)/dim.z);

        VolumeAtomic<uint8_t> sliceNoisy(tgt::vec3(dim.x, dim.y, 1));
        VolumeAtomic<uint8_t> sliceGT(tgt::vec3(dim.x, dim.y, 1));
        for(int y=0; y<dim.y; ++y) {
            for(int x=0; x<dim.x; ++x) {
                tgt::ivec3 p(x,y,z);

                bool inside = balls.inside(p);
                float val = inside ? insideBase : outsideBase;

                val += noiseDistr(input.randomEngine);
                val = tgt::clamp(val, 0.0f, 1.0f);
                sliceNoisy.setVoxelNormalized(val, x,y,0);

                sliceGT.setVoxelNormalized(inside ? 1.0f : 0.0f, x,y,0);
            }
        }
        input.outputVolumeNoisy->writeSlices(&sliceNoisy, z);
        input.outputVolumeGT->writeSlices(&sliceGT, z);
    }

    return {
        std::move(foregroundLabels),
        std::move(backgroundLabels),
        input.outputVolumeNoisy->getFileName(),
        input.outputVolumeGT->getFileName(),
    };
    //outputVolume will be destroyed and thus closed now.
}
void LargeTestDataGenerator::processComputeOutput(LargeTestDataGeneratorOutput output) {
    // outputVolume has been destroyed and thus closed by now.
    // So we can open it again (and use HDF5VolumeReader's implementation to read all the metadata with the file)
    const VolumeBase* volNoisy = HDF5VolumeReader().read(output.outputVolumeNoisyFilePath)->at(0);
    const VolumeBase* volGT = HDF5VolumeReader().read(output.outputVolumeGTFilePath)->at(0);

    outportNoisy_.setData(volNoisy);
    outportGT_.setData(volGT);

    foregroundLabelsPort_.setData(output.foregroundLabels.release());
    backgroundLabelsPort_.setData(output.backgroundLabels.release());
}

LargeTestDataGenerator::~LargeTestDataGenerator() {
}
VoreenSerializableObject* LargeTestDataGenerator::create() const {
    return new LargeTestDataGenerator();
}

} // namespace voreen
