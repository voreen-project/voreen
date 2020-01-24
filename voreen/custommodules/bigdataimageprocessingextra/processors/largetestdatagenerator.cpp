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
    , outport_(Port::OUTPORT, "largetestdatagenerator.outport", "Volume Output")
    , foregroundLabelsPort_(Port::OUTPORT, "largetestdatagenerator.foregroundLabelsPort", "Foreground Seeds Outport")
    , backgroundLabelsPort_(Port::OUTPORT, "largetestdatagenerator.backgroundLabelsPort", "Background Seeds Outport")
    , volumeDimensions_("volumeDimensions", "Volume dimensions", tgt::ivec3(1), tgt::ivec3(1), tgt::ivec3(std::numeric_limits<int>::max()))
    , outputVolumeFilePath_("outputVolumeFilePath", "Output Volume", "Path", "", "HDF5 (*.h5)", FileDialogProperty::SAVE_FILE, Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , scenario_("scenario", "Test Data Scenario")
{
    addPort(outport_);
    addPort(foregroundLabelsPort_);
    addPort(backgroundLabelsPort_);

    addProperty(volumeDimensions_);
    addProperty(outputVolumeFilePath_);
    addProperty(scenario_);
        scenario_.addOption("cells", "Cells", LargeTestDataGeneratorInput::CELLS);
        scenario_.addOption("vessels", "Vessels", LargeTestDataGeneratorInput::VESSELS);
        scenario_.selectByValue(LargeTestDataGeneratorInput::CELLS);

}

LargeTestDataGeneratorInput LargeTestDataGenerator::prepareComputeInput() {
    // Reset output volume to make sure it (and the hdf5filevolume) are not used any more
    outport_.setData(nullptr);

    const std::string volumeFilePath = outputVolumeFilePath_.get();
    const std::string volumeLocation = HDF5VolumeWriter::VOLUME_DATASET_NAME;

    tgt::svec3 dim = volumeDimensions_.get();
    if(char* cubeSizeOverride = std::getenv("CUBE_SIZE")) {
        try {
            dim = tgt::svec3(std::stoi(cubeSizeOverride));
            LINFO("Environment variable CUBE_SIZE overrides output size with " << dim);
        } catch (std::invalid_argument&) {
            LERROR("Invalid argument passt via CUBE_SIZE environment variable: " << cubeSizeOverride);
        }
    }

    float noiseRange = 0.1;
    if(char* noiseOverride = std::getenv("NOISE_RANGE")) {
        try {
            noiseRange = std::stof(noiseOverride);
            LINFO("Environment variable NOISE_RANGE overrides noise parameter with " << noiseRange);
        } catch (std::invalid_argument&) {
            LERROR("Invalid argument passt via CUBE_SIZE environment variable: " << noiseOverride);
        }
    }

    LargeTestDataGeneratorInput::random_engine_type randomEngine(randomDevice());
    if(char* randomEngineSeed = std::getenv("RANDOM_ENGINE_SEED")) {
        try {
            int seed = std::stoi(randomEngineSeed);
            LINFO("Environment variable NOISE_RANGE overrides noise parameter with " << seed);
            randomEngine.seed(seed);
        } catch (std::invalid_argument&) {
            LERROR("Invalid argument passt via RANDOM_ENGINE_SEED environment variable: " << randomEngineSeed);
        }
    }

    const std::string baseType = "uint8";
    const tgt::vec3 spacing = tgt::vec3::one;
    const tgt::vec3 offset = tgt::vec3::zero;
    const RealWorldMapping rwm = RealWorldMapping(1,0,"");
    const int deflateLevel = 1;

    if(volumeFilePath.empty()) {
        throw InvalidInputException("No volume file path specified!", InvalidInputException::S_ERROR);
    }

    std::unique_ptr<HDF5FileVolume> outputVolume = nullptr;
    try {
        outputVolume = std::unique_ptr<HDF5FileVolume>(HDF5FileVolume::createVolume(volumeFilePath, volumeLocation, baseType, dim, 1, true, deflateLevel, tgt::svec3(dim.x, dim.y, 1), false));
    } catch(tgt::IOException& e) {
        throw InvalidInputException("Could not create output volume.", InvalidInputException::S_ERROR);
    }


    outputVolume->writeSpacing(spacing);
    outputVolume->writeOffset(offset);
    outputVolume->writeRealWorldMapping(rwm);

    return LargeTestDataGeneratorInput(
        scenario_.getValue(),
        std::move(outputVolume),
        randomEngine,
        noiseRange
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
        for(int i=0; i < size; ++i) {
            int radius = radius_[i];
            tgt::ivec3 center = center_[i];
            if(tgt::distanceSq(p, center) < radius * radius) {
                return true;
            }
        }
        return false;
    }

private:
    std::vector<tgt::ivec3> center_;
    std::vector<int> radius_;
};
LargeTestDataGeneratorOutput LargeTestDataGenerator::compute(LargeTestDataGeneratorInput input, ProgressReporter& progress) const {
    tgtAssert(input.outputVolume, "No outputVolume");
    const tgt::svec3 dim = input.outputVolume->getDimensions();
    Balls balls{};

    std::unique_ptr<PointSegmentListGeometryVec3> foregroundLabels(new PointSegmentListGeometryVec3());
    std::unique_ptr<PointSegmentListGeometryVec3> backgroundLabels(new PointSegmentListGeometryVec3());

    std::uniform_int_distribution<> xDistr(0, dim.x-1);
    std::uniform_int_distribution<> yDistr(0, dim.y-1);
    std::uniform_int_distribution<> zDistr(0, dim.z-1);
    std::uniform_int_distribution<> rDistr(1, std::max(1UL, dim.x/10));

    std::normal_distribution<float> noiseDistr(0.0, input.noiseRange);

    size_t numObjects = 0;

    for(int i=0; i<10; ++i) {
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

        VolumeAtomic<uint8_t> slice(tgt::vec3(dim.x, dim.y, 1));
        for(int y=0; y<dim.y; ++y) {
            for(int x=0; x<dim.x; ++x) {
                tgt::ivec3 p(x,y,z);

                float val = balls.inside(p) ? insideBase : outsideBase;

                val += noiseDistr(input.randomEngine);
                slice.setVoxelNormalized(val, x,y,0);
            }
        }
        input.outputVolume->writeSlices(&slice, z);
    }

    return {
        std::move(foregroundLabels),
        std::move(backgroundLabels),
        input.outputVolume->getFileName(),
    };
    //outputVolume will be destroyed and thus closed now.
}
void LargeTestDataGenerator::processComputeOutput(LargeTestDataGeneratorOutput output) {
    // outputVolume has been destroyed and thus closed by now.
    // So we can open it again (and use HDF5VolumeReader's implementation to read all the metadata with the file)
    const VolumeBase* vol = HDF5VolumeReader().read(output.outputVolumeFilePath)->at(0);

    outport_.setData(vol);
    foregroundLabelsPort_.setData(output.foregroundLabels.release());
    backgroundLabelsPort_.setData(output.backgroundLabels.release());
}

LargeTestDataGenerator::~LargeTestDataGenerator() {
}
VoreenSerializableObject* LargeTestDataGenerator::create() const {
    return new LargeTestDataGenerator();
}

} // namespace voreen
