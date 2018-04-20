/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
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

#include "similaritydataviewer.h"

#include "tgt/immediatemode/immediatemode.h"
#include "voreen/core/voreenapplication.h"

#include "voreen/core/utils/statistics.h"

#include "../ports/conditions/portconditionensemble.h"


namespace voreen {

static int PIXEL_BRUSH_ENABLE_DISTANCE = 3;
static const tgt::vec2 NO_SELECTION(-2.0f); // needs to be any value not inside range [-1, 1]

const std::string SimilartyDataVolume::loggerCat_("voreen.viscontest2018.SimilartyDataVolume");

SimilartyDataVolume::SimilartyDataVolume()
    : RenderProcessor()
    , inport_(Port::INPORT, "ensembleinport", "Ensemble Data Input")
    , outport_(Port::OUTPORT, "volumehandle.volumehandle", "Volume Output")
    , similarityVolume_(nullptr)
    , similarityVolumeRepresentation_(nullptr)
    , similarityMethod_("similarityMethod", "Similarity Method")
{
    // Ports
    addPort(inport_);

    addPort(outport_);
    ON_CHANGE(outport_, SimilartyDataVolume, onInportChange);

    addProperty(similarityMethod_);
    similarityMethod_.addOption("variance", "Variance");
    similarityMethod_.addOption("minmax", "Min/Max Comparison");
    similarityMethod_.select("variance");

}

void SimilartyDataVolume::onInportChange() {
    if(!inport_.hasData()) return;

    tgt::ivec3 dimensions = inport_.getData()->getDimensions();
    try {
        similarityVolumeRepresentation_ = new VolumeRAM_Float(tgt::svec3(dimensions.x, dimensions.y, dimensions.z));
    } catch(std::bad_alloc e) {
        LERRORC("voreen.similaritydataviewer", "Failed to allocate similarity volume data memory");
        throw;
    }
    memset(similarityVolumeRepresentation_->voxel(), 0, similarityVolumeRepresentation_->getNumBytes());
    similarityVolume_ = new Volume(similarityVolumeRepresentation_, tgt::vec3::one, tgt::vec3::zero);
}

SimilartyDataVolume::~SimilartyDataVolume() {
}

Processor* SimilartyDataVolume::create() const {
    return new SimilartyDataVolume();
}

void SimilartyDataVolume::initialize() {
    RenderProcessor::initialize();
}

void SimilartyDataVolume::deinitialize() {

    RenderProcessor::deinitialize();
}

bool SimilartyDataVolume::isReady() const {
    return inport_.isReady();
}

void SimilartyDataVolume::process() {
    // TODO
}

void SimilartyDataVolume::onEvent(tgt::Event *e) {
    Processor::onEvent(e);
}

void SimilartyDataVolume::initSimilarityVolume() {
    if(!inport_.hasData()) return;

    const std::vector<std::vector<float>>& similarityData = inport_.getData()->getData();
    tgt::ivec3 dimensions = inport_.getData()->getDimensions();

    for(size_t z = 0; z < static_cast<size_t>(dimensions.z); z++) {
        for(size_t y = 0; y < static_cast<size_t>(dimensions.y); y++) {
            for(size_t x = 0; x < static_cast<size_t>(dimensions.x); x++) {
                size_t index = x + (y + z * dimensions.z) * dimensions.x;
                const std::vector<float> voxelData = similarityData[index];
                float similarityVoxelValue = 0;
                if(similarityMethod_.get() == "variance") {
                    similarityVoxelValue = calculateVariance(voxelData);
                }
                else if(similarityMethod_.get() == "minmax") {
                    similarityVoxelValue = calculateMinMaxDiff(voxelData);
                }
                similarityVolumeRepresentation_->voxel()[index] = similarityVoxelValue;
            }
        }
    }

    outport_.setData(similarityVolume_);
}

float SimilartyDataVolume::calculateVariance(const std::vector<float> voxelData) {
    Statistics statistics(true);

    for(float voxelRunValue : voxelData) {
        statistics.addSample(voxelRunValue);
    }

    return statistics.getVariance();
}

float SimilartyDataVolume::calculateMinMaxDiff(const std::vector<float> voxelData) {
    float min = std::numeric_limits<float>::max();
    float max = std::numeric_limits<float>::lowest();

    for(float voxelRunValue : voxelData) {
        min = std::min(min, voxelRunValue);
        max = std::max(max, voxelRunValue);
    }

    return std::abs(min - max);
}

} // namespace voreen
