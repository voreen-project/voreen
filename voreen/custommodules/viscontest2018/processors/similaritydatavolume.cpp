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

#include "similaritydatavolume.h"

#include "tgt/immediatemode/immediatemode.h"
#include "voreen/core/voreenapplication.h"

#include "voreen/core/utils/statistics.h"

#include "../ports/conditions/portconditionensemble.h"


namespace voreen {
    
const std::string SimilartyDataVolume::loggerCat_("voreen.viscontest2018.SimilartyDataVolume");

SimilartyDataVolume::SimilartyDataVolume()
    : RenderProcessor()
    , inport_(Port::INPORT, "ensembleinport", "Ensemble Data Input")
    , outport_(Port::OUTPORT, "volumehandle.volumehandle", "Volume Output")
    , similarityVolume_(nullptr)
    , similarityVolumeRepresentation_(nullptr)
    , similarityMethod_("similarityMethod", "Similarity Method")
    , group1_("similartyDataVolumeGroup1", "Group 1")
    , group2_("similartyDataVolumeGroup2", "Group 2")
{
    // Ports
    addPort(inport_);

    addPort(outport_);
    ON_CHANGE(outport_, SimilartyDataVolume, onInportChange);

    addProperty(similarityMethod_);
    similarityMethod_.addOption("variance", "Variance");
    similarityMethod_.addOption("minmax", "Min/Max Comparison");
    similarityMethod_.select("variance");

    addProperty(group1_);
    addProperty(group2_);


}

void SimilartyDataVolume::onInportChange() {
    if(!inport_.hasData()) return;

    tgt::ivec3 dimensions = inport_.getData()->getDimensions();
    try {
        similarityVolumeRepresentation_ = new VolumeRAM_Float(tgt::svec3(dimensions.x, dimensions.y, dimensions.z));
    } catch(std::bad_alloc e) {
        LERRORC("voreen.similaritydatavolume", "Failed to allocate similarity volume data memory");
        throw;
    }
    memset(similarityVolumeRepresentation_->voxel(), 0, similarityVolumeRepresentation_->getNumBytes());
    similarityVolume_ = new Volume(similarityVolumeRepresentation_, tgt::vec3::one, tgt::vec3::zero);

    for(std::string runName: inport_.getData()->getRuns()) {
        group1_.addRow(runName);
        group2_.addRow(runName);
    }
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
    tgt::svec3 dimensions = inport_.getData()->getDimensions();

    for(size_t z = 0; z < dimensions.z; z++) {
        for(size_t y = 0; y < dimensions.y; y++) {
            for(size_t x = 0; x < dimensions.x; x++) {
                size_t index = VolumeRAM_Float::calcPos(dimensions, x, y, z);
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

const std::vector<float> SimilartyDataVolume::applyGroupLogic(const std::vector<float>& rawVoxelData) {
    std::vector<float> modifiedVoxelData = rawVoxelData;
    // compare two groups
    if(!group1_.getSelectedRowIndices().empty() && !group2_.getSelectedRowIndices().empty()) {
        modifiedVoxelData.clear();
        float group1Avg;
        float group2Avg;

        Statistics statistics(true);
        for(int selectedIndex : group1_.getSelectedRowIndices()) {
            statistics.addSample(rawVoxelData.at(selectedIndex));
        }
        group1Avg = statistics.getMean();

        statistics.reset();
        for(int selectedIndex : group2_.getSelectedRowIndices()) {
            statistics.addSample(rawVoxelData.at(selectedIndex));
        }
        group2Avg = statistics.getMean();

        modifiedVoxelData.push_back(group1Avg);
        modifiedVoxelData.push_back(group2Avg);
    }
    // normal run filtering
    else if(!group1_.getSelectedRowIndices().empty() > 0) {
        modifiedVoxelData.clear();
        for(int selectedIndex : group1_.getSelectedRowIndices()) {
            modifiedVoxelData.push_back(rawVoxelData.at(selectedIndex));
        }
    }
    return modifiedVoxelData;
}

float SimilartyDataVolume::calculateVariance(const std::vector<float>& voxelData) {
    Statistics statistics(true);

    for(float voxelRunValue : voxelData) {
        statistics.addSample(voxelRunValue);
    }

    return statistics.getVariance();
}

float SimilartyDataVolume::calculateMinMaxDiff(const std::vector<float>& voxelData) {
    float min = std::numeric_limits<float>::max();
    float max = std::numeric_limits<float>::lowest();

    for(float voxelRunValue : voxelData) {
        min = std::min(min, voxelRunValue);
        max = std::max(max, voxelRunValue);
    }

    return std::abs(min - max);
}

} // namespace voreen
