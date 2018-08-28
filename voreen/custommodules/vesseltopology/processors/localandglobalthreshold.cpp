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

#include "localandglobalthreshold.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"

#include <memory>
#include <vector>
#include <algorithm>

using tgt::vec3;
using tgt::vec4;
using tgt::ivec3;

namespace voreen {

const std::string LocalAndGlobalThreshold::loggerCat_("voreen.vesseltopology.LocalAndGlobalThreshold");

LocalAndGlobalThreshold::LocalAndGlobalThreshold()
    : AsyncComputeProcessor()
    , localApplicationRange_("localApplicationRange", "Local threshold application range", 0, std::numeric_limits<float>::lowest(), std::numeric_limits<float>::max())
    , windowExtent_("windowExtent", "Thresholding window extent", 1, 1, 10)
    , windowSize_("windowSize", "Resulting thresholding window size", "")
    , inport_(Port::INPORT, "volumehandle.inport_")
    , outport_(Port::OUTPORT, "volumehandle.outport")
{
    addPort(inport_);
    addPort(outport_);

    addProperty(localApplicationRange_);
    addProperty(windowExtent_);
        ON_CHANGE_LAMBDA(windowExtent_, [this] () {
                std::string sizeStr = std::to_string(2*windowExtent_.get() + 1);
                windowSize_.set(sizeStr + "x" + sizeStr + "x" + sizeStr);
                });
    addProperty(windowSize_);
        windowSize_.setReadOnlyFlag(true);
}

LocalAndGlobalThreshold::~LocalAndGlobalThreshold() {
}

static uint8_t foreground = 0xff;
static uint8_t background = 0x00;

static float dist(float f1, float f2) {
    return std::abs(f1 - f2);
}

static std::pair<float, float> find_clusters(std::vector<float>& values) {
    tgtAssert(values.size() >= 2, "no values to find threshold for");
    // Let's go with k-means, k=2!
    const float min = *std::min_element(values.begin(), values.end());
    const float max = *std::max_element(values.begin(), values.end());
    const float threshold = (max - min) * 0.01;
    const int max_it = 100;

    float c1 = min;
    float c2 = max;
    float prevc1 = std::numeric_limits<float>::infinity();
    float prevc2 = std::numeric_limits<float>::infinity();
    for(int i=0; i < max_it && dist(c1, prevc1) > threshold && dist(c2, prevc2) > threshold; ++i) {
        float sum_c1 = 0;
        size_t num_in_c1 = 0;
        float sum_c2 = 0;
        size_t num_in_c2 = 0;
        for(float val : values) {
            if(dist(val, c1) < dist(val, c2)) {
                // val belongs to c1
                sum_c1 += val;
                ++num_in_c1;
            } else {
                // val belongs to c2
                sum_c2 += val;
                ++num_in_c2;
            }
        }
        c1 = sum_c1/num_in_c1;
        c2 = sum_c2/num_in_c2;
    }
    return std::make_pair(c1, c2);
}
static uint8_t classify(const VolumeRAM& input, float backgroundUpperBoundNormalized, float foregroundLowerBoundNormalized, uint8_t windowExtent, const tgt::ivec3& pos) {
    float val = input.getVoxelNormalized(pos);
    if(val < backgroundUpperBoundNormalized) {
        return background;
    } else if(val > foregroundLowerBoundNormalized) {
        return foreground;
        return background;
    } else {
        //return 0x7f;
        tgt::ivec3 dim = input.getDimensions();
        int extent = windowExtent;
        int xStart = std::max(pos.x - extent, 0);
        int xEnd   = std::min(pos.x + extent + 1, dim.x);
        int yStart = std::max(pos.y - extent, 0);
        int yEnd   = std::min(pos.y + extent + 1, dim.y);
        int zStart = std::max(pos.z - extent, 0);
        int zEnd   = std::min(pos.z + extent + 1, dim.z);
        std::vector<float> values;
        for(int z = zStart; z < zEnd; ++z) {
            for(int y = yStart; y < yEnd; ++y) {
                for(int x = xStart; x < xEnd; ++x) {
                    values.push_back(input.getVoxelNormalized(tgt::svec3(x, y, z)));
                }
            }
        }
        auto clusters = find_clusters(values);
        float threshold = (clusters.first + clusters.second)*0.5f;
        return (val > threshold) ? foreground : background;
    }
}

LAGTInput LocalAndGlobalThreshold::prepareComputeInput() {
    if(!inport_.hasData()) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    const VolumeBase& input = *inport_.getData();

    float backgroundUpperBoundNormalized = input.getRealWorldMapping().realWorldToNormalized(localApplicationRange_.get().x);
    float foregroundLowerBoundNormalized = input.getRealWorldMapping().realWorldToNormalized(localApplicationRange_.get().y);
    uint8_t windowExtent = windowExtent_.get();

    return LAGTInput(
        input,
        backgroundUpperBoundNormalized,
        foregroundLowerBoundNormalized,
        windowExtent
    );
}

LAGTOutput LocalAndGlobalThreshold::compute(LAGTInput input, ProgressReporter& progressReporter) const {
    const VolumeRAM* inputRAM = input.volume.getRepresentation<VolumeRAM>();
    tgtAssert(inputRAM, "No ram representation");

    std::unique_ptr<VolumeRAM_UInt8> output(new VolumeRAM_UInt8(input.volume.getDimensions()));
    tgt::ivec3 dim = input.volume.getDimensions();
    for(int z=0; z<dim.z; ++z) {
        progressReporter.setProgress(static_cast<float>(z)/dim.z);
        #pragma omp parallel for
        for(int y=0; y<dim.y; ++y) {
            for(int x=0; x<dim.x; ++x) {
                tgt::ivec3 pos(x,y,z);
                output->voxel(pos) = classify(*inputRAM, input.backgroundUpperBoundNormalized, input.foregroundLowerBoundNormalized, input.windowExtent, pos);
            }
        }
    }

    std::unique_ptr<Volume> vol(new Volume(output.release(), &input.volume));
    vol->setRealWorldMapping(RealWorldMapping());
    return LAGTOutput {
        std::move(vol)
    };
}
void LocalAndGlobalThreshold::processComputeOutput(LAGTOutput output) {
    outport_.setData(output.volume.release());
}

void LocalAndGlobalThreshold::adjustPropertiesToInput() {
    if(inport_.hasData()) {
        const VolumeMinMax* mm = inport_.getData()->getDerivedData<VolumeMinMax>();
        tgtAssert(mm, "No VolumeMinMax");

        localApplicationRange_.setMinValue(mm->getMin());
        localApplicationRange_.setMaxValue(mm->getMax());
    }
}

} // namespace voreen
