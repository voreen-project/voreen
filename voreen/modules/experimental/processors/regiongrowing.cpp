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

#include "regiongrowing.h"

#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"
#include "voreen/core/datastructures/volume/modality.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"
#include "modules/core/io/datvolumewriter.h"
#include "modules/core/io/datvolumereader.h"
#include "tgt/quadric.h"

#include <fstream>
#include <stack>

using tgt::vec3;
using tgt::vec4;
using tgt::ivec3;

namespace voreen {

const std::string RegionGrowingProcessor::loggerCat_("voreen.experimental.RegionGrowingProcessor");

RegionGrowingProcessor::RegionGrowingProcessor()
    : VolumeRenderer(),
      segmentationHandle_(0),
      strictness_("strictness", "strictness", 0.8f, 0.f, 65535.f),
      thresholdFilling_("thresholdFilling", "apply threshold on flood fill", false),
      thresholds_("thresholds", "Thresholds", tgt::vec2(0.f, 1.f), tgt::vec2(0.f), tgt::vec2(1.f)),
      thresholdL_("threshold.lower", "Lower Threshold", 0.f, 0.f, 1.f),
      thresholdU_("threshold.upper", "Upper Threshold", 1.f, 0.f, 1.f),
      fillCostFunction_("fillCostFunction", "Cost function"),
      adaptive_("adaptive", "use adaptive growing criteria", false),
      maxSeedDistance_("maxSeedDistance", "maximum distance to seed point", 0, 0, 999),
      startGrowing_("startGrowing", "Start Growing"),
      lastSegmentation_(0),
      currentHandleHash_(""),
      seedPoint_(Port::INPORT, "pointlist.seedpoints"),
      inport_(Port::INPORT, "volumehandle.volume"),
      outport_(Port::OUTPORT, "volumehandle.segmentation", "volumehandle.segmentation", segmentationHandle_)
{
    addPort(inport_);
    addPort(outport_);
    addPort(seedPoint_);

    addProperty(strictness_);
    addProperty(thresholdFilling_);
    //addProperty(thresholds_);
    addProperty(thresholdL_);
    addProperty(thresholdU_);

    fillCostFunction_.addOption("intensity", "intensity");
    fillCostFunction_.addOption("gradient-magnitude", "gradient magnitude");
    fillCostFunction_.addOption("weighted", "weighted");
    addProperty(fillCostFunction_);

    addProperty(adaptive_);
    addProperty(maxSeedDistance_);

    startGrowing_.onClick(MemberFunctionCallback<RegionGrowingProcessor>(this, &RegionGrowingProcessor::startGrowing));
    addProperty(startGrowing_);
}

RegionGrowingProcessor::~RegionGrowingProcessor() {
    delete segmentationHandle_;
    delete lastSegmentation_;
}

namespace {

// calculate mean, variance and standard deviation for some float values
class Stats {
public:
    Stats()
        : n_(0), mean_(0.f), m2_(0.f) {}

    void add(float x) {
        // on-line algorithm, see http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
        n_++;
        float delta = x - mean_;
        mean_ += delta / n_;
        m2_ += delta * (x - mean_);
    }

    float mean() const { return mean_; }
    float variance() const { return m2_ / (n_ - 1); }
    float stddev() const { return std::sqrt(variance()); }

protected:
    int n_;
    float mean_;
    float m2_;
};


// Calculate standard deviation of the intensity values from the 26 neighbors of the voxel at
// position pos. Note: position must not be a border voxel.
template<class T>
float neighborStandardDeviation(const tgt::ivec3& pos, T* dataset, float max, Stats& stats) {
    for (int z=-1; z <= 1; z++)
        for (int y=-1; y <= 1; y++)
            for (int x=-1; x <= 1; x++) {
                ivec3 p(x, y, z);
                if (p != ivec3::zero)
                    stats.add(dataset->voxel(pos + p) / max);
            }

    return stats.stddev();
}

// Calculate gradient magnitude. Note: position must not be a border voxel.
template<class T>
float gradientMagnitude(const tgt::ivec3& pos, T* dataset, float max) {
    float v0 = dataset->voxel(pos + ivec3(1, 0, 0)) / max;
    float v1 = dataset->voxel(pos + ivec3(0, 1, 0)) / max;
    float v2 = dataset->voxel(pos + ivec3(0, 0, 1)) / max;
    float v3 = dataset->voxel(pos - ivec3(1, 0, 0)) / max;
    float v4 = dataset->voxel(pos - ivec3(0, 1, 0)) / max;
    float v5 = dataset->voxel(pos - ivec3(0, 0, 1)) / max;
    vec3 gradient = vec3(v3 - v0, v4 - v1, v5 - v2);
    return length(gradient);
}

// Calculate standard deviation from the gradient magnitudes of the 26 neighbors of the voxel at
// position pos. Note: position must be at least 2 voxels away from the border.
template<class T>
float neighborStandardDeviationGradients(const tgt::ivec3& pos, T* dataset, float max, Stats& stats) {
    for (int z=-1; z <= 1; z++)
        for (int y=-1; y <= 1; y++)
            for (int x=-1; x <= 1; x++) {
                ivec3 p(x, y, z);
                if (p != ivec3::zero)
                    stats.add(gradientMagnitude(pos + p, dataset, max));
            }

    return stats.stddev();
}

enum FloodFillMode { FLOODFILL_INTENSITY, FLOODFILL_GRADMAG, FLOODFILL_WEIGHTED };

template<class T>
int floodFill(const ivec3& seed_pos, int segment, float lowerThreshold, float upperThreshold,
              float strictness, const T* dataset, VolumeRAM_UInt8* segvol, FloodFillMode mode,
              bool useThresholds, bool adaptive, float maxSeedDistance)
{
    ivec3 dims = dataset->getDimensions();
    float max = static_cast<float>(VolumeElement<typename T::VoxelType>::rangeMax());

    VolumeAtomic<bool> markedVoxels(dims);
    markedVoxels.clear();

    std::vector<ivec3> neighbors;
    for (int z=-1; z <= 1; z++)
        for (int y=-1; y <= 1; y++)
            for (int x=-1; x <= 1; x++) {
                ivec3 p(x, y, z);
                if (p != ivec3::zero)
                    neighbors.push_back(p);
            }


    std::stack<ivec3> voxelStack;
    voxelStack.push(seed_pos);
    for (size_t i=0; i < neighbors.size(); i++)
        voxelStack.push(seed_pos + neighbors[i]);

    Stats stats_value;
    Stats stats_gradmag;

    float seed_value = dataset->voxel(seed_pos) / max;
    float seed_stddev26 = neighborStandardDeviation(seed_pos, dataset, max, stats_value);

    float seed_gradmag = gradientMagnitude(seed_pos, dataset, max);
    float seed_gradmag_stddev26 = neighborStandardDeviationGradients(seed_pos, dataset, max, stats_gradmag);

    if (mode == FLOODFILL_INTENSITY || mode == FLOODFILL_WEIGHTED)
        std::cout << "seed value: " << seed_value << ", "
                  << "stddev 26: " << seed_stddev26 << std::endl;

    if (mode == FLOODFILL_GRADMAG || mode == FLOODFILL_WEIGHTED)
        std::cout << "seed_gradmag: " << seed_gradmag << ", "
                  << "stddev gradmag 26: " << seed_gradmag_stddev26 << std::endl;

    while (!voxelStack.empty()) {
        ivec3 pos = voxelStack.top();
        voxelStack.pop();

        if (pos.x < 2 || pos.x > dims.x - 3 ||
            pos.y < 2 || pos.y > dims.y - 3 ||
            pos.z < 2 || pos.z > dims.z - 3 ||
            markedVoxels.voxel(pos))
        {
            // on border or already visited
            continue;
        }

        float value = dataset->voxel(pos) / max;

        if (adaptive /* && not in initial neighbor */) {
            stats_value.add(value);
            seed_stddev26 = stats_value.stddev();

            if (mode == FLOODFILL_GRADMAG || mode == FLOODFILL_WEIGHTED) {
                stats_gradmag.add(gradientMagnitude(pos, dataset, max));
                seed_gradmag_stddev26 = stats_gradmag.stddev();
            }
        }


        if (useThresholds && (value < lowerThreshold || value > upperThreshold)) {
            // invalid value
            continue;
        }

        // Cost function: if less than 1 then voxel is within the region.
        //
        // Based on: Runzhen Huang, Kwan-Liu Ma. RGVis: Region growing based techniques for
        // volume visualization, 2003.
        float cost = 0.f;

        if (mode == FLOODFILL_INTENSITY)
            cost = fabs(value - seed_value) / (strictness * seed_stddev26);
        else if (mode == FLOODFILL_GRADMAG) {
            float gradmag = gradientMagnitude(pos, dataset, max);
            cost = fabs(gradmag - seed_gradmag) / (strictness * seed_gradmag_stddev26);
        }
        else if (mode == FLOODFILL_WEIGHTED) {
            float cost_a = fabs(value - seed_value) / (strictness * seed_stddev26);

            float gradmag = gradientMagnitude(pos, dataset, max);
            float cost_b = fabs(gradmag - seed_gradmag) / (strictness * seed_gradmag_stddev26);

            // weight p
            float p = (seed_gradmag_stddev26 / (seed_stddev26 + seed_gradmag_stddev26));
            cost = cost_a * p + cost_b * (1.f - p);
        }

        if (cost >= 1.f)
            continue;

        if (maxSeedDistance > 0 && tgt::distance(vec3(pos), vec3(seed_pos)) > maxSeedDistance)
            continue;

        // voxel is valid
        markedVoxels.voxel(pos) = true;

        // add neighbors to stack if not already visited
        for (size_t i=0; i < neighbors.size(); i++)
            if (!markedVoxels.voxel(pos + neighbors[i]))
                voxelStack.push(pos + neighbors[i]);
    }

    // now fill segmentation volume with all marked voxels
    float count = 0;
    for (size_t z=0; z < markedVoxels.getDimensions().z; z++) {
        for (size_t y=0; y < markedVoxels.getDimensions().y; y++) {
            for (size_t x=0; x < markedVoxels.getDimensions().x; x++) {
                if (markedVoxels.voxel(x, y, z)) {
                    uint8_t& v = segvol->voxel(x, y, z);
                    if (v == 0 || segment == 0) {
                        v = segment;
                        count++;
                    }
                }
            }
        }
    }

    return static_cast<int>(count);
}

} // namespace

void RegionGrowingProcessor::mark(const ivec3& seedpos, int segment) {

    if (!inport_.getData() || !segmentationHandle_)
        return;

    const VolumeRAM* vol = inport_.getData()->getRepresentation<VolumeRAM>();
    VolumeRAM_UInt8* segvol = dynamic_cast<VolumeRAM_UInt8*>(segmentationHandle_->getWritableRepresentation<VolumeRAM>());

    if (!vol || !segvol) {
        LERROR("no volume or segmentation");
        return;
    }

    float seedval = vol->getVoxelNormalized(seedpos);
    LINFO("seed pos " << seedpos << " with value " << seedval);

    if (seedval <= 0.f || seedpos == ivec3(0)) {
        LERROR("ignoring this seed value");
        return;
    }

    // copy seg volume to last volume to enable undo
    if (!lastSegmentation_ || lastSegmentation_->getDimensions() != segvol->getDimensions()) {
        delete lastSegmentation_;
        lastSegmentation_ = segvol->clone();
    } else {
        memcpy(lastSegmentation_->getData(), segvol->getData(), segvol->getNumBytes());
    }

    // start the flood fill
    //VolumeRAM_UInt8* s = dynamic_cast<VolumeRAM_UInt8*>(segmentation_->getRepresentation<VolumeRAM>());
    FloodFillMode mode = FLOODFILL_INTENSITY;
    if (fillCostFunction_.get() == "gradient-magnitude")
        mode = FLOODFILL_GRADMAG;
    else if (fillCostFunction_.get() == "weighted")
        mode = FLOODFILL_WEIGHTED;

    int count = 0;
    if (dynamic_cast<const VolumeRAM_UInt8*>(vol))
        count = floodFill(seedpos, segment, thresholdL_.get(), thresholdU_.get(), strictness_.get(),
                          static_cast<const VolumeRAM_UInt8*>(vol), segvol, mode, thresholdFilling_.get(), adaptive_.get(),
                          static_cast<float>(maxSeedDistance_.get()));
    else if (dynamic_cast<const VolumeRAM_UInt16*>(vol))
        count = floodFill(seedpos, segment, thresholdL_.get(), thresholdU_.get(), strictness_.get(),
                          static_cast<const VolumeRAM_UInt16*>(vol), segvol, mode, thresholdFilling_.get(), adaptive_.get(),
                          static_cast<float>(maxSeedDistance_.get()));

    LINFO("filled voxels: " << count);
}

void RegionGrowingProcessor::clearSegmentation() {
    if (segmentationHandle_) {
        VolumeRAM_UInt8* segvol = dynamic_cast<VolumeRAM_UInt8*>(segmentationHandle_->getWritableRepresentation<VolumeRAM>());
        if (segvol) {
            segvol->clear();
            volumeModified(segmentationHandle_);
        }
    }
}

void RegionGrowingProcessor::clearSegment(int segmentID) {
    if (segmentationHandle_) {
        VolumeRAM_UInt8* segvol = dynamic_cast<VolumeRAM_UInt8*>(segmentationHandle_->getWritableRepresentation<VolumeRAM>());
        if (segvol) {
            LINFO("clearing segment " << segmentID);
            for (size_t i=0; i < segvol->getNumVoxels(); i++) {
                uint8_t& v = segvol->voxel(i);
                if (v == segmentID)
                    v = 0;
            }
            volumeModified(segmentationHandle_);
        }
    }
}


void RegionGrowingProcessor::undoLastGrowing() {
    if (lastSegmentation_ && segmentationHandle_) {
        LINFO("undo");
        VolumeRAM_UInt8* segvol = dynamic_cast<VolumeRAM_UInt8*>(segmentationHandle_->getWritableRepresentation<VolumeRAM>());
        if (segvol && lastSegmentation_->getDimensions() == segvol->getDimensions()) {
            // copy last volume to seg volume
            memcpy(segvol->getData(), lastSegmentation_->getData(), lastSegmentation_->getNumBytes());
            volumeModified(segmentationHandle_);
            LINFO("undo finished");
        }
    }
}

void RegionGrowingProcessor::startGrowing() {
    if (seedPoint_.isConnected()) {
        if (segmentationHandle_) {

            PointSegmentListGeometryVec3* seedList = const_cast<PointSegmentListGeometryVec3*>(dynamic_cast<const PointSegmentListGeometryVec3*>(seedPoint_.getData()));
            if(seedList && seedList->getNumPoints() == 0) {
                LERROR("No seedpoints found.");
                return;
            } else if(!seedList) {
                const PointListGeometryVec3* seedListPoints = dynamic_cast<const PointListGeometryVec3*>(seedPoint_.getData());
                if(!seedListPoints) {
                    LERROR("Input of wrong type.");
                    return;
                }
                if(seedListPoints->getNumPoints() == 0) {
                    LERROR("No seedpoints found.");
                    return;
                }
                seedList = new PointSegmentListGeometryVec3();
                seedList->addSegment(seedListPoints->getData());
            }

            clearSegmentation();
            int segmentID = 128;
            for(int i = 0; i < seedList->getNumSegments(); i++) {
                ivec3 voxelpos = seedList->getSegment(i).front();
                mark(voxelpos, segmentID);
            }
            volumeModified(segmentationHandle_);
            outport_.invalidatePort();
        }
        else {
            LERROR("No segmentation volume.");
        }
    } else {
        LERROR("Unable to retrieve first-hit-points: No RenderStore coprocessor connected");
    }
}

void RegionGrowingProcessor::saveSegmentation(std::string filename) const {
    if (segmentationHandle_) {
        LINFO("saving segmentation volume to " << filename);
        DatVolumeWriter writer;
        writer.write(filename, segmentationHandle_);
    }
}

void RegionGrowingProcessor::loadSegmentation(std::string filename) {
    LINFO("Reading segmentation: " << filename);

    delete segmentationHandle_;
    segmentationHandle_ = 0;

    DatVolumeReader reader;
    VolumeList* volumeList = 0;
    try {
        volumeList = reader.read(filename);
    }
    catch (std::exception& /*e*/) {}

    if (volumeList && !volumeList->empty()) {
        VolumeBase* volHandle = volumeList->first();
        if (volHandle) {
            segmentationHandle_ = dynamic_cast<Volume*>(volHandle);
            if(!segmentationHandle_)
                delete volHandle;
            invalidate();
            delete volumeList;
            return;
        }
    }

    delete volumeList;

    LERROR("Failed reading segmentation: " << filename);

}

void RegionGrowingProcessor::saveSegment(int segmentID, std::string filename) {
    if (segmentationHandle_) {
        const VolumeRAM_UInt8* segvol = dynamic_cast<const VolumeRAM_UInt8*>(segmentationHandle_->getRepresentation<VolumeRAM>());
        if (segvol) {
            LINFO("saving segment " << segmentID << " to " << filename);
            VolumeRAM_UInt8* tmpvol = segvol->clone();

            ivec3 dims = tmpvol->getDimensions();
            for (int z=0; z < dims.z; z++) {
                for (int y=0; y < dims.y; y++) {
                    for (int x=0; x < dims.x; x++) {
                        // pfskel doesn't like voxel on the border ("object touching bounding
                        // box"), so zero them out here.
                        const int s = 4;
                        bool border = (x < s || y < s || z < s ||
                            x >= dims.x - s || y >= dims.y - s || z >= dims.z - s);
                        uint8_t& v = tmpvol->voxel(x, y, z);
                        if (v == segmentID && !border)
                            v = 1;
                        else
                            v = 0;
                    }
                }
            }
            DatVolumeWriter writer;
            writer.write(filename, new Volume(tmpvol, segmentationHandle_));
            delete tmpvol;
        }
    }
}



void RegionGrowingProcessor::process() {
    // if segmentation volume not present or volume has changed, create
    // segmentation volume with the same dimensions and spacing as the input volume
    if ((!segmentationHandle_ || inport_.hasChanged()) && inport_.getData()) {
        if(inport_.hasChanged()) {
            if(inport_.getData()->getHash() == currentHandleHash_)
                return;
            else
                currentHandleHash_ = inport_.getData()->getHash();
        }
        delete segmentationHandle_;
        VolumeRAM_UInt8* segVol = new VolumeRAM_UInt8(inport_.getData()->getRepresentation<VolumeRAM>()->getDimensions());
        //segVol->clear();
        segmentationHandle_ = new Volume(segVol, inport_.getData());
        clearSegmentation();
        outport_.setData(segmentationHandle_, false);
    }
}

void RegionGrowingProcessor::volumeModified(VolumeBase* v) {
    // Free the hardware volume, will be re-generated on the next repaint.
    v->removeRepresentation<VolumeGL>();
    v->getRepresentation<VolumeGL>();

    // hack to notify everybody that the volume was changed
    //MsgDistr.postMessage(new TemplateMessage<VolumeSetContainer*>("volumesetcontainer.clear", 0));
    invalidate();
}

FloatVec2Property& RegionGrowingProcessor::getThresholdProp() {
    return thresholds_;
}

BoolProperty& RegionGrowingProcessor::getThresholdFillingProp() {
    return thresholdFilling_;
}

FloatProperty& RegionGrowingProcessor::getStrictnessProp() {
    return strictness_;
}

BoolProperty& RegionGrowingProcessor::getAdaptiveProp() {
    return adaptive_;
}

IntProperty& RegionGrowingProcessor::getMaxSeedDistanceProp() {
    return maxSeedDistance_;
}

StringOptionProperty& RegionGrowingProcessor::getCostFunctionProp(){
    return fillCostFunction_;
}

} // namespace
