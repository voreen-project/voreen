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

#include "aortasegmentation.h"

#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"
#include "voreen/core/datastructures/volume/modality.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"

#include <queue>
#include <numeric>

using tgt::vec3;
using tgt::vec4;
using tgt::ivec3;

namespace voreen {

HeartSegmentationBorder::HeartSegmentationBorder(const tgt::vec3& aortaPos, const tgt::vec3& heartPos, float radius)
    : center_()
    , planeNormal_()
    , effectRadiusSquared_(radius*radius)
    , cylinderHeight_(tgt::distance(aortaPos, heartPos))
{
    tgtAssert(aortaPos != heartPos, "Invalid aorta/heart pos: They are equal");

    tgt::vec3 fromHeartToCenter = (aortaPos-heartPos)*0.5f;
    center_ = heartPos + fromHeartToCenter;
    planeNormal_ = tgt::normalize(fromHeartToCenter);
}

bool HeartSegmentationBorder::isInSphereOfInfluence(const tgt::vec3& pos) const {
    return tgt::distanceSq(pos, center_) < effectRadiusSquared_;
}
bool HeartSegmentationBorder::areOnOppositeSidesOfPlane(const tgt::vec3& p1, const tgt::vec3& p2) const {
    return tgt::dot((p1 - center_), planeNormal_) * tgt::dot((p2 - center_), planeNormal_) < 0;
}
bool HeartSegmentationBorder::isCrossedBy(const tgt::ivec3& before, const tgt::ivec3& after) const {

    /*
    if(!isInSphereOfInfluence(before)) {
        return false;
    }
    if(!isInSphereOfInfluence(after)) {
        return false;
    }

    return areOnOppositeSidesOfPlane(before, after);
    */
    tgt::vec3 p = tgt::vec3(after)-center_;
    float projection_length = tgt::dot(planeNormal_, p);
    if(2*tgt::abs(projection_length) > cylinderHeight_) {
        return false;
    }
    tgt::vec3 p_projected = p - planeNormal_ * projection_length;
    return tgt::lengthSq(p_projected) < effectRadiusSquared_;
}

const std::string AortaSegmentation::loggerCat_("voreen.vesseltopology.AortaSegmentation");

AortaSegmentation::AortaSegmentation()
    : Processor()
    , strictness_("strictness", "strictness", 0.8f, 0.f, 65535.f)
    , thresholdFilling_("thresholdFilling", "apply threshold on flood fill", false)
    , volumeThresholds_("volumeThresholds", "Thresholds (orig. Volume)", 0, std::numeric_limits<float>::lowest(), std::numeric_limits<float>::max())
    , vesselnessThresholds_("vesselnessThresholds", "Thresholds (Vesselness)", 0, std::numeric_limits<float>::lowest(), std::numeric_limits<float>::max())
    , separatingDiskRadius_("separatingDiskRadius", "Radius of separating Disk", 0, 0, std::numeric_limits<float>::max())
    , fillCostFunction_("fillCostFunction", "Cost function")
    , adaptive_("adaptive", "use adaptive growing criteria", false)
    , maxSteps_("maxSteps", "Maximum number of voxel steps to grow", 0, 0, std::numeric_limits<int>::max())
    , updateContinuously_("updateContinuously", "Update continuously", false)
    , startGrowing_("startGrowing", "Start Growing")
    , seedPoints_(Port::INPORT, "pointlist.seedpoints")
    , heartPoints_(Port::INPORT, "pointlist.heartPoints_")
    , volumeInport_(Port::INPORT, "volumehandle.volume")
    , vesselnessInport_(Port::INPORT, "volumehandle.vesselness")
    , outport_(Port::OUTPORT, "volumehandle.segmentation", "volumehandle.segmentation")
    , updateForced_(false)
{
    addPort(volumeInport_);
    addPort(vesselnessInport_);
    addPort(outport_);
    addPort(seedPoints_);
    addPort(heartPoints_);

    addProperty(strictness_);
    addProperty(thresholdFilling_);
    addProperty(volumeThresholds_);
    addProperty(vesselnessThresholds_);
    addProperty(separatingDiskRadius_);

    fillCostFunction_.addOption("intensity", "intensity");
    fillCostFunction_.addOption("gradient-magnitude", "gradient magnitude");
    fillCostFunction_.addOption("weighted", "weighted");
    addProperty(fillCostFunction_);

    addProperty(adaptive_);
    addProperty(maxSteps_);

    addProperty(updateContinuously_);
    addProperty(startGrowing_);
        ON_CHANGE_LAMBDA(startGrowing_, [this] () {
                updateForced_ = true;
                });
}

AortaSegmentation::~AortaSegmentation() {
}

namespace {

struct QueueVoxel {
    QueueVoxel(tgt::ivec3 pos, uint32_t step)
        : pos_(pos)
        , step_(step)
    {
    }
    tgt::ivec3 pos_;
    uint32_t step_;
};

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

class ThresholdCriterion {
public:
    ThresholdCriterion(const VolumeBase& volume, float lowerBound, float upperBound)
        : volume_(volume.getRepresentation<VolumeRAM>())
        , lowerBoundNormalized_(volume.getRealWorldMapping().realWorldToNormalized(lowerBound))
        , upperBoundNormalized_(volume.getRealWorldMapping().realWorldToNormalized(upperBound))
    {
    }
    bool isFullfilledAt(const tgt::ivec3& pos) const {
        float val = volume_->getVoxelNormalized(pos);
        return lowerBoundNormalized_ <= val && val <= upperBoundNormalized_;
    }

private:
    const VolumeRAM* volume_;
    float lowerBoundNormalized_;
    float upperBoundNormalized_;
};

// Calculate standard deviation of the intensity values from the 26 neighbors of the voxel at
// position pos. Note: position must not be a border voxel.
template<class T>
float neighborStandardDeviation(const tgt::ivec3& pos, T* dataset, float max, Stats& stats) {
    for(int z=-1; z <= 1; z++)
        for(int y=-1; y <= 1; y++)
            for(int x=-1; x <= 1; x++) {
                ivec3 p(x, y, z);
                if(p != ivec3::zero)
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
    for(int z=-1; z <= 1; z++)
        for(int y=-1; y <= 1; y++)
            for(int x=-1; x <= 1; x++) {
                ivec3 p(x, y, z);
                if(p != ivec3::zero)
                    stats.add(gradientMagnitude(pos + p, dataset, max));
            }

    return stats.stddev();
}


enum FloodFillMode { FLOODFILL_INTENSITY, FLOODFILL_GRADMAG, FLOODFILL_WEIGHTED };

template<class T>
int floodFill(const ivec3& seed_pos, int segment, const ThresholdCriterion& crit1, const ThresholdCriterion& crit2,
              float strictness, const T* dataset, VolumeRAM_UInt8& segvol, FloodFillMode mode,
              bool useThresholds, bool adaptive, int maxSteps, const HeartSegmentationBorder& border)
{
    ivec3 dims = dataset->getDimensions();
    float max = static_cast<float>(VolumeElement<typename T::VoxelType>::rangeMax());

    VolumeAtomic<bool> markedVoxels(dims);
    markedVoxels.clear();

    std::vector<ivec3> neighbors;
    for(int z=-1; z <= 1; z++)
        for(int y=-1; y <= 1; y++)
            for(int x=-1; x <= 1; x++) {
                ivec3 p(x, y, z);
                if(p != ivec3::zero)
                    neighbors.push_back(p);
            }


    std::queue<QueueVoxel> voxelQueue;
    voxelQueue.push(QueueVoxel(seed_pos, 0));
    /*
    for(size_t i=0; i < neighbors.size(); i++)
        voxelQueue.push(seed_pos + neighbors[i]);
        */

    Stats stats_value;
    Stats stats_gradmag;

    float seed_value = dataset->voxel(seed_pos) / max;
    float seed_stddev26 = neighborStandardDeviation(seed_pos, dataset, max, stats_value);

    float seed_gradmag = gradientMagnitude(seed_pos, dataset, max);
    float seed_gradmag_stddev26 = neighborStandardDeviationGradients(seed_pos, dataset, max, stats_gradmag);

    if(mode == FLOODFILL_INTENSITY || mode == FLOODFILL_WEIGHTED)
        std::cout << "seed value: " << seed_value << ", "
                  << "stddev 26: " << seed_stddev26 << std::endl;

    if(mode == FLOODFILL_GRADMAG || mode == FLOODFILL_WEIGHTED)
        std::cout << "seed_gradmag: " << seed_gradmag << ", "
                  << "stddev gradmag 26: " << seed_gradmag_stddev26 << std::endl;

    while(!voxelQueue.empty()) {
        QueueVoxel voxel = voxelQueue.front();
        ivec3 pos = voxel.pos_;
        uint32_t step = voxel.step_;
        voxelQueue.pop();

        if(pos.x < 2 || pos.x > dims.x - 3 ||
            pos.y < 2 || pos.y > dims.y - 3 ||
            pos.z < 2 || pos.z > dims.z - 3 ||
            markedVoxels.voxel(pos))
        {
            // on border or already visited
            continue;
        }

        float value = dataset->voxel(pos) / max;

        if(adaptive /* && not in initial neighbor */) {
            stats_value.add(value);
            seed_stddev26 = stats_value.stddev();

            if(mode == FLOODFILL_GRADMAG || mode == FLOODFILL_WEIGHTED) {
                stats_gradmag.add(gradientMagnitude(pos, dataset, max));
                seed_gradmag_stddev26 = stats_gradmag.stddev();
            }
        }


        if(useThresholds && (!crit1.isFullfilledAt(pos) || !crit2.isFullfilledAt(pos))) {
            // invalid value
            continue;
        }

        // Cost function: if less than 1 then voxel is within the region.
        //
        // Based on: Runzhen Huang, Kwan-Liu Ma. RGVis: Region growing based techniques for
        // volume visualization, 2003.
        float cost = 0.f;

        if(mode == FLOODFILL_INTENSITY)
            cost = fabs(value - seed_value) / (strictness * seed_stddev26);
        else if(mode == FLOODFILL_GRADMAG) {
            float gradmag = gradientMagnitude(pos, dataset, max);
            cost = fabs(gradmag - seed_gradmag) / (strictness * seed_gradmag_stddev26);
        }
        else if(mode == FLOODFILL_WEIGHTED) {
            float cost_a = fabs(value - seed_value) / (strictness * seed_stddev26);

            float gradmag = gradientMagnitude(pos, dataset, max);
            float cost_b = fabs(gradmag - seed_gradmag) / (strictness * seed_gradmag_stddev26);

            // weight p
            float p = (seed_gradmag_stddev26 / (seed_stddev26 + seed_gradmag_stddev26));
            cost = cost_a * p + cost_b * (1.f - p);
        }

        if(cost >= 1.f)
            continue;

        if(maxSteps > 0 && step >= maxSteps)
            continue;

        // voxel is valid
        markedVoxels.voxel(pos) = true;

        // add neighbors to stack if not already visited
        for(size_t i=0; i < neighbors.size(); i++) {
            tgt::ivec3 npos = pos + neighbors[i];
            if(!border.isCrossedBy(pos, npos) && !markedVoxels.voxel(pos + neighbors[i])) {
                voxelQueue.push(QueueVoxel(npos, step+1));
            }
        }
    }

    // now fill segmentation volume with all marked voxels
    float count = 0;
    for(size_t z=0; z < markedVoxels.getDimensions().z; z++) {
        for(size_t y=0; y < markedVoxels.getDimensions().y; y++) {
            for(size_t x=0; x < markedVoxels.getDimensions().x; x++) {
                if(markedVoxels.voxel(x, y, z)) {
                    uint8_t& v = segvol.voxel(x, y, z);
                    if(v == 0 || segment == 0) {
                        v = segment;
                        count++;
                    }
                }
            }
        }
    }

    return static_cast<int>(count);
}

} // anonymous namespace


void AortaSegmentation::mark(const ivec3& seedpos, int segment, const HeartSegmentationBorder& border, VolumeRAM_UInt8& segvol, const VolumeBase& volumeInput, const VolumeBase& vesselnessInput) {

    const VolumeRAM* vol = volumeInput.getRepresentation<VolumeRAM>();

    if(!vol) {
        LERROR("no volume or segmentation");
        return;
    }

    float seedval = vol->getVoxelNormalized(seedpos);
    LINFO("seed pos " << seedpos << " with value " << volumeInput.getRealWorldMapping().normalizedToRealWorld(seedval));

    if(seedval <= 0.f) {
        LERROR("ignoring this seed value");
        return;
    }

    // start the flood fill
    //VolumeRAM_UInt8* s = dynamic_cast<VolumeRAM_UInt8*>(segmentation_->getRepresentation<VolumeRAM>());
    FloodFillMode mode = FLOODFILL_INTENSITY;
    if(fillCostFunction_.get() == "gradient-magnitude")
        mode = FLOODFILL_GRADMAG;
    else if(fillCostFunction_.get() == "weighted")
        mode = FLOODFILL_WEIGHTED;

    ThresholdCriterion volCriterion(volumeInput, volumeThresholds_.get().x, volumeThresholds_.get().y);
    ThresholdCriterion vesselnessCriterion(vesselnessInput, vesselnessThresholds_.get().x, vesselnessThresholds_.get().y);

    int count = 0;
    if(const VolumeRAM_UInt8* volu8 = dynamic_cast<const VolumeRAM_UInt8*>(vol))
        count = floodFill(seedpos, segment, volCriterion, vesselnessCriterion, strictness_.get(),
                          volu8, segvol, mode, thresholdFilling_.get(), adaptive_.get(),
                          maxSteps_.get(), border);
    else if(const VolumeRAM_UInt16* volu16 = dynamic_cast<const VolumeRAM_UInt16*>(vol))
        count = floodFill(seedpos, segment, volCriterion, vesselnessCriterion, strictness_.get(),
                          volu16, segvol, mode, thresholdFilling_.get(), adaptive_.get(),
                          maxSteps_.get(), border);
    else if(const VolumeRAM_Float* volf = dynamic_cast<const VolumeRAM_Float*>(vol))
        count = floodFill(seedpos, segment, volCriterion, vesselnessCriterion, strictness_.get(),
                          volf, segvol, mode, thresholdFilling_.get(), adaptive_.get(),
                          maxSteps_.get(), border);

    LINFO("filled voxels: " << count);
}

static bool tryExtractPoints(const Geometry* geometry, std::vector<tgt::vec3>& output) {
    if(const PointSegmentListGeometryVec3* seedList = dynamic_cast<const PointSegmentListGeometryVec3*>(geometry)) {
        for(int i = 0; i < seedList->getNumSegments(); i++) {
            auto segment = seedList->getSegment(i);
            output.insert(output.end(), segment.begin(), segment.end());
        }
    } else if(const PointListGeometryVec3* seeds = dynamic_cast<const PointListGeometryVec3*>(geometry)) {
        output.insert(output.end(), seeds->begin(), seeds->end());
    } else {
        return false;
    }
    return true;
}

static bool tryExtractCenterPoint(const GeometryPort& g, tgt::vec3& output) {
    if(!g.hasData()) {
        return false;
    }
    std::vector<tgt::vec3> points;
    if(!tryExtractPoints(g.getData(), points)) {
        return false;
    }
    if(points.empty()) {
        return false;
    }
    output = std::accumulate(points.begin(), points.end(), tgt::vec3::zero) / static_cast<float>(points.size());
    return true;
}

void AortaSegmentation::process() {
    if(!updateContinuously_.get() && !updateForced_) {
        return;
    }
    updateForced_ = false;

    if(!volumeInport_.hasData() || !vesselnessInport_.hasData()) {
        LERROR("Missing input volume(s).");
        return;
    }
    const VolumeBase& volumeInput = *volumeInport_.getData();
    const VolumeBase& vesselnessInput = *vesselnessInport_.getData();
    if(volumeInput.getDimensions() != vesselnessInput.getDimensions()) {
        LERROR("Input volume dimension mismatch.");
        return;
    }

    tgt::vec3 avgSeedPoint;
    if(!tryExtractCenterPoint(seedPoints_, avgSeedPoint)) {
        LERROR("No seed points");
        return;
    }
    tgt::vec3 avgHeartPoint;
    if(!tryExtractCenterPoint(heartPoints_, avgHeartPoint)) {
        LERROR("No hearth points");
        return;
    }

    HeartSegmentationBorder border(avgSeedPoint, avgHeartPoint, separatingDiskRadius_.get());

    std::unique_ptr<VolumeRAM_UInt8> segVol(new VolumeRAM_UInt8(volumeInput.getDimensions()));
    int segmentID = 128;
    ivec3 voxelSeed = tgt::round(avgSeedPoint);
    mark(voxelSeed, segmentID, border, *segVol, volumeInput, vesselnessInput);

    outport_.setData(new Volume(segVol.release(), &volumeInput));
}

void AortaSegmentation::adjustPropertiesToInput() {
    if(volumeInport_.hasData()) {
        const VolumeMinMax* mm = volumeInport_.getData()->getDerivedData<VolumeMinMax>();
        tgtAssert(mm, "No VolumeMinMax");

        volumeThresholds_.setMinValue(mm->getMin());
        volumeThresholds_.setMaxValue(mm->getMax());

        maxSteps_.setMaxValue(3*tgt::hadd(volumeInport_.getData()->getDimensions()));
    }
    if(vesselnessInport_.hasData()) {
        const VolumeMinMax* mm = vesselnessInport_.getData()->getDerivedData<VolumeMinMax>();
        tgtAssert(mm, "No VolumeMinMax");

        vesselnessThresholds_.setMinValue(mm->getMin());
        vesselnessThresholds_.setMaxValue(mm->getMax());
    }

    tgt::vec3 avgSeedPoint;
    tgt::vec3 avgHeartPoint;
    if((seedPoints_.hasChanged() || heartPoints_.hasChanged())
            && tryExtractCenterPoint(seedPoints_, avgSeedPoint)
            && tryExtractCenterPoint(heartPoints_, avgHeartPoint)) {
       float halfDist = tgt::distance(avgSeedPoint, avgHeartPoint)/2;
       separatingDiskRadius_.setMaxValue(5*halfDist);
       //separatingDiskRadius_.set(halfDist);
    }
}

} // namespace
