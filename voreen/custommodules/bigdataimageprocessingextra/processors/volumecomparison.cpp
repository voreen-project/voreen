/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "volumecomparison.h"

#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/utils/stringutils.h"

#include "tgt/stopwatch.h"

#include <chrono>

namespace voreen {

const std::string VolumeComparison::loggerCat_("voreen.VolumeComparison");

VolumeComparison::VolumeComparison()
    : Processor()
    , firstSegmentationVolume_(Port::INPORT, "firstsegmentation", "First Segmentation Volume", false)
    , secondSegmentationVolume_(Port::INPORT, "secondsegmentation", "Second Segmentation Volume", false)
    , enabled_("enabled", "Enabled", true)
    , useClipRegion_("useClipRegion", "Use Clip Region", false)
    , clipRegion_("clipRegion", "Clip Region", tgt::IntBounds(tgt::ivec3(0), tgt::ivec3(1)), tgt::ivec3(0), tgt::ivec3(1))
    , binarizationThreshold_("binarizationThreshold", "Binarization Threshold", 0.5, 0.0, 1.0)
    , progressProperty_("progressProperty", "Quantification Progress")
    , csvSaveFile_("csvFileProp", "CSV Export Path", "CSV Export Path", ".", "Comma seperated values (*.csv)", FileDialogProperty::SAVE_FILE)
{
    addPort(firstSegmentationVolume_);
    addPort(secondSegmentationVolume_);

    addProperty(enabled_);

    addProperty(useClipRegion_);
    clipRegion_.setVisibleFlag(false);
    addProperty(clipRegion_);

    addProperty(progressProperty_);
    addProgressBar(&progressProperty_);

    addProperty(csvSaveFile_);

    ON_CHANGE(useClipRegion_, VolumeComparison, useClipRegionChanged);
}

Processor* VolumeComparison::create() const {
    return new VolumeComparison();
}

void VolumeComparison::initialize() {
    Processor::initialize();

    adjustPropertiesToInput();
}

template<typename V1>
struct QuantificationImpl {
template<typename V2>
static void quantificationImpl(const VolumeAtomic<V1>& slice1, const VolumeRAM& dynamicSlice, VolumeComparison::ScanSummary& globalSummary, tgt::svec3 llf, tgt::svec3 urb, RealWorldMapping rwm, float rwThreadhold) {
    const VolumeAtomic<V2>& slice2 = dynamic_cast<const VolumeAtomic<V2>&>(dynamicSlice);

    float threshold = rwm.realWorldToNormalized(rwThreadhold);

    VolumeComparison::ScanSummary summary;
    for (size_t y = llf.y; y <= urb.y; ++y) {
        size_t offset = y * (urb.x-llf.x);
        size_t xlo = offset + llf.x;
        size_t xhi = offset + urb.x+1;
        for (size_t x = xlo; x < xhi; ++x) {
            V1 v1g = slice1.voxel()[x];
            float v1 = getTypeAsFloat(v1g);

            V2 v2g = slice2.voxel()[x];
            float v2 = getTypeAsFloat(v2g);

            bool foregroundOne = v1 > threshold;
            bool foregroundTwo = v2 > threshold;

            if (foregroundOne && foregroundTwo)
                summary.numForegroundBoth_++;
            if (!foregroundOne && foregroundTwo)
                summary.numForegroundOnlyTwo_++;
            if (foregroundOne && !foregroundTwo)
                summary.numForegroundOnlyOne_++;
            if (!foregroundOne && !foregroundTwo)
                summary.numBackgroundBoth_++;

            float diff = v1 - v2;
            summary.sumOfVoxelDiffsSquared_ += diff*diff;
            summary.sumOfVoxelDiffsAbs_ += std::abs(diff);
        }
    }
    globalSummary.merge(summary);
}
};

template<typename V1>
static void quantificationIntermediate(const VolumeRAM& dynamicSlice, const VolumeRAM& slice2, VolumeComparison::ScanSummary& summary, tgt::svec3 llf, tgt::svec3 urb, RealWorldMapping rwm, float rwThreadhold) {
    const VolumeAtomic<V1>& slice1 = dynamic_cast<const VolumeAtomic<V1>&>(dynamicSlice);
    DISPATCH_FOR_BASETYPE(slice2.getBaseType(), QuantificationImpl<V1>::template quantificationImpl, slice1, slice2, summary, llf, urb, rwm, rwThreadhold);
}

static void quantification(const VolumeRAM& slice1, const VolumeRAM& slice2, VolumeComparison::ScanSummary& summary, tgt::svec3 llf, tgt::svec3 urb, RealWorldMapping rwm, float rwThreadhold) {
    DISPATCH_FOR_BASETYPE(slice1.getBaseType(), quantificationIntermediate, slice1, slice2, summary, llf, urb, rwm, rwThreadhold);
}

void VolumeComparison::process() {
    if(!enabled_.get()) {
        return;
    }

    setProgress(0.f);
    //quantificationPlot_.clear();

    // get volumes
    const VolumeBase* volume1Ptr = firstSegmentationVolume_.getData();
    const VolumeBase* volume2Ptr = secondSegmentationVolume_.getData();

    if(!volume1Ptr || !volume2Ptr) {
        return;
    }

    const VolumeBase& volume1 = *volume1Ptr;
    const VolumeBase& volume2 = *volume2Ptr;

    // check number of channels, dimensions, and data type of the volumes
    if (volume1.getNumChannels() != 1) {
        LERROR("First volume has more than one channel!");
        return;
    }
    if (volume2.getNumChannels() != 1) {
        LERROR("Second volume has more than one channel!");
        return;
    }

    tgt::svec3 dimensions = volume1.getDimensions();
    if (dimensions != volume2.getDimensions()) {
        LERROR("Volumes are not of the same size!");
        return;
    }

    RealWorldMapping rwm = volume1.getRealWorldMapping();
    if (rwm != volume2.getRealWorldMapping()) {
        LERROR("Real world mappings differ!");
        return;
    }

    ScanSummary summary{};

    // get the quantification dimensions
    tgt::svec3 llf, urb;
    llf = tgt::svec3::zero;

    urb = dimensions - tgt::svec3::one;

    // if the clip region is used, crop our bounds
    if (useClipRegion_.get()) {
        llf = tgt::max(llf, tgt::svec3(clipRegion_.get().getLLF()));
        urb = tgt::min(urb, tgt::svec3(clipRegion_.get().getURB()));
    }

    typedef std::chrono::steady_clock Clock;
    typedef std::chrono::time_point<Clock> TimePoint;

    TimePoint start = Clock::now();

#ifdef VRN_MODULE_OPENMP
#pragma omp parallel
#endif
    {
        VolumeComparison::ScanSummary threadSummary;

#ifdef VRN_MODULE_OPENMP
#pragma omp for
#endif
        for (size_t z = llf.z; z <= urb.z; ++z) {

            // we do not know how large a single slice is, so we only load one slice at a time for each volume
            std::unique_ptr<VolumeRAM> slice1(volume1.getSlice(z));
            std::unique_ptr<VolumeRAM> slice2(volume2.getSlice(z));
            tgtAssert(slice1, "No slice 1");
            tgtAssert(slice2, "No slice 2");

            // quantify slice according to the data type
            quantification(*slice1.get(), *slice2.get(), threadSummary, llf, urb, rwm, binarizationThreshold_.get());

            //setProgress(std::min(0.99f, static_cast<float>(z - llf.z) / static_cast<float>(urb.z - llf.z)));
        }

#ifdef VRN_MODULE_OPENMP
#pragma omp critical
#endif
        {
            summary.merge(threadSummary);
        }
    }

    setProgress(1.f);
    TimePoint end = Clock::now();

    auto time = end - start;
    auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(time);

    LINFO("Comparison required " << std::to_string(static_cast<float>(millis.count())/1000.0) << " seconds");

    std::stringstream header, line;
    size_t numVoxels = summary.totalNumberOfVoxels();
    float avgDiffAbs = summary.sumOfVoxelDiffsAbs_ / numVoxels;
    float variance = summary.sumOfVoxelDiffsSquared_ / (numVoxels * numVoxels);
    header << "Number_of_voxels,voxels_in_volume1,voxels_in_volume2,voxels_in_both,dice_score,sumDiffAbs,avgDiffAbs,sumDiffSquared,variance" << std::endl;
    line << numVoxels
        << "," << summary.numForegroundOnlyOne_
        << "," << summary.numForegroundOnlyTwo_
        << "," << summary.numForegroundBoth_
        << "," << summary.diceScore()
        << "," << summary.sumOfVoxelDiffsAbs_
        << "," << avgDiffAbs
        << "," << summary.sumOfVoxelDiffsSquared_
        << "," << variance;

    LINFO("Header: " << header.str());
    LINFO("Result: " << line.str());

    std::ofstream file(csvSaveFile_.get());
    file << header.str();
    file << line.str();
    file.close();
}

void VolumeComparison::adjustPropertiesToInput() {
    if (!isInitialized())
        return;

    // adjust clipping area to input
    const VolumeBase* firstVol = firstSegmentationVolume_.getData();
    const VolumeBase* secondVol = secondSegmentationVolume_.getData();
    const VolumeBase* someVol = firstVol ? firstVol : secondVol;
    if (firstVol && secondVol) {
        // check if dimensions are the same and adjust clipping area
        if (firstVol->getDimensions() == secondVol->getDimensions())
            clipRegion_.setMaxValue(tgt::ivec3(firstVol->getDimensions() - tgt::svec3::one));
    }
    else if (firstVol)
        clipRegion_.setMaxValue(tgt::ivec3(firstVol->getDimensions() - tgt::svec3::one));
    else if (secondVol)
        clipRegion_.setMaxValue(tgt::ivec3(secondVol->getDimensions() - tgt::svec3::one));

    if(someVol) {
        auto vmm = someVol->getDerivedData<VolumeMinMax>();
        binarizationThreshold_.setMinValue(vmm->getMin());
        binarizationThreshold_.setMaxValue(vmm->getMax());
    }
}

void VolumeComparison::useClipRegionChanged() {
    clipRegion_.setVisibleFlag(useClipRegion_.get());
}

VolumeComparison::ScanSummary::ScanSummary()
    : numForegroundBoth_(0)
    , numForegroundOnlyOne_(0)
    , numForegroundOnlyTwo_(0)
    , numBackgroundBoth_(0)
    , sumOfVoxelDiffsSquared_(0.0f)
    , sumOfVoxelDiffsAbs_(0.0f)
{
}

size_t VolumeComparison::ScanSummary::totalNumberOfVoxels() const {
    return numForegroundBoth_ + numForegroundOnlyOne_ + numForegroundOnlyTwo_ + numBackgroundBoth_;
}
float VolumeComparison::ScanSummary::diceScore() const {
    return 2.0f * numForegroundBoth_ / (numForegroundOnlyOne_ + numForegroundOnlyTwo_ + 2*numForegroundBoth_);
}
void VolumeComparison::ScanSummary::merge(ScanSummary& other) {
    numForegroundBoth_ += other.numForegroundBoth_;
    numForegroundOnlyOne_ += other.numForegroundOnlyOne_;
    numForegroundOnlyTwo_ += other.numForegroundOnlyTwo_;
    numBackgroundBoth_ += other.numBackgroundBoth_;
    sumOfVoxelDiffsSquared_ += other.sumOfVoxelDiffsSquared_;
    sumOfVoxelDiffsAbs_ += other.sumOfVoxelDiffsAbs_;
}

} // namespace
