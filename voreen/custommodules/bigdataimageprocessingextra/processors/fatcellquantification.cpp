/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "fatcellquantification.h"

#include "voreen/core/datastructures/volume/volumedisk.h"
#include "modules/bigdataimageprocessing/util/csvwriter.h"

#include <vector>

namespace voreen {

const std::string FatCellQuantification::loggerCat_("voreen.fatcellquantification");

FatCellQuantification::FatCellQuantification()
    : Processor()
    , firstSegmentationVolume_(Port::INPORT, "firstsegmentation", "First Segmentation Volume", false)
    //, secondSegmentationVolume_(Port::INPORT, "secondsegmentation", "Second (Reference) Segmentation Volume", false)
    //, useClipRegion_("useClipRegion", "Use Clip Region", false)
    //, clipRegion_("clipRegion", "Clip Region", tgt::IntBounds(tgt::ivec3(0), tgt::ivec3(1)), tgt::ivec3(0), tgt::ivec3(1))
    , csvSaveFile_("csvFileProp", "CSV Export Path", "CSV Export Path", ".", "Comma seperated values (*.csv)", FileDialogProperty::SAVE_FILE)
    //, saveToCsv_("savetocsv", "Save to CSV")
    //, quantificationPlot_(Port::OUTPORT, "QuantificationPlot", "Quantification Plot", true,  Processor::VALID)
{
    addPort(firstSegmentationVolume_);

    ON_CHANGE(firstSegmentationVolume_, FatCellQuantification, adjustToInputVolumes);

    addProperty(csvSaveFile_);
}

Processor* FatCellQuantification::create() const {
    return new FatCellQuantification();
}

bool FatCellQuantification::isReady() const {
    return isInitialized() && (firstSegmentationVolume_.isReady());
}

struct FatCellStats {
    FatCellStats()
        : id_(0)
        , numVoxels_(0)
        , voxelPosSum_(0.0)
        , centroid_(0.0)
        , surfaceDistanceSum_(0.0)
        , numSurfaceVoxels_(0)
        , avgDistToSurface_(0.0)
        , distSquareSum_(0.0)
        , distToSurfaceStddev_(0.0)
    {
    }
    uint32_t id_;
    uint32_t numVoxels_;
    tgt::vec3 voxelPosSum_;

    tgt::vec3 centroid_;

    float surfaceDistanceSum_;
    uint32_t numSurfaceVoxels_;

    float avgDistToSurface_;

    float distSquareSum_;

    float distToSurfaceStddev_;
};

void FatCellQuantification::process() {
    const std::string csvPath = csvSaveFile_.get();
    if (csvPath.empty()) {
        LERROR("No csv output path specified.");
        return;
    }
    const VolumeBase* volume1 = firstSegmentationVolume_.getData();
    if (!volume1) {
        LERROR("No input volume!");
        return;
    }
    const VolumeBase& volume = *volume1;
    if(volume.getFormat() != "uint32") {
        LERROR("Expected uint32 volume");
        return;
    }
    const VolumeRAM* volumeram = volume.getRepresentation<VolumeRAM>();
    tgtAssert(volumeram, "No volumeram");
    const VolumeAtomic<uint32_t>& vol = *dynamic_cast<const VolumeAtomic<uint32_t>*>(volumeram);

    const bool FOREGROUND = true;
    const bool BACKGROUND = false;
    auto dim = volume.getDimensions();
    tgt::ivec3 idim = dim;
    auto getVoxel = [&] (tgt::ivec3 p, bool outsideVolumeValue) {
        if (0 > p.x || p.x >= idim.x ||
            0 > p.y || p.y >= idim.y ||
            0 > p.z || p.z >= idim.z) {
            return outsideVolumeValue;
        }
        if(vol.voxel(p.x, p.y, p.z) > 0){
            return FOREGROUND;
        } else {
            return BACKGROUND;
        }
    };

    auto isSurfaceVoxel = [&] (tgt::ivec3 p) {
        return getVoxel(p, FOREGROUND) == FOREGROUND && (
               getVoxel(tgt::ivec3(p.x+1, p.y  , p.z  ), FOREGROUND) == BACKGROUND
            || getVoxel(tgt::ivec3(p.x-1, p.y  , p.z  ), FOREGROUND) == BACKGROUND
            || getVoxel(tgt::ivec3(p.x  , p.y+1, p.z  ), FOREGROUND) == BACKGROUND
            || getVoxel(tgt::ivec3(p.x  , p.y-1, p.z  ), FOREGROUND) == BACKGROUND
            || getVoxel(tgt::ivec3(p.x  , p.y  , p.z+1), FOREGROUND) == BACKGROUND
            || getVoxel(tgt::ivec3(p.x  , p.y  , p.z-1), FOREGROUND) == BACKGROUND);
    };

    auto transform = volume.getVoxelToWorldMatrix();
    std::vector<FatCellStats> statvec;
    for(size_t z = 0; z<dim.z; ++z) {
        for(size_t y = 0; y<dim.y; ++y) {
            for(size_t x = 0; x<dim.x; ++x) {
                uint32_t label = vol.voxel(x, y, z);
                if(label == 0) {
                    continue;
                }
                if(label >= statvec.size()) {
                    statvec.resize(label+1, FatCellStats());
                }

                tgt::vec3 pos = transform.transform(tgt::vec3(x,y,z));

                auto& stats = statvec.at(label-1);
                stats.id_ = label;
                stats.numVoxels_ += 1;
                stats.voxelPosSum_ += pos;
            }
        }
    }
    for(auto& stats : statvec) {
        if(stats.id_ != 0) {
            stats.centroid_ = stats.voxelPosSum_ / static_cast<float>(stats.numVoxels_);
        }
    }
    for(size_t z = 0; z<dim.z; ++z) {
        for(size_t y = 0; y<dim.y; ++y) {
            for(size_t x = 0; x<dim.x; ++x) {
                uint32_t label = vol.voxel(x, y, z);
                if(label == 0) {
                    continue;
                }
                auto& stats = statvec.at(label-1);

                tgt::ivec3 ipos(x,y,z);

                if(isSurfaceVoxel(ipos)) {
                    tgt::vec3 pos = transform.transform(tgt::vec3(ipos));
                    stats.surfaceDistanceSum_ += tgt::distance(pos, stats.centroid_);
                    stats.numSurfaceVoxels_ += 1;
                }
            }
        }
    }
    for(auto& stats : statvec) {
        if(stats.id_ != 0) {
            tgtAssert(stats.numSurfaceVoxels_ >= 1, "invalid number of surface voxels");
            stats.avgDistToSurface_ = stats.surfaceDistanceSum_ / stats.numSurfaceVoxels_;
        }
    }
    for(size_t z = 0; z<dim.z; ++z) {
        for(size_t y = 0; y<dim.y; ++y) {
            for(size_t x = 0; x<dim.x; ++x) {
                uint32_t label = vol.voxel(x, y, z);
                if(label == 0) {
                    continue;
                }
                auto& stats = statvec.at(label-1);

                tgt::ivec3 ipos(x,y,z);

                if(isSurfaceVoxel(ipos)) {
                    tgt::vec3 pos = transform.transform(tgt::vec3(ipos));
                    float distDiff = std::abs(tgt::distance(pos, stats.centroid_) - stats.avgDistToSurface_);
                    stats.distSquareSum_ += distDiff * distDiff;
                }
            }
        }
    }
    for(auto& stats : statvec) {
        if(stats.id_ != 0) {
            tgtAssert(stats.numSurfaceVoxels_ >= 2, "Invalid number of surface voxels (at least 2)");
            float variance = stats.distSquareSum_/(stats.numSurfaceVoxels_-1);
            stats.distToSurfaceStddev_ = std::sqrt(variance);
        }
    }

    CSVWriter<uint32_t, uint32_t, float, float, float, float, float, float> writer(csvPath);
    writer.writeHeader("label", "count", "volume", "center_x", "center_y", "center_z", "surface_dist_avg", "surface_dist_std_dev");
    float voxelVolume = tgt::hmul(volume.getSpacing());
    for(auto& stats : statvec) {
        if(stats.id_ != 0) {
            uint32_t count = stats.numVoxels_;
            float volume = count * voxelVolume;
            tgt::vec3 centroid = stats.centroid_;
            writer.write(stats.id_, count, volume, centroid.x, centroid.y, centroid.z, stats.avgDistToSurface_, stats.distToSurfaceStddev_);
        }
    }
}

void FatCellQuantification::adjustToInputVolumes() {
}

/*void SegmentationQuantification::exportToCSV() {
    std::ofstream file(csvSaveFile_.get());
    file << "Number_of_voxels,voxels_in_volume1,voxels_in_volume2,voxels_in_both" << std::endl;
    file << numVoxelsTotal_ << "," << numVoxelsInOne_ << "," << numVoxelsInTwo_ << "," << numVoxelsInBoth_ << std::endl;

    file.close();
}*/

} // namespace
