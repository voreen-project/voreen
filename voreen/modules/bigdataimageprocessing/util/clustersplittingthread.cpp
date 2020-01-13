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

#include "clustersplittingthread.h"

#include "voreen/core/datastructures/volume/volume.h"

#include "../operators/volumeoperatordistancetransform.h"
#include "../operators/volumeoperatorfastvolumecombine.h"
#include "../operators/volumeoperatorgradientdescent.h"
#include "../operators/volumeoperatorwatershed.h"

#include "voreen/core/datastructures/volume/operators/volumeoperatorgaussian.h"

#include "voreen/core/datastructures/geometry/pointlistgeometry.h"

namespace voreen {

ClusterSplittingThread::ClusterSplittingThread(NucleiClusterSplitting* processor, ConnectedComponentQueue& queue, boost::mutex& globalValueMutex, bool smoothForSeedDetection, size_t gaussianKernelSize, float gaussianSigma, bool maskBeforeSmoothing, bool expandMarkers, size_t minSeedRank)
    : BackgroundThread()
    , processor_(processor)
    , queue_(queue)
    , globalValueMutex_(globalValueMutex)
    , smoothForSeedDetection_(smoothForSeedDetection)
    , gaussianKernelSize_(gaussianKernelSize)
    , gaussianSigma_(gaussianSigma)
    , maskBeforeSmoothing_(maskBeforeSmoothing)
    , expandMarkers_(expandMarkers)
    , minSeedRank_(minSeedRank)
{ }

void ClusterSplittingThread::clearLocalData() {
    // clear list
    std::vector< std::vector<tgt::vec3> > emptyList;
    localResult_.setData(emptyList);
}

void ClusterSplittingThread::pushLocalData() {
    boost::mutex::scoped_lock lock(globalValueMutex_);

    // iterate over segments and add them to the processor's result
    for(auto& segment : localResult_.getData())
        processor_->nucleiSegmentList_.addSegment(segment);
}

void ClusterSplittingThread::handleInterruption() {
    pushLocalData();
}


void ClusterSplittingThread::threadMain() {

    clearLocalData();

    // we use this while-loop because at some point we will wait for the empty queue and get the interruptAndJoin signal, which will cause this thread to terminate as it is the only interruption point
    while(true) {

        // if the queue is empty the thread will be paused, which also provides our interruption point for termination
        ConnectedComponentQueue::ConnectedComponent component = queue_.pop();
        
        // create a foreground mask using the component id
        tgt::svec3 dim = component.image_->getDimensions();
        VolumeRAM_UInt16* foregroundRAM = component.image_->clone(); 
        for (size_t i = 0; i < foregroundRAM->getNumVoxels(); ++i) {
            if (static_cast<size_t>(component.labels_->voxel(i)) == component.id_)
                foregroundRAM->voxel(i) = 1;
            else
                foregroundRAM->voxel(i) = 0;
        } 

        // label-volume is not needed anymore -> free memory
        delete component.labels_;
        component.labels_ = 0;

        // create Volume objects for volume operators
        Volume* image = new Volume(component.image_, component.spacing_, component.offset_);
        Volume* foreground = new Volume(foregroundRAM, component.spacing_, component.offset_);

        // create masked volume
        Volume* maskedVolume = 0;
        if(smoothForSeedDetection_) {
            if (maskBeforeSmoothing_) {
                // first apply the masking, then smooth (border lies near background)
                Volume* tmpMaskImage = VolumeOperatorFastVolumeCombine::get(image, foreground)->apply(image, foreground, OP_MASK_A_BY_B, 0, 0, 1, 0, 0, 0);
                maskedVolume = VolumeOperatorGaussian::APPLY_OP(tmpMaskImage, gaussianKernelSize_, gaussianSigma_, 0);
                delete tmpMaskImage;
            }
            else {
                // first smooth, then apply masking
                Volume* filteredImage = VolumeOperatorGaussian::APPLY_OP(image, gaussianKernelSize_, gaussianSigma_, 0);
                maskedVolume = VolumeOperatorFastVolumeCombine::get(filteredImage, foreground)->apply(filteredImage, foreground, OP_MASK_A_BY_B, 0, 0, 1, 0, 0, 0);
                delete filteredImage;
            }
        }
        else {
            maskedVolume = VolumeOperatorFastVolumeCombine::get(image, foreground)->apply(image, foreground, OP_MASK_A_BY_B, 0, 0, 1, 0, 0, 0);
        }
        
        // original image not needed anymore
        delete image;
        
        Volume* detectedMarkers = VolumeOperatorGradientDescent::APPLY_OP(maskedVolume, 0);
        VolumeRAM_UInt16* markers = const_cast<VolumeRAM_UInt16*>(reinterpret_cast<const VolumeRAM_UInt16*>(detectedMarkers->getRepresentation<VolumeRAM>()));
        
        //collect seeds
        std::vector<SeedDescriptor> seeds;

        VRN_FOR_EACH_VOXEL(idx, tgt::svec3::zero, dim){
            uint16_t cellValue = markers->voxel(idx);
            if(cellValue != 0){
                SeedDescriptor se = {idx, FLT_MAX, cellValue};
                seeds.push_back(se);
            }
        }

        //relabel seeds and delete too small seeds
        uint16_t label = 1;
        PointListGeometryVec3* plg = new PointListGeometryVec3();
        for(size_t j = 0; j < seeds.size(); ++j) {
            if(seeds[j].identifier >= (uint16_t) minSeedRank_){
                seeds[j].identifier = label;
                markers->voxel(seeds[j].CenterOfMaxima) = label;
                plg->addPoint(tgt::vec3(seeds[j].CenterOfMaxima));
                label += 1;
            }else{
                markers->voxel(seeds[j].CenterOfMaxima) = 0;
                seeds.erase(seeds.begin() + j);
                j--;
            }
        }
        
        //Expand Markers!
        if(expandMarkers_){
            Volume* distanceTransform = VolumeOperatorEuclideanDistanceTransform::APPLY_OP(foreground, 0);
            const VolumeRAM_Float* dt = reinterpret_cast<const VolumeRAM_Float*>(distanceTransform->getRepresentation<VolumeRAM>());
            if(seeds.size() > 1) {
                for(size_t j = 0; j < seeds.size(); ++j){
                    float ClosestDistance = FLT_MAX;
                    for(size_t k = 0; k < seeds.size(); ++k){
                        if(k != j){
                            ClosestDistance = std::min(ClosestDistance, tgt::distance(tgt::vec3(seeds[j].CenterOfMaxima), tgt::vec3(seeds[k].CenterOfMaxima)));
                        }
                    }
                    seeds[j].DistanceToClosestMaximum = ClosestDistance;
                }

                uint16_t ColorIndex = 1;
                for(size_t j = 0; j < seeds.size(); ++j) {
                    float distanceToBackground = dt->voxel(seeds[j].CenterOfMaxima);
                    float actualDistance = floorf(distanceToBackground) == distanceToBackground ? (distanceToBackground - 1.0f) : floorf(distanceToBackground);
                    float minHeight = FLT_MAX;
                    int minRadius = static_cast<int>(std::min(floorf(0.5f*(seeds[j].DistanceToClosestMaximum) - 1), actualDistance));
                    for (tgt::ivec3 sub_idx = tgt::ivec3(-minRadius); sub_idx.z <= minRadius; ++sub_idx.z){
                        for (sub_idx.y = -minRadius; sub_idx.y <= minRadius; ++sub_idx.y){
                            for (sub_idx.x = -minRadius; sub_idx.x <= minRadius; ++sub_idx.x){
                                long long x = seeds[j].CenterOfMaxima.x + sub_idx.x;
                                long long y = seeds[j].CenterOfMaxima.y + sub_idx.y;
                                long long z = seeds[j].CenterOfMaxima.z + sub_idx.z;

                                tgt::svec3 index = tgt::clamp(tgt::svec3((size_t)x,(size_t)y,(size_t)z), tgt::svec3::zero, dim - tgt::svec3::one);
                                if(tgt::length(tgt::vec3(sub_idx)) <= minRadius){
                                    markers->voxel(index) = seeds[j].identifier;
                                }
                            }
                        }
                    }
                }
            }
            delete distanceTransform;
        }

        //segment!
        Volume* segmentedImage = VolumeOperatorWatershedTransform::get(maskedVolume, detectedMarkers, foreground)->apply(maskedVolume, detectedMarkers, foreground);

        // delete volume data to free memory
        delete maskedVolume;
        delete detectedMarkers;
        delete foreground;

        // for each label: get its voxels
        const VolumeRAM_UInt16* segmentationRAM = reinterpret_cast<const VolumeRAM_UInt16*>(segmentedImage->getRepresentation<VolumeRAM>());
        std::map<uint16_t, std::vector<tgt::svec3> > segments = std::map<uint16_t, std::vector<tgt::svec3> >();
        std::vector<tgt::vec3> centroids;
        for (size_t voxel_z=0; voxel_z<dim.z; voxel_z++) {
            for (size_t voxel_y=0; voxel_y<dim.y; voxel_y++) {
                for (size_t voxel_x=0; voxel_x<dim.x; voxel_x++) {
                    tgt::svec3 pos = tgt::ivec3(voxel_x, voxel_y, voxel_z);
                    uint16_t v = segmentationRAM->voxel(pos);
                    if (v) {    // do not insert a background segment
                        auto it = segments.find(v);
                        if(it == segments.end()) {
                            std::vector<tgt::svec3> positions = std::vector<tgt::svec3>();
                            positions.push_back(pos);
                            segments.insert(std::make_pair(v, positions));
                        } else {
                            it->second.push_back(pos);
                        }
                    }
                }
            }
        }
        
        // segmentation can now be deleted
        delete segmentedImage;

        // compute centroids of each label
        for(auto it = segments.begin(); it != segments.end(); it++) {
            tgt::dvec3 avg = tgt::dvec3(0.0);
            std::vector<tgt::svec3> positions = it->second;
    
            for(size_t i = 0; i < positions.size(); i++) {
                avg += tgt::dvec3(positions.at(i));
            }
            avg = avg / tgt::dvec3(positions.size());
    
            // add offset to voxel position
            centroids.push_back(tgt::vec3(avg) + tgt::vec3(component.llf_)); 
        }

        // add centroids to local result
        localResult_.addSegment(centroids);

        // increment the counter in the cluster splitting processor for progress feedback
        boost::mutex::scoped_lock lock(globalValueMutex_);
        processor_->threadComponentCounter_++;
        lock.unlock();

        // only push data at the end (i.e., when interrupted while waiting for the empty queue)
        //pushLocalData();
    }
}

}   // namespace
