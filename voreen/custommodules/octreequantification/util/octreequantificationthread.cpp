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

#include "octreequantificationthread.h"

#include "voreen/core/datastructures/octree/octreeutils.h"

namespace voreen {

OctreeQuantificationThread::OctreeQuantificationThread(OctreeSegmentationQuantification* quantificationProcessor, OctreeQuantificationNodeQueue& queue, 
            const VolumeOctree* octree, std::vector<const VolumeRAM*> rois, std::vector<tgt::vec3> coordinateRatios, 
            OctreeSegmentationQuantification::SegmentationInterpolationMode simMode, tgt::svec3 volumeURB, size_t numChannels, boost::mutex& globalValueMutex)
    : BackgroundThread()
    , quantificationProcessor_(quantificationProcessor)
    , queue_(queue)
    , volumeOctree_(octree)
    , rois_(rois)
    , coordinateRatios_(coordinateRatios)
    , simMode_(simMode)
    , volumeURB_(volumeURB)
    , numChannels_(numChannels)
    , globalValueMutex_(globalValueMutex)
    , localResults_(numChannels, rois.size())
{ }  

void OctreeQuantificationThread::clearLocalData() {
    localResults_.clear();
}

void OctreeQuantificationThread::pushLocalData() {
    boost::mutex::scoped_lock lock(globalValueMutex_);

    quantificationProcessor_->currentResults_.addResults(localResults_);
}

void OctreeQuantificationThread::handleInterruption() {
    pushLocalData();
}


void OctreeQuantificationThread::threadMain() {

    // variables needed during iteration over bricks and voxels
    tgt::svec3 brickDim = volumeOctree_->getBrickDim();
    size_t voxelsPerBrick = volumeOctree_->getNumVoxelsPerBrick();
    tgt::svec3 brickLLF, brickURB, pos, brickBufferPos;
    size_t brickBufferIndex; 

    std::vector<uint16_t> voxelValues;
    std::vector<float> roiProbabilities;
    std::vector<bool> inROIs;
    //size_t voxelAScatterValue, voxelBScatterValue, voxelCScatterValue;
    //float outerSegmentationProbability, innerStructureProbability;
    //bool isInOuterBoundary, isInInnerStructure;
    const uint16_t* nodeBrick = 0;
    size_t blockVoxelCounter = 0;

    clearLocalData();

    // we use this while-loop because at some point we will wait for the empty queue and get the interruptAndJoin signal, which will cause this thread to terminate as it is the only interruption point
    while(true) {

        // if the queue is empty the thread will be paused, which also provides our interruption point for termination
        std::pair<tgt::IntBounds, const VolumeOctreeNode*> quantificationPair = queue_.pop();
        
        tgt::IntBounds currentBounds = quantificationPair.first;
        const VolumeOctreeNode* currentNode = quantificationPair.second;

        if (currentNode->hasBrick()) {
            nodeBrick = volumeOctree_->getNodeBrick(currentNode);
        }
        else {
            // node is homogeneous -> take average voxel value for each channel
            voxelValues.clear();
            nodeBrick = 0;
            for (size_t c = 0; c < numChannels_; ++c) 
                voxelValues.push_back(currentNode->getAvgValue(c));
        }

        //TODO: could apply axis-aligned clipping here?
        brickLLF = tgt::svec3(currentBounds.getLLF());
        brickURB = tgt::min(tgt::svec3(currentBounds.getURB()), volumeURB_); // <- this makes sure that we do not count voxels which are added to fill up the block outside of the volume 
        
        for (pos = brickLLF; pos.z <= brickURB.z; ++pos.z) {
            for (pos.y = brickLLF.y; pos.y <= brickURB.y; ++pos.y) {
                for (pos.x = brickLLF.x; pos.x <= brickURB.x; ++pos.x) {

                    // clear temporary values
                    roiProbabilities.clear();
                    inROIs.clear();

                    if (nodeBrick) {
                        voxelValues.clear();
                        // read voxel values from brick
                        brickBufferPos = pos - brickLLF;
                        //brickBufferIndex = brickBufferPos.z * brickDim.y * brickDim.x + brickBufferPos.y * brickDim.x + brickBufferPos.x;
                        brickBufferIndex = cubicCoordToLinear(brickBufferPos, brickDim) * numChannels_;
                        for (size_t c = 0; c < numChannels_; ++c) 
                            voxelValues.push_back(nodeBrick[brickBufferIndex + c]);
                    }

                    // from now on we do the same, since we already have the voxel values

                    // get probabilities for each of the ROIS
                    if (simMode_ == OctreeSegmentationQuantification::SIM_NEAREST) {
                        // iterate over rois
                        for (size_t r = 0; r < rois_.size(); ++r) {
                            if (rois_.at(r))
                                roiProbabilities.push_back(rois_.at(r)->getVoxelNormalized(tgt::svec3(tgt::iround(tgt::vec3(pos) * coordinateRatios_.at(r)))));
                            else
                                roiProbabilities.push_back(0.f);
                        }
                    }
                    else if (simMode_ == OctreeSegmentationQuantification::SIM_LINEAR) {
                        // iterate over rois
                        for (size_t r = 0; r < rois_.size(); ++r) {
                            if (rois_.at(r))
                                roiProbabilities.push_back(rois_.at(r)->getVoxelNormalizedLinear(tgt::vec3(pos) * coordinateRatios_.at(r)));
                            else
                                roiProbabilities.push_back(0.f);
                        }
                    }
                    else {
                        // iterate over rois
                        for (size_t r = 0; r < rois_.size(); ++r) {
                            if (rois_.at(r))
                                roiProbabilities.push_back(rois_.at(r)->getVoxelNormalizedCubic(tgt::vec3(pos) * coordinateRatios_.at(r)));
                            else
                                roiProbabilities.push_back(0.f);
                        }
                    }

                    // determine boolean values
                    for (size_t i = 0; i < rois_.size(); ++i) {
                        if (roiProbabilities.at(i) >= 0.5f)
                            inROIs.push_back(true);
                        else
                            inROIs.push_back(false);
                    }

                    // add voxel to local results
                    localResults_.addVoxel(voxelValues, inROIs);
                    blockVoxelCounter++;
                }
            }
        }

        // be sure to release the brick after use
        if (nodeBrick) 
            volumeOctree_->releaseNodeBrick(currentNode);

        // give feedback about quantification status

        // directly write back our local data for progress feedback as well as prevention of data loss from interruption at the end
        boost::mutex::scoped_lock lock(globalValueMutex_);
        quantificationProcessor_->threadVoxelCounter_ += blockVoxelCounter;
        //lock.unlock();
        blockVoxelCounter = 0;
    }
}

}   // namespace
