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

#include "rwmultilabelloopfinalizer.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumelist.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "voreen/core/voreenapplication.h"

#include "tgt/vector.h"

namespace voreen {

const std::string RWMultiLabelLoopFinalizer::loggerCat_("voreen.RandomWalker.RWMultiLabelLoopFinalizer");

RWMultiLabelLoopFinalizer::RWMultiLabelLoopFinalizer()
    : VolumeProcessor(),
      probabilitiesInport_(Port::INPORT, "probabilities.in"),
      segmentationInport_(Port::INPORT, "segmentation.in"),
      probabilitiesCollectionOutport_(Port::OUTPORT, "probabilitiesCollection.out"),
      segmentationCollectionOutport_(Port::OUTPORT, "segmentationCollection.out"),
      segmentationOutport_(Port::OUTPORT, "segmentation.out"),
      segmentationProbOutport_(Port::OUTPORT, "segmentationProb.out"),
      loopOutport_(Port::OUTPORT, "loop.out"),
      segmentationCollection_(0),
      probabilitiesCollection_(0),
      segmentationVolume_(0),
      segmentationProbVolume_(0)
{
    addPort(probabilitiesInport_);
    addPort(segmentationInport_);
    addPort(segmentationOutport_);
    addPort(segmentationProbOutport_);
    addPort(probabilitiesCollectionOutport_);
    addPort(segmentationCollectionOutport_);

    loopOutport_.setLoopPort(true);
    addPort(loopOutport_);
}

RWMultiLabelLoopFinalizer::~RWMultiLabelLoopFinalizer() {
}

Processor* RWMultiLabelLoopFinalizer::create() const {
    return new RWMultiLabelLoopFinalizer();
}

void RWMultiLabelLoopFinalizer::initialize() {
    VolumeProcessor::initialize();
}

void RWMultiLabelLoopFinalizer::deinitialize() {
    VolumeProcessor::deinitialize();
}

bool RWMultiLabelLoopFinalizer::isReady() const {
    bool ready = true;
    ready &= segmentationInport_.isReady();
    ready &= probabilitiesInport_.isReady();
    ready &= (segmentationOutport_.isReady() || segmentationProbOutport_.isReady() ||
              probabilitiesCollectionOutport_.isReady() || segmentationCollectionOutport_.isReady());
    return ready;
}

void RWMultiLabelLoopFinalizer::process() {

    if (loopOutport_.getLoopIteration() == 0) {
        segmentationOutport_.setData(0);
        segmentationProbOutport_.setData(0);
        delete segmentationVolume_;
        delete segmentationProbVolume_;
        segmentationVolume_ = 0;
        segmentationProbVolume_ = 0;

        segmentationCollectionOutport_.setData(0);
        probabilitiesCollectionOutport_.setData(0);

        if (segmentationCollection_) {
            // note: deleting a volume(handle) automatically removes from a collection it as contained by
            while (!segmentationCollection_->empty())
                delete segmentationCollection_->first();
            delete segmentationCollection_;
        }

        if (probabilitiesCollection_) {
            while (!probabilitiesCollection_->empty())
                delete probabilitiesCollection_->first();
            delete probabilitiesCollection_;
        }

        segmentationCollection_ = new VolumeList();
        probabilitiesCollection_ = new VolumeList();

        if (probabilitiesInport_.hasData()) {
            tgt::ivec3 dim = probabilitiesInport_.getData()->getDimensions();
            segmentationVolume_ = new VolumeRAM_UInt8(dim);
            segmentationVolume_->clear();
            segmentationProbVolume_ = new VolumeRAM_UInt16(dim);
            segmentationProbVolume_->clear();
        }

        runtime_ = clock();
    }
    tgtAssert(segmentationCollection_, "no segmentation volume collection");
    tgtAssert(probabilitiesCollection_, "no probabilities volume collection");

    int iteration = loopOutport_.getLoopIteration();

    const VolumeRAM* curSegVolume = 0;
    const VolumeRAM* curProbVolume = 0;
    if (segmentationInport_.hasData()) {
        curSegVolume = segmentationInport_.getData()->getRepresentation<VolumeRAM>();
    }
    if (probabilitiesInport_.hasData()) {
        curProbVolume = probabilitiesInport_.getData()->getRepresentation<VolumeRAM>();
    }

    if (curSegVolume) {
        Volume* clonedHandle = new Volume(curSegVolume->clone(), probabilitiesInport_.getData());
        segmentationCollection_->add(clonedHandle);
    }
    if (curProbVolume) {
        Volume* clonedHandle = new Volume(curProbVolume->clone(), probabilitiesInport_.getData());
        probabilitiesCollection_->add(clonedHandle);
    }

    // update multi-label segmentation result
    if (curProbVolume) {
        tgtAssert(segmentationVolume_, "no segmentation volume");
        tgtAssert(segmentationProbVolume_, "no prob volume");
        tgtAssert(curProbVolume->getDimensions() == segmentationProbVolume_->getDimensions(), "dimension mis-match");
        tgtAssert(curProbVolume->getDimensions() == segmentationVolume_->getDimensions(), "dimension mis-match");

        for (size_t voxelID = 0; voxelID<curProbVolume->getNumVoxels(); voxelID++) {
            float curProbValue = curProbVolume->getVoxelNormalized(voxelID);
            if (curProbValue > segmentationProbVolume_->getVoxelNormalized(voxelID)) {
                segmentationVolume_->voxel(voxelID) = iteration;
                segmentationProbVolume_->setVoxelNormalized(curProbValue, voxelID);
            }
        }
    }

    // last iteration: output data
    if (iteration == loopOutport_.getNumLoopIterations()-1) {
        segmentationOutport_.setData(new Volume(segmentationVolume_, probabilitiesInport_.getData()));
        segmentationProbOutport_.setData(new Volume(segmentationProbVolume_, probabilitiesInport_.getData()));
        segmentationCollectionOutport_.setData(segmentationCollection_, false);
        probabilitiesCollectionOutport_.setData(probabilitiesCollection_, false);
        LINFO("Multi-label runtime: " << (static_cast<float>(clock() - runtime_)/CLOCKS_PER_SEC) << " sec");
    }
}

}   // namespace
