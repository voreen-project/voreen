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

#include "volumelistloopfinalizer.h"

#include "voreen/core/datastructures/imagesequence.h"
#include "tgt/gpucapabilities.h"

namespace voreen {

const std::string VolumeListLoopFinalizer::loggerCat_("voreen.volumelistloopfinalizer");

VolumeListLoopFinalizer::VolumeListLoopFinalizer()
    : Processor(),
      inport_(Port::INPORT, "rendering.in", "Image Input"),
      outport_(Port::OUTPORT, "imagesequence.out", "ImageSequence Output"),
      loopOutport_(Port::OUTPORT, "loop.out", "Loop Outport"),
      progress_("progress", "Progress"),
      volumeList_(nullptr)
{
    loopOutport_.setLoopPort(true);


    addPort(inport_);
    addPort(outport_);
    addPort(loopOutport_);

    addProgressBar(&progress_);
    addProperty(progress_);
}

Processor* VolumeListLoopFinalizer::create() const {
    return new VolumeListLoopFinalizer();
}

bool VolumeListLoopFinalizer::isReady() const {
    return (outport_.isReady());
}

void VolumeListLoopFinalizer::process() {
    // clear current list, if no valid data at inport
    if (!inport_.isReady()) {
        outport_.setData(nullptr);
        return;
    }

    int iteration = loopOutport_.getLoopIteration();
    // first iteration: clear previous volume list
    if (iteration == 0) {
        outport_.setData(nullptr);
        volumeList_.reset(new VolumeList());
    }

    setProgress(static_cast<float>(iteration)/std::max(1, loopOutport_.getNumLoopIterations() - 1));

    // add current volume to volumelist
    volumeList_->add(getVolumeFromInport());

    // last iteration: volumelist is complete => write it to outport
    if (iteration == loopOutport_.getNumLoopIterations()-1) {
        outport_.setData(volumeList_.get(), false);
    } else {
        loopOutport_.invalidatePort();
    }
}

VolumeBase* VolumeListLoopFinalizer::getVolumeFromInport() {
    const VolumeBase* inputVol = inport_.getData();
    tgtAssert(inputVol, "No input data");

    //TODO: maybe do something more intelligent than copying? what about large disk volumes?
    return new Volume(inputVol->getRepresentation<VolumeRAM>()->clone(), inputVol);
}

} // voreen namespace
