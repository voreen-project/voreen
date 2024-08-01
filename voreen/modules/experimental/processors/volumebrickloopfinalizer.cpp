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

#include "volumebrickloopfinalizer.h"

#include "tgt/logmanager.h"

namespace voreen {

using tgt::Texture;
using tgt::ivec3;

const std::string VolumeBrickLoopFinalizer::loggerCat_("voreen.experimental.VolumeBrickLoopFinalizer");

VolumeBrickLoopFinalizer::VolumeBrickLoopFinalizer()
    : Processor(),
      dimension_("dimension", "Dimensions", tgt::ivec3(256), tgt::ivec3(1), tgt::ivec3(4096)),
      inport_(Port::INPORT, "volumebrick.in"),
      outport_(Port::OUTPORT, "volume.out"),
      loopOutport_(Port::OUTPORT, "loop.out"),
      volume_(0)
{
    loopOutport_.setLoopPort(true);

    addProperty(dimension_);
    dimension_.setWidgetsEnabled(false);

    addPort(inport_);
    addPort(outport_);
    addPort(loopOutport_);
}

VolumeBrickLoopFinalizer::~VolumeBrickLoopFinalizer() {
    delete volume_;
}

Processor* VolumeBrickLoopFinalizer::create() const {
    return new VolumeBrickLoopFinalizer();
}

bool VolumeBrickLoopFinalizer::isReady() const {
    return (outport_.isReady());
}

void VolumeBrickLoopFinalizer::defineNewOffset(){
    int iter = loopOutport_.getLoopIteration();
    if (iter>0){
        //move offset ahead
        tgt::svec3 volOutSize = tgt::svec3(dimension_.get());
        currentOffset_.x += lastVolumeSize_.x;

        if(currentOffset_.x >= volOutSize.x){
            currentOffset_.x = 0;
            currentOffset_.y += lastVolumeSize_.y;
        }
        if(currentOffset_.y >= volOutSize.y){
            currentOffset_.x = 0;
            currentOffset_.y = 0;
            currentOffset_.z += lastVolumeSize_.z;
        }
    }
    else
        currentOffset_ = tgt::svec3(0,0,0);

    lastVolumeSize_ = inport_.getData()->getDimensions();
}

void VolumeBrickLoopFinalizer::process() {

    //if no valid data, clear
    if (!inport_.isReady()) {
        delete volume_;
        volume_ = 0;
        currentOffset_ = tgt::svec3(0,0,0);
        lastVolumeSize_ = tgt::svec3(0,0,0);
        outport_.setData(0);
        return;
    }

    // first iteration: clear previous volume and create new one
    if (loopOutport_.getLoopIteration() == 0) {
        VolumeBase* volumeHandle = const_cast<VolumeBase*>(outport_.getData());
        if(volumeHandle) {
            Volume* temp = dynamic_cast<Volume*>(volumeHandle);
            temp->deleteAllRepresentations();
            delete temp;
            outport_.setData(0);
        }

        if (inport_.hasData())
            volume_ = inport_.getData()->getRepresentation<VolumeRAM>()->createNew(tgt::svec3(dimension_.get()), true);
    }

    // add volume as part of our volume
    if (inport_.hasData()){
        defineNewOffset();
        volume_->setSubVolume(inport_.getData()->getRepresentation<VolumeRAM>(), currentOffset_);
    }

    // last iteration: volume operations are complete => write it to outport
    if (loopOutport_.getLoopIteration() == loopOutport_.getNumLoopIterations()-1) {
        Volume* vh = new Volume(volume_, tgt::vec3(1.0f), tgt::vec3(0.0f));//FIXME: metadata
        oldVolumePosition(vh); //FIXME
        outport_.setData(vh, false);
    } else {
        loopOutport_.invalidatePort();
    }

}

void VolumeBrickLoopFinalizer::initialize() {
    Processor::initialize();
}

void VolumeBrickLoopFinalizer::deinitialize() {
    Processor::deinitialize();
}

} // voreen namespace
