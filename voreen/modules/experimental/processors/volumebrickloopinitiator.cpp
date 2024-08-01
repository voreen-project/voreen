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

#include "volumebrickloopinitiator.h"
#include "tgt/glmath.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/volumegl.h"

namespace voreen {

const std::string VolumeBrickLoopInitiator::loggerCat_("voreen.experimental.VolumeBrickLoopInitiator");

VolumeBrickLoopInitiator::VolumeBrickLoopInitiator()
    : Processor(),
      volOutput_("volOutput", "Output Usage"),
      allocDirective_("allocDirective", "Bricking Directive"),
      //maxMemorySize_("maxMemorySize", "Max Output Size (MB)", 0, 0, 10000),
      divideBy_("divideBy", "Divide volume by certain number", 0, 0, 32),
      brickSize_("brickSize_", "Brick Size (Each Dim)"),
      dimension_("dimension", "Dimensions", tgt::ivec3(256), tgt::ivec3(1), tgt::ivec3(4096)),
      inport_(Port::INPORT, "volume.in"),
      outport_(Port::OUTPORT, "volumebrick.out"),
      loopInport_(Port::INPORT, "loop.in")
{
    volOutput_.addOption("CPU", "CPU (Main Memory)", CPU);
    volOutput_.addOption("GPU", "GPU (TextureGL)", GPU);
    volOutput_.select("CPU");
    addProperty(volOutput_);

    allocDirective_.addOption("MIN", "Minimum Number of Bricks", MIN_NUM_BRICK);
    //allocDirective_.addOption("DIV", "Divide by Number", DIV_BY_NUM);
    allocDirective_.addOption("SIZE", "Size of Bricks", SIZE_BRICK);
    allocDirective_.select("MIN");
    addProperty(allocDirective_);

    //addProperty(maxMemorySize_);
    addProperty(divideBy_);

    //Add values which have a power of two
    brickSize_.addOption("1", "1", 1);
    brickSize_.addOption("2", "2", 2);
    brickSize_.addOption("4", "4", 4);
    brickSize_.addOption("8", "8", 8);
    brickSize_.addOption("16", "16", 16);
    brickSize_.addOption("32", "32", 32);
    brickSize_.addOption("64", "64", 64);
    brickSize_.addOption("128", "128", 128);
    brickSize_.addOption("256", "256", 256);
    brickSize_.addOption("512", "512", 512);
    brickSize_.addOption("1024", "1024", 1024);
    brickSize_.select("32");
    addProperty(brickSize_);

    addProperty(dimension_);
    dimension_.setWidgetsEnabled(false);

    addPort(inport_);
    addPort(outport_);
    addPort(loopInport_);

    loopInport_.setLoopPort(true);
}

Processor* VolumeBrickLoopInitiator::create() const {
    return new VolumeBrickLoopInitiator();
}

bool VolumeBrickLoopInitiator::isReady() const {
    return (inport_.isReady() && outport_.isReady());
}

void VolumeBrickLoopInitiator::defineMaxAllocationSize(){
    if (inport_.hasData()){
        const VolumeRAM* vol = 0;
        currentBrickSize_ = inport_.getData()->getDimensions();
        try{
            vol = inport_.getData()->getRepresentation<VolumeRAM>();
            allFitsOnMainMem_ = true;
        }
        catch (std::bad_alloc&){
            vol = 0;
            allFitsOnMainMem_ = false;
            const VolumeDisk* diskRep = inport_.getData()->getRepresentation<VolumeDisk>();
            Volume* vh = new Volume(0, tgt::vec3(0.f), tgt::vec3(0.f));
            tgt::svec3 newDim = diskRep->getDimensions();
            size_t divide = 2;
            while(!vol){
                try{
                    newDim.z /= divide;
                    vh->addRepresentation(diskRep->getSubVolume(newDim));
                    vol = vh->getRepresentation<VolumeRAM>();
                }
                catch (std::bad_alloc&){
                    vol = 0;
                    divide++;
                }
            }
            currentBrickSize_ = newDim;
            delete vh;
        }
        if(volOutput_.getValue() == GPU){
            const VolumeGL* volGL = 0;
            try{
                volGL = inport_.getData()->getRepresentation<VolumeGL>();
                //volGL->generateTexture(vol);
                allFitsOnGpuMem_ = true;
            }
            catch (std::bad_alloc&){
                allFitsOnGpuMem_ = false;
                volGL = 0;
                tgt::svec3 newDim = vol->getDimensions();
                Volume* vh = new Volume(0, tgt::vec3(0.f), tgt::vec3(0.f));
                size_t divide = 2;
                while(!volGL){
                    try{
                        newDim.z /= divide;
                        vh->addRepresentation(vol->getSubVolume(newDim));
                        volGL = vh->getRepresentation<VolumeGL>();
                        ///TODO: VolumeGL doesn't upload texture properly as it should cast bad_alloc for very large volume.
                    }
                    catch (std::bad_alloc&){
                        volGL = 0;
                        divide++;
                    }
                }
                currentBrickSize_ = newDim;
                delete vh;
            }
        }
    }
}

tgt::svec3 VolumeBrickLoopInitiator::calculateNewOffset(int loopIteration, tgt::svec3 brickSize, tgt::svec3 dims, tgt::svec3 offset){
    if (loopIteration>0){
        //move offset ahead
        offset.x += brickSize.x;

        if(offset.x >= dims.x){
            offset.x = 0;
            offset.y += brickSize.y;
        }
        if(offset.y >= dims.y){
            offset.x = 0;
            offset.y = 0;
            offset.z += brickSize.z;
        }

        return offset;
    }
    else
        return tgt::svec3(0,0,0);
}

tgt::svec3 VolumeBrickLoopInitiator::calculateCorrectBrickSize(tgt::svec3 brickSize, tgt::svec3 dims, tgt::svec3 offset){
    if(offset.x + brickSize.x > dims.x)
        brickSize.x = dims.x - offset.x;
    if(offset.y + brickSize.y > dims.y)
        brickSize.y = dims.y - offset.y;
    if(offset.z + brickSize.z > dims.z)
        brickSize.z = dims.z - offset.z;

    return brickSize;
}

int VolumeBrickLoopInitiator::calculateLoopIterations(tgt::svec3 brickSize, tgt::svec3 dims, bool hasData){
    if (!hasData)
        return 0;
    else{
        //calculate number of bricks
        tgt::ivec3 bricksPerDim;
        bricksPerDim.x = tgt::iceil(static_cast<double>(dims.x)/static_cast<double>(brickSize.x));
        bricksPerDim.y = tgt::iceil(static_cast<double>(dims.y)/static_cast<double>(brickSize.y));
        bricksPerDim.z = tgt::iceil(static_cast<double>(dims.z)/static_cast<double>(brickSize.z));
        return bricksPerDim.x*bricksPerDim.y*bricksPerDim.z;
    }
}

void VolumeBrickLoopInitiator::createNewVolume(Volume* volumeHandle, tgt::svec3 size, tgt::svec3 offset){
    if(volumeHandle){
        delete volumeHandle;
        volumeHandle = 0;
    }

    if(allocDirective_.getValue() == MIN_NUM_BRICK){
        if(volOutput_.getValue() == GPU){
            if(allFitsOnMainMem_){
                if(!allFitsOnGpuMem_){
                    const VolumeRAM* vol = inport_.getData()->getRepresentation<VolumeRAM>();
                    volumeHandle = new Volume(vol->getSubVolume(size, offset), inport_.getData());
                }
            }
            else{
                const VolumeDisk* diskRep = inport_.getData()->getRepresentation<VolumeDisk>();
                volumeHandle = new Volume(diskRep->getSubVolume(size, offset), inport_.getData());
                volumeHandle->addRepresentation(volumeHandle->getWritableRepresentation<VolumeRAM>());
            }
        }
        else if(!allFitsOnMainMem_){
            const VolumeDisk* diskRep = inport_.getData()->getRepresentation<VolumeDisk>();
            volumeHandle = new Volume(diskRep->getSubVolume(size, offset), inport_.getData());
        }
    }
    else if(allocDirective_.getValue() == SIZE_BRICK) {
        const VolumeRAM* vol = inport_.getData()->getRepresentation<VolumeRAM>();
        volumeHandle = new Volume(vol->getSubVolume(size, offset), inport_.getData());
    }
}

void VolumeBrickLoopInitiator::process() {
    LGL_ERROR;

    if (!inport_.isReady()) {
        delete volHandle_;
        volHandle_ = 0;
        currentOffset_ = tgt::svec3(0,0,0);
        outport_.setData(0);
        return;
    }

    if(inport_.hasData()){
        if(inport_.hasChanged()){
            dimension_.set(inport_.getData()->getDimensions());

            if(allocDirective_.getValue() == MIN_NUM_BRICK)
                defineMaxAllocationSize();
            else if(allocDirective_.getValue() == SIZE_BRICK)
                currentBrickSize_ = tgt::svec3(brickSize_.getValue());

            loopInport_.setNumLoopIterations(calculateLoopIterations(currentBrickSize_, inport_.getData()->getDimensions(), inport_.hasData()));

            inport_.setValid();
        }

        currentOffset_ = calculateNewOffset(loopInport_.getLoopIteration(), currentBrickSize_, inport_.getData()->getDimensions(), currentOffset_);

        currentBrickSize_ = calculateCorrectBrickSize(currentBrickSize_, inport_.getData()->getDimensions(), currentOffset_);

        createNewVolume(volHandle_, currentBrickSize_, currentOffset_);
    }

    // put out result volume
    if(volHandle_)
        outport_.setData(volHandle_);
    else
        outport_.setData(const_cast<VolumeBase*>(inport_.getData()));
}

void VolumeBrickLoopInitiator::initialize() {
    Processor::initialize();

    allocDirective_.onChange(CallMemberAction<VolumeBrickLoopInitiator>(this, &VolumeBrickLoopInitiator::showAndHideProperties));

    showAndHideProperties();

    currentBrickSize_ = tgt::svec3(0,0,0);
    currentOffset_ = tgt::svec3(0,0,0);

    volHandle_ = 0;
}

void VolumeBrickLoopInitiator::deinitialize() {
    Processor::deinitialize();
}

void VolumeBrickLoopInitiator::showAndHideProperties(){
    divideBy_.setVisible(allocDirective_.getValue() == DIV_BY_NUM);
    brickSize_.setVisible(allocDirective_.getValue() == SIZE_BRICK);
}

} // voreen namespace
