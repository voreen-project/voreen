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

#include "volumeoctreelevelextractor.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/octree/volumeoctree.h"
#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"

namespace voreen {

const std::string VolumeOctreeLevelExtractor::loggerCat_("voreen.base.VolumeOctreeLevelExtractor");

VolumeOctreeLevelExtractor::VolumeOctreeLevelExtractor()
    : VolumeProcessor()
    , inport_(Port::INPORT, "volume.input", "Volume Input")
    , outport_(Port::OUTPORT, "volume.output", "Volume Output", false, Processor::VALID)
    , enableProcessing_("enabled", "Enable Processing", false)
    , mode_("mode", "Mode", Processor::VALID)
    //, channel_("channel", "Channel", 1, 1, 4)
    , maxEdgeLength_("edgelength", "Maximum Edge Length", Processor::VALID)
    , level_("level", "Octree Level", 0, 0, 7)
    , extractButton_("createvolume", "Create Volume")
    , automaticallyComputeOnChange_("automaticallycompute", "Automatically Update Output on Change", false, Processor::VALID)
    , buttonPressed_(false)
{
    inport_.onChange(MemberFunctionCallback<VolumeOctreeLevelExtractor>(this, &VolumeOctreeLevelExtractor::adjustPropertiesToInput));
    inport_.onChange(MemberFunctionCallback<VolumeOctreeLevelExtractor>(this, &VolumeOctreeLevelExtractor::adjustPropertyActivationState));
    enableProcessing_.onChange(MemberFunctionCallback<VolumeOctreeLevelExtractor>(this, &VolumeOctreeLevelExtractor::adjustPropertyActivationState));
    addPort(inport_);
    addPort(outport_);

    addProperty(enableProcessing_);
    enableProcessing_.onChange(MemberFunctionCallback<VolumeOctreeLevelExtractor>(this, &VolumeOctreeLevelExtractor::adjustPropertyActivationState));

    mode_.addOption("maxedge", "Select by Maximum Edge Length");
    mode_.addOption("octreelevel", "Select by Octree Level");
    mode_.onChange(MemberFunctionCallback<VolumeOctreeLevelExtractor>(this, &VolumeOctreeLevelExtractor::adjustPropertyActivationState));
    addProperty(mode_);

    //addProperty(channel_);

    maxEdgeLength_.addOption("64", "64 Voxels", 64);
    maxEdgeLength_.addOption("128", "128 Voxels", 128);
    maxEdgeLength_.addOption("256", "256 Voxels", 256);
    maxEdgeLength_.addOption("512", "512 Voxels", 512);
    maxEdgeLength_.onChange(MemberFunctionCallback<VolumeOctreeLevelExtractor>(this, &VolumeOctreeLevelExtractor::synchronizePropertiesOnChange));
    level_.onChange(MemberFunctionCallback<VolumeOctreeLevelExtractor>(this, &VolumeOctreeLevelExtractor::synchronizePropertiesOnChange));
    addProperty(maxEdgeLength_);
    addProperty(level_);

    addProperty(extractButton_);
    extractButton_.onClick(MemberFunctionCallback<VolumeOctreeLevelExtractor>(this, &VolumeOctreeLevelExtractor::buttonPressed));

    addProperty(automaticallyComputeOnChange_);
}

Processor* VolumeOctreeLevelExtractor::create() const {
    return new VolumeOctreeLevelExtractor();
}

void VolumeOctreeLevelExtractor::initialize() {
    VolumeProcessor::initialize();
    adjustPropertyActivationState();
}

void VolumeOctreeLevelExtractor::process() {
    tgtAssert(inport_.getData(), "No input volume");

    // processor has been invalidated -> delete output volume
    outport_.setData(0);

    // if the processor is not enabled: do nothing
    if (!enableProcessing_.get())
        return;

    // check if either the button has been pressed or the automatic computation is enabled
    if (!buttonPressed_ && !automaticallyComputeOnChange_.get())
        return;

    // set button press to false to prevent multiple computations
    buttonPressed_ = false;

    // get input volume and check if it has an octree representation
    const VolumeBase* inputVolume = inport_.getData();
    if (!inputVolume->hasRepresentation<VolumeOctree>()) {
        LERROR("Input volume does not have an octree representation");
        return;
    }

    // perform computation and set output volume to outport (outport should take ownership of the volume)
    const VolumeOctree* octree = inputVolume->getRepresentation<VolumeOctree>();
    VolumeRAM* levelVolume = octree->createVolume(octree->getNumLevels() - 1 - level_.get());
    tgt::vec3 ratio = tgt::vec3(inputVolume->getDimensions()) / tgt::vec3(levelVolume->getDimensions());

    Volume* result = new Volume(levelVolume, inputVolume);
    result->setSpacing(inputVolume->getSpacing() * ratio);
    outport_.setData(result);
}

void VolumeOctreeLevelExtractor::adjustPropertiesToInput() {
    if (!inport_.getData() || !inport_.getData()->hasRepresentation<VolumeOctree>())
        return;

    const VolumeOctree* octree = inport_.getData()->getRepresentation<VolumeOctree>();

    //channel_.setMaxValue(octree->getNumChannels());

    // get the absolute maximum level
    size_t maxLevel = octree->getActualTreeDepth() - 1;

    // get the brick size (power-of-two) and compute the maximum level for a volume <= 512^3
    size_t brickEdge = octree->getBrickDim().x;
    if (brickEdge >= 512)
        maxLevel = 0;
    else {
        //size_t maxVolumeEdge = tgt::max(octree->getVolumeDim());
        size_t s = 512;
        size_t maxLevel = std::min(maxLevel, static_cast<size_t>(tgt::ilog2(static_cast<int>(s) / static_cast<int>(brickEdge))));
    }

    // set the level property max value
    level_.setMaxValue(maxLevel);

    synchronizePropertiesOnChange();
}

void VolumeOctreeLevelExtractor::synchronizePropertiesOnChange() {
    if (!inport_.getData() || !inport_.getData()->hasRepresentation<VolumeOctree>())
        return;

    // get the brick size (power-of-two)
    const VolumeOctree* octree = inport_.getData()->getRepresentation<VolumeOctree>();
    size_t brickEdge = octree->getBrickDim().x;

    // depending on the mode, compute level_ from maxEdgeLength_ or vice versa (using the brick size and volume dimensions)
    if (mode_.getKey() == "octreelevel") {
        // get the current octree level and adjust the edge length property according to it
        size_t level = level_.get();
        size_t dim = static_cast<size_t>(std::pow(2.0, static_cast<double>(level)) * static_cast<double>(brickEdge));

        if (dim >= 512)
            maxEdgeLength_.selectByKey("512");
        else if (dim >= 256)
            maxEdgeLength_.selectByKey("256");
        else if (dim >= 128)
            maxEdgeLength_.selectByKey("128");
        else
            maxEdgeLength_.selectByKey("64");
    }
    else {
        int maxEdgeLength = maxEdgeLength_.getValue();
        int level = std::min(level_.getMaxValue(), tgt::ilog2(maxEdgeLength / static_cast<int>(brickEdge)));
        level_.set(std::min(level_.getMaxValue(), level));

        // check if this setting is even possible and if not switch
        size_t dim = static_cast<size_t>(std::pow(2.0, static_cast<double>(level_.get())) * static_cast<double>(brickEdge));
        if (dim < static_cast<size_t>(maxEdgeLength_.getValue())) {
            if (dim <= 64)
                maxEdgeLength_.selectByKey("64");
            else if (dim <= 128)
                maxEdgeLength_.selectByKey("128");
            else if (dim <= 256)
                maxEdgeLength_.selectByKey("256");
        }
        invalidate();
    }
}

void VolumeOctreeLevelExtractor::adjustPropertyActivationState() {
    // depending on the mode, setReadOnlyFlag at level_ or maxEdgeLength_
    maxEdgeLength_.setReadOnlyFlag(mode_.getKey() == "octreelevel");
    level_.setReadOnlyFlag(mode_.getKey() == "maxedge");

    // if processor is mot enabled or no input volume / input volume has no octree representation: disable the compute button
    if (!enableProcessing_.get() || !inport_.getData() || !inport_.getData()->hasRepresentation<VolumeOctree>()) {
        extractButton_.setReadOnlyFlag(true);
        buttonPressed_ = false;
    }
    else
        extractButton_.setReadOnlyFlag(false);

}

void VolumeOctreeLevelExtractor::buttonPressed() {
    buttonPressed_ = true;
    invalidate();
}

}   // namespace
