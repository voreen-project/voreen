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

#include "physicalclippinglinker.h"

namespace voreen {

const std::string PhysicalClippingLinker::loggerCat_("voreen.ensembleanalysis.PhysicalClippingLinker");

PhysicalClippingLinker::PhysicalClippingLinker()
    : Processor()
    , inport_(Port::INPORT, "volumeinport", "Volume Inport")
    , clipRegionVoxel_("clipRegionVoxel", "Clip Region Voxel", tgt::IntBounds(tgt::ivec3(0), tgt::ivec3(10000)), tgt::ivec3(0), tgt::ivec3(10000))
    , clipRegionPhysical_("clipRegionPhysical", "Clip Region Physical", tgt::Bounds(tgt::vec3(0), tgt::vec3(10000)), tgt::vec3(0), tgt::vec3(10000))
{
    addPort(inport_);
    ON_CHANGE(inport_, PhysicalClippingLinker, inputVolumeChanged);

    ON_CHANGE(clipRegionVoxel_, PhysicalClippingLinker, clipRegionVoxelChanged);
    ON_CHANGE(clipRegionPhysical_, PhysicalClippingLinker, clipRegionPhysicalChanged);
    addProperty(clipRegionVoxel_);
    addProperty(clipRegionPhysical_);
}

Processor* PhysicalClippingLinker::create() const {
    return new PhysicalClippingLinker();
}

void PhysicalClippingLinker::process() {
    // do nothing
}

void PhysicalClippingLinker::inputVolumeChanged() {
    if (!inport_.getData())
        return;

    // set max values
    clipRegionVoxel_.setMaxValue(tgt::ivec3(inport_.getData()->getDimensions()) - tgt::ivec3::one);
    clipRegionPhysical_.setMinValue(inport_.getData()->getLLF());
    clipRegionPhysical_.setMaxValue(inport_.getData()->getURB());
}

void PhysicalClippingLinker::clipRegionVoxelChanged() {
    if (!inport_.getData())
        return;

    tgt::vec3 llf = inport_.getData()->getVoxelToPhysicalMatrix() * tgt::vec3(clipRegionVoxel_.get().getLLF());
    tgt::vec3 urb = inport_.getData()->getVoxelToPhysicalMatrix() * tgt::vec3(clipRegionVoxel_.get().getURB());
    clipRegionPhysical_.set(tgt::Bounds(llf, urb));
}

void PhysicalClippingLinker::clipRegionPhysicalChanged() {
    if (!inport_.getData())
        return;

    tgt::ivec3 llf = inport_.getData()->getPhysicalToVoxelMatrix() * clipRegionPhysical_.get().getLLF();
    tgt::ivec3 urb = inport_.getData()->getPhysicalToVoxelMatrix() * clipRegionPhysical_.get().getURB();
    clipRegionVoxel_.set(tgt::IntBounds(llf, urb));
}

}   // namespace
