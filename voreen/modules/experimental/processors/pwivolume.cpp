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

#include "pwivolume.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"

using tgt::ivec3;
using tgt::svec3;
using tgt::vec3;

namespace voreen {

const std::string PWIVolume::loggerCat_("voreen.PWIVolume");

PWIVolume::PWIVolume()
    : CachingVolumeProcessor()
    , inport_(Port::INPORT, "volumes.in")
    , outport_(Port::OUTPORT, "volume.out")
    , plotOutport_(Port::OUTPORT, "plot.out")
{
    addPort(inport_);
    addPort(outport_);
    addPort(plotOutport_);
}

PWIVolume::~PWIVolume() {}

Processor* PWIVolume::create() const {
    return new PWIVolume();
}

void PWIVolume::process() {
    if(!inport_.getData())
        return;

    if(!outport_.hasData() || inport_.hasChanged())
        computePerfusion(inport_.getData());
}

void PWIVolume::computePerfusion(const VolumeList* volumeCollection) {
    if(volumeCollection->empty())
        return;

    const VolumeRAM* origin = volumeCollection->at(0)->getRepresentation<VolumeRAM>();
    if(!origin)
        return;

    VolumeRAM_Float* vol = new VolumeRAM_Float(origin->getDimensions());
    vol->clear();

    PlotData* newData = new PlotData(0,0);
    newData->reset(1, 1);
    newData->setColumnLabel(0, "Index");
    newData->setColumnLabel(1, "Size");

    for(size_t j = 1; j < volumeCollection->size(); j++) {
        float avgSum = 0.f;
        for(size_t i = 0; i < hmul(vol->getDimensions()); i++) {
            float diff = origin->getVoxelNormalized(i) - volumeCollection->at(j)->getRepresentation<VolumeRAM>()->getVoxelNormalized(i);
            vol->voxel(i) += diff;
            //avgSum += diff;
            avgSum += volumeCollection->at(j)->getRepresentation<VolumeRAM>()->getVoxelNormalized(i);
        }
        std::vector<PlotCellValue> pCellVector;
        pCellVector.push_back(PlotCellValue(static_cast<plot_t>(j)));
        pCellVector.push_back(PlotCellValue(avgSum / hmul(vol->getDimensions())));
        newData->insert(pCellVector);
    }

    Volume* result = new Volume(vol, volumeCollection->at(0)->getSpacing(), tgt::vec3(0.f));
    centerVolume(result);
    outport_.setData(result);
    plotOutport_.setData(newData);
}


} // namespace
