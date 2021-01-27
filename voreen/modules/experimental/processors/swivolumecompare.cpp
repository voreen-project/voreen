/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "swivolumecompare.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"
using tgt::ivec3;
using tgt::svec3;
using tgt::vec3;

namespace voreen {

const std::string SWIVolumeCompare::loggerCat_("voreen.SWIVolumeCompare");

SWIVolumeCompare::SWIVolumeCompare()
    : CachingVolumeProcessor()
    , roi_(Port::INPORT, "volume.mask")
    , swi_(Port::INPORT, "volume.swi")
    , outport_(Port::OUTPORT, "volume.out")
    //, plotOutport_(Port::OUTPORT, "plot.out")
    , thresh_("threshold", "SWI threshold", 0.5f, 0.f, 1.f)
    , ratio_("ratio", "Resulting Ratio", 0.f, 0.f, 1.f)
    , recalc_(true)
{
    addProperty(thresh_);
    addProperty(ratio_);

    addPort(roi_);
    addPort(swi_);
    addPort(outport_);
    //addPort(plotOutport_);

    thresh_.onChange(MemberFunctionCallback<SWIVolumeCompare>(this, &SWIVolumeCompare::enableRecalc));
    ratio_.setReadOnlyFlag(true);
}

SWIVolumeCompare::~SWIVolumeCompare() {}

Processor* SWIVolumeCompare::create() const {
    return new SWIVolumeCompare();
}

void SWIVolumeCompare::enableRecalc() {
    recalc_ = true;
}

void SWIVolumeCompare::process() {
    if(!roi_.getData() || !swi_.getData())
        return;

    if(!outport_.hasData() || roi_.hasChanged() || swi_.hasChanged())
        enableRecalc();

    if(recalc_)
        computeSWIDifference();
}

void SWIVolumeCompare::computeSWIDifference() {

    const VolumeRAM* swi = swi_.getData()->getRepresentation<VolumeRAM>();
    const VolumeRAM_UInt8* mask = roi_.getData()->getRepresentation<VolumeRAM_UInt8>();

    if(!swi || !mask)
        return;

    VolumeRAM_UInt8* vol = new VolumeRAM_UInt8(swi->getDimensions());
    vol->clear();

    int pos = 0;
    int roi = 0;

    float thresh = thresh_.get();

    for(size_t i = 0; i < hmul(roi_.getData()->getDimensions()); i++) {
        if(!mask->voxel(i))
            continue;

        std::cout << swi->getVoxelNormalized(i) << std::endl;
        if(swi->getVoxelNormalized(i) > thresh) {
            vol->voxel(i) = 255;
            pos++;
        }
        roi++;
    }

    ratio_.set((float)(pos) / roi);

    Volume* result = new Volume(vol, swi_.getData());
    centerVolume(result);
    outport_.setData(result);
    //plotOutport_.setData(newData);
}


} // namespace
