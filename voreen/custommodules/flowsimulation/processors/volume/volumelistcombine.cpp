/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "volumelistcombine.h"

#include "voreen/core/ports/conditions/portconditionvolumelist.h"

namespace voreen {

VolumeListCombine::VolumeListCombine()
    : Processor()
    , inport_(Port::INPORT, "volumelist.input", "Volume Input", false)
    , outport_(Port::OUTPORT, "volume.output", "Volume Output", false)
    , combineMethod_("combineMethod", "Combine Method")
{
    addPort(inport_);
    inport_.addCondition(new PortConditionVolumeListEnsemble());
    addPort(outport_);

    addProperty(combineMethod_);
    combineMethod_.addOption("mean", "Mean", MEAN);
}

VolumeListCombine::~VolumeListCombine() {}

Processor* VolumeListCombine::create() const {
    return new VolumeListCombine();
}

void VolumeListCombine::process() {
    const VolumeList* data = inport_.getData();
    if(!data || data->empty()) {
        outport_.setData(nullptr);
    } else if(data->size() == 1) {
        outport_.setData(data->first(), false);
    } else if(inport_.hasChanged()) {
        // Take the first as a reference (we already know we got an ensemble).
        const VolumeBase* reference = data->first();

        Volume* combined = reference->clone();
        VolumeRAM* combinedRepresentation = combined->getWritableRepresentation<VolumeRAM>();
        tgt::svec3 dim = combined->getDimensions();

        CombineMethod combineMethod = combineMethod_.getValue();
        for(size_t i=0; i<data->size(); i++) {
            VolumeRAMRepresentationLock representation(data->at(i));
            for(size_t z=0; z < dim.z; z++) {
                for(size_t y=0; y < dim.y; y++) {
                    for(size_t x=0; x < dim.x; x++) {
                        for(size_t channel=0; channel<combined->getNumChannels(); channel++) {

                            float part = 0.0f;

                            switch(combineMethod) {
                            case MEAN:
                                part = representation->getVoxelNormalized(x, y, z, channel);
                                part /= data->size();
                                break;
                            default:
                                tgtAssert(false, "unhandled combine method");
                            }

                            float value = combinedRepresentation->getVoxelNormalized(x, y, z, channel);
                            value += part;
                            combinedRepresentation->setVoxelNormalized(value, x, y, z, channel);
                        }
                    }
                }
            }
        }

        outport_.setData(combined, true);
    }
}

}   // namespace
