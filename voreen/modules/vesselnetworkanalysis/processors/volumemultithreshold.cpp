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

#include "volumemultithreshold.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"

#include "voreen/core/datastructures/callback/lambdacallback.h"

namespace voreen {

const std::string VolumeMultiThreshold::loggerCat_("voreen.vesselnetworkanalysis.volumemultithreshold");

VolumeMultiThreshold::VolumeMultiThreshold()
    : VolumeProcessor()
    , inport_(Port::INPORT, "volumethinning.inport", "Volume Input")
    , outport_(Port::OUTPORT, "volumethinning.outport", "Volume Output", false, Processor::VALID /* do not update if something connects to this port*/)
    , binarizationThresholdRange_("binarizationThresholdRange","Binarization Threshold Range", 0, 0, std::numeric_limits<float>::max(), 0, INT_MAX, Processor::VALID)
    , steps_("steps","Number of Output Volumes", 20, 0, 100, Processor::VALID)
    , start_("start","Compute",true)
    , startPressed_(false)
{
    addPort(inport_);
    addPort(outport_);
    addProperty(binarizationThresholdRange_);
    addProperty(steps_);
    addProperty(start_);

    ON_CHANGE_LAMBDA(start_, [this] () {
                startPressed_ = true;
            });
}

VolumeMultiThreshold::~VolumeMultiThreshold() {
}

static VolumeRAM* applyThreshold(const VolumeRAM& vol, float normalizedThreshold) {
    VolumeRAM* output = VolumeGeneratorUInt8().create(vol.getDimensions());

    for(size_t z = 0; z < vol.getDimensions().z; ++z) {
        for(size_t y = 0; y < vol.getDimensions().y; ++y) {
            for(size_t x = 0; x < vol.getDimensions().x; ++x) {
                tgt::svec3 pos(x,y,z);
                float value = vol.getVoxelNormalized(pos);
                output->setVoxelNormalized(value > normalizedThreshold ? 1.0f : 0.0f, pos);
            }
        }
    }
    return output;
}

void VolumeMultiThreshold::process() {
    const VolumeBase* invol = inport_.getData();
    if(!invol) {
        return;
    }
    if(inport_.hasChanged()) {
        const VolumeMinMax* vmm = invol->getDerivedData<VolumeMinMax>();
        binarizationThresholdRange_.setMinValue(vmm->getMin());
        binarizationThresholdRange_.setMaxValue(vmm->getMax());
        binarizationThresholdRange_.set(tgt::vec2((vmm->getMax() + vmm->getMin())/2));
    }
    if(!startPressed_) {
        return;
    }
    startPressed_ = false;
    const VolumeRAM* input = invol->getRepresentation<VolumeRAM>();
    tgtAssert(input, "No Volumeram");

    RealWorldMapping rwm = invol->getRealWorldMapping();
    float normalizedMinThreshold = rwm.realWorldToNormalized(binarizationThresholdRange_.get().x);
    float normalizedMaxThreshold = rwm.realWorldToNormalized(binarizationThresholdRange_.get().y);

    VolumeList* output = new VolumeList;
    for(int i = 0; i < steps_.get(); ++i) {
        float normalizedThreshold = normalizedMinThreshold + (normalizedMaxThreshold - normalizedMinThreshold) * i / (steps_.get() - 1);
        float rwThreshold = rwm.normalizedToRealWorld(normalizedThreshold);
        Volume* binVol = new Volume(applyThreshold(*input, normalizedThreshold), invol);
        // Hack: store the threshold used as offset in the RWM.
        RealWorldMapping rwmEncodedThreshold(1.0f, rwThreshold, "");
        binVol->setRealWorldMapping(rwmEncodedThreshold);
        output->add(binVol);
    }
    outport_.setData(output);
}
VoreenSerializableObject* VolumeMultiThreshold::create() const {
    return new VolumeMultiThreshold();
}

} // namespace voreen
