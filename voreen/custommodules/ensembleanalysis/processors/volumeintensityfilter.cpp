/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
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

#include "volumeintensityfilter.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeoperator.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"

namespace voreen {

VolumeIntensityFilter::VolumeIntensityFilter()
    : CachingVolumeProcessor()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false)
    , enableProcessingProp_("enableProcessing", "Enable", true)
    , intensityRange_("intensityRange", "Selected Intensity-Range")
    , continousUpdate_("continousUpdates", "Update continuously", true, Processor::VALID)
    , updateButton_("update", "Update")
    , forceUpdate_(true)
{
    addPort(inport_);
    addPort(outport_);

    enableProcessingProp_.onChange(MemberFunctionCallback<VolumeIntensityFilter>(this, &VolumeIntensityFilter::forceUpdate));
    intensityRange_.onChange(MemberFunctionCallback<VolumeIntensityFilter>(this, &VolumeIntensityFilter::forceUpdate));
    updateButton_.onChange(MemberFunctionCallback<VolumeIntensityFilter>(this, &VolumeIntensityFilter::maskVolume));
    inport_.onChange(MemberFunctionCallback<VolumeIntensityFilter>(this, &VolumeIntensityFilter::initVolumeRangeProperty));

    addProperty(enableProcessingProp_);
    addProperty(intensityRange_);
    addProperty(continousUpdate_);
    addProperty(updateButton_);
}

VolumeIntensityFilter::~VolumeIntensityFilter() {}

Processor* VolumeIntensityFilter::create() const {
    return new VolumeIntensityFilter();
}

void VolumeIntensityFilter::initVolumeRangeProperty() {
    if(!inport_.hasData()) return;

    const VolumeMinMax* minMax = inport_.getData()->getDerivedData<VolumeMinMax>();
    tgtAssert(minMax, "min max null");

    intensityRange_.setMinValue(minMax->getMin());
    intensityRange_.setMaxValue(minMax->getMax());
    intensityRange_.set(tgt::vec2(minMax->getMin(), minMax->getMax()));
    intensityRange_.updateWidgets();
}

void VolumeIntensityFilter::process() {
    if (!enableProcessingProp_.get()) {
        outport_.setData(const_cast<VolumeBase*>(inport_.getData()), false);
    }
    else if (forceUpdate_ || inport_.hasChanged()) {
        if(continousUpdate_.get())
            maskVolume();
    }
}

// private methods
//

void VolumeIntensityFilter::forceUpdate() {
    forceUpdate_ = true;
}

void VolumeIntensityFilter::maskVolume() {
    tgtAssert(inport_.hasData(), "Inport has not data");

    forceUpdate_ = false;
    const VolumeRAM_Float* vol = dynamic_cast<const VolumeRAM_Float*>(inport_.getData()->getRepresentation<VolumeRAM>());

    if (vol) {
        VolumeRAM_Float* v = 0;
        v = vol->clone();

#ifdef VRN_MODULE_OPENMP
        #pragma omp parallel for
#endif
        for(size_t i = 0; i < vol->getNumVoxels(); i++) {
            float voxelValue = vol->voxel(i);
            if(voxelValue < intensityRange_.get().x || voxelValue > intensityRange_.get().y) {
                v->voxel(i) = 0.0f;
            }
        }

        Volume* vh = new Volume(v, inport_.getData());
        outport_.setData(vh);
    }
    else {
        outport_.setData(0);
    }
}

}   // namespace
