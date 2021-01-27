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

#include "volumemasking.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeoperator.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"
#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"

namespace voreen {

VolumeMasking::VolumeMasking()
    : CachingVolumeProcessor()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false)
    , enableProcessingProp_("enableProcessing", "Enable")
    , maskFunction_("maskFunction", "Masking function")
    , passedVoxelAction_("passedVoxelAction", "Passed voxels")
    , continousUpdate_("continousUpdates", "Update continuously", true, Processor::VALID)
    , updateButton_("update", "Update")
    , forceUpdate_(true)
{
    addPort(inport_);
    addPort(outport_);

    passedVoxelAction_.addOption("passThrough", "Pass Through");
    passedVoxelAction_.addOption("maxIntensity", "Assign Max Intensity");
    passedVoxelAction_.addOption("alpha", "Assign Alpha Value");

    enableProcessingProp_.onChange(MemberFunctionCallback<VolumeMasking>(this, &VolumeMasking::forceUpdate));
    maskFunction_.onChange(MemberFunctionCallback<VolumeMasking>(this, &VolumeMasking::forceUpdate));
    passedVoxelAction_.onChange(MemberFunctionCallback<VolumeMasking>(this, &VolumeMasking::forceUpdate));
    updateButton_.onChange(MemberFunctionCallback<VolumeMasking>(this, &VolumeMasking::maskVolume));

    addProperty(enableProcessingProp_);
    addProperty(maskFunction_);
    addProperty(passedVoxelAction_);
    addProperty(continousUpdate_);
    addProperty(updateButton_);
}

VolumeMasking::~VolumeMasking() {}

Processor* VolumeMasking::create() const {
    return new VolumeMasking();
}

void VolumeMasking::beforeProcess() {
    maskFunction_.setVolume(inport_.getData());
    LGL_ERROR;
}

void VolumeMasking::process() {
    LGL_ERROR;

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

void VolumeMasking::forceUpdate() {
    forceUpdate_ = true;
}

void VolumeMasking::maskVolume() {
    tgtAssert(inport_.hasData(), "Inport has not data");

    forceUpdate_ = false;
    const VolumeRAM* vol = inport_.getData()->getRepresentation<VolumeRAM>();

    if (vol) {
        VolumeRAM* v = 0;
        v = vol->clone();

        tgt::Texture* maskTexture = maskFunction_.get()->getTexture();
        maskTexture->downloadTexture();
        const int maskTexDim = maskTexture->getDimensions().x;

        RealWorldMapping rwm = inport_.getData()->getRealWorldMapping();
        // apply masking
        if (passedVoxelAction_.isSelected("maxIntensity")) {
            for (size_t i=0; i<v->getNumVoxels(); i++) {
                float intensity = vol->getVoxelNormalized(i);
                intensity = rwm.normalizedToRealWorld(intensity);
                intensity = maskFunction_.get()->realWorldToNormalized(intensity);
                int index = static_cast<int>(intensity*(maskTexDim-1));
                index = tgt::clamp(index, 0, maskTexDim-1);

                int alpha = maskTexture->texel< tgt::Vector4<uint8_t> >(static_cast<size_t>(index)).a;
                v->setVoxelNormalized(alpha == 0 ? 0.f : 1.f, i);
                //progressBar_->setProgress(static_cast<float>(i) / static_cast<float>(v->getNumVoxels()));
            }
        }
        else if (passedVoxelAction_.isSelected("passThrough")) {
            for (size_t i=0; i<v->getNumVoxels(); i++) {
                float intensity = vol->getVoxelNormalized(i);
                intensity = rwm.normalizedToRealWorld(intensity);
                intensity = maskFunction_.get()->realWorldToNormalized(intensity);
                int index = static_cast<int>(intensity*(maskTexDim-1));
                index = tgt::clamp(index, 0, maskTexDim-1);

                int alpha = maskTexture->texel< tgt::Vector4<uint8_t> >(static_cast<size_t>(index)).a;
                if (alpha == 0)
                    v->setVoxelNormalized(0.f, i);
            }
        }
        else if (passedVoxelAction_.isSelected("alpha")) {
            for (size_t i=0; i<v->getNumVoxels(); i++) {
                float intensity = vol->getVoxelNormalized(i);
                intensity = rwm.normalizedToRealWorld(intensity);
                intensity = maskFunction_.get()->realWorldToNormalized(intensity);
                int index = static_cast<int>(intensity*(maskTexDim-1));
                index = tgt::clamp(index, 0, maskTexDim-1);

                int alpha = maskTexture->texel< tgt::Vector4<uint8_t> >(static_cast<size_t>(index)).a;
                v->setVoxelNormalized((float)alpha / 255.0f, i);
            }
        }
        else {
            LWARNING("Unknown voxel action: " << passedVoxelAction_.get());
        }

        Volume* vh = new Volume(v, inport_.getData());
        if (passedVoxelAction_.isSelected("alpha"))
            vh->setRealWorldMapping(RealWorldMapping());
        outport_.setData(vh);
    }
    else {
        outport_.setData(0);
    }
}

}   // namespace
