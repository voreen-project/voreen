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

#include "imagethreshold.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

#include "tgt/textureunit.h"

using tgt::TextureUnit;

namespace voreen {

ImageThreshold::ImageThreshold()
    : ImageProcessorBypassable("image/imagethreshold"),
    inport_(Port::INPORT, "image.inport", "Image Input"),
    outport_(Port::OUTPORT, "image.outport", "Image Output"),
    lowerThreshold_("lowerThreshold", "Lower Threshold", 0.f, 0.0f, 1.f),
    upperThreshold_("upperThreshold", "Upper Threshold", 1.f, 0.0f, 1.f),
    lowerMaskColor_("lowerMaskColor", "Lower Mask Color", tgt::Color(0.f,0.f,0.f,1.f)),
    upperMaskColor_("upperMaskColor", "Upper Mask Color", tgt::Color(0.f,0.f,0.f,1.f))

{
    addPort(inport_);
    addPort(outport_);

    addProperty(lowerThreshold_);
    ON_PROPERTY_CHANGE(lowerThreshold_,ImageThreshold,onLowerThresholdChange);
    addProperty(upperThreshold_);
    ON_PROPERTY_CHANGE(upperThreshold_,ImageThreshold,onUpperThresholdChange);
    addProperty(lowerMaskColor_);
    addProperty(upperMaskColor_);
}

Processor* ImageThreshold::create() const {
    return new ImageThreshold();
}

void ImageThreshold::process() {

    if (!enableSwitch_.get()){
        bypass(&inport_, &outport_);
        return;
    }

    outport_.activateTarget();
    outport_.clearTarget();

    TextureUnit colorUnit, depthUnit;
    inport_.bindTextures(colorUnit, depthUnit);

    // initialize shader
    program_->activate();
    setGlobalShaderParameters(program_);
    program_->setUniform("colorTex_", colorUnit.getUnitNumber());
    program_->setUniform("depthTex_", depthUnit.getUnitNumber());
    inport_.setTextureParameters(program_, "texParams_");
    program_->setUniform("lowerThreshold_", lowerThreshold_.get());
    program_->setUniform("upperThreshold_", upperThreshold_.get());
    program_->setUniform("lowerMaskColor_", lowerMaskColor_.get());
    program_->setUniform("upperMaskColor_", upperMaskColor_.get());

    renderQuad();

    program_->deactivate();
    outport_.deactivateTarget();
    TextureUnit::setZeroUnit();
    LGL_ERROR;
}

void ImageThreshold::onLowerThresholdChange() {
    if(lowerThreshold_.get() > upperThreshold_.get())
        upperThreshold_.set(lowerThreshold_.get());
}

void ImageThreshold::onUpperThresholdChange() {
    if(lowerThreshold_.get() > upperThreshold_.get())
        lowerThreshold_.set(upperThreshold_.get());
}

} // voreen namespace
