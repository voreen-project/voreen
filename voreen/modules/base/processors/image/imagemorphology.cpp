/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "imagemorphology.h"

#include "voreen/core/utils/voreenqualitymode.h"

#include "tgt/textureunit.h"

using tgt::TextureUnit;

namespace voreen {

ImageMorphology::ImageMorphology()
    : ImageProcessorBypassable("image/imagemorphology", true),
      inport_(Port::INPORT, "image.inport", "Image Input"),
      outport_(Port::OUTPORT, "image.outport", "Image Output"),
      kernelRadius_("kernelRadius", "Kernel Radius", 1, 1, 50),
      modeProp_("effectModeAsString", "Mode"),
      shapeProp_("shapeModeAsString", "Shape")
{
    modeProp_.addOption("dilation", "Dilation");
    modeProp_.addOption("erosion", "Erosion");
    modeProp_.select("dilation");

    shapeProp_.addOption("square", "Square");
    shapeProp_.addOption("sphere", "Sphere");
    shapeProp_.select("square");

    addProperty(modeProp_);
    addProperty(shapeProp_);
    addProperty(kernelRadius_);

    addPort(inport_);
    addPort(outport_);
}

Processor* ImageMorphology::create() const {
    return new ImageMorphology();
}

void ImageMorphology::process() {

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

    //get mode and pass it to shader
    int mode;

    if (modeProp_.isSelected("dilation"))
        mode = 0;
    else
        mode = 1;

    program_->setUniform("mode_", mode);

    //get shape and pass it to shader
    int shape;

    if (shapeProp_.isSelected("square"))
        shape = 0;
    else
        shape = 1;

    program_->setUniform("shape_", shape);

    //get Kernel Radius and pass it to the shader
    int kernelRadius = kernelRadius_.get();
    if(QualityMode.isInteractionMode() && interactionAdapt_.get())
        kernelRadius = std::max(1, (int)(kernelRadius * (1.f / interactionFactor_.get())));
    program_->setUniform("kernelRadius_", kernelRadius);

    inport_.setTextureParameters(program_, "textureParameters_");

    renderQuad();

    program_->deactivate();
    outport_.deactivateTarget();
    TextureUnit::setZeroUnit();
    LGL_ERROR;
}

} // voreen namespace
