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

#include "voreen/core/processors/imageprocessor.h"

namespace voreen {

const std::string ImageProcessor::loggerCat_("voreen.ImageProcessor");

ImageProcessor::ImageProcessor(const std::string& shaderFilename, bool enableCoarsenessAdaptation)
    : RenderProcessor()
    , shaderProp_("shader.program", "Image shader", shaderFilename.empty() ? "" : shaderFilename + ".frag", shaderFilename.empty() ? "" : "passthrough.vert","",
      Processor::INVALID_PROGRAM, Property::LOD_ADVANCED)
    , interactionAdapt_("iamode.adapt", "Adapt to interaction coarseness", false)
    , interactionFactor_("iamode.factor", "Interaction coarseness", 1, 1, 16, Processor::VALID)
    , program_(0)
    , shaderFilename_(shaderFilename)
{
    addProperty(shaderProp_);
    if(enableCoarsenessAdaptation) {
        addProperty(interactionAdapt_);
        addProperty(interactionFactor_);
        interactionAdapt_.setGroupID("coarseness");
        interactionFactor_.setGroupID("coarseness");
        setPropertyGroupGuiName("coarseness", "Interaction coarseness adaptation");
    }
}

ImageProcessor::~ImageProcessor() {
}

void ImageProcessor::initialize() {
    RenderProcessor::initialize();
    compile();
}

void ImageProcessor::compile() {
    program_ = 0;
    if(!shaderFilename_.empty()) {
        shaderProp_.setHeader(generateHeader());
        shaderProp_.rebuild();
        if(shaderProp_.hasValidShader())
            program_ = shaderProp_.getShader();
    }
}

void ImageProcessor::beforeProcess() {
    RenderProcessor::beforeProcess();

    // compile program if needed
    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM)
        compile();
    LGL_ERROR;
}

/**
* Read back depth buffer and determine min and max depth value.
*
* @param port the port to analyze the depth buffer
* @return tgt::vec2 with x = minDepth, y = maxDepth
*/
tgt::vec2 ImageProcessor::computeDepthRange(RenderPort* port) {
    port->getDepthTexture()->downloadTexture();
    float* pixels = (float*)port->getDepthTexture()->getCpuTextureData();
    float curDepth = *(pixels);
    float minDepth = curDepth;
    float maxDepth = curDepth;
    int numPixels = tgt::hmul(port->getSize());
    for (int i = 0; i < numPixels; i++) {
        curDepth = *(pixels++);
        minDepth = std::min(minDepth, curDepth);
        maxDepth = std::max(maxDepth, curDepth);
    }
    return tgt::vec2(minDepth, maxDepth);
}

} // voreen namespace
