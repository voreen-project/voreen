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

#include "preintegrationtablerenderer.h"

#include "tgt/textureunit.h"
#include "tgt/immediatemode/immediatemode.h"
#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"
#include "voreen/core/datastructures/transfunc/1d/preintegrationtable.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
using tgt::TextureUnit;

namespace voreen {

PreIntegrationTableRenderer::PreIntegrationTableRenderer()
    : ImageProcessor("copyimage")
    , outport_(Port::OUTPORT, "image.out", "Image Output", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , transFunc_("transferFunction", "Transfer Function")
    , piMode_("pinmode", "Pre-Integration Table Mode")
    , showGrid_("showGrid", "Show Grid", false)
    , gridDivisions_("gridDivisions", "Grid Divisions", 10, 1, 20)
{
    addPort(outport_);

    addProperty(transFunc_);
    addProperty(piMode_);
    piMode_.addOption("pre-integrated-gpu", "Pre-Integrated TF (GPU)");
    piMode_.addOption("pre-integrated-cpu", "Pre-integrated TF (CPU)");
    piMode_.addOption("pre-integrated-approximation", "Pre-integrated TF (Approximation)");

    addProperty(showGrid_);
    showGrid_.setGroupID("grid");
    addProperty(gridDivisions_);
    gridDivisions_.setGroupID("grid");
    setPropertyGroupGuiName("grid", "Grid");

    ON_CHANGE_LAMBDA(showGrid_, [&]{
        gridDivisions_.setVisibleFlag(showGrid_.get());
    });
    showGrid_.invalidate();
}

Processor* PreIntegrationTableRenderer::create() const {
    return new PreIntegrationTableRenderer();
}

std::string PreIntegrationTableRenderer::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string header = ImageProcessor::generateHeader(version);
    header += "#define NO_DEPTH_TEX\n";
    return header;
}

bool PreIntegrationTableRenderer::isReady() const {
    return (isInitialized() && outport_.isReady());
}

void PreIntegrationTableRenderer::process() {

    tgtAssert(outport_.isReady(), "Outport not ready");
    tgtAssert(program_, "Shader missing");

    // get preintegration table
    const tgt::Texture* tex = 0;
    if (piMode_.getKey() == "pre-integrated-gpu")
        tex = transFunc_.get()->getPreIntegrationTable(1.f, 0, false, true)->getTexture();
    else if (piMode_.getKey() == "pre-integrated-cpu")
        tex = transFunc_.get()->getPreIntegrationTable(1.f, 0, false, false)->getTexture();
    else if (piMode_.getKey() == "pre-integrated-approximation")
        tex = transFunc_.get()->getPreIntegrationTable(1.f, 0, true, false)->getTexture();

    if (!tex) {
        LERROR("No pre-integration texture!");
        return;
    }

    // render preintegration table to quad
    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (showGrid_.get()){
        int divisions = gridDivisions_.get();
        float size = 2.0f/divisions;
        tgt::vec4 colors[2] = {tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f), tgt::vec4(220.0f, 220.0f, 220.0f, 255.0f)*(1.0f/255.0f)};
        IMode.begin(tgt::ImmediateMode::TRIANGLES);
        for(int i = 0; i != divisions; i++){
            for(int j = 0; j != divisions; j++){
                IMode.color(colors[(i+j)%2]);
                IMode.vertex(i*size-1.0f, j*size-1.0f);
                IMode.vertex((i+1)*size-1.0f, j*size-1.0f);
                IMode.vertex(i*size-1.0f, (j+1)*size-1.0f);

                IMode.vertex((i+1)*size-1.0f, j*size-1.0f);
                IMode.vertex(i*size-1.0f, (j+1)*size-1.0f);
                IMode.vertex((i+1)*size-1.0f, (j+1)*size-1.0f);
            }

        }
        IMode.end();
        IMode.color(tgt::vec4::one);

        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    }

    program_->activate();
    // set common uniforms
    setGlobalShaderParameters(program_);

    // bind pre-integration table to tex unit
    glActiveTexture(GL_TEXTURE0);
    tex->bind();
    LGL_ERROR;

    // pass texture parameters to the shader
    program_->setUniform("colorTex_", 0);
    program_->setIgnoreUniformLocationError(true);
    program_->setUniform("texParams_.dimensions_", tgt::vec2(tex->getDimensions().xy()));
    program_->setUniform("texParams_.dimensionsRCP_", tgt::vec2(1.f) / tgt::vec2(tex->getDimensions().xy()));
    program_->setUniform("texParams_.matrix_", tgt::mat4::identity);
    program_->setIgnoreUniformLocationError(false);
    renderQuad();
    program_->deactivate();
    if (showGrid_.get()){
        glDisable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ZERO);
    }
    outport_.deactivateTarget();
    LGL_ERROR;
}

} // namespace voreen
