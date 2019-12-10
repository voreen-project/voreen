/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "canny.h"

#include "tgt/textureunit.h"

using tgt::TextureUnit;

namespace voreen {

Canny::Canny()
    : ImageProcessor("canny")
    , startThreshold_("startThreshold", "Start threshold", 1.0f, 0.0f, 10.0f)
    , runThreshold_("runThreshold", "Run threshold", 0.5f, 0.0f, 10.0f)
    , edgeColor_("edgeColor", "Edge color", tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f))
    , inport_(Port::INPORT, "image.inport")
    , outport_(Port::OUTPORT, "image.outport", "image.outport", true)
    , privatePort_(Port::OUTPORT, "image.privateport", "image.privateport", true)
{
    addProperty(startThreshold_);
    addProperty(runThreshold_);
    addProperty(edgeColor_);

    addPort(inport_);
    addPort(outport_);
    addPrivateRenderPort(&privatePort_);
}

Processor* Canny::create() const {
    return new Canny();
}

void Canny::process() {
    privatePort_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    TextureUnit colorUnit, depthUnit;
    inport_.bindTextures(colorUnit.getEnum(), depthUnit.getEnum());

    // initialize shader
    // the shader performs the non-maximum suppression
    program_->activate();
    setGlobalShaderParameters(program_);
    program_->setUniform("colorTex_", colorUnit.getUnitNumber());
    program_->setUniform("depthTex_", depthUnit.getUnitNumber());
    inport_.setTextureParameters(program_, "texParams_");
    program_->setUniform("startThreshold_", startThreshold_.get());

    renderQuad();

    program_->deactivate();
    LGL_ERROR;

    // read back texture and trace edges
    privatePort_.getColorTexture()->downloadTexture();
    tgt::Texture* tex = privatePort_.getColorTexture();
    if(!tex) return;
    bool* processed = new bool[privatePort_.getSize().x*privatePort_.getSize().y];
    for (int x = 0; x < privatePort_.getSize().x; ++x)
        for (int y = 0; y < privatePort_.getSize().y; ++y)
            processed[y*privatePort_.getSize().x+x] = false;
    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    tgt::ivec2 dimensions = privatePort_.getSize();
    MatStack.loadMatrix(tgt::mat4::createOrtho(0.0, dimensions.x, 0.0, dimensions.y, 0.0, 1.0));
    for (int x = 0; x < dimensions.x; ++x) {
        for (int y = 0; y < dimensions.y; ++y) {
            if (!processed[y*dimensions.x+x]) {
                tgt::Color curTexel = tex->texelAsFloat(x,y);
                if (curTexel.b > 0.0) {
                    processed[y*dimensions.x+x] = true;
                    // mark the current pixel
                    glBegin(GL_POINTS);
                    glColor4fv(edgeColor_.get().elem);
                    glVertex2i(x,y);
                    glEnd();
                    // forward trace the edge following the next vector
                    tgt::ivec2 nextCoords = tgt::ivec2(x,y) + tgt::ivec2((curTexel.xy()-tgt::vec2(0.5))*tgt::vec2(2.0));
                    while (nextCoords.x >= 0 && nextCoords.x < dimensions.x &&
                           nextCoords.y >= 0 && nextCoords.y < dimensions.y &&
                           !processed[nextCoords.y*dimensions.x+nextCoords.x]) {
                        processed[nextCoords.y*dimensions.x+nextCoords.x] = true;
                        tgt::Color nextTexel = tex->texelAsFloat(nextCoords.x,nextCoords.y);
                        if (nextTexel.b > 0.0 && nextTexel.a > runThreshold_.get()) {
                            // mark the neighbor pixel
                            glBegin(GL_POINTS);
                                glColor4fv(edgeColor_.get().elem);
                                glVertex2iv(nextCoords.elem);
                            glEnd();
                        }
                        nextCoords = nextCoords + tgt::ivec2((nextTexel.xy()-tgt::vec2(0.5))*tgt::vec2(2.0));
                    }
                    // backward trace the edge following the next vector
                    nextCoords = tgt::ivec2(x,y) - tgt::ivec2((curTexel.xy()-tgt::vec2(0.5))*tgt::vec2(2.0));
                    while (nextCoords.x >= 0 && nextCoords.x < dimensions.x &&
                           nextCoords.y >= 0 && nextCoords.y < dimensions.y &&
                           !processed[nextCoords.y*dimensions.x+nextCoords.x]) {
                        processed[nextCoords.y*dimensions.x+nextCoords.x] = true;
                        tgt::Color nextTexel = tex->texelAsFloat(nextCoords.x,nextCoords.y);
                        if (nextTexel.b > 0.0 && nextTexel.a > runThreshold_.get()) {
                            // mark the neighbor pixel
                            glBegin(GL_POINTS);
                                glColor4fv(edgeColor_.get().elem);
                                glVertex2iv(nextCoords.elem);
                            glEnd();
                        }
                        nextCoords = nextCoords - tgt::ivec2((nextTexel.xy()-tgt::vec2(0.5))*tgt::vec2(2.0));
                    }
                }
            }
        }
    }
    MatStack.popMatrix();
    outport_.deactivateTarget();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    privatePort_.getColorTexture()->setCpuTextureData(0,false);
    delete[] processed;
}

} // voreen namespace
