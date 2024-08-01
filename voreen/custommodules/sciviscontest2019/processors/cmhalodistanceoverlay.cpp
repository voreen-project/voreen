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

#include "cmhalodistanceoverlay.h"

#include "tgt/textureunit.h"
#include "tgt/event/mouseevent.h"
#include <sstream>

using tgt::TextureUnit;

namespace voreen {

CMHaloDistanceOverlay::CMHaloDistanceOverlay()
    : ImageProcessor("image/compositor")
    , imageInport_(Port::INPORT, "image.in", "Image Inpput")
    , inport_(Port::INPORT, "halohandle.in", "Halo Data Input")
    , outport_(Port::OUTPORT, "image.out", "Image Output")
    , camera_("camera", "Camera")
    , font_("font", "Font")
    , selectedHaloIDProp_("selectedHaloIDProp", "ID of focused halo", 257, -1, 1000000)
    , mouseOverHaloIDProp_("mouseOverHaloIDProp", "ID of hovered over halo", CMMergerTree::NO_HALO_ID, -1, 1000000, Processor::INVALID_RESULT, NumericProperty<int>::STATIC, Property::LOD_DEBUG)
    , shaderProp_("cmhdoshaderProp", "Shader", "cmhalodistance.frag", "cmhalodistance.vert", "cmhalodistance.geom" )
    , vbo_(0)
    , vao_(0)
{
    addPort(imageInport_);
    addPort(inport_);
    addPort(outport_);

    addProperty(font_);
    addProperty(camera_);
    addProperty(selectedHaloIDProp_);
    addProperty(mouseOverHaloIDProp_);
    addProperty(shaderProp_);

    ON_CHANGE(selectedHaloIDProp_, CMHaloDistanceOverlay, halosChanged);
    ON_CHANGE(mouseOverHaloIDProp_, CMHaloDistanceOverlay, halosChanged);
}

Processor* CMHaloDistanceOverlay::create() const {
    return new CMHaloDistanceOverlay();
}

void CMHaloDistanceOverlay::initialize() {
    ImageProcessor::initialize();
    setupBuffers();
    shaderProp_.rebuild();
}

void CMHaloDistanceOverlay::deinitialize() {
    ImageProcessor::deinitialize();
}
void CMHaloDistanceOverlay::setupBuffers() {
    if (vao_) glDeleteVertexArrays(1, &vao_);
    if (vbo_) glDeleteBuffers(1, &vbo_);

    glGenVertexArrays(1, &vao_);
    glGenBuffers(1, &vbo_);

    glBindVertexArray(vao_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferData(GL_ARRAY_BUFFER, 2*sizeof(tgt::vec3), NULL, GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(tgt::vec3), 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    LGL_ERROR;
}
void CMHaloDistanceOverlay::halosChanged() {
    if (!vao_ | !vbo_) {
        return;
    }
    std::vector<tgt::vec3> positions;
    if(!inport_.getData()) {
        return;
    }
    const CMHalo* selected = inport_.getData()->haloByID(selectedHaloIDProp_.get());
    const CMHalo* mouseOver = inport_.getData()->haloByID(mouseOverHaloIDProp_.get());
    if(!selected || !mouseOver || mouseOver == selected) {
        return;
    }
    //Update data
    positions.push_back(selected->pos);
    positions.push_back(mouseOver->pos);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferSubData(GL_ARRAY_BUFFER, 0, 2*sizeof(tgt::vec3), positions.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    LGL_ERROR;
}
static tgt::ivec2 toScreenPos(const tgt::mat4& projectionMatrix, const tgt::mat4& viewMatrix, const tgt::vec2& canvasSize, const tgt::vec3& pos) {
    tgt::vec4 centerPos = tgt::vec4(pos, 1);
    tgt::vec4 centerPosTransformed = projectionMatrix*viewMatrix*centerPos;
    return ((centerPosTransformed.xy()/centerPosTransformed.w*0.5f)+tgt::vec2(0.5,0.5))*canvasSize;
}

void CMHaloDistanceOverlay::process() {
    int lineWidth = 3.0f;
    if (!vao_ | !vbo_) {
        return;
    }
    //same code as in ImageOverlay
    //          |
    //          v
    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, imageInport_.getColorTexture()->getId());
    renderQuad();
    glBindTexture(GL_TEXTURE_2D, 0);
    //===========================

    //Lines
    int mID = mouseOverHaloIDProp_.get();
    int sID = selectedHaloIDProp_.get();
    if (mID == CMMergerTree::NO_HALO_ID || sID == CMMergerTree::NO_HALO_ID || mID == sID) {
        outport_.deactivateTarget();
        return;
    }

    glDisable(GL_DEPTH_TEST);
    glLineWidth(lineWidth);

    const tgt::mat4 projectionMatrix = camera_.get().getProjectionMatrix(outport_.getSize());
    const tgt::mat4 viewMatrix = camera_.get().getViewMatrix();
    tgt::Shader* shader = shaderProp_.getShader();
    shader->activate();
    shader->setUniform("projectionMatrix", projectionMatrix);
    shader->setUniform("viewMatrix", viewMatrix);

    glBindVertexArray(vao_);
    glDrawArrays(GL_LINES, 0, 2);
    glBindVertexArray(0);

    shader->deactivate();

    glLineWidth(1.0f);


    //Text
    const CMHalo* selected = inport_.getData()->haloByID(sID);
    const CMHalo* mouseOver = inport_.getData()->haloByID(mID);

	if (selected == nullptr || mouseOver == nullptr) {
		glEnable(GL_DEPTH_TEST);
		outport_.deactivateTarget();
		return;
	}

    tgt::ivec2 p1 = toScreenPos(projectionMatrix, viewMatrix, imageInport_.getSize(), selected->pos);
    tgt::ivec2 p2 = toScreenPos(projectionMatrix, viewMatrix, imageInport_.getSize(), mouseOver->pos);
    tgt::ivec2 p1p2 = p2-p1;
    tgt::ivec2 screenPos = (p1+p2)/2;

    std::stringstream ss;
    ss << std::fixed;
    ss << " ";
    ss << tgt::distance(selected->pos, mouseOver->pos);
    ss << " Mpc/h ";

    tgt::Font* f= font_.get();
    if(p1p2.x>=0) {
        f->setTextAlignment(tgt::Font::TextAlignment::MiddleRight);
        f->setLineWidth(screenPos.x);
        screenPos.x = 0;
    } else {
        f->setTextAlignment(tgt::Font::TextAlignment::MiddleLeft);
    }
    if(p1p2.y>=0) {
        f->setTextAlignment(tgt::Font::TextAlignment::TopCentered);
        screenPos.y += lineWidth;
    } else {
        f->setTextAlignment(tgt::Font::TextAlignment::BottomCentered);
        screenPos.y -= lineWidth;
    }
    //glBlendFuncSeparate(GL_ONE_MINUS_DST_COLOR, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    f->render(tgt::ivec3(screenPos,0), ss.str(), imageInport_.getSize());

    glEnable(GL_DEPTH_TEST);
    outport_.deactivateTarget();
    LGL_ERROR;
}

} // namespace voreen
