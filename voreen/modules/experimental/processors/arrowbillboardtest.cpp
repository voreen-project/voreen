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

#include "arrowbillboardtest.h"


using tgt::TextureUnit;

namespace voreen {

ArrowBillboardTest::ArrowBillboardTest()
    : RenderProcessor()
    , outport_(Port::OUTPORT, "image.out", "Image out", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , camera_("camera", "Camera")
    , cylinderShader_("cylinderShader", "Cylinder shader", "arrowbillboard_cylinder.frag", "arrowbillboard_cylinder.vert", "arrowbillboard_cylinder.geom")
    , capShader_("cylinderCapShader", "CylinderCap shader", "arrowbillboard_cap.frag", "arrowbillboard_cap.vert", "arrowbillboard_cap.geom")
{
    addPort(outport_);
    addProperty(camera_);
    addProperty(cylinderShader_);
    addProperty(capShader_);

    interactionHandler_ = new CameraInteractionHandler("cameraHandler", "Camera Handler", &camera_);
    addInteractionHandler(interactionHandler_);
}

ArrowBillboardTest::~ArrowBillboardTest()
{
    delete interactionHandler_;
}

Processor* ArrowBillboardTest::create() const {
    return new ArrowBillboardTest();
}

void ArrowBillboardTest::initialize() {
    RenderProcessor::initialize();
    cylinderShader_.rebuild();
    capShader_.rebuild();
}

void ArrowBillboardTest::deinitialize() {
    RenderProcessor::deinitialize();
}

void ArrowBillboardTest::process() {

    outport_.activateTarget();

    std::cout<<outport_.getSize()<<std::endl;
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    GLuint vbo;
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    std::vector<Vertex> vs;
    Vertex v;
    v.pos = tgt::vec3(2, 0, 0);
    v.direction = tgt::vec3(-4, 0, 0);
    v.color = tgt::vec3(1, 0, 0);
    vs.push_back(v);
    v.pos = tgt::vec3(0, 2, 0);
    v.direction = tgt::vec3(0, -4, 0);
    v.color = tgt::vec3(0, 1, 0);
    vs.push_back(v);


    glBufferData(GL_ARRAY_BUFFER, sizeof(Vertex)*vs.size(), vs.data(), GL_DYNAMIC_DRAW);
    //glDisable(GL_DEPTH_TEST);

    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (GLvoid*)offsetof(Vertex, pos));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (GLvoid*)offsetof(Vertex, direction));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (GLvoid*)offsetof(Vertex, color));


    tgt::mat4 M = camera_.get().getViewMatrix();
    tgt::mat4 P = camera_.get().getProjectionMatrix(outport_.getSize());


    tgt::Shader *shader = cylinderShader_.getShader();
    shader->activate();


    shader->setUniform("M_", M);
    shader->setUniform("P_", P);

    glDrawArrays(GL_POINTS, 0, vs.size());
    shader->deactivate();



    shader = capShader_.getShader();
    shader->activate();
    shader->setUniform("M_", M);
    shader->setUniform("P_", P);

    glDrawArrays(GL_POINTS, 0, vs.size());
    shader->deactivate();



    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);


    outport_.deactivateTarget();
    LGL_ERROR;
}



} // namespace voreen
