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

//header
#include "cosmologyarrowrenderer3d.h"
//texture handling
#include "tgt/texturemanager.h"
#include "tgt/textureunit.h"
#include "tgt/camera.h"
//shader support test
#include "tgt/gpucapabilities.h"
#include "voreen/core/utils/glsl.h"

#include "../datastructures/cmparticledata.h"
#include "../utils/cmmath.h"

#include <iostream> 
namespace voreen {

const std::string CosmologyArrowRenderer3D::loggerCat_("voreen.viscontest2015.CosmologyArrowRenderer3D");

CosmologyArrowRenderer3D::CosmologyArrowRenderer3D()
    : RenderProcessor()
    , inport_       (Port::INPORT, "particlehandle.output", "Particle Data Output")
    , renderOutport_(Port::OUTPORT, "renderOutport", "Arrow Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , timeStep_          ("timeStep",           "Time Step", 0, 0, 624)
	, useAlpha_			 ("useAlpha", "Use alpha", false)
	, alphaFactor_		 ("alphaFactor", "Particle alpha factor", 1.0f, 0.0f, 1.0f)
    , shaderProp_        ("arrow.prg",          "Arrow shader", "cmarrowshader.frag",  "cmarrowshader.vert",  "cmarrowshader.geom")
    , shaderPropFast_    ("arrow.prg2",          "Arrow shader", "cmarrowshader2.frag", "cmarrowshader2.vert", "cmarrowshader2.geom")
    , cameraProp_        ("cameraProp",         "Camera", tgt::Camera(tgt::vec3(0.0f, 0.0f, 50.f), tgt::vec3(0.0f, 0.0f, 0.0f), 
                                                          tgt::vec3(0.0f, 1.0f, 0.0f), 45.f,1.f,0.1f,500.f))
    , arrowColorMode_    ("arrowColorMode",     "Derive arrow color...")
    , arrowComponentMode_("arrowComponentMode", "Use from data...")
    , arrowSize_         ("arrowSize",          "Arrow Size",1.f,0.01f,300.f)
    , tfProp_            ("tfProp",             "Transfer Function")
    , xVoxelOffset_      ("xVoxelOffset",       "X Shift",0.5f,0.f,1.f)
    , yVoxelOffset_      ("yVoxelOffset",       "Y Shift",0.5f,0.f,1.f)
    , zVoxelOffset_      ("zVoxelOffset",       "Z Shift",0.5f,0.f,1.f)
    , reduceQuality_     ("reduceQuality",      "Reduce Quality", false)
    , adaptCamera_       ("adaptCamera",        "Should adapt camera", true)
    , universeDimensions_("universeDimensions", "Bounds in universe coords")
    , universeMatrix_    ("universeMatrix_",    "Matrix to real world", tgt::mat4::identity)
    , cameraHandler_(0)
    , vbo_(0)
{
    addPort(inport_); 
    addPort(renderOutport_); 

    addProperty(timeStep_);

    addProperty(arrowColorMode_);
        arrowColorMode_.addOption("fromtf",  "from Arrow Length",    CosmologyArrowRenderer3D::ARROW_LENGTH);
        arrowColorMode_.addOption("fromdir", "from Arrow Direction", CosmologyArrowRenderer3D::ARROW_DIRECTION);
        arrowColorMode_.setGroupID("arrow");
    addProperty(arrowSize_);
        arrowSize_.setGroupID("arrow");
    addProperty(tfProp_);
        tfProp_.setGroupID("arrow");
    addProperty(arrowComponentMode_);
        arrowComponentMode_.addOption("ac_all",    "all Dimensions",   CosmologyArrowRenderer3D::AC_ALL_DIM);
        arrowComponentMode_.addOption("ac_only_x", "only X Dimension", CosmologyArrowRenderer3D::AC_ONLY_X_DIM);
        arrowComponentMode_.addOption("ac_only_y", "only Y Dimension", CosmologyArrowRenderer3D::AC_ONLY_Y_DIM);
        arrowComponentMode_.addOption("ac_only_z", "only Z Dimension", CosmologyArrowRenderer3D::AC_ONLY_Z_DIM);
        arrowComponentMode_.setGroupID("arrow");
    setPropertyGroupGuiName("arrow","Arrow Color Settings");

	addProperty(useAlpha_);
	addProperty(alphaFactor_);

    addProperty(reduceQuality_);

    addProperty(xVoxelOffset_);
        xVoxelOffset_.setGroupID("pos");
    addProperty(yVoxelOffset_);
        yVoxelOffset_.setGroupID("pos");
    addProperty(zVoxelOffset_);
        zVoxelOffset_.setGroupID("pos");
    setPropertyGroupGuiName("pos","Arrow Position Settings");


    addProperty(cameraProp_);
    addProperty(shaderProp_);
    addProperty(shaderPropFast_);
    addProperty(adaptCamera_);

    addProperty(universeDimensions_);
    universeDimensions_.setReadOnlyFlag(true);
    universeDimensions_.setVisibleFlag(false);
    addProperty(universeMatrix_);
    universeMatrix_.setReadOnlyFlag(true);
    universeMatrix_.setVisibleFlag(false);

    universeMatrix_.setMaxValue(tgt::mat4(100000000.0f));
    universeMatrix_.setMinValue(tgt::mat4(-100000000.0f));

    cameraHandler_ = new CameraInteractionHandler("cameraHandler", "Camera Handler", &cameraProp_);
    addInteractionHandler(cameraHandler_);

    inport_.onNewData(MemberFunctionCallback<CosmologyArrowRenderer3D>(this, &CosmologyArrowRenderer3D::buildBuffers));
    ON_CHANGE(arrowColorMode_,CosmologyArrowRenderer3D,changePropertyVisibility);
    ON_CHANGE(arrowComponentMode_, CosmologyArrowRenderer3D, buildBuffers);
    ON_CHANGE(timeStep_, CosmologyArrowRenderer3D, buildBuffers);

}

CosmologyArrowRenderer3D::~CosmologyArrowRenderer3D() {
    delete cameraHandler_;
}

Processor* CosmologyArrowRenderer3D::create() const {
    return new CosmologyArrowRenderer3D();
}

void CosmologyArrowRenderer3D::initialize() {
    RenderProcessor::initialize();
    compile();
}

void CosmologyArrowRenderer3D::deinitialize() {
    RenderProcessor::deinitialize();
}

void CosmologyArrowRenderer3D::compile() {
    shaderProp_.setHeader(generateHeader());
    shaderProp_.rebuild();
    shaderPropFast_.setHeader(generateHeader());
    shaderPropFast_.rebuild();
}


void CosmologyArrowRenderer3D::buildBuffers() {
    if (!inport_.getData())
        return;

    CMParticleDataTimeSlice*  timeSlice         = inport_.getData()->sliceAtTimeStep(timeStep_.get());
    int                       numberOfParticles = timeSlice->getNumberOfParticles();
    std::vector<VertexLayout> verts;
    const CMParticle *        particles = timeSlice->startUsingParticles();
    maxLength_ = 0.0f;
    float timeStep_ = timeSlice->getTimeStep();
    
    const char * enabledState = inport_.getData()->getEnabledState();
    for(int i = 0; i != numberOfParticles; i++){
        CMParticle   particle     = particles[i];
        if (!enabledState[particle.ident])
            continue;
        VertexLayout vertexLayout;
        vertexLayout.pos = particle.pos;
        vertexLayout.dir = particle.vel;

        switch(arrowComponentMode_.getValue()) {
        case AC_ONLY_X_DIM:
            vertexLayout.dir.y = 0.f; vertexLayout.dir.z = 0.f;
            break;
        case AC_ONLY_Y_DIM:
            vertexLayout.dir.x = 0.f; vertexLayout.dir.z = 0.f;
            break;
        case AC_ONLY_Z_DIM:
            vertexLayout.dir.x = 0.f; vertexLayout.dir.y = 0.f;
            break;
        }

        verts.push_back(vertexLayout);
        maxLength_ = std::max(maxLength_, tgt::length(particle.vel));
    }
    timeSlice->finishedUsingParticles();

    glGenBuffers(1, &vbo_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(VertexLayout)*verts.size(), verts.data(), GL_STATIC_DRAW);
    vertexCount_ = (int)verts.size();

    tgt::Bounds bounds = timeSlice->getBounds();

    if (adaptCamera_.get())
        cameraProp_.adaptInteractionToScene(bounds);




    //bounds = CMtransformBounds(timeSlice->getNormalisationTransformation(), bounds);
    universeDimensions_.setMinValue(bounds.getLLF() - tgt::vec3::one);
    universeDimensions_.setMaxValue(bounds.getURB() + tgt::vec3::one);
    universeDimensions_.set(bounds);

    universeMatrix_.set(tgt::mat4::identity);
}


void CosmologyArrowRenderer3D::changePropertyVisibility() {
    tfProp_.setVisibleFlag(!arrowColorMode_.getValue());
}


static tgt::mat3 normalMatrix(tgt::mat4 M){
    // invert and transpose matrix
    tgt::mat4 Minv;
    M.invert(Minv);
    tgt::mat4 Mtinv = tgt::transpose(Minv);

    // get the 3x3 xyz matrix, because w is 0 for normals anyway
    tgt::mat3 Mtinv3;
    for(int x = 0; x != 3; x++){
        for(int y = 0; y != 3; y++){
            Mtinv3[x][y] = Mtinv[x][y];
        }
    }

    return Mtinv3;
}

//---------------------------------------------------------------------
//      process
//---------------------------------------------------------------------
void CosmologyArrowRenderer3D::process() {
    //activate outport
    renderOutport_.activateTarget();
    renderOutport_.clearTarget();

	bool       useAlpha = useAlpha_.get();

    //bind tftexture
    tgt::TextureUnit tfUnit;
    tgt::Texture* tex = tfProp_.get()->getTexture();
    tfUnit.activate();
    tex->bind();

    tgt::Camera cam = cameraProp_.get();
    tgt::mat4 projectionMatrix = cam.getProjectionMatrix(renderOutport_.getSize());
    tgt::mat4 normMatrix = inport_.getData()->sliceAtTimeStep(timeStep_.get())->getNormalisationTransformation();
    tgt::mat4 viewMatrix = cam.getViewMatrix();
    viewMatrix = viewMatrix*normMatrix;

    tgt::mat4 MVP = projectionMatrix*viewMatrix;

    tgt::mat3 MVPtinv3 = normalMatrix(MVP);

    tgt::Shader* program;
    if (reduceQuality_.get()){
        program = shaderPropFast_.getShader();
    }else{
        program = shaderProp_.getShader();
    }

    //activate shader and set uniforms
    program->activate();
	program->setUniform("alphaFactor_", alphaFactor_.get());
	program->setUniform("hight_",        arrowSize_.get());
	program->setUniform("MV_",           viewMatrix);
	program->setUniform("P_",            projectionMatrix);
	program->setUniform("MVtinv_",       normalMatrix(viewMatrix));

    if (!reduceQuality_.get()){
        program->setUniform("MVP_", MVP);
        program->setUniform("maxVelocity", maxLength_);
        program->setUniform("ColorTexture", tfUnit.getUnitNumber());
        program->setUniform("colorFromDir_", arrowColorMode_.getValue() == ARROW_LENGTH ? false : true);
        program->setUniform("Ptinv_", normalMatrix(projectionMatrix));
        program->setUniform("MVPtinv_", MVPtinv3);
    }

	if (useAlpha) {
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE);
		glDisable(GL_DEPTH_TEST);
	}
    
    // setup vertex buffers
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VertexLayout), (void*)offsetof(VertexLayout, pos));

    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(VertexLayout), (void*)offsetof(VertexLayout, dir));
    
    // Draw
    glDrawArrays(GL_POINTS, 0, vertexCount_);

    // clean up

	glDisable(GL_BLEND);
	glBlendFunc(GL_ONE, GL_ZERO);
	glEnable(GL_DEPTH_TEST);


    glBindVertexArray(0);
    glDeleteVertexArrays(1, &vao);
    program->deactivate();
    renderOutport_.deactivateTarget();
    LGL_ERROR;
}

} // namespace voreen
