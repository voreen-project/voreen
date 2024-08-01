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

//header file
#include "cosmologytracingpath.h"
//needed headers (used in process())
#include "tgt/textureunit.h"
#include "tgt/textureunit.h"
#include "../utils/cmmath.h"
#include <assert.h>
#include <random>

//we are in namespace voreen
namespace voreen {

CosmologyTracingPath::CosmologyTracingPath()
    : RenderProcessor()
    , outport_(Port::OUTPORT, "outport","Modified Image", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , inport_(Port::INPORT, "particlehandle.output", "Particle Data Output")
    , camera_("camera", "Camera")
    , shaderProp_("shaderProp", "ASDSAF", "cmtracingpath.frag", "cmtracingpath.vert")
    , adaptCamera_("adaptCamera", "Should adapt camera", true)
	, traceJumps_("traceJumps", "Trace Periodoc Particles", false)
    , transFunc_("transFunc", "Transfer function")
    , timeStep_("timeStep", "Time Step", 0.0f, 0.0f, 624.0f)
    , useAlpha_         ("useAlpha",           "Use alpha", false)
    , alphaFactor_      ("alphaFactor",        "Alpha Factor", 1.0f, 0.0f, 1.0f)
    , universeDimensions_("universeDimensions", "Bounds in universe coords")
    , universeMatrix_    ("universeMatrix_",    "Matrix to real world", tgt::mat4::identity)
	, baryonColor ("baryonColor", "Color of Baryons")
	, darkMatterColor("darkMatterColor", "Color of Darkmatter")
	, starColor("starColor", "Color of Star Particles")
	, windColor("windColor", "Color of Wind Particles")
	, gasColor("gasColor", "Color of Gas Particles")
	, agnColor("agnColor", "Color of AGN Particles")
    , vao_(0)
{

    cameraHandler_ = new CameraInteractionHandler("cameraHandler", "Camera", &camera_);
    addInteractionHandler(cameraHandler_);
    addPort(outport_);
    addPort(inport_);

    addProperty(transFunc_);
    addProperty(timeStep_);
    addProperty(camera_);
    addProperty(shaderProp_);
    addProperty(useAlpha_);
    addProperty(alphaFactor_);
    addProperty(adaptCamera_);
	addProperty(traceJumps_);
    addProperty(universeDimensions_);
	addProperty(baryonColor);
	addProperty(darkMatterColor);
	addProperty(starColor);
	addProperty(windColor);
	addProperty(gasColor);
	addProperty(agnColor);
    addProperty(universeMatrix_);

    universeMatrix_.setMaxValue(tgt::mat4(100000000.0f));
    universeMatrix_.setMinValue(tgt::mat4(-100000000.0f));

    inport_.onNewData(MemberFunctionCallback<CosmologyTracingPath>(this, &CosmologyTracingPath::invalidateBuffers));
    //ON_CHANGE(timeStep_, CosmologyTracingPath, invalidateBuffers);
	ON_CHANGE(traceJumps_, CosmologyTracingPath, invalidateBuffers);
    buffersInvalid_ = true;
}

CosmologyTracingPath::~CosmologyTracingPath(){
    delete cameraHandler_;
}

void CosmologyTracingPath::initialize() {
    RenderProcessor::initialize();

    shaderProp_.rebuild();
}

void CosmologyTracingPath::deinitialize() {
    RenderProcessor::deinitialize();
}

void CosmologyTracingPath::invalidateBuffers(){
    buffersInvalid_ = true;
}

namespace{
struct VertexLayout{
    tgt::vec3 pos;
    float time;
    GLint matrix;
	GLint type;
};
}

void CosmologyTracingPath::setupBuffers(){
    if (vao_){ 
        glDeleteVertexArrays(1, &vao_);
        glDeleteBuffers(1, &vbo_);
        glDeleteBuffers(1, &ibo_);
        vao_ = 0;
    }
    if(!inport_.getData()) {
        return;
    }

    const std::vector<CMParticleDataTimeSlice*>& timeslices = inport_.getData()->particleDataTimeSlices();

    /*if (inport_.getData()->isFiltered()){
        LERRORC("voreen.CosmologyTracingPath", "Can't use filtered ParticleData! Use a CosmologyParticleCompator before this processor");
        return;
    }*/


    std::vector<const CMParticle*> particles;
    std::vector<const int*> remappings;
    std::vector<float> oneByTimeStep;
    std::vector<float> matrixId;
    
    normMatricies_.clear();
    timeSliceBounds_.clear();

    //tgt::vec2 interval = timeStep_.get();

    for(auto slice: timeslices){
        float timeStep = slice->getTimeStep();
        //if (timeStep >= interval.x && timeStep <= interval.y){
            remappings.push_back(slice->getRemappingTable());
            particles.push_back(slice->startUsingParticles());
            oneByTimeStep.push_back(1.0f/slice->getTimeStep());

            matrixId.push_back(normMatricies_.size());
            normMatricies_.push_back(slice->getNormalisationTransformation());
            timeSliceBounds_.push_back(tgt::Bounds());
            //normMatricies_.push_back(tgt::mat4::identity);
            
        //}
        
    }
    //normMatricies_.resize(101);

    glGenVertexArrays(1, &vao_);
    glGenBuffers(1, &vbo_);
    glGenBuffers(1, &ibo_);

    std::vector<VertexLayout> verts;
    std::vector<int> indicies;
    
	const char* enableStates = inport_.getData()->getEnabledState();

    //tgt::Bounds bounds;
    
    int numOfParticles = timeslices[0]->getNumberOfParticles();

    for(int i = 0; i < numOfParticles; i++){
		bool crossedBoundary = false;
		std::vector<VertexLayout> tmpVerts;
		std::vector<int> tmpIndicies;
		tgt::vec3 lastPosition;
		
		if (!enableStates[i]) {
			continue;
		}

        for(int t = 0; t < timeslices.size(); t++){
           
            CMParticle p = particles[t][remappings[t][i]];
            //p.pos *= oneByTimeStep[t];
            p.pos = CMtransform(normMatricies_[t], p.pos);
			if (t == 0) {
				lastPosition = p.pos;
			}
			
			else {
				//particle jumps
				if (!traceJumps_.get() && tgt::lengthSq(lastPosition - p.pos) > 1024.0f){
					crossedBoundary = true;
					break;
				}
			}
            VertexLayout v;
            v.pos = p.pos;
            v.time = 1.0f/oneByTimeStep[t];
            v.matrix = matrixId[t];
			switch (p.mask & (256 + 128 + 64 + 32 + 2)) {
			case 0: 
				//dark matter
				v.type = 1;
				break;
			case 2: 
				//baryons
				v.type = 0;
				break;
			case 34: 
				//star
				v.type = 2;
				break;
			case 66: 
				//wind
				v.type = 3;
				break;
			case 130: 
				//gas
				v.type = 4;
				break;
			case 256: 
				//agn
				v.type = 5;
				break;
			}

			tmpIndicies.push_back((int)verts.size() + (int)tmpVerts.size());
			tmpVerts.push_back(v);
            //bounds.addPoint(p.pos);
            timeSliceBounds_[t].addPoint(p.pos);
        }
		if (!crossedBoundary || traceJumps_.get()) {
			verts.insert(verts.end(), tmpVerts.begin(), tmpVerts.end());
			indicies.insert(indicies.end(), tmpIndicies.begin(), tmpIndicies.end());
			indicies.push_back(0x7fffffff);
		}
    }

    for(auto slice: timeslices){
        slice->finishedUsingParticles();
    }
    
    tgt::Bounds bounds;
    for(auto b: timeSliceBounds_){
        bounds.addVolume(b);
    }

    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(VertexLayout)*verts.size(), verts.data(), GL_STATIC_DRAW);
    assert(glGetError() == GL_NO_ERROR);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int)*indicies.size(), indicies.data(), GL_STATIC_DRAW);
    assert(glGetError() == GL_NO_ERROR);

    drawCount_ = (int)indicies.size();
    if (adaptCamera_.get())
        camera_.adaptInteractionToScene(bounds);
    //transFunc_.get()->setDomain(tgt::vec2(0.0f, 1.0f));

    setupBounds();

	buffersInvalid_ = false;
}


void CosmologyTracingPath::setupBounds(){
    if(!inport_.getData()) {
        return;
    }

    const std::vector<CMParticleDataTimeSlice*>& timeslices = inport_.getData()->particleDataTimeSlices();

    tgtAssert(timeslices.size() == timeSliceBounds_.size(), "Sizes need to match");

    tgt::Bounds bounds;
    tgt::vec2 interval = timeStep_.get();
    for(int i = 0; i != timeslices.size(); i++){
        float timeStep = timeslices[i]->getTimeStep();
        if (timeStep >= interval.x && timeStep <= interval.y){
            bounds.addVolume(timeSliceBounds_[i]);
        }
    }

    universeDimensions_.setMinValue(bounds.getLLF() -tgt::vec3::one);
    universeDimensions_.setMaxValue(bounds.getURB() + tgt::vec3::one);
    universeDimensions_.set(bounds);

    universeMatrix_.set(tgt::mat4::identity);
}

void CosmologyTracingPath::process() {
    if (buffersInvalid_)
        setupBuffers();
    
    if (!vao_ || ! vbo_) return;

    
    outport_.activateTarget();
    outport_.clearTarget();
    LGL_ERROR;
    tgt::Shader* shader;
    shader = shaderProp_.getShader();
    
    shader->activate();
	//shader->rebuild();

    tgt::Camera cam = camera_.get();
    tgt::mat4 projectionMatrix = cam.getProjectionMatrix(outport_.getSize());
    tgt::mat4 viewMatrix = cam.getViewMatrix();

	tgt::vec4 typeColors[6];
	typeColors[0] = baryonColor.get();
	typeColors[1] = darkMatterColor.get();
	typeColors[2] = starColor.get();
	typeColors[3] = windColor.get();
	typeColors[4] = gasColor.get();
	typeColors[5] = agnColor.get();

    
    shader->setUniform("viewProjectionMatrix_", projectionMatrix*viewMatrix);
    shader->setUniform("interval_", timeStep_.get());
    shader->setUniform("alphaFactor_", alphaFactor_.get());
	shader->setUniform("typeColors_", typeColors, 6);
    //shader->setUniform("normMats_", normMatricies_.data(), normMatricies_.size());
    glUniformMatrix4fv(glGetUniformLocation(shader->getID(), "normMats_"), normMatricies_.size(), GL_TRUE, &(normMatricies_[0][0][0]));
    glBindVertexArray(vao_);
    assert(glGetError() == GL_NO_ERROR);
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
	glEnableVertexAttribArray(3);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VertexLayout), (void*)offsetof(VertexLayout, pos));
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(VertexLayout), (void*)offsetof(VertexLayout, time));
    glVertexAttribIPointer(2, 1, GL_INT, sizeof(VertexLayout), (void*)offsetof(VertexLayout, matrix));
	glVertexAttribIPointer(3, 1, GL_INT, sizeof(VertexLayout), (void*)offsetof(VertexLayout, type));
    //assert(glGetError() == GL_NO_ERROR);
    LGL_ERROR;
    
    transFunc_.get()->setDomain(timeStep_.get());
    transFunc_.invalidate();

    //std::cerr<<"Tracing path: " << drawCount_ << std::endl;
    tgt::TextureUnit transferUnit;
    if (transFunc_.get()) {
        transferUnit.activate();
        transFunc_.get()->getTexture()->bind();
        transFunc_.get()->setUniform(shader, "transferFunc_", "transferFuncTex_", transferUnit.getUnitNumber());
            
    }
    LGL_ERROR;

    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_);
    assert(glGetError() == GL_NO_ERROR);
    glPrimitiveRestartIndex(0x7fffffff);
    assert(glGetError() == GL_NO_ERROR);
    glEnable(GL_PRIMITIVE_RESTART);
    glFinish();
    if (useAlpha_.get()){
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE);
        glDisable(GL_DEPTH_TEST);
    }

	glDrawElements(GL_LINE_STRIP, drawCount_,  GL_UNSIGNED_INT, 0);
    glDisable(GL_PRIMITIVE_RESTART);
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glBlendFunc(GL_ONE, GL_ZERO);
    //printf("gl err: %x\n", glGetError());
    //assert(glGetError() == GL_NO_ERROR);

    // cleanup
    shader->deactivate();
    outport_.deactivateTarget();
    tgt::TextureUnit::setZeroUnit();
    LGL_ERROR;
}

} // namespace

