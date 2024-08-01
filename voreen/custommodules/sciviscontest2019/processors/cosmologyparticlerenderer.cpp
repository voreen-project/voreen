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

#include "cosmologyparticlerenderer.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "tgt/textureunit.h"
#include "tgt/textureunit.h"
#include "modules/staging/utils/simdraycaster/jobqueue.h"
#include "../utils/cmmath.h"

#include <assert.h>
#include <random>

namespace voreen {

CosmologyParticleRenderer::CosmologyParticleRenderer()
    : RenderProcessor()
    , outport_(Port::OUTPORT, "outport","Modified Image", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , inport_(Port::INPORT,   "particlehandle.output", "Particle Data Output")
    , useRelativeRadius_  ("useRelativeRadius_",  "Use relative radius")
    , radiusProp_         ("radiusProp",          "Radius" ,0.25f ,0.0f, 4.0f)                     ///< min and max value
    , relativRadiusProp_  ("relativRadiusProp_",  "Relative radius", 0.5f, 0.00f, 1.0f)
    , camera_             ("camera",              "Camera")
    , shaderProp_         ("shaderProp",          "Shader for scalar properties", "particlerenderer.frag", "particlerenderer.vert", "particlerenderer.geom" )
    , shaderPropDVC_      ("shaderPropDVC",       "Shader for vector properties", "particlerenderer_vec.frag", "particlerenderer_vec.vert", "particlerenderer_vec.geom" )
	, shaderPropPTP_	  ("shaderPropPTP",		  "Shader for vector properties", "particlerenderer_ptp.frag", "particlerenderer_ptp.vert", "particlerenderer_ptp.geom")
    , shaderPropMASS_     ("shaderPropMASS",      "Shader for vector properties", "particlerenderer_mass.frag", "particlerenderer_mass.vert", "particlerenderer_mass.geom")
    , transFunc_          ("transFunc",           "Transfer function")
    , massTransFunc_      ("massTransFunc",       "Pseudo-TF for density")
    , renderMode_         ("renderMode",          "Mode for rendering")
    , useAlpha_           ("useAlpha",            "Use alpha", false)
    , alphaFactor_        ("alphaFactor",         "Particle alpha factor", 1.0f, 0.0f, 1.0f)
    , massColor_          ("massColor",           "Color for mass visulations")
    , timeStep_           ("timeStep",            "Time Step", 0.0f, 0.0f, 624.0f)
    , parallelizeProp_    ("parallelizeProp_",    "Parallel buffer creation", true)
    , adaptCamera_        ("adaptCamera",         "Should adapt camera",      true)
    , universeDimensions_ ("universeDimensions",  "Bounds in universe coords")
    , universeMatrix_     ("universeMatrix_",     "Matrix to real world", tgt::mat4::identity)
    , visualizesDirection_("visualizesDirection", "Is visualizing a direction with colors")
    , usesTf_             ("usesTf",              "Is using transfunc for visualisation")
    , visualizesMass_     ("visualizesMass_",     "Is visualizing mass")
	, unitDisplayed_	  ("unit",				  "Unit Displayed")
    , vao_(0)
    , velValues_(nullptr)
    , phiValues_(nullptr)
	//, rhoValues_(nullptr)
	, muValues_(nullptr)
	, uuValues_(nullptr)
	//, hhValues_(nullptr)
{

    cameraHandler_ = new CameraInteractionHandler("cameraHandler", "Camera", &camera_, true, true);
	camera_.setPosition(tgt::vec3(100.f, 100.f, 100.f));
    buffersInvalid_ = true;
    addInteractionHandler(cameraHandler_);
    addPort(outport_);
    addPort(inport_);

    addProperty(timeStep_);

    addProperty(renderMode_);
    addProperty(useRelativeRadius_);
    addProperty(radiusProp_);
    addProperty(relativRadiusProp_);
    relativRadiusProp_.setStepping(0.0001f);
    addProperty(transFunc_);
    addProperty(massTransFunc_);

    addProperty(useAlpha_);
    addProperty(alphaFactor_);

    addProperty(massColor_);

    addProperty(camera_);
    addProperty(shaderProp_);
    addProperty(shaderPropDVC_);
    addProperty(shaderPropMASS_);

    addProperty(parallelizeProp_);
    addProperty(adaptCamera_);

    addProperty(unitDisplayed_);
		unitDisplayed_.setReadOnlyFlag(true);

    addProperty(universeDimensions_);
        universeDimensions_.setReadOnlyFlag(true);
        universeDimensions_.setVisibleFlag(false);
    addProperty(universeMatrix_);
        universeMatrix_.setReadOnlyFlag(true);
        universeMatrix_.setVisibleFlag(false);
    addProperty(visualizesDirection_);
        visualizesDirection_.setReadOnlyFlag(true);
        visualizesDirection_.setVisibleFlag(false);
    addProperty(usesTf_);
        usesTf_.setReadOnlyFlag(true);
        usesTf_.setVisibleFlag(false);
    addProperty(visualizesMass_);
        visualizesMass_.setReadOnlyFlag(true);
        visualizesMass_.setVisibleFlag(false);
    universeMatrix_.setMaxValue(tgt::mat4(100000000.0f));
    universeMatrix_.setMinValue(tgt::mat4(-100000000.0f));

    renderMode_.addOption("km/s", "Velocity transfer Func",      VELOCITY_TRANSFER_FUNC);
    renderMode_.addOption("directVelocityColor",  "Direct velocity Color" ,      DIRECT_VELOCITY_COLOR);
    renderMode_.addOption("phi",      "Potential Transfer Func" ,    PHI_TRANSFER_FUNC);
	renderMode_.addOption("uu", "Internal Energy Transfer Func", INT_EGY_TRANSFER_FUNC);
	renderMode_.addOption("mu", "Molecular Weight Transfer Func", MOL_WGT_TRANSFER_FUNC);
    renderMode_.addOption("mass",                 "Mass",                        MASS);
	//renderMode_.addOption("rho", "SPH Density", SPH_DEN_TRANSFER_FUNC);
	//renderMode_.addOption("hh", "SPH Softening Length", SPH_LEN_TRANSFER_FUNC);
	//renderMode_.addOption("particletype", "Particle Type", PARTICLE_TYPE);

    useRelativeRadius_.addOption("absolute", "Absolute radius", false);
    useRelativeRadius_.addOption("relative", "Relative radius", true);

    ON_CHANGE(renderMode_,        CosmologyParticleRenderer,   changedRenderMode);
    ON_CHANGE(useAlpha_,          CosmologyParticleRenderer,   changeUseAlpha);
    ON_CHANGE(useRelativeRadius_, CosmologyParticleRenderer,   changeUseRelativeRadius);
    ON_CHANGE(alphaFactor_,       CosmologyParticleRenderer,   changedAlphaFactor);
    ON_CHANGE(timeStep_,          CosmologyParticleRenderer,   invalidateBuffers);
    inport_.onNewData(MemberFunctionCallback<CosmologyParticleRenderer>(this, &CosmologyParticleRenderer::invalidateBuffers));

    changeUseAlpha();
    changeUseRelativeRadius();
}

const std::string CosmologyParticleRenderer::loggerCat_("voreen.CMParticleRenderer");

CosmologyParticleRenderer::~CosmologyParticleRenderer(){
    delete cameraHandler_;
}

void CosmologyParticleRenderer::initialize() {
    RenderProcessor::initialize();

    shaderProp_.rebuild();
    shaderPropDVC_.rebuild();
    shaderPropMASS_.rebuild();
}

void CosmologyParticleRenderer::deinitialize() {
    RenderProcessor::deinitialize();

	if (phiValues_ != NULL) {
		phiValues_->~Volume();
	}
	if (velValues_ != NULL) {
		velValues_->~Volume();
	}
	if (uuValues_ != NULL) {
		uuValues_->~Volume();
	}
	if (muValues_ != NULL) {
		muValues_->~Volume();
	}
	//if (hhValues_ != NULL) {
	//	hhValues_->~Volume();
	//}
	//if (rhoValues_ != NULL) {
	//	rhoValues_->~Volume();
	//}
}

void CosmologyParticleRenderer::changedRenderMode(){
    RenderMode rendermode = renderMode_.getValue();
	bool tf = rendermode == VELOCITY_TRANSFER_FUNC || rendermode == PHI_TRANSFER_FUNC || rendermode == INT_EGY_TRANSFER_FUNC || rendermode == MOL_WGT_TRANSFER_FUNC; // || rendermode == SPH_LEN_TRANSFER_FUNC || rendermode == SPH_DEN_TRANSFER_FUNC;
    bool dir = rendermode == DIRECT_VELOCITY_COLOR;
	bool type = rendermode == PARTICLE_TYPE;
    bool mass = rendermode == MASS;
    transFunc_.setVisibleFlag(tf);
    massColor_.setVisibleFlag(mass);
    massTransFunc_.setVisibleFlag(mass);

    visualizesDirection_.set(dir);
    usesTf_.set(tf);
    visualizesMass_.set(mass);

    if(transFunc_.get()) {
        switch (rendermode){
        case VELOCITY_TRANSFER_FUNC:
            //transFunc_.get()->setDomain(minmaxVel_);
			unitDisplayed_.set("km/s");
            transFunc_.setVolume(velValues_);
            //transFunc_.get()->setThreshold(tgt::vec2(0, 1));
            break;
        case PHI_TRANSFER_FUNC:
            //transFunc_.get()->setDomain(minmaxPhi_);
			unitDisplayed_.set("");
            transFunc_.setVolume(phiValues_);
            //transFunc_.get()->setThreshold(tgt::vec2(0, 1));
            break;
		case INT_EGY_TRANSFER_FUNC:
			//transFunc_.get()->setDomain(minmaxUu_);
			unitDisplayed_.set("(km/s)^2");
			transFunc_.setVolume(uuValues_);
			//transFunc_.get()->setThreshold(tgt::vec2(0, 1));
			break;
		case MOL_WGT_TRANSFER_FUNC:
			//transFunc_.get()->setDomain(minmaxMu_);			
			unitDisplayed_.set("");
			transFunc_.setVolume(muValues_);
			//transFunc_.get()->setThreshold(tgt::vec2(0, 1));
			break;
		//case SPH_DEN_TRANSFER_FUNC:
		//	//transFunc_.get()->setDomain(minmaxRho_);
		//	unitDisplayed_.set("h^2*Msolar/Mpc^3");
		//	transFunc_.setVolume(rhoValues_);
		//	//transFunc_.get()->setThreshold(tgt::vec2(0, 1));
		//	break;
		//case SPH_LEN_TRANSFER_FUNC:
		//	//transFunc_.get()->setDomain(minmaxHh_);
		//	unitDisplayed_.set("Mpc/h");
		//	transFunc_.setVolume(hhValues_);
		//	//transFunc_.get()->setThreshold(tgt::vec2(0, 1));
		//	break;
		default:
			unitDisplayed_.set("");
        }
    }
    transFunc_.invalidate();
}

void CosmologyParticleRenderer::changeUseRelativeRadius(){
    bool useRelative = useRelativeRadius_.getValue();
    radiusProp_.setVisibleFlag(!useRelative);
    relativRadiusProp_.setVisibleFlag(useRelative);
}

void CosmologyParticleRenderer::changeUseAlpha(){
    alphaFactor_.setVisibleFlag(useAlpha_.get());
}

void CosmologyParticleRenderer::changedAlphaFactor(){
    float alpha = alphaFactor_.get();
    if(massTransFunc_.get()) {
        massTransFunc_.get()->setDomain(tgt::vec2(0, 1/alpha));
        massTransFunc_.get()->setThreshold(tgt::vec2(0, 1));
    }
    massTransFunc_.invalidate();
}

void CosmologyParticleRenderer::invalidateBuffers(){
    buffersInvalid_ = true;

}

namespace{
struct Ranges{
    tgt::vec2 minMaxVel;
    tgt::vec2 minMaxPhi;
	tgt::vec2 minMaxUu;
	tgt::vec2 minMaxMu;
	//tgt::vec2 minMaxHh;
	//tgt::vec2 minMaxRho;
    tgt::Bounds bounds;
};

struct FillBufferData{
    const CMParticle* particles;
    int offset;
    int count;
    float oneByTimeStep;
    tgt::vec3* posVec;
    tgt::vec3* velVec;
    float* velScalar;
    float* phiScalar;
	float* uuScalar;
	float* muScalar;
	//float* hhScalar;
	//float* rhoScalar;
	int* pType;
    Ranges* ranges;

    int pc;
};
}
namespace{
inline float FastMin(float a, float b){
    return (a < b)?a:b;
}

inline float FastMax(float a, float b){
    return (a > b)?a:b;
}

struct FastBB{
    float minx;
    float miny;
    float minz;
    float maxx;
    float maxy;
    float maxz;
    FastBB(){
        minx = 100000;   
        miny = 100000;
        minz = 100000;
        maxx = -100000;   
        maxy = -100000;
        maxz = -100000;
    }

    inline void addPoint(float x, float y, float z){
        minx = FastMin(x, minx);
        miny = FastMin(y, miny);
        minz = FastMin(z, minz);

        maxx = FastMax(x, maxx);
        maxy = FastMax(y, maxy);
        maxz = FastMax(z, maxz);

    }
};
inline float fastLength(float x, float y, float z){
#if defined(WIN32) || defined(__linux__)
    __m128 v = _mm_set_ps(x, y, z, 0.0f);
    v = _mm_mul_ps(v, v);
    v = _mm_hadd_ps(v, v);
    v = _mm_hadd_ps(v, v);
    __m128 rt = _mm_rsqrt_ps(v);
    v = _mm_mul_ps(v, rt);
    //return v.m128_f32[0];
    return *((float*)&v);
#else
    return std::sqrt(x*x+y*y+z*z);
#endif
}


}


static void fillBuffers(FillBufferData d){
    tgt::vec2 minmaxVel = tgt::vec2(1000000, -1000000);
    tgt::vec2 minmaxPhi = tgt::vec2(1000000, -1000000);
	//tgt::vec2 minmaxRho = tgt::vec2(1000000, -1000000);
	tgt::vec2 minmaxUu = tgt::vec2(1000000, -1000000);
	tgt::vec2 minmaxMu = tgt::vec2(1000000, -1000000);
	//tgt::vec2 minmaxHh = tgt::vec2(1000000, -1000000);
    float oneByTimeStep = d.oneByTimeStep;
    //tgt::Bounds bounds;
    FastBB bb;
    for(int i = d.offset; i != d.offset+d.count; i++){
        assert(i >= 0);
        assert(i < d.pc);
        CMParticle p = d.particles[i];
        //Now in normalizationmatrix
        //p.pos *= oneByTimeStep;
        float velM = fastLength(p.vel.x, p.vel.y, p.vel.z);
        d.posVec[i] = p.pos;
        d.velVec[i] = p.vel;

        d.phiScalar[i] = p.phi;
        d.velScalar[i] = velM;
		//d.rhoScalar[i] = p.rho;
		d.uuScalar[i] = p.uu;
		d.muScalar[i] = p.mu;
		//d.hhScalar[i] = p.hh;
		//d.pType[i] = p.mask % 2;

        //bounds.addPoint(p.pos);
        bb.addPoint(p.pos.x, p.pos.y, p.pos.z);

        minmaxVel = tgt::vec2(FastMin(minmaxVel.x, velM),  FastMax(minmaxVel.y, velM));
        minmaxPhi = tgt::vec2(FastMin(minmaxPhi.x, p.phi), FastMax(minmaxPhi.y, p.phi));
		//minmaxRho = tgt::vec2(FastMin(minmaxRho.x, p.rho), FastMax(minmaxRho.y, p.rho));
		minmaxMu  = tgt::vec2(FastMin(minmaxMu.x, p.mu), FastMax(minmaxMu.y, p.mu));
		minmaxUu  = tgt::vec2(FastMin(minmaxUu.x, p.uu), FastMax(minmaxUu.y, p.uu));
		//minmaxHh = tgt::vec2(FastMin(minmaxHh.x, p.hh), FastMax(minmaxHh.y, p.hh));
    }
    d.ranges->minMaxVel = minmaxVel;
    d.ranges->minMaxPhi = minmaxPhi;
	//d.ranges->minMaxRho = minmaxRho;
	d.ranges->minMaxMu = minmaxMu;
	d.ranges->minMaxUu = minmaxUu;
	//d.ranges->minMaxHh = minmaxHh;
    d.ranges->bounds    = tgt::Bounds(tgt::vec3(bb.minx, bb.miny, bb.minz), tgt::vec3(bb.maxx, bb.maxy, bb.maxz));

}


double Benchtime(){
#ifdef WIN32
    LARGE_INTEGER counter;
    LARGE_INTEGER freq;
    QueryPerformanceCounter(&counter);
    QueryPerformanceFrequency(&freq);
    return (double)counter.QuadPart/(double)freq.QuadPart;
#elif __linux__
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double f =  1.0*tv.tv_sec+tv.tv_usec/1000000.0;
    return f;

#else
    return 0.0;
#endif
}

void CosmologyParticleRenderer::setupBuffers(){
    if (vao_){ 
        glDeleteVertexArrays(1, &vao_);
        glDeleteBuffers(MAX_RENDER_MODE, vbo_);
        glDeleteBuffers(1, &posvbo_);
        vao_ = 0;
    }

    if(!inport_.getData()) {
        return;
    }

    const CMParticleData *   particleData = inport_.getData();
    CMParticleDataTimeSlice* timeSlice    = particleData->sliceAtTimeStep(timeStep_.get());

    float oneByTimeStep = 1.0f/timeSlice->getTimeStep();
    glGenVertexArrays(1, &vao_);
    glGenBuffers(5, vbo_);
    glGenBuffers(1, &posvbo_);
    //LDEBUG("generted buffer " + std::to_string(vbo_));

	std::vector<float> phiScalar, velScalar, muScalar, uuScalar; //hhScalar, rhoScalar;
    std::vector<tgt::vec3> posVec, velVec;
	//std::vector<int> pType;

    minmaxVel_ = tgt::vec2(1000000, -1000000);
    minmaxPhi_ = tgt::vec2(1000000, -1000000);
	//minmaxRho_ = tgt::vec2(1000000, -1000000);
	minmaxMu_ = tgt::vec2(1000000, -1000000);
	minmaxUu_ = tgt::vec2(1000000, -1000000);
	//minmaxHh_ = tgt::vec2(1000000, -1000000);
    bounds_ = timeSlice->getBounds();
    tgt::Bounds bounds;
    double begin = Benchtime();
    const CMParticle* particles      = timeSlice->startUsingParticles();
    int               numOfParticles = timeSlice->getNumberOfParticles();
    if (!parallelizeProp_.get() || particleData->isFiltered() || true){
        const char* enabledState = particleData->getEnabledState();
        for(int i = 0; i != numOfParticles; i++){
            CMParticle p = particles[i];
            if (!enabledState[p.ident])
                continue;
            float velM = tgt::length(p.vel);
            posVec.push_back(p.pos);
            velVec.push_back(p.vel);

            phiScalar.push_back(p.phi);
			//rhoScalar.push_back(p.rho);
            velScalar.push_back(velM);
			muScalar.push_back(p.mu);
			uuScalar.push_back(p.uu);
			//hhScalar.push_back(p.hh);
			//pType.push_back(p.mask);
            bounds.addPoint(p.pos);

            minmaxVel_ = tgt::vec2(std::min(minmaxVel_.x, velM), std::max(minmaxVel_.y, velM));
            minmaxPhi_ = tgt::vec2(std::min(minmaxPhi_.x, p.phi), std::max(minmaxPhi_.y, p.phi));
			//minmaxRho_ = tgt::vec2(std::min(minmaxRho_.x, p.rho), std::max(minmaxRho_.y, p.rho));
			minmaxMu_ = tgt::vec2(std::min(minmaxMu_.x, p.mu), std::max(minmaxMu_.y, p.mu));
			minmaxUu_ = tgt::vec2(std::min(minmaxUu_.x, p.mu), std::max(minmaxUu_.y, p.uu));
			//minmaxHh_ = tgt::vec2(std::min(minmaxHh_.x, p.hh), std::max(minmaxHh_.y, p.hh));
        }
        //Calculate the histogram for the currently selected ghj
        tgt::svec3 dataVolumeDim = tgt::svec3(velVec.size(), 1, 1);
        VolumeRAM_Float* phivram  = new VolumeRAM_Float(dataVolumeDim);
		//VolumeRAM_Float* rhovram = new VolumeRAM_Float(dataVolumeDim);
        VolumeRAM_Float* velvram  = new VolumeRAM_Float(dataVolumeDim);
		VolumeRAM_Float* muvram = new VolumeRAM_Float(dataVolumeDim);
		VolumeRAM_Float* uuvram = new VolumeRAM_Float(dataVolumeDim);
		//VolumeRAM_Float* hhvram = new VolumeRAM_Float(dataVolumeDim);
		//VolumeRAM_Float* ptpvram = new VolumeRAM_Float(dataVolumeDim);
        /*for(size_t i = 0; i != velVec.size() ; ++i) {
            float val=0;
            switch(renderMode_.getValue()) {
                phivram->setVoxelNormalized
                case VELOCITY_TRANSFER_FUNC:
                    val = velScalar[i];
                    break;
                case PHI_TRANSFER_FUNC:
                    val = phiScalar[i];
                    break;
                case ACC_TRANSFER_FUNC:
                    val = accScalar[i];
                    break;
            }
            bodyvram->setVoxelNormalized(val, i, 0, 0);
        }*/

        memcpy(static_cast<float*>(phivram->getData()), phiScalar.data(), sizeof(float)*phiScalar.size());
		//memcpy(static_cast<float*>(rhovram->getData()), rhoScalar.data(), sizeof(float)*rhoScalar.size());
        memcpy(static_cast<float*>(velvram->getData()), velScalar.data(), sizeof(float)*velScalar.size());
		memcpy(static_cast<float*>(muvram->getData()), muScalar.data(), sizeof(float)*muScalar.size());
		memcpy(static_cast<float*>(uuvram->getData()), uuScalar.data(), sizeof(float)*uuScalar.size());
		//memcpy(static_cast<float*>(hhvram->getData()), hhScalar.data(), sizeof(float)*hhScalar.size());
		//memcpy(static_cast<float*>(ptpvram->getData()), pType.data(), sizeof(float)*pType.size());

        delete(velValues_);
        velValues_ = new Volume(velvram, tgt::vec3::zero, tgt::vec3::zero);

        delete(phiValues_);
        phiValues_ = new Volume(phivram, tgt::vec3::zero, tgt::vec3::zero);

		//delete(rhoValues_);
		//rhoValues_ = new Volume(rhovram, tgt::vec3::zero, tgt::vec3::zero);

		delete(muValues_);
		muValues_ = new Volume(muvram, tgt::vec3::zero, tgt::vec3::zero);

		delete(uuValues_);
		uuValues_ = new Volume(uuvram, tgt::vec3::zero, tgt::vec3::zero);

		//delete(hhValues_);
		//hhValues_ = new Volume(hhvram, tgt::vec3::zero, tgt::vec3::zero);

		//delete(ptpValues_);
		//ptpValues_ = new Volume(ptpvram, tgt::vec3::zero, tgt::vec3::zero);
    }else{
        const int JOBS = 16;
        Ranges ranges[JOBS];
        int jobSize = numOfParticles/JOBS;
        FillBufferData d;
        posVec.resize(numOfParticles);
        velVec.resize(numOfParticles);
        velScalar.resize(numOfParticles);
        phiScalar.resize(numOfParticles);
		//rhoScalar.resize(numOfParticles);
		muScalar.resize(numOfParticles);
		uuScalar.resize(numOfParticles);
		//hhScalar.resize(numOfParticles);
		//pType.resize(numOfParticles);

        d.posVec = posVec.data();
        d.velVec = velVec.data(); 
        d.velScalar = velScalar.data(); 
        d.phiScalar = phiScalar.data(); 
		//d.rhoScalar = rhoScalar.data();
		d.muScalar = muScalar.data();
		d.uuScalar = uuScalar.data();
		//d.hhScalar = hhScalar.data();
		//d.pType = pType.data();
        d.pc = numOfParticles;
        d.oneByTimeStep = oneByTimeStep;

        d.particles = particles;
        for(int i = 0; i != JOBS-1; i++){
            d.ranges = &(ranges[i]);
            d.offset = i*jobSize;
            d.count = jobSize;
            //std::cout << d.offset << " " << d.count << std::endl;
            globalJobQueue.addFunctionJob(fillBuffers, d);
        }
        d.ranges= &(ranges[JOBS-1]);
        d.offset = numOfParticles-(JOBS-1)*jobSize;
        d.count  = numOfParticles-d.offset;
        globalJobQueue.addFunctionJob(fillBuffers, d);
        globalJobQueue.waitAll();

        for(int i = 0; i != JOBS; i++){
            Ranges r = ranges[i];

            minmaxVel_ = tgt::vec2(std::min(minmaxVel_.x, r.minMaxVel.x), std::max(minmaxVel_.y, r.minMaxVel.y));
            minmaxPhi_ = tgt::vec2(std::min(minmaxPhi_.x, r.minMaxPhi.x), std::max(minmaxPhi_.y, r.minMaxPhi.y));
			//minmaxRho_ = tgt::vec2(std::min(minmaxRho_.x, r.minMaxRho.x), std::max(minmaxRho_.y, r.minMaxRho.y));
			minmaxMu_ = tgt::vec2(std::min(minmaxMu_.x, r.minMaxMu.x), std::max(minmaxMu_.y, r.minMaxMu.y));
			minmaxUu_ = tgt::vec2(std::min(minmaxUu_.x, r.minMaxUu.x), std::max(minmaxUu_.y, r.minMaxUu.y));
			//minmaxHh_ = tgt::vec2(std::min(minmaxHh_.x, r.minMaxHh.x), std::max(minmaxHh_.y, r.minMaxHh.y));
            bounds.addPoint(r.bounds.getLLF());
            bounds.addPoint(r.bounds.getURB());
        }
    }
    timeSlice->finishedUsingParticles();
    double time = Benchtime()-begin;
    LDEBUG("Buffer creation took " + std::to_string(time) + "s with" +  (parallelizeProp_.get()?"":"out") + " parallel building");

    glBindBuffer(GL_ARRAY_BUFFER, posvbo_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(tgt::vec3)*posVec.size(), posVec.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, vbo_[VELOCITY_TRANSFER_FUNC]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*velScalar.size(), velScalar.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, vbo_[DIRECT_VELOCITY_COLOR]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(tgt::vec3)*velVec.size(), velVec.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, vbo_[PHI_TRANSFER_FUNC]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*phiScalar.size(), phiScalar.data(), GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, vbo_[INT_EGY_TRANSFER_FUNC]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*uuScalar.size(), uuScalar.data(), GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, vbo_[MOL_WGT_TRANSFER_FUNC]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*muScalar.size(), muScalar.data(), GL_STATIC_DRAW);

	//glBindBuffer(GL_ARRAY_BUFFER, vbo_[SPH_DEN_TRANSFER_FUNC]);
	//glBufferData(GL_ARRAY_BUFFER, sizeof(float)*rhoScalar.size(), rhoScalar.data(), GL_STATIC_DRAW);

	//glBindBuffer(GL_ARRAY_BUFFER, vbo_[SPH_LEN_TRANSFER_FUNC]);
	//glBufferData(GL_ARRAY_BUFFER, sizeof(float)*hhScalar.size(), hhScalar.data(), GL_STATIC_DRAW);

	//glBindBuffer(GL_ARRAY_BUFFER, vbo_[PARTICLE_TYPE]);
	//glBufferData(GL_ARRAY_BUFFER, sizeof(float)*pType.size(), pType.data(), GL_STATIC_DRAW);

    particleCount_ = (int)posVec.size();
    buffersInvalid_ = false;
    if (adaptCamera_.get())
        camera_.adaptInteractionToScene(timeSlice->getBounds());

    bounds = CMtransformBounds(timeSlice->getNormalisationTransformation(), bounds);
    universeDimensions_.setMinValue(bounds.getLLF() - tgt::vec3::one);
    universeDimensions_.setMaxValue(bounds.getURB() + tgt::vec3::one);
    universeDimensions_.set(bounds);

    universeMatrix_.set(tgt::mat4::identity);
    changedRenderMode();
}




//#pragma GCC pop_options
void CosmologyParticleRenderer::process() {    
    if (buffersInvalid_)
        setupBuffers();

    if (!vao_ || ! vbo_) return;


    RenderMode rendermode = renderMode_.getValue();
    bool       useAlpha   = useAlpha_.get();
    bool       scalarMode = isScalarRenderMode(rendermode);
    float      radius;
    if (useRelativeRadius_.getValue()){
        radius = relativRadiusProp_.get()*tgt::length(bounds_.diagonal())*0.02f;
    }else{
        radius = radiusProp_.get();
    }

    
    outport_.activateTarget();
    outport_.clearTarget();

    tgt::Shader* shader;
	if (rendermode == MASS) {
		shader = shaderPropMASS_.getShader();
	//}else if (rendermode == PARTICLE_TYPE) {
	//	shader = shaderPropPTP_.getShader();
    }else if (scalarMode){
        shader = shaderProp_.getShader();
    }else{
        shader = shaderPropDVC_.getShader();
    }
    shader->activate();

    tgt::Camera cam            = camera_.get();
    tgt::mat4 projectionMatrix = cam.getProjectionMatrix(outport_.getSize());
    if(!inport_.getData()) {
        return;
    }
    tgt::mat4 normMatrix       = inport_.getData()->sliceAtTimeStep(timeStep_.get())->getNormalisationTransformation();
    tgt::mat4 viewMatrix       = cam.getViewMatrix();
    tgt::mat4 normAndViewMatrix = viewMatrix*normMatrix;
    
    shader->setUniform("projectionMatrix", projectionMatrix);
    shader->setUniform("viewMatrix", normAndViewMatrix);
    shader->setUniform("radius_", radius);
    shader->setUniform("alphaFactor_", alphaFactor_.get());
	GLenum gl_error;
    glBindVertexArray(vao_);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, posvbo_);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	gl_error = glGetError();
	if (gl_error) {
		std::cout << "GL Error: " << gl_error << "\n";
	}
    if (rendermode == MASS){
        shader->setUniform("color_", massColor_.get());
    } else if (scalarMode){
        glBindBuffer(GL_ARRAY_BUFFER, vbo_[rendermode]);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, 0);

        //TODO: code for domain!
        tgt::vec2 domain = transFunc_.get()->getDomain();
        //std::cerr << domain << std::endl;
        //shader->setUniform("scale_", domain.y-domain.x);
        //shader->setUniform("offset_", domain.x);

        tgt::TextureUnit transferUnit;
        if (transFunc_.get()) {
            transferUnit.activate();
            transFunc_.get()->getTexture()->bind();
            transFunc_.get()->setUniform(shader, "transferFunc_", "transferFuncTex_", transferUnit.getUnitNumber());  
        }
    }else{
        glBindBuffer(GL_ARRAY_BUFFER, vbo_[rendermode]);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
    }

    if (useAlpha){
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE);
        glDisable(GL_DEPTH_TEST);
    } 
    glDrawArrays(GL_POINTS, 0, particleCount_);
    glEnable(GL_DEPTH_TEST);
    glBindVertexArray(0);
	gl_error = glGetError();
	if (gl_error) {
		std::cout << "GL Error: " << gl_error << "\n";
	}
    glDisable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ZERO);
	gl_error = glGetError();
	if (gl_error) {
		std::cout << "GL Error: " << gl_error << "\n";
	}
	assert(glGetError() == GL_NO_ERROR);

    // cleanup
    shader->deactivate();
    outport_.deactivateTarget();
    tgt::TextureUnit::setZeroUnit();
    LGL_ERROR;
}

} // namespace