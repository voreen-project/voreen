/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2015 University of Muenster, Germany.                        *
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
#define M_PI 3.14f
//Bitflags for particle type
#define FLAG_BARYON     2  // 2^0, bit 0   Baryon 1, Dark matter 0
#define FLAG_STAR      32  // 2^6, bit 6   star   1
#define FLAG_WIND      64  // 2^7, bit 7   wind   1
#define FLAG_GAS      128  // 2^8, bit 8   gas    1
#define FLAG_AGN      256  // 2^9, bit 9   AGN    1

#include <iostream>
#include "cosmologyvolumeconverter.h"
#include "../utils/cmmath.h"
#include <math.h>
#include <algorithm>
#include <chrono>

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumegl.h"


namespace voreen{

/******************************************************************************
 *
 *                            Some Volume Builders
 *
 *****************************************************************************/
namespace{
struct NearestVolumeBuilder{
    tgt::Bounds      bounds;
    tgt::svec3       dim;
    VolumeRAM_Float* volumeData;
    MetaDataContainer metadata;

    inline void addPoint(tgt::vec3 pos, float intensity, float smoothingLength, voreen::CosmologyVolumeConverter::ParticleType, uint16_t mask, voreen::CosmologyVolumeConverter::ParticleProperty particleProperty, float time){
        tgt::vec3 pos_vx = (pos - bounds.getLLF()) / (bounds.getURB() - bounds.getLLF())*(tgt::vec3(dim)-tgt::vec3::one);
		assert(pos_vx.x >= 0 && pos_vx.y >= 0 && pos_vx.z >= 0);
		assert(pos_vx.x < dim.x && pos_vx.y < dim.y && pos_vx.z < dim.z);

		tgt::ivec3 nearest_vx = tgt::ivec3(tgt::floor(pos_vx+tgt::vec3(0.5f)));
        assert(nearest_vx.x < dim.x && nearest_vx.y < dim.y && nearest_vx.z < dim.z);

		volumeData->voxel(nearest_vx) += intensity;
    }
};

struct BackwardLinearInterpolationBuilder{
    tgt::Bounds      bounds;
    tgt::svec3       dim;
    VolumeRAM_Float* volumeData;
    MetaDataContainer metadata;

    inline void addPoint(tgt::vec3 pos, float intensity, float smoothingLength, voreen::CosmologyVolumeConverter::ParticleType, uint16_t mask, voreen::CosmologyVolumeConverter::ParticleProperty particleProperty, float time){
        tgt::vec3 pos_vx = (pos - bounds.getLLF()) / (bounds.getURB() - bounds.getLLF())*(tgt::vec3(dim)-tgt::vec3::one);
		assert(pos_vx.x >= 0 && pos_vx.y >= 0 && pos_vx.z >= 0);
		assert(pos_vx.x < dim.x && pos_vx.y < dim.y && pos_vx.z < dim.z);

		tgt::ivec3 llf_vx = tgt::ivec3(tgt::floor(pos_vx));

		tgt::vec3 w_llf_vx = pos_vx - tgt::vec3(llf_vx);
		tgt::vec3 w_urb_vx = tgt::vec3::one - w_llf_vx;

		volumeData->voxel(llf_vx + tgt::ivec3(0, 0, 0)) += intensity * w_urb_vx.x * w_urb_vx.y * w_urb_vx.z;
		volumeData->voxel(llf_vx + tgt::ivec3(0, 0, 1)) += intensity * w_urb_vx.x * w_urb_vx.y * w_llf_vx.z;
		volumeData->voxel(llf_vx + tgt::ivec3(0, 1, 0)) += intensity * w_urb_vx.x * w_llf_vx.y * w_urb_vx.z;
		volumeData->voxel(llf_vx + tgt::ivec3(0, 1, 1)) += intensity * w_urb_vx.x * w_llf_vx.y * w_llf_vx.z;

		volumeData->voxel(llf_vx + tgt::ivec3(1, 0, 0)) += intensity * w_llf_vx.x * w_urb_vx.y * w_urb_vx.z;
		volumeData->voxel(llf_vx + tgt::ivec3(1, 0, 1)) += intensity * w_llf_vx.x * w_urb_vx.y * w_llf_vx.z;
		volumeData->voxel(llf_vx + tgt::ivec3(1, 1, 0)) += intensity * w_llf_vx.x * w_llf_vx.y * w_urb_vx.z;
		volumeData->voxel(llf_vx + tgt::ivec3(1, 1, 1)) += intensity * w_llf_vx.x * w_llf_vx.y * w_llf_vx.z;
    }
};

struct WeightAdditionBuilder{
    tgt::Bounds      bounds;
    tgt::svec3       dim;
    VolumeRAM_Float* volumeData;
    MetaDataContainer metadata;

    inline void addPoint(tgt::vec3 pos, float intensity, float smoothingLength, voreen::CosmologyVolumeConverter::ParticleType, uint16_t mask, voreen::CosmologyVolumeConverter::ParticleProperty particleProperty, float time){
        tgt::vec3 pos_vx = (pos - bounds.getLLF()) / (bounds.getURB() - bounds.getLLF())*(tgt::vec3(dim)-tgt::vec3::one);
		assert(pos_vx.x > 0 && pos_vx.y > 0 && pos_vx.z > 0);
		assert(pos_vx.x < dim.x && pos_vx.y < dim.y && pos_vx.z < dim.z);

		tgt::ivec3 llf_vx = tgt::ivec3(tgt::floor(pos_vx));

		tgt::vec3 w_llf_vx = pos_vx - tgt::vec3(llf_vx);
		tgt::vec3 w_urb_vx = tgt::vec3::one - w_llf_vx;

        assert(w_llf_vx.x >= 0 && w_llf_vx.y >= 0 && w_llf_vx.z >= 0);

        // 1/12 makes sure we conserve the total mass of the particles
        // because we have 12 pairs w_urb_vx.d+w_llf_vx.d == 1
		volumeData->voxel(llf_vx + tgt::ivec3(0, 0, 0)) += intensity*(1.0f/12.0f)*(w_urb_vx.x + w_urb_vx.y + w_urb_vx.z);
		volumeData->voxel(llf_vx + tgt::ivec3(0, 0, 1)) += intensity*(1.0f/12.0f)*(w_urb_vx.x + w_urb_vx.y + w_llf_vx.z);
		volumeData->voxel(llf_vx + tgt::ivec3(0, 1, 0)) += intensity*(1.0f/12.0f)*(w_urb_vx.x + w_llf_vx.y + w_urb_vx.z);
		volumeData->voxel(llf_vx + tgt::ivec3(0, 1, 1)) += intensity*(1.0f/12.0f)*(w_urb_vx.x + w_llf_vx.y + w_llf_vx.z);

		volumeData->voxel(llf_vx + tgt::ivec3(1, 0, 0)) += intensity*(1.0f/12.0f)*(w_llf_vx.x + w_urb_vx.y + w_urb_vx.z);
		volumeData->voxel(llf_vx + tgt::ivec3(1, 0, 1)) += intensity*(1.0f/12.0f)*(w_llf_vx.x + w_urb_vx.y + w_llf_vx.z);
		volumeData->voxel(llf_vx + tgt::ivec3(1, 1, 0)) += intensity*(1.0f/12.0f)*(w_llf_vx.x + w_llf_vx.y + w_urb_vx.z);
		volumeData->voxel(llf_vx + tgt::ivec3(1, 1, 1)) += intensity*(1.0f/12.0f)*(w_llf_vx.x + w_llf_vx.y + w_llf_vx.z);
    }
};

struct SPHBuilder{
    tgt::Bounds      bounds;
    tgt::svec3       dim;
    VolumeRAM_Float* volumeData;
    MetaDataContainer metadata;




	inline float square(float value) {
		return value * value;
	}

	inline float cube(float value) {
		return value * value * value;
	}

	inline float powSix(float value) {
		return value * value * value * value * value * value;
	}

	float C = 495.0f / (32.0f * M_PI);

	inline float calculateKernelWeight(float distance, float smoothingLength, float constant) {
		float relativeDistance = distance / smoothingLength;

		return constant * powSix(1.0f - relativeDistance) * ((1.0f + 6.0f * relativeDistance + (35.0f / 3.0f) * square(relativeDistance)));

	}

    /*inline float calculateDistance(tgt::vec3 pos1, tgt::vec3 pos2){
        float distance = std::sqrt(std::pow((pos2.x - pos1.x), 2.0) + std::pow((pos2.y - pos1.y), 2.0) + std::pow((pos2.y - pos1.y), 2.0));

        return distance;
    }*/

    /*inline float calculateSPHValue(tgt::vec3 pos1, tgt::vec3 pos2, float smoothingLength, float intensity) {
        float sphValue = 0.0f;
        //float distance = calculateDistance(pos1, pos2);
		float distance = tgt::distance(pos1, pos2);
        float kernelWeigth = calculateKernelWeight(distance, smoothingLength);
        sphValue = intensity * kernelWeigth;
        return sphValue;
    }*/


    inline void addPoint(tgt::vec3 pos, float intensity, float smoothingLength, voreen::CosmologyVolumeConverter::ParticleType, uint16_t mask, voreen::CosmologyVolumeConverter::ParticleProperty particleProperty, float time){

		if (smoothingLength <= 0.0f) {
			return;
		}

        float timeStepScaleFactor = (float)(((1.0l - (1.0l / 201.0l)) / 625.0l) * (time + 1.0l)) + (1.0l / 201.0l);
        float temperatureConstant = 4.8f * std::pow(10, 5) * std::pow(timeStepScaleFactor, 3);

		tgt::vec3 pos_particle = (pos + 32.0f);

		float constant = C / cube(smoothingLength);
		float scalingFactor = dim.x / 64.0f;

        /*assert(pos_vx.x >= 0 && pos_vx.y >= 0 && pos_vx.z >= 0);
        assert(pos_vx.x < dim.x && pos_vx.y < dim.y && pos_vx.z < dim.z);*/

        tgt::ivec3 llf_vx = tgt::ivec3(tgt::ceil((pos_particle - smoothingLength) * scalingFactor));
        tgt::ivec3 urb_vx = tgt::ivec3(tgt::floor((pos_particle + smoothingLength) * scalingFactor));

        tgt::Bounds checkbounds = tgt::Bounds(tgt::vec3(0.0f,0.0f,0.0f), (dim - tgt::svec3(1,1,1)));

#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for shared(llf_vx, urb_vx, checkbounds, constant, smoothingLength, temperatureConstant)
#endif
		for (int k = llf_vx.z; k < urb_vx.z; ++k) {
			for (int j = llf_vx.y; j < urb_vx.y; ++j) {
#ifdef VRN_MODULE_OPENMP
#define omp simd
#endif
				for (int i = llf_vx.x; i < urb_vx.x; ++i) {
					tgt::vec3 voxelpos = tgt::ivec3(i, j, k);
					voxelpos = (voxelpos + tgt::vec3(0.5f, 0.5f, 0.5f)) / scalingFactor;
					if (i >= 0 && j >= 0 && k >= 0 && i < dim.x && j < dim.y && k < dim.z) {

						float distance = tgt::distance(pos_particle, voxelpos);

						if (distance > smoothingLength) {
							continue;
						}

						int index = i + dim.x * j + dim.x * dim.y * k;

						float kernelWeigth = calculateKernelWeight(distance, smoothingLength, constant);
						float value = kernelWeigth * intensity;

						volumeData->voxel(index) += (particleProperty == 9) ? value * temperatureConstant : value;

					}
				}
			}
		}
    }
};

struct AmountVolumeBuilder{
    tgt::Bounds      bounds;
    tgt::svec3       dim;
    VolumeRAM_Float* volumeData;
    MetaDataContainer metadata;

    uint16_t allBaryons_ = 2;
    uint16_t starBaryons_ = 34;
    uint16_t windBaryons_ = 66;
    uint16_t gasBaryons_ = 130;
    uint16_t agn_ = 256;


    inline void addPoint(tgt::vec3 pos, float intensity, float smoothingLength, voreen::CosmologyVolumeConverter::ParticleType particleType, uint16_t mask, voreen::CosmologyVolumeConverter::ParticleProperty particleProperty, float time){

        int b_ = false;

        switch (particleType) {
            case 00:
                b_ = (mask & allBaryons_) == allBaryons_;
                break;
            case 01:
                b_ = (mask & allBaryons_) == allBaryons_;
                b_ = !b_;
                break;
            case 2:
                b_ = (mask & windBaryons_) == windBaryons_;
                break;
            case 3:
                b_ = (mask & starBaryons_) == starBaryons_;
                break;
            case 4:
                b_ = (mask & gasBaryons_) == gasBaryons_;
                break;
            case 5:
                b_ = (mask & agn_) == agn_;
            case 6:
                b_ = true;
                break;
        }

        if(b_){
            tgt::vec3 pos_vx = (pos - bounds.getLLF()) / (bounds.getURB() - bounds.getLLF())*(tgt::vec3(dim)-tgt::vec3::one);
            assert(pos_vx.x >= 0 && pos_vx.y >= 0 && pos_vx.z >= 0);
            assert(pos_vx.x < dim.x && pos_vx.y < dim.y && pos_vx.z < dim.z);

            tgt::ivec3 nearest_vx = tgt::ivec3(tgt::floor(pos_vx+tgt::vec3(0.5f)));
            assert(nearest_vx.x < dim.x && nearest_vx.y < dim.y && nearest_vx.z < dim.z);


            volumeData->voxel(nearest_vx) += 1;
        }
    }
};


}

/******************************************************************************
 *
 *                              CosmologyVolumeConverter
 *
 *****************************************************************************/
CosmologyVolumeConverter::CosmologyVolumeConverter()
    : CachingVolumeProcessor()
    , outport_         (Port::OUTPORT, "volumehandle.output",   "Volume Output", false)
    , inport_          (Port::INPORT,  "particlehandle.output", "Particle Data Output")
    , volumeColor_     ("volumeColor",      "Volume Color", 0.5f)
    , timeStep_        ("timeStep",         "Time Step", 0.0f, 0.0f, 624.0f)
	, volumeDimensions_("volumeDimensions", "Volume Dimensions", tgt::ivec3(64), tgt::ivec3::zero, tgt::ivec3(1024))
    , spreadMode_      ("spreadMode",       "Mode for spreading of Mass")
    , particleProperty_("particleProperty", "Particle Property")
    , particleType_    ("particleType", "Particle Type")
    , propertyUnit_    ("propertyUnit",     "Unit of data displayed", "")
    , volume_(nullptr)
{
    addPort(outport_);
    addPort(inport_);

    addProperty(volumeColor_);
    addProperty(volumeDimensions_);
    addProperty(spreadMode_);
    addProperty(timeStep_);
    addProperty(propertyUnit_);

    addProperty(particleProperty_);
    addProperty(particleType_);
    
    spreadMode_.addOption("backwardLinearInterpolation", "Backward linear interpolation", BACKWARD_LINEAR_INTERPOLATION);
    spreadMode_.addOption("weightAddition"             , "Addition of edge-weights"     , WEIGHT_ADDITION);
    spreadMode_.addOption("nearestVoxel"               , "Nearest voxel "               , NEAREST_VOXEL);
    spreadMode_.addOption("SPH"                        , "SPH"                          , SPH);
    spreadMode_.addOption("AMOUNT"                     , "Amount"                       , AMOUNT);
    spreadMode_.addOption("VELOCITYFIELD"              , "Velocitifield"                , VELOCITYFIELD);

    particleProperty_.addOption("mass",          "Mass",                        MASS);
    particleProperty_.addOption("phi",           "Potential Energy",            PHI);
    particleProperty_.addOption("velocity",      "Magnitude of Velocity",       VELOCITY);
    particleProperty_.addOption("velocityx",     "Velocity in X Direction",     VELOCITY_X);
    particleProperty_.addOption("velocityy",     "Velocity in Y Direction",     VELOCITY_Y);
    particleProperty_.addOption("velocityz",     "Velocity in Z Direction",     VELOCITY_Z);
    particleProperty_.addOption("uu",            "internal energy",             UU);
    particleProperty_.addOption("temperature",   "temperature",                 TEMPERATURE);
    particleProperty_.addOption("sphDensity",    "sphDensity",                  SPHDensity);

    particleType_.addOption("baryon",         "Baryon",         BARYON);
    particleType_.addOption("dark matter",    "Dark matter",    DARK_MATTER);
    particleType_.addOption("wind",           "Wind",           WIND);
    particleType_.addOption("star",           "Star",           STAR);
    particleType_.addOption("gas",            "Gas",            GAS);
    particleType_.addOption("agn",            "AGN",            AGN);
    particleType_.addOption("all",            "ALL",            ALL);

    particleType_.setVisibleFlag(false);
    particleProperty_.setVisibleFlag(false);
    propertyUnit_.setVisibleFlag(false);

    ON_CHANGE(volumeDimensions_, CosmologyVolumeConverter, changedVolumeDimensions);
    ON_CHANGE(particleProperty_, CosmologyVolumeConverter, changedParticleProperty);
    ON_CHANGE(spreadMode_, CosmologyVolumeConverter, changedSpreadMode);
}

CosmologyVolumeConverter::~CosmologyVolumeConverter(){

}

Processor* CosmologyVolumeConverter::create() const{
    return new CosmologyVolumeConverter;
}

void CosmologyVolumeConverter::initialize() {
	changedVolumeDimensions();
}

void CosmologyVolumeConverter::changedVolumeDimensions(){
	//delete volume_;	
    volume_ = 0;
}
void CosmologyVolumeConverter::changedParticleProperty(){
    switch(particleProperty_.getValue()) {
        case MASS:
            propertyUnit_.set("Particles");
            break;
        case PHI:
            propertyUnit_.set("phi");
            break;
        case VELOCITY:
        case VELOCITY_X:
        case VELOCITY_Y:
        case VELOCITY_Z:
            propertyUnit_.set("km/s");
            break;
        case UU:
            propertyUnit_.set("(km/s)^2");
            break;
        default:
            propertyUnit_.set("");
    }
}

void CosmologyVolumeConverter::changedSpreadMode(){
    switch(spreadMode_.getValue()) {
        case BACKWARD_LINEAR_INTERPOLATION:
        case WEIGHT_ADDITION:
        case NEAREST_VOXEL:
        case SPH:
           // removeProperty(particleType_);
            propertyUnit_.setVisibleFlag(true);
            particleType_.setVisibleFlag(false);
            particleProperty_.setVisibleFlag(true);
            break;
        case AMOUNT:
            propertyUnit_.setVisibleFlag(false);
            particleProperty_.setVisibleFlag(false);
            particleType_.setVisibleFlag(true);
            break;
        case VELOCITYFIELD:
            propertyUnit_.setVisibleFlag(false);
            particleProperty_.setVisibleFlag(false);
            particleType_.setVisibleFlag(false);
            particleProperty_.set("Velocities");
        default:
            ;
    }
}

void CosmologyVolumeConverter::deinitialize() {
    //delete volume_;
    volume_ = 0;
}

void CosmologyVolumeConverter::process(){
    if (!inport_.getData()) return;

    float            color            = volumeColor_.get();
    SpreadMode       spreadMode       = spreadMode_.getValue();
    ParticleProperty particleProperty = particleProperty_.getValue();
    ParticleType     particleType     = particleType_.getValue();


    CMParticleDataTimeSlice* timeSlice     = inport_.getData()->sliceAtTimeStep(timeStep_.get());
    tgt::mat4                matrix        = timeSlice->getNormalisationTransformation();
    tgt::Bounds              bounds        = timeSlice->getBounds();
    int                      particleCount = timeSlice->getNumberOfParticles();
    tgt::Bounds              dataBounds    = CMtransformBounds(CMinvert(matrix), bounds);

    // Add some margin so we don't go out of range with particles exactly on the border
    dataBounds.addPoint(dataBounds.getLLF() - tgt::vec3(0.5f));
    dataBounds.addPoint(dataBounds.getURB() + tgt::vec3(0.5f));

    tgt::svec3 dim    = tgt::svec3(volumeDimensions_.get());
    tgt::vec3 spacing = bounds.diagonal()/tgt::vec3(dim);
    //tgt::vec3 offset  = -bounds.center();
    tgt::vec3 offset  = tgt::vec3::zero;
	
	VolumeRAM_Float* volumeData = new VolumeRAM_Float(dim, true);

    const CMParticle* particles = timeSlice->startUsingParticles();

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    // Zero volume
    std::fill_n((float*)volumeData->getData(), tgt::hmul(dim), 0.0f);
    const char* enabledState = inport_.getData()->getEnabledState();
#define buildVolumeExplicit(builder_type, intensity_expr, smoothingLength_expr, partyleType_expr, partyleMask_expr, partyleProperty_expr, time_expr)\
{\
    builder_type builder;\
    builder.dim = dim;\
    builder.bounds = dataBounds;\
    builder.volumeData = volumeData;\
    for (int i = 0; i != particleCount; i++){\
        CMParticle p = particles[i];\
        if (!enabledState[p.ident])\
            continue;\
        tgt::vec3 pos = p.pos;\
        builder.addPoint(pos, (intensity_expr), (smoothingLength_expr), (partyleType_expr), (partyleMask_expr), (partyleProperty_expr), (time_expr));\
    }\
}

#define buildVolume(expr, expr2, expr3, expr4, expr5, expr6)\
{\
    if (spreadMode == NEAREST_VOXEL){\
        buildVolumeExplicit(NearestVolumeBuilder, (expr), (expr2), (expr3), (expr4), (expr5), (expr6));\
    }else if (spreadMode == BACKWARD_LINEAR_INTERPOLATION){\
        buildVolumeExplicit(BackwardLinearInterpolationBuilder, (expr), (expr2), (expr3), (expr4), (expr5), (expr6));\
    }else if (spreadMode == WEIGHT_ADDITION){\
        buildVolumeExplicit(WeightAdditionBuilder, (expr), (expr2), (expr3), (expr4), (expr5), (expr6));\
    }else if (spreadMode == SPH){\
        buildVolumeExplicit(SPHBuilder, (expr), (expr2), (expr3), (expr4), (expr5), (expr6));\
    }else if (spreadMode == AMOUNT){\
        buildVolumeExplicit(AmountVolumeBuilder, (expr), (expr2), (expr3), (expr4), (expr5), (expr6));\
    }\
}


    switch (particleProperty) {
        case MASS: buildVolume(p.mass, 0.0f, particleType, p.mask, particleProperty, timeStep_.get());
            break;
        case PHI: buildVolume(p.phi, p.hh, particleType, p.mask, particleProperty, timeStep_.get());
            break;
        case VELOCITY: buildVolume(tgt::length(p.vel), 0.0f, particleType, p.mask, particleProperty, timeStep_.get());
            break;
        case VELOCITY_X: buildVolume(p.vel.x, 0.0f, particleType, p.mask, particleProperty, timeStep_.get());
            break;
        case VELOCITY_Y: buildVolume(p.vel.y, 0.0f, particleType, p.mask, particleProperty, timeStep_.get());
            break;
        case VELOCITY_Z: buildVolume(p.vel.z, 0.0f, particleType, p.mask, particleProperty, timeStep_.get());
            break;
        case UU: buildVolume(p.uu, p.hh, particleType, p.mask, particleProperty, timeStep_.get());
            break;
        case TEMPERATURE: buildVolume(p.uu, p.hh, particleType, p.mask, particleProperty, timeStep_.get());
            break;
        case SPHDensity: buildVolume(p.rho, p.hh, particleType, p.mask, particleProperty, timeStep_.get());
            break;
        case Velocities: {
            tgt::Bounds bounds;
            tgt::svec3 dim;
            MetaDataContainer metadata;
            VolumeRAM_3xFloat *volumeData2 = new VolumeRAM_3xFloat(dim, true);
            for (int i = 0; i != particleCount; i++) {
                \
        CMParticle p = particles[i];\
        if (!enabledState[p.ident])\

                    continue;
                        \

                tgt::vec3 pos = p.pos;\


                tgt::vec3 pos_vx = (pos - bounds.getLLF()) / (bounds.getURB() - bounds.getLLF()) *
                                   (tgt::vec3(dim) - tgt::vec3::one);
                assert(pos_vx.x >= 0 && pos_vx.y >= 0 && pos_vx.z >= 0);
                assert(pos_vx.x < dim.x && pos_vx.y < dim.y && pos_vx.z < dim.z);
                tgt::ivec3 nearest_vx = tgt::ivec3(tgt::floor(pos_vx + tgt::vec3(0.5f)));
                assert(nearest_vx.x < dim.x && nearest_vx.y < dim.y && nearest_vx.z < dim.z);

                volumeData2->voxel(nearest_vx) += p.vel;


            }
            volume_ = new Volume(volumeData2, spacing, offset);
            volume_->setTimestep(timeStep_.get());
            volume_->getMetaDataContainer().addMetaData("Scalar", new StringMetaData("uu"));
            volume_->getMetaDataContainer().addMetaData("simulated_time", new FloatMetaData(timeStep_.get()));
            break;
        }
    }

    timeSlice->finishedUsingParticles();

	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	std::cout << "Execution time: " << (duration / 1000.0f);


    if (particleProperty != Velocities){
        volume_ = new Volume(volumeData, spacing, offset);
        volume_->setTimestep(timeStep_.get());
        volume_->getMetaDataContainer().addMetaData("Scalar", new StringMetaData("uu"));
        volume_->getMetaDataContainer().addMetaData("simulated_time", new FloatMetaData(timeStep_.get()));
    }
    outport_.setData(volume_);
}
}