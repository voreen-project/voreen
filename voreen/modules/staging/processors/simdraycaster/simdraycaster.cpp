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

#include "simdraycaster.h"

#include "../../utils/simdraycaster/jobqueue.h"
#include "../../utils/simdraycaster/memory.h"
#include "../../utils/simdraycaster/brickedvolume.h"

#include "tgt/cpucapabilities.h"

#include "voreen/core/utils/backgroundthread.h"
#include "voreen/core/datastructures/transfunc/1d/1dkeys/transfunc1dkeys.h"
#include "voreen/core/datastructures/transfunc/1d/preintegrationtable.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/voreenapplication.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

using namespace voreen;
using namespace tgt;


namespace voreen {
// Parameter for SIMDRaycaster raycasting function
namespace simdraycaster{

struct RayCasterParameters{
    int offset;
    int linelength;
    int count;
    int stride;
    vec4 *entryBuffer;
    vec4 *exitBuffer;
    int bytesPerVoxel;
    vec4* tf;
    int tfsize;
    float tfscale;
    float tfoffset;
    uint32_t* img;
    ivec3 volDim;
    const void* vol;
    float stepsize;
    const size_t* offsetx;
    const size_t* offsety;
    const size_t* offsetz;
};

#define RAYCASTER_FUNCTION(name) void name(RayCasterParameters p)

    /**
     * A function that does the raycasting.
     */
    typedef void (*RayCasterFunction)(RayCasterParameters p);

// Declaration of all raycasting functions
#ifdef SIMD_SSE41
    RAYCASTER_FUNCTION(raycaster_uint8_t_sse41);
    RAYCASTER_FUNCTION(raycaster_uint16_t_sse41);
    RAYCASTER_FUNCTION(raycaster_uint8_t_bricked_sse41);
    RAYCASTER_FUNCTION(raycaster_uint16_t_bricked_sse41);
    RAYCASTER_FUNCTION(raycaster_uint8_t_preintegration_sse41);
    RAYCASTER_FUNCTION(raycaster_uint16_t_preintegration_sse41);
    RAYCASTER_FUNCTION(raycaster_uint8_t_bricked_preintegration_sse41);
    RAYCASTER_FUNCTION(raycaster_uint16_t_bricked_preintegration_sse41);
    RAYCASTER_FUNCTION(raycaster_uint8_t_mip_sse41);
    RAYCASTER_FUNCTION(raycaster_uint16_t_mip_sse41);
    RAYCASTER_FUNCTION(raycaster_uint8_t_bricked_mip_sse41);
    RAYCASTER_FUNCTION(raycaster_uint16_t_bricked_mip_sse41);
#endif 
#ifdef SIMD_SSE3
    RAYCASTER_FUNCTION(raycaster_uint8_t_sse3);
    RAYCASTER_FUNCTION(raycaster_uint16_t_sse3);
    RAYCASTER_FUNCTION(raycaster_uint8_t_bricked_sse3);
    RAYCASTER_FUNCTION(raycaster_uint16_t_bricked_sse3);
    RAYCASTER_FUNCTION(raycaster_uint8_t_preintegration_sse3);
    RAYCASTER_FUNCTION(raycaster_uint16_t_preintegration_sse3);
    RAYCASTER_FUNCTION(raycaster_uint8_t_bricked_preintegration_sse3);
    RAYCASTER_FUNCTION(raycaster_uint16_t_bricked_preintegration_sse3);
    RAYCASTER_FUNCTION(raycaster_uint8_t_mip_sse3);
    RAYCASTER_FUNCTION(raycaster_uint16_t_mip_sse3);
    RAYCASTER_FUNCTION(raycaster_uint8_t_bricked_mip_sse3);
    RAYCASTER_FUNCTION(raycaster_uint16_t_bricked_mip_sse3);
#endif

    enum SSELevel{
        SSE41,
        SSE3
    };

    /**
     * A Raycasting function with it's properties
     */
    struct RayCasterConf{
        /**
         * How many byte every voxel hast.
         *
         * 1 => uint8_t
         * 2 => uint16_t
         */
        int bytes_per_voxel;
        /**
         * Does the raycaster use bricked �volumes
         * or usual  (3D-Array) volumes
         */
        bool bricked;
        /**
         * Does the raycaster use preintegration,
         * instead of 1D Transfer functions
         */
        bool preintegration;
        /**
         * Does the raycaster do MIP rendering instead of DVR
         */
        bool mip;

        SSELevel sselevel;

        /**
         * The function, that does the specified raycasting
         */
        RayCasterFunction function;
    };

    /**
     * List of all Raycasters with their properties
     */
    RayCasterConf raycasters[] = {
    // [bytes_per_voxel]
    //      [bricked]
    //             [preintegraton]
    //                    [mip]
    //                           [sselevel]
#ifdef SIMD_SSE41
        {1, false, false, false, SSE41, raycaster_uint8_t_sse41},
        {2, false, false, false, SSE41, raycaster_uint16_t_sse41},
        {1, true,  false, false, SSE41, raycaster_uint8_t_bricked_sse41},
        {2, true,  false, false, SSE41, raycaster_uint16_t_bricked_sse41},
        {1, false, true,  false, SSE41, raycaster_uint8_t_preintegration_sse41},
        {2, false, true,  false, SSE41, raycaster_uint16_t_preintegration_sse41},
        {1, true,  true,  false, SSE41, raycaster_uint8_t_bricked_preintegration_sse41},
        {2, true,  true,  false, SSE41, raycaster_uint16_t_bricked_preintegration_sse41},
        {1, false, false, true,  SSE41, raycaster_uint8_t_mip_sse41},
        {2, false, false, true,  SSE41, raycaster_uint16_t_mip_sse41},
        {1, true,  false, true,  SSE41, raycaster_uint8_t_bricked_mip_sse41},
        {2, true,  false, true,  SSE41, raycaster_uint16_t_bricked_mip_sse41},
#endif
#ifdef SIMD_SSE3
        {1, false, false, false, SSE3, raycaster_uint8_t_sse3},
        {2, false, false, false, SSE3, raycaster_uint16_t_sse3},
        {1, true,  false, false, SSE3, raycaster_uint8_t_bricked_sse3},
        {2, true,  false, false, SSE3, raycaster_uint16_t_bricked_sse3},
        {1, false, true,  false, SSE3, raycaster_uint8_t_preintegration_sse3},
        {2, false, true,  false, SSE3, raycaster_uint16_t_preintegration_sse3},
        {1, true,  true,  false, SSE3, raycaster_uint8_t_bricked_preintegration_sse3},
        {2, true,  true,  false, SSE3, raycaster_uint16_t_bricked_preintegration_sse3},
        {1, false, false, true,  SSE3, raycaster_uint8_t_mip_sse3},
        {2, false, false, true,  SSE3, raycaster_uint16_t_mip_sse3},
        {1, true,  false, true,  SSE3, raycaster_uint8_t_bricked_mip_sse3},
        {2, true,  false, true,  SSE3, raycaster_uint16_t_bricked_mip_sse3},
#endif
    };

    /**
     * Find the raycaster with the specified properties if it is found
     * in the raycasters table.
     */
    RayCasterFunction GetRayCasterFunction(int bytes_per_voxel, bool bricked, bool preintegration,
                                           bool mip, SSELevel sselevel){
        for(int i = 0; i != sizeof(raycasters)/sizeof(raycasters[0]); i++){
            RayCasterConf c = raycasters[i];
            if (c.bytes_per_voxel != bytes_per_voxel) continue;
            if (c.preintegration  != preintegration)  continue;
            if (c.bricked         != bricked)         continue;
            if (c.mip             != mip)             continue;
            if (c.sselevel        != sselevel)        continue;

            return c.function;
        }
        return 0;
    }
}
}

using namespace voreen::simdraycaster;

voreen::SIMDRayCaster::SIMDRayCaster()
    : performanceInfoPort_(Port::OUTPORT, "performanceInfoPort", "Performance info port")
    , transferFunc_("transferFunction", "Transfer Function", Processor::INVALID_RESULT)
    , camera_("camera", "Camera", tgt::Camera(vec3(0.f, 0.f, 3.5f), vec3(0.f, 0.f, 0.f), vec3(0.f, 1.f, 0.f)), true)
    , enableBrickedVolume_("brickedVolume", "Use Bricked Volome", false, Processor::INVALID_RESULT)
    , preintegration_("preintegration", "Use Preintegration", false, Processor::INVALID_RESULT)
    , mipRendering_("miprendering", "Use Mip Rendering")
    , pixelBlockConfiguration_("blockSize", "Pixel group size", tgt::vec2(16), tgt::vec2(1), tgt::vec2(512))
    , brickSize_("brickSize", "Brick size", 8, 1, 32)

{
        volumeInport_.showTextureAccessProperties(true);

        addPort(performanceInfoPort_);
        addProperty(transferFunc_);
        addProperty(camera_);
        
        
        mipRendering_.addOption("mip", "MIP", 1);
        mipRendering_.addOption("dvr", "DVR", 0);
        addProperty(mipRendering_);
        mipRendering_.setGroupID("options");
        addProperty(preintegration_);
        preintegration_.setGroupID("options");
        setPropertyGroupGuiName("options", "Raycasting settings");

        addProperty(pixelBlockConfiguration_);

        addProperty(enableBrickedVolume_);
        enableBrickedVolume_.setGroupID("bricked");
        addProperty(brickSize_);
        brickSize_.setGroupID("bricked");
        setPropertyGroupGuiName("bricked", "Bricked volume");


        // reseting of performance metrics
        volumeInport_.onChange(MemberFunctionCallback<PerformanceMetric>(&performanceMetric_, &PerformanceMetric::clearRuns));
        enableBrickedVolume_.onChange(MemberFunctionCallback<PerformanceMetric>(&performanceMetric_, &PerformanceMetric::clearRuns));
        preintegration_.onChange(MemberFunctionCallback<PerformanceMetric>(&performanceMetric_, &PerformanceMetric::clearRuns));
        mipRendering_.onChange(MemberFunctionCallback<PerformanceMetric>(&performanceMetric_, &PerformanceMetric::clearRuns));
        pixelBlockConfiguration_.onChange(MemberFunctionCallback<PerformanceMetric>(&performanceMetric_, &PerformanceMetric::clearRuns));
        brickSize_.onChange(MemberFunctionCallback<PerformanceMetric>(&performanceMetric_, &PerformanceMetric::clearRuns));

        // invalidation with new volume
        ON_CHANGE_LAMBDA(brickSize_, [this]{
            delete brickedVolume_;
            brickedVolume_ = 0;
        });

        // property visibility
        ON_CHANGE_LAMBDA(enableBrickedVolume_, [this]{
            brickSize_.setVisibleFlag(enableBrickedVolume_.get());
        });
        ON_CHANGE_LAMBDA(mipRendering_, [this]{
            preintegration_.setVisibleFlag(!mipRendering_.getValue());
        });

        brickSize_.setVisibleFlag(enableBrickedVolume_.get());
        preintegration_.setVisibleFlag(!mipRendering_.getValue());

        brickedVolume_ = 0;
}

voreen::SIMDRayCaster::~SIMDRayCaster(){
    if (brickedVolume_)
        delete brickedVolume_;
}

void voreen::SIMDRayCaster::initialize() {
    RenderProcessor::initialize();

    glGenTextures(1, &resultTexure_);
    glBindTexture(GL_TEXTURE_2D, resultTexure_);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    volumeInport_.onNewData(MemberFunctionCallback<SIMDRayCaster>(this, &SIMDRayCaster::adaptToNewVolume));
}
void voreen::SIMDRayCaster::deinitialize() {
    RenderProcessor::deinitialize();
    if (resultTexure_){
        glDeleteTextures(1, &resultTexure_);
    }
}



void voreen::SIMDRayCaster::process(){
    transferFunc_.setVolume(volumeInport_.getData());

    tgtAssert(CPUCapabilities::getRef().hasSSE3(), "SIMDRayCaster needs at least SSE3");
#if defined SIMD_SSE3 && defined SIMD_SSE41
    SSELevel sselevel = CPUCapabilities::getRef().hasSSE41()?SSE41:SSE3;
#elif defined SIMD_SSE3
    SSELevel sselevel = SSE3;
#elif defined SIMD_SSE41
    SSELevel sselevel = SSE41;
#else
    #error "SIMDRaycaster needs sse"
#endif
    bool bricked        = enableBrickedVolume_.get();
    bool preintegration = preintegration_.get();
    bool mip            = mipRendering_.getValue();

    transferFunc_.setVolume(volumeInport_.getData());

    if (mip) {
        // mip does not use preintegration (or even any integration at all)
        preintegration = false;
    }

    ivec2 size     = entryPort_.getSize();
    float stepsize = 1.0f/samplingRate_.get();

    const VolumeRAM* volume          = volumeInport_.getData()->getRepresentation<VolumeRAM>();
    tgt::svec3       volDim          = volume->getDimensions();
    size_t           bytes_per_voxel = volume->getBytesPerVoxel();
    int              pixelCount      = size.x*size.y;
    auto             vf              = volume->getFormat();

    // Allocate memory for raycasting
    uint32_t * img         = aligned_malloc<uint32_t>(pixelCount);
    tgt::vec4* entryBuffer = aligned_malloc<vec4>(pixelCount);
    tgt::vec4* exitBuffer  = aligned_malloc<vec4>(pixelCount);

    // Get data from entry and exit ports
    entryPort_.getColorTexture()->downloadTextureToBuffer(GL_RGBA, GL_FLOAT,
        reinterpret_cast<GLubyte*>(entryBuffer), sizeof(vec4)*pixelCount);
    exitPort_.getColorTexture()->downloadTextureToBuffer(GL_RGBA, GL_FLOAT,
        reinterpret_cast<GLubyte*>(exitBuffer), sizeof(vec4)*pixelCount);

    

    if (bricked && !brickedVolume_){
        brickedVolume_ = createBrickedVolume(volume, brickSize_.get());
    }

    size_t tfsize;
    vec4* tf; // NOT ALIGNED!
    if (!preintegration){
        tf = getTransFunc1D(&tfsize);
        if (!mip){
            // in case of mip rendering we don't need
            // opacityCorrection and alpha premultiplicatoin
            precalcOpacityCorrectionAndAlphaPremultiplication(tf, tfsize, volDim, stepsize);
        }
    } else /* => if (preintegration) */ {
        tf = getPreintegrateTransFunc(&tfsize, stepsize/tgt::max(volDim));
        premultiplyAlpha(tf, tfsize*tfsize);
    }

    RayCasterFunction function = GetRayCasterFunction((int)bytes_per_voxel, bricked, preintegration, mip, sselevel);

    // distribute works into jobs for multiple threads
    JobQueue* queue = &globalJobQueue;

    RayCasterParameters job;
    job.bytesPerVoxel = (int)bytes_per_voxel;
    job.entryBuffer   = entryBuffer;
    job.exitBuffer    = exitBuffer;

    RealWorldMapping rwm = volumeInport_.getData()->getRealWorldMapping();
    tgt::vec2 domain = transferFunc_.get()->getDomain();
    job.tf = tf;
    
    job.tfsize = (int)tfsize;

    // Combine real world mapping and domain into a single affine mapping
    job.tfoffset = (rwm.getOffset()-domain.x)/(domain.y-domain.x)*(job.tfsize-1);
    job.tfscale = rwm.getScale()/(domain.y-domain.x);

    job.img = img;
    job.volDim = volDim;
    if (!bricked){
        job.vol = (void*)volume->getData();
    }else{
        job.vol     = brickedVolume_->getData();
        job.offsetx = brickedVolume_->getOffsetX();
        job.offsety = brickedVolume_->getOffsetY();
        job.offsetz = brickedVolume_->getOffsetZ();
    }

    job.stepsize = stepsize;

    ivec2 blocksize = pixelBlockConfiguration_.get();

    
    performanceMetric_.beginRun();
    // run jobs
    for(int x = 0; x < size.x; x += blocksize.x){
        for(int y = 0; y < size.y; y += blocksize.y){
            int lenx       = std::min(x+blocksize.x, size.x)-x;
            int leny       = std::min(y+blocksize.y, size.y)-y;
            job.offset     = x + y * size.x;
            job.linelength = lenx;
            job.stride     = size.x-lenx;
            job.count      = leny;
            queue->addFunctionJob(function, job);
        }
    }
    // wait for all jobs to finish
    queue->waitAll();

    performanceMetric_.endRun();

    // copy image to outport
    blitImage(img, size);


    // free allocated memory
    aligned_free(img);
    aligned_free(entryBuffer);
    aligned_free(exitBuffer);
    aligned_free(tf);
}

tgt::vec4* voreen::SIMDRayCaster::getTransFunc1D(size_t *size){
    Texture *tftex  = transferFunc_.get()->getTexture();
    size_t   tfsize = tftex->getDimensions().x;
    vec4*    tf     = aligned_malloc<vec4>(tfsize);
    
    tftex->downloadTextureToBuffer(GL_RGBA, GL_FLOAT,
        reinterpret_cast<GLubyte*>(tf), sizeof(vec4)*tfsize);

    *size = tfsize;
    return tf;
}

tgt::vec4* voreen::SIMDRayCaster::getPreintegrateTransFunc(size_t *size, float samplingrate){
    TransFunc1DKeys *          tr      = const_cast<TransFunc1DKeys*>(dynamic_cast<const TransFunc1DKeys*>(transferFunc_.get()));
    const PreIntegrationTable* pr      = tr->getPreIntegrationTable(samplingrate, 0, true, false);
    size_t                     ptxsize = pr->getDimension();
    vec4 *                     ptx     = aligned_malloc<vec4>(ptxsize*ptxsize);

    memcpy(ptx, pr->getTable(), sizeof(vec4) * ptxsize * ptxsize);

    *size = ptxsize;
    return ptx;
}

void voreen::SIMDRayCaster::premultiplyAlpha(tgt::vec4* data, size_t elements){
    for(size_t i = 0; i != elements; i++){
        vec4 val = data[i];
        if (val.a < 0.0f) val.a = 0.0f;
        if (val.a > 1.0f) val.a = 1.0f;

        val.xyz() *= val.a;
        data[i] = val;
    }
}

void voreen::SIMDRayCaster::precalcOpacityCorrectionAndAlphaPremultiplication(vec4* data, size_t elements,
                                                     tgt::svec3 volumeDimension,
                                                     float stepsize){
    // replicate the way opacity correction works in the
    // gpu raycasters
    //
    // This calculates the opacity correction for a stepsize
    // of 0.005 * the maximum dimension of the volume
    const float SAMPLING_BASE_INTERVAL_RCP = 200.0f;
    float normalizedStepSize = stepsize/(tgt::max(volumeDimension));
    float correctionExponent = normalizedStepSize * SAMPLING_BASE_INTERVAL_RCP;
    for(size_t i = 0; i != elements; i++){
        vec4 val = data[i];
        if (val.a < 0.0f) val.a = 0.0f;
        if (val.a > 1.0f) val.a = 1.0f;

        val.a = 1.0f - pow(1.0f - val.a, correctionExponent);
        val.xyz() *= val.a;
        data[i] = val;
    }
}

void voreen::SIMDRayCaster::blitImage(std::uint32_t *img, ivec2 size){
    outport1_.activateTarget();
    outport1_.clearTarget();
    tgt::Shader::deactivate();

    // Load buffer in Texture
    glBindTexture(GL_TEXTURE_2D, resultTexure_);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, size.x, size.y, 0, GL_RGBA, GL_UNSIGNED_BYTE, img);

    // Render to framebuffer
    renderQuad();

    outport1_.deactivateTarget();
}

BrickedVolumeBase * voreen::SIMDRayCaster::createBrickedVolume(const VolumeRAM * vol, int size){
    tgt::svec3 volDim = vol->getDimensions();
    size_t bytes_per_voxel = vol->getBytesPerVoxel();
    if (bytes_per_voxel == 1){
        return new BrickedVolume<uint8_t>((uint8_t*)vol->getData(), volDim, size);
    }
    if (bytes_per_voxel == 2){
        return new BrickedVolume<uint16_t>((uint16_t*)vol->getData(), volDim, size);
    }

    return 0;
}

void voreen::SIMDRayCaster::adaptToNewVolume(){
    if (brickedVolume_){
        delete brickedVolume_;
    }
    brickedVolume_ = 0;
    camera_.adaptInteractionToScene(volumeInport_.getData()->getBoundingBox().getBoundingBox());
}

void voreen::SIMDRayCaster::afterProcess(){
    RenderProcessor::afterProcess();
    performanceInfoPort_.setData(performanceMetric_.getTextInfo());   
}


bool voreen::SIMDRayCaster::isReady() const{
    if (!volumeInport_.hasData()) return false;
    if (!entryPort_.hasData())    return false;
    if (!exitPort_.hasData())     return false;
    if (!outport1_.isConnected()) return false;
    return true;
}
