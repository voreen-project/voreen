/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include <cmath>
#include <ctime>
#include <queue>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "modules/opencl/processors/raycaster_cl.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/properties/buttonproperty.h"

#include "tgt/glmath.h"
#include "tgt/gpucapabilities.h"

#include "modules/opencl/openclmodule.h"

using tgt::vec3;
using tgt::mat4;

namespace voreen {

using namespace cl;

const std::string RaycasterCL::loggerCat_("voreen.RaycasterCL");

RaycasterCL::RaycasterCL()
    : OpenCLProcessor(),
      volumePort_(Port::INPORT, "volumehandle.volumehandle", "Volume Input", false, VolumeRenderer::INVALID_PROGRAM),
      entryPort_(Port::INPORT, "image.entrypoints", "Entry-points Input", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN),
      exitPort_(Port::INPORT, "image.exitpoints", "Exit-points Input", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN),
      // TODO: depth formats like GL_DEPTH_COMPONENT24 are not supported by OpenCL, investigate. Depth values are not calculated for now. FL
      outport_(Port::OUTPORT, "image.output", "Image Output", true, Processor::INVALID_PROGRAM, RenderPort::RENDERSIZE_RECEIVER),
      samplingRate_("samplingRate", "Sampling Rate", 4.f, 0.01f, 20.0f),
      transferFunc_("transferFunction", "Transfer Function"),
      openclProp_("openclprog", "OpenCL Program", VoreenApplication::app()->getBasePath() + "/modules/opencl/cl/raycaster.cl"),
      opencl_(0),
      context_(0),
      queue_(0),
      volumeTex_(0),
      tfData_(0),
      entryTexCol_(0),
      exitTexCol_(0),
      outCol_(0)
{
    addPort(volumePort_);
    addPort(entryPort_);
    addPort(exitPort_);

    outport_.onSizeReceiveChange<RaycasterCL>(this, &RaycasterCL::portResized);
    addPort(outport_);

    addProperty(samplingRate_);
    transferFunc_.onChange(MemberFunctionCallback<RaycasterCL>(this, &RaycasterCL::updateTF));
    addProperty(transferFunc_);
    addProperty(openclProp_);

}

RaycasterCL::~RaycasterCL() {
}

Processor* RaycasterCL::create() const {
    return new RaycasterCL();
}

void RaycasterCL::initialize() {
    OpenCLProcessor::initialize();

    portResized();
    updateTF();
}

void RaycasterCL::deinitialize() {
    OpenCLProcessor::deinitialize();
}

void RaycasterCL::initializeCL() {
    opencl_ = OpenCLModule::getInstance()->getOpenCL();
    context_ = OpenCLModule::getInstance()->getCLContext();
    queue_ = OpenCLModule::getInstance()->getCLCommandQueue();
}

void RaycasterCL::deinitializeCL() {
    delete volumeTex_;
    volumeTex_ = 0;
    delete tfData_;
    tfData_ = 0;
    delete entryTexCol_;
    delete exitTexCol_;
    delete outCol_;
    //delete entryTexDepth_;
    //delete exitTexDepth_;
    //delete outDepth_;
    entryTexCol_ = nullptr;
    exitTexCol_ = nullptr;
    outCol_ = nullptr;

    // Reset pointers
    opencl_ = nullptr;
    context_ = nullptr;
    queue_ = nullptr;
}

void RaycasterCL::beforeProcess() {
    VolumeRenderer::beforeProcess();

    // compile program if needed
    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM) {
        std::ostringstream clDefines;
        clDefines << " -cl-fast-relaxed-math -cl-mad-enable";

        // set include path for modules
        clDefines << " -I" << tgt::FileSystem::cleanupPath(VoreenApplication::app()->getBasePath() + "/modules/opencl/cl/");

        openclProp_.setDefines(clDefines.str());
        openclProp_.rebuild();

        delete volumeTex_;
        volumeTex_ = new ImageObject3D(context_, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, volumePort_.getData()->getRepresentation<VolumeRAM>());
        transferFunc_.setVolume(volumePort_.getData());
    }

    LGL_ERROR;

    if (getInvalidationLevel() >= Processor::INVALID_RESULT) {
        if(!entryTexCol_) {
            entryTexCol_ = new SharedTexture(context_, CL_MEM_READ_ONLY, entryPort_.getColorTexture());
            exitTexCol_  = new SharedTexture(context_, CL_MEM_READ_ONLY, exitPort_.getColorTexture());
            outCol_      = new SharedTexture(context_, CL_MEM_WRITE_ONLY, outport_.getColorTexture());
            //entryTexDepth_ = new SharedTexture(context_, CL_MEM_READ_ONLY, entryPort_.getDepthTexture());
            //exitTexDepth_  = new SharedTexture(context_, CL_MEM_READ_ONLY, exitPort_.getDepthTexture());
            //outDepth_      = new SharedTexture(context_, CL_MEM_READ_ONLY, outport_.getDepthTexture());
        }
    }

}

void RaycasterCL::updateTF() {
    if(!context_)
        return;

    tgt::Texture* tfTex = transferFunc_.get()->getTexture();
    delete tfData_;
    ImageFormat imgf(CL_RGBA, CL_UNORM_INT8);
    tfData_ = new ImageObject2D(context_, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, imgf, tfTex->getWidth(), 1, tfTex->getWidth() * tfTex->getBpt(), transferFunc_.get()->getTexture()->getCpuTextureData());
}


void RaycasterCL::process() {

    LGL_ERROR;

    Kernel* kernel = openclProp_.getProgram()->getKernel("raycast");

    if (kernel) {

        // TODO: Does this really have to be done every time to see the image on the canvas? FL
        outport_.activateTarget();
        outport_.deactivateTarget();

        LGL_ERROR;
        glFinish();

        tgt::svec3 dims = tgt::svec3(outport_.getColorTexture()->getDimensions());
        float samplingStepSize = 1.f / (tgt::min(volumePort_.getData()->getRepresentation<VolumeRAM>()->getDimensions()) * samplingRate_.get());

        kernel->setArg(0, volumeTex_);
        kernel->setArg(1, tfData_);
        kernel->setArg(2, entryTexCol_);
        kernel->setArg(3, exitTexCol_);
        kernel->setArg(4, outCol_);
        //kernel->setArg(5, entryTexDepth_);
        //kernel->setArg(6, exitTexDepth_);
        //kernel->setArg(7, outDepth_);
        kernel->setArg(5, samplingStepSize);
        kernel->setArg(6, static_cast<tgt::ivec2>(outport_.getColorTexture()->getDimensions().xy()));

        queue_->enqueueAcquireGLObject(entryTexCol_);
        queue_->enqueueAcquireGLObject(exitTexCol_);
        queue_->enqueueAcquireGLObject(outCol_);

        //queue_->enqueueAcquireGLObject(&entryTexDepth_);
        //queue_->enqueueAcquireGLObject(&exitTexDepth_);
        //queue_->enqueueAcquireGLObject(&outDepth_);

        queue_->enqueue(kernel, dims.xy() + tgt::svec2(8 - dims.x % 8, 8 - dims.y % 8), tgt::ivec2(8));

        queue_->enqueueReleaseGLObject(entryTexCol_);
        queue_->enqueueReleaseGLObject(exitTexCol_);
        queue_->enqueueReleaseGLObject(outCol_);

        //queue_->enqueueReleaseGLObject(&entryTexDepth_);
        //queue_->enqueueReleaseGLObject(&exitTexDepth_);
        //queue_->enqueueReleaseGLObject(&outDepth_);
        queue_->finish();
    }
    else {
        LERROR("Kernel 'raycast' not found");
        return;
    }
}

bool RaycasterCL::isReady() const {

    if (!OpenCLProcessor::isReady())
        return false;

    //check if all inports are connected
    if(!volumePort_.isReady() || !entryPort_.isReady() || !exitPort_.isReady())
        return false;

    if(!outport_.isReady())
        return false;

    return true;
}

void RaycasterCL::portResized() {
    // it is important to delete the CL memory objects before deleting the associated OpenGL textures from the fbo
    delete entryTexCol_; entryTexCol_ = 0;
    delete exitTexCol_;  exitTexCol_  = 0;
    delete outCol_;      outCol_      = 0;

    entryPort_.requestSize(outport_.getReceivedSize());
    exitPort_.requestSize(outport_.getReceivedSize());
}

} // namespace voreen
