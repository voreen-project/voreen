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

#include "grayscale_cl.h"

#include "modules/opencl/openclmodule.h"
#include "voreen/core/voreenapplication.h"

namespace voreen {

using namespace cl;

GrayscaleCL::GrayscaleCL()
    : OpenCLProcessor()
    , saturation_("saturation", "Saturation", 0.0f)
    , inport_(Port::INPORT, "inport", "Image Input")
    , outport_(Port::OUTPORT, "outport", "Image Output")
    , prog_(nullptr)
{
    // register properties and ports:
    addProperty(saturation_);

    addPort(inport_);
    addPort(outport_);
}

Processor* GrayscaleCL::create() const {
    return new GrayscaleCL();
}

bool GrayscaleCL::isDeviceChangeSupported() const {
    return true;
}

void GrayscaleCL::initializeCL() {
    prog_ = new Program(OpenCLModule::getInstance()->getCLContext());
    prog_->loadSource(VoreenApplication::app()->getModulePath("opencl") + "/cl/grayscale.cl");
    prog_->build(OpenCLModule::getInstance()->getCLDevice());
}

void GrayscaleCL::deinitializeCL() {
    delete prog_;
    prog_ = nullptr;
}

void GrayscaleCL::process() {
    if (prog_) {
        Kernel* k = prog_->getKernel("gr");
        if(k) {
            const bool sharing = OpenCLModule::getInstance()->getGLSharing();

            cl::Context* context = OpenCLModule::getInstance()->getCLContext();
            cl::CommandQueue* commandQueue = OpenCLModule::getInstance()->getCLCommandQueue();
            tgtAssert(context, "No OpenCL context");
            tgtAssert(commandQueue, "No OpenCL command queue");

            glFinish();

            std::unique_ptr<cl::MemoryObject> in;
            std::unique_ptr<cl::MemoryObject> out;
            if(sharing) {
                in.reset(new SharedTexture(context, CL_MEM_READ_ONLY, inport_.getColorTexture()));
                out.reset(new SharedTexture(context, CL_MEM_WRITE_ONLY, outport_.getColorTexture()));
            } else {
                inport_.getColorTexture()->downloadTexture();
                in.reset(new cl::ImageObject2D(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, inport_.getColorTexture()));
                out.reset(new cl::ImageObject2D(context, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, outport_.getColorTexture()));
            }

            k->setArg(0, in.get());
            k->setArg(1, out.get());
            k->setArg(2, saturation_.get());

            if(sharing) {
                commandQueue->enqueueAcquireGLObject(in.get());
                commandQueue->enqueueAcquireGLObject(out.get());
            }
            commandQueue->enqueue(k, inport_.getSize());
            if(sharing) {
                commandQueue->enqueueReleaseGLObject(in.get());
                commandQueue->enqueueReleaseGLObject(out.get());
            } else {
                cl::ImageObject2D* outPtr = dynamic_cast<cl::ImageObject2D*>(out.get());
                tgtAssert(outPtr, "Expected image ImageObject2D");
                commandQueue->enqueueReadImage(*outPtr, outport_.getColorTexture(), true);
                outport_.getColorTexture()->uploadTexture();
            }

            commandQueue->finish();

            outport_.validateResult();
            outport_.invalidatePort();
        } else {
            LWARNING("Kernel \"gr\" not found");
        }
    }
}

} // namespace
