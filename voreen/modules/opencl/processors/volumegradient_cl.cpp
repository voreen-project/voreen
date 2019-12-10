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

#include "volumegradient_cl.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "voreen/core/processors/volumerenderer.h"

#include "modules/opencl/openclmodule.h"
#include "voreen/core/voreenapplication.h"

namespace voreen {

using namespace cl;

const std::string VolumeGradientCL::loggerCat_("voreen.VolumeGradientCL");

VolumeGradientCL::VolumeGradientCL()
    : OpenCLProcessor(),
    technique_("technique", "Technique", VolumeRenderer::INVALID_PROGRAM),
    outputType_("output_type", "Output Datatype"),
    copyIntensityChannel_("copyIntensityChannel", "Copy intensity channel", false),
    inport_(Port::INPORT, "volumehandle.input", "Volume Input", false, VolumeRenderer::INVALID_PROGRAM),
    outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false),
    opencl_(0),
    context_(0),
    queue_(0),
    prog_(0),
    volumeTex_(0)
{
    technique_.addOption("central-differences", "Central differences");
    technique_.addOption("filtered-central", "Filtered central differences");
    //technique_.addOption("sobel", "Sobel");
    //technique_.addOption("linear-regression", "Linear regression");

    outputType_.addOption("uint8", "8bit signed integer");
    outputType_.addOption("uint16", "16bit signed integer");
    outputType_.addOption("float", "32bit float");

    addProperty(technique_);
    addProperty(outputType_);
    addProperty(copyIntensityChannel_);

    addPort(inport_);
    addPort(outport_);
}

Processor* VolumeGradientCL::create() const {
    return new VolumeGradientCL();
}

void VolumeGradientCL::initializeCL() {

    opencl_ = OpenCLModule::getInstance()->getOpenCL();
    context_ = OpenCLModule::getInstance()->getCLContext();
    queue_ = OpenCLModule::getInstance()->getCLCommandQueue();

    std::string kernelFile = VoreenApplication::app()->getModulePath("opencl") + "/cl/gradient.cl";
    LINFO("Loading program " << kernelFile);
    prog_ = new Program(OpenCLModule::getInstance()->getCLContext());
    prog_->loadSource(kernelFile);
}

void VolumeGradientCL::deinitializeCL() {
    delete prog_;
    prog_ = 0;
}

void VolumeGradientCL::beforeProcess() {
    VolumeProcessor::beforeProcess();

    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM) {
        std::ostringstream clDefines;
        clDefines << " -cl-fast-relaxed-math -cl-mad-enable";

        // set include path for modules
        // TODO:  replace slashes with double-backslashes
        std::string inc = VoreenApplication::app()->getBasePath() + "/modules/opencl/cl/";
        clDefines << " -I" << inc << " ";

        if(technique_.get() == "central-differences")
            clDefines << " -DGRADIENT_CENTRAL ";
        else if(technique_.get() == "filtered-central")
            clDefines << " -DGRADIENT_FILTERED ";
        prog_->setBuildOptions(clDefines.str());

        bool success = prog_->build(OpenCLModule::getInstance()->getCLDevice());
        if (!success)
            throw VoreenException("Unable to build program");

        if (!prog_->getKernel("gradient"))
            throw VoreenException("Kernel 'gradient' not found");

        delete volumeTex_;
        volumeTex_ = new ImageObject3D(context_, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, inport_.getData()->getRepresentation<VolumeRAM>());
    }
}

void VolumeGradientCL::process() {

    tgtAssert(prog_, "No program");

    const VolumeRAM* inputVolume = inport_.getData()->getRepresentation<VolumeRAM>();
    VolumeRAM* outputVolume = 0;

    Kernel* kernel = prog_->getKernel("gradient");
    if (kernel) {
        glFinish();

        delete outport_.getData();
        outport_.setData(0);

        if(outputType_.get() == "uint8") {
            if (copyIntensityChannel_.get())
                outputVolume = new VolumeRAM_4xUInt8(inputVolume->getDimensions());
            else
                outputVolume = new VolumeRAM_3xUInt8(inputVolume->getDimensions());
        } else if(outputType_.get() == "uint16") {
            if (copyIntensityChannel_.get())
                outputVolume = new VolumeRAM_4xUInt16(inputVolume->getDimensions());
            else
                outputVolume = new VolumeRAM_3xUInt16(inputVolume->getDimensions());
        } else if(outputType_.get() == "float") {
            if (copyIntensityChannel_.get())
                outputVolume = new VolumeRAM_4xFloat(inputVolume->getDimensions());
            else
                outputVolume = new VolumeRAM_3xFloat(inputVolume->getDimensions());
        }

        VolumeWriteBuffer volBuffer(context_, outputVolume);

        kernel->setArg(0, volumeTex_);
        // At the moment, the info-Buffer and the data-Buffer have to be passed to the kernel as separate arguments.
        // They must be re-combined within the kernel manually before using the write_volume functions.
        kernel->setArg(1, volBuffer.infoBuffer_);
        kernel->setArg(2, volBuffer.dataBuffer_);
        kernel->setArg(3, copyIntensityChannel_.get() ? 1 : 0);

        queue_->enqueue(kernel, inputVolume->getDimensions());
        queue_->enqueueReadBuffer(volBuffer.dataBuffer_, outputVolume->getData(), true);

        queue_->finish();
    }
    else {
        LERROR("Kernel 'gradient' not found");
    }

    // assign computed volume to outport
    if (outputVolume)
        outport_.setData(new Volume(outputVolume, inport_.getData()));
    else
        outport_.setData(0);

}

}   // namespace
