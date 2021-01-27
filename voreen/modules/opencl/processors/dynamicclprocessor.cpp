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

#include <algorithm>

#include "modules/opencl/processors/dynamicclprocessor.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/datastructures/transfunc/1d/1dkeys/transfunc1dkeys.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/matrixproperty.h"
#include "voreen/core/processors/processorwidget.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

#include "modules/opencl/openclmodule.h"

namespace voreen {

using namespace cl;

const std::string DynamicCLProcessor::loggerCat_("voreen.DynamicCLProcessor");
std::vector<std::string> DynamicCLProcessor::storeModsIgn_ = std::vector<std::string>();
std::vector<std::string> DynamicCLProcessor::storeMods_ = std::vector<std::string>();
std::vector<std::string> DynamicCLProcessor::clTypesIgn_ = std::vector<std::string>();
std::vector<std::string> DynamicCLProcessor::clTypes_ = std::vector<std::string>();
std::vector<std::string> DynamicCLProcessor::tagVals_ = std::vector<std::string>();
std::vector<std::string> DynamicCLProcessor::numericTagVals_ = std::vector<std::string>();

std::map<std::string, DynamicCLProcessor::AttributeType> DynamicCLProcessor::propertyMap_ = std::map<std::string, DynamicCLProcessor::AttributeType>();

DynamicCLProcessor::DynamicCLProcessor()
    : OpenCLProcessor()
    , openclProp_("openclprog", "OpenCL Program", VoreenApplication::app()->getModulePath("opencl") + "/cl/simple.cl")
    , opencl_(nullptr)
    , context_(nullptr)
    , queue_(nullptr)
    , curArgs_(std::vector<ArgInfo>())
    , curSharedTexs_(std::vector<SharedTexture*>())
    , curVolTexs_(std::vector<ImageObject3D*>())
    , curTFs_(std::vector<ImageObject2D*>())
    , regenerateSharedTextures_(false)
    , portArgMap_(std::map<Port*, ArgInfo*>())
    , curWorkDimensionSource_(0)
    , wasDeserialized_(false)
{
    if(clTypes_.size() == 0)
        setupKeywords();

    std::ostringstream clDefines;
    clDefines << " -cl-fast-relaxed-math -cl-mad-enable";

    // set include path for modules
    clDefines << " -I" << VoreenApplication::app()->getBasePath() + "/modules/opencl/cl/";

    openclProp_.setDefines(clDefines.str());
    openclProp_.onChange(MemberFunctionCallback<DynamicCLProcessor>(this, &DynamicCLProcessor::buildAndInit));
}

DynamicCLProcessor::~DynamicCLProcessor() {
    for (Port* port : dynamicPorts_) {
        port->deinitialize();
        delete port;
    }
    dynamicPorts_.clear();

    for (Property* property : dynamicProperties_)
        delete property;
    dynamicProperties_.clear();
}

Processor* DynamicCLProcessor::create() const {
    return new DynamicCLProcessor();
}

void DynamicCLProcessor::initialize() {

    openclProp_.initialize();

    // In case we deserialized the program, there is no need to parse the program, since
    // it already has been performed in order to deserialization the connections between the ports.
    // Is has to, however, if we create a fresh instance of this processor to initialize the default source.
    if (!wasDeserialized_) {
        openclProp_.resetProgramSource(); // set's programSource
        parseProgram();
    }

    // Initialize processor afterwards, which will initialize the ports and properties having been added while parsing.
    OpenCLProcessor::initialize();
}

void DynamicCLProcessor::deinitialize() {
    OpenCLProcessor::deinitialize();

    openclProp_.deinitialize();
}

void DynamicCLProcessor::initializeCL() {
    opencl_ = OpenCLModule::getInstance()->getOpenCL();
    context_ = OpenCLModule::getInstance()->getCLContext();
    queue_ = OpenCLModule::getInstance()->getCLCommandQueue();

    buildProgram();
}

void DynamicCLProcessor::deinitializeCL() {
    clearSharedTextures();
    clearTransferFunctions();
    clearVolumeRepresentations();

    // Reset pointers.
    opencl_ = nullptr;
    context_ = nullptr;
    queue_ = nullptr;
}

bool DynamicCLProcessor::isDeviceChangeSupported() const {
    return true;
}

void DynamicCLProcessor::buildProgram() {
    if(!opencl_)
        return;

    openclProp_.rebuild();
    getProcessorWidget()->updateFromProcessor();
}

void DynamicCLProcessor::parseProgram() {

    std::string source = openclProp_.get().programSource_;
    if (source.empty() && openclProp_.getProgram()) {
        source = openclProp_.getProgram()->getInfo<std::string>(CL_PROGRAM_SOURCE);
    }

    if (source.empty()) {
        LWARNING("No Program Source available.");
        return;
    }

    LINFO("Parsing program source...");
    parseArguments(source);
    assignAttributes();
    processAnnotations();
    setInvalidationLevels();
    setWorkDimensions();

    LINFO("Building ports and properties...");
    removeObsoletePortsAndProperties();
    refreshPortsAndProperties();
}

void DynamicCLProcessor::buildAndInit() {

    // Try to build program first.
    try {
        buildProgram();
    }
    catch (const VoreenException& e) {
        LERROR(e.what());
    }

    // Even if build failed, parsing the code could still work.
    std::vector<ArgInfo> tmpArgs = curArgs_;
    try {
        parseProgram();
    }
    catch (const std::exception& e) {
        LERROR(e.what());
        LERROR("Ports and properties could not be updated.");
        curArgs_ = tmpArgs;
        return;
    }

    portResized();
    invalidate(Processor::INVALID_PROGRAM);
}

void DynamicCLProcessor::beforeProcess() {
    OpenCLProcessor::beforeProcess();

    // handle TFs as pure OpenCL image objects as 1d texture sharing is not supported (it seems)
    // There is currently no possibility for the processor to react directly to port connection changes.
    // Therefore, we have to regenerate the shared textures for each rendering at the moment.
    // TODO: Do this only if the invalidation level of the individual argument has been reached.
    if(getInvalidationLevel() >= Processor::INVALID_RESULT) {
        generateTransferFunctions();
        generateSharedTextures();
    }

    if(getInvalidationLevel() >= Processor::INVALID_PROGRAM) {
        generateVolumeRepresentations();
    }

    LGL_ERROR;
}

void DynamicCLProcessor::clearTransferFunctions() {
    for (size_t i = 0; i < curTFs_.size(); i++)
        delete curTFs_.at(i);
    curTFs_.clear();
}

void DynamicCLProcessor::generateTransferFunctions() {
    clearTransferFunctions();

    for(size_t i = 0; i < curArgs_.size(); i++) {
        if(curArgs_.at(i).attType_ == TRANSFUNCPROPERTY) {
            ImageObject2D* tfImg = 0;
            TransFunc1DKeysProperty* tfProp = (TransFunc1DKeysProperty*)(getProperty(curArgs_.at(i).values_["name"]));
            tgt::Texture* tfTex = tfProp->get()->getTexture();
            if(tfTex) {
                ImageFormat imgf(CL_RGBA, CL_UNORM_INT8);
                tfImg = new ImageObject2D(context_, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, imgf, tfTex->getWidth(), 1, tfTex->getWidth() * tfTex->getBpt(), tfProp->get()->getTexture()->getCpuTextureData());
            }
            curTFs_.push_back(tfImg);
        }
    }
}

void DynamicCLProcessor::clearVolumeRepresentations() {
    for (size_t i = 0; i < curVolTexs_.size(); i++)
        delete curVolTexs_.at(i);

    for (size_t i = 0; i < curVolBuffers_.size(); i++)
        delete curVolBuffers_.at(i);

    curVolTexs_.clear();
    curVolBuffers_.clear();
}

void DynamicCLProcessor::generateVolumeRepresentations() {
    clearVolumeRepresentations();

    // TODO: Copying spacing and other stuff from the first volume inport, if it exists and is ready, to volume outports is only a heuristic.  It would be better
    // to provide a mechanism in which in- and outports can be assigend to each other via comment tags in the cl code.
    VolumeBase* inHandle = 0;

    for(size_t i = 0; i < curArgs_.size(); i++) {
        if(curArgs_.at(i).attType_ == IMAGE3D && curArgs_.at(i).att_ == INPORT) {
            ImageObject3D* volObj = 0;
            const VolumeBase* handle  = ((VolumePort*)(getPort(curArgs_.at(i).values_["name"])))->getData();
            if(handle) {
                volObj = new ImageObject3D(context_, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, handle->getRepresentation<VolumeRAM>());
                if(!inHandle)
                    inHandle = const_cast<VolumeBase*>(handle);
            }
            curVolTexs_.push_back(volObj);
        }
    }

    for(size_t i = 0; i < curArgs_.size(); i++) {
        if(curArgs_.at(i).attType_ == IMAGE3D && curArgs_.at(i).att_ == OUTPORT) {

            if (curArgs_.at(i).values_.find("volumedimensions") == curArgs_.at(i).values_.end()) {
                LWARNING("At argument '" << curArgs_.at(i).values_["name"] << "': Please specify the target volume dimensions with the tag 'volumedimensions' using a value of the form XxYxZ.");
                LWARNING("Using default value 32x32x32.");
                curArgs_.at(i).values_["volumedimensions"] = "32x32x32";
            }
            if (curArgs_.at(i).values_.find("volumetype") == curArgs_.at(i).values_.end()) {
                LWARNING("At argument '" << curArgs_.at(i).values_["name"] << "': Please specify the target volume type with the tag 'volumetype', e.g. 'volumetype=4xUInt8.");
                LWARNING("Defaulting to UInt8.");
                curArgs_.at(i).values_["volumetype"] = "UInt8";
            }

            VolumeRAM* vol = generateVolumeFromTags(curArgs_.at(i).values_["volumedimensions"], curArgs_.at(i).values_["volumetype"]);
            VolumeWriteBuffer* vwb = 0;

            if(vol) {
                VolumePort* vp = (VolumePort*)(getPort(curArgs_.at(i).values_["name"]));
                vp->setData(nullptr);
                vwb = new VolumeWriteBuffer(context_, vol);
                Volume* vh;
                if(inHandle)
                    vh = new Volume(vol, inHandle);
                else {
                    vh = new Volume(vol, tgt::vec3(1.f), tgt::vec3(0.0f));
                    oldVolumePosition(vh);
                }
                vp->setData(vh);
            }

            curVolBuffers_.push_back(vwb);
        }
    }
}

VolumeRAM* DynamicCLProcessor::generateVolumeFromTags(const std::string& dims, const std::string& type) {

    size_t posX = dims.find_first_of('x', 0);
    if(posX == std::string::npos)
        posX = dims.find_first_of('X', 0);
    if(posX == std::string::npos) {
        LERROR("Please use the tag 'volumedimensions' with a value of the form XxYxZ.");
        return 0;
    }

    int valX, valY, valZ;
    std::stringstream sX;
    sX << dims.substr(0, posX);
    sX >> valX;

    if(sX.rdstate() & std::ifstream::failbit) {
        LERROR("Could not read value for tag 'volumedimensions'.");
        return 0;
    }

    size_t posY = dims.find('x', posX + 1);
    if(posY == std::string::npos)
        posY = dims.find('X', posX + 1);
    if(posY == std::string::npos) {
        LERROR("Please use the tag 'volumedimensions' with a value of the form XxYxZ.");
        return 0;
    }

    std::stringstream sY;
    sY << dims.substr(posX + 1, posY - (posX + 1));
    sY >> valY;

    if(sY.rdstate() & std::ifstream::failbit) {
        LERROR("Could not read value for tag 'volumedimensions'.");
        return 0;
    }

    std::stringstream sZ;
    sZ << dims.substr(posY + 1);
    sZ >> valZ;

    if(sZ.rdstate() & std::ifstream::failbit) {
        LERROR("Could not read value for tag 'volumedimensions'.");
        return 0;
    }

    tgt::ivec3 volumeDims = tgt::ivec3(valX, valY, valZ);
    std::string lowerType(type);
    std::transform(lowerType.begin(), lowerType.end(), lowerType.begin(), ::tolower);

    // TODO: signed types
    if(lowerType == "uint8")
        return new VolumeRAM_UInt8(volumeDims);
    else if(lowerType == "uint16")
        return new VolumeRAM_UInt16(volumeDims);
    //else if(lowerType == "uint32")
        //return new VolumeRAM_UInt32(volumeDims);
    else if(lowerType == "3xuint8")
        return new VolumeRAM_3xUInt8(volumeDims);
    else if(lowerType == "3xuint16")
        return new VolumeRAM_3xUInt8(volumeDims);
    else if(lowerType == "4xuint8")
        return new VolumeRAM_4xUInt8(volumeDims);
    else if(lowerType == "4xuint16")
        return new VolumeRAM_4xUInt16(volumeDims);
    else if(lowerType == "float")
        return new VolumeRAM_Float(volumeDims);
    else if(lowerType == "4xfloat")
        return new VolumeRAM_4xFloat(volumeDims);
    else {
        LERROR("Volume type '" << type << "' is currently not supported.");
        return 0;
    }
}

void DynamicCLProcessor::clearSharedTextures() {
    for (size_t i = 0; i < curSharedTexs_.size(); i++)
        delete curSharedTexs_.at(i);
    curSharedTexs_.clear();
}

void DynamicCLProcessor::generateSharedTextures() {
    clearSharedTextures();

    for(size_t i = 0; i < curArgs_.size(); i++) {
        if(curArgs_.at(i).attType_ == IMAGE2D) {
            SharedTexture* shrTex = 0;
            RenderPort* rp = (RenderPort*)(getPort(curArgs_.at(i).values_["name"]));
            tgt::Texture* colTex = 0;

            // FIXME: there is a voreen bug where render targets are destroyed when a renderport is disconnected.  After reconnection, the renderport has to be ready (which it
            // isn't due to the missing renderport).  Thus, the render processor is not ready, and therefore it's beforeProcess method will not be called.  However, this method reinitializes the port's RenderTarget...
            if(rp && rp->isConnected() && !rp->getRenderTarget())
                LWARNING("RenderTarget was not (re-)initialized for port '" << rp->getID() << "'.");
            else
                colTex = rp->getColorTexture();

            if(colTex)
                shrTex = new SharedTexture(context_, (curArgs_.at(i).att_ == INPORT) ? CL_MEM_READ_ONLY : CL_MEM_WRITE_ONLY, colTex);

            curSharedTexs_.push_back(shrTex);
        }
    }
}

// TODO: ATM, any outport resize will lead to size propagation from the first connected outport to all inports; all outports will be also be resized to that size.
// This should be configurable using tags in the cl code.
void DynamicCLProcessor::portResized() {
    // apparently, it is important to delete the CL memory objects before deleting the associated OpenGL textures from the fbo
    clearSharedTextures();

    regenerateSharedTextures_ = true;

    tgt::ivec2 size = tgt::ivec2(-1);

    for(size_t i = 0; i < curArgs_.size(); i++) {
        if(curArgs_.at(i).attType_ == IMAGE2D && curArgs_.at(i).att_ == OUTPORT) {
            RenderPort* rp = (RenderPort*)(getPort(curArgs_.at(i).values_["name"]));
            if(rp->getReceivedSize() != tgt::ivec2(-1)) {
                size = rp->getReceivedSize();
                break;
            }
        }
    }

    if(size == tgt::ivec2(-1))
        return;

    for(size_t i = 0; i < curArgs_.size(); i++) {
        if(curArgs_.at(i).attType_ == IMAGE2D) {
            RenderPort* rp = (RenderPort*)(getPort(curArgs_.at(i).values_["name"]));
            if(curArgs_.at(i).att_ == INPORT)
                rp->requestSize(size);
            else if(curArgs_.at(i).att_ == OUTPORT) {
                RenderSizeReceiveProperty* srp = rp->getSizeReceiveProperty();
                if(srp)
                    srp->set(size);
            }
        }
    }
}

void DynamicCLProcessor::process() {

    LGL_ERROR;

    if (!openclProp_.getProgram()) {
        LWARNING("No valid program present. Skipping execution.");
        return;
    }

    if(openclProp_.getProgram()->getCurrentKernels().empty()) {
        LWARNING("No valid kernels present. Skipping execution.");
        return;
    }

    if(openclProp_.getProgram()->getCurrentKernels().size() > 1) {
        LERROR("Multiple kernels are currently not supported. Skipping execution.");
        return;
    }

    if(curWorkDimensionSource_ == 0) {
        LWARNING("No work dimension for kernel was specified.  Skipping execution.");
        return;
    }

    Kernel* kernel = (*openclProp_.getProgram()->getCurrentKernels().begin()).second;

    if (kernel) {

        for(size_t i = 0; i < curArgs_.size(); i++) {
            if(curArgs_.at(i).att_ == OUTPORT && curArgs_.at(i).attType_ == IMAGE2D && ((RenderPort*)(getPort(curArgs_.at(i).values_["name"])))->hasRenderTarget()) {
                ((RenderPort*)(getPort(curArgs_.at(i).values_["name"])))->activateTarget();
                glClear(GL_COLOR_BUFFER_BIT);
                ((RenderPort*)(getPort(curArgs_.at(i).values_["name"])))->deactivateTarget();
            }
        }

        LGL_ERROR;
        glFinish();

        unsigned int volTexCounter = 0;
        unsigned int volBufCounter = 0;
        unsigned int shrCounter = 0;
        unsigned int tfCounter  = 0;

        for(size_t i = 0; i < curArgs_.size(); i++) {
            if (curArgs_.at(i).attType_ == IMAGE2D) {
                SharedTexture* shrTex = curSharedTexs_.at(shrCounter);
                if(shrTex)
                    kernel->setArg(static_cast<cl_uint>(i), curSharedTexs_.at(shrCounter));
                else {
                    LWARNING("Image port '" << curArgs_.at(i).values_["name"] << "' not connected, delaying execution.");
                    return;
                }

                shrCounter++;
            } else if(curArgs_.at(i).attType_ == IMAGE3D && curArgs_.at(i).att_ == INPORT) {
                ImageObject3D* volTex = curVolTexs_.at(volTexCounter);
                if(volTex)
                    kernel->setArg(static_cast<cl_uint>(i), volTex);
                else {
                    LWARNING("Volume port '" << curArgs_.at(i).values_["name"] << "' not connected, delaying execution.");
                    return;
                }
                volTexCounter++;
            } else if(curArgs_.at(i).attType_ == IMAGE3D && curArgs_.at(i).att_ == OUTPORT) {
                VolumeWriteBuffer* volBuf = curVolBuffers_.at(volBufCounter);
                if(volBuf) {
                    kernel->setArg(static_cast<cl_uint>(i), volBuf->infoBuffer_);
                    kernel->setArg(static_cast<cl_uint>(i+1), volBuf->dataBuffer_);
                    i++;
                }
                else {
                    LWARNING("Volume port '" << curArgs_.at(i).values_["name"] << "' not initialized with Volume, delaying execution.");
                    return;
                }
                volBufCounter++;
            } else if(curArgs_.at(i).attType_ == UNSUPPORTED && curArgs_.at(i).att_ == UNKNOWN) {
                // do nothing
            } else {
                Property* curProp = getProperty(curArgs_.at(i).values_["name"]);

                if(dynamic_cast<FloatProperty*>(curProp))
                    kernel->setArg(static_cast<cl_uint>(i), dynamic_cast<FloatProperty*>(curProp)->get());
                else if(dynamic_cast<FloatVec2Property*>(curProp))
                    kernel->setArg(static_cast<cl_uint>(i), dynamic_cast<FloatVec2Property*>(curProp)->get());
                else if(dynamic_cast<FloatVec3Property*>(curProp))
                    kernel->setArg(static_cast<cl_uint>(i), dynamic_cast<FloatVec3Property*>(curProp)->get());
                else if(dynamic_cast<FloatVec4Property*>(curProp))
                    kernel->setArg(static_cast<cl_uint>(i), dynamic_cast<FloatVec4Property*>(curProp)->get());
                else if(dynamic_cast<IntProperty*>(curProp))
                    kernel->setArg(static_cast<cl_uint>(i), dynamic_cast<IntProperty*>(curProp)->get());
                else if(dynamic_cast<IntVec2Property*>(curProp))
                    kernel->setArg(static_cast<cl_uint>(i), dynamic_cast<IntVec2Property*>(curProp)->get());
                else if(dynamic_cast<IntVec3Property*>(curProp))
                    kernel->setArg(static_cast<cl_uint>(i), dynamic_cast<IntVec3Property*>(curProp)->get());
                else if(dynamic_cast<IntVec4Property*>(curProp))
                    kernel->setArg(static_cast<cl_uint>(i), dynamic_cast<IntVec4Property*>(curProp)->get());
                else if(dynamic_cast<FloatMat4Property*>(curProp))
                    kernel->setArg(static_cast<cl_uint>(i), dynamic_cast<FloatMat4Property*>(curProp)->get());
                else if(dynamic_cast<TransFunc1DKeysProperty*>(curProp)) {
                    ImageObject2D* tf = curTFs_.at(tfCounter);
                    if(tf)
                        kernel->setArg(static_cast<cl_uint>(i), tf);
                    else
                        LWARNING("Transfer function '" << curArgs_.at(i).values_["name"] << "' not initialized.");
                    tfCounter++;
                }
                // matrices are currently not supported by OpenCL natively.  We can use the type float16 for 4x4 matrices.  Other matrices cannot be passed to the kernel right now.
                //else if(dynamic_cast<FloatMat2Property*>(curProp))
                    //kernel->setArg(i, dynamic_cast<FloatMat2Property*>(curProp)->get());
                //else if(dynamic_cast<FloatMat3Property*>(curProp))
                    //kernel->setArg(i, dynamic_cast<FloatMat3Property*>(curProp)->get());
                //else if(dynamic_cast<ColorProperty*>(curProp))
                    //kernel->setArg(i, dynamic_cast<ColorProperty*>(curProp)->get());
                else
                    LERROR("Type for '" << curArgs_.at(i).values_["name"] << "' can currently not be passed to the kernel.");
            }
        }

        for(size_t i = 0; i < curSharedTexs_.size(); i++) {
            SharedTexture* shrTex = curSharedTexs_.at(i);
            if(shrTex)
                queue_->enqueueAcquireGLObject(curSharedTexs_.at(i));
        }

        std::string wdName = curWorkDimensionSource_->values_["name"];

        if(curWorkDimensionSource_->att_ == INPORT || curWorkDimensionSource_->att_ == OUTPORT) {
            if(curWorkDimensionSource_->attType_ == IMAGE2D) {
                RenderPort* rp = (RenderPort*)(getPort(wdName));
                if(!rp->isConnected() || !rp->getRenderTarget()) {
                    LWARNING("Port '" << wdName << "' is not connected or has no render target, no working dimensions available. Skipping execution.");
                    return;
                }
                queue_->enqueue(kernel, rp->getSize());
            }
            else if(curWorkDimensionSource_->attType_ == IMAGE3D) {
                const VolumeBase* h = ((VolumePort*)getPort(wdName))->getData();
                if(!h) {
                    LWARNING("Port '" << wdName << "' has no volume, no working dimensions available. Skipping execution.");
                    return;
                }
                queue_->enqueue(kernel, h->getDimensions());
            }
            else {
                LWARNING("Type of port '" << wdName << "' cannot be used as work dimension source. Skipping execution.");
                return;
            }
        }
        else if(curWorkDimensionSource_->att_ == PROPERTY) {
            if(curWorkDimensionSource_->attType_ == INTPROPERTY)
                queue_->enqueue(kernel, ((IntProperty*)getProperty(wdName))->get());
            else if(curWorkDimensionSource_->attType_ == INT2PROPERTY)
                queue_->enqueue(kernel, ((IntVec2Property*)getProperty(wdName))->get());
            else if(curWorkDimensionSource_->attType_ == INT3PROPERTY)
                queue_->enqueue(kernel, ((IntVec3Property*)getProperty(wdName))->get());
            else {
                LWARNING("Type of property '" << wdName << "' cannot be used as work dimension source. Skipping execution.");
                return;
            }
        }

        // read into volumes if necessary
        for(size_t i = 0; i < curVolBuffers_.size(); i++) {
            VolumeWriteBuffer* vwb = curVolBuffers_.at(i);
            if(vwb)
                vwb->readBackData(queue_);
                //queue_->enqueueReadBuffer(vwb->dataBuffer_, vwb->volume_->getData(), true);
        }

        for(size_t i = 0; i < curSharedTexs_.size(); i++) {
            SharedTexture* shrTex = curSharedTexs_.at(i);
            if(shrTex)
                queue_->enqueueReleaseGLObject(curSharedTexs_.at(i));
        }

        queue_->finish();
    }
    else {
        LERROR("Kernel not found");
        return;
    }
}

void DynamicCLProcessor::serialize(Serializer& s) const {
    const bool usePointerContentSerialization = s.getUsePointerContentSerialization();
    try {

        std::vector<const Property*> clVec;
        clVec.push_back(&openclProp_);

        // serialize properties
        s.setUsePointerContentSerialization(true);
        s.serialize("PriorityProperty", clVec, "Property");

    } catch (SerializationException &e) {
        LWARNING(std::string("OpenCL program serialization failed: ") + e.what());
    }
    s.setUsePointerContentSerialization(usePointerContentSerialization);

    OpenCLProcessor::serialize(s);
}

void DynamicCLProcessor::deserialize(Deserializer& s) {
    const bool usePointerContentSerialization = s.getUsePointerContentSerialization();
    try {
        std::vector<Property*> clVec;
        clVec.push_back(&openclProp_);

        // deserialize property
        s.setUsePointerContentSerialization(true);
        s.deserialize("PriorityProperty", clVec, "Property");
    } catch (SerializationNoSuchDataException) {
        s.removeLastError();
    }
    s.setUsePointerContentSerialization(usePointerContentSerialization);

    // Parse the program before deserialization because it will create
    // properties and ports. The latter needs to exist before connections can be deserialized.
    parseProgram();
    OpenCLProcessor::deserialize(s);

    // Tell the initialization procedure we have deserialized data.
    wasDeserialized_ = true;
}

void DynamicCLProcessor::parseArguments(const std::string& programSource) {

    curArgs_.clear();

    size_t kernelPos   = programSource.find("__kernel");
    size_t argPosBegin = programSource.find("(", kernelPos);
    size_t argPosEnd   = programSource.find(")", kernelPos);
    std::string args   = programSource.substr(argPosBegin + 1, argPosEnd - argPosBegin - 1);

    // delete ordinary comments
    size_t commentPos = args.find("//");
    while(commentPos != std::string::npos) {
        if(commentPos + 2 < args.size() && args.at(commentPos + 2) == '$') {
            commentPos = args.find("//", commentPos + 2);
            continue;
        }
        // delete comment until end of line, but keep newline token until annotations have been parsed
        size_t lineEnd = args.find("\n", commentPos);
        args.erase(commentPos, lineEnd - commentPos);
        commentPos = args.find("//");
    }

    // parse annotations
    size_t annotationPos = args.find("//$");
    while(annotationPos != std::string::npos) {
        size_t lineEnd = args.find("\n", annotationPos);
        std::string curAnt = args.substr(annotationPos + 3, lineEnd - annotationPos - 3);
        curArgs_.push_back(ArgInfo(curAnt));
        args.erase(annotationPos, lineEnd - annotationPos + 1);
        annotationPos = args.find("//$");
    }

    // get rid of all remaining whitespace
    args.erase(std::remove(args.begin(), args.end(), '\t'), args.end());
    args.erase(std::remove(args.begin(), args.end(), '\n'), args.end());
    args.erase(std::remove(args.begin(), args.end(), '\r'), args.end());
    args.erase(std::remove(args.begin(), args.end(), '\f'), args.end());
    args.erase(std::remove(args.begin(), args.end(), '\v'), args.end());

    // parse actual arguments
    bool done = false;
    std::vector<std::string> tmpArgs;
    std::set<std::string> names;

    while(!done) {
        size_t argEnd = args.find(",");
        if(argEnd == std::string::npos) {
            done = true,
            argEnd = args.find_last_not_of(' ')+1;
        }

        std::string curArg = args.substr(0, argEnd);
        curArg = curArg.substr(args.find_first_not_of(' '), args.find_last_not_of(' ') + 1);

        std::string name = curArg.substr(curArg.find_last_of(' ')+1);

        // Check for duplicates.
        if(names.find(name) != names.end()) {
            throw VoreenException("Found duplicate argument name: " + name);
        }
        names.insert(name);

        tmpArgs.push_back(curArg);
        if(!done)
            args = args.substr(argEnd + 1);
    }

    if(tmpArgs.size() != curArgs_.size()) {
        //LWARNING("Number of annotations is not equal to number of kernel arguments.");
        throw VoreenException("Number of annotations is not equal to number of kernel arguments.");
    }
    else {
        for(size_t i = 0; i < curArgs_.size(); i++)
            curArgs_.at(i).arg_ = tmpArgs.at(i);
    }

}

void DynamicCLProcessor::assignAttributes() {

    for(size_t i = 0; i < curArgs_.size(); i++) {
        ArgInfo* curArg = &(curArgs_.at(i));
        std::string curMod  = "";
        std::string curType = "";

        for(size_t j = 0; j < storeModsIgn_.size(); j++) {
            // remove storage class modifiers from arguments that are irrelevant for the network
            size_t modPos = curArg->arg_.find(storeModsIgn_.at(j));
            if(modPos != std::string::npos)
                curArg->arg_.erase(modPos, storeModsIgn_.at(j).size());
        }

        for(size_t j = 0; j < storeMods_.size(); j++) {
            // remove storage class modifiers from arguments
            size_t modPos = curArg->arg_.find(storeMods_.at(j));
            if(modPos != std::string::npos) {
                curMod = storeMods_.at(j);
                curArg->arg_.erase(modPos, storeMods_.at(j).size());
                break; // there can be only one significant modifier
            }
        }

        for(size_t j = 0; j < clTypesIgn_.size(); j++) {
            // remove storage class modifiers from arguments that are irrelevant for the network
            size_t typePos = curArg->arg_.find(clTypesIgn_.at(j));
            if(typePos != std::string::npos)
                throw VoreenException(std::string("Kernel argument '").append(curArg->arg_).append("': unsupported argument type '").append(clTypesIgn_.at(j)));
        }

        for(size_t j = 0; j < clTypes_.size(); j++) {
            // remove type names from arguments
            size_t typePos = curArg->arg_.find(clTypes_.at(j));
            if(typePos != std::string::npos) {
                curType = clTypes_.at(j);
                curArg->arg_.erase(typePos, clTypes_.at(j).size());
                break; // there can be only one type
            }
        }

        // remove remaining whitespace, which should leave us with only the argument name
        curArg->arg_.erase(std::remove(curArg->arg_.begin(), curArg->arg_.end(), ' '), curArg->arg_.end());
        curArg->registerAttribute(curMod, curType);
    }
}

void DynamicCLProcessor::processAnnotations() {
    for(size_t i = 0; i < curArgs_.size(); i++) {
        ArgInfo* curArg = &(curArgs_.at(i));
        std::string curAnt = curArg->ant_;
        std::replace(curAnt.begin(), curAnt.end(), ',', ' ');
        std::replace(curAnt.begin(), curAnt.end(), ';', ' ');

        size_t antPos = curAnt.find("@");
        while(antPos != std::string::npos) {
            size_t nextAnt = curAnt.find("@", antPos + 1);
            std::string curVal = curAnt.substr(antPos, nextAnt - antPos);

            size_t eqPos = curVal.find("=");
            if(eqPos == std::string::npos)
                throw VoreenException(std::string("Annotation without assignment at '").append(curVal).append("'"));

            std::string tag = curVal.substr(1, eqPos - 1);

            size_t valPos = curVal.find_first_not_of(" ", eqPos + 1);
            if(valPos == std::string::npos)
                throw VoreenException(std::string("Assignment of empty value at '").append(curVal).append("'"));

            std::string value;
            if(curVal.at(valPos) == '"') {
                size_t closeQuotePos = curVal.find("\"", valPos + 1);
                if(closeQuotePos == std::string::npos)
                    throw VoreenException(std::string("No closing quotes at '").append(curVal).append("'"));
                else
                    value = curVal.substr(valPos + 1, closeQuotePos - valPos - 1);
            }
            else
                value = curVal.substr(valPos, curVal.find_first_of(" ", valPos + 1) - valPos);

            tag.erase(std::remove(tag.begin(), tag.end(), ' '), tag.end());
            //value.erase(std::remove(value.begin(), value.end(), ' '), value.end());

            if(std::find(tagVals_.begin(), tagVals_.end(), tag) != tagVals_.end())
                curArg->registerAnnotation(tag, value);
            else
                throw VoreenException(std::string("Unknown annotation tag '").append(tag).append("'"));

            antPos = nextAnt;
        }

        // if no name was specified via the @name tag, take the argument name from the code as the name
        if(curArgs_.at(i).values_["name"].empty())
            curArgs_.at(i).values_["name"] = curArgs_.at(i).arg_;
        if(curArgs_.at(i).values_["label"].empty())
            curArgs_.at(i).values_["label"] = curArgs_.at(i).values_["name"];
    }

}

void DynamicCLProcessor::setWorkDimensions() {

    // look for the argument which will tell the kernel its workload.  If no such argument was specified with the appropriate tag, use the first
    // image or volume, if present, and otherwise use the first int(2|3) argument.
    curWorkDimensionSource_ = 0;

    for(size_t i = 0; i < curArgs_.size(); i++) {
        if(!curWorkDimensionSource_ && curArgs_[i].values_["workdimssource"] == "true") {
            if(curArgs_[i].attType_ == IMAGE2D && curArgs_[i].att_ == INPORT) {
                LWARNING("Currently, image inports can not be used as the kernel working dimension source, please use outports instead.");
                break;
            }
            curWorkDimensionSource_ = &(curArgs_[i]);
            LINFO("Using argument '" << curArgs_[i].values_["name"] << "' as kernel dimensions source");
        }
        else if(curArgs_[i].values_["workdimssource"] == "true")
            LWARNING("Please specify only one argument as the kernel working dimension source.");
    }

    if(!curWorkDimensionSource_) {
        LINFO("Working dimensions source not specified.");
        for(size_t i = 0; i < curArgs_.size(); i++) {
            if((curArgs_[i].attType_ == IMAGE2D && curArgs_[i].att_ == OUTPORT) || curArgs_[i].attType_ == IMAGE3D) {
                LINFO("Using argument '" << curArgs_[i].values_["name"] << "' as source");
                curWorkDimensionSource_ = &(curArgs_[i]);
                break;
            }
        }
    }

    if(!curWorkDimensionSource_) {
        for(size_t i = 0; i < curArgs_.size(); i++) {
            if(curArgs_[i].attType_ == INTPROPERTY || curArgs_[i].attType_ == INT2PROPERTY || curArgs_[i].attType_ == INT3PROPERTY) {
                LINFO("Using argument '" << curArgs_[i].values_["name"] << "' as source");
                curWorkDimensionSource_ = &(curArgs_[i]);
                break;
            }
        }
    }
}

void DynamicCLProcessor::refreshPortsAndProperties() {

    for(size_t i = 0; i < curArgs_.size(); i++) {
        ArgInfo curArg = curArgs_.at(i);

        if(curArg.att_ == INPORT) {
            Port* p = getPort(curArg.values_["name"]);
            tgtAssert(!p || !p->isOutport(), "Port was expected to be inport");
            if(!p)
                addNewInport(curArg);
        } else if(curArg.att_ == OUTPORT) {
            Port* p = getPort(curArg.values_["name"]);
            tgtAssert(!p || p->isOutport(), "Port was expected to be outport");
            if(!p)
                addNewOutport(curArg);
        } else if(curArg.att_ == PROPERTY) {
            Property* p = getProperty(curArg.values_["name"]);
            if(p)
                updateProperty(p, curArg);
            else
                addProperty(generateNewProperty(curArg));
        }
    }

    // inform the observers that ports and properties might have changed
    notifyPortsChanged();
    notifyPropertiesChanged();
}

void DynamicCLProcessor::removeObsoletePortsAndProperties() {

    // Handle properties.
    std::vector<Property*> obsoleteProperties;

    for(Property* p : dynamicProperties_) {

        bool found = false;
        for(size_t j = 0; j < curArgs_.size(); j++) {
            if (curArgs_.at(j).att_ == PROPERTY && getProperty(curArgs_.at(j).values_["name"]) == p) {
                found = true;
                break;
            }
        }
        if(!found)
            obsoleteProperties.push_back(p);
    }

    if (!obsoleteProperties.empty()) {
        for (Property* p : obsoleteProperties) {
            removeProperty(p);
            dynamicProperties_.erase(std::find(dynamicProperties_.begin(), dynamicProperties_.end(), p));
            delete p;
        }

        notifyPropertiesChanged();
    }

    // Handle ports.
    std::vector<Port*> obsoletePorts;

    for(Port* p : dynamicPorts_) {
        bool found = false;
        for(size_t j = 0; j < curArgs_.size(); j++) {
            if(curArgs_.at(j).att_ != PROPERTY && getPort(curArgs_.at(j).values_["name"]) == p) {
                if ((p->isOutport() && curArgs_.at(j).att_ == OUTPORT) || (p->isInport() && curArgs_.at(j).att_ == INPORT)) {
                    found = true;
                    break;
                }
            }
        }
        if(!found)
            obsoletePorts.push_back(p);
    }

    if (!obsoletePorts.empty()) {

        for (size_t i = 0; i < obsoletePorts.size(); i++) {
            Port* p = obsoletePorts.at(i);
            deinitializePort(p);
            removePort(p);
            portArgMap_.erase(p);
            dynamicPorts_.erase(std::find(dynamicPorts_.begin(), dynamicPorts_.end(), p));
        }

        notifyPortsChanged();

        for (size_t i = 0; i < obsoletePorts.size(); i++)
            delete obsoletePorts.at(i);
    }
}

void DynamicCLProcessor::addNewInport(ArgInfo& arg) {
    if(arg.attType_ == IMAGE2D) {
        LINFO("Adding image inport '" + arg.values_["name"] + "'");
        RenderPort* imgInport = new RenderPort(Port::INPORT, arg.values_["name"], arg.values_["name"], false, arg.invLevel_, RenderPort::RENDERSIZE_ORIGIN, GL_RGBA16);
        dynamicPorts_.push_back(imgInport);
        addPort(imgInport);
        initializePort(imgInport);
        portArgMap_[imgInport] = &arg;
    } else if(arg.attType_ == IMAGE3D) {
        LINFO("Adding volume inport '" + arg.values_["name"] + "'");
        VolumePort* volInport = new VolumePort(Port::INPORT, arg.values_["name"], arg.values_["name"], false, arg.invLevel_);
        dynamicPorts_.push_back(volInport);
        addPort(volInport);
        initializePort(volInport);
        portArgMap_[volInport] = &arg;
    }
}

void DynamicCLProcessor::addNewOutport(ArgInfo& arg) {
    if(arg.attType_ == IMAGE2D) {
        LINFO("Adding image outport '" + arg.values_["name"] + "'");
        RenderPort* imgOutport = new RenderPort(Port::OUTPORT, arg.values_["name"], arg.values_["name"], true, arg.invLevel_, RenderPort::RENDERSIZE_RECEIVER, GL_RGBA16);
        dynamicPorts_.push_back(imgOutport);
        addPort(imgOutport);
        initializePort(imgOutport);
        imgOutport->onSizeReceiveChange<DynamicCLProcessor>(this, &DynamicCLProcessor::portResized);
        portArgMap_[imgOutport] = &arg;
    } else if(arg.attType_ == IMAGE3D) {
        LINFO("Adding volume outport '" + arg.values_["name"] + "'");
        VolumePort* volOutport = new VolumePort(Port::OUTPORT, arg.values_["name"], arg.values_["name"], true, arg.invLevel_);
        dynamicPorts_.push_back(volOutport);
        addPort(volOutport);
        initializePort(volOutport);
        portArgMap_[volOutport] = &arg;
    }
}

std::map<std::string, float> DynamicCLProcessor::getNumericPropertyValues(ArgInfo& arg) {

    std::map<std::string, float> rangeValues;

    for(size_t i = 0; i < numericTagVals_.size(); i++) {
        // value was not specified: default value will be used
        if(arg.values_.find(numericTagVals_.at(i)) == arg.values_.end())
            continue;

        std::stringstream s;
        float val;
        s << arg.values_[numericTagVals_.at(i)];
        s >> val;

        if(s.rdstate() & std::ifstream::failbit)
            throw VoreenException(std::string("Could not read numeric annotation value '").append(numericTagVals_.at(i)).append("' for argument '").append(arg.values_["name"]).append("'"));
        else
            rangeValues[numericTagVals_.at(i)] = val;
    }

    return rangeValues;
}

Property* DynamicCLProcessor::generateNewProperty(ArgInfo& arg) {

    LINFO("Adding property '" + arg.values_["name"] + "'");
    Property* p = 0;

    if(   arg.attType_ == FLOATPROPERTY || arg.attType_ == FLOAT2PROPERTY || arg.attType_ == FLOAT3PROPERTY || arg.attType_ == FLOAT4PROPERTY
       || arg.attType_ == INTPROPERTY   || arg.attType_ == INT2PROPERTY   || arg.attType_ == INT3PROPERTY   || arg.attType_ == INT4PROPERTY
       || arg.attType_ == MAT2PROPERTY  || arg.attType_ == MAT3PROPERTY   || arg.attType_ == MAT4PROPERTY) {

        std::map<std::string, float> rangeValues = getNumericPropertyValues(arg);

        // fill in default values if they were not specified
        if (rangeValues.find("min") == rangeValues.end()) rangeValues["min"] = 0.f;
        if (rangeValues.find("max") == rangeValues.end()) rangeValues["max"] = rangeValues["min"] + 1.f;
        if (rangeValues.find("value") == rangeValues.end()) rangeValues["value"] = rangeValues["min"];

        switch(arg.attType_) {
            case FLOATPROPERTY:  p =  new FloatProperty(arg.values_["name"], arg.values_["label"], rangeValues["value"], rangeValues["min"], rangeValues["max"], arg.invLevel_); break;
            case FLOAT2PROPERTY: p = new FloatVec2Property(arg.values_["name"], arg.values_["label"], tgt::vec2(rangeValues["value"]), tgt::vec2(rangeValues["min"]), tgt::vec2(rangeValues["max"]), arg.invLevel_); break;
            case FLOAT3PROPERTY: p = new FloatVec3Property(arg.values_["name"], arg.values_["label"], tgt::vec3(rangeValues["value"]), tgt::vec3(rangeValues["min"]), tgt::vec3(rangeValues["max"]), arg.invLevel_); break;
            case FLOAT4PROPERTY: p = new FloatVec4Property(arg.values_["name"], arg.values_["label"], tgt::vec4(rangeValues["value"]), tgt::vec4(rangeValues["min"]), tgt::vec4(rangeValues["max"]), arg.invLevel_); break;

            case INTPROPERTY:  p =  new IntProperty(arg.values_["name"], arg.values_["label"], tgt::iround(rangeValues["value"]), tgt::iround(rangeValues["min"]), tgt::iround(rangeValues["max"]), arg.invLevel_); break;
            case INT2PROPERTY: p = new IntVec2Property(arg.values_["name"], arg.values_["label"], tgt::ivec2(tgt::iround(rangeValues["value"])),
                                                       tgt::ivec2(tgt::iround(rangeValues["min"])), tgt::ivec2(tgt::iround(rangeValues["max"])), arg.invLevel_);
                               break;
            case INT3PROPERTY: p = new IntVec3Property(arg.values_["name"], arg.values_["label"], tgt::ivec3(tgt::iround(rangeValues["value"])),
                                                       tgt::ivec3(tgt::iround(rangeValues["min"])), tgt::ivec3(tgt::iround(rangeValues["max"])), arg.invLevel_);
                               break;
            case INT4PROPERTY: p = new IntVec4Property(arg.values_["name"], arg.values_["label"], tgt::ivec4(tgt::iround(rangeValues["value"])),
                                                       tgt::ivec4(tgt::iround(rangeValues["min"])), tgt::ivec4(tgt::iround(rangeValues["max"])), arg.invLevel_);
                               break;
            case MAT2PROPERTY: p = new FloatMat2Property(arg.values_["name"], arg.values_["label"], tgt::mat2(rangeValues["value"]), tgt::mat2(rangeValues["min"]), tgt::mat2(rangeValues["max"]), arg.invLevel_); break;
            case MAT3PROPERTY: p = new FloatMat3Property(arg.values_["name"], arg.values_["label"], tgt::mat3(rangeValues["value"]), tgt::mat3(rangeValues["min"]), tgt::mat3(rangeValues["max"]), arg.invLevel_); break;
            case MAT4PROPERTY: p = new FloatMat4Property(arg.values_["name"], arg.values_["label"], tgt::mat4(rangeValues["value"]), tgt::mat4(rangeValues["min"]), tgt::mat4(rangeValues["max"]), arg.invLevel_); break;

            default: break;
        }
    }

    //if(arg.attType_ == BOOLPROPERTY) {
        //bool val = true;
        //if(arg.values_.find("value") == arg.values_.end() && arg.values_["value"] == "false")
            //val = false;
        //p = new BoolProperty(arg.values_["name"], arg.values_["label"], val);
    //}

    if(arg.attType_ == TRANSFUNCPROPERTY) {
        p = new TransFunc1DKeysProperty(arg.values_["name"], arg.values_["label"], arg.invLevel_);
        ((TransFunc1DKeysProperty*)p)->set1DKeys(new TransFunc1DKeys());
    }

    // Keep track about generated properties.
    if (p) {
        // If the property has been added after the initialization of the processor, it needs to be initialized explicitly.
        if (isInitialized())
            p->initialize();
        dynamicProperties_.push_back(p);
    }

    return p;
}

void DynamicCLProcessor::updateProperty(Property* p, ArgInfo& arg) {

    LINFO("Updating property '" + arg.values_["name"] + "'");
    bool handleAsNew = false;

    p->setInvalidationLevel(arg.invLevel_);

    if(   arg.attType_ == FLOATPROPERTY || arg.attType_ == FLOAT2PROPERTY || arg.attType_ == FLOAT3PROPERTY || arg.attType_ == FLOAT4PROPERTY
       || arg.attType_ == INTPROPERTY   || arg.attType_ == INT2PROPERTY   || arg.attType_ == INT3PROPERTY   || arg.attType_ == INT4PROPERTY
       || arg.attType_ == MAT2PROPERTY  || arg.attType_ == MAT3PROPERTY   || arg.attType_ == MAT4PROPERTY) {

        std::map<std::string, float> rangeValues = getNumericPropertyValues(arg);

        switch(arg.attType_) {
            case FLOATPROPERTY:
                if(FloatProperty* pp = dynamic_cast<FloatProperty*>(p)) {
                    if(rangeValues.find("value") != rangeValues.end()) pp->set(rangeValues["value"]);
                    if(rangeValues.find("min")   != rangeValues.end()) pp->setMinValue(rangeValues["min"]);
                    if(rangeValues.find("max")   != rangeValues.end()) pp->setMaxValue(rangeValues["max"]);
                } else // previously used name is now attached to property of different kind
                    handleAsNew = true;
                break;
            case FLOAT2PROPERTY:
                if(FloatVec2Property* pp = dynamic_cast<FloatVec2Property*>(p)) {
                    if(rangeValues.find("value") != rangeValues.end()) pp->set(tgt::vec2(rangeValues["value"]));
                    if(rangeValues.find("min")   != rangeValues.end()) pp->setMinValue(tgt::vec2(rangeValues["min"]));
                    if(rangeValues.find("max")   != rangeValues.end()) pp->setMaxValue(tgt::vec2(rangeValues["max"]));
                } else // previously used name is now attached to property of different kind
                    handleAsNew = true;
                break;
            case FLOAT3PROPERTY:
                if(FloatVec3Property* pp = dynamic_cast<FloatVec3Property*>(p)) {
                    if(rangeValues.find("value") != rangeValues.end()) pp->set(tgt::vec3(rangeValues["value"]));
                    if(rangeValues.find("min")   != rangeValues.end()) pp->setMinValue(tgt::vec3(rangeValues["min"]));
                    if(rangeValues.find("max")   != rangeValues.end()) pp->setMaxValue(tgt::vec3(rangeValues["max"]));
                } else // previously used name is now attached to property of different kind
                    handleAsNew = true;
                break;
            case FLOAT4PROPERTY:
                if(FloatVec4Property* pp = dynamic_cast<FloatVec4Property*>(p)) {
                    if(rangeValues.find("value") != rangeValues.end()) pp->set(tgt::vec4(rangeValues["value"]));
                    if(rangeValues.find("min")   != rangeValues.end()) pp->setMinValue(tgt::vec4(rangeValues["min"]));
                    if(rangeValues.find("max")   != rangeValues.end()) pp->setMaxValue(tgt::vec4(rangeValues["max"]));
                } else // previously used name is now attached to property of different kind
                    handleAsNew = true;
                break;

            case INTPROPERTY:
                if(IntProperty* pp = dynamic_cast<IntProperty*>(p)) {
                    if (rangeValues.find("value") != rangeValues.end()) pp->set(tgt::iround(rangeValues["value"]));
                    if (rangeValues.find("min")   != rangeValues.end()) pp->setMinValue(tgt::iround(rangeValues["min"]));
                    if (rangeValues.find("max")   != rangeValues.end()) pp->setMaxValue(tgt::iround(rangeValues["max"]));
                } else // previously used name is now attached to property of different kind
                    handleAsNew = true;
                break;
            case INT2PROPERTY:
                if(IntVec2Property* pp = dynamic_cast<IntVec2Property*>(p)) {
                    if(rangeValues.find("value") != rangeValues.end()) pp->set(tgt::ivec2(tgt::iround(rangeValues["value"])));
                    if(rangeValues.find("min")   != rangeValues.end()) pp->setMinValue(tgt::ivec2(tgt::iround(rangeValues["min"])));
                    if(rangeValues.find("max")   != rangeValues.end()) pp->setMaxValue(tgt::ivec2(tgt::iround(rangeValues["max"])));
                } else // previously used name is now attached to property of different kind
                    handleAsNew = true;
                break;
            case INT3PROPERTY:
                if(IntVec3Property* pp = dynamic_cast<IntVec3Property*>(p)) {
                    if(rangeValues.find("value") != rangeValues.end()) pp->set(tgt::ivec3(tgt::iround(rangeValues["value"])));
                    if(rangeValues.find("min")   != rangeValues.end()) pp->setMinValue(tgt::ivec3(tgt::iround(rangeValues["min"])));
                    if(rangeValues.find("max")   != rangeValues.end()) pp->setMaxValue(tgt::ivec3(tgt::iround(rangeValues["max"])));
                } else // previously used name is now attached to property of different kind
                    handleAsNew = true;
                break;
            case INT4PROPERTY:
                if(IntVec4Property* pp = dynamic_cast<IntVec4Property*>(p)) {
                    if(rangeValues.find("value") != rangeValues.end()) pp->set(tgt::ivec4(tgt::iround(rangeValues["value"])));
                    if(rangeValues.find("min")   != rangeValues.end()) pp->setMinValue(tgt::ivec4(tgt::iround(rangeValues["min"])));
                    if(rangeValues.find("max")   != rangeValues.end()) pp->setMaxValue(tgt::ivec4(tgt::iround(rangeValues["max"])));
                } else // previously used name is now attached to property of different kind
                    handleAsNew = true;
                break;

            case MAT2PROPERTY:
                if(FloatMat2Property* pp = dynamic_cast<FloatMat2Property*>(p)) {
                    if(rangeValues.find("value") != rangeValues.end()) pp->set(tgt::mat2(rangeValues["value"]));
                    if(rangeValues.find("min")   != rangeValues.end()) pp->setMinValue(tgt::mat2(rangeValues["min"]));
                    if(rangeValues.find("max")   != rangeValues.end()) pp->setMaxValue(tgt::mat2(rangeValues["max"]));
                } else // previously used name is now attached to property of different kind
                    handleAsNew = true;
                break;
            case MAT3PROPERTY:
                if(FloatMat3Property* pp = dynamic_cast<FloatMat3Property*>(p)) {
                    if(rangeValues.find("value") != rangeValues.end()) pp->set(tgt::mat3(rangeValues["value"]));
                    if(rangeValues.find("min")   != rangeValues.end()) pp->setMinValue(tgt::mat3(rangeValues["min"]));
                    if(rangeValues.find("max")   != rangeValues.end()) pp->setMaxValue(tgt::mat3(rangeValues["max"]));
                } else // previously used name is now attached to property of different kind
                    handleAsNew = true;
                break;
            case MAT4PROPERTY:
                if(FloatMat4Property* pp = dynamic_cast<FloatMat4Property*>(p)) {
                    if(rangeValues.find("value") != rangeValues.end()) pp->set(tgt::mat4(rangeValues["value"]));
                    if(rangeValues.find("min")   != rangeValues.end()) pp->setMinValue(tgt::mat4(rangeValues["min"]));
                    if(rangeValues.find("max")   != rangeValues.end()) pp->setMaxValue(tgt::mat4(rangeValues["max"]));
                } else // previously used name is now attached to property of different kind
                    handleAsNew = true;
                break;

            default: break;
        }
    }

    //if(arg.attType_ == BOOLPROPERTY) {
        //if(BoolProperty* pp = dynamic_cast<BoolProperty*>(p)) {
            //if(arg.values_["value"] == "false")
                //pp->set(false);
            //else
                //pp->set(true);
        //} else // previously used name is now attached to property of different kind
            //handleAsNew = true;
    //}
    //
    if(arg.attType_ == TRANSFUNCPROPERTY) {
        if(dynamic_cast<TransFunc1DKeysProperty*>(p))
            return;
        else
            handleAsNew = true;
    }

    if(handleAsNew) {
        removeProperty(p);
        dynamicProperties_.erase(std::find(dynamicProperties_.begin(), dynamicProperties_.end(), p));
        delete p;
        addProperty(generateNewProperty(arg));
    }
}

void DynamicCLProcessor::setInvalidationLevels() {

    for(size_t i = 0; i < curArgs_.size(); i++) {
        if (curArgs_.at(i).values_.find("invalidationlevel") == curArgs_.at(i).values_.end()) {
            curArgs_.at(i).invLevel_ = Processor::INVALID_RESULT;
            continue;

        } else {
            std::string value = curArgs_.at(i).values_["invalidationlevel"];

            if(value == "valid")
                curArgs_.at(i).invLevel_ = Processor::VALID;
            else if(value == "invalid_result")
                curArgs_.at(i).invLevel_ = Processor::INVALID_RESULT;
            else if(value == "invalid_parameters")
                curArgs_.at(i).invLevel_ = Processor::INVALID_PARAMETERS;
            else if(value == "invalid_program")
                curArgs_.at(i).invLevel_ = Processor::INVALID_PROGRAM;
            else if(value == "invalid_ports")
                curArgs_.at(i).invLevel_ = Processor::INVALID_PORTS;
            else if(value == "invalid_processor")
                curArgs_.at(i).invLevel_ = Processor::INVALID_PROCESSOR;
            else {
                LERROR("Unknown invalidation level specified in argument '" << curArgs_.at(i).values_["name"] << "', defaulting to INVALID_RESULT.");
                curArgs_.at(i).invLevel_ = Processor::INVALID_RESULT;
            }
        }
    }
}

void DynamicCLProcessor::ArgInfo::registerAttribute(const std::string& mod, const std::string& type) {
    if(type == "") {
        LWARNING("Kernel argument " << arg_ << ": unsupported argument type " << type);
        return;
    }
    else if(type == "image2d_t") {
        if(mod == "write_only")
            att_ = OUTPORT;
        else
            att_ = INPORT;
        attType_ = IMAGE2D;
    } else if(type == "image3d_t") {
        if(mod == "write_only") {
            LWARNING("Currently, only read_only 3d images are supported. Please use mod_volumewrite.cl and the volume_t type instead.");
            return;
        }
        att_ = INPORT;
        attType_ = IMAGE3D;
    } else if(type == "volume_t*" || type == "volume_t *") {
        if(mod == "read_only") {
            LWARNING("Please use image3d_t instead of volume_t for read-only arguments.");
            return;
        }
        att_ = OUTPORT;
        attType_ = IMAGE3D;
    } else if(type == "void*" || type == "void *") {
        att_ = UNKNOWN;
        attType_ = UNSUPPORTED;
    } else { // other datatypes that can be mapped to a property
        att_ = PROPERTY;
        attType_ = propertyMap_[type];
    }
}

void DynamicCLProcessor::ArgInfo::registerAnnotation(const std::string& tag, const std::string& value) {
    if(value == "") {
        LWARNING("Empty value assigned for tag '" << tag << "'");
        return;
    }

    // if the transferfunction tag was used, change the attribute from port to property
    if(tag == "transferfunction" && value == "true" && attType_ == IMAGE2D) {
        if(att_ == OUTPORT) {
            LERROR("Write-only transfer function not supported in argument '" << values_["name"] << "'.");
            return;
        }
        att_ = PROPERTY;
        attType_ = TRANSFUNCPROPERTY;
    }

    values_[tag] = value;
}

void DynamicCLProcessor::setupKeywords() {

    // the following commented types are reserved, but not yet available
    //clTypesIgn_.push_back("float2x3"); clTypesIgn_.push_back("float2x4");
    //clTypesIgn_.push_back("float2x8"); clTypesIgn_.push_back("float2x16"); clTypesIgn_.push_back("float3x2");
    //clTypesIgn_.push_back("float3x4"); clTypesIgn_.push_back("float3x8"); clTypesIgn_.push_back("float3x16"); clTypesIgn_.push_back("float4x2");
    //clTypesIgn_.push_back("float4x3"); clTypesIgn_.push_back("float4x8"); clTypesIgn_.push_back("float4x16");
    //clTypesIgn_.push_back("float8x2"); clTypesIgn_.push_back("float8x3"); clTypesIgn_.push_back("float8x4"); clTypesIgn_.push_back("float8x8");
    //clTypesIgn_.push_back("float8x16"); clTypesIgn_.push_back("float16x2"); clTypesIgn_.push_back("float16x3"); clTypesIgn_.push_back("float16x4");
    //clTypesIgn_.push_back("float16x8"); clTypesIgn_.push_back("float16x16"); clTypesIgn_.push_back("double2x2"); clTypesIgn_.push_back("double2x3");
    //clTypesIgn_.push_back("double2x4"); clTypesIgn_.push_back("double2x8"); clTypesIgn_.push_back("double2x16"); clTypesIgn_.push_back("double3x2");
    //clTypesIgn_.push_back("double3x3"); clTypesIgn_.push_back("double3x4"); clTypesIgn_.push_back("double3x8"); clTypesIgn_.push_back("double3x16");
    //clTypesIgn_.push_back("double4x2"); clTypesIgn_.push_back("double4x3"); clTypesIgn_.push_back("double4x4"); clTypesIgn_.push_back("double4x8");
    //clTypesIgn_.push_back("double4x16"); clTypesIgn_.push_back("double8x2"); clTypesIgn_.push_back("double8x3"); clTypesIgn_.push_back("double8x4");
    //clTypesIgn_.push_back("double8x8"); clTypesIgn_.push_back("double8x16"); clTypesIgn_.push_back("double16x2"); clTypesIgn_.push_back("double16x3");
    //clTypesIgn_.push_back("double16x4"); clTypesIgn_.push_back("double16x8"); clTypesIgn_.push_back("double16x16");
    clTypesIgn_.push_back("unsigned char"); clTypesIgn_.push_back("unsigned short"); clTypesIgn_.push_back("unsigned int"); clTypesIgn_.push_back("unsigned long");
    clTypesIgn_.push_back("bool2"); clTypesIgn_.push_back("bool3"); clTypesIgn_.push_back("bool4"); clTypesIgn_.push_back("bool8"); clTypesIgn_.push_back("bool16");
    clTypesIgn_.push_back("quad2"); clTypesIgn_.push_back("quad3"); clTypesIgn_.push_back("quad4"); clTypesIgn_.push_back("quad8"); clTypesIgn_.push_back("quad16");
    clTypesIgn_.push_back("char4"); clTypesIgn_.push_back("char8"); clTypesIgn_.push_back("char16"); clTypesIgn_.push_back("uchar2");
    clTypesIgn_.push_back("uchar3"); clTypesIgn_.push_back("uchar4"); clTypesIgn_.push_back("uchar8"); clTypesIgn_.push_back("uchar16");
    clTypesIgn_.push_back("short2"); clTypesIgn_.push_back("short3"); clTypesIgn_.push_back("short4"); clTypesIgn_.push_back("short8");
    clTypesIgn_.push_back("short16"); clTypesIgn_.push_back("ushort2"); clTypesIgn_.push_back("ushort3"); clTypesIgn_.push_back("ushort4");
    clTypesIgn_.push_back("ushort8"); clTypesIgn_.push_back("ushort16"); clTypesIgn_.push_back("int8"); clTypesIgn_.push_back("int16"); clTypesIgn_.push_back("uint2");
    clTypesIgn_.push_back("uint3"); clTypesIgn_.push_back("uint4"); clTypesIgn_.push_back("uint8"); clTypesIgn_.push_back("uint16");
    clTypesIgn_.push_back("long2"); clTypesIgn_.push_back("long3"); clTypesIgn_.push_back("long4"); clTypesIgn_.push_back("long8");
    clTypesIgn_.push_back("long16"); clTypesIgn_.push_back("ulong2"); clTypesIgn_.push_back("ulong3"); clTypesIgn_.push_back("ulong4");
    clTypesIgn_.push_back("ulong8"); clTypesIgn_.push_back("ulong16"); clTypesIgn_.push_back("float8");
    clTypesIgn_.push_back("double2"); clTypesIgn_.push_back("double3"); clTypesIgn_.push_back("double4"); clTypesIgn_.push_back("double8");
    clTypesIgn_.push_back("double16"); clTypesIgn_.push_back("half2"); clTypesIgn_.push_back("half3"); clTypesIgn_.push_back("half4");
    clTypesIgn_.push_back("char2"); clTypesIgn_.push_back("char3"); clTypesIgn_.push_back("half8"); clTypesIgn_.push_back("half16");
    clTypesIgn_.push_back("uchar"); clTypesIgn_.push_back("ushort"); clTypesIgn_.push_back("uint"); clTypesIgn_.push_back("ulong");
    clTypesIgn_.push_back("half"); clTypesIgn_.push_back("double"); clTypesIgn_.push_back("quad");clTypesIgn_.push_back("ptrdiff_t");
    clTypesIgn_.push_back("intptr_t"); clTypesIgn_.push_back("uintptr_t"); clTypesIgn_.push_back("sampler_t"); clTypesIgn_.push_back("event_t");

    // the following commented types are reserved, but not yet available
    //clTypes_.push_back("float2x2");
    //clTypes_.push_back("float3x3");
    //clTypes_.push_back("float4x4");
    clTypes_.push_back("int2");
    clTypes_.push_back("int3");
    clTypes_.push_back("int4");
    clTypes_.push_back("float2");
    clTypes_.push_back("float3");   // only available in OpenCL 1.1
    clTypes_.push_back("float4");
    clTypes_.push_back("float16");  // 4x4 matrix
    clTypes_.push_back("bool");
    clTypes_.push_back("int");
    clTypes_.push_back("float");
    clTypes_.push_back("image2d_t");
    clTypes_.push_back("image3d_t");
    clTypes_.push_back("volume_t*");
    clTypes_.push_back("volume_t *");
    clTypes_.push_back("void*");
    clTypes_.push_back("void *");

    //propertyMap_["float2x2"] = MAT2PROPERTY;
    //propertyMap_["float3x3"] = MAT3PROPERTY;
    //propertyMap_["float4x4"] = MAT4PROPERTY;
    propertyMap_["int2"] = INT2PROPERTY;
    propertyMap_["int3"] = INT3PROPERTY;
    propertyMap_["int4"] = INT4PROPERTY;
    propertyMap_["float2"] = FLOAT2PROPERTY;
    propertyMap_["float3"] = FLOAT3PROPERTY;
    propertyMap_["float4"] = FLOAT4PROPERTY;
    propertyMap_["float16"] = MAT4PROPERTY;
    propertyMap_["bool"] = BOOLPROPERTY;
    propertyMap_["int"] = INTPROPERTY;
    propertyMap_["float"] = FLOATPROPERTY;
    propertyMap_[""] = UNSUPPORTED;

    storeModsIgn_.push_back("__global");
    storeModsIgn_.push_back("global");
    storeModsIgn_.push_back("__local");
    storeModsIgn_.push_back("local");
    storeModsIgn_.push_back("__constant");
    storeModsIgn_.push_back("constant");
    storeModsIgn_.push_back("__private");
    storeModsIgn_.push_back("private");

    storeMods_.push_back("__read_only");
    storeMods_.push_back("read_only");
    storeMods_.push_back("__write_only");
    storeMods_.push_back("write_only");

    tagVals_.push_back("min");
    tagVals_.push_back("max");
    tagVals_.push_back("value");

    numericTagVals_ = tagVals_;

    tagVals_.push_back("name");
    tagVals_.push_back("label");
    tagVals_.push_back("transferfunction");
    tagVals_.push_back("workdimssource");
    tagVals_.push_back("invalidationlevel");
    tagVals_.push_back("volumedimensions");
    tagVals_.push_back("volumetype");
}

} // namespace voreen
