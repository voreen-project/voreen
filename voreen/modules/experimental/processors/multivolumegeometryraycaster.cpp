/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
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

#include "multivolumegeometryraycaster.h"

#include "tgt/textureunit.h"
#include "tgt/glmath.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "voreen/core/utils/glsl.h"
#include "voreen/core/utils/classificationmodes.h"
#include "voreen/core/utils/voreenqualitymode.h"
#include "voreen/core/datastructures/geometry/geometrysequence.h"

#include <sstream>

#ifdef GL_ATOMIC_COUNTER_BUFFER //disable compilation for old gl headers

#ifndef GL_SHADER_STORAGE_BUFFER
#define GL_SHADER_STORAGE_BUFFER          0x90D2
#endif

#ifndef GL_SHADER_STORAGE_BARRIER_BIT
#define GL_SHADER_STORAGE_BARRIER_BIT     0x2000
#endif

using tgt::vec3;
using tgt::TextureUnit;

namespace voreen {

MultiVolumeGeometryRaycaster::MultiVolumeGeometryRaycaster()
    : VolumeRaycaster(true, true)
    , volumeInport2_(Port::INPORT, "volume2", "Volume2 Input", false, Processor::INVALID_PROGRAM)
    , volumeInport3_(Port::INPORT, "volume3", "Volume3 Input", false, Processor::INVALID_PROGRAM)
    , volumeInport4_(Port::INPORT, "volume4", "Volume4 Input", false, Processor::INVALID_PROGRAM)
    , pgPort1_(Port::INPORT, "pg1", "Proxy-Geometry 1")
    , pgPort2_(Port::INPORT, "pg2", "Proxy-Geometry 2")
    , pgPort3_(Port::INPORT, "pg3", "Proxy-Geometry 3")
    , pgPort4_(Port::INPORT, "pg4", "Proxy-Geometry 4")
    , geometryPort_(Port::INPORT, "geometry", "Geometry", true)
    //, outport1_(Port::OUTPORT, "output2", "Image Output 2", true, Processor::INVALID_PROGRAM, RenderPort::RENDERSIZE_DEFAULT, GL_RGBA16F_ARB)
    //, outport2_(Port::OUTPORT, "output3", "Image Output 3", true, Processor::INVALID_PROGRAM, RenderPort::RENDERSIZE_DEFAULT, GL_RGBA16F_ARB)
    , tmpPort_(Port::OUTPORT, "tmp", "Temp", true, Processor::INVALID_PROGRAM, RenderPort::RENDERSIZE_DEFAULT, GL_RGBA16F)
    , headPort_(Port::OUTPORT, "head", "Head", true, Processor::INVALID_PROGRAM, RenderPort::RENDERSIZE_DEFAULT, GL_R32UI)
    , rcShaderProp_("raycast_shader", "Raycasting Shader", "oit_rc.frag", "passthrough.vert")
    , oitShaderProp_("oit_shader", "OIT Shader", "oit.frag", "oit.vert")
    , oitClearShaderProp_("oitClear_shader", "OIT Clear Shader", "oit_clear.frag", "oit_clear.vert")
    , oitSortShaderProp_("oitSort_shader", "OIT Sort Shader", "oit_sort.frag", "oit_clear.vert")
    , sortingAlgorithm_("sortingAlgo", "Sorting Algorithm", Processor::INVALID_PROGRAM)
    , numPages_("numPages", "Num Pages", 4, 1, 32)
    , maxDepthComplexity_("maxDepthComplexity", "Max Depth Complexity", 4, 2, 32, Processor::INVALID_PROGRAM)
    , benchmark_("benchmark", "Benchmark", false)
    , adjustSize_("adjustSize", "Adjust Size", false)
    , storageBufferSizeProp_("storageBufferSize", "Storage Buffer Size", 7000000, 1000000, 50000000)
    , classificationMode2_("classification2", "Classification 2", Processor::INVALID_PROGRAM)
    , classificationMode3_("classification3", "Classification 3", Processor::INVALID_PROGRAM)
    , classificationMode4_("classification4", "Classification 4", Processor::INVALID_PROGRAM)
    , atomicBuffer_(0)
    , atomicBufferSize_(0)
    , atomicBufferInitialized_(false)
    , storageBuffer_(0)
    , storageBufferSize_(0)
    , targetStorageBufferSize_(7000000)
    , storageBufferInitialized_(false)
    , proxyGeometryBuffer_(0)
    , shadeMode1_("shading1", "Shading 1", Processor::INVALID_PROGRAM)
    , shadeMode2_("shading2", "Shading 2", Processor::INVALID_PROGRAM)
    , shadeMode3_("shading3", "Shading 3", Processor::INVALID_PROGRAM)
    , shadeMode4_("shading4", "Shading 4", Processor::INVALID_PROGRAM)
    , clippingGradientDepth_("clippingGradientWidth", "Cliiping Gradient Depth", 0.0001, 0.0f, 1.0f)
    , transferFunc1_("transferFunction1", "Transfer Function 1")
    , transferFunc2_("transferFunction2", "Transfer Function 2")
    , transferFunc3_("transferFunction3", "Transfer Function 3")
    , transferFunc4_("transferFunction4", "Transfer Function 4")
    , camera_("camera", "Camera", tgt::Camera(vec3(0.f, 0.f, 3.5f), vec3(0.f, 0.f, 0.f), vec3(0.f, 1.f, 0.f)))
{
    // ports
    volumeInport2_.addCondition(new PortConditionVolumeTypeGL());
    volumeInport3_.addCondition(new PortConditionVolumeTypeGL());
    volumeInport4_.addCondition(new PortConditionVolumeTypeGL());
    volumeInport_.showTextureAccessProperties(true);
    volumeInport2_.showTextureAccessProperties(true);
    volumeInport3_.showTextureAccessProperties(true);
    volumeInport4_.showTextureAccessProperties(true);
    addPort(pgPort1_);
    addPort(volumeInport2_);
    addPort(pgPort2_);
    addPort(volumeInport3_);
    addPort(pgPort3_);
    addPort(volumeInport4_);
    addPort(pgPort4_);

    addPort(geometryPort_);

    //addPort(outport1_);
    //addPort(outport2_);
    tmpPort_.setDeinitializeOnDisconnect(false);
    addPort(tmpPort_);
    headPort_.setDeinitializeOnDisconnect(false);
    addPort(headPort_);

    // shader property
    addProperty(oitClearShaderProp_);
    addProperty(oitShaderProp_);
    addProperty(oitSortShaderProp_);
    sortingAlgorithm_.addOption("BUBBLE_SORT", "Bubble Sort");
    sortingAlgorithm_.addOption("SELECTION_SORT", "Selection Sort");
    sortingAlgorithm_.addOption("INSERTION_SORT", "Insertion Sort");
    addProperty(sortingAlgorithm_);
    addProperty(rcShaderProp_);
    addProperty(numPages_);
    addProperty(maxDepthComplexity_);
    addProperty(benchmark_);
    //addProperty(adjustSize_);
    addProperty(storageBufferSizeProp_);

    addProperty(classificationMode_);
    ClassificationModes::fillProperty(&classificationMode2_);
    addProperty(classificationMode2_);
    ClassificationModes::fillProperty(&classificationMode3_);
    addProperty(classificationMode3_);
    ClassificationModes::fillProperty(&classificationMode4_);
    addProperty(classificationMode4_);

    // tf properties
    addProperty(transferFunc1_);
    addProperty(transferFunc2_);
    addProperty(transferFunc3_);
    addProperty(transferFunc4_);

    cameraHandler_ = new CameraInteractionHandler("cameraHandler", "Camera", &camera_);
    cameraHandler_->setEnabled(false); // not needed most of the time
    addInteractionHandler(cameraHandler_);
    addProperty(camera_);

    // shading properties
    addProperty(gradientMode_);
    GLSL::fillShadingModesProperty(shadeMode1_);
    addProperty(shadeMode1_);
    GLSL::fillShadingModesProperty(shadeMode2_);
    addProperty(shadeMode2_);
    GLSL::fillShadingModesProperty(shadeMode3_);
    addProperty(shadeMode3_);
    GLSL::fillShadingModesProperty(shadeMode4_);
    addProperty(shadeMode4_);

    addProperty(clippingGradientDepth_);

    addProperty(isoValue_);

    // lighting properties
    addProperty(lightPosition_);
    addProperty(lightAmbient_);
    addProperty(lightDiffuse_);
    addProperty(lightSpecular_);
    addProperty(materialShininess_);
    addProperty(applyLightAttenuation_);
    addProperty(lightAttenuation_);

    // assign lighting properties to property group
    lightPosition_.setGroupID("lighting");
    lightAmbient_.setGroupID("lighting");
    lightDiffuse_.setGroupID("lighting");
    lightSpecular_.setGroupID("lighting");
    materialShininess_.setGroupID("lighting");
    applyLightAttenuation_.setGroupID("lighting");
    lightAttenuation_.setGroupID("lighting");
    setPropertyGroupGuiName("lighting", "Lighting Parameters");

    // listen to changes of properties that influence the GUI state (i.e. visibility of other props)
    classificationMode_.onChange(MemberFunctionCallback<MultiVolumeGeometryRaycaster>(this, &MultiVolumeGeometryRaycaster::adjustPropertyVisibilities));
    shadeMode_.onChange(MemberFunctionCallback<MultiVolumeGeometryRaycaster>(this, &MultiVolumeGeometryRaycaster::adjustPropertyVisibilities));
    compositingMode1_.onChange(MemberFunctionCallback<MultiVolumeGeometryRaycaster>(this, &MultiVolumeGeometryRaycaster::adjustPropertyVisibilities));
    compositingMode2_.onChange(MemberFunctionCallback<MultiVolumeGeometryRaycaster>(this, &MultiVolumeGeometryRaycaster::adjustPropertyVisibilities));
    compositingMode3_.onChange(MemberFunctionCallback<MultiVolumeGeometryRaycaster>(this, &MultiVolumeGeometryRaycaster::adjustPropertyVisibilities));
    applyLightAttenuation_.onChange(MemberFunctionCallback<MultiVolumeGeometryRaycaster>(this, &MultiVolumeGeometryRaycaster::adjustPropertyVisibilities));
}

MultiVolumeGeometryRaycaster::~MultiVolumeGeometryRaycaster() {
    delete cameraHandler_;
}

Processor* MultiVolumeGeometryRaycaster::create() const {
    return new MultiVolumeGeometryRaycaster();
}

void MultiVolumeGeometryRaycaster::initializeAtomicBuffer(int numAtomics) {
    if(atomicBufferInitialized_ && (numAtomics == atomicBufferSize_))
        return;

    if(atomicBufferInitialized_) {
        glDeleteBuffers(1, &atomicBuffer_);
        atomicBufferInitialized_ = false;
    }

    // declare and generate a buffer object name
    glGenBuffers(1, &atomicBuffer_);
    // bind the buffer and define its initial storage capacity
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, atomicBuffer_);
    glBufferData(GL_ATOMIC_COUNTER_BUFFER, sizeof(GLuint) * numAtomics, NULL, GL_STATIC_DRAW);
    // unbind the buffer
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, 0);

    atomicBufferInitialized_ = true;
    atomicBufferSize_ = numAtomics;

    invalidate(INVALID_PROGRAM);
}

void MultiVolumeGeometryRaycaster::initializeProxyGeometryBuffer(size_t s) {
    //if(storageBufferInitialized_ && (s == storageBufferSize_))
        //return;

    //if(storageBufferInitialized_) {
        //glDeleteBuffers(1, &proxyGeometryBuffer_);
        //storageBufferInitialized_ = false;
    //}

    glGenBuffers(1, &proxyGeometryBuffer_);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, proxyGeometryBuffer_);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(ProxyGeometryGPU) * s, NULL, GL_STATIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    //storageBufferInitialized_ = true;
    //storageBufferSize_ = s;

    invalidate(INVALID_PROGRAM);
}

void MultiVolumeGeometryRaycaster::initializeStorageBuffer(size_t s) {
    if(storageBufferInitialized_ && (s == storageBufferSize_))
        return;

    if(storageBufferInitialized_) {
        glDeleteBuffers(1, &storageBuffer_);
        storageBufferInitialized_ = false;
    }

    glGenBuffers(1, &storageBuffer_);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, storageBuffer_);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(LinkedListStructGPU) * s, NULL, GL_STATIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    storageBufferInitialized_ = true;
    storageBufferSize_ = s;

    invalidate(INVALID_PROGRAM);
}

void MultiVolumeGeometryRaycaster::initialize() {
    //TODO: check for OpenGL extensions

    VolumeRaycaster::initialize();
    compile();

    //portGroup_.initialize();
    //portGroup_.addPort(outport_);
    //portGroup_.addPort(outport1_);
    //portGroup_.addPort(outport2_);
    //portGroup_.deactivateTargets();
    LGL_ERROR;

    GLint retVal;
    glGetIntegerv(GL_MAX_FRAGMENT_ATOMIC_COUNTERS, &retVal);
    LINFO("GL_MAX_FRAGMENT_ATOMIC_COUNTERS: " << retVal);
    glGetIntegerv(GL_MAX_FRAGMENT_ATOMIC_COUNTER_BUFFERS, &retVal);
    LINFO("GL_MAX_FRAGMENT_ATOMIC_COUNTER_BUFFERS: " << retVal);

    initializeAtomicBuffer(numPages_.get());
    initializeStorageBuffer(storageBufferSizeProp_.get());
    //initializeStorageBuffer(targetStorageBufferSize_);
    initializeProxyGeometryBuffer(256); //TODO

    adjustPropertyVisibilities();
}

void MultiVolumeGeometryRaycaster::deinitialize() {
    //portGroup_.removePort(outport_);
    //portGroup_.removePort(outport1_);
    //portGroup_.removePort(outport2_);
    //portGroup_.deinitialize();
    LGL_ERROR;

    if(atomicBufferInitialized_) {
        glDeleteBuffers(1, &atomicBuffer_);
        atomicBufferInitialized_ = false;
    }

    if(storageBufferInitialized_) {
        glDeleteBuffers(1, &storageBuffer_);
        storageBufferInitialized_ = false;
    }

    VolumeRaycaster::deinitialize();
}

void MultiVolumeGeometryRaycaster::compile() {
    std::string header = generateHeader();

    rcShaderProp_.setHeader(header);
    rcShaderProp_.rebuild();

    oitClearShaderProp_.setHeader(header);
    oitClearShaderProp_.rebuild();

    oitSortShaderProp_.setHeader(header);
    oitSortShaderProp_.rebuild();

    oitShaderProp_.setHeader(header);
    oitShaderProp_.rebuild();
}

bool MultiVolumeGeometryRaycaster::isReady() const {
    //if(!pgPort1_.isReady() || !volumeInport1_.isReady())
        //return false;

    //if(!tmpPort_.isReady() || !headPort_.isReady())
        //return false;

    //check if at least one outport is connected:
    if (!outport1_.isReady()/* && !outport1_.isReady() && !outport2_.isReady()*/)
        return false;

    return true;
}

void MultiVolumeGeometryRaycaster::beforeProcess() {
    VolumeRaycaster::beforeProcess();

    initializeAtomicBuffer(numPages_.get());
    //initializeStorageBuffer(targetStorageBufferSize_);
    initializeStorageBuffer(storageBufferSizeProp_.get());

    // compile program if needed
    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM) {
        PROFILING_BLOCK("compile");
        compile();
    }
    LGL_ERROR;

    tmpPort_.resize(outport1_.getSize());
    headPort_.resize(outport1_.getSize());

    transferFunc1_.setVolume(volumeInport_.getData());
    transferFunc2_.setVolume(volumeInport2_.getData());
    transferFunc3_.setVolume(volumeInport3_.getData());
    transferFunc4_.setVolume(volumeInport4_.getData());

    if(volumeInport_.hasChanged() || volumeInport2_.hasChanged() || volumeInport3_.hasChanged() || volumeInport4_.hasChanged()) {
        tgt::Bounds b;
        if(volumeInport_.hasData())
            b.addVolume(volumeInport_.getData()->getBoundingBox().getBoundingBox());
        if(volumeInport2_.hasData())
            b.addVolume(volumeInport2_.getData()->getBoundingBox().getBoundingBox());
        if(volumeInport3_.hasData())
            b.addVolume(volumeInport3_.getData()->getBoundingBox().getBoundingBox());
        if(volumeInport4_.hasData())
            b.addVolume(volumeInport4_.getData()->getBoundingBox().getBoundingBox());
        if(length(b.diagonal()) != 0.0)
            camera_.adaptInteractionToScene(b);
    }
}

void MultiVolumeGeometryRaycaster::clearDatastructures() {
    // reset atomic counter:
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, atomicBuffer_);
    //GLuint clear = 1; // start at 1, 0 is for clear pixels
    GLuint* clear = new GLuint[atomicBufferSize_];
    for(int i=0; i<atomicBufferSize_; i++)
        clear[i] = 1;
    glBufferSubData(GL_ATOMIC_COUNTER_BUFFER, 0 , sizeof(GLuint) * atomicBufferSize_, clear);
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, 0);
    delete[] clear;

    tmpPort_.activateTarget();

    // activate shader:
    tgt::Shader* shader = oitClearShaderProp_.getShader();
    shader->activate();

    glBindImageTexture(3 /* unit */, headPort_.getColorTexture()->getId(), 0, false, 0, GL_READ_WRITE, GL_R32UI);

    renderQuad();

    shader->deactivate();
}

void MultiVolumeGeometryRaycaster::renderGeometries(const std::vector<ProxyGeometry>& pgs) {
    tmpPort_.activateTarget();
    glClearDepth(1.0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    LGL_ERROR;

    // set modelview and projection matrices
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadMatrix(camera_.get().getProjectionMatrix(outport1_.getSize()));
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadMatrix(camera_.get().getViewMatrix());
    LGL_ERROR;

    // activate shader:
    tgt::Shader* shader = oitShaderProp_.getShader();
    shader->activate();

    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, atomicBuffer_);
    glBindBufferBase(GL_ATOMIC_COUNTER_BUFFER, 0, atomicBuffer_);
    glBindImageTexture(3 /* unit */, headPort_.getColorTexture()->getId(), 0, false, 0, GL_READ_WRITE, GL_R32UI);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, storageBuffer_);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, proxyGeometryBuffer_);

    // render proxy-geometries:
    for(size_t i=0; i<pgs.size(); i++) {
        shader->setUniform("proxyGeometryId_", (int)(i+1));
        pgs[i].g_->render();
    }

    // render other geometries:
    shader->setUniform("proxyGeometryId_", 0);
    if(geometryPort_.isReady()) {
        for(size_t i = 0; i < geometryPort_.getAllData().size(); i++)
            geometryPort_.getAllData().at(i)->render();
    }

    shader->deactivate();
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, 0);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    // reset matrices:
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();
    LGL_ERROR;

    tmpPort_.deactivateTarget();
}

void MultiVolumeGeometryRaycaster::sortFragments() {
    tmpPort_.activateTarget();
    glClearDepth(1.0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    LGL_ERROR;

    // activate shader:
    tgt::Shader* shader = oitSortShaderProp_.getShader();
    shader->activate();

    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, atomicBuffer_);
    glBindBufferBase(GL_ATOMIC_COUNTER_BUFFER, 0, atomicBuffer_);
    glBindImageTexture(3 /* unit */, headPort_.getColorTexture()->getId(), 0, false, 0, GL_READ_WRITE, GL_R32UI);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, storageBuffer_);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, proxyGeometryBuffer_);

    tgt::Camera cam = camera_.get();
    setGlobalShaderParameters(shader, &cam);

    renderQuad();

    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, 0);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    shader->deactivate();

    LGL_ERROR;

    tmpPort_.deactivateTarget();
}

bool MultiVolumeGeometryRaycaster::adjustStorageBufferSize() {
    // read result:
    GLuint* counter = new GLuint[atomicBufferSize_];
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, atomicBuffer_);
    glGetBufferSubData(GL_ATOMIC_COUNTER_BUFFER, 0, sizeof(GLuint)*atomicBufferSize_, counter);
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, 0);
    std::string out;
    size_t totalFragments = 0;
    size_t pageSize = storageBufferSize_ / atomicBufferSize_;

    size_t maxPageFill = 0;
    size_t minPageFill = pageSize;
    for(int i=0; i<atomicBufferSize_; i++) {
        size_t numFragmentsPage = counter[i] - 1;
        totalFragments += numFragmentsPage;
        out += itos(numFragmentsPage) + " ";

        maxPageFill = std::max(maxPageFill, numFragmentsPage);
        minPageFill = std::min(minPageFill, numFragmentsPage);
    }
    delete[] counter;

    return false; //TODO
    //return true;
    //LINFO(out << "Fragments (" << totalFragments << " total)");

    if(maxPageFill > pageSize) {
        targetStorageBufferSize_ = static_cast<size_t>(maxPageFill * 1.2f * atomicBufferSize_);
        LWARNING("Increasing page size! " << pageSize);
        return false;
    }
    else if(maxPageFill < (0.7f * pageSize)) {
        targetStorageBufferSize_ = static_cast<size_t>(0.8f * pageSize * atomicBufferSize_);
        LWARNING("decreasing page size! " << pageSize);
        return true; // no need to render twice
    }
    else
        return false;
}

//TODO: code duplication with MultiVolumeRaycaster
float MultiVolumeGeometryRaycaster::getVoxelSamplingStepSize(const VolumeBase* vol, float worldSamplingStepSize) {
    return tgt::min(worldSamplingStepSize * vol->getPhysicalToTextureMatrix().getScalingPart());
}

void MultiVolumeGeometryRaycaster::performRaycasting() {
    // initialize shader
    tgt::Shader* raycastPrg = rcShaderProp_.getShader();
    raycastPrg->activate();

    // vector containing the volumes to bind; is passed to bindVolumes()
    std::vector<VolumeStruct> volumeTextures;
    std::vector<const VolumeBase*> volumeHandles;

    // bind volumes
    TextureUnit volUnit1, volUnit2, volUnit3, volUnit4;
    if (volumeInport_.isReady() && pgPort1_.isReady()) {
        volumeTextures.push_back(VolumeStruct(
                    volumeInport_.getData(),
                    &volUnit1,
                    "volume1_","volumeStruct1_",
                    volumeInport_.getTextureClampModeProperty().getValue(),
                    tgt::vec4(volumeInport_.getTextureBorderIntensityProperty().get()),
                    volumeInport_.getTextureFilterModeProperty().getValue())
                );
        volumeHandles.push_back(volumeInport_.getData());
    }
    if (volumeInport2_.isReady() && pgPort2_.isReady()) {
        volumeTextures.push_back(VolumeStruct(
                    volumeInport2_.getData(),
                    &volUnit2,
                    "volume2_","volumeStruct2_",
                    volumeInport2_.getTextureClampModeProperty().getValue(),
                    tgt::vec4(volumeInport2_.getTextureBorderIntensityProperty().get()),
                    volumeInport2_.getTextureFilterModeProperty().getValue())
                );
        volumeHandles.push_back(volumeInport2_.getData());
    }
    if (volumeInport3_.isReady() && pgPort3_.isReady()) {
        volumeTextures.push_back(VolumeStruct(
                    volumeInport3_.getData(),
                    &volUnit3,
                    "volume3_","volumeStruct3_",
                    volumeInport3_.getTextureClampModeProperty().getValue(),
                    tgt::vec4(volumeInport3_.getTextureBorderIntensityProperty().get()),
                    volumeInport3_.getTextureFilterModeProperty().getValue())
                );
        volumeHandles.push_back(volumeInport3_.getData());
    }
    if (volumeInport4_.isReady() && pgPort4_.isReady()) {
        volumeTextures.push_back(VolumeStruct(
                    volumeInport4_.getData(),
                    &volUnit4,
                    "volume4_","volumeStruct4_",
                    volumeInport4_.getTextureClampModeProperty().getValue(),
                    tgt::vec4(volumeInport4_.getTextureBorderIntensityProperty().get()),
                    volumeInport4_.getTextureFilterModeProperty().getValue())
                );
        volumeHandles.push_back(volumeInport4_.getData());
    }

    //portGroup_.activateTargets(); //TODO
    //portGroup_.clearTargets();
    outport1_.activateTarget();
    glClearDepth(1.0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    LGL_ERROR;
    LGL_ERROR;

    glBindImageTexture(3 /* unit */, headPort_.getColorTexture()->getId(), 0, false, 0, GL_READ_WRITE, GL_R32UI);
    //glBindImageTexture(3 [> unit <], headPort_.getColorTexture()->getId(), 0, false, 0, GL_READ_ONLY, GL_R32UI); //TODO

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, storageBuffer_);

    // set common uniforms used by all shaders
    tgt::Camera cam = camera_.get();
    setGlobalShaderParameters(raycastPrg, &cam);
    raycastPrg->setUniform("near_", cam.getNearDist());
    raycastPrg->setUniform("far_", cam.getFarDist());
    // bind the volumes and pass the necessary information to the shader
    bindVolumes(raycastPrg, volumeTextures, &cam, lightPosition_.get());

    // determine ray step length in world coords
    float samplingStepSizeWorld = 0.0f;
    if (volumeTextures.size() > 0) {
        float voxelSizeWorld = 999.f;
        float voxelSizeTexture = 999.f;
        for(size_t i=0; i<volumeHandles.size(); ++i) {
            const VolumeBase* volume = volumeHandles[i];
            tgtAssert(volume, "No volume");
            tgt::ivec3 volDim = volume->getDimensions();
            tgt::vec3 cubeSizeWorld = volume->getCubeSize() * volume->getPhysicalToWorldMatrix().getScalingPart();

            float tVoxelSizeWorld = tgt::max(cubeSizeWorld / tgt::vec3(volDim));
            if (tVoxelSizeWorld < voxelSizeWorld) {
                voxelSizeWorld = tVoxelSizeWorld;
                voxelSizeTexture = tgt::max(1.f / tgt::vec3(volDim));
            }
        }

        samplingStepSizeWorld = voxelSizeWorld / samplingRate_.get();
        float samplingStepSizeTexture = voxelSizeTexture / samplingRate_.get();

        if (QualityMode.isInteractionMode()) {
            samplingStepSizeWorld /= interactionQuality_.get();
            samplingStepSizeTexture /= interactionQuality_.get();
        }

        raycastPrg->setUniform("samplingStepSize_", samplingStepSizeWorld);
        if (compositingMode1_.isSelected("dvr") )// ||
            //(compositingMode1_.isSelected("dvr") && outport1_.isConnected()) ||
            //(compositingMode2_.isSelected("dvr") && outport2_.isConnected()) ) {
            // adapts the compositing of the multivolume RC to the one of the singlevolume RC (see below).
            raycastPrg->setUniform("mvOpacityCorrectionFactor_", samplingStepSizeTexture / samplingStepSizeWorld);
        //}
        LGL_ERROR;
    }
    LGL_ERROR;

    // bind transfer functions
    TextureUnit transferUnit1, transferUnit2, transferUnit3, transferUnit4;
    transferUnit1.activate();
    if (transferFunc1_.get() && volumeInport_.getData())
        ClassificationModes::bindTexture(classificationMode_.get(), transferFunc1_.get(), getVoxelSamplingStepSize(volumeInport_.getData(), samplingStepSizeWorld));

    transferUnit2.activate();
    if (transferFunc2_.get() && volumeInport2_.getData())
        ClassificationModes::bindTexture(classificationMode2_.get(), transferFunc2_.get(), getVoxelSamplingStepSize(volumeInport2_.getData(), samplingStepSizeWorld));

    transferUnit3.activate();
    if (transferFunc3_.get() && volumeInport3_.getData())
        ClassificationModes::bindTexture(classificationMode3_.get(), transferFunc3_.get(), getVoxelSamplingStepSize(volumeInport3_.getData(), samplingStepSizeWorld));

    transferUnit4.activate();
    if (transferFunc4_.get() && volumeInport4_.getData())
        ClassificationModes::bindTexture(classificationMode4_.get(), transferFunc4_.get(), getVoxelSamplingStepSize(volumeInport4_.getData(), samplingStepSizeWorld));


    raycastPrg->setUniform("clippingGradientDepth_", clippingGradientDepth_.get());

    if (compositingMode1_.get() ==  "iso" ||
        compositingMode1_.get() == "iso" ||
        compositingMode2_.get() == "iso")
        raycastPrg->setUniform("isoValue_", isoValue_.get());

    if(volumeInport_.isReady() && pgPort1_.isReady() && ClassificationModes::usesTransferFunction(classificationMode_.get()) ) {
        transferFunc1_.get()->setUniform(raycastPrg, "transferFunc1_", "transferFuncTex1_", transferUnit1.getUnitNumber());
    }
    if(volumeInport2_.isReady() && ClassificationModes::usesTransferFunction(classificationMode_.get())) {
        transferFunc2_.get()->setUniform(raycastPrg, "transferFunc2_", "transferFuncTex2_", transferUnit2.getUnitNumber());
    }
    if(volumeInport3_.isReady() && ClassificationModes::usesTransferFunction(classificationMode_.get())) {
        transferFunc3_.get()->setUniform(raycastPrg, "transferFunc3_", "transferFuncTex3_", transferUnit3.getUnitNumber());
    }
    if(volumeInport4_.isReady() && ClassificationModes::usesTransferFunction(classificationMode_.get())) {
        transferFunc4_.get()->setUniform(raycastPrg, "transferFunc4_", "transferFuncTex4_", transferUnit4.getUnitNumber());
    }

    renderQuad(); // perform the actual rendering by drawing a screen-aligned quad

    raycastPrg->deactivate();
    //portGroup_.deactivateTargets(); //TODO
    outport1_.deactivateTarget();
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    glActiveTexture(GL_TEXTURE0);
    LGL_ERROR;
}

void MultiVolumeGeometryRaycaster::gatherProxyGeometriesRecursive(const Geometry* g, std::vector<ProxyGeometry>& pgs, int volumeId, tgt::mat4 transformation, const tgt::mat4& physicalToTextureMatrix) {
    const GeometrySequence* gs = dynamic_cast<const GeometrySequence*>(g);
    if(gs) {
        for(size_t i=0; i<gs->getNumGeometries(); i++) {
            gatherProxyGeometriesRecursive(gs->getGeometry(i), pgs, volumeId, transformation * gs->getTransformationMatrix(), physicalToTextureMatrix);
        }
    }
    else {
        ProxyGeometry pg;
        pg.textureToWorldMatrix_ = transformation * g->getTransformationMatrix();
        tgt::mat4 pToW = transformation * g->getTransformationMatrix() * physicalToTextureMatrix;
        tgt::mat4 wToP = tgt::mat4::identity;
        pToW.invert(wToP);
        pg.cameraPositionPhysical_ = wToP * camera_.get().getPosition();
        pg.lightPositionPhysical_ = wToP * lightPosition_.get().xyz();
        pg.g_ = g;
        pg.volumeId_ = volumeId;
        pgs.push_back(pg);
    }
}

void MultiVolumeGeometryRaycaster::updateProxyGeometryBuffer(const std::vector<ProxyGeometry>& pgs) {
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, proxyGeometryBuffer_);
    std::vector<ProxyGeometryGPU> pgsGPU;

    pgsGPU.push_back(ProxyGeometryGPU()); //0 = no PG
    for(size_t i=0; i<pgs.size(); i++) {
        ProxyGeometryGPU pgGPU;
        pgGPU.textureToWorldMatrix_ = pgs[i].textureToWorldMatrix_;
        pgGPU.textureToWorldMatrix_.invert(pgGPU.worldToTextureMatrix_);

        // transpose matrices for OpenGL-compatibility:
        pgGPU.textureToWorldMatrix_ = transpose(pgGPU.textureToWorldMatrix_);
        pgGPU.worldToTextureMatrix_ = transpose(pgGPU.worldToTextureMatrix_);

        pgGPU.cameraPositionPhysical_ = pgs[i].cameraPositionPhysical_;
        pgGPU.lightPositionPhysical_ = pgs[i].lightPositionPhysical_;

        pgGPU.volumeId_ = pgs[i].volumeId_;
        pgsGPU.push_back(pgGPU);
    }

    glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0 , sizeof(ProxyGeometryGPU) * pgsGPU.size(), &pgsGPU[0]);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
}

void MultiVolumeGeometryRaycaster::process() {
    // gather all proxygeometries:
    std::vector<ProxyGeometry> pgs;

    if (volumeInport_.isReady() && pgPort1_.isReady()) {
        gatherProxyGeometriesRecursive(pgPort1_.getData(), pgs, 1, tgt::mat4::identity, volumeInport_.getData()->getPhysicalToTextureMatrix());
    }
    if (volumeInport2_.isReady() && pgPort2_.isReady()) {
        gatherProxyGeometriesRecursive(pgPort2_.getData(), pgs, 2, tgt::mat4::identity, volumeInport2_.getData()->getPhysicalToTextureMatrix());
    }
    if (volumeInport3_.isReady() && pgPort3_.isReady()) {
        gatherProxyGeometriesRecursive(pgPort3_.getData(), pgs, 3, tgt::mat4::identity, volumeInport3_.getData()->getPhysicalToTextureMatrix());
    }
    if (volumeInport4_.isReady() && pgPort4_.isReady()) {
        gatherProxyGeometriesRecursive(pgPort4_.getData(), pgs, 4, tgt::mat4::identity, volumeInport4_.getData()->getPhysicalToTextureMatrix());
    }

    updateProxyGeometryBuffer(pgs);

    unsigned int queryID[4];
    glGenQueries(4, queryID);

    // issue the first query
    // Records the time only after all previous
    // commands have been completed
    glQueryCounter(queryID[0], GL_TIMESTAMP);

    clearDatastructures();

    glQueryCounter(queryID[1], GL_TIMESTAMP);

    renderGeometries(pgs);

    sortFragments();

    glQueryCounter(queryID[2], GL_TIMESTAMP);

    performRaycasting();

    glQueryCounter(queryID[3], GL_TIMESTAMP);

    //TODO: reading atomic counters is SLOOOOOOW!
    //if(adjustSize_.get() && adjustStorageBufferSize()) {
        //TODO: prevent processor from being flagged as valid in afterProcess()
    //}

    if(!benchmark_.get()) {
        glDeleteQueries(4, queryID);
        return;
    }
    // wait until the results are available
    GLint stopTimerAvailable = 0;
    while (!stopTimerAvailable)
        glGetQueryObjectiv(queryID[3], GL_QUERY_RESULT_AVAILABLE, &stopTimerAvailable);

    GLuint64 startTime, afterClearTime, afterGeomTime, stopTime;
    // get query results
    glGetQueryObjectui64v(queryID[0], GL_QUERY_RESULT, &startTime);
    glGetQueryObjectui64v(queryID[1], GL_QUERY_RESULT, &afterClearTime);
    glGetQueryObjectui64v(queryID[2], GL_QUERY_RESULT, &afterGeomTime);
    glGetQueryObjectui64v(queryID[3], GL_QUERY_RESULT, &stopTime);
    glDeleteQueries(4, queryID);

    printf("Time spent on the GPU: %f ms\n", (stopTime - startTime) / 1000000.0);
    printf("Clear: %f ms\n", (afterClearTime - startTime) / 1000000.0);
    printf("Geometry: %f ms\n", (afterGeomTime - afterClearTime) / 1000000.0);
    printf("Raycast: %f ms\n", (stopTime - afterGeomTime) / 1000000.0);

}

std::string MultiVolumeGeometryRaycaster::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string headerSource = VolumeRaycaster::generateHeader();

    if(atomicBufferInitialized_) {
        //if(atomicBufferSize_ > 1) //TODO?
            headerSource += "#define NUM_PAGES " + itos(atomicBufferSize_) + "\n";
    }
    else
        headerSource += "#define NUM_PAGES 1\n";

    headerSource += "#define MAX_DEPTHCOMPLEXITY " + itos(maxDepthComplexity_.get()) + "\n";

    if(storageBufferInitialized_)
        headerSource += "#define PAGE_SIZE " + itos(storageBufferSize_ / atomicBufferSize_) + "\n";
    else
        headerSource += "#define PAGE_SIZE 1\n";

    headerSource += "#define " + sortingAlgorithm_.get() + "\n";

    if(volumeInport_.isReady())
        headerSource += "#define VOLUME_1_ACTIVE\n";
    if(volumeInport2_.isReady())
        headerSource += "#define VOLUME_2_ACTIVE\n";
    if(volumeInport3_.isReady())
        headerSource += "#define VOLUME_3_ACTIVE\n";
    if(volumeInport4_.isReady())
        headerSource += "#define VOLUME_4_ACTIVE\n";

    headerSource += ClassificationModes::getShaderDefineSamplerType(classificationMode_.get(), transferFunc1_.get(), "TF_SAMPLER_TYPE_1");
    headerSource += ClassificationModes::getShaderDefineSamplerType(classificationMode2_.get(), transferFunc2_.get(), "TF_SAMPLER_TYPE_2");
    headerSource += ClassificationModes::getShaderDefineSamplerType(classificationMode3_.get(), transferFunc3_.get(), "TF_SAMPLER_TYPE_3");
    headerSource += ClassificationModes::getShaderDefineSamplerType(classificationMode4_.get(), transferFunc4_.get(), "TF_SAMPLER_TYPE_4");

    headerSource += ClassificationModes::getShaderDefineFunction(classificationMode_.get(), "RC_APPLY_CLASSIFICATION");
    headerSource += ClassificationModes::getShaderDefineFunction(classificationMode2_.get(), "RC_APPLY_CLASSIFICATION2");
    headerSource += ClassificationModes::getShaderDefineFunction(classificationMode3_.get(), "RC_APPLY_CLASSIFICATION3");
    headerSource += ClassificationModes::getShaderDefineFunction(classificationMode4_.get(), "RC_APPLY_CLASSIFICATION4");

    // configure shading mode
    headerSource += GLSL::getShaderDefine(shadeMode1_.get(), "APPLY_SHADING_1");
    headerSource += GLSL::getShaderDefine(shadeMode2_.get(), "APPLY_SHADING_2");
    headerSource += GLSL::getShaderDefine(shadeMode3_.get(), "APPLY_SHADING_3");
    headerSource += GLSL::getShaderDefine(shadeMode4_.get(), "APPLY_SHADING_4");

    // DVR opacity correction function adapting the MV compositing to the SVRC compositing,
    // used by the compositing macros below.
    // The adaption is necessary, because the multivolume RC samples in world space
    // instead of in texture space. Due to differing sampling base intervals, we would otherwise
    // still get correct compositing results, but the compositing would slightly differ from
    // the one performed by the SingleVolumeRaycaster.
    headerSource += "uniform float mvOpacityCorrectionFactor_;\n";
    headerSource += "vec4 mvOpacityCorrection(in vec4 color) {\n";
    headerSource += "  return vec4(color.rgb, 1.0 - pow(1.0-color.a, mvOpacityCorrectionFactor_));\n";
    headerSource += "}\n";

    //headerSource += outport1_.generateHeader(rcShaderProp_.getShader());
    return headerSource;
}

void MultiVolumeGeometryRaycaster::adjustPropertyVisibilities() {
    bool useLighting = !shadeMode1_.isSelected("none") |
                       !shadeMode2_.isSelected("none") |
                       !shadeMode3_.isSelected("none") |
                       !shadeMode4_.isSelected("none");
    setPropertyGroupVisible("lighting", useLighting);

    bool useIsovalue = (compositingMode1_.isSelected("iso")  ||
        compositingMode1_.isSelected("iso") ||
        compositingMode2_.isSelected("iso")   );
    isoValue_.setVisibleFlag(useIsovalue);

    lightAttenuation_.setVisibleFlag(applyLightAttenuation_.get());
}

} // namespace

#endif //GL_ATOMIC_COUNTER_BUFFER
