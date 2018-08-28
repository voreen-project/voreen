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

#include "multisliceviewer.h"
#include "voreen/core/utils/glsl.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"

#include "tgt/glmath.h"
#include "tgt/textureunit.h"
#include "tgt/immediatemode/immediatemode.h"

using tgt::ivec2;
using tgt::vec3;
using tgt::mat4;
using tgt::Texture;
using tgt::TextureUnit;

namespace voreen {

MultiSliceViewer::SliceCache::SliceCache(int cacheSize, float samplingRate) : cacheSize_(cacheSize), samplingRate_(samplingRate) {
}

MultiSliceViewer::SliceCache::~SliceCache() {
}

SliceTexture* MultiSliceViewer::SliceCache::getSlice(tgt::plane pl, const VolumeBase* vh) const {
    std::string volumeId = vh->getHash();
    for(std::list<CacheEntry>::iterator it=slices_.begin(); it != slices_.end(); it++) {
        if( ((*it).plane_.toVec4() == pl.toVec4()) &&
            ((*it).volumeId_ == volumeId) ) {
            // move to front:
            CacheEntry ce = *it;
            slices_.erase(it);
            slices_.push_front(ce);

            LDEBUGC("", "Returning cached slice!");
            return ce.slice_;
        }
    }
    SliceTexture* sl = SliceHelper::getVolumeSlice(vh, pl, samplingRate_);
    setSlice(sl, pl, vh);
    LDEBUGC("", "Returning new slice!");
    return sl;
}

void MultiSliceViewer::SliceCache::setSlice(SliceTexture* sl, tgt::plane pl, const VolumeBase* vh) const {
    slices_.push_front(CacheEntry(pl, vh->getHash(), sl));

    // Limit size:
    while(static_cast<int>(slices_.size()) > cacheSize_) {
        delete slices_.back().slice_;
        slices_.pop_back();
    }
}

void MultiSliceViewer::SliceCache::setCacheSize(int cacheSize) {
    cacheSize_ = cacheSize;
    while(static_cast<int>(slices_.size()) > cacheSize_) {
        delete slices_.back().slice_;
        slices_.pop_back();
    }
}

void MultiSliceViewer::SliceCache::setSamplingRate(float samplingRate) {
    if(samplingRate != samplingRate_) {
        clear();
        samplingRate_ = samplingRate;
    }
}

void MultiSliceViewer::SliceCache::clear() {
    while(!slices_.empty()) {
        delete slices_.back().slice_;
        slices_.pop_back();
    }
}

//-----------------------------------------------------------------------------

//TODO: Also in multislicerenderer.cpp
void fillBlendingProp2(StringOptionProperty& sop) {
    sop.addOption("maximumAlpha", "Maximum Alpha");
    sop.addOption("minimumAlpha", "Minimum Alpha");
    sop.addOption("replace", "Replace");
    sop.addOption("deactivate", "deactivate");
    sop.addOption("add", "Add");
    sop.addOption("subtract", "Subtract");
    sop.addOption("multiply", "Multiply");
    sop.addOption("divide", "Divide");
    sop.addOption("difference", "Difference");
    sop.addOption("decal", "Decal");
    sop.addOption("over", "Over");
    sop.addOption("setHue", "Set Hue");
    sop.addOption("blend", "Blend");
}

const std::string MultiSliceViewer::loggerCat_("voreen.MultiSliceViewer");

MultiSliceViewer::MultiSliceViewer()
    : VolumeRenderer()
    , sh_(0)
    , eepShader_(0)
    , sliceCache_(20, 2.0f)
    , mainInport_("mainInport", "Main Volume")
    , restrictToMainVolume_("restrictToMainVolume", "Restrict Rendering to Main Volume", false)
    , sliceAlignment_("sliceAlignmentProp", "Slice Alignment")
    , xSliceIndexProp_("xSliceIndex", "X Slice Number", 0, 0, 10000)
    , ySliceIndexProp_("ySliceIndex", "Y Slice Number", 0, 0, 10000)
    , zSliceIndexProp_("zSliceIndex", "Z Slice Number", 0, 0, 10000)
    , camera_("camera", "Camera")
    , alignCameraButton_("alignCameraButton", "Align Camera")
    , renderCrosshair_("renderCrosshair", "Render Crosshair")
    , crosshairColor_("crosshairColor", "Crosshair Color", tgt::Color(0.f, 1.f, 0.f, 1.f))
    , crosshairWidth_("crosshairWidth", "Crosshair Width", 1.f, 0.1f, 10.f)
    , crosshairRadius_("crosshairRadius", "Crosshair Radius", 15.f, 1.f, 50.f)
    , grabbedX_(false)
    , grabbedY_(false)
    , plane_("plane", "Plane", vec3(1.0f, 0.0f, 0.0f), vec3(-5.0f), vec3(5.0f), Processor::VALID)
    , planeDist_("planeDist", "Plane Distance", 0.0f, -1000.0f, 1000.0f, Processor::VALID)
    , transferFunc1_("transferFunction1", "Transfer function 1")
    , transferFunc2_("transferFunction2", "Transfer function 2")
    , transferFunc3_("transferFunction3", "Transfer function 3")
    , transferFunc4_("transferFunction4", "Transfer function 4")
    , blendingMode1_("blendingMode1", "Blending 1", Processor::INVALID_PROGRAM)
    , blendingMode2_("blendingMode2", "Blending 2", Processor::INVALID_PROGRAM)
    , blendingMode3_("blendingMode3", "Blending 3", Processor::INVALID_PROGRAM)
    , blendingMode4_("blendingMode4", "Blending 4", Processor::INVALID_PROGRAM)
    , blendingFactor1_("blendingFactor1", "Blending Factor 1", 0.5f, 0.0f, 1.0f)
    , blendingFactor2_("blendingFactor2", "Blending Factor 2", 0.5f, 0.0f, 1.0f)
    , blendingFactor3_("blendingFactor3", "Blending Factor 3", 0.5f, 0.0f, 1.0f)
    , blendingFactor4_("blendingFactor4", "Blending Factor 4", 0.5f, 0.0f, 1.0f)
    , texFilterMode1_("textureFilterMode1_", "Texture Filtering 1")
    , texFilterMode2_("textureFilterMode2_", "Texture Filtering 2")
    , texFilterMode3_("textureFilterMode3_", "Texture Filtering 3")
    , texFilterMode4_("textureFilterMode4_", "Texture Filtering 4")
    , texClampMode1_("textureClampMode1_", "Texture Clamp 1")
    , texClampMode2_("textureClampMode2_", "Texture Clamp 2")
    , texClampMode3_("textureClampMode3_", "Texture Clamp 3")
    , texClampMode4_("textureClampMode4_", "Texture Clamp 4")
    , texBorderIntensity_("textureBorderIntensity", "Texture Border Intensity", 0.f)
    , mwheelCycleHandler_("mouseWheelHandler", "Slice Cycling", &zSliceIndexProp_)
    , interactionHandler_("interactionHandler", "Camera Interaction", &camera_)
    , texMode1_("textureMode1", "Texture Mode 1", Processor::INVALID_PROGRAM)
    , texMode2_("textureMode2", "Texture Mode 2", Processor::INVALID_PROGRAM)
    , texMode3_("textureMode3", "Texture Mode 3", Processor::INVALID_PROGRAM)
    , texMode4_("textureMode4", "Texture Mode 4", Processor::INVALID_PROGRAM)
    , samplingRate_("samplingRate", "Sampling Rate", 2.0f, 1.0f, 20.0f)
    , cacheSize_("cacheSize", "# Cached 2D Textures", 20, 1, 1000)
    , inport1_(Port::INPORT, "volume1", "Volume 1", false, Processor::INVALID_PROGRAM)
    , inport2_(Port::INPORT, "volume2", "Volume 2", false, Processor::INVALID_PROGRAM)
    , inport3_(Port::INPORT, "volume3", "Volume 3", false, Processor::INVALID_PROGRAM)
    , inport4_(Port::INPORT, "volume4", "Volume 4", false, Processor::INVALID_PROGRAM)
    , geomPort_(Port::OUTPORT, "geometry")
    , sliceIndexTextPort_(Port::OUTPORT, "text")
    , intensityTextPort_(Port::OUTPORT, "textIntensity")
    , outport_(Port::OUTPORT, "output", "Rendering", true, INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , entryPort_(Port::OUTPORT, "entrypoints", "entrypoints", false, INVALID_RESULT, RenderPort::RENDERSIZE_DEFAULT, GL_RGBA16F)
{
    addInteractionHandler(mwheelCycleHandler_);
    addInteractionHandler(interactionHandler_);

    mainInport_.addOption("1", "p1", 1);
    mainInport_.addOption("2", "p2", 2);
    mainInport_.addOption("3", "p3", 3);
    mainInport_.addOption("4", "p4", 4);
    mainInport_.onChange( MemberFunctionCallback<MultiSliceViewer>(this, &MultiSliceViewer::updateSliceProperties) );
    addProperty(mainInport_);
    addProperty(restrictToMainVolume_);

    sliceAlignment_.addOption("xy-plane", "XY-Plane (axial)", XY_PLANE);
    sliceAlignment_.addOption("xz-plane", "XZ-Plane (coronal)", XZ_PLANE);
    sliceAlignment_.addOption("yz-plane", "YZ-Plane (sagittal)", YZ_PLANE);
    sliceAlignment_.onChange( MemberFunctionCallback<MultiSliceViewer>(this, &MultiSliceViewer::updateSliceProperties) );
    addProperty(sliceAlignment_);

    addProperty(xSliceIndexProp_);
    addProperty(ySliceIndexProp_);
    addProperty(zSliceIndexProp_);

    addProperty(camera_);

    addProperty(alignCameraButton_);
    alignCameraButton_.onChange( MemberFunctionCallback<MultiSliceViewer>(this, &MultiSliceViewer::alignCamera) );

    addProperty(renderCrosshair_);
    renderCrosshair_.onChange( MemberFunctionCallback<MultiSliceViewer>(this, &MultiSliceViewer::adjustPropertyVisibility) );
    addProperty(crosshairColor_);
    addProperty(crosshairWidth_);
    addProperty(crosshairRadius_);

    addProperty(plane_);
    addProperty(planeDist_);

    addProperty(transferFunc1_);
    addProperty(transferFunc2_);
    addProperty(transferFunc3_);
    addProperty(transferFunc4_);

    fillBlendingProp2(blendingMode1_);
    addProperty(blendingMode1_);
    addProperty(blendingFactor1_);

    fillBlendingProp2(blendingMode2_);
    addProperty(blendingMode2_);
    addProperty(blendingFactor2_);

    fillBlendingProp2(blendingMode3_);
    addProperty(blendingMode3_);
    addProperty(blendingFactor3_);

    fillBlendingProp2(blendingMode4_);
    addProperty(blendingMode4_);
    addProperty(blendingFactor4_);

    // volume texture filtering
    texFilterMode1_.addOption("nearest", "Nearest",  GL_NEAREST);
    texFilterMode1_.addOption("linear",  "Linear",   GL_LINEAR);
    texFilterMode1_.selectByKey("linear");
    addProperty(texFilterMode1_);
    texFilterMode2_.addOption("nearest", "Nearest",  GL_NEAREST);
    texFilterMode2_.addOption("linear",  "Linear",   GL_LINEAR);
    texFilterMode2_.selectByKey("linear");
    addProperty(texFilterMode2_);
    texFilterMode3_.addOption("nearest", "Nearest",  GL_NEAREST);
    texFilterMode3_.addOption("linear",  "Linear",   GL_LINEAR);
    texFilterMode3_.selectByKey("linear");
    addProperty(texFilterMode3_);
    texFilterMode4_.addOption("nearest", "Nearest",  GL_NEAREST);
    texFilterMode4_.addOption("linear",  "Linear",   GL_LINEAR);
    texFilterMode4_.selectByKey("linear");
    addProperty(texFilterMode4_);

    // volume texture clamping
    texClampMode1_.addOption("clamp-to-edge",   "Clamp to Edge",     GL_CLAMP_TO_EDGE);
    texClampMode1_.addOption("clamp-to-border", "Clamp to Border",   GL_CLAMP_TO_BORDER);
    texClampMode1_.selectByKey("clamp-to-edge");
    addProperty(texClampMode1_);
    texClampMode2_.addOption("clamp-to-edge",   "Clamp to Edge",     GL_CLAMP_TO_EDGE);
    texClampMode2_.addOption("clamp-to-border", "Clamp to Border",   GL_CLAMP_TO_BORDER);
    texClampMode2_.selectByKey("clamp-to-edge");
    addProperty(texClampMode2_);
    texClampMode3_.addOption("clamp-to-edge",   "Clamp to Edge",     GL_CLAMP_TO_EDGE);
    texClampMode3_.addOption("clamp-to-border", "Clamp to Border",   GL_CLAMP_TO_BORDER);
    texClampMode3_.selectByKey("clamp-to-edge");
    addProperty(texClampMode3_);
    texClampMode4_.addOption("clamp-to-edge",   "Clamp to Edge",     GL_CLAMP_TO_EDGE);
    texClampMode4_.addOption("clamp-to-border", "Clamp to Border",   GL_CLAMP_TO_BORDER);
    texClampMode4_.selectByKey("clamp-to-edge");
    addProperty(texClampMode4_);
    addProperty(texBorderIntensity_);

    texMode1_.addOption("2d-texture", "2D Textures", TEXTURE_2D);
    texMode1_.addOption("3d-texture", "3D Texture", TEXTURE_3D);
    texMode1_.selectByKey("3d-texture");
    addProperty(texMode1_);
    texMode2_.addOption("2d-texture", "2D Textures", TEXTURE_2D);
    texMode2_.addOption("3d-texture", "3D Texture", TEXTURE_3D);
    texMode2_.selectByKey("3d-texture");
    addProperty(texMode2_);
    texMode3_.addOption("2d-texture", "2D Textures", TEXTURE_2D);
    texMode3_.addOption("3d-texture", "3D Texture", TEXTURE_3D);
    texMode3_.selectByKey("3d-texture");
    addProperty(texMode3_);
    texMode4_.addOption("2d-texture", "2D Textures", TEXTURE_2D);
    texMode4_.addOption("3d-texture", "3D Texture", TEXTURE_3D);
    texMode4_.selectByKey("3d-texture");
    addProperty(texMode4_);
    addProperty(samplingRate_);
    addProperty(cacheSize_);

    addPort(inport1_);
    addPort(inport2_);
    addPort(inport3_);
    addPort(inport4_);

    addPort(outport_);
    addPrivateRenderPort(entryPort_);
    addPort(geomPort_);
    addPort(sliceIndexTextPort_);
    addPort(intensityTextPort_);
}

MultiSliceViewer::~MultiSliceViewer() {
}

void MultiSliceViewer::adjustPropertyVisibility() {
    crosshairColor_.setVisibleFlag(renderCrosshair_.get());
    crosshairWidth_.setVisibleFlag(renderCrosshair_.get());
    crosshairRadius_.setVisibleFlag(renderCrosshair_.get());
}

void MultiSliceViewer::initialize() {
    VolumeRenderer::initialize();

    sh_ = ShdrMgr.loadSeparate("passthrough.vert", "multisliceviewer.frag", generateHeader(), false);
    if (!sh_)
        throw VoreenException("Failed to load shaders");

    eepShader_ = ShdrMgr.load("eep_simple", VolumeRenderer::generateHeader(), false);
    if (!eepShader_)
        throw VoreenException("Failed to load shaders");
}

void MultiSliceViewer::deinitialize() {
    ShdrMgr.dispose(sh_);
    sh_ = 0;
    LGL_ERROR;

    ShdrMgr.dispose(eepShader_);
    eepShader_ = 0;
    LGL_ERROR;

    sliceCache_.clear();

    VolumeRenderer::deinitialize();
}

void MultiSliceViewer::compile() {
    sh_->setHeaders(generateHeader());
    sh_->rebuild();
}

const VolumePort* MultiSliceViewer::getMainInport() const {
    switch(mainInport_.getValue()) {
        case 1: return &inport1_;
        case 2: return &inport2_;
        case 3: return &inport3_;
        case 4: return &inport4_;
        default: return 0;
    }
}

const std::vector<const VolumeBase*> MultiSliceViewer::getSecondaryVolumes() const {
    std::vector<const VolumeBase*> sv;

    if((mainInport_.getValue() != 1) && inport1_.isReady())
        sv.push_back(inport1_.getData());

    if((mainInport_.getValue() != 2) && inport2_.isReady())
        sv.push_back(inport2_.getData());

    if((mainInport_.getValue() != 3) && inport3_.isReady())
        sv.push_back(inport3_.getData());

    if((mainInport_.getValue() != 4) && inport4_.isReady())
        sv.push_back(inport4_.getData());

    return sv;
}

std::string MultiSliceViewer::getProcessorInfo() const {
    return "";
}

bool MultiSliceViewer::isReady() const {
    if(getMainInport()->isReady() && outport_.isReady())
        return true;
    else
        return false;
}

const IntProperty* MultiSliceViewer::getSliceIndexProperty() const {
    switch(sliceAlignment_.getValue()) {
        case YZ_PLANE: return &xSliceIndexProp_;
        case XZ_PLANE: return &ySliceIndexProp_;
        case XY_PLANE: return &zSliceIndexProp_;
        default:
            {
                tgtAssert(false, "should not get here!");
                return &xSliceIndexProp_;
            }
    }
}

int MultiSliceViewer::getSliceIndex() const {
    switch(sliceAlignment_.getValue()) {
        case YZ_PLANE: return xSliceIndexProp_.get();
        case XZ_PLANE: return ySliceIndexProp_.get();
        case XY_PLANE: return zSliceIndexProp_.get();
        default:
            {
                return xSliceIndexProp_.get();
                tgtAssert(false, "should not get here!");
            }
    }
}

void MultiSliceViewer::process() {
    sliceCache_.setCacheSize(cacheSize_.get());
    sliceCache_.setSamplingRate(samplingRate_.get());
    LGL_ERROR;

    if (getMainInport()->hasChanged()) {
        updateSliceProperties();  // validate the currently set values and adjust them if necessary
        //TODO: Check if re-alignment is really necessary
    }

    // compile program if needed
    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM)
        compile();
    LGL_ERROR;

    const VolumeBase* volh = getMainInport()->getData();

    // Calculate geometry:
    TriangleMeshGeometryNormal* slice = SliceHelper::getSliceGeometry(volh, sliceAlignment_.getValue(), (float)getSliceIndex(), true, restrictToMainVolume_.get() ? std::vector<const VolumeBase*>() : getSecondaryVolumes());
    geomPort_.setData(slice);

    //calculate plane equation:
    tgt::mat4 m = slice->getTransformationMatrix();
    tgt::plane p(m * slice->getTriangle(0).v_[0].pos_, m * slice->getTriangle(0).v_[1].pos_, m * slice->getTriangle(0).v_[2].pos_);

    plane_.set(p.n);
    planeDist_.set(p.d);

    // render entry points:
    entryPort_.activateTarget();

    LGL_ERROR;
    glClearDepth(1.0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    LGL_ERROR;

    // set modelview and projection matrices
    tgt::Camera cam = camera_.get();
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.loadMatrix(cam.getProjectionMatrix(outport_.getSize()));

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.loadMatrix(cam.getViewMatrix());

    MatStack.multMatrix(slice->getTransformationMatrix());

    // activate shader program
    eepShader_->activate();
    setGlobalShaderParameters(eepShader_, &cam);
    eepShader_->setUniform("useTextureCoordinates_", false);
    //mat4 viewToWorldMatrix = mat4::identity;
    //cam.getViewMatrix().invert(viewToWorldMatrix);
    //eepShader_->setUniform("inverseViewMatrix_", viewToWorldMatrix);
    LGL_ERROR;

    slice->render();
    LGL_ERROR;

    entryPort_.deactivateTarget();

    // deactivate shader program
    eepShader_->deactivate();

    // restore OpenGL state
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();
    LGL_ERROR;

    TextureUnit::setZeroUnit();
    LGL_ERROR;
    //----------------------------------------------------------------------------------
    // Perform the actual rendering:

    // bind transfer function
    TextureUnit transferUnit, transferUnit2, transferUnit3, transferUnit4;
    transferUnit.activate();
    if (transferFunc1_.get())
        transferFunc1_.get()->getTexture()->bind();

    transferUnit2.activate();
    if (transferFunc2_.get())
        transferFunc2_.get()->getTexture()->bind();

    transferUnit3.activate();
    if (transferFunc3_.get())
        transferFunc3_.get()->getTexture()->bind();

    transferUnit4.activate();
    if (transferFunc4_.get())
        transferFunc4_.get()->getTexture()->bind();

    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    LGL_ERROR;

    transferFunc1_.setVolume(inport1_.getData());
    transferFunc2_.setVolume(inport2_.getData());
    transferFunc3_.setVolume(inport3_.getData());
    transferFunc4_.setVolume(inport4_.getData());

    TextureUnit entryUnit, entryDepthUnit;
    // bind entry params
    entryPort_.bindTextures(entryUnit.getEnum(), entryDepthUnit.getEnum());
    LGL_ERROR;

    // initialize shader
    sh_->activate();

    // vector containing the volumes to bind; is passed to bindVolumes()
    std::vector<VolumeStruct> volumeTextures;

    TextureUnit volUnit1, volUnit2, volUnit3, volUnit4;
    // add volumes
    if(inport1_.isReady()) {
        if(texMode1_.getValue() == TEXTURE_3D) {
            volumeTextures.push_back(VolumeStruct(
                        inport1_.getData(),
                        &volUnit1,
                        "volume1_","volumeStruct1_",
                        texClampMode1_.getValue(),
                        tgt::vec4(texBorderIntensity_.get()),
                        texFilterMode1_.getValue())
                    );
        }
        else {
            SliceTexture* slice = sliceCache_.getSlice(p, inport1_.getData());

            sh_->setIgnoreUniformLocationError(true); // we do not need all of the information, so we ignore missing uniforms
            GLSL::setUniform(sh_, "volume1_", "volumeStruct1_", slice, &volUnit1);
            sh_->setIgnoreUniformLocationError(false);

            sh_->setUniform("volumeStruct1_.worldToTextureMatrix_", slice->getWorldToTextureMatrix());
            volUnit1.activate();
            slice->bind();
            LGL_ERROR;
        }
    }
    if(inport2_.isReady()) {
        if(texMode2_.getValue() == TEXTURE_3D) {
            volumeTextures.push_back(VolumeStruct(
                        inport2_.getData(),
                        &volUnit2,
                        "volume2_","volumeStruct2_",
                        texClampMode2_.getValue(),
                        tgt::vec4(texBorderIntensity_.get()),
                        texFilterMode2_.getValue())
                    );
        }
        else {
            volUnit2.activate();
            SliceTexture* slice = sliceCache_.getSlice(p, inport2_.getData());

            sh_->setIgnoreUniformLocationError(true); // we do not need all of the information, so we ignore missing uniforms
            GLSL::setUniform(sh_, "volume2_", "volumeStruct2_", slice, &volUnit2);
            sh_->setIgnoreUniformLocationError(false);

            sh_->setUniform("volumeStruct2_.worldToTextureMatrix_", slice->getWorldToTextureMatrix());
            slice->bind();
            LGL_ERROR;
        }
    }
    if(inport3_.isReady()) {
        if(texMode3_.getValue() == TEXTURE_3D) {
            volumeTextures.push_back(VolumeStruct(
                        inport3_.getData(),
                        &volUnit3,
                        "volume3_","volumeStruct3_",
                        texClampMode3_.getValue(),
                        tgt::vec4(texBorderIntensity_.get()),
                        texFilterMode3_.getValue())
                    );
        }
        else {
            volUnit3.activate();
            SliceTexture* slice = sliceCache_.getSlice(p, inport3_.getData());

            sh_->setIgnoreUniformLocationError(true); // we do not need all of the information, so we ignore missing uniforms
            GLSL::setUniform(sh_, "volume3_", "volumeStruct3_", slice, &volUnit3);
            sh_->setIgnoreUniformLocationError(false);

            sh_->setUniform("volumeStruct3_.worldToTextureMatrix_", slice->getWorldToTextureMatrix());
            slice->bind();
            LGL_ERROR;
        }
    }
    if(inport4_.isReady()) {
        if(texMode4_.getValue() == TEXTURE_3D) {
            volumeTextures.push_back(VolumeStruct(
                        inport4_.getData(),
                        &volUnit4,
                        "volume4_","volumeStruct4_",
                        texClampMode4_.getValue(),
                        tgt::vec4(texBorderIntensity_.get()),
                        texFilterMode4_.getValue())
                    );
        }
        else {
            volUnit4.activate();
            SliceTexture* slice = sliceCache_.getSlice(p, inport4_.getData());

            sh_->setIgnoreUniformLocationError(true); // we do not need all of the information, so we ignore missing uniforms
            GLSL::setUniform(sh_, "volume4_", "volumeStruct4_", slice, &volUnit4);
            sh_->setIgnoreUniformLocationError(false);

            sh_->setUniform("volumeStruct4_.worldToTextureMatrix_", slice->getWorldToTextureMatrix());
            slice->bind();
            LGL_ERROR;
        }
    }
    LGL_ERROR;

    // set common uniforms used by all shaders
    setGlobalShaderParameters(sh_, &cam);

    // bind the volumes and pass the necessary information to the shader
    bindVolumes(sh_, volumeTextures);

    // pass the remaining uniforms to the shader
    sh_->setUniform("entryPoints_", entryUnit.getUnitNumber());
    sh_->setUniform("entryPointsDepth_",entryDepthUnit.getUnitNumber());
    entryPort_.setTextureParameters(sh_, "entryParameters_");

    if(inport1_.isReady()) {
        transferFunc1_.get()->setUniform(sh_, "transferFunc1_", "transferFuncTex1_", transferUnit.getUnitNumber());

        sh_->setIgnoreUniformLocationError(true);
        sh_->setUniform("blendingFactor1_", blendingFactor1_.get());
        sh_->setIgnoreUniformLocationError(false);
    }
    if(inport2_.isReady()) {
        transferFunc2_.get()->setUniform(sh_, "transferFunc2_", "transferFuncTex2_", transferUnit2.getUnitNumber());

        sh_->setIgnoreUniformLocationError(true);
        sh_->setUniform("blendingFactor2_", blendingFactor2_.get());
        sh_->setIgnoreUniformLocationError(false);
    }
    if(inport3_.isReady()) {
        transferFunc3_.get()->setUniform(sh_, "transferFunc3_", "transferFuncTex3_", transferUnit3.getUnitNumber());

        sh_->setIgnoreUniformLocationError(true);
        sh_->setUniform("blendingFactor3_", blendingFactor3_.get());
        sh_->setIgnoreUniformLocationError(false);
    }
    if(inport4_.isReady()) {
        transferFunc4_.get()->setUniform(sh_, "transferFunc4_", "transferFuncTex4_", transferUnit4.getUnitNumber());

        sh_->setIgnoreUniformLocationError(true);
        sh_->setUniform("blendingFactor4_", blendingFactor4_.get());
        sh_->setIgnoreUniformLocationError(false);
    }

    LGL_ERROR;

    renderQuad();
    LGL_ERROR;

    sh_->deactivate();
    LGL_ERROR;

    glActiveTexture(GL_TEXTURE0);
    LGL_ERROR;

    if(renderCrosshair_.get()) {
        glLineWidth(crosshairWidth_.get());
        glDepthFunc(GL_ALWAYS);
        glEnable(GL_BLEND);

        vec3 pVoxel = vec3((float)xSliceIndexProp_.get(), (float)ySliceIndexProp_.get(), (float)zSliceIndexProp_.get()) + vec3(0.5f);
        vec3 p = volh->getVoxelToWorldMatrix() * pVoxel;
        tgt::vec2 c = camera_.get().project(outport_.getSize(), p).xy();

        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.pushMatrix();
        tgt::vec2 ss = outport_.getSize();
        MatStack.loadMatrix(tgt::mat4::createOrtho(0, ss.x, 0, ss.y, 0, 1));

        float circleRadius = crosshairRadius_.get();

        IMode.begin(tgt::ImmediateMode::LINES);
        IMode.color(crosshairColor_.get());
        IMode.vertex(0.0f, c.y);
        IMode.vertex(c.x-circleRadius, c.y);

        IMode.vertex(c.x+circleRadius, c.y);
        IMode.vertex(ss.x, c.y);

        IMode.vertex(c.x, 0.0f);
        IMode.vertex(c.x, c.y-circleRadius);

        IMode.vertex(c.x, c.y+circleRadius);
        IMode.vertex(c.x, ss.y);
        IMode.color(tgt::vec4::one);
        IMode.end();

        IMode.begin(tgt::ImmediateMode::LINE_LOOP);
        int numSegments = 100;
        IMode.color(crosshairColor_.get());
        for(int i=0; i<numSegments; i++) {
            float angle = ((float)i / (float) numSegments) * tgt::PIf * 2.0f;
            IMode.vertex(c.x + (sinf(angle) * circleRadius), c.y + (cosf(angle) * circleRadius));
        }
        IMode.color(tgt::vec4::one);
        IMode.end();

        glDisable(GL_BLEND);
        glDepthFunc(GL_LESS);
        glLineWidth(1.f);
        MatStack.popMatrix();
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        LGL_ERROR;
    }

    outport_.deactivateTarget();

    //Generate slice index text output:
    std::stringstream strstrIndex;
    strstrIndex << getSliceIndex() << "/" << getSliceIndexProperty()->getMaxValue();
    sliceIndexTextPort_.setData(strstrIndex.str());

    //Generate intensity text output:
    std::string strIntensity;
    tgt::vec3 posWorld = volh->getVoxelToWorldMatrix() * (vec3((float)xSliceIndexProp_.get(), (float)ySliceIndexProp_.get(), (float)zSliceIndexProp_.get()) + vec3(0.5f));
    if(inport1_.isReady()) {
        RealWorldMapping rwm = inport1_.getData()->getRealWorldMapping();
        vec3 posVoxel = inport1_.getData()->getWorldToVoxelMatrix() * posWorld;
        if (inport1_.getData()->hasRepresentation<VolumeRAM>())
            strIntensity += inport1_.getData()->getRepresentation<VolumeRAM>()->getVoxelValueAsString(posVoxel, &rwm);
        else
            ; // TODO
    }
    if(inport2_.isReady()) {
        if(!strIntensity.empty())
            strIntensity += "\n";

        RealWorldMapping rwm = inport2_.getData()->getRealWorldMapping();
        vec3 posVoxel = inport2_.getData()->getWorldToVoxelMatrix() * posWorld;
        if (inport2_.getData()->hasRepresentation<VolumeRAM>())
            strIntensity += inport2_.getData()->getRepresentation<VolumeRAM>()->getVoxelValueAsString(posVoxel, &rwm);
        else
            ; // TODO
    }
    if(inport3_.isReady()) {
        if(!strIntensity.empty())
            strIntensity += "\n";

        RealWorldMapping rwm = inport3_.getData()->getRealWorldMapping();
        vec3 posVoxel = inport3_.getData()->getWorldToVoxelMatrix() * posWorld;
        if (inport3_.getData()->hasRepresentation<VolumeRAM>())
            strIntensity += inport3_.getData()->getRepresentation<VolumeRAM>()->getVoxelValueAsString(posVoxel, &rwm);
        else
            ; // TODO
    }
    if(inport4_.isReady()) {
        if(!strIntensity.empty())
            strIntensity += "\n";

        RealWorldMapping rwm = inport4_.getData()->getRealWorldMapping();
        vec3 posVoxel = inport4_.getData()->getWorldToVoxelMatrix() * posWorld;
        if (inport4_.getData()->hasRepresentation<VolumeRAM>())
           strIntensity += inport4_.getData()->getRepresentation<VolumeRAM>()->getVoxelValueAsString(posVoxel, &rwm);
        else
            ; // TODO
    }
    intensityTextPort_.setData(strIntensity);
}

void MultiSliceViewer::updateSliceProperties() {
    if (!getMainInport()->getData())
        return;

    tgt::ivec3 volumeDim = getMainInport()->getData()->getDimensions();
    tgtAssert(sliceAlignment_.getValue() >= 0 && sliceAlignment_.getValue() <= 2, "Invalid alignment value");

    xSliceIndexProp_.setMaxValue(volumeDim.x - 1);
    ySliceIndexProp_.setMaxValue(volumeDim.y - 1);
    zSliceIndexProp_.setMaxValue(volumeDim.z - 1);
    switch(sliceAlignment_.getValue()) {
        case YZ_PLANE: mwheelCycleHandler_.setProperty(&xSliceIndexProp_);
                       break;
        case XZ_PLANE: mwheelCycleHandler_.setProperty(&ySliceIndexProp_);
                       break;
        case XY_PLANE: mwheelCycleHandler_.setProperty(&zSliceIndexProp_);
                       break;
        default: tgtAssert(false, "should not get here!");
    }

    alignCamera();
}

void MultiSliceViewer::alignCamera() {
    const VolumeBase* volh = getMainInport()->getData();

    if(!volh)
        return;

    vec3 urb = volh->getURB();
    vec3 llf = volh->getLLF();

    vec3 center = (urb + llf) / 2.0f;
    vec3 xVec = vec3(llf.x, center.y, center.z);
    vec3 yVec = vec3(center.x, llf.y, center.z);
    vec3 zVec = vec3(center.x, center.y, llf.z);

    tgt::mat4 m = volh->getPhysicalToWorldMatrix();

    center = m * center;
    xVec = m * xVec;
    yVec = m * yVec;
    zVec = m * zVec;

    vec3 mainDir;
    vec3 upDir;
    vec3 strafeDir;

    switch(sliceAlignment_.getValue()) {
        case YZ_PLANE: {
                           mainDir = xVec - center;
                           upDir = -(zVec - center);
                           strafeDir = -(yVec - center);
                       }
                       break;
        case XZ_PLANE: {
                           mainDir = yVec - center;
                           upDir = -(zVec - center);
                           strafeDir = -(xVec - center);
                       }
                       break;
        case XY_PLANE: {
                           mainDir = zVec - center;
                           upDir = (yVec - center);
                           strafeDir = (xVec - center);
                       }
                       break;
        default: tgtAssert(false, "should not get here!");
    }

    //setup camera:
    vec3 pos = center + (mainDir * 1.5f);
    vec3 up = normalize(upDir);

    //setup frustum:
    float size = std::max(length(strafeDir), length(upDir));
    float farDist = length(mainDir) * 3.0f;
    tgt::Frustum f(-size, +size, -size, +size, 0.01f, farDist);

    bool updateCam = false;

    vec3 look = normalize(center - pos);
    tgt::Camera currentCam = camera_.get();

    if(dot(look, currentCam.getLook()) < 0.99f || ( (tgt::length(pos - currentCam.getPosition()) / farDist) > 0.001f) )
        updateCam = true;

    //TODO: more checks, make configureable?

    if(updateCam) {
        tgt::Camera cam(pos, center, up);
        cam.setProjectionMode(tgt::Camera::ORTHOGRAPHIC);

        cam.setFrustum(f);

        camera_.set(cam);
    }
}

void MultiSliceViewer::setIntProperty(IntProperty& p, int value) {
    if(value < p.getMinValue())
        p.set(p.getMinValue());
    else if(value > p.getMaxValue())
        p.set(p.getMaxValue());
    else
        p.set(value);
}

void MultiSliceViewer::onEvent(tgt::Event* e) {
    tgt::MouseEvent* me = dynamic_cast<tgt::MouseEvent*>(e);

    if(me) {
        tgt::ivec2 coord = me->coord();
        coord.y = me->viewport().y - coord.y;

        if(renderCrosshair_.get()) {
            if((me->action() == tgt::MouseEvent::PRESSED) && (me->button() == tgt::MouseEvent::MOUSE_BUTTON_LEFT)) {
                const VolumeBase* volh = getMainInport()->getData();
                if(volh) {
                    vec3 pVoxel = vec3((float)xSliceIndexProp_.get(), (float)ySliceIndexProp_.get(), (float)zSliceIndexProp_.get()) + vec3(0.5f);
                    vec3 p = volh->getVoxelToWorldMatrix() * pVoxel;
                    tgt::vec2 c = camera_.get().project(outport_.getSize(), p).xy();

                    if(distance(tgt::vec2(coord), c) < crosshairRadius_.get()) {
                        grabbedX_ = true;
                        grabbedY_ = true;
                    }
                    else if(abs(coord.x - c.x) < 5)
                        grabbedY_ = true;
                    else if(abs(coord.y - c.y) < 5)
                        grabbedX_ = true;

                    if(grabbedX_ || grabbedY_)
                        e->accept();
                }
            }
            else if((me->action() == tgt::MouseEvent::RELEASED) && (me->button() == tgt::MouseEvent::MOUSE_BUTTON_LEFT)) {
                if(grabbedX_ || grabbedY_) {
                    grabbedX_ = false;
                    grabbedY_ = false;
                    e->accept();
                    return;
                }
            }
            else if(me->action() == tgt::MouseEvent::MOTION) {
                if(grabbedX_ || grabbedY_) {
                    const VolumeBase* volh = getMainInport()->getData();
                    if(volh) {
                        tgt::line3 line = camera_.get().getViewRay(me->viewport(), me->coord());
                        vec3 pWorld = line.getFromParam(0.5f);
                        vec3 pVoxel = volh->getWorldToVoxelMatrix() * pWorld;

                        switch(sliceAlignment_.getValue()) {
                            case YZ_PLANE: {
                                               if(grabbedY_)
                                                   setIntProperty(ySliceIndexProp_, tgt::iround(pVoxel.y));
                                               if(grabbedX_)
                                                   setIntProperty(zSliceIndexProp_, tgt::iround(pVoxel.z));
                                           }
                                           break;
                            case XZ_PLANE: {
                                               if(grabbedY_)
                                                   setIntProperty(xSliceIndexProp_, tgt::iround(pVoxel.x));
                                               if(grabbedX_)
                                                   setIntProperty(zSliceIndexProp_, tgt::iround(pVoxel.z));
                                           }
                                           break;
                            case XY_PLANE: {
                                               if(grabbedY_)
                                                   setIntProperty(xSliceIndexProp_, tgt::iround(pVoxel.x));
                                               if(grabbedX_)
                                                   setIntProperty(ySliceIndexProp_, tgt::iround(pVoxel.y));
                                           }
                                           break;
                            default: tgtAssert(false, "should not get here!");
                        }
                    }
                    e->accept();
                    return;
                }
            }
            else if(me->action() == tgt::MouseEvent::EXIT) {
                //grabbedX_ = false;
                //grabbedY_ = false;
            }
        }

        if(me->action() == tgt::MouseEvent::DOUBLECLICK) {
            const VolumeBase* volh = getMainInport()->getData();
            if(volh) {
                tgt::line3 line = camera_.get().getViewRay(me->viewport(), me->coord());
                vec3 pWorld = line.getFromParam(0.5f);
                vec3 pVoxel = volh->getWorldToVoxelMatrix() * pWorld;

                switch(sliceAlignment_.getValue()) {
                    case YZ_PLANE: {
                                       setIntProperty(ySliceIndexProp_, tgt::iround(pVoxel.y));
                                       setIntProperty(zSliceIndexProp_, tgt::iround(pVoxel.z));
                                   }
                                   break;
                    case XZ_PLANE: {
                                       setIntProperty(xSliceIndexProp_, tgt::iround(pVoxel.x));
                                       setIntProperty(zSliceIndexProp_, tgt::iround(pVoxel.z));
                                   }
                                   break;
                    case XY_PLANE: {
                                       setIntProperty(xSliceIndexProp_, tgt::iround(pVoxel.x));
                                       setIntProperty(ySliceIndexProp_, tgt::iround(pVoxel.y));
                                   }
                                   break;
                    default: tgtAssert(false, "should not get here!");
                }
            }
            e->accept();
            return;
        }
        lastMousePos_ = coord;
    }

    VolumeRenderer::onEvent(e);
}

std::string MultiSliceViewer::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string headerSource = VolumeRenderer::generateHeader();

    if(inport1_.isReady()) {
        headerSource += "#define VOLUME_1_ACTIVE\n";
        if(texMode1_.getValue() == TEXTURE_3D)
            headerSource += "#define VOLUME_1_3D\n";
        else
            headerSource += "#define VOLUME_1_2D\n";
    }
    if(inport2_.isReady()) {
        headerSource += "#define VOLUME_2_ACTIVE\n";
        if(texMode2_.getValue() == TEXTURE_3D)
            headerSource += "#define VOLUME_2_3D\n";
        else
            headerSource += "#define VOLUME_2_2D\n";
    }
    if(inport3_.isReady()) {
        headerSource += "#define VOLUME_3_ACTIVE\n";
        if(texMode3_.getValue() == TEXTURE_3D)
            headerSource += "#define VOLUME_3_3D\n";
        else
            headerSource += "#define VOLUME_3_2D\n";
    }
    if(inport4_.isReady()) {
        headerSource += "#define VOLUME_4_ACTIVE\n";
        if(texMode4_.getValue() == TEXTURE_3D)
            headerSource += "#define VOLUME_4_3D\n";
        else
            headerSource += "#define VOLUME_4_2D\n";
    }

    headerSource += "#define TF_SAMPLER_TYPE_1 " + transferFunc1_.get()->getSamplerType() + "\n";
    headerSource += "#define TF_SAMPLER_TYPE_2 " + transferFunc2_.get()->getSamplerType() + "\n";
    headerSource += "#define TF_SAMPLER_TYPE_3 " + transferFunc3_.get()->getSamplerType() + "\n";
    headerSource += "#define TF_SAMPLER_TYPE_4 " + transferFunc4_.get()->getSamplerType() + "\n";

    headerSource += "#define BLENDINGMODE1(result, color, factor) " + blendingMode1_.get() + "(result, color, factor);\n";
    headerSource += "#define BLENDINGMODE2(result, color, factor) " + blendingMode2_.get() + "(result, color, factor);\n";
    headerSource += "#define BLENDINGMODE3(result, color, factor) " + blendingMode3_.get() + "(result, color, factor);\n";
    headerSource += "#define BLENDINGMODE4(result, color, factor) " + blendingMode4_.get() + "(result, color, factor);\n";

    return headerSource;
}

} // namespace voreen
