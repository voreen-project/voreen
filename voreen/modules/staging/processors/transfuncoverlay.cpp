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

#include "transfuncoverlay.h"

#include "tgt/textureunit.h"
#include "tgt/immediatemode/immediatemode.h"

#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"
#include "voreen/core/datastructures/transfunc/1d/preintegrationtable.h"

#include "voreen/core/utils/stringutils.h"

#include <sstream>

using tgt::TextureUnit;

namespace voreen {

TransFuncOverlay::TransFuncOverlay()
    : ImageProcessor("image/compositor")
    , imageInport_(Port::INPORT, "image.in", "Image Inpput")
    , privatePort_(Port::OUTPORT, "private", "Privateport")
    , outport_(Port::OUTPORT, "image.out", "Image Output")
    , fontProp_("voreen.fontprop", "Font:")
    , transferFunc_("transferFunction", "Transfer Function:")
    , renderPreIntegrationTable_("renderPreIntegrationTable", "Render Pre-Integration Table", false)
    , renderOverlay_("renderOverlay", "Render Overlay", true)
    , usePixelCoordinates_("usePixelCoordinates", "Move/Resize Overlay by ")
    , overlayBottomLeft_("overlayBottomLeft", "Overlay Bottom Left", tgt::ivec2(10), tgt::ivec2(-4096), tgt::ivec2(4096))
    , overlayDimensions_("overlayDimensions", "Overlay Dimensions", tgt::ivec2(100), tgt::ivec2(0), tgt::ivec2(4096))
    , overlayBottomLeftRelative_("overlayBottomLeftRelative", "Overlay Bottom Left (Relative)", tgt::vec2(0.05f), tgt::vec2(-1.f), tgt::vec2(1.f))
    , overlayDimensionsRelative_("overlayDimensionsRelative", "Overlay Dimensions (Relative)", tgt::vec2(0.25f), tgt::vec2(0.f), tgt::vec2(1.f))
    , overlayOpacity_("overlayOpacity", "Overlay Opacity", 1.f)
    , fontColor_("fontColor", "Font Color", tgt::Color(1.f, 1.f, 1.f, 1.f))
    , renderRangeProp_("renderRangeProp","Render Range", true)
    , valueRangeProp_("valueRangeProp", "Value Range (Thresholds)", tgt::vec2(0.f, 1.f), tgt::vec2(-10000.f), tgt::vec2(10000.f))
    , deriveRangeFromTransFunc_("deriveRangeFromTransFunc", "Derive range from transfer function", false)
    , rangePrecisionProp_("rangePrecisionProp", "Range Precision")
    , tfUnit_("tfUnit", "Unit (max 20 Chars)", "")
    , renderBorder_("renderBorder", "Render Border", true)
    , borderWidth_("borderWidth", "Border Width", 2.f, 0.1f, 5.f)
    , borderColor_("borderColor", "Border Color", tgt::Color(0.f, 0.f, 0.f, 1.f))
    , copyShader_(0)
{
    addPort(imageInport_);
    addPrivateRenderPort(&privatePort_);
    addPort(outport_);

    addProperty(fontProp_);
    fontProp_.setVisibleFlag(false);
    addProperty(transferFunc_);
    addProperty(renderPreIntegrationTable_);

    addProperty(renderOverlay_);
    renderOverlay_.setGroupID("general");
    addProperty(overlayOpacity_);
    overlayOpacity_.setGroupID("general");
    addProperty(usePixelCoordinates_);
    usePixelCoordinates_.setGroupID("general");
    usePixelCoordinates_.addOption("true", "Pixel Coordinates", true);
    usePixelCoordinates_.addOption("false", "Relative Position", false);
    usePixelCoordinates_.select("true");
    ON_CHANGE(usePixelCoordinates_, TransFuncOverlay, onChangeUsePixelCoordinates);
    setPropertyGroupGuiName("general", "General");

    addProperty(overlayBottomLeft_);
    overlayBottomLeft_.setGroupID("overlayAbs");
    addProperty(overlayDimensions_);
    overlayDimensions_.setGroupID("overlayAbs");

    addProperty(overlayBottomLeftRelative_);
    overlayBottomLeftRelative_.setGroupID("overlayPer");
    addProperty(overlayDimensionsRelative_);
    overlayDimensionsRelative_.setGroupID("overlayPer");

    setPropertyGroupGuiName("overlayAbs", "Overlay Options (absolute)");
    setPropertyGroupGuiName("overlayPer", "Overlay Options (relative)");

    addProperty(fontColor_);
    fontColor_.setGroupID("overlay settings");
    addProperty(tfUnit_);
    tfUnit_.setGroupID("overlay settings");
    addProperty(renderRangeProp_);
    renderRangeProp_.setGroupID("overlay settings");
    addProperty(deriveRangeFromTransFunc_);
    deriveRangeFromTransFunc_.setGroupID("overlay settings");
    addProperty(valueRangeProp_);
    valueRangeProp_.setGroupID("overlay settings");
    addProperty(rangePrecisionProp_);
        rangePrecisionProp_.addOption("default","Default",-1);
        rangePrecisionProp_.addOption("zero","Zero", 0);
        rangePrecisionProp_.addOption("one","One", 1);
        rangePrecisionProp_.addOption("two","Two", 2);
        rangePrecisionProp_.setGroupID("overlay settings");
    addProperty(renderBorder_);
    renderBorder_.setGroupID("overlay settings");
    addProperty(borderWidth_);
    borderWidth_.setGroupID("overlay settings");
    addProperty(borderColor_);
    borderColor_.setGroupID("overlay settings");
    setPropertyGroupGuiName("overlay settings", "Overlay Settings");

    onChangeUsePixelCoordinates();
    ON_CHANGE(deriveRangeFromTransFunc_, TransFuncOverlay, onChangeDeriveRangeFromTransFunc);
}

Processor* TransFuncOverlay::create() const {
    return new TransFuncOverlay();
}

void TransFuncOverlay::initialize() {
    ImageProcessor::initialize();

    copyShader_ = ShdrMgr.loadSeparate("passthrough.vert", "copyimage.frag",
        ImageProcessor::generateHeader(), false);
}

void TransFuncOverlay::deinitialize() {
    ShdrMgr.dispose(copyShader_);
    copyShader_ = 0;

    ImageProcessor::deinitialize();
}

bool TransFuncOverlay::isReady() const {
    return (imageInport_.isReady() && outport_.isReady());
}

std::string TransFuncOverlay::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string header = ImageProcessor::generateHeader(version);
    header += "#define MODE_ALPHA_BLENDING_B_OVER_A\n";
    //header += "#define MODE_WEIGHTED_AVERAGE\n";
    return header;
}

void TransFuncOverlay::beforeProcess() {
    ImageProcessor::beforeProcess();
    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM)
        compile();
}

void TransFuncOverlay::process() {

    tgtAssert(outport_.isReady(), "Outport not ready");
    tgtAssert(imageInport_.isReady(), "Inport not ready");
    tgtAssert(program_ && copyShader_, "Shader missing");

    //same code as in ImageOverlay
    //          |
    //          v
    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // bind input image to tex unit
    TextureUnit imageUnit, imageUnitDepth;
    imageInport_.bindTextures(imageUnit.getEnum(), imageUnitDepth.getEnum());

    // 1. copy input image to outport
    copyShader_->activate();
    setGlobalShaderParameters(copyShader_);
    imageInport_.setTextureParameters(copyShader_, "texParams_");
    copyShader_->setUniform("useTexcoords", false);
    copyShader_->setUniform("colorTex_", imageUnit.getUnitNumber());
    copyShader_->setUniform("depthTex_", imageUnitDepth.getUnitNumber());
    renderQuad();
    copyShader_->deactivate();
    outport_.deactivateTarget();
    LGL_ERROR;


    tgt::ivec2 overlayDim, overlayBL;
    getCurrentOverlayPixelCoordinates(overlayBL, overlayDim);

    //render overlay
    // Slightly hackish: Resize the private port every process, because Renderprocessor
    // overwrites the size of all private ports, when a resize event happens.
    privatePort_.resize(overlayDim);
    privatePort_.activateTarget();
    //glPushAttrib(GL_ALL_ATTRIB_BITS);
    glClearColor(fontColor_.get().r, fontColor_.get().g, fontColor_.get().b, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //render the transfer function texture
    const tgt::Texture* tfTex = 0;
    if (renderPreIntegrationTable_.get())
        tfTex = transferFunc_.get()->getPreIntegrationTable(1.0f / (41.0f * 2.0f))->getTexture();
    else
        tfTex = transferFunc_.get()->getTexture();

    tgtAssert(tfTex, "No transfer function texture");

    IMode.color(1.f, 1.f, 1.f, 1.f);
    glDisable(GL_DEPTH_TEST);
    IMode.begin(tgt::ImmediateMode::QUADS);
        IMode.vertex(-0.8f, -0.9f);
        IMode.vertex(-0.5f, -0.9f);
        IMode.vertex(-0.5f, 0.7f);
        IMode.vertex(-0.8f, 0.7f);
    IMode.end();
    IMode.color(0.8f, 0.8f, 0.8f, 1.f);
    IMode.begin(tgt::ImmediateMode::QUADS);
        IMode.vertex(-0.8f, -0.9f);
        IMode.vertex(-0.65f, -0.9f);
        IMode.vertex(-0.65f, -0.58f);
        IMode.vertex(-0.8f, -0.58f);

        IMode.vertex(-0.65f, -0.58f);
        IMode.vertex(-0.5f, -0.58f);
        IMode.vertex(-0.5f, -0.26f);
        IMode.vertex(-0.65f, -0.26f);

        IMode.vertex(-0.8f, -0.26f);
        IMode.vertex(-0.65f, -0.26f);
        IMode.vertex(-0.65f, 0.06f);
        IMode.vertex(-0.8f, 0.06f);

        IMode.vertex(-0.65f, 0.06f);
        IMode.vertex(-0.5f, 0.06f);
        IMode.vertex(-0.5f, 0.38f);
        IMode.vertex(-0.65f, 0.38f);

        IMode.vertex(-0.8f, 0.38f);
        IMode.vertex(-0.65f, 0.38f);
        IMode.vertex(-0.65f, 0.7f);
        IMode.vertex(-0.8f, 0.7f);
    IMode.end();
    IMode.color(1.f, 1.f, 1.f, 1.f);

    if (renderPreIntegrationTable_.get()){
        IMode.setTextureMode(tgt::ImmediateMode::TexturMode::TEX2D);
    }
    else{
        IMode.setTextureMode(tgt::ImmediateMode::TexturMode::TEX1D);
    }
    tfTex->bind();
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    IMode.begin(tgt::ImmediateMode::QUADS);
    if (renderPreIntegrationTable_.get()) {
        IMode.texcoord(0.f, 0.f); IMode.vertex(-0.8f, -0.9f);
        IMode.texcoord(1.f, 0.f); IMode.vertex(-0.5f, -0.9f);
        IMode.texcoord(1.f, 1.f); IMode.vertex(-0.5f, 0.7f);
        IMode.texcoord(0.f, 1.f); IMode.vertex(-0.8f, 0.7f);
    }
    else {
        IMode.texcoord(0.f); IMode.vertex(-0.8f, -0.9f);
        IMode.texcoord(0.f); IMode.vertex(-0.5f, -0.9f);
        IMode.texcoord(1.f); IMode.vertex(-0.5f, 0.7f);
        IMode.texcoord(1.f); IMode.vertex(-0.8f, 0.7f);
    }
    IMode.end();
    //glBlendFuncSeparate(GL_ONE,GL_ZERO,GL_ONE,GL_ZERO);
    glBlendFunc(GL_ONE, GL_ZERO);
    glDisable(GL_BLEND);
    //glEnable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    //render fonts
    tgt::ivec2 privatePortSize = privatePort_.getSize();
    tgt::vec2 range;
    if (deriveRangeFromTransFunc_.get()){
        range = transferFunc_.get()->getDomain();
    }
    else{
        range = valueRangeProp_.get();
    }

    int fontSize = privatePort_.getSize().y / 8;
    fontProp_.get()->setFontColor(fontColor_.get()*tgt::vec4(1.0f, 1.0f, 1.0f, overlayOpacity_.get()));
    fontProp_.get()->setFontSize(fontSize);
    fontProp_.get()->setTextAlignment(tgt::Font::MiddleLeft);
    fontProp_.get()->setLineWidth(0.0f);
    tgt::vec2 textSize = fontProp_.get()->getSize(tgt::vec3(0, privatePort_.getSize().y*0.925f, 0), tfUnit_.get(), privatePortSize);
    float offset = 0;
    if (textSize.x < privatePort_.getSize().x*0.35f){
        offset = (privatePort_.getSize().x*0.35f - textSize.x)*0.5f;
    }
    // render unit always
    fontProp_.get()->render(tgt::vec3(offset, privatePort_.getSize().y*0.925f, 0), tfUnit_.get(), privatePortSize);
    fontProp_.get()->setLineWidth(0.0f);
    fontProp_.get()->setTextAlignment(tgt::Font::MiddleLeft);

    //render range only if activated
    if(renderRangeProp_.get()) {
        std::stringstream strstr;
        if(rangePrecisionProp_.getValue() == -1)
            strstr << range.x; //better formatting than ftos...
        else
            strstr << ftos(range.x, rangePrecisionProp_.getValue());
        fontProp_.get()->render(tgt::vec3(privatePort_.getSize().x*0.3f, privatePort_.getSize().y*0.10f, 0), strstr.str(), privatePortSize);
        strstr.clear();
        strstr.str("");
        if(rangePrecisionProp_.getValue() == -1)
            strstr << (range.y - range.x) / 2.f + range.x; //better formatting than ftos...
        else
            strstr << ftos((range.y - range.x) / 2.f + range.x, rangePrecisionProp_.getValue());
        fontProp_.get()->render(tgt::vec3(privatePort_.getSize().x*0.3f, privatePort_.getSize().y*0.45f, 0), strstr.str(), privatePortSize);
        strstr.clear();
        strstr.str("");
        if(rangePrecisionProp_.getValue() == -1)
            strstr << range.y; //better formatting than ftos...
        else
            strstr << ftos(range.y, rangePrecisionProp_.getValue());
        fontProp_.get()->render(tgt::vec3(privatePort_.getSize().x*0.3f, privatePort_.getSize().y*0.80f, 0), strstr.str(), privatePortSize);
    }


    // render border around overlay
    if (renderBorder_.get()) {
        //glPushAttrib(GL_ALL_ATTRIB_BITS);
        IMode.color(borderColor_.get().r, borderColor_.get().g, borderColor_.get().b, borderColor_.get().a*overlayOpacity_.get());
        glEnable(GL_LINE_SMOOTH);
        glLineWidth(borderWidth_.get());
        glDepthFunc(GL_ALWAYS);
        IMode.setTextureMode(tgt::ImmediateMode::TexturMode::TEXNONE);
        IMode.begin(tgt::ImmediateMode::LINE_STRIP);
            IMode.vertex(-0.8f, -0.9f);
            IMode.vertex(-0.5f, -0.9f);
            IMode.vertex(-0.5f, 0.7f);
            IMode.vertex(-0.8f, 0.7f);
            IMode.vertex(-0.8f, -0.9f);
        IMode.end();
        glDepthFunc(GL_LESS);
        glDisable(GL_LINE_SMOOTH);
    }
    LGL_ERROR;
    privatePort_.deactivateTarget();

    // 2. render overlay over copied input image (using compositor shader)
    // check, if overlay dims are greater zero
    bool dimensionsValid = ((usePixelCoordinates_.getValue() && tgt::hand(tgt::greaterThan(overlayDimensions_.get(), tgt::ivec2(0)))) ||
        (!usePixelCoordinates_.getValue() && tgt::hand(tgt::greaterThan(overlayDimensionsRelative_.get(), tgt::vec2(0.f)))));

    outport_.activateTarget();
    if (renderOverlay_.get() && /*overlayInport_.isReady() &&*/ dimensionsValid) {
        // bind overlay to tex unit
        TextureUnit overlayUnit;
        tgt::Texture* overlayTex = privatePort_.getColorTexture();//overlayInport_.getColorTexture();
        tgtAssert(overlayTex, "No overlay texture");
        overlayUnit.activate();
        overlayTex->bind();

        program_->activate();
        setGlobalShaderParameters(program_);

        // image texture parameters
        imageInport_.setTextureParameters(program_, "textureParameters0_");
        program_->setUniform("colorTex0_", imageUnit.getUnitNumber());
        program_->setUniform("depthTex0_", imageUnitDepth.getUnitNumber());
        program_->setUniform("colorTex1_", overlayUnit.getUnitNumber());
        program_->setUniform("depthTex1_", 0);
        //program_->setUniform("weightingFactor_", 1.f-overlayOpacity_.get());

        tgt::vec2 outportDim = tgt::vec2(outport_.getSize());


        // overlay texture matrix mapping from normalized frag coords (outport) to overlay tex coords
        tgt::mat4 overlayTexCoordMatrix = tgt::mat4::identity;
        overlayTexCoordMatrix *= tgt::mat4::createScale(tgt::vec3(outportDim / tgt::vec2(overlayDim), 0.f));
        overlayTexCoordMatrix *= tgt::mat4::createTranslation(-tgt::vec3(tgt::vec2(overlayBL) / outportDim, 1.f));

        // overlay texture parameters
        bool oldIgnoreError = program_->getIgnoreUniformLocationError();
        program_->setIgnoreUniformLocationError(true);
        program_->setUniform("textureParameters1_.dimensions_", tgt::vec2(overlayDim));
        program_->setUniform("textureParameters1_.dimensionsRCP_", tgt::vec2(1.f) / tgt::vec2(overlayDim));
        program_->setUniform("textureParameters1_.matrix_", overlayTexCoordMatrix);
        program_->setIgnoreUniformLocationError(oldIgnoreError);
        LGL_ERROR;

        // render overlay at specified position and size
        tgt::vec2 bl = 2.f*tgt::vec2(overlayBL) / outportDim - 1.f;
        tgt::vec2 dim = 2.f*tgt::vec2(overlayDim) / outportDim;
        glDepthFunc(GL_ALWAYS);
        IMode.begin(tgt::ImmediateMode::QUADS);
            IMode.vertex(bl.x, bl.y);
            IMode.vertex(bl.x + dim.x, bl.y);
            IMode.vertex(bl.x + dim.x, bl.y + dim.y);
            IMode.vertex(bl.x, bl.y + dim.y);
        IMode.end();
        glDepthFunc(GL_LESS);
        program_->deactivate();
        LGL_ERROR;
    }
    IMode.color(tgt::vec4::one);
    IMode.setTextureMode(tgt::ImmediateMode::TEXNONE);
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    outport_.deactivateTarget();
    glLineWidth(1.0f);
    TextureUnit::setZeroUnit();
    LGL_ERROR;
}

void TransFuncOverlay::onChangeUsePixelCoordinates(){
    if (usePixelCoordinates_.getValue()){
        setPropertyGroupVisible("overlayAbs", true);
        setPropertyGroupVisible("overlayPer", false);
    }
    else {
        setPropertyGroupVisible("overlayAbs", false);
        setPropertyGroupVisible("overlayPer", true);
    }
}

void TransFuncOverlay::onChangeDeriveRangeFromTransFunc(){
    valueRangeProp_.setVisibleFlag(!deriveRangeFromTransFunc_.get());
}

void TransFuncOverlay::getCurrentOverlayPixelCoordinates(tgt::ivec2& bl, tgt::ivec2& dim) const {
    tgt::ivec2 outportDim = outport_.getSize();
    if (usePixelCoordinates_.getValue()) {
        bl = overlayBottomLeft_.get();
        dim = overlayDimensions_.get();
    }
    else {
        bl = tgt::ivec2(tgt::round(overlayBottomLeftRelative_.get() * tgt::vec2(outportDim)));
        dim = tgt::ivec2(tgt::round(overlayDimensionsRelative_.get() * tgt::vec2(outportDim)));
    }
    dim = tgt::max(dim, tgt::ivec2::two /* one => TEXTURE_1D according to tgt::Texture */);
}
} // namespace voreen
