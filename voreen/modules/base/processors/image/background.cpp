/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "background.h"
#include "voreen/core/voreenapplication.h"
#include "tgt/immediatemode/immediatemode.h"

#include "tgt/texturemanager.h"
#include "tgt/gpucapabilities.h"
#include "tgt/textureunit.h"

using tgt::TextureUnit;

namespace voreen {

namespace {

GLubyte* blur(GLubyte* image, int size) {
    GLubyte* n = new GLubyte[size * size];

    for (int y = 0; y < size; y++) {
        for (int x = 0; x < size; x++) {
            int pos = x + y * size;
            n[pos] = image[pos] / 5;

            if (x == 0)
                n[pos] += image[pos + (size - 1)] / 5;
            else
                n[pos] += image[pos - 1] / 5;

            if (y == 0)
                n[pos] += image[pos + (size - 1) * size] / 5;
            else
                n[pos] += image[pos - size] / 5;

            if (x == size-1)
                n[pos] += image[y * size] / 5;
            else
                n[pos] += image[pos + 1] / 5;

            if (y == size-1)
                n[pos] += image[x] / 5;
            else
                n[pos] +=image[pos + size] / 5;
        }
    }
    delete[] image;
    return n;
}

GLubyte* resize(GLubyte* image, int size) {
    GLubyte* n = new GLubyte[4 * size * size];

    for (int y = 0; y < size*2; y++)
        for (int x = 0; x < size*2; x++)
            n[x + y*size*2] = image[(x / 2) + (y / 2) * size];

    delete[] image;

    return blur(blur(n, size*2), size*2);
}

GLubyte* tile(GLubyte* image, int size) {
    GLubyte* n = new GLubyte[4 * size * size];

    for (int y = 0; y < size*2; y++)
        for (int x = 0; x < size*2; x++)
            n[x + y*size*2] = image[(x % size) + (y % size) * size];

    delete[] image;

    return n;
}

} // namespace

Background::Background()
    : ImageProcessorBypassable("image/background")
    , firstcolor_("color1", "First Color", tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f))
    , secondcolor_("color2", "Second Color", tgt::vec4(0.392f, 0.392f, 0.392f, 1.0f)) // 100/255
    , angle_("angle", "Angle", 60, 0, 359)
    , tex_(0)
    , textureLoaded_(false)
    , filename_("texture", "Texture", "Select texture",
                VoreenApplication::app()->getCoreResourcePath("textures"), "Image Files (*.jpg *.png *.bmp)",
                FileDialogProperty::OPEN_FILE, Processor::INVALID_RESULT)
    , tile_("repeat", "Repeat Background", 1.0f, 0.f, 100.f)
    , randomClouds_("randomClouds","Random?",true)
    , modeProp_("backgroundModeAsString", "Type")
    , blendMode_("blendMode", "Blend Mode", Processor::INVALID_PROGRAM)
    , inport_(Port::INPORT, "image.input", "Image Input")
    , outport_(Port::OUTPORT, "image.output", "Image Output")
    , privatePort_(Port::OUTPORT, "image.tmp", "image.tmp", false)
    , textureInvalid_(true)
{
    modeProp_.addOption("monochrome", "Monochrome");
    modeProp_.addOption("gradient",   "Gradient");
    modeProp_.addOption("radial",     "Radial");
    modeProp_.addOption("texture",    "Texture");
    modeProp_.addOption("cloud",      "Cloud");
    modeProp_.select("gradient");
    modeProp_.onChange(MemberFunctionCallback<Background>(this, &Background::onBackgroundModeChanged));
    onBackgroundModeChanged(); // setup initial property visibilty

    addProperty(modeProp_);
    addProperty(firstcolor_);
    addProperty(secondcolor_);
    firstcolor_.onChange(MemberFunctionCallback<Background>(this, &Background::invalidateTexture));
    secondcolor_.onChange(MemberFunctionCallback<Background>(this, &Background::invalidateTexture));

#ifdef VRN_MODULE_DEVIL
    addProperty(filename_);
    filename_.onChange(MemberFunctionCallback<Background>(this, &Background::invalidateTexture));
#endif
    addProperty(angle_);
    addProperty(tile_);
    addProperty(randomClouds_);

    blendMode_.addOption("alpha-blending",      "Alpha Blending (Default Setting)",   "BLEND_MODE_ALPHA");
    blendMode_.addOption("add",                 "Add (Use with Black Backgrounds)",     "BLEND_MODE_ADD");
    blendMode_.selectByKey("alpha-blending");
    addProperty(blendMode_);

    addPort(inport_);
    addPort(outport_);
    addPrivateRenderPort(&privatePort_);
}

Processor* Background::create() const {
    return new Background();
}

void Background::initialize() {
    ImageProcessorBypassable::initialize();
    loadTexture();
}

void Background::deinitialize() {
    if (tex_) {
        if (textureLoaded_)
            TexMgr.dispose(tex_);
        else
            delete tex_;
    }
    LGL_ERROR;
    tex_ = 0;

    ImageProcessorBypassable::deinitialize();
}

bool Background::isReady() const {
    return outport_.isReady();
}

void Background::process() {
    //bypass the processor
    if (!enableSwitch_.get()){
        bypass(&inport_, &outport_);
        return;
    }

    if(!shaderProp_.hasValidShader()) {
        return;
    }

    if (!inport_.isReady())
        outport_.activateTarget();
    else
        privatePort_.activateTarget();

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();

    glClearDepth(1.0);

    glDisable(GL_DEPTH_TEST);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // first the background
    renderBackground();

    if (!inport_.isReady())
        outport_.deactivateTarget();
    else
        privatePort_.deactivateTarget();

    glEnable(GL_DEPTH_TEST);

    // leave if there's nothing to render above
    if (inport_.isReady()) {
        outport_.activateTarget();

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // use the shader to to blend the original picture over the background and to keep the
        // depth values

        TextureUnit colorUnit0, depthUnit0, colorUnit1, depthUnit1;
        privatePort_.bindTextures(colorUnit0.getEnum(), depthUnit0.getEnum());
        inport_.bindTextures(colorUnit1.getEnum(), depthUnit1.getEnum());

        // initialize shader
        program_->activate();
        setGlobalShaderParameters(program_);
        program_->setUniform("colorTex0_", colorUnit0.getUnitNumber());
        program_->setUniform("depthTex0_", depthUnit0.getUnitNumber());
        program_->setUniform("colorTex1_", colorUnit1.getUnitNumber());
        program_->setUniform("depthTex1_", depthUnit1.getUnitNumber());
        privatePort_.setTextureParameters(program_, "textureParameters0_");
        inport_.setTextureParameters(program_, "textureParameters1_");

        renderQuad();

        outport_.deactivateTarget();

        program_->deactivate();
        TextureUnit::setZeroUnit();
    }
    LGL_ERROR;
}

std::string Background::generateHeader(const tgt::GpuCapabilities::GlVersion* version /*= 0*/) {
    std::string header = ImageProcessorBypassable::generateHeader(version);
    header += "#define " + blendMode_.getValue() + "\n";
    return header;
}

void Background::renderBackground() {
    if (modeProp_.isSelected("monochrome")) {
        glClearColor(firstcolor_.get().r, firstcolor_.get().g,
                    firstcolor_.get().b, firstcolor_.get().a);
        glClear(GL_COLOR_BUFFER_BIT);
        glClearColor(0.f,0.f,0.f,0.f);
        LGL_ERROR;
    }
    else if (modeProp_.isSelected("gradient")) {
        // linear gradient
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.pushMatrix();
        MatStack.loadIdentity();
        MatStack.rotate(static_cast<float>(angle_.get()), 0.0f, 0.0f, 1.0f);

        // when you rotate the texture, you need to scale it.
        // otherwise the edges don't cover the complete screen
        // magic number: 0.8284271247461900976033774484194f = sqrt(8)-2;
        MatStack.scale(1.0f + (45 - abs(angle_.get() % 90 - 45)) / 45.0f*0.8284271247461900976033774484194f,
            1.0f + (45 - abs(angle_.get() % 90 - 45)) / 45.0f*0.8284271247461900976033774484194f, 1.0f);


        IMode.begin(tgt::ImmediateMode::QUADS);
        IMode.color(firstcolor_.get());
        IMode.vertex(-1.f, -1.f);
        IMode.vertex( 1.f, -1.f);
        IMode.color(secondcolor_.get());
        IMode.vertex( 1.f, 1.f);
        IMode.vertex(-1.f, 1.f);
        IMode.end();
        IMode.color(1.f, 1.f, 1.f, 1.f);
        MatStack.popMatrix();
        LGL_ERROR;

    }
    else if (modeProp_.isSelected("radial")) {

        if (!tex_ || textureInvalid_)
            loadTexture();

        if (tex_) {
            tex_->bind();
            tgt::ImmediateMode::TexturMode texmode = static_cast<tgt::ImmediateMode::TexturMode>(tex_->getType());
            IMode.setTextureMode(texmode);

            MatStack.pushMatrix();
            MatStack.scale(1.44f, 1.44f, 0.0f);

            IMode.begin(tgt::ImmediateMode::QUADS);
            IMode.texcoord(0.0f,  0.0f);
            IMode.vertex( -1.0f, -1.0f);
            IMode.texcoord(tile_.get(), 0.0f);
            IMode.vertex( 1.0f, -1.0f);
            IMode.texcoord(tile_.get(), tile_.get());
            IMode.vertex( 1.0f, 1.0f);
            IMode.texcoord(0.0f, tile_.get());
            IMode.vertex( -1.0f, 1.0f);
            IMode.end();

            MatStack.popMatrix();

            IMode.setTextureMode(tgt::ImmediateMode::TEXNONE);
            LGL_ERROR;
        }
        else {
            LWARNING("No texture");
        }
    }
    else if (modeProp_.isSelected("texture")) {

        if (!textureLoaded_ || textureInvalid_)
            loadTexture();

        if (tex_) {
            tex_->bind();
            tgt::ImmediateMode::TexturMode texmode = static_cast<tgt::ImmediateMode::TexturMode>(tex_->getType());
            IMode.setTextureMode(texmode);

            IMode.begin(tgt::ImmediateMode::QUADS);
            IMode.texcoord(0.0f, 0.0f);
            IMode.vertex(-1.0f, -1.0f);
            IMode.texcoord(tile_.get(), 0.0f);
            IMode.vertex( 1.0f, -1.0f);
            IMode.texcoord(tile_.get(), tile_.get());
            IMode.vertex( 1.0f, 1.0f);
            IMode.texcoord(0.0f, tile_.get());
            IMode.vertex(-1.0f, 1.0f);
            IMode.end();

            IMode.setTextureMode(tgt::ImmediateMode::TEXNONE);
            LGL_ERROR;
        }
    }
    else if (modeProp_.isSelected("cloud")) {

        // cloud texture
        if (!tex_ || textureInvalid_)
            loadTexture();

        if (tex_) {
            tex_->bind();
            tgt::ImmediateMode::TexturMode texmode = static_cast<tgt::ImmediateMode::TexturMode>(tex_->getType());
            IMode.setTextureMode(texmode);

            IMode.begin(tgt::ImmediateMode::QUADS);
            IMode.texcoord(0.0f, 0.0f);
            IMode.vertex(-1.0f, -1.0f);
            IMode.texcoord(tile_.get(), 0.0f);
            IMode.vertex( 1.0f, -1.0f);
            IMode.texcoord(tile_.get(), tile_.get());
            IMode.vertex( 1.0f, 1.0f);
            IMode.texcoord(0.0f, tile_.get());
            IMode.vertex(-1.0f, 1.0f);
            IMode.end();

            IMode.setTextureMode(tgt::ImmediateMode::TEXNONE);
            LGL_ERROR;
        }
        else {
            LWARNING("No texture");
        }
    }
    else {
        LWARNING("Unknown background mode");
    }
}

void Background::onBackgroundModeChanged() {

    if (modeProp_.isSelected("monochrome")) {
        firstcolor_.setVisibleFlag(true);
        secondcolor_.setVisibleFlag(false);
        angle_.setVisibleFlag(false);
        filename_.setVisibleFlag(false);
        tile_.setVisibleFlag(false);
        randomClouds_.setVisibleFlag(false);
    }
    else if (modeProp_.isSelected("gradient")) {
        firstcolor_.setVisibleFlag(true);
        secondcolor_.setVisibleFlag(true);
        angle_.setVisibleFlag(true);
        filename_.setVisibleFlag(false);
        tile_.setVisibleFlag(false);
        randomClouds_.setVisibleFlag(false);
    }
    else if (modeProp_.isSelected("radial")) {
        tile_.set(1.f);
        firstcolor_.setVisibleFlag(true);
        secondcolor_.setVisibleFlag(true);
        angle_.setVisibleFlag(false);
        filename_.setVisibleFlag(false);
        tile_.setVisibleFlag(false);
        randomClouds_.setVisibleFlag(false);
    }
    else if (modeProp_.isSelected("cloud")) {
        firstcolor_.setVisibleFlag(false);
        secondcolor_.setVisibleFlag(false);
        angle_.setVisibleFlag(false);
        filename_.setVisibleFlag(false);
        tile_.setVisibleFlag(true);
        randomClouds_.setVisibleFlag(true);
    }
    else if (modeProp_.isSelected("texture")) {
        firstcolor_.setVisibleFlag(false);
        secondcolor_.setVisibleFlag(false);
        angle_.setVisibleFlag(false);
        filename_.setVisibleFlag(true);
        tile_.setVisibleFlag(true);
        randomClouds_.setVisibleFlag(false);
    }

    textureInvalid_ = true;
    invalidate();
}

void Background::loadTexture() {
    // only use OpenGL after processor is initialized
    if (!tgt::TextureManager::isInited())
        return;

    // is a texture already loaded? -> then delete
    if (tex_) {
        if (textureLoaded_) {
            TexMgr.dispose(tex_);
            textureLoaded_ = false;
        }
        else
            delete tex_;

        LGL_ERROR;
        tex_ = 0;
    }

    // create Texture
    if (modeProp_.isSelected("radial")) {
        createRadialTexture();
    }
    else if (modeProp_.isSelected("texture")) {
        if (!filename_.get().empty()) {
            tex_ = TexMgr.load(filename_.get(), tgt::Texture::LINEAR, false, false, true, true);
        }
        if (tex_)
            textureLoaded_ = true;
    }
    else if (modeProp_.isSelected("cloud")) {
        createCloudTexture();
    }
    textureInvalid_ = false;

    // adapt texture's data type to input texture, if necessary
    /*tgt::Texture* portTex = inport_.hasData() ? inport_.getData()->getColorTexture() : 0;
    if (tex_ && portTex && (tex_->getDataType() != portTex->getDataType()) ) {
        delete[] tex_->getPixelData();
        tex_->setPixelData(0);
        GLubyte* convBuffer = tex_->downloadTextureToBuffer(tex_->getFormat(), portTex->getDataType());
        LGL_ERROR;
        tex_->setPixelData(convBuffer);
        tex_->setDataType(portTex->getDataType());
        tex_->setInternalFormat(portTex->getInternalFormat());
        tex_->uploadTexture();
        LGL_ERROR;
    }*/
}

void Background::createCloudTexture() {
    int imgx = 256;
    int imgy = 256;
    int size = 32;

    //random?
    if(!randomClouds_.get())
        srand(0);

    // creates 4 blured and resized noise-datas
    GLubyte* tex1 = new GLubyte[size * size];

    for (int k = 0; k < size*size; k++)
        tex1[k] = rand() % 256;

    tex1 = tile(tile(tile(blur(tex1, size), size), size * 2), size * 4);
    GLubyte* tex2 = new GLubyte[size * size];

    for (int k = 0; k < size*size; k++)
        tex2[k] = rand() % 256;

    tex2 = tile(tile(resize(blur(tex2, size), size), size * 2), size * 4);
    GLubyte* tex3 = new GLubyte[size * size];

    for (int k = 0; k < size*size; k++) {
        tex3[k] = rand() % 256;
    }

    tex3 = tile(resize(resize(blur(tex3, size), size), size * 2), size * 4);
    GLubyte* tex4 = new GLubyte[size * size];

    for (int k = 0; k < size*size; k++)
        tex4[k] = rand() % 256;

    tex4 = resize(resize(resize(blur(tex4, size), size), size * 2), size * 4);

    // creates the final cloud data with the 4 textures from above
    GLubyte* tex5 = new GLubyte[imgx * imgy];

    for (int y = 0; y < imgy; y++) {
        for (int x = 0; x < imgx; x++) {
            int pos = x + imgx * y;
            // mix textures and don't use all values, so that there are free areas in the sky
            float c = (tex4[pos] + tex3[pos] / 2 + tex2[pos] / 4 + tex1[pos] / 8) / 2 - 110.0f;
            if (c < 0)
                c = 0;
            // exp-function to calc the final value for a nicer look
            tex5[pos] = (GLubyte)(255 - (pow(0.95f, c) * 255));
        }
    }

    tex5 = blur(tex5, imgx);

    tex_ = new tgt::Texture(tgt::ivec3(imgx, imgy, 1));
    tex_->alloc(true);

    for (int y = 0; y < imgy; y++) {
        for (int x = 0; x < imgx; x++) {
            tgt::vec4 cloudCol = tgt::vec4(1.f);
            tgt::vec4 backgroundCol = tgt::vec4(0.5f, 0.5f, 1.f, 1.f);
            tgt::vec4 mixColor = tgt::mix(backgroundCol, cloudCol, tgt::clamp(tex5[x + imgx * y] / 255.f, 0.f, 1.f));
            tex_->texel<tgt::col4>(x,y) = tgt::col4(tgt::iround(mixColor*255.f));
        }
    }

    delete[] tex1;
    delete[] tex2;
    delete[] tex3;
    delete[] tex4;
    delete[] tex5;

    tex_->uploadTexture();
    LGL_ERROR;
}

void Background::createRadialTexture() {
    int imgx = 256;
    int imgy = 256;
    tex_ = new tgt::Texture(tgt::ivec3(imgx, imgy, 1));
    tex_->alloc(true);

    // create pixel data
    float r;
    for (int y = 0; y < imgy; y++) {
        for (int x = 0; x < imgx; x++) {
            // calculate radius (Pythagoras)
            r = static_cast<float>((x - imgx / 2) * (x - imgx / 2) + (y - imgy / 2) * (y - imgy / 2));
            r = sqrt(r);
            // norm to half size of smaller value
            r = 2 * r / ((imgx < imgy) ? imgx : imgy);
            // mix both colors according to radius
            r = tgt::clamp(r, 0.f, 1.f);
            tgt::vec4 mixedColor = tgt::mix(firstcolor_.get(), secondcolor_.get(), r);
            tex_->texel<tgt::col4>(x,y) = tgt::col4(tgt::iround(mixedColor * 255.f));;
        }
    }

    tex_->uploadTexture();
    LGL_ERROR;
}

void Background::createEmptyTexture() {
    int imgx = 256;
    int imgy = 256;
    tex_ = new tgt::Texture(tgt::ivec3(imgx, imgy, 1));
    GLubyte* buffer = tex_->alloc(true);

    memset(buffer,0,tex_->getSizeOnCPU());

    tex_->uploadTexture();
    LGL_ERROR;
}

void Background::invalidateTexture() {
    textureInvalid_ = true;
}

} // namespace
