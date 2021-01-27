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

#include "voreen/qt/mainwindow/menuentities/rendertargetviewermenuentity/rendertargetviewerpainter.h"
#include "voreen/qt/mainwindow/menuentities/rendertargetviewermenuentity/rendertargetviewer.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/glsl.h"
#include "voreen/core/ports/renderport.h"
#include "voreen/core/datastructures/rendertarget/rendertarget.h"

#include "tgt/glcanvas.h"
#include "tgt/framebufferobject.h"
#include "tgt/exception.h"
#include "tgt/shadermanager.h"
#include "tgt/font.h"
#include "tgt/immediatemode/immediatemode.h"

namespace {

    int fontSize = 12;

tgt::vec3 rgbToHsv(tgt::vec3 rgb) {
    tgt::vec3 hsv;
    float min = std::min(rgb.r, std::min(rgb.g, rgb.b));
    float max = std::max(rgb.r, std::max(rgb.g, rgb.b));
    float diff = max - min;
    if(max == min)
        hsv[0] = 0;
    else if(max == rgb.r)
        hsv[0] = 60 * (0 + rgb.g - rgb.b / diff);
    else if(max == rgb.g)
        hsv[0] = 60 * (2 + rgb.b - rgb.r / diff);
    else if(max == rgb.b)
        hsv[0] = 60 * (4 + rgb.r - rgb.g / diff);
    if(hsv[0] < 0)
        hsv[0] += 360;
    if(max == 0)
        hsv[1] = 0;
    else
        hsv[1] = diff / max;
    hsv[2] = max;
    return hsv;
}

QString internalTypeToString(GLint internalType) {
    switch (internalType) {
        case GL_RGB: return "GL_RGB";
        case GL_RGBA: return "GL_RGBA";
        case GL_RGBA16: return "GL_RGBA16";
        case GL_R32F: return "GL_R32F";
        case GL_R32UI: return "GL_R32UI";
#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
        case GL_RGB16F_ARB: return "GL_RGB16F_ARB";
        case GL_RGBA16F_ARB: return "GL_RGBA16F_ARB";
        case GL_RGBA32F_ARB: return "GL_RGBA32F_ARB";
#else
        case GL_RGB16F: return "GL_RGB16F";
        case GL_RGBA16F: return "GL_RGBA16F";
        case GL_RGBA32F: return "GL_RGBA32F";
#endif
        case GL_DEPTH_COMPONENT16: return "GL_DEPTH_COMPONENT16";
        case GL_DEPTH_COMPONENT24: return "GL_DEPTH_COMPONENT24";
        case GL_DEPTH_COMPONENT32: return "GL_DEPTH_COMPONENT32";
        default:
            return "Unknown type";
    }
}

}

namespace voreen {

const std::string RenderTargetViewerPainter::loggerCat_ = "voreen.RenderTargetViewerPainter";

RenderTargetViewerPainter::RenderTargetViewerPainter(tgt::GLCanvas* canvas, RenderTargetViewer* viewer)
    : Painter(canvas)
    , renderTargetViewer_(viewer)
    , colorProgram_(0)
    , inversecolorProgram_(0)
    , fbo_(0)
    , colorTex_(0)
    , fontTex_(0)
    , quadVbo_(0)
    , quadVao_(0)
    , font_(0)
{
    tgtAssert(viewer, "no render target viewer");
}

RenderTargetViewerPainter::~RenderTargetViewerPainter() {
    delete fbo_;
    delete font_;
    delete colorTex_;
    delete fontTex_;
    ShdrMgr.dispose(colorProgram_);
    ShdrMgr.dispose(inversecolorProgram_);

    glDeleteBuffers(1, &quadVbo_);
    glDeleteVertexArrays(1, &quadVao_);
}

void RenderTargetViewerPainter::sizeChanged(const tgt::ivec2& size) {
    glViewport(0,0, size.x, size.y);
}

void RenderTargetViewerPainter::initialize(){
    Painter::initialize();
    tgtAssert(canvas_, "no canvas");
    glViewport(0,0, canvas_->getPhysicalSize().x, canvas_->getPhysicalSize().y);

    try {
        colorProgram_ = ShdrMgr.loadSeparate("passthrough.vert", "rendertargetviewer/color.frag", GLSL::generateStandardShaderHeader(), false);
        inversecolorProgram_ = ShdrMgr.loadSeparate("passthrough.vert", "rendertargetviewer/inversecolor.frag", GLSL::generateStandardShaderHeader(), false);
    }
    catch(tgt::Exception& e) {
        LERROR(e.what());
    }

    fbo_ = new tgt::FramebufferObject();

    tgt::ivec3 size(canvas_->getPhysicalSize().x, canvas_->getPhysicalSize().y, 1);

    colorTex_ = new tgt::Texture(size, GL_RGBA, GL_RGBA16, GL_UNSIGNED_SHORT, tgt::Texture::LINEAR, tgt::Texture::CLAMP_TO_EDGE);
    colorTex_->uploadTexture();

    fontTex_ = new tgt::Texture(size, GL_RGBA, GL_RGBA, GL_UNSIGNED_SHORT, tgt::Texture::LINEAR, tgt::Texture::CLAMP_TO_EDGE);
    fontTex_->uploadTexture();

    glGenBuffers(1, &quadVbo_);
    glGenVertexArrays(1, &quadVao_);

    font_ = new tgt::Font(VoreenApplication::app()->getFontPath("VeraMono.ttf"), fontSize);
}

void RenderTargetViewerPainter::paint() {
    glDisable(GL_DEPTH_TEST);

    glClearColor(0.7f, 0.7f, 0.7f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);

    LGL_ERROR;

    //return if we have no network
    if (!renderTargetViewer_->evaluator_) {
        glEnable(GL_DEPTH_TEST);
        glClearColor(0.f,0.f,0.f,0.f);
        glColor4f(1.f, 1.f, 1.f, 1.f);
        return;
    }

    std::vector<RenderPort*> renderPorts = renderTargetViewer_->collectRenderPorts();

    // update window title
    int memsize = 0;
    for (size_t i = 0; i < renderPorts.size(); ++i) {
        RenderTarget* rt = renderPorts[i]->getRenderTarget();
        if (rt->getColorTexture())
            memsize += rt->getColorTexture()->getSizeOnGPU();
        if (rt->getDepthTexture())
            memsize += rt->getDepthTexture()->getSizeOnGPU();
    }
    memsize /= 1024 * 1024;
    QString title = QObject::tr("Render Target Viewer: %1 Render Targets (%2 mb)").arg(renderPorts.size()).arg(memsize);
    if (renderTargetViewer_->parentWidget())
        renderTargetViewer_->parentWidget()->setWindowTitle(title);

    // render port contents
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.multMatrix(tgt::mat4::createOrtho(0,canvas_->getPhysicalSize().x,0,canvas_->getPhysicalSize().y,-1,1));

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.loadIdentity();

    //return if we have no render ports
    if (renderPorts.empty()) {
        renderFont(13, 13, "No rendertargets available.");
        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.popMatrix();
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.popMatrix();
        glEnable(GL_DEPTH_TEST);
        glClearColor(0.f,0.f,0.f,0.f);
        glColor4f(1.f, 1.f, 1.f, 1.f);
        return;
    }

    //render on large port, or all ports
    if (renderTargetViewer_->maximizeOnePort_) {
        if(renderTargetViewer_->selectedRenderPortIndex_ >= 0 && renderTargetViewer_->selectedRenderPortIndex_ < (int)renderPorts.size()) {
            paintPort(renderPorts[renderTargetViewer_->selectedRenderPortIndex_], renderTargetViewer_->selectedRenderPortIndex_,
                      canvas_->getPhysicalSize().x, canvas_->getPhysicalSize().y);
        }
    }
    else {
        // update port grid
        int countX = static_cast<int>(std::ceil(std::sqrt(static_cast<float>(renderPorts.size()))));
        int countY = static_cast<int>(std::ceil(static_cast<float>(renderPorts.size()) / countX));

        int pixelPerPortWidth = static_cast<int>(std::ceil(static_cast<float>(canvas_->getPhysicalSize().x) / countX));
        int pixelPerPortHeight = static_cast<int>(std::ceil(static_cast<float>(canvas_->getPhysicalSize().y) / countY));

        for (int y = 0; y < countY; ++y) {
            for (int x = 0; x < countX; ++x) {
                int index = (countX*y)+x;
                if (index >= static_cast<int>(renderPorts.size()))
                    break;

                MatStack.pushMatrix();
                MatStack.translate(pixelPerPortWidth * x, pixelPerPortHeight * (countY - 1 - y), 0.0);
                paintPort(renderPorts[index], static_cast<int>(index), pixelPerPortWidth, pixelPerPortHeight);
                MatStack.popMatrix();
            }
        }
    }

    //clean up
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();

    glEnable(GL_DEPTH_TEST);
    glClearColor(0.f,0.f,0.f,0.f);
    glColor4f(1.f, 1.f, 1.f, 1.f);
}

void RenderTargetViewerPainter::paintPort(RenderPort* rp, int index, int requestedWidth, int requestedHeight) {

    if (!rp || !rp->getRenderTarget())
        return;

    RenderTarget* rt = rp->getRenderTarget();

    // set default display type if nothing has been set before
    if(renderTargetViewer_->showType_[index] == 0)
        renderTargetViewer_->showType_[index] = RenderTargetViewer::Color |
                                                RenderTargetViewer::R |
                                                RenderTargetViewer::G |
                                                RenderTargetViewer::B |
                                                RenderTargetViewer::A |
                                                RenderTargetViewer::BackgroundBlack;

    // calculate scaling
    if (renderTargetViewer_->keepAspectRatioACT_->isChecked()) {
        float scaleX = 1.0f;
        float scaleY = 1.0f;
        tgt::ivec2 size = rt->getSize();

        float aspectRatioRT = (float)size.x / size.y;
        if(requestedHeight * aspectRatioRT > requestedWidth) {
            scaleY = (requestedWidth / aspectRatioRT) / requestedHeight;
        }
        else {
            scaleX = (requestedHeight * aspectRatioRT) / requestedWidth;
        }

        // Requested size must be greater zero.
        requestedWidth = std::max<int>(requestedWidth * scaleX, 1);
        requestedHeight = std::max<int>(requestedHeight * scaleY, 1);
    }

    MatStack.pushMatrix();
    MatStack.loadIdentity();

    GLfloat currentDepthValue = 0.f;
    tgt::Vector4<GLfloat> currentTexel = tgt::Vector4<GLfloat>(1.f);
    float currentValueScale = 1.f;

    if (renderTargetViewer_->maximizeOnePort_) {

        MatStack.scale(renderTargetViewer_->zoomScale_, renderTargetViewer_->zoomScale_, 1.0);
        MatStack.translate(renderTargetViewer_->zoomTranslateX_, renderTargetViewer_->zoomTranslateY_, 0);

        // fetch depthBuffer value under cursor
        renderTargetToTexture(rt, RenderTargetViewer::Depth, colorTex_, requestedWidth, requestedHeight);
        {
            tgt::GLFramebufferObjectGuard guard(fbo_);
            fbo_->attachTexture(colorTex_);
            fbo_->isComplete();
            GLfloat tmpDepthValue;
            glReadPixels(renderTargetViewer_->mouseX_, renderTargetViewer_->mouseY_, 1, 1, GL_RED, GL_FLOAT, &tmpDepthValue);
            currentDepthValue = pow(tmpDepthValue, 50); //WTF? y?
        }

        // fetch texel value under cursor
        renderTargetToTexture(rt, RenderTargetViewer::Color | RenderTargetViewer::R | RenderTargetViewer::G | RenderTargetViewer::B |
            RenderTargetViewer::A | RenderTargetViewer::BackgroundBlack, colorTex_, requestedWidth, requestedHeight);

        if (renderTargetViewer_->mouseIsInside_)
        {
            tgt::GLFramebufferObjectGuard guard(fbo_);
            fbo_->attachTexture(rp->getColorTexture());
            fbo_->isComplete();

            tgt::vec2 p = (tgt::vec2(renderTargetViewer_->mouseX_, renderTargetViewer_->mouseY_) / tgt::vec2(requestedWidth, requestedHeight));
            p /= renderTargetViewer_->zoomScale_;
            p -= tgt::vec2(renderTargetViewer_->zoomTranslateX_, renderTargetViewer_->zoomTranslateY_) / tgt::vec2(requestedWidth, requestedHeight);

            if (p.x < 1.f && p.x > 0.f && p.y < 1.f && p.y > 0.f) {
                tgt::ivec2 origPos = p * tgt::vec2(rp->getColorTexture()->getDimensions().xy());

                // The following are the internal formats used in RenderTarget.  If other formats are added in that class, this will have to be extended.
                if (rp->getColorTexture()->getGLInternalFormat() == GL_RGB) {
                    tgt::Vector3<GLubyte> val;
                    glReadPixels(origPos.x, origPos.y, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, val.elem);
                    currentTexel = tgt::vec4(tgt::vec3(val), 255.f);
                    currentValueScale = 255.f;
                }
                else if (rp->getColorTexture()->getGLInternalFormat() == GL_RGBA) {
                    tgt::Vector4<GLubyte> val;
                    glReadPixels(origPos.x, origPos.y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, val.elem);
                    currentTexel = tgt::vec4(val);
                    currentValueScale = 255.f;
                }
                else if (rp->getColorTexture()->getGLInternalFormat() == GL_RGBA16) {
                    tgt::Vector4<GLushort> val;
                    glReadPixels(origPos.x, origPos.y, 1, 1, GL_RGBA, GL_UNSIGNED_SHORT, val.elem);
                    currentTexel = tgt::vec4(val);
                    currentValueScale = 65535.f;
                }
                else if (rp->getColorTexture()->getGLInternalFormat() == GL_R32UI) {
                    GLuint val;
                    glReadPixels(origPos.x, origPos.y, 1, 1, GL_RED_INTEGER, GL_UNSIGNED_INT, &val);
                    currentTexel = tgt::vec4(val, 0.f, 0.f, 0.f);
                    currentValueScale = 4294967295.f;
                }
                else {  // all float formats are treated equally
                    glReadPixels(origPos.x, origPos.y, 1, 1, GL_RGBA, GL_FLOAT, currentTexel.elem);
                    currentValueScale = 1.f;
                }
            }

            fbo_->detachAll();
        }
    }

    // render current type
    if(!renderTargetViewer_->maximizeOnePort_ || renderTargetViewer_->showType_[index] != (RenderTargetViewer::Color|RenderTargetViewer::R|
                                                                                           RenderTargetViewer::G|RenderTargetViewer::B|
                                                                                           RenderTargetViewer::A|RenderTargetViewer::BackgroundBlack))
        renderTargetToTexture(rt, renderTargetViewer_->showType_[index], colorTex_, requestedWidth, requestedHeight);

    if((renderTargetViewer_->showType_[index] & RenderTargetViewer::Depth) == RenderTargetViewer::Depth) {
        renderTargetToTexture(rt, RenderTargetViewer::Depth, colorTex_, requestedWidth, requestedHeight);

        //used to color scale the rendering
        GLfloat minDepth = 1;
        GLfloat maxDepth = 0;

        {
            tgt::GLFramebufferObjectGuard guard(fbo_);
            fbo_->attachTexture(colorTex_);
            fbo_->isComplete();
            size_t size = colorTex_->getWidth() * colorTex_->getHeight();
            GLfloat* texelBuffer = new GLfloat[size];
            glReadPixels(0, 0, colorTex_->getWidth(), colorTex_->getHeight(), GL_RED, GL_FLOAT, texelBuffer);

            for (size_t i = 0; i < size; i++) {
                if (texelBuffer[i] >= 0)
                    minDepth = std::min(minDepth, texelBuffer[i]);
                maxDepth = std::max(maxDepth, texelBuffer[i]);
            }
            maxDepth = std::max(maxDepth, minDepth);
            delete[] texelBuffer;
        }

        renderTargetToTexture(rt, RenderTargetViewer::Depth, colorTex_, requestedWidth, requestedHeight, minDepth, maxDepth);
    }

    MatStack.popMatrix();

    renderInfosToFontTexture(rt, requestedWidth, requestedHeight, currentDepthValue, currentTexel, currentValueScale);

    // combine color-/font-texture by shader
    paintCombinedTextures(requestedWidth, requestedHeight);

    // draw red line bounds around selected render target
    if (!renderTargetViewer_->maximizeOnePort_ && renderTargetViewer_->mouseIsInside_ && index == renderTargetViewer_->selectedRenderPortIndex_)
        paintOutline(requestedWidth, requestedHeight);
}

void RenderTargetViewerPainter::renderTargetToTexture(RenderTarget* rt, unsigned int showType, tgt::Texture* texture,
                                                      int requestedWidth, int requestedHeight,  float depthMin, float depthMax) {

    tgt::Shader* shaderProgram = colorProgram_;
    if (!shaderProgram)
        return;

    tgt::ivec3 size(requestedWidth, requestedHeight, 1);

    glActiveTexture(GL_TEXTURE0);
    // resize textures if necessary
    if (size != texture->getDimensions()) {
        texture->updateDimensions(size,true);
    }

    // activate fbo
    tgt::GLFramebufferObjectGuard guard(fbo_);
    fbo_->attachTexture(texture);
    fbo_->isComplete();

    glClearColor(0.7, 0.7, 0.7, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);
    glClearColor(0, 0, 0, 0);

    tgtAssert(rt, "No render target");

    if ((showType & RenderTargetViewer::Depth) == RenderTargetViewer::Depth) {
        rt->bindDepthTexture();
    }
    else {
        rt->bindColorTexture();
    }

    tgtAssert(shaderProgram, "No shader");
    shaderProgram->activate();
    shaderProgram->setUniform("tex_", 0);
    LGL_ERROR;
    shaderProgram->setUniform("enableColorR_", (showType & RenderTargetViewer::R) == RenderTargetViewer::R);
    shaderProgram->setUniform("enableColorG_", (showType & RenderTargetViewer::G) == RenderTargetViewer::G);
    shaderProgram->setUniform("enableColorB_", (showType & RenderTargetViewer::B) == RenderTargetViewer::B);
    shaderProgram->setUniform("enableColorA_", (showType & RenderTargetViewer::A) == RenderTargetViewer::A);
    shaderProgram->setUniform("showAlpha_", (showType & RenderTargetViewer::Alpha) == RenderTargetViewer::Alpha);
    shaderProgram->setUniform("showDepth_", (showType & RenderTargetViewer::Depth) == RenderTargetViewer::Depth);
    shaderProgram->setUniform("showHue_", (showType & RenderTargetViewer::H) == RenderTargetViewer::H);
    shaderProgram->setUniform("showSaturation_", (showType & RenderTargetViewer::S) == RenderTargetViewer::S);
    shaderProgram->setUniform("showValue_", (showType & RenderTargetViewer::V) == RenderTargetViewer::V);
    shaderProgram->setUniform("minDepth_", (float)depthMin);
    shaderProgram->setUniform("maxDepth_", (float)depthMax);
    shaderProgram->setUniform("enableBackgroundCheckerboardPattern_", (showType & RenderTargetViewer::CheckerboardPattern) == RenderTargetViewer::CheckerboardPattern);
    shaderProgram->setUniform("checkerboardPatternDimension_", 0.1f);
    float bgColor = ((showType & RenderTargetViewer::BackgroundWhite) == RenderTargetViewer::BackgroundWhite) ? 1.0f : 0.0f;
    shaderProgram->setUniform("backgroundColorR_", bgColor);
    shaderProgram->setUniform("backgroundColorB_", bgColor);
    shaderProgram->setUniform("backgroundColorG_", bgColor);

    LGL_ERROR;
    shaderProgram->setIgnoreUniformLocationError(true);
    shaderProgram->setUniform("texParameters_.dimensions_", tgt::vec2(rt->getSize()));
    shaderProgram->setUniform("texParameters_.dimensionsRCP_", 1.f / tgt::vec2(rt->getSize()));
    shaderProgram->setUniform("texParameters_.matrix_", tgt::mat4::identity);
    shaderProgram->setIgnoreUniformLocationError(false);
    LGL_ERROR;

    IMode.setMatstackUniforms(shaderProgram);
    renderQuad(requestedWidth, requestedHeight);

    shaderProgram->deactivate();
    LGL_ERROR;
}

void RenderTargetViewerPainter::renderInfosToFontTexture(RenderTarget* rt, int requestedWidth, int requestedHeight,
                                                         float currentDepthValue, tgt::vec4 currentTexel, float currentValueScale) {
    MatStack.pushMatrix();

    tgt::ivec3 size(requestedWidth, requestedHeight, 1);

    MatStack.loadIdentity();
    LGL_ERROR;

    glActiveTexture(GL_TEXTURE0);
    // resize texture if necessary
    if (size != fontTex_->getDimensions()) {
        fontTex_->updateDimensions(size,true);
    }
    LGL_ERROR;

    // activate fbo
    tgt::GLFramebufferObjectGuard guard(fbo_);
    fbo_->attachTexture(fontTex_);
    fbo_->isComplete();

    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT);
    glColor4f(1.0, 1.0, 1.0, 1.0);
    LGL_ERROR;

    if (renderTargetViewer_->showInfosACT_->isChecked()) {

        // render fonts
        QString colorStr = "";
        QString depthStr = "";
        QString sizeStr = QObject::tr("%1x%2").arg(rt->getSize().x).arg(rt->getSize().y);
        int deltaY = 15; // line height

        if (rt->getColorTexture())
            colorStr += internalTypeToString(rt->getColorTexture()->getGLInternalFormat());
        else
            colorStr += "No color texture";

        if (rt->getDepthTexture())
            depthStr += internalTypeToString(rt->getDepthTexture()->getGLInternalFormat());
        else
            depthStr += "No depth texture";

        int mem = 0;
        if (rt->getColorTexture())
            mem += rt->getColorTexture()->getSizeOnGPU();
        if (rt->getDepthTexture())
            mem += rt->getDepthTexture()->getSizeOnGPU();
        mem /= 1024;

        QString memStr = QObject::tr("%1 kb").arg(mem);
        QString numUpdateStr = QObject::tr("@%1").arg(rt->getNumUpdates());
        LGL_ERROR;

        if (renderTargetViewer_->mouseIsInside_ && renderTargetViewer_->maximizeOnePort_ && renderTargetViewer_->showInfosDetailsACT_->isChecked()){
            int textureCoordX = static_cast<int>((float)(renderTargetViewer_->zoomOffsetX_ + (renderTargetViewer_->mouseX_ - renderTargetViewer_->zoomOffsetX_) /
                renderTargetViewer_->zoomScale_) * rt->getSize().x / requestedWidth);
            int textureCoordY = static_cast<int>((float)(renderTargetViewer_->zoomOffsetY_ + (renderTargetViewer_->mouseY_ - renderTargetViewer_->zoomOffsetY_) /
                renderTargetViewer_->zoomScale_) * rt->getSize().y / requestedHeight);

            if (textureCoordX < rt->getSize().x && textureCoordY < rt->getSize().y) {
                int numLines = 7;
                int offsetX = 3;
                int offsetY = requestedHeight - numLines * deltaY - 3;

                renderFont(offsetX, offsetY, QObject::tr("Y: %1").arg(textureCoordY));

                offsetY += deltaY;
                renderFont(offsetX, offsetY, QObject::tr("X: %1").arg(textureCoordX));

                offsetY += deltaY;
                renderFont(offsetX, offsetY, QObject::tr("D: %1").arg(currentDepthValue));

                offsetY += deltaY;
                std::stringstream ssA;
                ssA << "A: " << currentTexel[3] << " " << std::setprecision(2) << std::fixed << currentTexel[3] / currentValueScale;
                renderFont(offsetX, offsetY, QObject::tr(ssA.str().c_str()));

                tgt::vec3 hsv = rgbToHsv(currentTexel.xyz() / currentValueScale);

                offsetY += deltaY;
                std::stringstream ssB;
                ssB << "B: " << currentTexel[2] << " " << std::setprecision(2) << std::fixed << currentTexel[2] / currentValueScale;
                ssB << " " << "V: " << std::fixed << hsv[2];
                renderFont(offsetX, offsetY, QObject::tr(ssB.str().c_str()));

                offsetY += deltaY;
                std::stringstream ssG;
                ssG << "G: " << currentTexel[1] << " " << std::setprecision(2) << std::fixed << currentTexel[1] / currentValueScale;
                ssG << " " << "S: " << std::fixed << hsv[1];
                renderFont(offsetX, offsetY, QObject::tr(ssG.str().c_str()));

                offsetY += deltaY;
                std::stringstream ssR;
                ssR << "R: " << currentTexel[0] << " " << std::setprecision(2) << std::fixed << currentTexel[0] / currentValueScale;
                ssR << " " << "H: " << std::setfill('0') << std::setw(3) << (int)hsv[0];
                renderFont(offsetX, offsetY, QObject::tr(ssR.str().c_str()));
            }
        }

        LGL_ERROR;

        int offsetX = 3;
        int offsetY = 5;

        renderFont(offsetX, offsetY, QObject::tr(rt->getDebugLabel().c_str()));
        offsetY += deltaY;
        renderFont(offsetX, offsetY, depthStr);
        offsetY += deltaY;
        renderFont(offsetX, offsetY, colorStr);
        offsetY += deltaY;
        renderFont(offsetX, offsetY, sizeStr);
        offsetY += deltaY;
        renderFont(offsetX, offsetY, memStr);
        offsetY += deltaY;
        renderFont(offsetX, offsetY, numUpdateStr);
        LGL_ERROR;
    }

    MatStack.popMatrix();
    glColor4f(1, 1, 1, 1);

    LGL_ERROR;
}

void RenderTargetViewerPainter::paintCombinedTextures(int requestedWidth, int requestedHeight) {
    tgt::Shader* shaderProgram = inversecolorProgram_;
    if (!shaderProgram)
        return;

    glActiveTexture(GL_TEXTURE0);
    LGL_ERROR;
    colorTex_->bind();
    LGL_ERROR;
    glActiveTexture(GL_TEXTURE1);
    LGL_ERROR;
    fontTex_->bind();
    LGL_ERROR;
    glActiveTexture(GL_TEXTURE0);
    LGL_ERROR;

    tgtAssert(shaderProgram, "No shader");
    shaderProgram->activate();
    LGL_ERROR;

    shaderProgram->setUniform("backTex_", 0);
    shaderProgram->setUniform("frontTex_", 1);
    shaderProgram->setUniform("threshold_", 0.15f);
    LGL_ERROR;
    shaderProgram->setIgnoreUniformLocationError(true);
    shaderProgram->setUniform("texParameters_.dimensions_", tgt::vec2(colorTex_->getDimensions().xy()));
    shaderProgram->setUniform("texParameters_.dimensionsRCP_", 1.f / tgt::vec2(colorTex_->getDimensions().xy()));
    shaderProgram->setUniform("texParameters_.matrix_", tgt::mat4::identity);
    shaderProgram->setIgnoreUniformLocationError(false);
    LGL_ERROR;
    IMode.setMatstackUniforms(shaderProgram);
    renderQuad(requestedWidth, requestedHeight);
    LGL_ERROR;
    shaderProgram->deactivate();
    LGL_ERROR;
}

void RenderTargetViewerPainter::paintOutline(int requestedWidth, int requestedHeight) {
    IMode.begin(tgt::ImmediateMode::LINE_LOOP);
    IMode.color(tgt::vec4(1.0f, 0.0f, 0.0f, 1.0f));
    IMode.vertex(tgt::vec2(0, 0));
    IMode.vertex(tgt::vec2(requestedWidth-1, 0));
    IMode.vertex(tgt::vec2(requestedWidth-1, requestedHeight-1));
    IMode.vertex(tgt::vec2(0, requestedHeight-1));
    IMode.end();
    IMode.color(tgt::vec4::one);
    LGL_ERROR;
}

void RenderTargetViewerPainter::renderQuad(int requestedWidth, int requestedHeight) {
    glBindVertexArray(quadVao_);
    glBindBuffer(GL_ARRAY_BUFFER, quadVbo_);

    // position
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, false, sizeof(tgt::vec2)*2, 0);
    // texcoord
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, false, sizeof(tgt::vec2)*2, (GLvoid*)sizeof(tgt::vec2));

    tgt::vec2 buffer[] = {
        tgt::vec2(0,0), tgt::vec2(0,0),
        tgt::vec2(requestedWidth,0), tgt::vec2(1,0),
        tgt::vec2(requestedWidth,requestedHeight), tgt::vec2(1,1),
        tgt::vec2(0,requestedHeight), tgt::vec2(0,1),
        tgt::vec2(0,0), tgt::vec2(0,0)
    };
    glBufferData(GL_ARRAY_BUFFER, sizeof(buffer), buffer, GL_DYNAMIC_DRAW);

    glDrawArrays(GL_TRIANGLE_STRIP , 0, 5);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void RenderTargetViewerPainter::renderFont(float x, float y, QString text) {
    font_->render(tgt::vec3(x, y, 0.0), text.toStdString(), canvas_->getPhysicalSize());
    LGL_ERROR;
}

} //namespace
