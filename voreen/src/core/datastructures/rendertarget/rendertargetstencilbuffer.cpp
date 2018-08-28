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

#include "voreen/core/datastructures/rendertarget/rendertargetstencilbuffer.h"

#include "tgt/logmanager.h"
#include "tgt/gpucapabilities.h"

using tgt::Texture;

namespace voreen {

const std::string RenderTargetStencilBuffer::loggerCat_ = "voreen.RenderTargetStencilBuffer";

RenderTargetStencilBuffer::RenderTargetStencilBuffer()
    : RenderTarget()
{ }

void RenderTargetStencilBuffer::initialize(GLint internalColorFormat, GLint internalDepthStencilFormat) {

    if (fbo_)
        return;

    if (!GpuCaps.isNpotSupported() && !GpuCaps.areTextureRectanglesSupported()) {
        LWARNING("Neither non-power-of-two textures nor texture rectangles seem to be supported!");
    }

    tgt::ivec3 size(128, 128, 1);

    switch(internalColorFormat) {
        case GL_RGB:
            colorTex_ = new Texture(size, GL_RGB, GL_RGB, GL_UNSIGNED_BYTE, Texture::LINEAR, Texture::CLAMP_TO_EDGE);
            break;
        case GL_RGBA:
            colorTex_ = new Texture(size, GL_RGBA, GL_RGBA, GL_UNSIGNED_BYTE, Texture::LINEAR, Texture::CLAMP_TO_EDGE);
            break;
        case GL_RGBA16:
            colorTex_ = new Texture(size, GL_RGBA, GL_RGBA16, GL_UNSIGNED_SHORT, Texture::LINEAR, Texture::CLAMP_TO_EDGE);
            break;
        case GL_RGB16F:
            colorTex_ = new Texture(size, GL_RGB, GL_RGB16F, GL_FLOAT, Texture::LINEAR, Texture::CLAMP_TO_EDGE);
            break;
        case GL_RGBA16F:
            colorTex_ = new Texture(size, GL_RGBA, GL_RGBA16F, GL_FLOAT, Texture::LINEAR, Texture::CLAMP_TO_EDGE);
            break;
        case GL_RGBA32F:
            colorTex_ = new Texture(size, GL_RGBA, GL_RGBA32F, GL_FLOAT, Texture::LINEAR, Texture::CLAMP_TO_EDGE);
            break;
        case GL_RED:
            colorTex_ = new Texture(size, GL_RED, GL_RED, GL_UNSIGNED_BYTE, Texture::LINEAR, Texture::CLAMP_TO_EDGE);
            break;
        case GL_R32F:
            colorTex_ = new Texture(size, GL_RED, GL_R32F, GL_FLOAT, Texture::LINEAR, Texture::CLAMP_TO_EDGE);
            break;
        case GL_R32UI:
            colorTex_ = new Texture(size, GL_RED_INTEGER, GL_R32UI, GL_UNSIGNED_INT , Texture::LINEAR, Texture::CLAMP_TO_EDGE);
            break;

        default:
            LERROR("Unknown internal format!");
    }
    colorTex_->uploadTexture();

    switch(internalDepthStencilFormat) {
        case GL_DEPTH24_STENCIL8:
            depthTex_ = new Texture(size, GL_DEPTH_STENCIL, GL_DEPTH24_STENCIL8, GL_UNSIGNED_INT_24_8, Texture::LINEAR, Texture::CLAMP_TO_EDGE);
            break;
#ifdef GL_DEPTH32F_STENCIL8
        case GL_DEPTH32F_STENCIL8:
            depthTex_ = new Texture(size, GL_DEPTH_STENCIL, GL_DEPTH32F_STENCIL8, GL_FLOAT_32_UNSIGNED_INT_24_8_REV, Texture::LINEAR, Texture::CLAMP_TO_EDGE);
            break;
#endif
        default:
            depthTex_ = new Texture(size, GL_DEPTH_STENCIL, GL_DEPTH24_STENCIL8, GL_UNSIGNED_INT_24_8, Texture::LINEAR, Texture::CLAMP_TO_EDGE);
            LERROR("Unknown internal depth format!");
    }
    depthTex_->uploadTexture();

    fbo_ = new tgt::FramebufferObject();
    if (!fbo_) {
        LERROR("Failed to initialize framebuffer object!");
        return;
    }
    fbo_->activate();

    fbo_->attachTexture(colorTex_);
    fbo_->isComplete();

    fbo_->attachTexture(depthTex_, GL_DEPTH_STENCIL_ATTACHMENT);
    fbo_->isComplete();

    fbo_->deactivate();
}

//void RenderTarget::bindDepthTexture() {
//    tgtAssert(depthTex_, "No depth texture available!");
//    if(depthTex_)
//        depthTex_->bind();
//}

//void RenderTarget::bindDepthTexture(GLint texUnit, GLint filterMode/* = GL_LINEAR*/, GLint wrapMode /*= GL_CLAMP_TO_EDGE*/, tgt::vec4 borderColor /*= tgt::vec4(0.f)*/) {
//    glActiveTexture(texUnit);
//    bindDepthTexture();

    // texture filtering
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, filterMode);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, filterMode);
//    LGL_ERROR;

    // texture wrapping
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrapMode);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrapMode);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, wrapMode);
//    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, static_cast<tgt::Vector4<GLfloat> >(borderColor).elem);
//    LGL_ERROR;
//}

}   // namespace
