/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2024 University of Muenster, Germany,           *
 * Department of Computer Science.                                    *
 *                                                                    *
 * This file is part of the tgt library. This library is free         *
 * software; you can redistribute it and/or modify it under the terms *
 * of the GNU Lesser General Public License version 2.1 as published  *
 * by the Free Software Foundation.                                   *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       *
 * GNU Lesser General Public License for more details.                *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License in the file "LICENSE.txt" along with this library.         *
 * If not, see <http://www.gnu.org/licenses/>.                        *
 *                                                                    *
 **********************************************************************/

#include "tgt/framebufferobject.h"
#include "tgt/glcontextmanager.h"
#include "tgt/logmanager.h"

namespace tgt {

const std::string FramebufferObject::loggerCat_("tgt.FramebufferObject");

FramebufferObject::FramebufferObject()
  : id_(0)
{
    generateId();
}

FramebufferObject::~FramebufferObject()
{
    glDeleteFramebuffers(1, &id_);
}

void FramebufferObject::activate()
{
    glBindFramebuffer(GL_FRAMEBUFFER, id_);
}

void FramebufferObject::deactivate()
{
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void FramebufferObject::attachTexture(Texture* texture, GLenum attachment, int mipLevel, int zSlice)
{
    switch(texture->getType()) {
        case GL_TEXTURE_1D:
            glFramebufferTexture1D(GL_FRAMEBUFFER, attachment, GL_TEXTURE_1D, texture->getId(), mipLevel);
            break;
        case GL_TEXTURE_3D:
            glFramebufferTexture3D(GL_FRAMEBUFFER, attachment, GL_TEXTURE_3D, texture->getId(), mipLevel, zSlice);
            break;
        case GL_TEXTURE_2D_ARRAY:
            glFramebufferTextureLayer(GL_FRAMEBUFFER, attachment, texture->getId(), mipLevel, zSlice);
            break;
        default: //GL_TEXTURE_2D, GL_TEXTURE_RECTANGLE
            glFramebufferTexture2D(GL_FRAMEBUFFER, attachment, texture->getType(), texture->getId(), mipLevel);
            break;
    }
    attachedTextures_[attachment] = texture;
}

Texture* FramebufferObject::getTextureAtAttachment(GLenum attachment) {
    std::map<GLenum, Texture*>::iterator iter = attachedTextures_.find(attachment);
    if( iter != attachedTextures_.end() ) {
        return attachedTextures_[attachment];
    }
    else
        return 0;
}

void FramebufferObject::detachTexture(GLenum attachment) {
    std::map<GLenum, Texture*>::iterator iter = attachedTextures_.find(attachment);
    if( iter != attachedTextures_.end() ) {
        attachedTextures_.erase(iter);
    }
    else {
        LWARNING("Trying to detach unknown texture!");
    }

    glFramebufferTexture2D(GL_FRAMEBUFFER, attachment, GL_TEXTURE_2D, 0, 0);
}

void FramebufferObject::detachAll() {
    while(!attachedTextures_.empty()) {
        detachTexture(attachedTextures_.begin()->first);
    }
}

bool FramebufferObject::isComplete()
{
    bool complete = false;

    GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
    switch(status) {
        case GL_FRAMEBUFFER_COMPLETE:
            complete = true;
            break;
        case GL_FRAMEBUFFER_UNDEFINED:
            LERROR("GL_FRAMEBUFFER_UNDEFINED");
            break;
        case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT:
            LERROR("GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT");
            break;
        case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT:
            LERROR("GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT");
            break;
        case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER:
            LERROR("GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER");
            break;
        case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER:
            LERROR("GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER");
            break;
        case GL_FRAMEBUFFER_UNSUPPORTED:
            LERROR("GL_FRAMEBUFFER_UNSUPPORTED");
            break;
        case GL_FRAMEBUFFER_INCOMPLETE_MULTISAMPLE:
            LERROR("GL_FRAMEBUFFER_INCOMPLETE_MULTISAMPLE");
            break;
        case GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS:
            LERROR("GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS");
            break;
        default:
            LERROR("Unknown error!");
    }

    return complete;
}

bool FramebufferObject::isActive() const {
    return ((getActiveObject() == id_) && (id_ != 0));
}

GLuint FramebufferObject::getActiveObject() {
    GLint fbo;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &fbo);
    return static_cast<GLuint>(fbo);
}

GLuint FramebufferObject::generateId() {
    id_ = 0;
    glGenFramebuffers(1, &id_);
    return id_;
}

GLFramebufferObjectGuard::GLFramebufferObjectGuard(FramebufferObject* fbo)
{
    context_ = GLContextMgr.getActiveContext();
    id_ = FramebufferObject::getActiveObject();
    if (fbo)
        fbo->activate();
}

GLFramebufferObjectGuard::~GLFramebufferObjectGuard() {
    tgtAssert(GLContextMgr.getActiveContext() == context_, "FBO was reactivated in wrong context");
    glBindFramebuffer(GL_FRAMEBUFFER, id_);
}

} // namespace
