/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2018 University of Muenster, Germany,           *
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

#ifndef TGT_FRAMEBUFFEROBJECT_H
#define TGT_FRAMEBUFFEROBJECT_H

#include "tgt/texture.h"
#include "tgt/types.h"

#include <map>

namespace tgt {

class TGT_API FramebufferObject {
public:
    FramebufferObject();
    virtual ~FramebufferObject();

    /**
     * Activates the FBO.
     *
     * Caution (1): It is still necessary to call glDrawBuffers if the FBO contains several attachments.
     * Caution (2): FBOs are not shared among contexts!!!
     */
    void activate();

    /**
     * Deactivates the FBO by binding the default FBO (0).
     */
    static void deactivate();

    /**
     * Determines, whether this FBO is currently bound.
     */
    bool isActive() const;

    /**
     * Determines whether the FBO is complete.
     */
    static bool isComplete();

    /**
     * Returns the currently active FBO name.
     */
    static GLuint getActiveObject();

    /// Bind a texture to the "attachment" point of this FBO
    void attachTexture(Texture* texture,
                       GLenum attachment = GL_COLOR_ATTACHMENT0,
                       int mipLevel      = 0,
                       int zSlice        = 0);

    void detachTexture(GLenum attachment);

    void detachAll();

    Texture* getTextureAtAttachment(GLenum attachment);

protected:
    GLuint generateId();

    GLuint id_;
    std::map<GLenum, Texture*> attachedTextures_;

    static const std::string loggerCat_; ///< category used in logging
};

class GLContextBase;

/**
 * This guard ensures that a given FBO is active at the time of instantiation
 * and the previously active FBO will be made active again at the time of destruction.
 */
class TGT_API GLFramebufferObjectGuard {
public:
    GLFramebufferObjectGuard(FramebufferObject* fbo = nullptr);
    ~GLFramebufferObjectGuard();
private:
    GLContextBase* context_; // context is used for debugging purposes
    GLint id_;
};

} // namespace tgt

#endif // TGT_FRAMEBUFFEROBJECT_H
