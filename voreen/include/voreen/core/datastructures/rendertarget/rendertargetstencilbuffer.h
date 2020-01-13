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

#ifndef VRN_RENDERTARGETSTENCILBUFFER_H
#define VRN_RENDERTARGETSTENCILBUFFER_H

#include "rendertarget.h"

namespace voreen {

class VRN_CORE_API RenderTargetStencilBuffer : public RenderTarget {
public:
    RenderTargetStencilBuffer();

    /**
     * Initializes the RenserTargetStencilBuffer with the specified format.
     *
     * @param internalColorFormat the format of the color texture
     * @param internalDepthStencilFormat must be GL_DEPTH24_STENCIL8 or GL_DEPTH32F_STENCIL8
     */
    virtual void initialize(GLint internalColorFormat = GL_RGBA16, GLint internalDepthStencilFormat = GL_DEPTH24_STENCIL8);

    //void activateTarget(const std::string& debugLabel = "");
    //void deactivateTarget();
    //bool isActive() const;

    //void bindColorTexture(GLint texUnit, GLint filterMode = GL_LINEAR, GLint wrapMode = GL_CLAMP_TO_EDGE, tgt::vec4 borderColor = tgt::vec4(0.f));

    //void bindColorTexture();

    //const tgt::Texture* getColorTexture() const { return colorTex_; }
    //tgt::Texture* getColorTexture() { return colorTex_; }

    //Get the color at position pos. This method activates the RenderTarget!
    //tgt::vec4 getColorAtPos(tgt::ivec2 pos);

    //tgt::ivec2 getSize() const;
    //void resize(tgt::ivec2 newsize);

    //these are just for debugging purposes:
    //std::string getDebugLabel() { return debugLabel_; }
    //void setDebugLabel(const std::string& debugLabel) { debugLabel_ = debugLabel; }

    //void increaseNumUpdates() { numUpdates_++; }
    //int getNumUpdates() { return numUpdates_; }

    virtual bool hasStencilBuffer() const { return true; }

protected:

    //void bindDepthTexture();
    //void bindDepthTexture(GLint texUnit, GLint filterMode = GL_LINEAR, GLint wrapMode = GL_CLAMP_TO_EDGE, tgt::vec4 borderColor = tgt::vec4(0.f));

    //tgt::Texture* getDepthTexture() { return depthTex_; }
    //const tgt::Texture* getDepthTexture() const { return depthTex_; }

    //tgt::FramebufferObject* fbo_;

    //tgt::Texture* colorTex_;
    //tgt::Texture* depthTex_;

    //std::string debugLabel_;
    //int numUpdates_;

    static const std::string loggerCat_;
};

}   // namespace

#endif  // VRN_RENDERTARGET_H
