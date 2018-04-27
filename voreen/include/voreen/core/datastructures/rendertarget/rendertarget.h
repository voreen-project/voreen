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

#ifndef VRN_RENDERTARGET_H
#define VRN_RENDERTARGET_H

#include "voreen/core/voreencoreapi.h"
#include "tgt/framebufferobject.h"

namespace voreen {

class VRN_CORE_API RenderTarget {
public:
    RenderTarget();
    virtual ~RenderTarget();

    virtual void initialize(GLint internalColorFormat = GL_RGBA16, GLint internalDepthFormat = GL_DEPTH_COMPONENT24);
    void deinitialize();

    void activateTarget(const std::string& debugLabel = "");
    void deactivateTarget();
    bool isActive() const;

    void bindColorTexture(GLint texUnit, GLint filterMode = GL_LINEAR, GLint wrapMode = GL_CLAMP_TO_EDGE, tgt::vec4 borderColor = tgt::vec4(0.f));
    void bindDepthTexture(GLint texUnit, GLint filterMode = GL_LINEAR, GLint wrapMode = GL_CLAMP_TO_EDGE, tgt::vec4 borderColor = tgt::vec4(0.f));

    void bindColorTexture();
    void bindDepthTexture();

    const tgt::Texture* getColorTexture() const { return colorTex_; }
    virtual const tgt::Texture* getDepthTexture() const { return depthTex_; }

    tgt::Texture* getColorTexture() { return colorTex_; }
    tgt::Texture* getDepthTexture() { return depthTex_; }

    /** 
     * Get the color at the specified position.
     * This method activates the target, reads the color and will deactivate the target.
     */
    tgt::vec4 getColorAtPos(tgt::ivec2 pos);

    tgt::ivec2 getSize() const;
    void resize(tgt::ivec2 newsize);

    //these are just for debugging purposes:
    std::string getDebugLabel() { return debugLabel_; }
    void setDebugLabel(const std::string& debugLabel) { debugLabel_ = debugLabel; }

    void increaseNumUpdates() { numUpdates_++; }
    int getNumUpdates() { return numUpdates_; }

    virtual bool hasStencilBuffer() const { return false; }

protected:

    tgt::FramebufferObject* fbo_;

    tgt::Texture* colorTex_;
    tgt::Texture* depthTex_;

    std::string debugLabel_;
    int numUpdates_;

    static const std::string loggerCat_;
};

}   // namespace

#endif  // VRN_RENDERTARGET_H
