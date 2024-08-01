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

#ifndef VRN_RENDERTARGETVIEWERPAINTER_H
#define VRN_RENDERTARGETVIEWERPAINTER_H

#include "tgt/painter.h"
#include "tgt/texture.h"

#include "voreen/qt/voreenqtapi.h"

#include <QFont>
#include <QString>


namespace tgt {
    class FramebufferObject;
    class Shader;
    class Font;
}

namespace voreen {

    class RenderTargetViewer;
    class RenderPort;
    class RenderTarget;

/**
 * Painter of the render target viewer
 */
class VRN_QT_API RenderTargetViewerPainter : public tgt::Painter {

        //to access paint
        friend class QtCanvas;

public:
    /// Constructor
    RenderTargetViewerPainter(tgt::GLCanvas* canvas, RenderTargetViewer* viewer);
    virtual ~RenderTargetViewerPainter();

    //----------------------------------------
    //  functions to override
    //----------------------------------------
protected:
    /** @override */
    virtual void paint();
    /** @override */
    virtual void sizeChanged(const tgt::ivec2&);
    /** @override */
    virtual void initialize();

private:
    //----------------------------------------
    //  paint helper
    //----------------------------------------
    /**
     * Renders on render port
     */
    void paintPort(RenderPort* rp, int index, int requestedWidth, int requestedHeight);
    /**
     * Renders the informations (showType) of the render target (rt) into the texture (texture).
     */
    void renderTargetToTexture(RenderTarget* rt, unsigned int showType, tgt::Texture* texture, int requestedWidth, int requestedHeight,
                               float depthMin = 0.f, float depthMax = 1.f);
    /**
     * Render informations into font texture.
     */
    void renderInfosToFontTexture(RenderTarget* rt, int requestedWidth, int requestedHeight, float currentDepthValue, tgt::vec4 currentTexel, float currentValueScale);
    /**
     * Renders the font. Uses tgt or Qt as defined.
     */
    void renderFont(float x, float y, QString text);
    /**
     * Combines color and font texture.
     */
    void paintCombinedTextures(int requestedWidth, int requestedHeight);
    /**
     * Draws an outline to highlight the selected port.
     */
    void paintOutline(int requestedWidth, int requestedHeight);
    /**
     * Renders a full screen quad.
     */
    void renderQuad(int requestedWidth, int requestedHeight);

    //----------------------------------------
    //  member
    //----------------------------------------
    RenderTargetViewer* renderTargetViewer_; ///< associated viewer
    //tgt::QtCanvas* canvas_;                ///< canvas defined in super class

    tgt::Shader* colorProgram_;
    tgt::Shader* inversecolorProgram_;

    tgt::FramebufferObject* fbo_;
    tgt::Texture* colorTex_;
    tgt::Texture* fontTex_;

    GLuint quadVbo_;
    GLuint quadVao_;

    tgt::Font* font_;

    static const std::string loggerCat_;
};

} // namespace

#endif // VRN_RENDERTARGETVIEWERPAINTER_H
