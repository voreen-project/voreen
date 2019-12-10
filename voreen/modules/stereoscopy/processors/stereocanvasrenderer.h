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

#ifndef VRN_STEREOCANVASRENDERER_H
#define VRN_STEREOCANVASRENDERER_H

//super class
#include "modules/core/processors/output/canvasrenderer.h"

// properties
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/interaction/camerainteractionhandler.h"

namespace tgt {
class GLCanvas;
class Texture;
}

namespace voreen {

/**
 * Supports several stereo modes: side by side, autostereoscopic, anaglyph, quad-buffered.
 * To be used as drop-in replacement for the standard canvas.
 */
class VRN_CORE_API StereoCanvasRenderer : public CanvasRenderer {
    friend class StereoCanvasRendererWidget;
public:
    /** Constructor */
    StereoCanvasRenderer();
    /** Destructor */
    ~StereoCanvasRenderer();
    //--------------------------------------
    //  basic processor functions
    //--------------------------------------
    virtual Processor* create() const        { return new StereoCanvasRenderer();}
    virtual std::string getClassName() const { return "StereoCanvas";            }
    virtual std::string getCategory() const  { return "Output";                  }
    virtual CodeState getCodeState() const   { return CODE_STATE_TESTING;        }
protected:
    virtual void setDescriptions() {
        setDescription("Supports several stereo modes: side-by-side, auto-stereoscopic, anaglyph and quad-buffered. "
                       "To be used as drop-in replacement for the standard canvas. "
                       "<p>See: <a href=\"http://voreen.uni-muenster.de/?q=stereoscopy\" >voreen.uni-muenster.de/?q=stereoscopy</a></p>");
    }
    virtual void initialize();
    virtual void deinitialize();
    virtual void onEvent(tgt::Event* e);
    virtual void process();

    //--------------------------------------
    //  resize functions
    //--------------------------------------
public:
    virtual void canvasResized(tgt::ivec2 newsize);
    virtual void resizeCanvas(tgt::ivec2 newsize);
    virtual void setCanvas(tgt::GLCanvas* canvas);
    void resizeAllPorts(tgt::ivec2 newsize);

    //--------------------------------------
    //  get/renderTo texture functions
    //--------------------------------------
    virtual const tgt::Texture* getImageColorTexture() const;
    virtual tgt::Texture* getImageColorTexture();
    virtual const tgt::Texture* getImageDepthTexture() const;
    virtual tgt::Texture* getImageDepthTexture();
    virtual bool renderToImage(const std::string &filename);
    virtual bool renderToImage(const std::string &filename, tgt::ivec2 dimensions);

protected:
    /** different types of supported stereo modes */
    enum StereoMode {
        NO_STEREO_MODE = 0,                 ///< no stereo, 'this' can be used as a normal canvas
        SPLITSCREEN_STEREO_MODE = 1,        ///< left eye and right eye are rendered side by side
        QUADBUFFER_STEREO_MODE = 2,         ///< this function is not supported yet
        ANAGLYPH_STEREO_MODE = 3,           ///< left and right eye are color coded in the same image
        AUTOSTEREOSCOPIC_STEREO_MODE = 4    ///< used for autostereoscopic displays. The frustum will be modified
    };
    friend class OptionProperty<StereoMode>;

    /** different types of handling the incomming image */
    enum NextExpectedImageMode {
        NORMAL = 0, ///< camera is in middle eye. Camera will be set to left eye and mode is changed to LEFT_1
        LEFT_1 = 1, ///< camera is in left eye. Camera will be set to right eye and mode is changed to RIGHT_1
        RIGHT_1 = 2,///< camera is in right eye. If lastRunWasInInteractionMode_ equals interactionMode() the camera will stay right and mode is changed to RIGHT_2.
                        // Otherwise camera is set left and mode is changed to LEFT_2.
        LEFT_2 = 3, ///< camera is in left eye. If lastRunWasInInteractionMode_ equals interactionMode() the camera will stay right and mode is changed to LEFT_1.
                        // Otherwise camera is set right and mode is changed to RIGHT_1.
        RIGHT_2 = 4 ///< camera is in right eye. Camera will be set to left eye and mode is changed to LEFT_2
    };
    friend class OptionProperty<NextExpectedImageMode>;

    /** different anaglyph color settings */
    enum AnaglyphMode {
        RED_GREEN = 0,      ///< red/green eye encoding
        RED_BLUE = 1,       ///< red/blue eye encoding
        RED_CYAN = 2        ///< red/cyan eye encoding. This mode retains gray-scale values
    };
    friend class OptionProperty<AnaglyphMode>;

    /** different autostereoscopic modes */
    enum AutostereoscopicMode {
        VERTICAL_INTERLEAVED = 0,      ///< vertical lines
        HORIZONTAL_INTERLEAVED = 1,    ///< horizontal lines
        CHECKER_INTERLEAVED = 2        ///< checkerboard lines
    };
    friend class OptionProperty<AutostereoscopicMode>;


    //--------------------------------------
    //  ports
    //--------------------------------------
    RenderPort tmpPort_;        ///< (private) port to stash the first (eye) image while waiting for the second
    RenderPort sideBySidePort_; ///< (private) port to store the last rendered two (eye) images side by side
    RenderPort finalPort_;      ///< (private) port to store the last completely rendered image

    //--------------------------------------
    //  properties
    //--------------------------------------
        //camera settings
    CameraProperty cameraProp_;                     ///< the camera
    void invalidateSideBySidePort();
    void invalidateFinalPort();
        //stereoscopic method
    OptionProperty<StereoMode> stereoModeProp_;     ///< the actual stereo mode
        void stereoModeOnChange();
    OptionProperty<AnaglyphMode> anaglyphModeProp_; ///< the actual anaglyph mode
    OptionProperty<AutostereoscopicMode> autostereoscopicModeProp_; ///< the current autostereoscopic mode
    BoolProperty eyeInvertProp_;                    ///< should eyes been inverted?
    BoolProperty calibrateDisplayProp_;             ///< show calibration texture?
        //stereo settings
    FloatProperty eyeSeparationProp_;               ///< eye separation in cm
        void eyeSeparationOnChange();
    FloatProperty focalLengthProp_;                 ///< focal length == real-world distance to monitor in cm
        void focalLengthOnChange();
    FloatProperty focalWidthProp_;                  ///< focal width == real-world display width in cm
        void focalWidthOnChange();
    FloatProperty relativeFocalLengthProp_;         ///< stereo focal length relative to distance between near and far plane (for non-real world use)
        void relativeFocalLengthPropOnChange();
    BoolProperty useRealWorldFrustumProp_;          ///< use real world frustum for stereo rendering or scale it to match internal frustum
        void useRealWorldFrustumPropOnChange();
    OptionProperty<tgt::Camera::StereoAxisMode> stereoAxisModeProp_;    ///< the stereo camera axis mode
        void stereoAxisModeOnChange();
        //events
    EventProperty<StereoCanvasRenderer> mouseMoveEventProp_;    ///< handle viewportin split screen mode
        void mouseMove(tgt::MouseEvent* e);

private:
    //--------------------------------------
    //  copy functions
    //--------------------------------------
        /**
         * Main function controlling the stereo rendering.
         * It is called in 'process' if the StereoMode is not NO_STEREO_MODE.
         * @see NextExpectedImageMode
         */
    void processStereo();
        /**
         * Copies the texture of the finalPort into the canvas buffer.
         */
    void copyIntoCanvas();
        /**
         * Copies the texture of 'input' into the 'output'.
         * @param input renderport which texture is been copied. Normally the Inport.
         * @param output renderport which the texture is copied to. Normally the tmpPort.
         */
    void copyIntoPort(RenderPort* input, RenderPort* output);
        /**
         * Copies the two textures into the sideBySidePort. If eyes are inverted, left and right will be changed automatically.
         * Will call 'copyIntoFinal' afterwards'
         * @param colorLeft color texture of the left eye.
         * @param colorRight color tecture of right eye.
         * @param depthLeft depth texture of the left eye.
         * @param depthRight depth texture of the right eye.
         */
    void copyIntoSideBySide(tgt::Texture* colorLeft, tgt::Texture* colorRight, tgt::Texture* depthLeft = 0, tgt::Texture* depthRight = 0);
        /**
         * Copies the two textures of the renderports into the sideBySidePort.
         * @overload
         * @param left renderport containing the texture of the left eye.
         * @param right renderport containing the texture of the right eye.
         */
    void copyIntoSideBySide(RenderPort* left, RenderPort* right);

        /**
         * Copies the sideBySidePort into the finalPort according to the current stereo mode.
         */
    void copyIntoFinal();

    //--------------------------------------
    //  render functions
    //--------------------------------------
    void renderAnaglyph();
    void renderAutostereoscopic();
    void renderSplitScreen();

    void renderFullScreen();
    //--------------------------------------
    //  members
    //--------------------------------------
    //textures
    tgt::Texture* calibrationRightTexture_;             ///< calibration texture for the right eye
    tgt::Texture* calibrationLeftTexture_;              ///< calibration texture for the left eye
    //shaders
    tgt::Shader* sideBySideShader_;         ///< shader for copying two images into side by side (copyimagesplitscreen.frag)
    tgt::Shader* anaglyphShader_;           ///< shader for anaglyph rendering (anaglyph.frag)
    tgt::Shader* autostereoscopicShader_;   ///< shader for autostereoscopic rendering (autostereoscopic.frag)
    tgt::Shader* copyTextureShader_;        ///< simple shader for copying two textures
    //last camera setting
    tgt::Camera::ProjectionMode previousCameraProjectionMode_;  ///< projection mode of the camera before changed to FRUSTUM
    //process stereo member
    NextExpectedImageMode nextExpectedImage_;   ///< next expected image mode
    bool lastRunWasInInteractionMode_;          ///< last image was rendered in interaction mode?
    bool w8ingOn2Eye_;                          ///< are we waiting on second image? (used for head tracking)

    //buffers used for rendering fullscreen
    GLuint fullScreenVertexArrayObjectID_;
    GLuint fullScreenVertexBufferObjectID_;

    static const std::string loggerCat_;    ///< category used in logging
};

} // namespace voreen

#endif // VRN_STEREOCANVASRENDERER_H
