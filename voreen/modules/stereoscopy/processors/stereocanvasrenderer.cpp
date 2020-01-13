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

#include "stereocanvasrenderer.h"

#include "voreen/core/utils/voreenpainter.h"
#include "voreen/core/utils/stringutils.h"
#include "voreen/core/processors/processorwidget.h"
#include "voreen/core/processors/processorwidgetfactory.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/voreenqualitymode.h"

#include "tgt/glcanvas.h"
#include "tgt/texturemanager.h"
#include "tgt/gpucapabilities.h"
#include "tgt/filesystem.h"
#include "tgt/camera.h"
#include "tgt/glcontextmanager.h"

namespace voreen {

const std::string StereoCanvasRenderer::loggerCat_("voreen.stereoscopy.StereoCanvasRenderer");

//------------------------------------------------------------------------------------------------------------------
//      basic processor functions
//------------------------------------------------------------------------------------------------------------------
StereoCanvasRenderer::StereoCanvasRenderer()
    : CanvasRenderer()
    //ports
    , tmpPort_(Port::OUTPORT, "private.tmp", "TmpPort")
    , sideBySidePort_(Port::OUTPORT, "private.sidebyside", "StoragePort")
    , finalPort_(Port::OUTPORT, "private.final", "FinalPort")
    //properties
        //camera settings
    , cameraProp_("cameraProp","Camera", tgt::Camera(tgt::vec3(0.0f, 0.0f, 50.f), tgt::vec3(0.0f, 0.0f, 0.0f),
                                                 tgt::vec3(0.0f, 1.0f, 0.0f), 45.f,1.f,0.1f,500.f))
    //stereoscopic method
    , stereoModeProp_("stereoModeProp","Stereoscopic Method")
    , anaglyphModeProp_("anaglyphModeProp","Anaglyph Color")
    , autostereoscopicModeProp_("autostereoscopicModeProp","Autostereoscopic Interleaving")
    , eyeInvertProp_("eyeInvertProp","Invert Eyes", false)
    , calibrateDisplayProp_("calibrateDisplay", "Show Calibration", false)
        //frustum settings
    , stereoAxisModeProp_("stereoAxisModeProp","Stereo Axis Mode")
    , eyeSeparationProp_("eyeSeparationProp","Eye Separation (cm)", 6.5f, 0.0f, 12.0f)
    , focalLengthProp_("focalLengthProp","Display Distance (cm)", 60.f, 10.f, 300.0f)
    , focalWidthProp_("focalWidthProp","Display Width (cm)", 45.f, 1.f, 200.0f)
    , relativeFocalLengthProp_("relativeFocalLengthProp","Focal Plane Factor", 0.5f, 0.001f, 1.0f)
    , useRealWorldFrustumProp_("useRealWorldFrustum", "Use Real-World Frustum Size", false)
        //events
    , mouseMoveEventProp_("mouseEvent.move", "Move Event", this, &StereoCanvasRenderer::mouseMove,
      tgt::MouseEvent::MOUSE_BUTTON_NONE, tgt::MouseEvent::MOTION | tgt::MouseEvent::CLICK | tgt::MouseEvent::ENTER_EXIT, tgt::MouseEvent::MODIFIER_NONE)
    //set member to default value
    , calibrationRightTexture_(0), calibrationLeftTexture_(0)
    , sideBySideShader_(0), anaglyphShader_(0), autostereoscopicShader_(0), copyTextureShader_(0)
    , previousCameraProjectionMode_(tgt::Camera::PERSPECTIVE)
    , nextExpectedImage_(NORMAL)
    , lastRunWasInInteractionMode_(false)
    , w8ingOn2Eye_(false)
    , fullScreenVertexArrayObjectID_(0)
    , fullScreenVertexBufferObjectID_(0)
{
    /*
    // private ports are not added here but in initialize() to avoid context problems
    addPrivateRenderPort(tmpPort_);
    addPrivateRenderPort(sideBySidePort_);
    addPrivateRenderPort(finalPort_);
    */

    //event handler for splitscreen
    addEventProperty(mouseMoveEventProp_);

    // stereoscopic method
    addProperty(stereoModeProp_);
        stereoModeProp_.addOption("nostereo","No Stereo",StereoCanvasRenderer::NO_STEREO_MODE);
        stereoModeProp_.addOption("splitscreen","Side by Side",StereoCanvasRenderer::SPLITSCREEN_STEREO_MODE);
        stereoModeProp_.addOption("autostereoscopic","Auto-Stereoscopic",StereoCanvasRenderer::AUTOSTEREOSCOPIC_STEREO_MODE);
        stereoModeProp_.addOption("anaglyph","Anaglyph",StereoCanvasRenderer::ANAGLYPH_STEREO_MODE);
        stereoModeProp_.addOption("quadbuffer","Quadbuffer",StereoCanvasRenderer::QUADBUFFER_STEREO_MODE);
        stereoModeProp_.onChange(MemberFunctionCallback<StereoCanvasRenderer>(this, &StereoCanvasRenderer::stereoModeOnChange));
        stereoModeProp_.setGroupID("stereoscopic method");
    addProperty(anaglyphModeProp_);
        anaglyphModeProp_.addOption("redcyan","Red - Cyan",StereoCanvasRenderer::RED_CYAN);
        anaglyphModeProp_.addOption("redblue","Red - Blue",StereoCanvasRenderer::RED_BLUE);
        anaglyphModeProp_.addOption("redgreen","Red - Green",StereoCanvasRenderer::RED_GREEN);
        anaglyphModeProp_.onChange(MemberFunctionCallback<StereoCanvasRenderer>(this, &StereoCanvasRenderer::invalidateFinalPort));
        anaglyphModeProp_.setGroupID("stereoscopic method");
    addProperty(autostereoscopicModeProp_);
        autostereoscopicModeProp_.addOption("vertical","Vertica-Interleaved",StereoCanvasRenderer::VERTICAL_INTERLEAVED);
        autostereoscopicModeProp_.addOption("horizontal","Horizontal-Interleaved",StereoCanvasRenderer::HORIZONTAL_INTERLEAVED);
        autostereoscopicModeProp_.addOption("checker","Checker-Interleaved",StereoCanvasRenderer::CHECKER_INTERLEAVED);
        autostereoscopicModeProp_.onChange(MemberFunctionCallback<StereoCanvasRenderer>(this, &StereoCanvasRenderer::invalidateFinalPort));
        autostereoscopicModeProp_.setGroupID("stereoscopic method");
    addProperty(eyeInvertProp_);
        eyeInvertProp_.onChange(MemberFunctionCallback<StereoCanvasRenderer>(this, &StereoCanvasRenderer::invalidateSideBySidePort));
        eyeInvertProp_.setGroupID("stereoscopic method");
    addProperty(calibrateDisplayProp_);
        calibrateDisplayProp_.onChange(MemberFunctionCallback<StereoCanvasRenderer>(this, &StereoCanvasRenderer::invalidateSideBySidePort));
        calibrateDisplayProp_.setGroupID("stereoscopic method");
    setPropertyGroupGuiName("stereoscopic method", "Stereoscopic Method");

    //frustum settings
    addProperty(stereoAxisModeProp_);
        stereoAxisModeProp_.addOption("onAxis", "On Axis", tgt::Camera::ON_AXIS);
        stereoAxisModeProp_.addOption("onAxisHMD","On Axis (HMD)",tgt::Camera::ON_AXIS_HMD);
        stereoAxisModeProp_.onChange(MemberFunctionCallback<StereoCanvasRenderer>(this, &StereoCanvasRenderer::stereoAxisModeOnChange));
        stereoAxisModeProp_.setGroupID("frustum settings");
    addProperty(eyeSeparationProp_);
        eyeSeparationProp_.onChange(MemberFunctionCallback<StereoCanvasRenderer>(this, &StereoCanvasRenderer::eyeSeparationOnChange));
        eyeSeparationProp_.setGroupID("frustum settings");
    addProperty(useRealWorldFrustumProp_);
        useRealWorldFrustumProp_.onChange(MemberFunctionCallback<StereoCanvasRenderer>(this, &StereoCanvasRenderer::useRealWorldFrustumPropOnChange));
        useRealWorldFrustumProp_.setGroupID("frustum settings");
    addProperty(focalLengthProp_);
        focalLengthProp_.onChange(MemberFunctionCallback<StereoCanvasRenderer>(this, &StereoCanvasRenderer::focalLengthOnChange));
        focalLengthProp_.setGroupID("frustum settings");
    addProperty(focalWidthProp_);
        focalWidthProp_.onChange(MemberFunctionCallback<StereoCanvasRenderer>(this, &StereoCanvasRenderer::focalWidthOnChange));
        focalWidthProp_.setGroupID("frustum settings");
    addProperty(relativeFocalLengthProp_);
        relativeFocalLengthProp_.onChange(MemberFunctionCallback<StereoCanvasRenderer>(this, &StereoCanvasRenderer::relativeFocalLengthPropOnChange));
        relativeFocalLengthProp_.setNumDecimals(3);
        relativeFocalLengthProp_.setGroupID("frustum settings");
    setPropertyGroupGuiName("frustum settings", "Frustum Settings");

        // camera settings
    addProperty(cameraProp_);

    focalLengthOnChange();
    focalWidthOnChange();
    relativeFocalLengthPropOnChange();
    useRealWorldFrustumPropOnChange();
}

StereoCanvasRenderer::~StereoCanvasRenderer() {
}

void StereoCanvasRenderer::initialize() {
    CanvasRenderer::initialize();

    tgt::GLConditionalContextStateGuard guard(canvas_ != nullptr, canvas_);

    // port settings (normal render port added by CanvasRenderer) ----------------
    // evil hack: we add the private ports here to avoid context problems, since we now have set our own context
    // since we have called CanvasRenderer::initialize() first (to obtain our context),
    // which calls RenderProcessor::initialize() that initializes private render ports, we have to do this manually now
    addPrivateRenderPort(tmpPort_);
    tmpPort_.initialize();
    addPrivateRenderPort(sideBySidePort_);
    sideBySidePort_.initialize();
    addPrivateRenderPort(finalPort_);
    finalPort_.initialize();
    adjustRenderOutportSizes();
    // ------- finished port settings

    // create geometry for fullscreen triangle
    glGenVertexArrays(1, &fullScreenVertexArrayObjectID_);
    glGenBuffers(1, &fullScreenVertexBufferObjectID_);
    float defaultFullScreenBuffer[] = {
        -1.f, -1.f, //vertex 1
         0.f,  0.f, //tex 1
         3.f, -1.f, //vertex 2
         2.f,  0.f, //tex 2
        -1.f,  3.f, //vertex 3
         0.f,  2.f  //tex 3
    };
    glBindVertexArray(fullScreenVertexArrayObjectID_);
    glBindBuffer(GL_ARRAY_BUFFER, fullScreenVertexBufferObjectID_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(defaultFullScreenBuffer), defaultFullScreenBuffer, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, false, sizeof(float)*4, 0); // vertex position
    glVertexAttribPointer(1, 2, GL_FLOAT, false, sizeof(float)*4, (char*)NULL + sizeof(float)*2); // texture position
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    //load shaders
    sideBySideShader_ = ShdrMgr.loadSeparate("stereo.vert", "copyimagesplitscreen.frag", generateHeader(), false);
    anaglyphShader_ = ShdrMgr.loadSeparate("stereo.vert", "anaglyph.frag", generateHeader(), false);
    autostereoscopicShader_ = ShdrMgr.loadSeparate("stereo.vert", "autostereoscopic.frag", generateHeader(), false);
    copyTextureShader_ = ShdrMgr.loadSeparate("stereo.vert", "copytexture.frag", generateHeader(), false);
    //load textures
    calibrationRightTexture_ = TexMgr.load(VoreenApplication::app()->getCoreResourcePath("textures/stereocalibrationR.png"));
    calibrationLeftTexture_ = TexMgr.load(VoreenApplication::app()->getCoreResourcePath("textures/stereocalibrationL.png"));
    //prepare camera (set ProjectionMode to FRUSTUM)
    tgt::Camera* cam = const_cast<tgt::Camera*>(&cameraProp_.get());
    previousCameraProjectionMode_ = cam->getProjectionMode();
    cam->setProjectionMode(tgt::Camera::FRUSTUM);
    cameraProp_.invalidate();
    //get properties right
    stereoModeOnChange();
}

void StereoCanvasRenderer::deinitialize() {
    tgt::GLConditionalContextStateGuard guard(canvas_ != nullptr, canvas_);
    
    //delete buffers
    glDeleteBuffers(1, &fullScreenVertexBufferObjectID_);
    glDeleteVertexArrays(1, &fullScreenVertexArrayObjectID_);
    //delete shaders
    ShdrMgr.dispose(sideBySideShader_);
    sideBySideShader_ = 0;
    ShdrMgr.dispose(anaglyphShader_);
    anaglyphShader_ = 0;
    ShdrMgr.dispose(autostereoscopicShader_);
    autostereoscopicShader_ = 0;
    ShdrMgr.dispose(copyTextureShader_);
    copyTextureShader_ = 0;
    //delete textures
    if (calibrationRightTexture_)
       TexMgr.dispose(calibrationRightTexture_);
    calibrationRightTexture_ = 0;
    if (calibrationLeftTexture_)
       TexMgr.dispose(calibrationLeftTexture_);
    calibrationLeftTexture_ = 0;
    //set camera projection mode back to previous mode
    tgt::Camera* cam = const_cast<tgt::Camera*>(&cameraProp_.get());
    cam->setProjectionMode(previousCameraProjectionMode_);

    CanvasRenderer::deinitialize();
}

void StereoCanvasRenderer::process() {
    if (!canvas_)
        return;

    // Context sensitive code.
    tgt::GLContextStateGuard guard(canvas_);

    //process based on stereo mode
    switch(stereoModeProp_.getValue()){
    case NO_STEREO_MODE:
        //set camera eye to middle and process like a normal canvas
        if(cameraProp_.setStereoEyeMode(tgt::Camera::EYE_MIDDLE))
            cameraProp_.invalidate();
        CanvasRenderer::process();
        break;
    case SPLITSCREEN_STEREO_MODE:
    case ANAGLYPH_STEREO_MODE:
    case AUTOSTEREOSCOPIC_STEREO_MODE:
    case QUADBUFFER_STEREO_MODE:
        if(calibrateDisplayProp_.get()){ // show calibration textures
            if (!calibrationRightTexture_ || !calibrationLeftTexture_) {
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
                return;
            }
            copyIntoSideBySide(calibrationLeftTexture_,calibrationRightTexture_);
            copyIntoCanvas();
        } else {
            if (inport_.isReady()) {
                if(inport_.hasChanged() || !sideBySidePort_.hasValidResult())
                    processStereo();
                else {
                    if(!finalPort_.hasValidResult())
                        copyIntoFinal();
                    copyIntoCanvas();
                }
            } else { // not ready (show error texture)
                if (!errorTex_) {
                    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
                    return;
                }
                copyIntoSideBySide(errorTex_,errorTex_);
                copyIntoCanvas();
            }
        }
        break;
    default:
        LERROR("StereoCanvas: Unknown StereoMode!!!");
        break;
    }
#ifdef VRN_MODULE_HEADTRACKING
     //head tracking
    getTrackingUpdate();
#endif
    LGL_ERROR;
}

//------------------------------------------------------------------------------------------------------------------
//      resize functions
//------------------------------------------------------------------------------------------------------------------
void StereoCanvasRenderer::resizeAllPorts(tgt::ivec2 newsize) {
    //check, if the call came from inside or outside i.e. context change or not
    tgt::GLContextStateGuard guard(canvas_);
    tgt::ivec2 halfNewSize(newsize.x/2,newsize.y);
    tgt::ivec2 doubleNewSize(newsize.x*2,newsize.y);
    switch(stereoModeProp_.getValue()){
    case NO_STEREO_MODE:
        inport_.requestSize(newsize);
        tmpPort_.resize(newsize);
        sideBySidePort_.resize(newsize);
        break;
    case ANAGLYPH_STEREO_MODE:
    case AUTOSTEREOSCOPIC_STEREO_MODE:
    case QUADBUFFER_STEREO_MODE:
        inport_.requestSize(newsize);
        tmpPort_.resize(newsize);
        sideBySidePort_.resize(doubleNewSize);
        break;
    case SPLITSCREEN_STEREO_MODE:
        inport_.requestSize(halfNewSize);
        tmpPort_.resize(halfNewSize);
        sideBySidePort_.resize(newsize);
        break;
    default:
        break;
    }
    finalPort_.resize(newsize);
    //reset camera
    nextExpectedImage_ = NORMAL;
    invalidate();
}

void StereoCanvasRenderer::canvasResized(tgt::ivec2 newsize) {
    if (canvas_) {
        tgt::GLContextStateGuard guard(canvas_);
        resizeAllPorts(newsize);
        canvasSize_.set(newsize);
    }
}

void StereoCanvasRenderer::resizeCanvas(tgt::ivec2 newsize) {
    if (!tgt::hand(tgt::greaterThanEqual(newsize, tgt::ivec2(canvasSize_.getMinValue()))) && tgt::hand(tgt::lessThanEqual(newsize, canvasSize_.getMaxValue()))) {
        LWARNING("Invalid canvas dimensions: " << newsize << ". Ignoring.");
    }

    if (!canvas_)
        return;

    if (getProcessorWidget() && getProcessorWidget()->getSize() != newsize) {
        if (getProcessorWidget()->isVisible()) {
            getProcessorWidget()->setSize(newsize.x, newsize.y);
        }
        newsize = getProcessorWidget()->getSize();
    }

    if (newsize != inport_.getSize()) {
        //check, if the call came from inside or outside i.e. context change or not
        tgt::GLContextStateGuard guard(canvas_);
        glViewport(0, 0, static_cast<GLint>(newsize.x), static_cast<GLint>(newsize.y));
        resizeAllPorts(newsize);
    }
    canvasSize_.set(newsize);
}

void StereoCanvasRenderer::setCanvas(tgt::GLCanvas* canvas) {
    if (canvas == canvas_)
        return;

    //remove from old canvas:
    if (canvas_) {
        tgt::EventHandler* eh = canvas_->getEventHandler();
        if (eh) {
            eh->removeListener(this);
        }
    }
    canvas_ = canvas;
    //register at new canvas:
    if (canvas_) {
        tgt::EventHandler* eh = canvas_->getEventHandler();
        if (eh) {
            eh->addListenerToFront(this);
        }
        resizeAllPorts(canvas->getSize());
    }
    invalidate();
}

//------------------------------------------------------------------------------------------------------------------
//      get/render to texture functions
//------------------------------------------------------------------------------------------------------------------
const tgt::Texture* StereoCanvasRenderer::getImageColorTexture() const {
    switch(stereoModeProp_.getValue()){
        case NO_STEREO_MODE:
            return CanvasRenderer::getImageColorTexture();
            break;
        case ANAGLYPH_STEREO_MODE:
        case AUTOSTEREOSCOPIC_STEREO_MODE:
        case SPLITSCREEN_STEREO_MODE:
        case QUADBUFFER_STEREO_MODE:
            if (finalPort_.hasRenderTarget())
                return finalPort_.getColorTexture();
            else
                return 0;
            break;
        default:
            tgtAssert(false,"Unknown StereoMode!");
            return 0;
    }
}

tgt::Texture* StereoCanvasRenderer::getImageColorTexture(){
    return const_cast<tgt::Texture*>(static_cast<const StereoCanvasRenderer*>(this)->getImageColorTexture());
}

const tgt::Texture* StereoCanvasRenderer::getImageDepthTexture() const{
    switch(stereoModeProp_.getValue()){
        case NO_STEREO_MODE:
            return CanvasRenderer::getImageDepthTexture();
            break;
        case ANAGLYPH_STEREO_MODE:
        case AUTOSTEREOSCOPIC_STEREO_MODE:
        case SPLITSCREEN_STEREO_MODE:
        case QUADBUFFER_STEREO_MODE:
            if (finalPort_.hasRenderTarget())
                return finalPort_.getDepthTexture();
            else
                return 0;
            break;
            break;
        default:
            tgtAssert(false,"Unknown StereoMode!");
            return 0;
    }
}

tgt::Texture* StereoCanvasRenderer::getImageDepthTexture() {
    return const_cast<tgt::Texture*>(static_cast<const StereoCanvasRenderer*>(this)->getImageDepthTexture());
}

bool StereoCanvasRenderer::renderToImage(const std::string &filename) {
    switch(stereoModeProp_.getValue()){
        case NO_STEREO_MODE:
            return CanvasRenderer::renderToImage(filename);
            break;
        case SPLITSCREEN_STEREO_MODE:
        case ANAGLYPH_STEREO_MODE:
        case AUTOSTEREOSCOPIC_STEREO_MODE:
        case QUADBUFFER_STEREO_MODE:
            if (!canvas_) {
                LWARNING("StereoCanvasRenderer::renderToImage(): no canvas assigned");
                return false;
            }

            if (!finalPort_.hasRenderTarget()) {
                LWARNING("StereoCanvasRenderer::renderToImage(): storagePort has no data");
                return false;
            }

            renderToImageError_.clear();

            try {
                finalPort_.saveToImage(filename);
                LINFO("Saved rendering " << finalPort_.getSize() << " to file: " << tgt::FileSystem::cleanupPath(renderToImageFilename_));
            }
            catch (std::bad_alloc& /*e*/) {
                LERROR("Exception in StereoCanvasRenderer::renderToImage(): bad allocation (" << getID() << ")");
                renderToImageError_ = "Not enough system memory (bad allocation)";
            }
            catch (VoreenException& e) {
                LERROR(e.what());
                renderToImageError_ = std::string(e.what());
            }
            catch (std::exception& e) {
                LERROR("Exception in StereoCanvasRenderer::renderToImage(): " << e.what() << " (" << getID() << ")");
                renderToImageError_ = std::string(e.what());
            }

            return (renderToImageError_.empty());
            break;
        default:
            tgtAssert(false,"Unknown StereoMode!");
            return false;
    }
}
bool StereoCanvasRenderer::renderToImage(const std::string &filename, tgt::ivec2 dimensions){
    bool success; tgt::ivec2 oldDimensions;
    switch(stereoModeProp_.getValue()){
        case NO_STEREO_MODE:
            return CanvasRenderer::renderToImage(filename,dimensions);
            break;
        case SPLITSCREEN_STEREO_MODE:
        case ANAGLYPH_STEREO_MODE:
        case AUTOSTEREOSCOPIC_STEREO_MODE:
        case QUADBUFFER_STEREO_MODE:
        {
            if (!canvas_) {
                LWARNING("StereoCanvasRenderer::renderToImage(): no canvas assigned");
                return false;
            }

            if (!finalPort_.hasRenderTarget()) {
                LWARNING("StereoCanvasRenderer::renderToImage(): storagePort has no data");
                return false;
            }

            oldDimensions = finalPort_.getSize();
            // resize texture container to desired image dimensions and propagate change
            tgt::GLContextStateGuard guard(canvas_);
            resizeAllPorts(dimensions);
            canvas_->repaint();
            canvas_->repaint();
            canvas_->repaint();
            // render with adjusted viewport size
            success = renderToImage(filename);
            // reset texture container dimensions from canvas size
            resizeAllPorts(oldDimensions);
            return success;
            break;
        }
        default:
            tgtAssert(false,"Unknown StereoMode!");
            return false;
    }
}

//------------------------------------------------------------------------------------------------------------------
//      on change and event functions
//------------------------------------------------------------------------------------------------------------------
void StereoCanvasRenderer::stereoModeOnChange() {
    nextExpectedImage_ = NORMAL;
    canvasResized(canvasSize_.get());
    //switch based on StereoMode
    switch(stereoModeProp_.getValue()) {
    case NO_STEREO_MODE:
        anaglyphModeProp_.setVisibleFlag(false);
        autostereoscopicModeProp_.setVisibleFlag(false);
        stereoAxisModeProp_.setReadOnlyFlag(true);
        eyeInvertProp_.set(false);
        eyeInvertProp_.setReadOnlyFlag(true);
        calibrateDisplayProp_.set(false);
        calibrateDisplayProp_.setReadOnlyFlag(true);
        eyeSeparationProp_.setReadOnlyFlag(true);
        focalLengthProp_.setReadOnlyFlag(true);
        focalWidthProp_.setReadOnlyFlag(true);
        relativeFocalLengthProp_.setReadOnlyFlag(true);
        useRealWorldFrustumProp_.setReadOnlyFlag(true);
        stereoAxisModeProp_.setReadOnlyFlag(true);
        break;
    case QUADBUFFER_STEREO_MODE:
    case SPLITSCREEN_STEREO_MODE:
        anaglyphModeProp_.setVisibleFlag(false);
        autostereoscopicModeProp_.setVisibleFlag(false);
        stereoAxisModeProp_.setReadOnlyFlag(false);
        eyeInvertProp_.setReadOnlyFlag(false);
        calibrateDisplayProp_.setReadOnlyFlag(false);
        eyeSeparationProp_.setReadOnlyFlag(false);
        useRealWorldFrustumProp_.setReadOnlyFlag(false);
        focalLengthProp_.setReadOnlyFlag(!useRealWorldFrustumProp_.get());
        focalWidthProp_.setReadOnlyFlag(!useRealWorldFrustumProp_.get());
        relativeFocalLengthProp_.setReadOnlyFlag(useRealWorldFrustumProp_.get());
        stereoAxisModeProp_.setReadOnlyFlag(false);
        break;
    case AUTOSTEREOSCOPIC_STEREO_MODE:
        anaglyphModeProp_.setVisibleFlag(false);
        autostereoscopicModeProp_.setVisibleFlag(true);
        stereoAxisModeProp_.setReadOnlyFlag(false);
        eyeInvertProp_.setReadOnlyFlag(false);
        calibrateDisplayProp_.setReadOnlyFlag(false);
        eyeSeparationProp_.setReadOnlyFlag(false);
        useRealWorldFrustumProp_.setReadOnlyFlag(false);
        focalLengthProp_.setReadOnlyFlag(!useRealWorldFrustumProp_.get());
        focalWidthProp_.setReadOnlyFlag(!useRealWorldFrustumProp_.get());
        relativeFocalLengthProp_.setReadOnlyFlag(!!useRealWorldFrustumProp_.get());
        stereoAxisModeProp_.setReadOnlyFlag(false);
        break;
    case ANAGLYPH_STEREO_MODE:
        anaglyphModeProp_.setVisibleFlag(true);
        autostereoscopicModeProp_.setVisibleFlag(false);
        stereoAxisModeProp_.setReadOnlyFlag(false);
        eyeInvertProp_.setReadOnlyFlag(false);
        calibrateDisplayProp_.setReadOnlyFlag(false);
        eyeSeparationProp_.setReadOnlyFlag(false);
        useRealWorldFrustumProp_.setReadOnlyFlag(false);
        focalLengthProp_.setReadOnlyFlag(!useRealWorldFrustumProp_.get());
        focalWidthProp_.setReadOnlyFlag(!useRealWorldFrustumProp_.get());
        relativeFocalLengthProp_.setReadOnlyFlag(!!useRealWorldFrustumProp_.get());
        stereoAxisModeProp_.setReadOnlyFlag(false);
        break;
    default:
        tgtAssert(false,"Unknown StereoMode!");
        break;
    }
    finalPort_.invalidateResult(); //force new rendering
}

void StereoCanvasRenderer::invalidateSideBySidePort() {
    sideBySidePort_.invalidateResult();
    invalidate();
}

void StereoCanvasRenderer::invalidateFinalPort() {
    finalPort_.invalidateResult();
    invalidate();
}

void StereoCanvasRenderer::eyeSeparationOnChange() {
    // convert to mm
    if(cameraProp_.setStereoEyeSeparation(eyeSeparationProp_.get() * 10.f))
        cameraProp_.invalidate();
}

void StereoCanvasRenderer::stereoAxisModeOnChange() {
    if(cameraProp_.setStereoAxisMode(stereoAxisModeProp_.getValue()))
        cameraProp_.invalidate();
}

void StereoCanvasRenderer::focalLengthOnChange() {
    // convert to mm
    if(cameraProp_.setStereoFocalLength(focalLengthProp_.get() * 10.f))
        cameraProp_.invalidate();
}

void StereoCanvasRenderer::focalWidthOnChange() {
    // convert to mm
    if(cameraProp_.setStereoWidth(focalWidthProp_.get() * 10.f))
        cameraProp_.invalidate();
}

void StereoCanvasRenderer::relativeFocalLengthPropOnChange() {
    float foc = relativeFocalLengthProp_.get();
    float maxDist = relativeFocalLengthProp_.getMaxValue();
    float minDist = relativeFocalLengthProp_.getMinValue();
    float maxLog = log(maxDist);
    float minLog = log(minDist);
    float scale = (maxLog - minLog) / (maxDist - minDist);
    foc = std::exp(minLog + scale * (foc - minDist));
    if(cameraProp_.setStereoRelativeFocalLength(foc))
        cameraProp_.invalidate();
}

void StereoCanvasRenderer::useRealWorldFrustumPropOnChange() {
    focalLengthProp_.setReadOnlyFlag(!useRealWorldFrustumProp_.get());
    focalWidthProp_.setReadOnlyFlag(!useRealWorldFrustumProp_.get());
    relativeFocalLengthProp_.setReadOnlyFlag(!!useRealWorldFrustumProp_.get());
    if(cameraProp_.setUseRealWorldFrustum(useRealWorldFrustumProp_.get()))
        cameraProp_.invalidate();
}

void StereoCanvasRenderer::onEvent(tgt::Event* e) {
    if (!canvas_)
        return;

    //canvas_->getGLFocus(); events must be in shared context, not canvas context
    tgt::MouseEvent* me = dynamic_cast<tgt::MouseEvent*>(e);
    //pass, if no mouseevent and stereo mode is split screen
    if (!me || mouseMoveEventProp_.accepts(me) || stereoModeProp_.getValue() != SPLITSCREEN_STEREO_MODE) {
        RenderProcessor::onEvent(e);
        return;
    }

    //set right viewport for mouse events in split screen mode
    tgt::ivec2 view = me->viewport();
    switch(stereoModeProp_.getValue()){
        case NO_STEREO_MODE:
        case AUTOSTEREOSCOPIC_STEREO_MODE:
        case ANAGLYPH_STEREO_MODE:
        case QUADBUFFER_STEREO_MODE:
            //should not get here
            tgtAssert(false,"StereoCanvasRenderer::onEvent: Unexpected stereo mode");
            break;
        case SPLITSCREEN_STEREO_MODE:
            view.x /= 2;
            if (me->x() < (me->viewport().x / 2)) {
                tgt::MouseEvent newme(me->x(), me->y(), me->action(), me->modifiers(), me->button(), view);
                newme.ignore();  // accepted is set to true by default
                inport_.distributeEvent(&newme);
                if (newme.isAccepted())
                    me->accept();
            }
            else {
                tgt::MouseEvent newme(me->x() - (me->viewport().x / 2), me->y(), me->action(), me->modifiers(), me->button(), view);
                newme.ignore();  // accepted is set to true by default
                inport_.distributeEvent(&newme);
                if (newme.isAccepted())
                    me->accept();
            }
            break;
        default:
            //should not get here
            tgtAssert(false,"StereoCanvasRenderer::onEvent: Unexpected stereo mode");
            break;
    }
}

void StereoCanvasRenderer::mouseMove(tgt::MouseEvent* e) {
    tgt::ivec2 view;
    tgt::MouseEvent* me;
    switch(stereoModeProp_.getValue()){
        case NO_STEREO_MODE:
        case ANAGLYPH_STEREO_MODE:
        case AUTOSTEREOSCOPIC_STEREO_MODE:
        case QUADBUFFER_STEREO_MODE:
            break;
        case SPLITSCREEN_STEREO_MODE:
            e->accept();
            view = e->viewport();
            view.x /= 2;
            me = new tgt::MouseEvent(e->x() % (e->viewport().x/2), e->y() , tgt::MouseEvent::MOTION, e->modifiers(), e->button(), view);
            me->ignore();
            inport_.distributeEvent(me);
            if(me->isAccepted())
                e->accept();
            delete me;
            break;
        default:
            break;
    }
}

//------------------------------------------------------------------------------------------------------------------
//      copy functions
//------------------------------------------------------------------------------------------------------------------
void StereoCanvasRenderer::processStereo() {
    if(nextExpectedImage_ == NORMAL){
        nextExpectedImage_ = LEFT_1;
        if (cameraProp_.setStereoEyeMode(tgt::Camera::EYE_LEFT))
            cameraProp_.invalidate();
    }
    else{
        switch(nextExpectedImage_){
        case LEFT_1:
            copyIntoCanvas();
            copyIntoPort(&inport_,&tmpPort_);
            lastRunWasInInteractionMode_ = QualityMode.isInteractionMode();
            nextExpectedImage_ = RIGHT_1;
            if (cameraProp_.setStereoEyeMode(tgt::Camera::EYE_RIGHT))
                cameraProp_.invalidate();
            w8ingOn2Eye_ = true;
            break;
        case RIGHT_1:
            if(QualityMode.isInteractionMode() == lastRunWasInInteractionMode_){
                copyIntoSideBySide(&tmpPort_,&inport_);
                copyIntoCanvas();
                w8ingOn2Eye_ = false;
                nextExpectedImage_ = RIGHT_2;
            }
            else { //interactionMode toggeled
                copyIntoPort(&inport_,&tmpPort_);
                lastRunWasInInteractionMode_ = QualityMode.isInteractionMode();
                nextExpectedImage_ = LEFT_2;
                copyIntoCanvas();
                if(cameraProp_.setStereoEyeMode(tgt::Camera::EYE_LEFT))
                    cameraProp_.invalidate();
                w8ingOn2Eye_ = true;
            }
            break;
        case RIGHT_2:
            copyIntoCanvas();
            copyIntoPort(&inport_,&tmpPort_);
            lastRunWasInInteractionMode_ = QualityMode.isInteractionMode();
            nextExpectedImage_ = LEFT_2;
            if(cameraProp_.setStereoEyeMode(tgt::Camera::EYE_LEFT))
                cameraProp_.invalidate();
            w8ingOn2Eye_ = true;
            break;
        case LEFT_2:
            if (QualityMode.isInteractionMode() == lastRunWasInInteractionMode_) {
                copyIntoSideBySide(&inport_,&tmpPort_);
                copyIntoCanvas();
                nextExpectedImage_ = LEFT_1;
                w8ingOn2Eye_ = false;
            }
            else { //interactionMode toggeled
                copyIntoPort(&inport_,&tmpPort_);
                lastRunWasInInteractionMode_ = QualityMode.isInteractionMode();
                nextExpectedImage_ = RIGHT_1;
                copyIntoCanvas();
                if(cameraProp_.setStereoEyeMode(tgt::Camera::EYE_RIGHT))
                    cameraProp_.invalidate();
                w8ingOn2Eye_ = true;
            }
            break;
        default:
            LERROR("StereoCanvas: Unknown NextExpectedImageMode!!!");
            break;
        }
    }
}

void StereoCanvasRenderer::copyIntoPort(RenderPort* input, RenderPort* output){
    output->activateTarget();
        // activate shader
        copyTextureShader_->activate();
        // bind input textures
        input->bindTextures(GL_TEXTURE0, GL_TEXTURE1);
        // pass texture parameters to the shader
        copyTextureShader_->setUniform("colorTex_", 0);
        copyTextureShader_->setUniform("depthTex_", 1);
        copyTextureShader_->setUniform("useQuadBuffer_", false);
        LGL_ERROR;
        // execute the shader
        renderFullScreen();
        copyTextureShader_->deactivate();
    output->deactivateTarget();
    glActiveTexture(GL_TEXTURE0); //default voreen settings
    LGL_ERROR;
}

void StereoCanvasRenderer::copyIntoSideBySide(tgt::Texture* colorLeft, tgt::Texture* colorRight, tgt::Texture* depthLeft, tgt::Texture* depthRight){
    //invert eye textures
    if (eyeInvertProp_.get()){
        tgt::Texture* help = colorLeft;
        colorLeft = colorRight;
        colorRight = help;
        help = depthLeft;
        depthLeft = depthRight;
        depthRight = help;
    }

    sideBySidePort_.activateTarget();

    // activate shader
    sideBySideShader_->activate();

    // pass texture parameters to the shader
    glActiveTexture(GL_TEXTURE0);
    colorLeft->bind();
    sideBySideShader_->setUniform("colorTexLeft_", 0);
    if(depthLeft) {
        glActiveTexture(GL_TEXTURE1);
        depthLeft->bind();
        sideBySideShader_->setUniform("depthTexLeft_", 1);
        sideBySideShader_->setUniform("useDepthTexLeft_",true);
    } else {
        sideBySideShader_->setUniform("useDepthTexLeft_",false);
    }
    glActiveTexture(GL_TEXTURE2);
    colorRight->bind();
    sideBySideShader_->setUniform("colorTexRight_", 2);
    if(depthRight) {
        glActiveTexture(GL_TEXTURE3);
        depthRight->bind();
        sideBySideShader_->setUniform("depthTexRight_", 3);
        sideBySideShader_->setUniform("useDepthTexRight_",true);
    } else {
        sideBySideShader_->setUniform("useDepthTexRight_",false);
    }

    renderFullScreen();

    sideBySideShader_->deactivate();
    sideBySidePort_.deactivateTarget();
    sideBySidePort_.validateResult();

    glActiveTexture(GL_TEXTURE0); //default voreen settings
    LGL_ERROR;
    //update final result
    copyIntoFinal();
}

void StereoCanvasRenderer::copyIntoSideBySide(RenderPort* left, RenderPort* right){
    copyIntoSideBySide(left->getColorTexture(), right->getColorTexture(),
                    left->getDepthTexture(), right->getDepthTexture());
}

void StereoCanvasRenderer::copyIntoFinal() {
    finalPort_.activateTarget();
    //switch based on StereoMode
    switch(stereoModeProp_.getValue()) {
    case NO_STEREO_MODE:
        //shouldn't get here
        tgtAssert(false,"copyIntoCanvas does not support NO_STEREO_MODE!!!");
        finalPort_.deactivateTarget();
        break;
    case AUTOSTEREOSCOPIC_STEREO_MODE:
        renderAutostereoscopic();
        break;
    case ANAGLYPH_STEREO_MODE:
        renderAnaglyph();
        break;
    case SPLITSCREEN_STEREO_MODE:
    case QUADBUFFER_STEREO_MODE:
        renderSplitScreen();
        break;
    default:
        //shouldn't get here
        tgtAssert(false,"Unknown StereoMode!!!");
        break;
    }
    finalPort_.deactivateTarget();
    glActiveTexture(GL_TEXTURE0); //default voreen settings
    LGL_ERROR;
}

void StereoCanvasRenderer::copyIntoCanvas(){
    //copied from CanvasRenderer
    tgt::GLContextStateGuard guard(canvas_);

    glViewport(0, 0, canvas_->getSize().x, canvas_->getSize().y);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //switch based on StereoMode
    switch(stereoModeProp_.getValue()) {
    case NO_STEREO_MODE:
        //shouldn't get here
        tgtAssert(false,"copyIntoCanvas does not support NO_STEREO_MODE!!!");
        break;
    case AUTOSTEREOSCOPIC_STEREO_MODE:
    case ANAGLYPH_STEREO_MODE:
    case SPLITSCREEN_STEREO_MODE:
        // activate shader
        copyTextureShader_->activate();
        // manually pass the viewport dimensions to the shader
        finalPort_.bindTextures(GL_TEXTURE0, GL_TEXTURE1);
        copyTextureShader_->setUniform("colorTex_", 0);
        copyTextureShader_->setUniform("depthTex_", 1);
        copyTextureShader_->setUniform("useQuadBuffer_", false);
        LGL_ERROR;
        renderFullScreen();
        copyTextureShader_->deactivate();
        if (canvas_->getAutoFlush())
            canvas_->finalizeRendering();
        break;
    case QUADBUFFER_STEREO_MODE:
        //render left
        glDrawBuffer(GL_BACK_LEFT);
        copyTextureShader_->activate();
        finalPort_.bindTextures(GL_TEXTURE0, GL_TEXTURE1);
        copyTextureShader_->setUniform("colorTex_", 0);
        copyTextureShader_->setUniform("depthTex_", 1);
        copyTextureShader_->setUniform("useQuadBuffer_", true);
        copyTextureShader_->setUniform("useLeftQuadBuffer_", true);
        LGL_ERROR;
        renderFullScreen();
        //right buffer
        glDrawBuffer(GL_BACK_RIGHT);
        copyTextureShader_->setUniform("useLeftQuadBuffer_", false);
        renderFullScreen();
        copyTextureShader_->deactivate();
        //set back to normel
        glDrawBuffer(GL_BACK);
        if (canvas_->getAutoFlush())
            canvas_->finalizeRendering();
        break;
    default:
        //shouldn't get here
        tgtAssert(false,"Unknown StereoMode!!!");
        break;
    }

    glActiveTexture(GL_TEXTURE0); //default voreen settings
    LGL_ERROR;
}

//------------------------------------------------------------------------------------------------------------------
//      render functions
//------------------------------------------------------------------------------------------------------------------
void StereoCanvasRenderer::renderAnaglyph() {
    // activate shader
    anaglyphShader_->activate();
    // manually pass the viewport dimensions to the shader
    anaglyphShader_->setUniform("screenDimRCP_", tgt::vec2(1.f) / static_cast<tgt::vec2>(canvas_->getSize()));
    sideBySidePort_.bindTextures(GL_TEXTURE0, GL_TEXTURE1);
    anaglyphShader_->setUniform("colorTex_", 0);
    anaglyphShader_->setUniform("depthTex_", 1);
    anaglyphShader_->setUniform("colorCode_", static_cast<int>(anaglyphModeProp_.getValue()));
    LGL_ERROR;
    renderFullScreen();
    anaglyphShader_->deactivate();
}

void StereoCanvasRenderer::renderAutostereoscopic() {
    // activate shader
    autostereoscopicShader_->activate();
    // manually pass the viewport dimensions to the shader
    autostereoscopicShader_->setUniform("screenDimRCP_", tgt::vec2(1.f) / static_cast<tgt::vec2>(canvas_->getSize()));
    sideBySidePort_.bindTextures(GL_TEXTURE0, GL_TEXTURE1);
    autostereoscopicShader_->setUniform("colorTex_", 0);
    autostereoscopicShader_->setUniform("depthTex_", 1);
    autostereoscopicShader_->setUniform("interleaveCode_", static_cast<int>(autostereoscopicModeProp_.getValue()));
    LGL_ERROR;
    renderFullScreen();
    autostereoscopicShader_->deactivate();
}

void StereoCanvasRenderer::renderSplitScreen() {
    // activate shader
    copyTextureShader_->activate();
    sideBySidePort_.bindTextures(GL_TEXTURE0, GL_TEXTURE1);
    copyTextureShader_->setUniform("colorTex_", 0);
    copyTextureShader_->setUniform("depthTex_", 1);
    copyTextureShader_->setUniform("useQuadBuffer_", false);
    LGL_ERROR;
    renderFullScreen();
    copyTextureShader_->deactivate();
}

void StereoCanvasRenderer::renderFullScreen() {
    LGL_ERROR;
    glDepthFunc(GL_ALWAYS);
    LGL_ERROR;
    glBindVertexArray(fullScreenVertexArrayObjectID_);
    LGL_ERROR;
    glDrawArrays(GL_TRIANGLES, 0, 3);
    LGL_ERROR;
    glDepthFunc(GL_LESS);
}

} // namespace voreen
