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

#include "canvasrenderer.h"

#include "voreen/core/utils/voreenpainter.h"
#include "voreen/core/utils/stringutils.h"
#include "voreen/core/processors/processorwidget.h"
#include "voreen/core/processors/processorwidgetfactory.h"
#include "voreen/core/voreenapplication.h"

#include "tgt/glcanvas.h"
#include "tgt/texturemanager.h"
#include "tgt/gpucapabilities.h"
#include "tgt/filesystem.h"
#include "tgt/glcontextmanager.h"

#ifdef VRN_MODULE_DEVIL
#include "modules/devil/devilmodule.h"
#endif

#include "voreen/core/utils/voreenqualitymode.h"

namespace voreen {

const std::string CanvasRenderer::loggerCat_("voreen.core.CanvasRenderer");

CanvasRenderer::CanvasRenderer()
    : RenderProcessor()
    , showCanvas_("showCanvas", "Show Canvas", true)
    , canvasSize_("canvasSize", "Canvas Size", tgt::ivec2(256), tgt::ivec2(32), tgt::ivec2(2 << 12), Processor::VALID)
    , showCursor_("showCursor", "Show Cursor", true)
    , showFullScreen_("showFullScreen", "Fullscreen (F11)", false)
    , screenshotFilename_("screenshotFilename", "File", "Select file...","", "*.*", FileDialogProperty::SAVE_FILE, Processor::INVALID_PATH)
    , saveScreenshotButton_("saveScreenshot", "Save Screenshot")
    , canvas_(0)
    , shader_(0)
    , inport_(Port::INPORT, "image.input", "Image Input", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN)
    , errorTex_(0)
    , renderToImage_(false)
    , renderToImageFilename_("")
{
    addPort(inport_);
    addProperty(canvasSize_);

    addProperty(showCanvas_);
    addProperty(showCursor_);
    addProperty(showFullScreen_);
    showCanvas_.setGroupID("set");
    showCursor_.setGroupID("set");
    showFullScreen_.setGroupID("set");
    setPropertyGroupGuiName("set", "Settings");

    addProperty(screenshotFilename_);
    addProperty(saveScreenshotButton_);
    screenshotFilename_.setGroupID("screenshot");
    saveScreenshotButton_.setGroupID("screenshot");
    setPropertyGroupGuiName("screenshot", "Screenshot");

    canvasSize_.onChange(MemberFunctionCallback<CanvasRenderer>(this, &CanvasRenderer::sizePropChanged));

    showCanvas_.onChange(MemberFunctionCallback<CanvasRenderer>(this, &CanvasRenderer::boolPropertyChanged));
    showCursor_.onChange(MemberFunctionCallback<CanvasRenderer>(this, &CanvasRenderer::boolPropertyChanged));
    showFullScreen_.onChange(MemberFunctionCallback<CanvasRenderer>(this, &CanvasRenderer::boolPropertyChanged));
    saveScreenshotButton_.onChange(MemberFunctionCallback<CanvasRenderer>(this, &CanvasRenderer::saveScreenshotClicked));
}

CanvasRenderer::~CanvasRenderer() {}

Processor* CanvasRenderer::create() const {
    return new CanvasRenderer();
}

void CanvasRenderer::process() {
    if (!canvas_)
        return;

    tgt::GLContextStateGuard guard(canvas_);

    tgt::ivec2 screenSize = canvas_->getPhysicalSize();

    glViewport(0, 0, screenSize.x, screenSize.y);

    // render inport to image, if renderToImage flag has been set
    if (renderToImage_) {
        if (inport_.isReady()) {
            // render inport to image, if renderToImage flag has been set
            try {
                renderInportToImage(renderToImageFilename_);
                LINFO("Saved rendering with dimensions " << inport_.getSize() << " to file: " << tgt::FileSystem::cleanupPath(renderToImageFilename_));
            }
            catch (std::bad_alloc& /*e*/) {
                LERROR("Exception in CanvasRenderer::renderInportToImage(): bad allocation (" << getID() << ")");
                renderToImageError_ = "Not enough system memory (bad allocation)";
            }
            catch (VoreenException& e) {
                LERROR(e.what());
                renderToImageError_ = std::string(e.what());
            }
            catch (std::exception& e) {
                LERROR("Exception in CanvasRenderer::renderInportToImage(): " << e.what() << " (" << getID() << ")");
                renderToImageError_ = std::string(e.what());
            }
            renderToImage_ = false;
        }
    }
    else {
        // activate shader
        shader_->activate();

        // set common uniforms
        setGlobalShaderParameters(shader_);

        // manually pass the viewport dimensions to the shader,
        // since setGlobalShaderParameters() expects a render outport, which we do not have
        shader_->setIgnoreUniformLocationError(true);
        shader_->setUniform("screenDim_", tgt::vec2(screenSize));
        shader_->setUniform("screenDimRCP_", 1.f / tgt::vec2(screenSize));
        shader_->setUniform("useTexcoords", false);
        // pass texture parameters to the shader
        shader_->setUniform("colorTex_", 0);
        shader_->setUniform("depthTex_", 1);
        shader_->setIgnoreUniformLocationError(false);
        LGL_ERROR;

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        if (inport_.isReady()) {
            // bind input textures
            inport_.bindTextures(GL_TEXTURE0, GL_TEXTURE1);
            inport_.setTextureParameters(shader_, "texParams_");
            renderQuad();
            glActiveTexture(GL_TEXTURE0);
            LGL_ERROR;
        } else if (errorTex_) {
            glActiveTexture(GL_TEXTURE0);
            errorTex_->bind();
            errorTex_->enable();
            // set texture parameter manually
            shader_->setIgnoreUniformLocationError(true);
            shader_->setUniform("texParams_.dimensions_", tgt::vec2(screenSize));
            shader_->setUniform("texParams_.dimensionsRCP_", tgt::vec2(1.f) / tgt::vec2(screenSize));
            shader_->setUniform("texParams_.matrix_", tgt::mat4::identity);
            shader_->setIgnoreUniformLocationError(false);

            renderQuad();
            LGL_ERROR;
            errorTex_->disable();
        }
        shader_->deactivate();
    }
}

bool CanvasRenderer::isReady() const {
    return true;
}

bool CanvasRenderer::isValid() const {
    return false;
}

void CanvasRenderer::initialize() {
    RenderProcessor::initialize();

#ifdef VRN_MODULE_DEVIL
    std::string imageExtensions = "*." + strJoin(DevILModule::getSupportedWriteExtensions(), ";; *.");
    screenshotFilename_.setFileFilter(imageExtensions);
#endif

    if (getProcessorWidget())
        getProcessorWidget()->updateFromProcessor();

    shader_ = ShdrMgr.loadSeparate("passthrough.vert", "copyimage.frag", generateHeader(), false);

    errorTex_ = TexMgr.load(VoreenApplication::app()->getCoreResourcePath("textures/error.tga"));
    canvasSize_.setMaxValue(tgt::ivec2(GpuCaps.getMaxTextureSize()));
}

void CanvasRenderer::deinitialize()  {
    ShdrMgr.dispose(shader_);
    shader_ = 0;
    if (errorTex_)
        TexMgr.dispose(errorTex_);
    errorTex_ = 0;
    LGL_ERROR;

    //NOTE: this function deletes the canvas
    RenderProcessor::deinitialize();
    //LGL_ERROR; no glerror check, since the active context has been deleted
    setCanvas(0);
}

bool CanvasRenderer::isEndProcessor() const {
    return true;
}

void CanvasRenderer::invalidate(int inv) {
    // We bypass Processor::invalidate() in order to prevent a needless double
    // processing of the network.
    PropertyOwner::invalidate(inv);
    //check, if is ready (isInitialized)
    bool ready = isReady();
    if(ready && (processorState_ == PROCESSOR_STATE_NOT_READY)){
        processorState_ = PROCESSOR_STATE_READY;
        notifyStateChanged();
    } else if(!ready && (processorState_ == PROCESSOR_STATE_READY)) {
        processorState_ = PROCESSOR_STATE_NOT_READY;
        notifyStateChanged();
    }

    if(canvas_)
        canvas_->update();

    //take screenshot, if path has been changed
    if(inv == Processor::INVALID_PATH && isInitialized())
        saveScreenshotClicked();
}

void CanvasRenderer::setCanvas(tgt::GLCanvas* canvas) {
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
        inport_.requestSize(canvas_->getPhysicalSize());
    }

    invalidate();
}

tgt::GLCanvas* CanvasRenderer::getCanvas() const {
    return canvas_;
}

bool CanvasRenderer::renderToImage(const std::string &filename) {
    if (!canvas_) {
        LWARNING("CanvasRenderer::renderToImage(): no canvas assigned");
        renderToImageError_ = "No canvas assigned";
        return false;
    }

    if (!inport_.hasRenderTarget()) {
        LWARNING("CanvasRenderer::renderToImage(): inport has no data");
        renderToImageError_ = "No rendering";
        return false;
    }

    // enable render-to-file on next process
    renderToImageFilename_ = filename;
    renderToImage_ = true;
    renderToImageError_.clear();

    // request high quality mode for screenshots
    VoreenQualityMode::getRef().requestQualityModeForced(VoreenQualityMode::RQ_HIGH, this);

    // force rendering pass
    canvas_->repaint();

    // set back quality mode
    VoreenQualityMode::getRef().requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);

    return (renderToImageError_.empty());
}

bool CanvasRenderer::renderToImage(const std::string &filename, tgt::ivec2 dimensions) {
    if (!canvas_) {
        LWARNING("CanvasRenderer::renderToImage(): no canvas assigned");
        renderToImageError_ = "No canvas assigned";
        return false;
    }

    if (!inport_.hasRenderTarget()) {
        LWARNING("CanvasRenderer::renderToImage(): inport has no data");
        renderToImageError_ = "No rendering";
        return false;
    }

    tgt::ivec2 oldDimensions = inport_.getSize();

    // resize texture container to desired image dimensions and propagate change
    if (oldDimensions != dimensions) {
        tgt::GLContextStateGuard guard(canvas_);
        inport_.requestSize(dimensions);
    }

    // render with adjusted viewport size
    bool success = renderToImage(filename);

    // reset texture container dimensions from canvas size
    if (oldDimensions != dimensions) {
        tgt::GLContextStateGuard guard(canvas_);
        inport_.requestSize(oldDimensions);
    }

    return success;
}

void CanvasRenderer::onEvent(tgt::Event* e) {
    tgt::GLContextStateGuard guard; // Activate Main Context.
    Processor::onEvent(e);
}

void CanvasRenderer::renderInportToImage(const std::string& filename) {
    if (!inport_.hasRenderTarget())
        throw VoreenException("No rendering");

    inport_.saveToImage(filename);
}

const tgt::Texture* CanvasRenderer::getImageColorTexture() const {
    if (inport_.hasRenderTarget())
        return inport_.getColorTexture();
    else
        return 0;
}

tgt::Texture* CanvasRenderer::getImageColorTexture() {
    if (inport_.hasRenderTarget())
        return inport_.getColorTexture();
    else
        return 0;
}

const tgt::Texture* CanvasRenderer::getImageDepthTexture() const {
    if (inport_.hasRenderTarget())
        return inport_.getDepthTexture();
    else
        return 0;
}

tgt::Texture* CanvasRenderer::getImageDepthTexture() {
    if (inport_.hasRenderTarget())
        return inport_.getDepthTexture();
    else
        return 0;
}

std::string CanvasRenderer::getRenderToImageError() const {
    return renderToImageError_;
}

tgt::ivec2 CanvasRenderer::getCanvasSize() const {
    return canvasSize_.get();
}

void CanvasRenderer::resizeCanvas(tgt::ivec2 newsize) {
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
        tgt::GLContextStateGuard guard(canvas_);
        glViewport(0, 0, static_cast<GLint>(newsize.x), static_cast<GLint>(newsize.y));
        inport_.requestSize(canvas_->getPhysicalSize());
    }

    canvasSize_.set(newsize);
}

void CanvasRenderer::sizePropChanged() {
    resizeCanvas(canvasSize_.get());
}


void CanvasRenderer::boolPropertyChanged() {
    if(getProcessorWidget())
        getProcessorWidget()->updateFromProcessor();
}

void CanvasRenderer::saveScreenshotClicked() {
    if (screenshotFilename_.get().empty()) {
        LWARNING("No screenshot filename set");
        return;
    }

    try {
        renderToImage(screenshotFilename_.get());
        //LINFO("Saved rendering with dimensions " << inport_.getSize() << " to file: " << screenshotFilename_.get());
    }
    catch (VoreenException& e) {
        LERROR("Failed to write screenshot to file '" << screenshotFilename_.get() << "': " << e.what());
    }

}

} // namespace voreen
