
/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2014 University of Muenster, Germany.                        *
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

#include "cminfooverlay.h"
#include "voreen/core/interaction/camerainteractionhandler.h"
#include "voreen/core/voreenapplication.h"

#include "tgt/texturemanager.h"
#include "tgt/gpucapabilities.h"
#include "tgt/glmath.h"
#include "tgt/textureunit.h"
#include "tgt/immediatemode/immediatemode.h"

#include "tgt/matrixstack.h"

namespace voreen {

static const std::string ICON_FILE = "custommodules/viscontest2015/resources/infoicon.png";

const std::string CMInfoOverlay::loggerCat_("voreen.CMInfoOverlay");

CMInfoOverlay::CMInfoOverlay()
    : ImageProcessor("image/compositor") //determines .frag shader
    //ports
    , inport_(Port::INPORT, "image.input", "Image Input")
    , outport_(Port::OUTPORT, "image.output", "Image Output")
    //basic
    , enableProp_("enableProp","Enable",true)
    //Info pic selection
    , infoImageFileProp_("infoImageFileProp", "Info image", "Select info image", "", "*.png", FileDialogProperty::OPEN_FILE)
    //position
    , positionProp_("positionProp", "Icon position", TOP_RIGHT)
    , sizeProp_("sizeProp", "Icon size", 0.08f, 0.0f, 1.f)
    , alphaProp_("alphaProp_", "Alpha for icon overlay", 0.5f, 0.0f, 1.f)
    , activatedProp_("activatedProp","Activated",false)
    , iconTexture_(nullptr)
    , infoTexture_(nullptr)
{
    //ports
    addPort(inport_);
    addPort(outport_);
    //basic
    addProperty(enableProp_);
    //image
    addProperty(infoImageFileProp_);
    //position
    addProperty(positionProp_);
    positionProp_.addOption("top_left", "Top left", TOP_LEFT);
    positionProp_.addOption("top_right", "Top right", TOP_RIGHT);
    positionProp_.addOption("bottom_left", "Bottom left", BOTTOM_LEFT);
    positionProp_.addOption("bottom_right", "Bottom right", BOTTOM_RIGHT);
    //size
    addProperty(sizeProp_);
    //addProperty(alphaProp_);

    addProperty(activatedProp_);

    ON_CHANGE(infoImageFileProp_, CMInfoOverlay, tryLoadInfoImageTex);
}

CMInfoOverlay::~CMInfoOverlay() {
}

Processor* CMInfoOverlay::create() const {
    return new CMInfoOverlay();
}

void CMInfoOverlay::initialize() {
    ImageProcessor::initialize();

    tryLoadIconTex();
}

void CMInfoOverlay::deinitialize() {
    //delete texture
    if(iconTexture_) {
        TexMgr.dispose(iconTexture_);
    }
    if(infoTexture_) {
        TexMgr.dispose(infoTexture_);
    }

    ImageProcessor::deinitialize();
}

void CMInfoOverlay::tryLoadIconTex() {
    if(iconTexture_) {
        TexMgr.dispose(iconTexture_);
    }
    std::string iconPath = VoreenApplication::app()->getBasePath(ICON_FILE);
    try {
        iconTexture_ = TexMgr.load(iconPath);
    } catch(...) {
        LERROR("Could not load icon texture " + iconPath);
    }
    if(iconTexture_) {
        iconTexture_->setFilter(tgt::Texture::Filter::LINEAR);
    }
}
void CMInfoOverlay::tryLoadInfoImageTex() {
    if(infoTexture_) {
        TexMgr.dispose(infoTexture_);
    }
    std::string infoImagePath = infoImageFileProp_.get();
    try {
        infoTexture_ = TexMgr.load(infoImagePath);
    } catch(...) {
        LWARNING("Could not load info texture " + infoImagePath);
    }
    if(infoTexture_) {
        infoTexture_->setFilter(tgt::Texture::Filter::LINEAR);
    }
}
tgt::vec4 CMInfoOverlay::getIconBounds() {
    tgt::vec4 bounds;

    const tgt::ivec2 screenSize = outport_.getSize();
    const float outAspectRatio = ((float)screenSize.x)/screenSize.y;
    float sizeY = sizeProp_.get();
    float sizeX = sizeY/outAspectRatio;
    float spacingX = sizeX*0.5f;
    float spacingY = sizeY*0.5f;

    switch(positionProp_.getValue()) {
        case TOP_LEFT:
            bounds.x = -1.0f + spacingX;
            bounds.y =  1.0f - spacingY - sizeY;
            break;
        case TOP_RIGHT:
            bounds.x =  1.0f - spacingX - sizeX;
            bounds.y =  1.0f - spacingY - sizeY;
            break;
        case BOTTOM_LEFT:
            bounds.x = -1.0f + spacingX;
            bounds.y = -1.0f + spacingY;
            break;
        case BOTTOM_RIGHT:
            bounds.x =  1.0f - spacingX - sizeX;
            bounds.y = -1.0f + spacingY;
            break;
    }
    bounds.z = bounds.x+sizeX;
    bounds.w = bounds.y+sizeY;
    return bounds;
}
void CMInfoOverlay::onEvent(tgt::Event *ev){
    tgt::MouseEvent *me = dynamic_cast<tgt::MouseEvent*>(ev);
    if (!me || !enableProp_.get()){
        ImageProcessor::onEvent(ev);
        return;
    }
    const tgt::ivec2 screenSize = outport_.getSize();
    const tgt::vec4 b = getIconBounds();
    tgt::ivec2 coordi = me->coord();
    tgt::vec2 coord = tgt::vec2(coordi)/tgt::vec2(screenSize) * tgt::vec2::two - tgt::vec2::one;
    coord.y *= -1;
    //Transform model coordinates further
    bool hit = coord.x >= b.x && coord.x <= b.z && coord.y >= b.y && coord.y <= b.w;
    bool activated = activatedProp_.get();
    if (me->button() == tgt::MouseEvent::MOUSE_BUTTON_LEFT && (me->action() == tgt::MouseEvent::RELEASED||me->action() == tgt::MouseEvent::EXIT) && (hit || activated)){
        activatedProp_.set(!activated);
        ev->accept();
    }
    ImageProcessor::onEvent(ev);
}

void CMInfoOverlay::process() {
    tgtAssert(outport_.isReady(), "Outport not ready");
    tgtAssert(inport_.isReady(), "Inport not ready");

    glDisable(GL_DEPTH_TEST);

    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    const tgt::ivec2 screenSize = outport_.getSize();
    const float outAspectRatio = ((float)screenSize.x)/screenSize.y;

    //First: Render inport image
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, inport_.getColorTexture()->getId());
    renderQuad();
    glBindTexture(GL_TEXTURE_2D, 0);
    LGL_ERROR;

    if(enableProp_.get()) {
        if(activatedProp_.get()) {
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            if(!infoTexture_) {
                LWARNING("Info texture is not loaded.");
            } else {
                const tgt::ivec3 texDim = infoTexture_->getDimensions();
                const float texAspectRatio = ((float)texDim.x)/texDim.y;
                float xmin, xmax, ymin, ymax;
                if(texAspectRatio>=outAspectRatio) {
                    xmin = -1.0f;
                    xmax = -xmin;
                    ymin = -outAspectRatio/texAspectRatio;
                    ymax = -ymin;
                } else {
                    xmin = -texAspectRatio/outAspectRatio;
                    xmax = -xmin;
                    ymin = -1.0f;
                    ymax = -ymin;
                }
                IMode.setTextureMode(tgt::ImmediateMode::TEX2D);
                glBindTexture(GL_TEXTURE_2D, infoTexture_->getId());
                IMode.begin(tgt::ImmediateMode::TRIANGLES);{
                    IMode.texcoord(0,0); IMode.vertex(xmin, ymin);
                    IMode.texcoord(0,1); IMode.vertex(xmin, ymax);
                    IMode.texcoord(1,1); IMode.vertex(xmax, ymax);

                    IMode.texcoord(1,1); IMode.vertex(xmax, ymax);
                    IMode.texcoord(1,0); IMode.vertex(xmax, ymin);
                    IMode.texcoord(0,0); IMode.vertex(xmin, ymin);
                }IMode.end();
                glBindTexture(GL_TEXTURE_2D, 0);
                IMode.setTextureMode(tgt::ImmediateMode::TEXNONE);
            }
        } else {
            if(!iconTexture_) {
                LWARNING("Icon texture is not loaded.");
            } else {
                const tgt::vec4 b = getIconBounds();

                glEnable(GL_BLEND);
                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

                IMode.setTextureMode(tgt::ImmediateMode::TEX2D);
                glBindTexture(GL_TEXTURE_2D, iconTexture_->getId());
                IMode.begin(tgt::ImmediateMode::TRIANGLES);{
                    IMode.texcoord(0,0); IMode.vertex(b.x, b.y);
                    IMode.texcoord(0,1); IMode.vertex(b.x, b.w);
                    IMode.texcoord(1,1); IMode.vertex(b.z, b.w);

                    IMode.texcoord(1,1); IMode.vertex(b.z, b.w);
                    IMode.texcoord(1,0); IMode.vertex(b.z, b.y);
                    IMode.texcoord(0,0); IMode.vertex(b.x, b.y);
                }IMode.end();
                glBindTexture(GL_TEXTURE_2D, 0);
                IMode.setTextureMode(tgt::ImmediateMode::TEXNONE);

                glBlendFunc(GL_ONE, GL_ZERO);
                glDisable(GL_BLEND);
            }
        }
    }
    glEnable(GL_DEPTH_TEST);
    outport_.deactivateTarget();
    LGL_ERROR;
}


} // namespace
