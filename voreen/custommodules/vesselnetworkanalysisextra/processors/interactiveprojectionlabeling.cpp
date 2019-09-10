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

#include "interactiveprojectionlabeling.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include <iostream>
#include "tgt/textureunit.h"
#include "tgt/init.h"
#include "tgt/immediatemode/immediatemode.h"

namespace voreen {

const std::string InteractiveProjectionLabeling::loggerCat_("voreen.vesselnetworkanalysisextra.interactiveprojectionlabeling");

template<class Vec>
struct PolyLinePoint {
    Vec pos_;
    float d_;
};
template<class Vec>
struct PolyLine {
    PolyLine(std::vector<Vec>& points)
        : points_()
    {
        tgtAssert(!points.empty(), "Points must not be empty!")
        float len = 0.0;
        for(int i=0; i<((int)points.size())-1; ++i) {
            len += tgt::distance(points[i], points[i+1]);
        }
        float norm_len = 0.0;
        for(int i=0; i<((int)points.size())-1; ++i) {
            points_.push_back( PolyLinePoint<Vec> {
                    points[i],
                    norm_len
                    });
            norm_len += tgt::distance(points[i], points[i+1])/len;
        }
        points_.push_back( PolyLinePoint<Vec> {
                points.back(),
                norm_len
                });
    }

    std::vector<PolyLinePoint<Vec>> points_;

    tgt::vec2 interpolate(float d) {
        tgtAssert(0.0 <= d && d <= 1.0, "Invalid interpolation parameter");
        if(points_.size() == 1) {
            return points_[0].pos_;
        }
        int i;
        for(i=0; i<((int)points_.size())-1; ++i) {
            if(d <= points_[i+1].d_) {
                break;
            }
        }
        auto& p1 = points_[i];
        auto& p2 = points_[i+1];
        float alpha = (d - p1.d_) / (p2.d_- p1.d_);
        return p1.pos_ * (1-alpha) + alpha * p2.pos_;
    }
};

void InteractiveProjectionLabeling::projectionEvent(tgt::MouseEvent* e) {
}
void InteractiveProjectionLabeling::overlayEvent(tgt::MouseEvent* e) {
    auto button = e->button();
    if(e->modifiers() != tgt::Event::CTRL
            || ((button & (tgt::MouseEvent::MOUSE_BUTTON_LEFT | tgt::MouseEvent::MOUSE_BUTTON_RIGHT)) == 0))
             {
        return;
    }
    tgt::ivec2 coords = e->coord();
    tgt::ivec2 viewport = e->viewport();

    coords.y = viewport.y - coords.y - 1;
    auto mouse = tgt::vec2(coords)/tgt::vec2(viewport);

    boost::optional<int> nearest = boost::none;
    int i = 0;
    const float MOUSE_DRAG_DIST = 0.02;
    for(auto& p : displayLine_) {
        float dist = tgt::distance(p, mouse);
        if (dist < MOUSE_DRAG_DIST && (!nearest || dist < tgt::distance(displayLine_[i], mouse))) {
            nearest = i;
        }
        ++i;
    }
    if(nearest) {
        if(e->action() == tgt::MouseEvent::PRESSED && button == tgt::MouseEvent::MOUSE_BUTTON_RIGHT) {
            displayLine_.erase(displayLine_.begin() + *nearest);
        } else {
            displayLine_[*nearest] = mouse;
        }
    } else {
        if(e->action() == tgt::MouseEvent::PRESSED && button == tgt::MouseEvent::MOUSE_BUTTON_LEFT) {
            displayLine_.push_back(mouse);
        }
    }

    invalidate();
    e->accept();
}

void InteractiveProjectionLabeling::onPortEvent(tgt::Event* e, Port* port) {
    if(tgt::MouseEvent* me = dynamic_cast<tgt::MouseEvent*>(e)) {
        if(port == &overlayOutput_) {
            overlayEvent(me);
        }
        else if(port == &projectionOutput_) {
            projectionEvent(me);

            // Definitely consume events for this port.
            me->accept();
        }
    }
    if(!e->isAccepted()) {
        RenderProcessor::onPortEvent(e, port);
    }
}

InteractiveProjectionLabeling::InteractiveProjectionLabeling()
    : RenderProcessor()
    , inport_(Port::INPORT, "interactiveprojectionlabeling.inport", "Volume Input")
    , labelVolume_(Port::OUTPORT, "interactiveprojectionlabeling.labelVolume", "Labels Output")
    , overlayInput_(Port::INPORT, "interactiveprojectionlabeling.overlayinput", "Overlay Input", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN)
    , overlayOutput_(Port::OUTPORT, "interactiveprojectionlabeling.overlayoutput", "Overlay (3D)", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , projectionOutput_(Port::OUTPORT, "interactiveprojectionlabeling.projectionoutput", "Projection (2D)", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , fhp_(Port::INPORT, "interactiveprojectionlabeling.fhp", "First hit points", false)
    , lhp_(Port::INPORT, "interactiveprojectionlabeling.lhp", "Last hit points", false)
    , outputVolume_(boost::none)
    , copyShader_(nullptr)
    , displayLine_()
    , interpolationValue_("interpolationvalue", "interp", 0.0, 0.0, 1.0)
{
    addPort(inport_);
    addPort(labelVolume_);
    addPort(overlayInput_);
    addPort(overlayOutput_);
    addPort(projectionOutput_);
    addPort(fhp_);
    addPort(lhp_);

    overlayOutput_.onSizeReceiveChange<InteractiveProjectionLabeling>(this, &InteractiveProjectionLabeling::updateSizes);

    addProperty(interpolationValue_);
}

void InteractiveProjectionLabeling::updateSizes() {
    overlayInput_.requestSize(overlayOutput_.getReceivedSize());
}


void InteractiveProjectionLabeling::initialize() {
    RenderProcessor::initialize();

    copyShader_ = ShdrMgr.loadSeparate("passthrough.vert", "copyimage.frag", RenderProcessor::generateHeader(), false);
}

void InteractiveProjectionLabeling::deinitialize() {
    ShdrMgr.dispose(copyShader_);
    copyShader_ = nullptr;

    RenderProcessor::deinitialize();
}
void InteractiveProjectionLabeling::renderOverlay() {
    overlayOutput_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    tgtAssert(copyShader_, "Shader missing");

    tgt::TextureUnit imageUnit, imageUnitDepth;
    overlayInput_.bindTextures(imageUnit.getEnum(), imageUnitDepth.getEnum());

    copyShader_->activate();
    setGlobalShaderParameters(copyShader_);
    overlayInput_.setTextureParameters(copyShader_, "texParams_");
    copyShader_->setUniform("colorTex_", imageUnit.getUnitNumber());
    copyShader_->setUniform("depthTex_", imageUnitDepth.getUnitNumber());
    renderQuad();
    copyShader_->deactivate();
    LGL_ERROR;

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.ortho(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

    IMode.color(tgt::vec3(1.0, 1.0, 0.0));
    IMode.begin(tgt::ImmediateMode::LINE_STRIP);
    for(auto& p : displayLine_) {
        IMode.vertex(p);
    }
    IMode.end();

    glPointSize(5.0);
    IMode.begin(tgt::ImmediateMode::POINTS);
    for(auto& p : displayLine_) {
        IMode.vertex(p);
    }
    IMode.end();

    IMode.color(tgt::vec3(1.0, 0.0, 0.0));
    if(!displayLine_.empty()) {
        IMode.begin(tgt::ImmediateMode::POINTS);
        IMode.vertex(PolyLine<tgt::vec2>(displayLine_).interpolate(interpolationValue_.get()));
        IMode.end();
    }

    MatStack.popMatrix();

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    IMode.color(tgt::vec4(1.0));
    glPointSize(1.0);

    overlayOutput_.deactivateTarget();
}
void InteractiveProjectionLabeling::renderProjection() {
    projectionOutput_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    projectionOutput_.deactivateTarget();
}

void InteractiveProjectionLabeling::withOutputVolume(std::function<void(LZ4SliceVolume<uint8_t>&)> func) {
    if(outputVolume_) {
        labelVolume_.setData(nullptr);
        func(*outputVolume_);
        labelVolume_.setData(LZ4SliceVolume<uint8_t>::open(outputVolume_->getFilePath()).toVolume().release());
    }
}

InteractiveProjectionLabeling::~InteractiveProjectionLabeling() {
}

void InteractiveProjectionLabeling::process() {
    renderOverlay();
    renderProjection();
    // initialize shader

    //auto* vol = new VolumeAtomic<uint8_t>(xz_.labels().getDimensions());
    //for(int i=0; i<tgt::hmul(xz_.labels().getDimensions()); ++i) {
    //    vol->voxel(i) = xz_.labels().voxel(i);
    //}
    //labelVolume_.setData(new Volume(vol, inport_.getData()));
}


void InteractiveProjectionLabeling::adjustPropertiesToInput() {
    labelVolume_.setData(nullptr);
    if(!inport_.hasData()) {
        return;
    }
    const auto& vol = *inport_.getData();
    auto dim = vol.getDimensions();

    LZ4SliceVolumeBuilder<uint8_t> builder(
            VoreenApplication::app()->getUniqueTmpFilePath("." + LZ4SliceVolumeBase::FILE_EXTENSION),
            LZ4SliceVolumeMetadata::fromVolume(vol).withRealWorldMapping(RealWorldMapping::createDenormalizingMapping<uint8_t>())
            );
    builder.fill(0);
    outputVolume_ = std::move(builder).finalize();
    withOutputVolume([&] (LZ4SliceVolume<uint8_t>& vol) { });
}
VoreenSerializableObject* InteractiveProjectionLabeling::create() const {
    return new InteractiveProjectionLabeling();
}

} // namespace voreen
