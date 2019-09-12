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

LabelGuard::LabelGuard(LabelProjection& labelProjection)
    : labelProjection_(labelProjection)
{
}
LabelGuard::~LabelGuard() {
    labelProjection_.ensureTexturesPresent();
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    labelProjection_.projectionTexture_->uploadTexture();
}
float& LabelGuard::at(tgt::svec2 p) {
    return labelProjection_.projection_.voxel(p.x, p.y, 0);
}

LabelProjection::LabelProjection(tgt::svec2 dimensions)
    : projection_(tgt::svec3(dimensions, 1))
    , projectionTexture_(boost::none)
{
    projection_.clear();
}
void LabelProjection::ensureTexturesPresent() {
    if(!projectionTexture_) {
        projectionTexture_ = tgt::Texture(projection_.getDimensions(), GL_RED, GL_RED, GL_FLOAT, tgt::Texture::LINEAR, tgt::Texture::CLAMP_TO_EDGE, (GLubyte*) projection_.voxel(), false);
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        projectionTexture_->uploadTexture();
    }
}
void LabelProjection::bindTexture() {
    ensureTexturesPresent();
    projectionTexture_->bind();
}

template<class Vec>
struct PolyLinePoint {
    Vec pos_;
    float d_;
};
template<class Vec>
struct PolyLine {
    PolyLine(std::deque<Vec>& points)
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

    std::deque<PolyLinePoint<Vec>> points_;

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

struct Line {
    Line(tgt::vec2 p1, tgt::vec2 p2)
        : p1_(p1)
        , p2_(p2)
    {
    }
    tgt::vec2 p1_;
    tgt::vec2 p2_;

    float len() {
        return tgt::distance(p1_, p2_);
    }
    float dist(tgt::vec2 q) {
        auto parallel = p1_-p2_;
        auto parallel_norm = tgt::normalize(parallel);

        float along = tgt::dot(parallel_norm, q - p2_);
        if(0 > along || along > len()) {
            return std::numeric_limits<float>::infinity();
        }

        tgt::vec2 orthogonal(parallel.y, -parallel.x);
        auto orth_norm = tgt::normalize(orthogonal);
        float dist = tgt::dot(orth_norm, p1_-q);
        return tgt::abs(dist);
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
        if(e->action() == tgt::MouseEvent::RELEASED && button == tgt::MouseEvent::MOUSE_BUTTON_RIGHT) {
            displayLine_.erase(displayLine_.begin() + *nearest);
        } else {
            displayLine_[*nearest] = mouse;
        }
    } else if(e->action() == tgt::MouseEvent::RELEASED && button == tgt::MouseEvent::MOUSE_BUTTON_LEFT) {
        int insert_index = -1;
        for(int i=0; i<((int)displayLine_.size())-1; ++i) {
            Line line(displayLine_[i], displayLine_[i+1]);
            float dist = line.dist(mouse);
            if(dist < MOUSE_DRAG_DIST) {
                insert_index = i;
            }
        }
        if(insert_index != -1) {
            displayLine_.insert(displayLine_.begin() + insert_index+1, mouse);
        } else {
            if(displayLine_.empty() || tgt::distance(displayLine_.front(), mouse) < tgt::distance(displayLine_.back(), mouse)) {
                displayLine_.push_front(mouse);
            } else {
                displayLine_.push_back(mouse);
            }
        }
    }

    updateProjection();

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
    , interpolationValue_("interpolationvalue", "interp", 0.0, 0.0, 1.0)
    , outputVolume_(boost::none)
    , copyShader_(nullptr)
    , projectionShader_("shader", "Shader", "interactiveprojectionlabeling.frag", "oit_passthrough.vert")
    , displayLine_()
    , projection_(boost::none)
{
    addPort(inport_);
    addPort(labelVolume_);
    addPort(overlayInput_);
    addPort(overlayOutput_);
    addPort(projectionOutput_);
    addPort(fhp_);
    addPort(lhp_);

    overlayOutput_.onSizeReceiveChange<InteractiveProjectionLabeling>(this, &InteractiveProjectionLabeling::updateSizes);

    addProperty(projectionShader_);
    addProperty(interpolationValue_);
}

void InteractiveProjectionLabeling::updateSizes() {
    overlayInput_.requestSize(overlayOutput_.getReceivedSize());
    updateProjection();
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
    auto program = projectionShader_.getShader();
    if(!program || !program->isLinked()) {
        LGL_ERROR;
        LERROR("Shader not compiled!");
        return;
    }

    projectionOutput_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if(projection_) {
        auto& p = *projection_;

        tgt::TextureUnit unit;
        unit.activate();
        p.bindTexture();

        program->activate();
        //program->setUniform("dimensions_", tgt::ivec3(vol.getDimensions()));
        //program->setUniform("realToProjectedMat_", p.realToProjected());
        //program->setUniform("projectedToRealMat_", p.projectedToReal());
        //program->setUniform("projectionRange_", projectionRange);
        //program->setUniform("volumeTex_", volUnit.getUnitNumber());
        program->setUniform("tex_", unit.getUnitNumber());

        glDepthFunc(GL_ALWAYS);
        renderQuad();
        glDepthFunc(GL_LESS);

        program->deactivate();
        glActiveTexture(GL_TEXTURE0);
    }
    projectionOutput_.deactivateTarget();
    LGL_ERROR;
}

void InteractiveProjectionLabeling::withOutputVolume(std::function<void(LZ4SliceVolume<uint8_t>&)> func) {
    if(outputVolume_) {
        labelVolume_.setData(nullptr);
        func(*outputVolume_);
        labelVolume_.setData(LZ4SliceVolume<uint8_t>::open(outputVolume_->getFilePath()).toVolume().release());
    }
}

void InteractiveProjectionLabeling::updateProjection() {
    if(!tgt::isInitedGL() || !fhp_.getColorTexture() || !lhp_.getColorTexture()) {
        return;
    }

    if(displayLine_.empty()) {
        return;
    }

    if(!inport_.hasData()) {
        return;
    }
    const auto& vol = *inport_.getData();

    const auto volram_ptr = vol.getRepresentation<VolumeRAM>();
    const auto& volram = *volram_ptr;

    auto dim = overlayOutput_.getReceivedSize();
    projection_ = LabelProjection(dim);
    auto proj = projection_->projection_mut();

    VolumeAtomic<tgt::vec4> front((tgt::vec4*)fhp_.getColorTexture()->downloadTextureToBuffer(GL_RGBA, GL_FLOAT), tgt::svec3(fhp_.getSize(),1));
    VolumeAtomic<tgt::vec4> back((tgt::vec4*)lhp_.getColorTexture()->downloadTextureToBuffer(GL_RGBA, GL_FLOAT), tgt::svec3(lhp_.getSize(),1));

    auto line = PolyLine<tgt::vec2>(displayLine_);

    auto tex_to_vox = vol.getTextureToVoxelMatrix();
    for(int x = 0; x < dim.x; ++x) {
        float d = ((float)x)/(dim.x-1);
        auto p = line.interpolate(d);

        tgt::vec3 normalized_query(p, 0);
        tgt::vec4 front_pos = tex_to_vox * front.getVoxelLinear(normalized_query * tgt::vec3(front.getDimensions()));
        tgt::vec4 back_pos = tex_to_vox * back.getVoxelLinear(normalized_query * tgt::vec3(back.getDimensions()));

        for(int y = 0; y < dim.y; ++y) {
            float alpha = ((float)y)/(dim.y-1);

            auto query_pos = front_pos * (1-alpha) + alpha * back_pos;

            float val = volram.getVoxelNormalizedLinear(query_pos.xyz());

            proj.at(tgt::svec2(x, y)) = val;
        }
    }
}

InteractiveProjectionLabeling::~InteractiveProjectionLabeling() {
}

void InteractiveProjectionLabeling::process() {
    if (getInvalidationLevel() == INVALID_PROGRAM) {
        //std::string header;
        //header += "#define UNLABELED " + std::to_string(LabelProjection::UNLABELED) + "\n";
        //header += "#define BACKGROUND " + std::to_string(LabelProjection::BACKGROUND) + "\n";
        //header += "#define FOREGROUND " + std::to_string(LabelProjection::FOREGROUND) + "\n";
        //header += "#define SUGGESTED_FOREGROUND " + std::to_string(LabelProjection::SUGGESTED_FOREGROUND) + "\n";
        //header += "#define INCONSISTENT " + std::to_string(LabelProjection::INCONSISTENT) + "\n";
        //projectionShader_.setHeader(header);
        projectionShader_.rebuild();
    }

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

    updateProjection();
}
VoreenSerializableObject* InteractiveProjectionLabeling::create() const {
    return new InteractiveProjectionLabeling();
}

} // namespace voreen
