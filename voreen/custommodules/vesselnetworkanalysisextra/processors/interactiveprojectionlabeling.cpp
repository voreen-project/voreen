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
tgt::vec2& LabelGuard::at(tgt::svec2 p) {
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
        projectionTexture_ = tgt::Texture(projection_.getDimensions(), GL_RG, GL_RG, GL_FLOAT, tgt::Texture::LINEAR, tgt::Texture::CLAMP_TO_EDGE, (GLubyte*) projection_.voxel(), false);
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

    tgt::vec2 interpolate(float d) const {
        tgtAssert(0.0 <= d && d <= 1.0, "Invalid interpolation parameter");
        if(points_.size() == 1) {
            return points_[0].pos_;
        }
        int i=0;
        while(d > points_[i+1].d_ && i < ((int)points_.size())-2) {
            ++i;
        }
        auto& p1 = points_[i];
        auto& p2 = points_[i+1];
        if(p1.d_ == p2.d_) {
            return p1.pos_;
        }
        float alpha = (d - p1.d_) / (p2.d_- p1.d_);
        tgt::vec2 res = p1.pos_ * (1-alpha) + alpha * p2.pos_;
        tgtAssert(std::isfinite(res.x) && std::isfinite(res.y), "Invalid interpolation result");
        return res;
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

static tgt::vec2 projectionDepthRange(const VolumeBase& vol, const VolumeAtomic<tgt::vec4>& front, const VolumeAtomic<tgt::vec4>& back, PolyLine<tgt::vec2>& line, tgt::vec3 camera) {

    auto dim = vol.getDimensions();

    auto tex_to_world = vol.getTextureToWorldMatrix();

    float min_dist = std::numeric_limits<float>::infinity();
    float max_dist = 0.0f;
    for(int x = 0; x < dim.x; ++x) {
        float d = ((float)x)/(dim.x-1);
        auto p = line.interpolate(d);

        tgt::vec3 normalized_query(p, 0);
        tgt::vec4 front_pos = front.getVoxelLinear(normalized_query * tgt::vec3(front.getDimensions()));
        tgt::vec4 back_pos = back.getVoxelLinear(normalized_query * tgt::vec3(back.getDimensions()));

        if(front_pos.a > 0) {
            min_dist = std::min(min_dist, tgt::distance(camera, (tex_to_world * front_pos).xyz()));
        }
        if(back_pos.a > 0) {
            max_dist = std::max(max_dist, tgt::distance(camera, (tex_to_world * back_pos).xyz()));
        }
    }

    return tgt::vec2(min_dist, max_dist);
}


#define MOUSE_INTERACTION_DIST 0.02f
static void handleLineEvent(std::deque<tgt::vec2>& points, tgt::MouseEvent* e) {
    auto button = e->button();
    if((button & (tgt::MouseEvent::MOUSE_BUTTON_LEFT | tgt::MouseEvent::MOUSE_BUTTON_RIGHT)) == 0) {
        return;
    }

    tgt::ivec2 coords = e->coord();
    tgt::ivec2 viewport = e->viewport();

    coords.y = viewport.y - coords.y - 1;
    auto mouse = tgt::vec2(coords)/tgt::vec2(viewport);

    boost::optional<int> nearest = boost::none;
    int i = 0;
    for(auto& p : points) {
        float dist = tgt::distance(p, mouse);
        if (dist < MOUSE_INTERACTION_DIST && (!nearest || dist < tgt::distance(points[i], mouse))) {
            nearest = i;
        }
        ++i;
    }
    if(nearest) {
        if(e->action() == tgt::MouseEvent::RELEASED && button == tgt::MouseEvent::MOUSE_BUTTON_RIGHT) {
            points.erase(points.begin() + *nearest);
        } else {
            points[*nearest] = mouse;
        }
    } else if(e->action() == tgt::MouseEvent::RELEASED && button == tgt::MouseEvent::MOUSE_BUTTON_LEFT) {
        int insert_index = -1;
        for(int i=0; i<((int)points.size())-1; ++i) {
            Line line(points[i], points[i+1]);
            float dist = line.dist(mouse);
            if(dist < MOUSE_INTERACTION_DIST) {
                insert_index = i;
            }
        }
        if(insert_index != -1) {
            points.insert(points.begin() + insert_index+1, mouse);
        } else {
            if(points.empty() || tgt::distance(points.front(), mouse) < tgt::distance(points.back(), mouse)) {
                points.push_front(mouse);
            } else {
                points.push_back(mouse);
            }
        }
    }
    e->accept();
}
void InteractiveProjectionLabeling::projectionEvent(tgt::MouseEvent* e) {
    auto button = e->button();
    if((button & (tgt::MouseEvent::MOUSE_BUTTON_LEFT | tgt::MouseEvent::MOUSE_BUTTON_RIGHT)) == 0) {
        return;
    }

    if(e->modifiers() == tgt::Event::CTRL) {
        handleLineEvent(projectionLabels_.lowerBackground_, e);
        std::sort(projectionLabels_.lowerBackground_.begin(), projectionLabels_.lowerBackground_.end(), [] (const tgt::vec2& p1, const tgt::vec2& p2) {
                return p1.x < p2.x;
                });
    } else if(e->modifiers() == tgt::Event::SHIFT) {
        handleLineEvent(projectionLabels_.upperBackground_, e);
        std::sort(projectionLabels_.upperBackground_.begin(), projectionLabels_.upperBackground_.end(), [] (const tgt::vec2& p1, const tgt::vec2& p2) {
                return p1.x < p2.x;
                });
    } else if(e->modifiers() == tgt::Event::MODIFIER_NONE) {
        handleLineEvent(projectionLabels_.foreground_, e);
        std::sort(projectionLabels_.foreground_.begin(), projectionLabels_.foreground_.end(), [] (const tgt::vec2& p1, const tgt::vec2& p2) {
                return p1.x < p2.x;
                });
    } else {
        return;
    }
    projectionLabelsModified_ = true;

    invalidate();
}

void InteractiveProjectionLabeling::overlayEvent(tgt::MouseEvent* e) {
    auto button = e->button();

    if(e->modifiers() != tgt::Event::CTRL
            || ((button & (tgt::MouseEvent::MOUSE_BUTTON_LEFT | tgt::MouseEvent::MOUSE_BUTTON_RIGHT)) == 0))
             {
        return;
    }

    handleLineEvent(displayLine_, e);

    updateProjection();

    invalidate();
}

void InteractiveProjectionLabeling::onPortEvent(tgt::Event* e, Port* port) {
    if(tgt::MouseEvent* me = dynamic_cast<tgt::MouseEvent*>(e)) {
        if(port == &overlayOutput_) {
            overlayEvent(me);
            if(!displayLine_.empty() && state_ == FREE) {
                state_ = LABELING;
                projectionLabelsModified_ = false;
            }
            if(state_ == LABELING) {
                me->accept();
            }
        }
        else if(port == &projectionOutput_) {
            projectionEvent(me);

            // Definitely consume events for this port.
            me->accept();
        }
    } else if(tgt::KeyEvent* ke = dynamic_cast<tgt::KeyEvent*>(e)) {
        switch(ke->keyCode()) {
            case tgt::KeyEvent::K_ESCAPE: {
                displayLine_.clear();
                projectionLabels_.foreground_.clear();
                projectionLabels_.lowerBackground_.clear();
                projectionLabels_.upperBackground_.clear();
                projection_ = boost::none;
                state_ = FREE;
                ke->accept();
                invalidate();
                break;
            }
            case tgt::KeyEvent::K_SPACE: {
                if(state_ == LABELING && !projectionLabels_.foreground_.empty()) {
                    finishProjection();
                    displayLine_.clear();
                    projectionLabels_.foreground_.clear();
                    projectionLabels_.lowerBackground_.clear();
                    projectionLabels_.upperBackground_.clear();
                    state_ = FREE;
                    invalidate();
                }
                ke->accept();
                break;
            }
            default:;
        }
    }
    if(!e->isAccepted() && state_ == FREE) {
        RenderProcessor::onPortEvent(e, port);
    }
}

void InteractiveProjectionLabeling::finishProjection() {

    if(!inport_.hasData()) {
        return;
    }
    const auto& vol = *inport_.getData();

    const int NUM_SAMPLES = 100;

    PolyLine<tgt::vec2> displayLine(displayLine_);

    auto maybe_front = getFhp();
    auto maybe_back = getLhp();
    if(!maybe_front || !maybe_back) {
        return;
    }
    auto& front = *maybe_front;
    auto& back = *maybe_back;

    tgt::vec3 camera = camera_.get().getPosition();

    auto tex_to_world = vol.getTextureToWorldMatrix();

    auto minmax = projectionDepthRange(vol, front, back, displayLine, camera);
    float min_dist = minmax.x;
    float max_dist = minmax.y;

    auto project3D = [&] (const PolyLine<tgt::vec2>& projectionLine) {
        std::vector<tgt::vec3> segment;
        for(int i=0; i<NUM_SAMPLES; ++i) {
            float projection_d = static_cast<float>(i)/(NUM_SAMPLES-1);

            tgt::vec2 projectionPoint = projectionLine.interpolate(projection_d);
            float normalized_depth = projectionPoint.y;

            float depth = normalized_depth * (max_dist - min_dist) + min_dist;

            float display_d = projectionPoint.x;

            tgt::vec2 display_point = displayLine.interpolate(display_d);

            tgt::vec3 normalized_query(display_point, 0.0);
            tgt::vec4 front_pos = front.getVoxelLinear(normalized_query * tgt::vec3(front.getDimensions()));
            tgt::vec4 back_pos = back.getVoxelLinear(normalized_query * tgt::vec3(front.getDimensions()));

            tgt::vec4 front_world = tex_to_world * front_pos;
            tgt::vec4 back_world = tex_to_world * back_pos;

            tgt::vec3 view_dir = tgt::normalize(back_world.xyz() - front_world.xyz());

            tgt::vec3 point = camera + depth*view_dir;
            segment.push_back(point);
        }
        return segment;
    };

    if(!projectionLabels_.foreground_.empty()) {
        foregroundLabelLines_.addSegment(project3D(projectionLabels_.foreground_));
    }

    if(!projectionLabels_.lowerBackground_.empty()) {
        backgroundLabelLines_.addSegment(project3D(projectionLabels_.lowerBackground_));
    }
    if(!projectionLabels_.upperBackground_.empty()) {
        backgroundLabelLines_.addSegment(project3D(projectionLabels_.upperBackground_));
    }
}
InteractiveProjectionLabeling::InteractiveProjectionLabeling()
    : RenderProcessor()
    , inport_(Port::INPORT, "interactiveprojectionlabeling.inport", "Volume Input")
    //, labelVolume_(Port::OUTPORT, "interactiveprojectionlabeling.labelVolume", "Labels Output")
    , foregroundLabelGeometry_(Port::OUTPORT, "interactiveprojectionlabeling.foregroundLabelGeometry", "Foreground Labels Output")
    , backgroundLabelGeometry_(Port::OUTPORT, "interactiveprojectionlabeling.backgroundLabelGeometry", "Background Labels Output")
    , overlayOutput_(Port::OUTPORT, "interactiveprojectionlabeling.overlayoutput", "Overlay (3D)", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , projectionOutput_(Port::OUTPORT, "interactiveprojectionlabeling.projectionoutput", "Projection (2D)", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , fhp_(Port::INPORT, "interactiveprojectionlabeling.fhp", "First hit points", false)
    , lhp_(Port::INPORT, "interactiveprojectionlabeling.lhp", "Last hit points", false)
    , camera_("camera", "Camera")
    , initializationMode_("initializationMode", "Initialization Mode")
    , outputVolume_(boost::none)
    , projectionShader_("shader", "Shader", "interactiveprojectionlabeling.frag", "oit_passthrough.vert")
    , displayLine_()
    , projection_(boost::none)
    , projectionLabels_()
    , projectionLabelsModified_(false)
    , foregroundLabelLines_()
    , backgroundLabelLines_()
    , state_(FREE)
{
    addPort(inport_);
    //addPort(labelVolume_);
    addPort(foregroundLabelGeometry_);
    addPort(backgroundLabelGeometry_);
    addPort(overlayOutput_);
    addPort(projectionOutput_);
    addPort(fhp_);
    addPort(lhp_);

    overlayOutput_.onSizeReceiveChange<InteractiveProjectionLabeling>(this, &InteractiveProjectionLabeling::updateSizes);

    addProperty(projectionShader_);
    addProperty(camera_);
    addProperty(initializationMode_);
        initializationMode_.addOption("none", "None", NONE);
        initializationMode_.addOption("brightlumen", "Bright Lumen", BRIGHT_LUMEN);
        ON_CHANGE(initializationMode_, InteractiveProjectionLabeling, initializeProjectionLabels);
    //initializationMode_.addOption("brightwall", "Bright Wall", BRIGHT_WALL);
}

void InteractiveProjectionLabeling::updateSizes() {
    updateProjection();
}


void InteractiveProjectionLabeling::initialize() {
    RenderProcessor::initialize();
}

void InteractiveProjectionLabeling::deinitialize() {
    RenderProcessor::deinitialize();
}
static void renderLine(const std::deque<tgt::vec2>& points, tgt::vec3 color) {
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.ortho(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

    IMode.color(color);
    IMode.begin(tgt::ImmediateMode::LINE_STRIP);
    for(auto& p : points) {
        IMode.vertex(p);
    }
    IMode.end();

    glPointSize(5.0);
    IMode.begin(tgt::ImmediateMode::POINTS);
    for(auto& p : points) {
        IMode.vertex(p);
    }
    IMode.end();

    MatStack.popMatrix();

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    IMode.color(tgt::vec4(1.0));
    glPointSize(1.0);
}
void InteractiveProjectionLabeling::renderOverlay() {
    overlayOutput_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    renderLine(displayLine_, tgt::vec3(1.0, 0.0, 0.0));

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

    renderLine(projectionLabels_.lowerBackground_, tgt::vec3(0.0, 1.0, 0.0));
    renderLine(projectionLabels_.foreground_, tgt::vec3(1.0, 0.0, 0.0));
    renderLine(projectionLabels_.upperBackground_, tgt::vec3(0.0, 1.0, 0.0));

    projectionOutput_.deactivateTarget();
    LGL_ERROR;
}

void InteractiveProjectionLabeling::withOutputVolume(std::function<void(LZ4SliceVolume<uint8_t>&)> func) {
    if(outputVolume_) {
        //labelVolume_.setData(nullptr);
        func(*outputVolume_);
        //labelVolume_.setData(LZ4SliceVolume<uint8_t>::open(outputVolume_->getFilePath()).toVolume().release());
    }
}

boost::optional<VolumeAtomic<tgt::vec4>> InteractiveProjectionLabeling::getFhp() {
    if(!tgt::isInitedGL() || !fhp_.getColorTexture()) {
        return boost::none;
    }
    return VolumeAtomic<tgt::vec4>((tgt::vec4*)fhp_.getColorTexture()->downloadTextureToBuffer(GL_RGBA, GL_FLOAT), tgt::svec3(fhp_.getSize(),1));
}
boost::optional<VolumeAtomic<tgt::vec4>> InteractiveProjectionLabeling::getLhp() {
    if(!tgt::isInitedGL() || !lhp_.getColorTexture()) {
        return boost::none;
    }
    return VolumeAtomic<tgt::vec4>((tgt::vec4*)lhp_.getColorTexture()->downloadTextureToBuffer(GL_RGBA, GL_FLOAT), tgt::svec3(lhp_.getSize(),1));
}

static std::vector<int> maxPath(const VolumeAtomic<float>& img) {
    VolumeAtomic<int> paths(img.getDimensions());
    tgt::ivec3 idim = img.getDimensions();
    std::vector<float> global_cost(idim.x, 0.0);
    const int MAX_NEIGHBOR_OFFSET = 3;
    for(int y=0; y < idim.y; ++y) {
        std::vector<float> next_global_cost(idim.x, 0.0);
        for(int x=0; x < idim.x; ++x) {
            int best_i = 0;
            float best_val = 0.0;
            for(int d=std::max(0, x-MAX_NEIGHBOR_OFFSET); d<std::min(idim.x, x+MAX_NEIGHBOR_OFFSET+1); ++d) {
                float val = img.voxel(d, y, 0) + global_cost.at(d);
                if(val > best_val) {
                    best_val = val;
                    best_i = d;
                }
            }
            paths.voxel(x, y, 0) = best_i;
            next_global_cost[x] = best_val;
        }
        global_cost = next_global_cost;
    }

    std::vector<int> path;
    float best_begin_val = 0.0;
    int best_begin_i = 0;
    for(int x=0; x < idim.x; ++x) {
        if(global_cost[x] > best_begin_val) {
            best_begin_val = global_cost[x];
            best_begin_i = x;
        }
    }
    int x = best_begin_i;
    path.push_back(x);
    for(int y=idim.y-2; y >= 0; --y) {
        x = paths.voxel(x, y, 0);
        path.push_back(x);
    }
    std::reverse(path.begin(), path.end());
    return path;
}
static void initBrightLumen(const LabelProjection& proj, ProjectionLabels& labels) {
    const auto& orig = proj.projection();
    tgt::ivec3 idim = orig.getDimensions();
    VolumeAtomic<float> top_gradients(tgt::svec3(idim.y, idim.x, idim.z));
    VolumeAtomic<float> bottom_gradients(tgt::svec3(idim.y, idim.x, idim.z));

    for(int y=0; y < idim.y; ++y) {
        for(int x=0; x < idim.x; ++x) {
            tgt::vec2 left = orig.voxel(x, std::max(0, y-1), 0);
            tgt::vec2 right = orig.voxel(x, std::min(y+1, idim.y-1), 0);
            float diff;
            if(left.y > 0.0 && right.y > 0.0) {
                diff = left.x - right.x;
            } else {
                diff = 0.0;
            }
            top_gradients.voxel(y, x, 0) = std::max(0.0f, diff);
            bottom_gradients.voxel(y, x, 0) = std::max(0.0f, -diff);
        }
    }

    auto bottom_path = maxPath(bottom_gradients);
    auto top_path = maxPath(top_gradients);

    tgtAssert(bottom_path.size() == top_path.size(), "Path size mismatch");

    for(int x=0; x < top_path.size(); ++x) {
        float x_pos = static_cast<float>(x)/(idim.x-1);
        float y_top = static_cast<float>(top_path.at(x))/(idim.y-1);
        float y_bottom = static_cast<float>(bottom_path.at(x))/(idim.y-1);

        float width = y_top - y_bottom;
        float center = (y_top + y_bottom)/2;
        labels.foreground_.emplace_back(x_pos, center);
        labels.lowerBackground_.emplace_back(x_pos, center-width);
        labels.upperBackground_.emplace_back(x_pos, center+width);
    }
}
void InteractiveProjectionLabeling::initializeProjectionLabels() {
    projectionLabels_.foreground_.clear();
    projectionLabels_.lowerBackground_.clear();
    projectionLabels_.upperBackground_.clear();
    switch(initializationMode_.getValue()) {
        case BRIGHT_LUMEN:
            {
                initBrightLumen(*projection_, projectionLabels_);
                break;
            }
        default:;
    }
    projectionLabelsModified_ = false;
}


void InteractiveProjectionLabeling::updateProjection() {
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

    auto maybe_front = getFhp();
    auto maybe_back = getLhp();
    if(!maybe_front || !maybe_back) {
        return;
    }
    auto& front = *maybe_front;
    auto& back = *maybe_back;

    tgt::vec3 camera = camera_.get().getPosition();

    auto line = PolyLine<tgt::vec2>(displayLine_);
    auto tex_to_world = vol.getTextureToWorldMatrix();

    auto minmax = projectionDepthRange(vol, front, back, line, camera);
    float min_dist = minmax.x;
    float max_dist = minmax.y;

    auto world_to_vox = vol.getWorldToVoxelMatrix();

    tgt::vec3 dimf = vol.getDimensions();
    for(int x = 0; x < dim.x; ++x) {
        float d = ((float)x)/(dim.x-1);
        auto p = line.interpolate(d);

        tgt::vec3 normalized_query(p, 0);
        tgt::vec4 front_pos = front.getVoxelLinear(normalized_query * tgt::vec3(front.getDimensions()));
        tgt::vec4 back_pos = back.getVoxelLinear(normalized_query * tgt::vec3(back.getDimensions()));

        tgt::vec4 front_world = tex_to_world * front_pos;
        tgt::vec4 back_world = tex_to_world * back_pos;

        tgt::vec3 view_dir = tgt::normalize(back_world.xyz() - front_world.xyz());

        for(int y = 0; y < dim.y; ++y) {
            float alpha = ((float)y)/(dim.y-1);
            float alpha_rw = max_dist * alpha + (1.0 - alpha) * min_dist;

            tgt::vec4 query_pos_rw(view_dir * alpha_rw + camera, 1.0);
            tgt::vec3 query_pos = (world_to_vox * query_pos_rw).xyz();

            tgt::vec2 val;
            if(tgt::hor(tgt::greaterThan(query_pos, dimf)) || tgt::hor(tgt::lessThan(query_pos, tgt::vec3::zero))) {
                val = tgt::vec2(0.0, 0.0);
            } else {
                val = tgt::vec2(volram.getVoxelNormalizedLinear(query_pos), 1.0);
            }

            proj.at(tgt::svec2(x, y)) = val;
        }
    }

    if(!projectionLabelsModified_) {
        initializeProjectionLabels();
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

    foregroundLabelGeometry_.setData(&foregroundLabelLines_, false);
    backgroundLabelGeometry_.setData(&backgroundLabelLines_, false);
    // initialize shader

    //auto* vol = new VolumeAtomic<uint8_t>(xz_.labels().getDimensions());
    //for(int i=0; i<tgt::hmul(xz_.labels().getDimensions()); ++i) {
    //    vol->voxel(i) = xz_.labels().voxel(i);
    //}
    //labelVolume_.setData(new Volume(vol, inport_.getData()));
}


void InteractiveProjectionLabeling::adjustPropertiesToInput() {
    //labelVolume_.setData(nullptr);
    if(!inport_.hasData()) {
        return;
    }
    const auto& vol = *inport_.getData();
    auto dim = vol.getDimensions();

    //LZ4SliceVolumeBuilder<uint8_t> builder(
    //        VoreenApplication::app()->getUniqueTmpFilePath("." + LZ4SliceVolumeBase::FILE_EXTENSION),
    //        LZ4SliceVolumeMetadata::fromVolume(vol).withRealWorldMapping(RealWorldMapping::createDenormalizingMapping<uint8_t>())
    //        );
    //builder.fill(0);
    //outputVolume_ = std::move(builder).finalize();
    //withOutputVolume([&] (LZ4SliceVolume<uint8_t>& vol) { });

    updateProjection();
}
VoreenSerializableObject* InteractiveProjectionLabeling::create() const {
    return new InteractiveProjectionLabeling();
}

} // namespace voreen
