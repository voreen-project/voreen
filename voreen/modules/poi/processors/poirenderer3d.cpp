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

#include "poirenderer3d.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "tgt/immediatemode/immediatemode.h"
#include "tgt/font.h"
#include "voreen/core/voreenapplication.h"
namespace voreen {

namespace{
    struct VertexLayout3D{
        tgt::vec3 pos;
        tgt::vec3 color;
        float selected;
    };
}   

POIRenderer3d::POIRenderer3d()
    : RenderProcessor()
    , renderPort_(Port::OUTPORT, "renderport", "Renderport", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , volumeInport_(Port::INPORT, "volumeInport", "Volume port")
    , shader_("shader", "Shader", "poipoint3d.frag", "poipoint3d.vert", "poipoint3d.geom", Processor::INVALID_PROGRAM, Property::LOD_DEBUG)
    , cpPort_(Port::INPORT, "cpPort", "Coprocessors", false)
    , camera_("camera", "Camera")
    , radius_("radius", "Radius", 1.0, 0.0, 100.0)
    , autoadjustRadius_("autoadjustradius", "Autoadjust radius", false)
    , radiusDivieder_("radiusDivider", "Radius divider", 100.0, 1.0, 1000.0)
    , distanceMeasureColor_("distanceMeasureColor", "Color for distance measurement")
    , selectionMode_("selectionMode", "Selection Mode")
{
    addPort(cpPort_);
    addPort(renderPort_);
    addPort(volumeInport_);
    
    addProperty(camera_);
        camera_.setGroupID("rendering");
    addProperty(autoadjustRadius_);
        autoadjustRadius_.setGroupID("rendering");
        autoadjustRadius_.setReadOnlyFlag(true);
    addProperty(radius_);
        radius_.setGroupID("rendering");
        radius_.setStepping(0.001f);
        radius_.setNumDecimals(4);
    addProperty(radiusDivieder_);
        radiusDivieder_.setGroupID("rendering");
        radiusDivieder_.setReadOnlyFlag(true);
    setPropertyGroupGuiName("rendering", "Rendering");

    addProperty(selectionMode_);
    selectionMode_.addOption("add",     "Add to selection",      POIStorage::ADD_TO_SELECTION);
    selectionMode_.addOption("replace", "Replace selection",     POIStorage::REPLACE_SELECTION);
    selectionMode_.addOption("toggle",  "Toggle selection",      POIStorage::TOGGLE_SELECTION);
    selectionMode_.addOption("remove",  "Remove from selection", POIStorage::REMOVE_SELECTION);

    addProperty(distanceMeasureColor_);
        distanceMeasureColor_.setGroupID("info");
    setPropertyGroupGuiName("info", "Point information");


    addProperty(shader_);
    shader_.setGroupID("shaders");
    setPropertyGroupGuiName("shaders", "Shaders");


    cameraHandler_ = std::unique_ptr<CameraInteractionHandler>
        (new CameraInteractionHandler("cameraHandler", "Camera Handler", &camera_, true));

    mouseMoveEventProp_ = std::unique_ptr<EventProperty<POIRenderer3d> >(new EventProperty<POIRenderer3d>("mouseMoveEventProp", "Move", this,
        &POIRenderer3d::mouseMoveEvent,
        tgt::MouseEvent::MOUSE_BUTTON_NONE,
        tgt::MouseEvent::ACTION_ALL,
        tgt::Event::MODIFIER_NONE, false, true));

    mouseClickEventProp_ = std::unique_ptr<EventProperty<POIRenderer3d> >(new EventProperty<POIRenderer3d>("mouseClickEventProp", "Click",
        this, &POIRenderer3d::mouseClickEvent,
        static_cast<tgt::MouseEvent::MouseButtons>(tgt::MouseEvent::MOUSE_BUTTON_LEFT),
        tgt::MouseEvent::CLICK, tgt::Event::MODIFIER_NONE, true));

    mousePressEventProp_ = std::unique_ptr<EventProperty<POIRenderer3d> >(new EventProperty<POIRenderer3d>("mousePressEventProp1", "Click",
       this, &POIRenderer3d::startSelecting,
       static_cast<tgt::MouseEvent::MouseButtons>(tgt::MouseEvent::MOUSE_BUTTON_LEFT),
       tgt::MouseEvent::PRESSED, tgt::Event::ALT, true));

    mouseReleasedEventProp_ = std::unique_ptr<EventProperty<POIRenderer3d> >(new EventProperty<POIRenderer3d>("mouseReleaseEventProp2", "Click",
        this, &POIRenderer3d::endSelecting,
        static_cast<tgt::MouseEvent::MouseButtons>(tgt::MouseEvent::MOUSE_BUTTON_LEFT),
        tgt::MouseEvent::RELEASED, tgt::Event::MODIFIER_NONE, true));
    

    addEventProperty(mouseMoveEventProp_.get());
    addEventProperty(mouseClickEventProp_.get());
    addEventProperty(mousePressEventProp_.get());
    addEventProperty(mouseReleasedEventProp_.get());

    shouldRebuildBuffers_ = true;
    ON_CHANGE_LAMBDA(cpPort_, [this]{
        shouldRebuildBuffers_ = true;
    });

    ON_CHANGE_LAMBDA(volumeInport_, [this]{
        bool hasData = volumeInport_.getData() != nullptr;
        if (!hasData && autoadjustRadius_.get())
            autoadjustRadius_.set(false);
        autoadjustRadius_.setReadOnlyFlag(!hasData);
        radius_.setReadOnlyFlag(autoadjustRadius_.get());
        autodjustRadius();
    });

    ON_CHANGE_LAMBDA(autoadjustRadius_, [this]{
        radius_.setReadOnlyFlag(autoadjustRadius_.get());
        radiusDivieder_.setReadOnlyFlag(!autoadjustRadius_.get());
        autodjustRadius();
    });

    ON_CHANGE_LAMBDA(radiusDivieder_, [this]{
        autodjustRadius();
    });
    addInteractionHandler(cameraHandler_.get());

    selecting_ = false;
    vao_ = 0;
    vbo_ = 0;

}

Processor* POIRenderer3d::create() const {
    return new POIRenderer3d();
}

std::string POIRenderer3d::getClassName() const {
    return "POIRenderer3d";
}

std::string POIRenderer3d::getCategory() const {
    return "Points of Interest";
}

void POIRenderer3d::setDescriptions() {
    setDescription("The POIRenderer3d renders points as 3d spheres. It is useful to blend it over a 3d "
                   "rendering of the volume which is beeing annotated with points.");
    autoadjustRadius_.setDescription("If a volume is supplied this enables calculating the radius of the spheres from its diameter.");
    radius_.setDescription("The radius of each sphere in world space.");
    radiusDivieder_.setDescription("If the radius is set to adjust to a supplied volume, the dimension is the volumes diameter divided by this property.");
    distanceMeasureColor_.setDescription("If a point is active and the mouse hovers "
                                         "over a different point a tool to show this distance is drawn. This is its color.");
}

void POIRenderer3d::process() {

    if (shouldRebuildBuffers_ || true){
        buildBuffers();
        shouldRebuildBuffers_ = false;
    }

    POIStorage* storage = cpPort_.getConnectedProcessor();

    tgt::Shader* shader;
    tgt::mat4 MV = camera_.get().getViewMatrix();
    tgt::mat4 P = camera_.get().getProjectionMatrix(renderPort_.getSize());

    // render image
    shader = shader_.getShader();
    shader->activate();
    shader->setUniform("MV_", MV, false);
    shader->setUniform("P_", P, false);
    shader->setUniform("radius_", radius_.get());
    renderPort_.activateTarget();
    renderPort_.clearTarget();
    renderPOIs();
    shader->deactivate();

    // render distance measure
    POIPointID mouseoverid = storage->getMouseOverPoint();
    POIPointID activeid = POI_NO_SUCH_POINT;
    if (storage->selectedPointCount() == 1) {
        activeid = storage->getSelectedPoints().at(0).id_;
    }
    if (mouseoverid != POI_NO_SUCH_POINT && activeid != POI_NO_SUCH_POINT){
        glDisable(GL_DEPTH_TEST);
        tgt::vec3 mouseover_world = storage->getPointById(mouseoverid).position_;
        tgt::vec3 active_world = storage->getPointById(activeid).position_;

        tgt::mat4 inv;
        (P*MV).invert(inv);

        tgt::vec4 viewdir4 = inv*tgt::vec4(0, 0, 1, 1);
        tgt::vec3 viewdir = viewdir4.xyz()/viewdir4.w;
        tgt::vec3 orth = normalize(tgt::cross(viewdir, mouseover_world-active_world));

        MatStack.pushMatrix();
        MatStack.loadIdentity();
        MatStack.multMatrix(P*MV);

        IMode.begin(tgt::ImmediateMode::LINES);
        IMode.color(distanceMeasureColor_.get());
        IMode.vertex(active_world);
        IMode.vertex(mouseover_world);

        IMode.vertex(active_world+2*radius_.get()*orth);
        IMode.vertex(active_world-2*radius_.get()*orth);

        IMode.vertex(mouseover_world+1.5f*radius_.get()*orth);
        IMode.vertex(mouseover_world-1.5f*radius_.get()*orth);
        IMode.color(tgt::vec4::one);
        IMode.end();


        MatStack.popMatrix();
        glEnable(GL_DEPTH_TEST);
    }

    shader->deactivate();
    if (selecting_){
        renderSelection(selectionBegin_, selectionEnd_);
    }
    renderPort_.deactivateTarget();
    
}

void POIRenderer3d::initialize() {
    RenderProcessor::initialize();
    shader_.rebuild();
}

bool POIRenderer3d::isReady() const
{
    return cpPort_.isReady() && renderPort_.isReady();
}

void POIRenderer3d::buildBuffers()
{
    if (vbo_){
        glDeleteBuffers(1, &vbo_);
        vbo_ = 0;
    }
    if (vao_){
        glDeleteVertexArrays(1, &vao_);
        vao_ = 0;
    }

    POIStorage* storage = cpPort_.getConnectedProcessor();
    POIPointID mouseoverid = storage->getMouseOverPoint();
    const std::vector<POIPoint>& points = storage->getPoints();
    glGenVertexArrays(1, &vao_);
    glBindVertexArray(vao_);

    std::vector<VertexLayout3D> verts;
    for(auto p: points){
        POIGroup g = storage->getGroup(p.group_);
        if (!g.enabled_)
            continue;

        VertexLayout3D v;
        v.pos = p.position_;
        v.color = g.color_;
        v.selected = p.selected_ ? 1.0f : 0.0f;
        v.selected += p.id_ == mouseoverid ? 2.0f : 0.0f;
        verts.push_back(v);
    }

    glGenBuffers(1, &vbo_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(VertexLayout3D)*verts.size(), verts.data(), GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VertexLayout3D), (void*)offsetof(VertexLayout3D, pos));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(VertexLayout3D), (void*)offsetof(VertexLayout3D, color));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(VertexLayout3D), (void*)offsetof(VertexLayout3D, selected));
    glBindVertexArray(0);

    primitiveCount_  = static_cast<int>(verts.size());
}

void POIRenderer3d::renderPOIs()
{
    glBindVertexArray(vao_);
    
    glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(primitiveCount_));
    glBindVertexArray(0);
}

void POIRenderer3d::mouseMoveEvent(tgt::MouseEvent* mouseEve)
{
    tgt::ivec2 pos = tgt::ivec2(mouseEve->x(), renderPort_.getSize().y-mouseEve->y());
    
    POIPointID id = idAtPosition(pos);
    POIStorage* storage = cpPort_.getConnectedProcessor();
    storage->setMouseOverPoint(id);
    if (selecting_){
        selectionEnd_ = tgt::vec2(pos);
    }
    mouseEve->accept();
}

void POIRenderer3d::mouseClickEvent(tgt::MouseEvent* mouseEve)
{
    tgt::ivec2 pos = tgt::ivec2(mouseEve->x(), renderPort_.getSize().y-mouseEve->y());
    POIPointID id = idAtPosition(pos);
    if (id != POI_NO_SUCH_POINT){
        POIStorage* storage = cpPort_.getConnectedProcessor();
        POIPoint p = storage->getPointById(id);
        p.selected_ = true;
        storage->setPointById(p);
    }

}

static tgt::dvec3 transformVec3(tgt::dmat4 M, tgt::dvec3 v) {
    tgt::vec4 p = M*tgt::vec4(v, 1.0);
    return p.xyz()/p.w;
}

static double intersection(tgt::dvec3 origin, tgt::dvec3 dir, tgt::dvec3 pos, double radius){
    tgt::dvec3 diff =origin-pos;
    double ff = tgt::dot(diff, dir);
    double sq = ff*ff-tgt::dot(diff, diff)+radius*radius;
    if (sq < 0)
        return -1;
    sq = sqrt(sq);
    double minf = -ff-sq;
    double maxf = -ff+sq;
    if (minf < 0 && maxf > 0)
        return 0.0001;
    return minf;
}

voreen::POIPointID POIRenderer3d::idAtPosition(tgt::ivec2 pos)
{
    POIStorage* storage = cpPort_.getConnectedProcessor();

    tgt::dvec3 origin;
    tgt::dvec3 dir;

    tgt::Camera cam = camera_.get();

    tgt::dmat4 MVP = cam.getProjectionMatrix(renderPort_.getSize())*cam.getViewMatrix();
    tgt::dmat4 MVPi;
    MVP.invert(MVPi);

    tgt::dmat4 ortho = tgt::mat4::createOrtho(0, 1.0f*renderPort_.getSize().x, 0, 1.0f*renderPort_.getSize().y, -1, 1);
    tgt::dmat4 orthoI;
    ortho.invert(orthoI);

    tgt::dmat4 transform = MVPi*ortho;
    tgt::dvec3 f = transformVec3(orthoI, tgt::vec3(pos, 0.0));

    tgt::dvec3 begin = transformVec3(transform, tgt::vec3(pos, -1.0));
    tgt::dvec3 end = transformVec3(transform, tgt::vec3(pos, 1.0));
    
    origin = begin;
    dir = tgt::normalize(end-begin);

    POIPointID id = POI_NO_SUCH_POINT;
    double dist = -10000000000;
    float radius = radius_.get();
    for(auto p : storage->getPoints()){
        double intersect = intersection(origin, dir, p.position_, radius);
        if (intersect > 0 && intersect > dist){
            dist = intersect;
            id = p.id_;
        }
    }
    return id;
}

void POIRenderer3d::autodjustRadius()
{
    const VolumeBase* v = volumeInport_.getData();
    if (autoadjustRadius_.get() && v){
        float len = tgt::length(v->getBoundingBox(true).getBoundingBox().diagonal());
        radius_.set(len/radiusDivieder_.get());
    }
}

void POIRenderer3d::endSelecting( tgt::MouseEvent* ev )
{
    POIStorage* storage = cpPort_.getConnectedProcessor();
    if (selecting_){
        selecting_ = false;
        tgt::Bounds b(tgt::vec3(selectionBegin_, -10000.0f), tgt::vec3(selectionEnd_, 10000.0f));
        tgt::vec2 llf = b.getLLF().xy();
        tgt::vec2 urb = b.getURB().xy();
        tgt::Camera cam =camera_.get();

        tgt::mat4 MVP = cam.getProjectionMatrix(renderPort_.getSize())*cam.getViewMatrix();
        tgt::mat4 ortho = tgt::mat4::createOrtho(0, 1.0f*renderPort_.getSize().x, 0, 1.0f*renderPort_.getSize().y, -1, 1);
        tgt::mat4 orthoI;
        ortho.invert(orthoI);

        tgt::mat4 transform = orthoI*MVP;
        std::vector<POIPoint> selectedPoints;
        for(auto p: storage->getPoints()){
            tgt::vec4 pos4Proj = transform*tgt::vec4(p.position_, 1.0f);
            tgt::vec3 posProj = pos4Proj.xyz()/pos4Proj.w;
            if (b.containsPoint(posProj)){
                selectedPoints.push_back(p);
            }
        }
        storage->selectPoints(selectedPoints, selectionMode_.getValue());
    }
}

void POIRenderer3d::startSelecting( tgt::MouseEvent* ev )
{
    tgt::ivec2 pos = tgt::ivec2(ev->x(), renderPort_.getSize().y-ev->y());
    if (ev->modifiers() != mousePressEventProp_->getModifier()){
        return;
    }

    POIPointID id = idAtPosition(pos);
    if (id == POI_NO_SUCH_POINT){
        selecting_ = true;
        selectionBegin_ = tgt::vec2(pos);
        selectionEnd_   = tgt::vec2(pos);
    }
}

void POIRenderer3d::renderSelection( tgt::vec2 begin, tgt::vec2 end )
{
    MatStack.pushMatrix();
    tgt::Camera cam = camera_.get();

    MatStack.loadIdentity();
    MatStack.ortho(0, 1.0f*renderPort_.getSize().x, 0, 1.0f*renderPort_.getSize().y, -1, 1);

    glDisable(GL_DEPTH_TEST);
    IMode.begin(tgt::ImmediateMode::LINE_LOOP);
        IMode.vertex(begin);
        IMode.vertex(begin.x,end.y);
        IMode.vertex(end.x,end.y);
        IMode.vertex(end.x,begin.y);
    IMode.end();
    IMode.color(tgt::vec4::one);
    glEnable(GL_DEPTH_TEST);
    MatStack.popMatrix();
}
} // namespace
