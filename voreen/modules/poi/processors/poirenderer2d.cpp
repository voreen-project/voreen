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

#include "poirenderer2d.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/utils/stringutils.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "tgt/textureunit.h"
#include "tgt/immediatemode/immediatemode.h"
#include "tgt/event/mouseevent.h"
#include "voreen/core/voreenapplication.h"

#include <cmath>
#include <cstdlib>
namespace voreen {
namespace{

struct VertexLayout
{
    tgt::vec3 pos;
    tgt::vec3 color;
    float active;
    float id;
};
int popcnt(int p){
#ifdef __linux__ 
    return __builtin_popcount(p);
#elif _WIN32
    return __popcnt(p);
#else
    int c = 0;
    while(p != 0){
        c += p&1;
        p>>=1;
    }
    return c;
#endif
}
};
const std::string POIRenderer2d::loggerCat_("voreen.POIRenderer2d");

POIRenderer2d::POIRenderer2d()
    : ImageProcessor("image/background")
    , inport_            (Port::INPORT,  "inport",            "Slice Inport")
    , volumeInport_ (Port::INPORT, "volumeInport", "Volume Inport")
    , cpPort_(Port::INPORT, "cpPort", "Coprocessors", false)
    , renderOutport_     (Port::OUTPORT, "renderOutport",     "Volume Outport")
    , radius_("radius", "Radius in Pixel", 0.5f, 0.0f, 100.0f)
    , showArrows_("showArrows", "Show Arrows", false)
    , arrowRange_("arrowRange", "Max slice distance for arrow", 10, 0, 10000)
    , shader_("shader", "Shader", "poipoint2d.frag", "poipoint2d.vert", "poipoint2d.geom", Processor::INVALID_PROGRAM, Property::LOD_DEBUG)
    , arrowShader_("arrowShader", "Arrow Shader", "poipoint2darrow.frag", "poipoint2darrow.vert", "poipoint2darrow.geom", Processor::INVALID_PROGRAM, Property::LOD_DEBUG)
    , showIds_("showIds", "Show Ids", false)
    , fontProp_("fontProp", "Font for ids")
    , distanceMeasureColor_("distanceMeasureColor", "Color for distance measurement")
    , selectionMode_("selectionMode", "Selection Mode")
    , selectBackgroud_("selectBackgroud", "Select points in other slices.")
    , linkedMousePositionInSlice_("linkedMousePositionInSlice","\"Mouse Position\" from SliceViewer",tgt::ivec3(-1),tgt::ivec3(-1),tgt::ivec3(std::numeric_limits<int>::max()-1),
                                  Processor::INVALID_RESULT,NumericProperty<tgt::ivec3>::STATIC,Property::LOD_DEBUG)
    , linkedSliceAlignment_("linkedSliceAlignment", "\"Slice Alignment\" from SliceViewer",Processor::INVALID_RESULT,false,Property::LOD_DEBUG)
    , linkedSliceIndex_("linkedSliceIndex", "\"Slice Number\" from SliceViewer", 0, 0, 10000,Processor::INVALID_RESULT,NumericProperty<int>::STATIC, Property::LOD_DEBUG)
    , linkedPickingMatrix_("linkedPickingMatrix", "\"Picking Matrix\" from SliceViewer", tgt::mat4::createIdentity(), tgt::mat4(-1e6f), tgt::mat4(1e6f),
                            Processor::INVALID_RESULT, NumericProperty<tgt::mat4>::STATIC,Property::LOD_DEBUG)
{
    addPort(inport_);
    addPort(volumeInport_);
    addPort(cpPort_);
    addPort(renderOutport_);

    addProperty(radius_);
    
    addProperty(showArrows_);
        showArrows_.setGroupID("arrows");
    addProperty(arrowRange_);
        arrowRange_.setGroupID("arrows");
    setPropertyGroupGuiName("arrows", "Arrows to other slices");

    
    addProperty(showIds_);
        showIds_.setGroupID("info");
    addProperty(fontProp_);
        fontProp_.setGroupID("info");
    setPropertyGroupGuiName("info", "Point information");
    
    addProperty(distanceMeasureColor_);
        distanceMeasureColor_.setGroupID("distance");
    setPropertyGroupGuiName("distance", "Distance mesurement");


    addProperty(selectionMode_);
        selectionMode_.addOption("add",     "Add to selection",      POIStorage::ADD_TO_SELECTION);
        selectionMode_.addOption("replace", "Replace selection",     POIStorage::REPLACE_SELECTION);
        selectionMode_.addOption("toggle",  "Toggle selection",      POIStorage::TOGGLE_SELECTION);
        selectionMode_.addOption("remove",  "Remove from selection", POIStorage::REMOVE_SELECTION);
        selectionMode_.setGroupID("selection");   
    addProperty(selectBackgroud_);
        selectBackgroud_.setGroupID("selection");   
    setPropertyGroupGuiName("selection", "Selection");

    addProperty(shader_);
        shader_.setGroupID("shaders");
    addProperty(arrowShader_);
        arrowShader_.setGroupID("shaders");

    addProperty(linkedMousePositionInSlice_);
        linkedMousePositionInSlice_.setGroupID("linked");
    addProperty(linkedSliceAlignment_);
        linkedSliceAlignment_.addOption("xy-plane", "XY-Plane (axial)", XY_PLANE);
        linkedSliceAlignment_.addOption("xz-plane", "XZ-Plane (coronal)", XZ_PLANE);
        linkedSliceAlignment_.addOption("yz-plane", "YZ-Plane (sagittal)", YZ_PLANE);
        linkedSliceAlignment_.setGroupID("linked");
    addProperty(linkedSliceIndex_);
        linkedSliceIndex_.setGroupID("linked");
    addProperty(linkedPickingMatrix_);
        linkedPickingMatrix_.setReadOnlyFlag(true);
        linkedPickingMatrix_.setGroupID("linked");
    setPropertyGroupGuiName("linked", "To be linked");

    mouseEventPress_ = std::unique_ptr<EventProperty<POIRenderer2d>>(new EventProperty<POIRenderer2d>("mouseEvent.cursorPositionPress1", "Cursor Position Press",
        this, &POIRenderer2d::mouseLocalization,
        static_cast<tgt::MouseEvent::MouseButtons>(tgt::MouseEvent::MOUSE_BUTTON_LEFT),
        tgt::MouseEvent::PRESSED | tgt::MouseEvent::RELEASED | tgt::MouseEvent::MOTION, tgt::Event::MODIFIER_NONE, true));

    mouseEventSet_ = std::unique_ptr<EventProperty<POIRenderer2d>>( new EventProperty<POIRenderer2d>("mouseEvent.set4", "Cursor Position Press",
        this, &POIRenderer2d::setPoint,
        static_cast<tgt::MouseEvent::MouseButtons>(tgt::MouseEvent::MOUSE_BUTTON_LEFT),
        tgt::MouseEvent::RELEASED, tgt::Event::SHIFT, true));

    removeEvent_ = std::unique_ptr<EventProperty<POIRenderer2d>>(new EventProperty<POIRenderer2d>("keyEvent.remove", "Remove button event", this, &POIRenderer2d::removePoint, tgt::KeyEvent::K_DELETE,
        tgt::KeyEvent::MODIFIER_NONE, false, true));

    addEventProperty(mouseEventPress_.get());
    addEventProperty(mouseEventSet_.get());
    addEventProperty(removeEvent_.get());

    shouldRebuildPointsPerSlice_ = true;
    ON_CHANGE_LAMBDA(cpPort_, [this]{
        shouldRebuildPointsPerSlice_ = true;
    });

    ON_CHANGE_LAMBDA(linkedSliceAlignment_, [this]{
        shouldRebuildPointsPerSlice_ = true;
    });

    ON_CHANGE_LAMBDA(volumeInport_, [this]{
        shouldRebuildPointsPerSlice_ = true;
    });


    ON_CHANGE_LAMBDA(linkedMousePositionInSlice_, [this]{
        POIStorage* storage = cpPort_.getConnectedProcessor();
        if (!storage) return;
        POIPointID id = FindPOI(linkedMousePositionInSlice_.get());
        storage->setMouseOverPoint(id);
    });

    isDragging_ = false;
    isSelecting_ = false;
}

POIRenderer2d::~POIRenderer2d() {
}

bool POIRenderer2d::isReady() const{
    return inport_.isReady() && renderOutport_.isReady() && volumeInport_.isReady();
}

void POIRenderer2d::process(){
    renderOutport_.activateTarget();
    glBindTexture(GL_TEXTURE_2D, inport_.getRenderTarget()->getColorTexture()->getId());
    renderQuad();
    
    if (isDragging_){
        POIStorage* storage = cpPort_.getConnectedProcessor();
        
        tgt::vec3 mousepos_voxel = linkedMousePositionInSlice_.get();
        tgt::vec3 pos_world = (volumeInport_.getData()->getVoxelToWorldMatrix()*tgt::vec4(mousepos_voxel-dragOffset_, 1.0)).xyz();
        POIPoint point = storage->getPointById(draggedPoint_);
        point.position_ = pos_world;
        storage->setPointById(point);
    }
    
    if (shouldRebuildPointsPerSlice_){
        buildPointsPerSlice();
        shouldRebuildPointsPerSlice_ = false;
    }
    renderPOIs();
     
    if (isSelecting_){
        tgt::mat4 screenToNDC = tgt::mat4::createOrtho(0.0f, 1.0f*renderOutport_.getSize().x, 0.0f, 1.0f*renderOutport_.getSize().y, -10000.0f, 10000.0f);
        tgt::mat4 voxelToScreen =  voxelToScreenMatrix();
        tgt::vec3 mousepos_voxel = linkedMousePositionInSlice_.get();
        tgt::Bounds selection(selectionEnd_voxel_, selectionStart_voxel_);
        renderSelection(selection);
    }
    renderOutport_.deactivateTarget(); 

}

void POIRenderer2d::setDescriptions(){
    setDescription("This is the primary processor for interacting with point data and it "
                   "shows points in a slice as colored circles. It is an overlay over a SliceViewer. ");
    radius_.setDescription(" The radius of the circles drawn on the screen");

    showArrows_.setDescription(" If this is active, points in other slices are shown as up and down "
                               "pointing arrows. The direction is the direction in which the mousewheel needs " 
                               "to be turned to scroll to the slice of the corresponding point in usual computer "
                               "configurations");
    arrowRange_.setDescription("This sets in how many of the next slices the "
                               "arrows are drawn. Points in more distant slices are not visible in any way");

    showIds_.setDescription("This allows to draw a number over every circle with its id number.");
    fontProp_.setDescription("The font used to draw the ids over the circles.");
    distanceMeasureColor_.setDescription("If a point is active and the mouse hovers "
                                         "over a different point a tool to show this distance is drawn. This is its color.");
    selectionMode_.setDescription("Lets the used change the mode on selection. Points can be toggled, added to the selection "
                                  "or replace the selection.");

    linkedMousePositionInSlice_.setDescription("Needs to be linked to mousePositionInSlice in the input sliceviewer.");
    linkedSliceAlignment_.setDescription("Needs to be linked to sliceAlignment in the input sliceviewer.");
    linkedSliceIndex_.setDescription("Needs to be linked to sliceIndex in the input sliceviewer.");
    linkedPickingMatrix_.setDescription("Needs to be linked to pickingMatrix in the input sliceviewer.");
}

void POIRenderer2d::setRenderState(){
    glDisable(GL_DEPTH_TEST);

}

void POIRenderer2d::resetRenderState(){
    glEnable(GL_DEPTH_TEST);
}

void POIRenderer2d::createSliceBaseVectors(tgt::vec3 &basex, tgt::vec3 &basey){
    tgt::mat4 voxelToWorldMatrix = volumeInport_.getData()->getVoxelToWorldMatrix();
    switch (linkedSliceAlignment_.getValue())
    {
    case XY_PLANE:
        basex = tgt::vec3(1, 0, 0);
        basey = tgt::vec3(0, 1, 0);
        break;
    case XZ_PLANE:
        basex = tgt::vec3(1, 0, 0);
        basey = tgt::vec3(0, 0, 1);
        break;
    case YZ_PLANE:
        basex = tgt::vec3(0, 1, 0);
        basey = tgt::vec3(0, 0, 1);
        break;
    default:
        LERROR("Unknown slice alignment!");
        return;
    }
    basex = (voxelToWorldMatrix*tgt::vec4(basex, 0.0)).xyz();
    basey = (voxelToWorldMatrix*tgt::vec4(basey, 0.0)).xyz();
}

void POIRenderer2d::mouseLocalization(tgt::MouseEvent* ev){
    POIStorage* storage = cpPort_.getConnectedProcessor();
    tgt::mat4 screenToVoxelMatrix = linkedPickingMatrix_.get();

    float screenHeight = static_cast<float>(renderOutport_.getSize().y);

    tgt::vec3 mousepos_voxel = (screenToVoxelMatrix*
        tgt::vec4(static_cast<float>(ev->x()), screenHeight-static_cast<float>(ev->y()), static_cast<float>(linkedSliceIndex_.get()), 1.0f)).xyz();
    POIPointID hitid = FindPOI(mousepos_voxel);;
    POIPoint hitpoint;
    if (hitid != POI_NO_SUCH_POINT) hitpoint = storage->getPointById(hitid);


    storage->setMouseOverPoint(hitid);
    tgt::mat4 worldtovoxel = volumeInport_.getData()->getWorldToVoxelMatrix();

    if (isSelecting_){
        float layers = 0.5f+ (selectBackgroud_.get()? static_cast<float>(arrowRange_.get()): 0.0f);
        tgt::vec3 mousepos_layer_voxel = (screenToVoxelMatrix*
            tgt::vec4(static_cast<float>(ev->x()), screenHeight-static_cast<float>(ev->y()), static_cast<float>(linkedSliceIndex_.get())+layers, 1.0f)).xyz();
        selectionEnd_voxel_ = mousepos_layer_voxel;
    }

    if (storage->isInteractionEnabled()){
        
        if (ev->action() == tgt::MouseEvent::PRESSED && ev->button() == tgt::MouseEvent::MOUSE_BUTTON_LEFT 
            && hitid != POI_NO_SUCH_POINT){
            // Selected a Point
            POIPoint p = storage->getPointById(hitid);
            if (!isDragging_){
                std::vector<POIPoint> selectedPoints;
                selectedPoints.push_back(p);
                storage->selectPoints(selectedPoints, selectionMode_.getValue());
            }
            isDragging_ = true;
            draggedPoint_ = hitid;
            dragOffset_ = mousepos_voxel-(worldtovoxel*tgt::vec4(hitpoint.position_, 1.0)).xyz();
            ev->accept();
        }else if (ev->action() == tgt::MouseEvent::PRESSED && ev->button() == tgt::MouseEvent::MOUSE_BUTTON_LEFT && !isSelecting_){
            // Start selection of area
            isSelecting_ = true;
            float layers = 0.5f+ (selectBackgroud_.get()? static_cast<float>(arrowRange_.get()): 0.0f);
            tgt::vec3 mousepos_layer_voxel = (screenToVoxelMatrix*
                tgt::vec4(static_cast<float>(ev->x()), screenHeight-static_cast<float>(ev->y()), static_cast<float>(linkedSliceIndex_.get())-layers, 1.0f)).xyz();
            selectionStart_voxel_ = mousepos_layer_voxel;
            selectionEnd_voxel_ = mousepos_layer_voxel;
        }else if(ev->action() == tgt::MouseEvent::RELEASED){
            // 
            isDragging_ = false;
            if (isSelecting_){
                float layers = 0.5f+(selectBackgroud_.get()? static_cast<float>(arrowRange_.get()): 0.0f);
                tgt::vec3 mousepos_layer_voxel = (screenToVoxelMatrix*
                    tgt::vec4(static_cast<float>(ev->x()), screenHeight-static_cast<float>(ev->y()), static_cast<float>(linkedSliceIndex_.get())+layers, 1.0f)).xyz();

                isSelecting_ = false;
                tgt::Bounds selectionRect(mousepos_layer_voxel, selectionStart_voxel_);
                auto points = storage->getPoints();
                std::vector<POIPoint> selectedPoints;
                for(auto p : points){
                    tgt::vec3 pos_voxel = (worldtovoxel*tgt::vec4(p.position_, 1.0)).xyz();
                    if (selectionRect.containsPoint(pos_voxel)){
                        selectedPoints.push_back(p);
                    }
                }
                storage->selectPoints(selectedPoints, selectionMode_.getValue());
            }
        }
    }
}

void POIRenderer2d::setPoint(tgt::MouseEvent* ev)
{
    POIStorage* storage = cpPort_.getConnectedProcessor();
    if (!storage->isInteractionEnabled())
        return;

    const VolumeBase * vol = volumeInport_.getData();
    if (!vol)
        return;

    tgt::ivec2 mouse_screen = ev->coord();
    mouse_screen.y = renderOutport_.getSize().y-mouse_screen.y;
    
    int sliceIndex = linkedSliceIndex_.get();

    tgt::mat4 screenToVoxelMatrix = linkedPickingMatrix_.get();
    tgt::mat4 voxelToWorldMatrix = volumeInport_.getData()->getVoxelToWorldMatrix();
    
    tgt::vec3 mouse_voxel = (screenToVoxelMatrix * tgt::vec4(mouse_screen, static_cast<float>(sliceIndex), 1.f)).xyz();
    tgt::vec3 mouse_world = (voxelToWorldMatrix * tgt::vec4(mouse_voxel, 1.f)).xyz();
    tgt::ivec3 dim = vol->getDimensions();

    POIGroupID activeGroup = storage->getActiveGroup();
    if (activeGroup == POI_NO_SUCH_GROUP){
        //VoreenApplication::app()->showMessageBox("No group selected", "A group needs to be selected to place a new point!", true);
        LERROR("No Group selected");
    }else if (mouse_voxel.x >= 0 && mouse_voxel.y >= 0 && mouse_voxel.z > 0 &&
        mouse_voxel.x < dim.x && mouse_voxel.y < dim.y && mouse_voxel.z < dim.z){
        storage->addPoint(mouse_world, activeGroup);
    }
    ev->accept();
}

void POIRenderer2d::buildPointsPerSlice()
{
    const VolumeBase * vol = volumeInport_.getData();
    if (!vol)
        return;
    POIStorage* storage = cpPort_.getConnectedProcessor();
    if (!storage)
        return;

    int comp = 0;
    switch (linkedSliceAlignment_.getValue())
    {
    case XY_PLANE:
        comp = 2;
        break;
    case XZ_PLANE:
        comp = 1;
        break;
    case YZ_PLANE:
        comp = 0;
        break;
    }

    pointsPerSlice_.clear();
    pointsPerSlice_.resize(vol->getDimensions()[comp]);

    tgt::vec4 worldToVoxel = vol->getWorldToVoxelMatrix()[comp];
    
    for(POIPoint p: storage->getPoints()){
        // FIXME: Does this work for translated world positions? Maybe not
        int slice = static_cast<int>(std::floor(dot(worldToVoxel, tgt::vec4(p.position_, 1.0))+0.5f));
        if (slice >= 0 && slice < (int)pointsPerSlice_.size()){
            pointsPerSlice_[slice].push_back(p);
        }
    }
}

void POIRenderer2d::renderPOIs()
{
    setRenderState();

    tgt::vec3 basex, basey;
    POIStorage* storage = cpPort_.getConnectedProcessor();

    createSliceBaseVectors(basex, basey);
    int sliceIndex = linkedSliceIndex_.get();

    float radius = getRadiusVoxel();
    
    tgt::Shader* shader = shader_.getShader();
    shader->activate();

    tgt::mat4 screenToNDC = tgt::mat4::createOrtho(0.0f, 1.0f*renderOutport_.getSize().x, 0.0f, 1.0f*renderOutport_.getSize().y, -10000.0f, 10000.0f);
    tgt::mat4 worldToScreen =  voxelToScreenMatrix()*
                  volumeInport_.getData()->getWorldToVoxelMatrix();

    shader->setUniform("radius_", radius_.get());
    shader->setUniform("worldToScreen_", worldToScreen);
    shader->setUniform("screenToNDC_", screenToNDC);
    shader->setUniform("selectionid_", static_cast<float>(storage->getMouseOverPoint()));
    
    drawPOISlice(sliceIndex);

    shader->deactivate();

    if (showArrows_.get()){
        shader = arrowShader_.getShader();
        //shader = shader_.getShader();
        shader->activate();
        shader->setUniform("radius_", radius_.get());
        shader->setUniform("worldToScreen_", worldToScreen);
        shader->setUniform("screenToNDC_", screenToNDC);
        shader->setUniform("sign_", -1.0f);
        shader->setUniform("selectionid_", static_cast<float>(storage->getMouseOverPoint()));

        int range = arrowRange_.get();
        for(int i = std::max((int)(sliceIndex-range), (int)0); i != std::min((int)(sliceIndex+range), (int)pointsPerSlice_.size()); i++){
            if (i == sliceIndex)
            {
                shader->setUniform("sign_", 1.0f);
                continue;
            }
            drawPOISlice(i);
        }
    }
    shader->deactivate();

    // render line 
    POIPointID active_id = POI_NO_SUCH_POINT;
    if (storage->selectedPointCount() == 1){
        active_id = storage->getSelectedPoints().at(0).id_;
    }
    POIPointID mouseover_id = storage->getMouseOverPoint();
    POIPoint active_point;
    POIPoint mouseover_point;
    if (active_id != POI_NO_SUCH_POINT) active_point = storage->getPointById(active_id);
    if (mouseover_id != POI_NO_SUCH_POINT) mouseover_point = storage->getPointById(mouseover_id);

    if (active_id != POI_NO_SUCH_POINT && mouseover_id != POI_NO_SUCH_POINT){
        MatStack.pushMatrix();
        MatStack.loadMatrix(screenToNDC);
        tgt::vec2 active_screen = (worldToScreen*tgt::vec4(active_point.position_, 1.0)).xy();
        tgt::vec2 mouseover_screen = (worldToScreen*tgt::vec4(mouseover_point.position_, 1.0)).xy();

        IMode.begin(tgt::ImmediateMode::LINES);
        IMode.color(distanceMeasureColor_.get());
        IMode.vertex(active_screen);
        IMode.vertex(mouseover_screen);

        tgt::vec2 dir = active_screen-mouseover_screen;
        tgt::vec2 orth = tgt::normalize(tgt::vec2(dir.y, -dir.x));

        IMode.vertex(active_screen+orth*radius_.get()*1.5f);
        IMode.vertex(active_screen-orth*radius_.get()*1.5f);

        IMode.vertex(mouseover_screen+orth*radius_.get()*1.5f);
        IMode.vertex(mouseover_screen-orth*radius_.get()*1.5f);

        IMode.end();
        IMode.color(tgt::vec4::one);
        MatStack.popMatrix();
    }
    resetRenderState();

    if (showIds_.get()){
        glDisable(GL_DEPTH_TEST);
        tgt::Font *font = fontProp_.get();
        font->setTextAlignment(tgt::Font::MiddleLeft);
        for(POIPoint p: pointsPerSlice_[sliceIndex]){
            if (!storage->getGroup(p.group_).enabled_) continue;
            std::stringstream ss;
            ss << p.id_;
            tgt::vec3 pos = (worldToScreen*tgt::vec4(p.position_, 1.0)).xyz();
            pos.z = 0;
            tgt::vec2 size = font->getSize(pos, ss.str(), renderOutport_.getSize());
            font->render(pos-tgt::vec3(size.x/2.0f, 0.0f, 0.0f), ss.str(), renderOutport_.getSize());

        }
        glEnable(GL_DEPTH_TEST);
    }
}

void POIRenderer2d::initialize() {
    ImageProcessor::initialize();
    shader_.rebuild();
    arrowShader_.rebuild();
}

tgt::mat4 POIRenderer2d::voxelToScreenMatrix() const
{
    tgt::mat4 pickingMatrix = linkedPickingMatrix_.get();
    tgt::mat4 invPickingMatrix;
    pickingMatrix.invert(invPickingMatrix);
    return invPickingMatrix;
}

void POIRenderer2d::removePoint(tgt::KeyEvent* ev)
{
    POIStorage *storage = cpPort_.getConnectedProcessor();
    if (storage->isInteractionEnabled())
        storage->removeSelectedPoints();
}

void POIRenderer2d::drawPOISlice(int slicenum)
{
    POIStorage* storage = cpPort_.getConnectedProcessor();
    std::vector<VertexLayout> verts;
    //POIPointID active = storage->getActivePoint();
    for(POIPoint p: pointsPerSlice_[slicenum]){
        POIGroup g = storage->getGroup(p.group_);
        if (!g.enabled_)
            continue;
        VertexLayout v;
        v.active = p.selected_ ? 1.0f : 0.0f;
        v.pos    = p.position_;
        v.color  = g.color_;
        v.id     = static_cast<float>(p.id_);
        verts.push_back(v);
    }

    GLuint vao;
    GLuint vbo;
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(VertexLayout)*verts.size(), verts.data(), GL_STREAM_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VertexLayout), (void*)offsetof(VertexLayout, pos));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(VertexLayout), (void*)offsetof(VertexLayout, color));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(VertexLayout), (void*)offsetof(VertexLayout, active));
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(VertexLayout), (void*)offsetof(VertexLayout, id));

    glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(verts.size()));
    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(1, &vbo);
}

float POIRenderer2d::getRadiusVoxel()
{
    tgt::mat4 worldToVoxel = volumeInport_.getData()->getWorldToVoxelMatrix();
    float len = pow(tgt::length(worldToVoxel.getScalingPart()), 1.0f/3.0f);
    return len*radius_.get();
}

bool POIRenderer2d::getRadiusPixel()
{
    tgt::mat4 worldtovoxel = volumeInport_.getData()->getWorldToVoxelMatrix();
    tgt::mat4 voxeltoscreen = voxelToScreenMatrix();
    tgt::mat4 worldtoscreen = voxeltoscreen*worldtovoxel;
    float len = pow(tgt::length(worldtoscreen.getScalingPart()), 1.0f/3.0f);
    return len*radius_.get();
}


voreen::POIPointID POIRenderer2d::FindPOI(tgt::vec3 mousepos_voxel)
{
    int sliceIndex = linkedSliceIndex_.get();

    if (!volumeInport_.getData()) return POI_NO_SUCH_POINT;
    if (sliceIndex >= (int)pointsPerSlice_.size()) return POI_NO_SUCH_POINT;

    float radius = radius_.get();
    tgt::mat4 worldtovoxel = volumeInport_.getData()->getWorldToVoxelMatrix();
    tgt::mat4 voxeltoscreen = voxelToScreenMatrix();
    tgt::mat4 worldtoscreen = voxeltoscreen*worldtovoxel;
    tgt::vec2 mousepos = (voxeltoscreen*tgt::vec4(mousepos_voxel, 1.0)).xy();


    for(POIPoint p: pointsPerSlice_[sliceIndex]){
        tgt::vec2 screenpos = (worldtoscreen*tgt::vec4(p.position_, 1.0)).xy();
        if (tgt::distance(screenpos, mousepos) < radius){
            return p.id_;
        }
    }
    return POI_NO_SUCH_POINT;
}

void POIRenderer2d::renderSelection( tgt::Bounds selection_voxel )
{
    tgt::mat4 screenToNDC = tgt::mat4::createOrtho(0.0f, 1.0f*renderOutport_.getSize().x, 0.0f, 1.0f*renderOutport_.getSize().y, -10000.0f, 10000.0f);

    MatStack.pushMatrix();

    MatStack.loadMatrix(screenToNDC*voxelToScreenMatrix());

    tgt::vec3 llf = selection_voxel.getLLF();
    tgt::vec3 urb = selection_voxel.getURB();

    IMode.begin(tgt::ImmediateMode::LINES);
    
    for(int i = 0; i != 8; i++){
        for(int j = 0; j != 8; j++){
            int hamming_dist = popcnt(j^i); // this is not portable
            if (hamming_dist == 1){
                tgt::vec3 a, b;
                for(int n = 0; n != 3; n++){
                    a[n] = (i&(1<<n))?llf[n]:urb[n];
                    b[n] = (j&(1<<n))?llf[n]:urb[n];
                }
                IMode.vertex(a);
                IMode.vertex(b);
            }
        }
    }
    IMode.end();
    IMode.color(tgt::vec4::one);
    MatStack.popMatrix();
}
} 
