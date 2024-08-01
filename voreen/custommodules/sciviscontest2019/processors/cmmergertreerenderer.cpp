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

//header file
#include "cmmergertreerenderer.h"
#include "tgt/texturemanager.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/properties/eventproperty.h"

#include <deque>
#include <vector>
#include <assert.h>

//we are in namespace voreen
namespace voreen {

CMMergerTreeRenderer::CMMergerTreeRenderer()
    : RenderProcessor()
    , outport_(Port::OUTPORT, "outport","Modified Image", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , inport_(Port::INPORT, "haloport.input", "Halo Data Input")
    , radiusProp_("radiusProp", "Radius", 0.5f, 0.f, 10.f, Processor::INVALID_RESULT, NumericProperty<float>::STATIC, Property::LOD_DEBUG)
    , lineWidthProp_("lineWidthProp", "Line width", 1.0f, 0.f, 5.0f, Processor::INVALID_RESULT, NumericProperty<float>::STATIC, Property::LOD_DEBUG)
    , spacingProp_("spacingProp", "Spacing", 2.0f, 0.f, 10.f, Processor::INVALID_RESULT, NumericProperty<float>::STATIC, Property::LOD_DEBUG)
    , useTexturesProp_("useTexturesProp", "Enable Textures", false)
    , vertShaderProp_("vertShaderProp", "Shaders for vertices", "mergertreevertrenderer.frag", "mergertreevertrenderer.vert", "mergertreevertrenderer.geom", Processor::INVALID_PROGRAM, Property::LOD_DEBUG)
    , vertSelectionShaderProp_("vertSelectionShaderProp", "SelectionShaders for vertices", "mergertreevertselectionrenderer.frag", "mergertreevertselectionrenderer.vert", "mergertreevertselectionrenderer.geom", Processor::INVALID_PROGRAM, Property::LOD_DEBUG)
    , linkShaderProp_("linkShaderProp", "Shaders for links", "mergertreelinkrenderer.frag", "mergertreelinkrenderer.vert", "mergertreelinkrenderer.geom", Processor::INVALID_PROGRAM, Property::LOD_DEBUG)
    , selectedHaloIDProp_("selectedHaloIDProp", "ID of focused halo", 257, -1, 1000000)
    , mouseOverHaloIDProp_("mouseOverHaloIDProp", "ID of hovered over halo", CMMergerTree::NO_HALO_ID, -1, 1000000, Processor::INVALID_RESULT, NumericProperty<int>::STATIC, Property::LOD_DEBUG)
    , timeStep_("timeStep", "time of selected halo", 0.0f, 0.0f, 624.0f)
    , center_("center", "Center", tgt::vec2(0.0f, 0.0f), tgt::vec2(-1000.0f, -1000.0f), tgt::vec2(1000.0f, 1000.0f), Processor::INVALID_RESULT, NumericProperty<tgt::vec2>::STATIC, Property::LOD_DEBUG)
    , windowWidth_("windowWidth", "Width of window", 100.0f, 0.001f, 1000.0f, Processor::INVALID_RESULT, NumericProperty<float>::STATIC, Property::LOD_DEBUG)
    , animationDurationProp_("animationDurationProp", "Duration of camera movement in ms", 1000, 1, 10000, Processor::INVALID_RESULT, NumericProperty<int>::STATIC, Property::LOD_DEBUG)
    , lastOutportWidth_(100)
    , haloTexture_(nullptr)
    , linkTexture_(nullptr)
    , linkTextureZoomLevel_("linkTextureZoomLevel", "Zoom level of link texture", 1.5f, 0.1f, 10.0f, Processor::INVALID_RESULT, NumericProperty<float>::STATIC, Property::LOD_DEBUG)
    , vertColorProp_("vertColor","Color of halos:",tgt::vec4(0.8f))
    , newVertColorProp_("newVertColor","Color of orphan halos:",tgt::vec4(0.0f, 0.8f, 0.8f, 1.0f))
    , selectedVertColorProp_("selectedVertColor","Color of selected halo:",tgt::vec4(0.0f, 0.8f, 0.0f, 1.0f))
    , mouseOverVertColorProp_("mouseOverVertColor","Color of hovered over halo:",tgt::vec4(0.8f, 0.8f, 0.0f, 1.0f))
    , onPathVertColorProp_("onPathVertColor","Color of halo on path:",tgt::vec4(0.8f))
    , linkColorProp_("linkColor","Color of links:",tgt::vec4(0.8f))
    , pathLinkColorProp_("pathLinkColor","Color of path links:",tgt::vec4(0.0f, 0.8f, 0.0f, 1.0f))
    , propertyAnimator_(60)
    , vao_(0)
    , selectionVao_(0)
    , vertVbo_(0)
    , linkIbo_(0)
    , linkCount_(0)
    , targetSelectionLock_(false)
{
    selectionManager_ = new CMSelectionManager("selectionHandler", "Selection", tgt::ivec2(256, 256), &selectedHaloIDProp_, &mouseOverHaloIDProp_, CMMergerTree::NO_HALO_ID);
    addInteractionHandler(selectionManager_);
    cameraHandler_ = new Camera2DInteractionHandler("cameraHandler", "Camera", &center_, &windowWidth_, &lastOutportWidth_);
    addInteractionHandler(cameraHandler_);
    addPort(outport_);
    addPort(inport_);

    //key events
    keyPressEvent_ = new EventProperty<CMMergerTreeRenderer>("CMMergerTreeRenderer.moveToZero", "CMMergerTreeRenderer moveToZero", this, &CMMergerTreeRenderer::processKeyEvent,
        tgt::KeyEvent::K_LAST,
        tgt::Event::MODIFIER_NONE, true, true);
    addEventProperty(keyPressEvent_);


    //register properties

    addProperty(radiusProp_);
    addProperty(lineWidthProp_);
    addProperty(spacingProp_);
    addProperty(useTexturesProp_);

    addProperty(selectedHaloIDProp_);
    addProperty(mouseOverHaloIDProp_);
    addProperty(timeStep_);

    addProperty(linkShaderProp_);
    addProperty(vertShaderProp_);
    addProperty(vertSelectionShaderProp_);
    addProperty(linkTextureZoomLevel_);
    addProperty(vertColorProp_);
    addProperty(newVertColorProp_);
    addProperty(selectedVertColorProp_);
    addProperty(mouseOverVertColorProp_);
    addProperty(onPathVertColorProp_);
    addProperty(linkColorProp_);
    addProperty(pathLinkColorProp_);

    addProperty(windowWidth_);
    addProperty(center_);
    addProperty(animationDurationProp_);

    inport_.onNewData(MemberFunctionCallback<CMMergerTreeRenderer>(this, &CMMergerTreeRenderer::setupBuffers));
    ON_CHANGE(selectedHaloIDProp_, CMMergerTreeRenderer, selectedHaloChanged);
    ON_CHANGE(mouseOverHaloIDProp_, CMMergerTreeRenderer, setupBuffers);
    ON_CHANGE(timeStep_, CMMergerTreeRenderer, timeStepChanged);
    ON_CHANGE(spacingProp_, CMMergerTreeRenderer, setupBuffers);
    ON_CHANGE(useTexturesProp_, CMMergerTreeRenderer, useTexturesChanged);

}

CMMergerTreeRenderer::~CMMergerTreeRenderer(){
    delete cameraHandler_;
    delete selectionManager_;
    delete keyPressEvent_;
    if(haloTexture_) {
        TexMgr.dispose(haloTexture_);
    }
    if(linkTexture_) {
        TexMgr.dispose(linkTexture_);
    }
}

void CMMergerTreeRenderer::initialize() {
    // call superclass function first
    RenderProcessor::initialize();

    selectionManager_->initialize();
    //Fill texture with zeros
    selectionManager_->resize(outport_.getSize());
    selectionManager_->activate();
    selectionManager_->clear();
    selectionManager_->deactivate();

    // load shaders
    useTexturesChanged(); //implicit
    vertSelectionShaderProp_.rebuild();

    haloTexture_ = TexMgr.load(VoreenApplication::app()->getBasePath("custommodules/viscontest2015/resources/halotexture.png"));
    haloTexture_->setFilter(tgt::Texture::Filter::MIPMAP);

    linkTexture_ = TexMgr.load(VoreenApplication::app()->getBasePath("custommodules/viscontest2015/resources/linktexture.png"));
    linkTexture_->setFilter(tgt::Texture::Filter::MIPMAP);
    linkTexture_->setWrapping(tgt::Texture::Wrapping::REPEAT);
    //Synchronize selectedHaloIDProp_ and timeStep_
    //selectedHaloChanged();
}

void CMMergerTreeRenderer::deinitialize() {
    deleteBuffers();
    // call superclass function last
    selectionManager_->deinitialize();
    RenderProcessor::deinitialize();
}

#define VERTEX_FLAG_NONE        0u
#define VERTEX_FLAG_SELECTED   (1u<<0)
#define VERTEX_FLAG_MOUSE_OVER (1u<<1)
#define VERTEX_FLAG_ON_PATH    (1u<<2)
#define VERTEX_FLAG_ORPHAN     (1u<<3)
float CMMergerTreeRenderer::addAncestors(const CMHalo* const currentHalo, std::deque<const CMHalo*>::const_reverse_iterator& path, const std::deque<const CMHalo*>::const_reverse_iterator& pathEnd, std::vector<GLushort>& vertlinks, float& yspace, float xpos, GLushort descVertPos) {
    const float spacing = spacingProp_.get();
    const tgt::vec3 orphanColor(1.0f, 0.0f, 0.0f);
    const tgt::vec3 parentedColor(1.0f, 1.0f, 0.0f);
    GLuint vertex_flags = VERTEX_FLAG_NONE;
    if(currentHalo->parentID == CMMergerTree::NO_HALO_ID) {
        vertex_flags |= VERTEX_FLAG_ORPHAN;
    }
    if(currentHalo->ID == selectedHaloIDProp_.get()) {
        vertex_flags |= VERTEX_FLAG_SELECTED;
    }
    if(currentHalo->ID == mouseOverHaloIDProp_.get()) {
        vertex_flags |= VERTEX_FLAG_MOUSE_OVER;
    }
    if(path != pathEnd && currentHalo->ID == (*path)->ID) {
        vertex_flags |= VERTEX_FLAG_ON_PATH;
        ++path;
    }

    //Secure position in verts array to know the index beforehand
    currentVerts_.emplace_back(tgt::vec3(xpos, yspace, 0), vertex_flags, currentHalo->ID);
    GLushort ownVertPos = currentVerts_.size() - 1;
    if(descVertPos != ownVertPos) {
        vertlinks.push_back(ownVertPos);
        vertlinks.push_back(descVertPos);
    }

    if(currentHalo->parentID == CMMergerTree::NO_HALO_ID) {
        yspace += spacing;
        return yspace - spacing;
    } else {
        float firstParentVertPos, lastParentVertPos;
        const CMHalo* parent = currentHalo->parent();
        firstParentVertPos = lastParentVertPos = addAncestors(parent, path, pathEnd, vertlinks, yspace, xpos-spacing, ownVertPos);
        parent = parent->spouse();
        while(parent) {
            lastParentVertPos = addAncestors(parent, path, pathEnd, vertlinks, yspace, xpos-spacing, ownVertPos);
            parent = parent->spouse();
        }
        float ownYPos = (firstParentVertPos + lastParentVertPos)*0.5f;
        currentVerts_[ownVertPos].pos.y = ownYPos;
        return ownYPos;
    }
}
void CMMergerTreeRenderer::calculateHaloPath() {
    const CMHalo* selectedHalo;
    try {
        selectedHalo = getSelectedHalo();
    } catch(...) {
        return;
    }
    // The selected halo is part of the path already. Nothing to be done here.
    if(std::find(selectedPath_.cbegin(), selectedPath_.cend(), selectedHalo)!=selectedPath_.cend()) {
        return;
    }
    selectedPath_.clear();
    selectedPath_.push_back(selectedHalo);
    const CMHalo* parent = selectedHalo->parent();
    while(parent) {
        const CMHalo* current = parent;
        while(current) {
            if(current->mass > parent->mass) {
                parent = current;
            }
            current = current->spouse();
        }
        selectedPath_.push_front(parent);
        parent = parent->parent();
    }
    const CMHalo* child = selectedHalo->descendant();
    while(child) {
        selectedPath_.push_back(child);
        child = child->descendant();
    }
}

void CMMergerTreeRenderer::deleteBuffers(){
    if (vao_) glDeleteVertexArrays(1, &vao_);
    if (selectionVao_) glDeleteVertexArrays(1, &selectionVao_);
    if (vertVbo_) glDeleteBuffers(1, &vertVbo_);
    if (linkIbo_) glDeleteBuffers(1, &linkIbo_);
    vao_     = 0;
    selectionVao_ = 0;
    vertVbo_ = 0;
    linkIbo_ = 0;
}
void CMMergerTreeRenderer::setupBuffers(){
    deleteBuffers();
    currentVerts_.clear();

    const CMMergerTree* tree = inport_.getData();
    if(!tree) {
        return;
    }

    const CMHalo* selectedHalo = tree->haloByID(selectedHaloIDProp_.get());
    if(!selectedHalo) {
        return;
    }
    calculateHaloPath();


    std::vector<GLushort> vertlinks;
    float yspace = 0.0f;
    auto pathRBegin = selectedPath_.crbegin();
    auto pathREnd = selectedPath_.crend();
    addAncestors(selectedPath_.back(), pathRBegin, pathREnd, vertlinks, yspace, 0.0f);

    linkCount_ = vertlinks.size();

    //Generate and setup buffers
    glGenBuffers(1, &vertVbo_);
    glBindBuffer(GL_ARRAY_BUFFER, vertVbo_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(VertexLayout)*currentVerts_.size(), currentVerts_.data(), GL_STATIC_DRAW);

    glGenBuffers(1, &linkIbo_);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, linkIbo_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLushort)*vertlinks.size(), vertlinks.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    //Generate and setup buffers vao_
    glGenVertexArrays(1, &vao_);
    glBindVertexArray(vao_);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VertexLayout), (void*)offsetof(VertexLayout, pos));

    glEnableVertexAttribArray(1);
    glVertexAttribIPointer(1, 1, GL_UNSIGNED_INT, sizeof(VertexLayout), (void*)offsetof(VertexLayout, flags));

    glBindVertexArray(0);


    //Generate and setup buffers selectionVao_
    glGenVertexArrays(1, &selectionVao_);
    glBindVertexArray(selectionVao_);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VertexLayout), (void*)offsetof(VertexLayout, pos));

    glEnableVertexAttribArray(1);
    glVertexAttribIPointer(1, 1, GL_UNSIGNED_INT, sizeof(VertexLayout), (void*)offsetof(VertexLayout, flags));

    glEnableVertexAttribArray(2);
    glVertexAttribIPointer(2, 1, GL_UNSIGNED_INT, sizeof(VertexLayout), (void*)offsetof(VertexLayout, haloID));

    glBindVertexArray(0);


    glBindBuffer(GL_ARRAY_BUFFER, 0);
}
const CMHalo* CMMergerTreeRenderer::getSelectedHalo() {
    const CMMergerTree* tree = inport_.getData();
    if(!tree) {
        throw "no tree";
    }
    const CMHalo* selectedHalo = tree->haloByID(selectedHaloIDProp_.get());
    if(!selectedHalo) {
        throw "no halo";
    }
    return selectedHalo;
}

void CMMergerTreeRenderer::useTexturesChanged() {
    std::string header = generateHeader();
    if(useTexturesProp_.get()) {
        header += "#define USE_TEXTURES\n";
        linkTextureZoomLevel_.setVisibleFlag(true);
    } else {
        linkTextureZoomLevel_.setVisibleFlag(false);
    }
    vertShaderProp_.setHeader(header);
    vertShaderProp_.rebuild();

    linkShaderProp_.setHeader(header);
    linkShaderProp_.rebuild();
}
void CMMergerTreeRenderer::selectedHaloChanged() {
    const CMHalo* selectedHalo;
    try {
        selectedHalo = getSelectedHalo();
    } catch(...) {
        deleteBuffers();
        return;
    }
    setupBuffers();

    tgt::vec2 newCenter;
    for(VertexLayout& v : currentVerts_) {
        if((int)v.haloID == selectedHalo->ID) {
            newCenter = v.pos.xy();
        }
    }

    targetSelectionLock_ = true;
	float timeStep = ((selectedHalo->scale - (1.0l / 201.0l)) / ((1.0l - (1.0l / 201.0l)) / 625.0l)) - 1.0f;
    timeStep_.set(timeStep);
    //timeStep_.set(selectedHalo->scale); 
    targetSelectionLock_ = false;
    moveCameraTo(newCenter);
}

void CMMergerTreeRenderer::timeStepChanged() {
    if(targetSelectionLock_) {
        return;
    }
    const CMHalo* selectedHalo;
    try {
        selectedHalo = getSelectedHalo();
    } catch(...) {
        deleteBuffers();
        return;
    }
    float targetTimeStep = timeStep_.get();
    auto it = selectedPath_.cbegin();
    const CMHalo* bestCenterHalo = *it++;
    for(; it != selectedPath_.cend(); ++it) {
        if(bestCenterHalo->scale + (*it)->scale > 2*targetTimeStep) {
            break;
        }
        bestCenterHalo = *it;
    }
    //Will call selectedHaloChanged and thus prepare do the rest of the work
    selectedHaloIDProp_.set(bestCenterHalo->ID);
}

void CMMergerTreeRenderer::moveCameraTo(const tgt::vec2& dest) {
    float animationCenterWidth = std::max(tgt::distance(center_.get(), dest), windowWidth_.get());
    std::vector<CMPropertyInterpolator*> props;
    props.push_back(new CMNumericPropertyInterpolator<tgt::vec2>(CMPropertyInterpolator::SMOOTH, &center_, dest));
    props.push_back(new CMTwoStepNumericPropertyInterpolator<float>(CMPropertyInterpolator::LINEAR, &windowWidth_, animationCenterWidth, windowWidth_.get()));
    propertyAnimator_.animate(animationDurationProp_.get(), props);
}
void CMMergerTreeRenderer::processKeyEvent(tgt::KeyEvent* e) {
    static const float kbspeed = 0.01f;
    tgt::vec2 p;
    float newTime;
    bool accepted = true;
    switch(e->keyCode()) {
        case tgt::KeyEvent::K_SPACE:
            selectedHaloChanged();
            break;
        case tgt::KeyEvent::K_J:
            p = center_.get();
            p.y-=kbspeed*windowWidth_.get();
            center_.set(p);
            break;
        case tgt::KeyEvent::K_K:
            p = center_.get();
            p.y+=kbspeed*windowWidth_.get();
            center_.set(p);
            break;
        case tgt::KeyEvent::K_H:
            p = center_.get();
            p.x-=kbspeed*windowWidth_.get();
            center_.set(p);
            break;
        case tgt::KeyEvent::K_L:
            p = center_.get();
            p.x+=kbspeed*windowWidth_.get();
            center_.set(p);
            break;
        case tgt::KeyEvent::K_N:
            newTime = timeStep_.get()+1.00f;
            if(e->pressed() && newTime <= 624.0f){
                timeStep_.set(newTime);
            }
            break;
        case tgt::KeyEvent::K_P:
            newTime = timeStep_.get()-1.00f;
            if(e->pressed() && newTime >= 0.0f){
                timeStep_.set(newTime);
            }
            break;
        default:
            accepted = false;
    }
    if(accepted) {
        e->accept();
    }
}

void CMMergerTreeRenderer::renderLinks(tgt::mat4& projectionMatrix, tgt::mat4& viewMatrix) {
    tgt::Shader* linkshader = linkShaderProp_.getShader();
    linkshader->activate();

    linkshader->setUniform("projectionMatrix", projectionMatrix);
    linkshader->setUniform("viewMatrix", viewMatrix);
    if(useTexturesProp_.get()) {
        linkshader->setUniform("radius_", radiusProp_.get());
    }
    linkshader->setUniform("lineWidth_", lineWidthProp_.get());

    linkshader->setUniform("colorOrdinary_", linkColorProp_.get());
    linkshader->setUniform("colorPath_", pathLinkColorProp_.get());

    if(useTexturesProp_.get()) {
        linkshader->setUniform("textureZoom_", linkTextureZoomLevel_.get());
        linkTexture_->bind();
        linkTexture_->enable();
    }

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, linkIbo_);
    glDrawElements(GL_LINES, linkCount_, GL_UNSIGNED_SHORT, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    if(useTexturesProp_.get()) {
        linkTexture_->disable();
    }
    linkshader->deactivate();
    LGL_ERROR;
}
void CMMergerTreeRenderer::renderVertIDs(tgt::mat4& projectionMatrix, tgt::mat4& viewMatrix, tgt::Shader* vertshader) {
    vertshader->activate();

    vertshader->setUniform("projectionMatrix", projectionMatrix);
    vertshader->setUniform("viewMatrix", viewMatrix);
    vertshader->setUniform("radius_", radiusProp_.get());

    if(useTexturesProp_.get()) {
        haloTexture_->bind();
        haloTexture_->enable();
    }

    glDrawArrays(GL_POINTS, 0, currentVerts_.size());

    if(useTexturesProp_.get()) {
        haloTexture_->disable();
    }
    vertshader->deactivate();
    LGL_ERROR;
}

void CMMergerTreeRenderer::renderVerts(tgt::mat4& projectionMatrix, tgt::mat4& viewMatrix, tgt::Shader* vertshader) {
    vertshader->activate();

    vertshader->setUniform("projectionMatrix", projectionMatrix);
    vertshader->setUniform("viewMatrix", viewMatrix);
    vertshader->setUniform("radius_", radiusProp_.get());

    vertshader->setUniform("colorOrdinary_", vertColorProp_.get());
    vertshader->setUniform("colorSelected_", selectedVertColorProp_.get());
    vertshader->setUniform("colorMouseOver_", mouseOverVertColorProp_.get());
    vertshader->setUniform("colorOrphan_", newVertColorProp_.get());
    vertshader->setUniform("colorOnPath_", onPathVertColorProp_.get());

    if(useTexturesProp_.get()) {
        haloTexture_->bind();
        haloTexture_->enable();
    }

    glDrawArrays(GL_POINTS, 0, currentVerts_.size());

    if(useTexturesProp_.get()) {
        haloTexture_->disable();
    }
    vertshader->deactivate();
    LGL_ERROR;
}

void CMMergerTreeRenderer::process() {
    outport_.activateTarget();
    outport_.clearTarget();
    const tgt::ivec2 outportSize = outport_.getSize();
    lastOutportWidth_ = outportSize.x;
    if (!vao_ || !vertVbo_ || !linkIbo_) {
        outport_.deactivateTarget();
        return;
    }

    tgt::mat4 projectionMatrix = tgt::Matrix4f::identity;
    projectionMatrix.t00 = 2.0f/windowWidth_.get();
    projectionMatrix.t11 = 2.0f*outportSize.x/(outportSize.y*windowWidth_.get());

    tgt::mat4 viewMatrix = tgt::Matrix4f::identity;
    viewMatrix.t03 = -center_.get().x;
    viewMatrix.t13 = -center_.get().y;
    //viewMatrix.t23 = 0.5f;

    glEnable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glBindVertexArray(vao_);
    //Render Verts first to avoid intel mesa bug?
    renderLinks(projectionMatrix, viewMatrix);
    renderVerts(projectionMatrix, viewMatrix, vertShaderProp_.getShader());
    glBindVertexArray(0);

    glBlendFunc(GL_ONE, GL_ZERO);

    outport_.deactivateTarget();

    //Render into selection frame buffer
    glDisable(GL_BLEND);
    glDisable(GL_DEPTH_TEST);
    selectionManager_->activate();
    selectionManager_->resize(outportSize);
    selectionManager_->clear();
    glBindVertexArray(selectionVao_);
    renderVertIDs(projectionMatrix, viewMatrix, vertSelectionShaderProp_.getShader());
    glBindVertexArray(0);
    selectionManager_->deactivate();

    glEnable(GL_DEPTH_TEST);
    LGL_ERROR;
}

} // namespace
