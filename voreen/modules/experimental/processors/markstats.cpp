/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "markstats.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/utils/stringutils.h"

#include "tgt/textureunit.h"
#include "tgt/immediatemode/immediatemode.h"
#include "tgt/event/mouseevent.h"

namespace voreen {

const std::string MarkStats::loggerCat_("voreen.MarkStats");

MarkStats::MarkStats()
    : ImageProcessor("image/background")
    , inport_            (Port::INPORT,  "inport",            "Slice Inport")
    , renderOutport_     (Port::OUTPORT, "renderOutport",     "Volume Outport")
    , textOutport_       (Port::OUTPORT, "textOutport",       "Text Outport")
    , geometryOutport_   (Port::OUTPORT, "geometryOutport",   "Geometry Outport")
    , segmentationId_("segementationId", "Segmentation ID", 1, 0, std::numeric_limits<uint16_t>::max())
    , brushAlpha_("brushAlpha", "Brush alpha", 1.0f, 0.0f, 1.0f)
    , linkedMousePositionInSlice_("linkedMousePositionInSlice","\"Mouse Position\" from SliceViewer",tgt::ivec3(-1),tgt::ivec3(-1),tgt::ivec3(INT_MAX),
                                  Processor::INVALID_RESULT,NumericProperty<tgt::ivec3>::STATIC,Property::LOD_DEBUG)
    , linkedSliceAlignment_("linkedSliceAlignment", "\"Slice Alignment\" from SliceViewer",Processor::INVALID_RESULT,false,Property::LOD_DEBUG)
    , linkedSliceIndex_("linkedSliceIndex", "\"Slice Number\" from SliceViewer", 0, 0, 10000,Processor::INVALID_RESULT,NumericProperty<int>::STATIC, Property::LOD_DEBUG)
    , linkedPickingMatrix_("linkedPickingMatrix", "\"Picking Matrix\" from SliceViewer", tgt::mat4::createIdentity(), tgt::mat4(-1e6f), tgt::mat4(1e6f),
                            Processor::INVALID_RESULT, NumericProperty<tgt::mat4>::STATIC,Property::LOD_DEBUG)
{
    addPort(inport_);

    addPort(renderOutport_);
    addPort(textOutport_);
    addPort(geometryOutport_);

    addProperty(segmentationId_);

    addProperty(brushAlpha_);
    brushAlpha_.setGroupID("brush");
    setPropertyGroupGuiName("brush", "Brush");

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

    mouseEventPress_ = new EventProperty<MarkStats>("mouseEvent.cursorPositionPress", "Cursor Position Press",
        this, &MarkStats::mouseLocalization,
        static_cast<tgt::MouseEvent::MouseButtons>(tgt::MouseEvent::MOUSE_BUTTON_LEFT),
        tgt::MouseEvent::PRESSED | tgt::MouseEvent::RELEASED, tgt::Event::MODIFIER_NONE, true);


    addEventProperty(mouseEventPress_);




    colors_.push_back(tgt::vec4(1.0f, 0.0f, 0.0f, 1.0f));
    colors_.push_back(tgt::vec4(0.0f, 1.0f, 0.0f, 1.0f));
    colors_.push_back(tgt::vec4(0.0f, 0.0f, 1.0f, 1.0f));
    colors_.push_back(tgt::vec4(1.0f, 1.0f, 0.0f, 1.0f));
    colors_.push_back(tgt::vec4(1.0f, 0.0f, 1.0f, 1.0f));
    colors_.push_back(tgt::vec4(0.0f, 1.0f, 1.0f, 1.0f));
    colors_.push_back(tgt::vec4(1.0f, 0.5f, 0.5f, 1.0f));
    colors_.push_back(tgt::vec4(0.5f, 1.0f, 0.5f, 1.0f));
    colors_.push_back(tgt::vec4(0.5f, 0.5f, 1.0f, 1.0f));
}

MarkStats::~MarkStats() {
    delete mouseEventPress_;
}

bool MarkStats::isReady() const{
    return inport_.isReady();
}


void MarkStats::process(){
    renderOutport_.activateTarget();
    glBindTexture(GL_TEXTURE_2D, inport_.getRenderTarget()->getColorTexture()->getId());
    renderQuad();

    setRenderState();

    IMode.begin(tgt::ImmediateMode::QUADS);
    for(size_t i = 0; i != positions_.size(); i++){
        MarkPoint p = positions_[i];
        IMode.color(colors_[(p.id_-1)%colors_.size()]);
        IMode.vertex(tgt::vec3(p.position_)+tgt::vec3(-0.5f,  0.5f, 0));
        IMode.vertex(tgt::vec3(p.position_)+tgt::vec3(-0.5f, -0.5f, 0));
        IMode.vertex(tgt::vec3(p.position_)+tgt::vec3( 0.5f, -0.5f, 0));
        IMode.vertex(tgt::vec3(p.position_)+tgt::vec3( 0.5f,  0.5f, 0));
    }

    IMode.end();

    drawBrush();


    resetRenderState();
    if (shouldSetOutport_){
        shouldSetOutport_ = false;
    }

    renderOutport_.deactivateTarget();
}

void MarkStats::setDescriptions(){
    setDescription("TODO.");
}

void MarkStats::setRenderState(){
    MatStack.pushMatrix();
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    tgt::mat4 pickingMatrix = linkedPickingMatrix_.get();
    tgt::mat4 invPickingMatrix;
    pickingMatrix.invert(invPickingMatrix);

    MatStack.multMatrix(tgt::mat4::createOrtho(0.0f, 1.0f*renderOutport_.getSize().x, 0.0f, 1.0f*renderOutport_.getSize().y, -100.0f, 100.0f));
    MatStack.multMatrix(invPickingMatrix);

    tgt::vec4 color;
    int id = segmentationId_.get();
    if (id == 0){
        color = tgt::vec4(1.0);
    }else{
        color = colors_[(id - 1) % colors_.size()];
    }
    color.a = brushAlpha_.get();

    IMode.color(color);
}

void MarkStats::resetRenderState(){
    IMode.color(tgt::vec4::one);
    glBlendFunc(GL_ONE, GL_ZERO);
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    MatStack.popMatrix();
}

void MarkStats::drawBrush(){

    tgt::ivec3 basex;
    tgt::ivec3 basey;
    tgt::ivec3 mouseVolume = linkedMousePositionInSlice_.get();

    createSliceBaseVectors(basex, basey);

    IMode.begin(tgt::ImmediateMode::QUADS);
        IMode.vertex(tgt::vec3(mouseVolume)+tgt::vec3(-0.5f,  0.5f, 0));
        IMode.vertex(tgt::vec3(mouseVolume)+tgt::vec3(-0.5f, -0.5f, 0));
        IMode.vertex(tgt::vec3(mouseVolume)+tgt::vec3( 0.5f, -0.5f, 0));
        IMode.vertex(tgt::vec3(mouseVolume)+tgt::vec3( 0.5f,  0.5f, 0));
    IMode.end();
}


void MarkStats::createSliceBaseVectors(tgt::ivec3 &basex, tgt::ivec3 &basey){
    switch (linkedSliceAlignment_.getValue())
    {
    case XY_PLANE:
        basex = tgt::ivec3(1, 0, 0);
        basey = tgt::ivec3(0, 1, 0);
        break;
    case XZ_PLANE:
        basex = tgt::ivec3(1, 0, 0);
        basey = tgt::ivec3(0, 0, 1);
        break;
    case YZ_PLANE:
        basex = tgt::ivec3(0, 1, 0);
        basey = tgt::ivec3(0, 0, 1);
        break;
    default:
        LERROR("Unknown slice alignment!");
        return;
    }
}

void MarkStats::mouseLocalization(tgt::MouseEvent* ev){
    if (ev->action() == tgt::MouseEvent::PRESSED && ev->button() == tgt::MouseEvent::MOUSE_BUTTON_LEFT){
        /*isPainting_ = true;
        lastMousePos_ = linkedMousePositionInSlice_.get();
        currentLine_.clear();
        currentLine_.push_back(lastMousePos_);*/
    }else if(ev->action() == tgt::MouseEvent::RELEASED){
        MarkPoint p;
        p.position_ = linkedMousePositionInSlice_.get();
        p.id_ = segmentationId_.get();
        positions_.push_back(p);
        shouldSetOutport_ = true;
    }
}

void MarkStats::pickSegmentationId(tgt::MouseEvent* ev){
    /*tgt::ivec3 pos = linkedMousePositionInSlice_.get();
    tgt::ivec3 dim = volumeInport_.getData()->getDimensions();

    if (tgt::hor(tgt::lessThan(pos, tgt::ivec3::zero)) || tgt::hor(tgt::greaterThanEqual(pos, dim))){
        return;
    }

    segmentationId_.set(volumeram_->voxel(pos));*/
}


} // namespace voreen
