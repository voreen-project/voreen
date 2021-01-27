/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "manualsegmentation.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/utils/stringutils.h"

#include "tgt/textureunit.h"
#include "tgt/immediatemode/immediatemode.h"
#include "tgt/event/mouseevent.h"

namespace voreen {

const std::string ManualSegmentation::loggerCat_("voreen.ManualSegmentation");

ManualSegmentation::ManualSegmentation()
    : ImageProcessor("image/background")
    , inport_                (Port::INPORT,  "inport",            "Slice Inport")
    , renderOutport_         (Port::OUTPORT, "renderOutport",     "Slice Outport")
    , storageCoprocessorPort_(Port::INPORT,  "coprocessorPort",   "Coprocessor Port")
    , segmentationId_("segementationId", "Segmentation ID", 1, 0, 999)
    , brushRadius_("brushSize",     "Brush size",  1, 1, 30)
    , brushAlpha_("brushAlpha",   "Brush alpha", 1.0f, 0.0f, 1.0f)
    , delayedDraw_("delayedDraw", "Draw delayed", true)
    , usingEraser_("usingEraser", "Use Eraser", false)
    , doubleclickLatancy_("doubleclickLatancy", "Latancy for double click(in ms)", 300, 0, 1000)
    , linkedMousePositionInSlice_("linkedMousePositionInSlice","\"Mouse Position\" from SliceViewer",tgt::ivec3(-1),tgt::ivec3(-1),tgt::ivec3(INT_MAX),
                                  Processor::INVALID_RESULT,NumericProperty<tgt::ivec3>::STATIC,Property::LOD_DEBUG)
    , linkedSliceAlignment_("linkedSliceAlignment", "\"Slice Alignment\" from SliceViewer",Processor::INVALID_RESULT,false,Property::LOD_DEBUG)
    , linkedSliceIndex_("linkedSliceIndex", "\"Slice Number\" from SliceViewer", 0, 0, 10000,Processor::INVALID_RESULT,NumericProperty<int>::STATIC, Property::LOD_DEBUG)
    , linkedPickingMatrix_("linkedPickingMatrix", "\"Picking Matrix\" from SliceViewer", tgt::mat4::createIdentity(), tgt::mat4(-1e6f), tgt::mat4(1e6f),
                            Processor::INVALID_RESULT, NumericProperty<tgt::mat4>::STATIC,Property::LOD_DEBUG)
    , isPainting_(false)
    , justFinishedPainting_(false)
    , shouldCreateBrush_(true)
{

    addPort(inport_);
    addPort(renderOutport_);
    addPort(storageCoprocessorPort_);


    addProperty(segmentationId_);

    addProperty(brushRadius_);
    brushRadius_.setGroupID("brush");
    addProperty(brushAlpha_);
    brushAlpha_.setGroupID("brush");
    addProperty(delayedDraw_);
    delayedDraw_.setGroupID("brush");
    addProperty(usingEraser_);
    usingEraser_.setGroupID("brush");
    setPropertyGroupGuiName("brush", "Brush");


    addProperty(doubleclickLatancy_);
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

    tgt::Event::Modifier ALL = static_cast<tgt::Event::Modifier>(tgt::Event::SHIFT | tgt::Event::CTRL | tgt::Event::ALT |
        tgt::Event::META | tgt::Event::NUM |tgt::Event::CAPS |tgt::Event::MODE);

    mouseEventPress_ = new EventProperty<ManualSegmentation>("mouseEvent.cursorPositionPress", "Cursor Position Press",
        this, &ManualSegmentation::mouseLocalization,
        static_cast<tgt::MouseEvent::MouseButtons>(tgt::MouseEvent::MOUSE_BUTTON_LEFT),
        tgt::MouseEvent::PRESSED | tgt::MouseEvent::RELEASED, ALL, true);

    keyEventToggleErase_ = new EventProperty<ManualSegmentation>("keyEvent.toggleEraser", "Toggle Eraser",
        this, &ManualSegmentation::toggleEraser, tgt::KeyEvent::K_E);

    mouseEventPick_ = new EventProperty<ManualSegmentation>("mouseEvent.picking", "Cursor Position Press",
        this, &ManualSegmentation::pickSegmentationId,
        static_cast<tgt::MouseEvent::MouseButtons>(tgt::MouseEvent::MOUSE_BUTTON_LEFT),
        tgt::MouseEvent::PRESSED | tgt::MouseEvent::RELEASED, tgt::Event::SHIFT, true);


    mouseEventStopPainting_ = new EventProperty<ManualSegmentation>("mouseEvent.stoppainting", "Cursor Position Press",
        this, &ManualSegmentation::stopPainting,
        static_cast<tgt::MouseEvent::MouseButtons>(tgt::MouseEvent::MOUSE_BUTTON_RIGHT),
        tgt::MouseEvent::PRESSED | tgt::MouseEvent::RELEASED, tgt::Event::MODIFIER_NONE, true);

    keyEventBlockPainting_ = new EventProperty<ManualSegmentation>("keyEvent.blockPainting", "Block Painting",
        this, &ManualSegmentation::blockPainting, tgt::KeyEvent::K_W);

    addEventProperty(mouseEventPress_);
    addEventProperty(keyEventToggleErase_);
    addEventProperty(mouseEventStopPainting_);
    addEventProperty(mouseEventPick_);
    addEventProperty(keyEventBlockPainting_);

    ON_CHANGE(brushRadius_, ManualSegmentation, invalidateBrush);
    ON_CHANGE(usingEraser_, ManualSegmentation, usingEraserChanged);

    blockPainting_ = false;
}

ManualSegmentation::~ManualSegmentation() {
    delete mouseEventPress_;
    delete mouseEventPick_;
    delete mouseEventStopPainting_;
    delete keyEventToggleErase_;
    delete keyEventBlockPainting_;
}


void ManualSegmentation::initialize() {
    ImageProcessor::initialize();

}

void ManualSegmentation::process(){
    renderOutport_.activateTarget();
    glBindTexture(GL_TEXTURE_2D, inport_.getRenderTarget()->getColorTexture()->getId());
    renderQuad();

    if(shouldCreateBrush_){
        createBrush();
        shouldCreateBrush_ = false;
    }
    setRenderState();
    drawBrush();
    if (isPainting_){
        tgt::ivec3 linkedMousePos = linkedMousePositionInSlice_.get();
        if (linkedMousePos != tgt::ivec3(-1)){
            // Mouse inside of slice
            currentLines_.back().line_.push_back(linkedMousePositionInSlice_.get());
        }else{
            // start a new line
            Line l;
            l.isErasing_ = usingEraser_.get();
            currentLines_.push_back(l);
        }
        if (delayedDraw_.get())
            drawLines();
    }
    resetRenderState();

    if (justFinishedPainting_ || (!delayedDraw_.get())){
        applyBrush();
        if (isPainting_ && !delayedDraw_.get()){
            tgt::ivec3 linkedMousePos = linkedMousePositionInSlice_.get();
            currentLines_.clear();
            Line l;
            l.isErasing_ = usingEraser_.get();
            if (linkedMousePos != tgt::ivec3(-1)){
                l.line_.push_back(linkedMousePos);
            }
            currentLines_.push_back(l);

            storageCoprocessorPort_.getConnectedProcessor()->invalidateSegmentation();
        }
        else if (justFinishedPainting_)
            storageCoprocessorPort_.getConnectedProcessor()->invalidateSegmentation();

        justFinishedPainting_ = false;
    }


    renderOutport_.deactivateTarget();
}

void ManualSegmentation::setDescriptions(){
    setDescription(
        "This processor allows to paint a manual segmentation.\n"
        "It is used as a overlay for a SliceViewer and need mouseport, pickingmatrix, sliceindex and slicealignment linked from it<br/>"
        "A ManualSegmentationStorage need to be connected on the side with the volume that should be segmented.<br/>"
        "While the mouse is hold down you can paint the current id.<br/>"
        "A eraser can be enabled and disable by pressing the 'E' key.<br/>"
        "Shift+Click for taking the segmentation from the current position and using it as the current id"
        );

    segmentationId_.setDescription("The id of the segmented object used for the next painting operation.");
    brushRadius_.setDescription("The radius of the brush for painting.");
    brushAlpha_.setDescription("This is the alpha of the brush and does not affect the segmentation.");
    delayedDraw_.setDescription(
        "If this is active the output is only updated after a line stroke and not while painting.<br/>"
        "While painting only a small line representing the stroke is shown.<br/>"
        "This is much faster and useful when not to much precision is needed.");

    inport_.setDescription("A port a sliceviewer as input.");
    renderOutport_.setDescription("Sliceviewer from input with the brushs overlayed.");
    storageCoprocessorPort_.setDescription("The ManualSegmentationStorage processor that stores the segmentation worked on.");
    usingEraser_.setDescription("Selects if the eraser is used at the moment. Can be toggle by pressing the 'E' button.");
}

void ManualSegmentation::setRenderState(){
    MatStack.pushMatrix();
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    tgt::mat4 pickingMatrix = linkedPickingMatrix_.get();
    tgt::mat4 invPickingMatrix;
    pickingMatrix.invert(invPickingMatrix);

    MatStack.multMatrix(tgt::mat4::createOrtho(0.0f, 1.0f*renderOutport_.getSize().x, 0.0f, 1.0f*renderOutport_.getSize().y, -10000.0f, 10000.0f));
    MatStack.multMatrix(invPickingMatrix);


}

void ManualSegmentation::resetRenderState(){
    IMode.color(tgt::vec4::one);
    glBlendFunc(GL_ONE, GL_ZERO);
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    MatStack.popMatrix();
}

void ManualSegmentation::drawBrush(){
    int brushsize = brushRadius_.get();
    int realbrushsize = brushsize*2+1;

    tgt::ivec3 basex;
    tgt::ivec3 basey;
    tgt::ivec3 mouseVolume = linkedMousePositionInSlice_.get();
    if (mouseVolume == tgt::ivec3(-1)){
        return;
    }

    createSliceBaseVectors(basex, basey);
    IMode.color(getColorForId(segmentationId_.get(), usingEraser_.get()));

    IMode.begin(tgt::ImmediateMode::QUADS);
    for(int y = -brushsize; y <= brushsize; y++){
        for(int x = -brushsize; x <= brushsize; x++){
            tgt::ivec3 pos = mouseVolume+y*basey+x*basex;
            if (brush_[x+brushsize+realbrushsize*(y+brushsize)] == 0){
                continue;
            }
            IMode.vertex(tgt::vec3(pos)-0.5f*tgt::vec3(basex)+0.5f*tgt::vec3(basey));//+tgt::vec3(-0.5f,  0.5f, 0));
            IMode.vertex(tgt::vec3(pos)-0.5f*tgt::vec3(basex)-0.5f*tgt::vec3(basey));//+tgt::vec3(-0.5f, -0.5f, 0));
            IMode.vertex(tgt::vec3(pos)+0.5f*tgt::vec3(basex)-0.5f*tgt::vec3(basey));//+tgt::vec3( 0.5f, -0.5f, 0));
            IMode.vertex(tgt::vec3(pos)+0.5f*tgt::vec3(basex)+0.5f*tgt::vec3(basey));//+tgt::vec3( 0.5f,  0.5f, 0));
        }
    }
    IMode.end();
}

void ManualSegmentation::drawLines(){
    for(auto line: currentLines_){
        IMode.begin(tgt::ImmediateMode::LINE_STRIP);
        IMode.color(getColorForId(segmentationId_.get(), line.isErasing_));
        for(auto point: line.line_){
            IMode.vertex(point);
        }
        IMode.end();
    }
}

void ManualSegmentation::createSliceBaseVectors(tgt::ivec3 &basex, tgt::ivec3 &basey){
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

void ManualSegmentation::mouseLocalization(tgt::MouseEvent* ev){
    //if (ev->action() == tgt::MouseEvent::CLICK)
    if (ev->action() == tgt::MouseEvent::PRESSED && ev->button() == tgt::MouseEvent::MOUSE_BUTTON_LEFT && ev->modifiers() == tgt::Event::MODIFIER_NONE && !blockPainting_){
        isPainting_ = true;
        Line l;
        currentLines_.clear();
        l.isErasing_ = usingEraser_.get();
        currentLines_.push_back(l);
        currentLines_.back().line_.push_back(linkedMousePositionInSlice_.get());
        invalidate();
    }else if(ev->action() == tgt::MouseEvent::RELEASED){
        isPainting_ = false;
        justFinishedPainting_ = true;
        invalidate();

    }
}



void ManualSegmentation::toggleEraser(tgt::KeyEvent* ev){
    if (!ev->pressed())
        usingEraser_.set(!usingEraser_.get());
    ev->accept();
}


void ManualSegmentation::pickSegmentationId(tgt::MouseEvent* ev){
    VolumeRAM_UInt16 * segVolume = storageCoprocessorPort_.getConnectedProcessor()->getSegmentationVolume();
    tgt::ivec3 pos = linkedMousePositionInSlice_.get();
    tgt::ivec3 dim = segVolume->getDimensions();

    if (tgt::hor(tgt::lessThan(pos, tgt::ivec3::zero)) || tgt::hor(tgt::greaterThanEqual(pos, dim))){
        return;
    }

    segmentationId_.set(segVolume->voxel(pos));
}

void ManualSegmentation::stopPainting(tgt::MouseEvent* ev){
    isPainting_ = false;
    currentLines_.clear();
}


void ManualSegmentation::blockPainting(tgt::KeyEvent* ev)
{
    if (ev->pressed()){
        blockPainting_ = true;
        if (isPainting_){
            justFinishedPainting_ = true;
            isPainting_ = false;
            invalidate();
        }
    }else{
        blockPainting_ = false;
    }
}



void ManualSegmentation::invalidateBrush(){
    shouldCreateBrush_ = true;
}


void ManualSegmentation::createBrush(){
    int brushsize = brushRadius_.get();
    int realBrushSize = 2*brushsize+1;
    brush_.resize(realBrushSize*realBrushSize);

    for(int y = 0; y != realBrushSize; y++){
        for(int x = 0; x != realBrushSize; x++){
            int dx = x-brushsize;
            int dy = y-brushsize;
            int d = dx*dx+dy*dy;
            brush_[realBrushSize*y+x] = (d < brushsize*brushsize ? 0xff : 0);
        }
    }
}

void ManualSegmentation::applyBrush(){

    if (currentLines_.empty())
        return;

    for(auto line: currentLines_){
        uint16_t id = line.isErasing_ ? 0: segmentationId_.get();
        storageCoprocessorPort_.getConnectedProcessor()->updateHighestId(id);
        if (line.line_.size() >= 2){
            for(size_t i = 0; i < line.line_.size()-1; i++){
                fillLine(line.line_[i], line.line_[i+1], id);
            }
        }else if (line.line_.size() == 1){
            fillLine(line.line_[0], line.line_[0], id);
        }
    }
    currentLines_.clear();
    storageCoprocessorPort_.getConnectedProcessor()->invalidateSegmentation();
}

void ManualSegmentation::fillLine(tgt::vec3 start, tgt::vec3 end, uint16_t id){
    VolumeRAM_UInt16 * segVolume = storageCoprocessorPort_.getConnectedProcessor()->getSegmentationVolume();
    if(!segVolume)
        return;
    float len = std::max(1.0f, tgt::distance(start, end));
    tgt::ivec3 dim = tgt::ivec3(segVolume->getDimensions());
    int brushsize = brushRadius_.get();
    int realbrushsize = brushsize*2+1;
    //tgt::ivec3 mouseVolume = linkedMousePositionInSlice_.get();

    tgt::ivec3 basex;
    tgt::ivec3 basey;

    createSliceBaseVectors(basex, basey);


    for(float t = 0; t <= 1; t+=0.9f/len){
        tgt::ivec3 p = tgt::ivec3(tgt::round(tgt::mix(start, end, t)));
        for(int y = -brushsize; y <= brushsize; y++){
            for(int x = -brushsize; x <= brushsize; x++){
                tgt::ivec3 pos = p+y*basey+x*basex;
                if (tgt::hor(tgt::lessThan(pos, tgt::ivec3::zero)) || tgt::hor(tgt::greaterThanEqual(pos, dim))){
                    continue;
                }
                if (brush_[x+brushsize+realbrushsize*(y+brushsize)] == 0){
                    continue;
                }
                segVolume->voxel(pos) = id;
            }
        }
    }
}

tgt::vec4 ManualSegmentation::getColorForId(int id, bool usingEraser)
{
    if (usingEraser)
        id = 0;
    tgt::vec4 color = storageCoprocessorPort_.getConnectedProcessor()->getColorForId(id);
    color.a = brushAlpha_.get();

    return color;
}

void ManualSegmentation::usingEraserChanged()
{
    if (isPainting_){
        tgt::ivec3 mousePos = linkedMousePositionInSlice_.get();
        // finish current line
        currentLines_.back().line_.push_back(mousePos);

        // start a new line with the new eraser setting
        Line l;
        l.isErasing_ = usingEraser_.get();
        currentLines_.push_back(l);
        currentLines_.back().line_.push_back(mousePos);
    }
}

} // namespace voreen
