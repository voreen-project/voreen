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

#include "radarglyphrenderer2d.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

RadarGlyphRenderer2D::RadarGlyphRenderer2D()
    : RadarGlyphRendererBase()
    , sliceIndex_("sliceIndex", "Slice Index:",0,0,10000)
    , boundaryColorProp_("boundaryColorProp","Boundary Color:",tgt::vec4(1.f))
    , mwheelCycleHandler_("mouseWheelHandler", "Slice Cycling", &sliceIndex_)
    , mwheelZoomHandler_("zoomHandler", "Slice Zoom", &zoomProp_, tgt::MouseEvent::CTRL)
    , zoomProp_("zoomProp", "Zoom:",1.f,0.01f,1.f)
    , xOffsetProp_("xOffsetProp", "X Offset:",0.f,-1.f,1.f)
    , yOffsetProp_("yOffsetProp", "Y Offset:",0.f,-1.f,1.f)
{
    //set render modi
    renderModeProp_.addOption("normal_radar_glyphs_2d_xy", "Normal RGs XY-Plane", NORMAL_RADAR_GLYPHS_2D_XY);
    renderModeProp_.addOption("normal_radar_glyphs_2d_xz", "Normal RGs XZ-Plane", NORMAL_RADAR_GLYPHS_2D_XZ);
    renderModeProp_.addOption("normal_radar_glyphs_2d_yz", "Normal RGs YZ-Plane", NORMAL_RADAR_GLYPHS_2D_YZ);

    addProperty(sliceIndex_);
        sliceIndex_.setGroupID("View Settings");
        ON_PROPERTY_CHANGE(sliceIndex_,RadarGlyphRenderer2D,sliceIndexOnChange);
    addProperty(zoomProp_);
        zoomProp_.setGroupID("View Settings");
    addProperty(xOffsetProp_);
        xOffsetProp_.setGroupID("View Settings");
    addProperty(yOffsetProp_);
        yOffsetProp_.setGroupID("View Settings");
    addProperty(boundaryColorProp_);
        boundaryColorProp_.setGroupID("View Settings");
    setPropertyGroupGuiName("View Settings","View Settings");

    mouseEventShift_ = new EventProperty<RadarGlyphRenderer2D>("mouseEvent.Shift", "Slice Shift",
        this, &RadarGlyphRenderer2D::shiftEvent,
        tgt::MouseEvent::MOUSE_BUTTON_LEFT,
        tgt::MouseEvent::PRESSED | tgt::MouseEvent::MOTION, tgt::Event::SHIFT);

    addEventProperty(mouseEventShift_);
    addInteractionHandler(mwheelCycleHandler_);
    addInteractionHandler(mwheelZoomHandler_);
}

RadarGlyphRenderer2D::~RadarGlyphRenderer2D() {
    delete mouseEventShift_;
}

//-----------------------------------------------------------------------------------
//      on change
//-----------------------------------------------------------------------------------
void RadarGlyphRenderer2D::renderModeOnChange(){
    if(vcInport_.isReady()){
        inportOnChange();
    }
}

void RadarGlyphRenderer2D::sliceIndexOnChange(){
    if(vcInport_.isReady()){
        calcGlyphsNew_ = true;
        invalidate();
    }
}

void RadarGlyphRenderer2D::inportOnChange() {
    RadarGlyphRendererBase::inportOnChange();
    const VolumeList* vc = vcInport_.getData();
    tgt::ivec3 dim(vc->first()->getDimensions());
    switch(renderModeProp_.getValue()){
        case NORMAL_RADAR_GLYPHS_2D_YZ:
            sliceIndex_.setMaxValue(dim.x-1);
            x_count_ = vc->first()->getDimensions().z;
            y_count_ = vc->first()->getDimensions().y;
            break;
        case NORMAL_RADAR_GLYPHS_2D_XZ:
            sliceIndex_.setMaxValue(dim.y-1);
            x_count_ = vc->first()->getDimensions().x;
            y_count_ = vc->first()->getDimensions().z;
            break;
        case NORMAL_RADAR_GLYPHS_2D_XY:
            sliceIndex_.setMaxValue(dim.z-1);
            x_count_ = vc->first()->getDimensions().x;
            y_count_ = vc->first()->getDimensions().y;
            break;
    }

    float x_pixel_ = (float)(imgOutport_.getSize().x) / x_count_;
    float y_pixel_ = (float)(imgOutport_.getSize().y) / y_count_;

    if(x_pixel_ <= y_pixel_) {
        x_offset_ = 0.f;
        y_offset_ = ((float)imgOutport_.getSize().y- (x_pixel_ * (float)(y_count_)))/2.f;
        glyphMaxDiameter_ = x_pixel_;
    } else {
        x_offset_ = ((float)imgOutport_.getSize().x- (y_pixel_ * (float)(x_count_)))/2.f;
        y_offset_ = 0.f;
        glyphMaxDiameter_ = y_pixel_;
    }
    numberOfGlyphs_ = static_cast<GLint>(x_count_*y_count_);
}


//-----------------------------------------------------------------------------------
//      process
//-----------------------------------------------------------------------------------
void RadarGlyphRenderer2D::process() {    //TODO volume not 3xfloat?
    //return, if inport is not ready or size == 0
    if (!vcInport_.isReady() || vcInport_.getData()->empty())
        return;

    bool hasChanged = vcInport_.hasChanged();
    const VolumeList* vc = vcInport_.getData();

    if(hasChanged) {
        //return, if all dimensions are not the same
        tgt::svec3 dim = vc->at(0)->getDimensions();
        if(tgt::hmul(dim) == 0) {
            LWARNING("Volume is empty!!!");
            return;
        }
        for(size_t i = 1; i < vc->size(); ++i){
            if(vc->at(i)->getDimensions() != dim){
                LWARNING("Dimensions in the collection are not the same!!!");
                return;
            }
        }
        inportOnChange();
    }

    //set sliceIndex right
    if(hasChanged){
        inportOnChange();
    }

    //test, if program is compiled
    if(!rgProgram_){
        LERROR("No 2D RG Shader!!!");
        return;
    }

    if(calcGlyphsNew_){
        updateGlyphVBO();
    }
    /*if(showIndexProp_.get() && calcIndexNew_){
        updateIndexVBO();
    }*/

    //activate outport
    imgOutport_.activateTarget();

    //bind tftexture
    tgt::TextureUnit tfUnit;
    tgt::Texture* tex = tfProp_.get()->getTexture();
    tfUnit.activate();
    tex->bind();
    LGL_ERROR;

    //get sizes
    size_t vcSize = vc->size();
    size_t volSize = vc->at(0)->getNumVoxels();

    //set projection modelview matrix
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.loadIdentity();
    glOrtho(0.0,imgOutport_.getSize().x,0,imgOutport_.getSize().y,1.0,-1.0);
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.loadIdentity();
    //clear buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //presettings
    MatStack.scale(1.f/zoomProp_.get(),1.f/zoomProp_.get(),1.f/zoomProp_.get());
    MatStack.translate(xOffsetProp_.get()*imgOutport_.getSize().x,yOffsetProp_.get()*imgOutport_.getSize().y,0.f);

    //draw border
    glColor4fv(boundaryColorProp_.get().elem);
    glBegin(GL_LINE_LOOP);
        glVertex3f(x_offset_,y_offset_,0.f);
        glVertex3f(x_offset_,imgOutport_.getSize().y-y_offset_,0.f);
        glVertex3f(imgOutport_.getSize().x-x_offset_,imgOutport_.getSize().y-y_offset_,0.f);
        glVertex3f(imgOutport_.getSize().x-x_offset_,y_offset_,0.f);
    glEnd();
    glColor4f(1.f,1.f,1.f,1.f);

    //render glyphs
    rgProgram_->activate();
    rgProgram_->setUniform("ColorTexture", tfUnit.getUnitNumber());
    rgProgram_->setUniform("iteration_", interpolationStepProp_.get());
    rgProgram_->setUniform("maximum_", maxLengthProp_.get());
    rgProgram_->setUniform("factor_", glyphMaxDiameter_/static_cast<float>(vc->size()*2)); //radius / size
    switch(renderModeProp_.getValue()){
        case NORMAL_RADAR_GLYPHS_2D_XY:
        case NORMAL_RADAR_GLYPHS_2D_YZ:
        case NORMAL_RADAR_GLYPHS_2D_XZ:
            renderRadarGlyphs();
            break;
    }
    rgProgram_->deactivate();

    //render point data
    /*if(!pointData_.empty()){
        pdProgram_->activate();
        pdProgram_->setUniform("ColorTexture", tfUnit.getUnitNumber());
        renderPointData();
        pdProgram_->deactivate();
    }
    pointData_.clear();*/


    //render active index
    if(showIndexProp_.get()){
        bbProgram_->activate();
        bbProgram_->setUniform("radius",glyphMaxDiameter_/16.f * sphereSizeProp_.get()/zoomProp_.get());
        bbProgram_->setUniform("factor_",(activeIndexProp_.get()+1)*glyphMaxDiameter_/(static_cast<float>(vc->size())*2.f));
        glDisable(GL_DEPTH_TEST);// otherwise the index will not be rendered BUG?
        renderActiveIndex();
        glEnable(GL_DEPTH_TEST);
        bbProgram_->deactivate();
     }

    //reset projection modelview matrix
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);

    //deactivate outport
    imgOutport_.deactivateTarget();
}

//-----------------------------------------------------------------------------------
//      fill VBOs
//-----------------------------------------------------------------------------------
void RadarGlyphRenderer2D::fillGlyphVBO() {
    const VolumeList* vc = vcInport_.getData();

    GLint size = sizeof(GLfloat)*7;
    GLfloat* vertexElement = new GLfloat[7](); //dir(x,y,z), center(x,y,z), step
    GLint offset = 0;

    for(size_t i = 0; i < y_count_; i++){
        for(size_t j = 0; j < x_count_; j++){
            setProgress(static_cast<float>(i*x_count_+j)/(y_count_*x_count_-1));
            glBindBuffer(GL_ARRAY_BUFFER, glyphVBO_); //should not be necessary, but otherwise setProgress is buggy
            tgt::vec3 center(x_offset_+(j+0.5f)*glyphMaxDiameter_,y_offset_+(i+0.5f)*glyphMaxDiameter_,0.f);
            vertexElement[0] = 0.f;
            vertexElement[1] = 0.f;
            vertexElement[2] = 0.f;
            vertexElement[3] = center.x;
            vertexElement[4] = center.y;
            vertexElement[5] = center.z;
            vertexElement[6] = 0.f;
            glBufferSubData(GL_ARRAY_BUFFER, offset, size, vertexElement);
            offset += size;
            for(size_t k = 0; k < vc->size(); k=k+volumeStepProp_.get()){
                tgt::vec3 dis= vc->at(k)->getRepresentation<VolumeRAM_3xFloat>()->voxel(permute(sliceIndex_.get(),i,j));
                switch(renderModeProp_.getValue()){
                case NORMAL_RADAR_GLYPHS_2D_YZ:
                    dis.x = dis.z;
                    dis.z = 0.f;
                    break;
                case NORMAL_RADAR_GLYPHS_2D_XZ:
                    dis.y = dis.z;
                    dis.z = 0.f;
                    break;
                case NORMAL_RADAR_GLYPHS_2D_XY:
                    dis.z = 0.f;
                    break;
                }
                vertexElement[0] = dis.x;
                vertexElement[1] = dis.y;
                vertexElement[2] = dis.z;
                vertexElement[3] = center.x;
                vertexElement[4] = center.y;
                vertexElement[5] = center.z;
                vertexElement[6] = static_cast<GLfloat>(k+1);
                glBufferSubData(GL_ARRAY_BUFFER, offset, size, vertexElement);
                offset += size;
            }
        }
    }
    delete[] vertexElement;
}

/*void RadarGlyphRenderer2D::fillIndexVBO() {
    const VolumeCollection* vc = vcInport_.getData();

    GLint size = sizeof(GLfloat)*3;
    GLfloat* vertexElement = new GLfloat[3](); //point
    GLint offset = 0;

    int k = activeIndexProp_.get();
    float factor = (k+1)*glyphMaxDiameter_/(static_cast<float>(vc->size())*2.f);

     for(int i = 0; i < y_count_; i++){
        for(int j = 0; j < x_count_; j++){
            tgt::vec3 center(x_offset_+(j+0.5f)*glyphMaxDiameter_,y_offset_+(i+0.5f)*glyphMaxDiameter_,0.f);
            tgt::vec3 dis= vc->at(k)->getRepresentation<VolumeRAM_3xFloat>()->voxel(permute(sliceIndex_.get(),i,j));
            switch(renderModeProp_.getValue()){
            case NORMAL_RADAR_GLYPHS_2D_YZ:
                dis.x = dis.z;
                dis.z = 0.f;
                break;
            case NORMAL_RADAR_GLYPHS_2D_XZ:
                dis.y = dis.z;
                dis.z = 0.f;
                break;
            case NORMAL_RADAR_GLYPHS_2D_XY:
                dis.z = 0.f;
                break;
            }
            if(tgt::length(dis) > 0.f)
                dis = tgt::normalize(dis);
            tgt::vec3 point = center + (dis*factor);
            vertexElement[0] = point.x;
            vertexElement[1] = point.y;
            vertexElement[2] = point.z;
            glBufferSubData(GL_ARRAY_BUFFER, offset, size, vertexElement);
            offset += size;
        }
    }
    delete[] vertexElement;
}*/

tgt::ivec3 RadarGlyphRenderer2D::permute(int slice, int y, int x){
    switch(renderModeProp_.getValue()){
    case NORMAL_RADAR_GLYPHS_2D_YZ:
        return tgt::ivec3(slice,y,x);
        break;
    case NORMAL_RADAR_GLYPHS_2D_XZ:
        return tgt::ivec3(x,slice,y);
        break;
    case NORMAL_RADAR_GLYPHS_2D_XY:
        return tgt::ivec3(x,y,slice);
        break;
    }
    tgtAssert(false,"Shouldn't get here!!!");
    return tgt::ivec3(0);
}

//-----------------------------------------------------------------------------------
//      events
//-----------------------------------------------------------------------------------
void RadarGlyphRenderer2D::shiftEvent(tgt::MouseEvent* e) {
    e->ignore();
    if (!imgOutport_.isReady())
        return;

    if (e->action() == tgt::MouseEvent::PRESSED) {
        mousePosition_ = e->coord();
        return;
    }

    tgt::vec2 mouseOffset = tgt::vec2(e->coord() - mousePosition_)/tgt::vec2(imgOutport_.getSize());

    xOffsetProp_.set(tgt::clamp(xOffsetProp_.get()+mouseOffset.x*zoomProp_.get(),-1.f,1.f));
    yOffsetProp_.set(tgt::clamp(yOffsetProp_.get()-mouseOffset.y*zoomProp_.get(),-1.f,1.f));
    mousePosition_ = e->coord();

    e->accept();
}

}   // namespace
