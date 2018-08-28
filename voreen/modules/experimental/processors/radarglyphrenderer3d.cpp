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

#include "radarglyphrenderer3d.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/interaction/camerainteractionhandler.h"
#include "tgt/glmath.h"

namespace voreen {

RadarGlyphRenderer3D::RadarGlyphRenderer3D()
    : RadarGlyphRendererBase()
    , cameraProp_("cameraProp","Camera", tgt::Camera(tgt::vec3(0.f, 0.f, 3.5f), tgt::vec3(0.f, 0.f, 0.f), tgt::vec3(0.f, 1.f, 0.f)))
{
    //set modi
    renderModeProp_.addOption("normal_radar_glyphs_3d", "Normal RGs", NORMAL_RADAR_GLYPHS_3D);

    addProperty(cameraProp_);
        cameraProp_.setGroupID("View Settings");
    setPropertyGroupGuiName("View Settings","View Settings");

    cameraHandler_ = new CameraInteractionHandler("cameraHandler", "Camera", &cameraProp_);
    addInteractionHandler(cameraHandler_);
}

RadarGlyphRenderer3D::~RadarGlyphRenderer3D() {
    delete cameraHandler_;
}

//-----------------------------------------------------------------------------------
//      on change
//-----------------------------------------------------------------------------------
void RadarGlyphRenderer3D::renderModeOnChange(){
}

void RadarGlyphRenderer3D::inportOnChange(){
    RadarGlyphRendererBase::inportOnChange();
    cameraProp_.adaptInteractionToScene(vcInport_.getData()->at(0)->getBoundingBox().getBoundingBox());
}

//-----------------------------------------------------------------------------------
//      process
//-----------------------------------------------------------------------------------
void RadarGlyphRenderer3D::process() {    //TODO volume not 3xfloat?
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

    //test, if program is compiled
    if (!rgProgram_){
        LERROR("No 3D RG Shader!!!");
        return;
    }

    if(calcGlyphsNew_){
        updateGlyphVBO();
    }
    /*if(showIndexProp_.get() && calcIndexNew_){
        updateIndexVBO();
    }*/

    //bind tftexture
    tgt::TextureUnit tfUnit;
    tgt::Texture* tex = tfProp_.get()->getTexture();
    tfUnit.activate();
    tex->bind();

    //activate outport
    imgOutport_.activateTarget();

    //get sizes
    size_t vcSize = vc->size();
    size_t volSize = vc->at(0)->getNumVoxels();
    float glyphMaxRadius = tgt::min(vc->at(0)->getSpacing())/2.f;

    //set projection modelview matrix
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.loadMatrix(cameraProp_.get().getProjectionMatrix(imgOutport_.getSize()));
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.loadMatrix(cameraProp_.get().getViewMatrix());
    LGL_ERROR;
    //clear buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //render glyphs
    rgProgram_->activate();
    rgProgram_->setUniform("ColorTexture", tfUnit.getUnitNumber());
    rgProgram_->setUniform("iteration_", interpolationStepProp_.get());
    rgProgram_->setUniform("maximum_", maxLengthProp_.get());
    rgProgram_->setUniform("factor_", glyphMaxRadius/static_cast<float>(vc->size())); //radius / size
    switch(renderModeProp_.getValue()){
    case NORMAL_RADAR_GLYPHS_3D:
        renderRadarGlyphs();;
        break;
    }
    rgProgram_->deactivate();

    tfUnit.setZeroUnit();

    //render point data
    /*if(!pointData_.empty()){
        pdProgram_->activate();
        pdProgram_->setUniform("ColorTexture", tfUnit.getUnitNumber());
        renderPointData();
        pdProgram_->deactivate();
        pointData_.clear();
    }*/

    //render bill boards
    if(showIndexProp_.get()) {
        bbProgram_->activate();
        bbProgram_->setUniform("radius",glyphMaxRadius/8.f * sphereSizeProp_.get());
        bbProgram_->setUniform("factor_",(activeIndexProp_.get()+1)*glyphMaxRadius/(static_cast<float>(vc->size())));
        renderActiveIndex();
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
//      calculate glyphs
//-----------------------------------------------------------------------------------
void RadarGlyphRenderer3D::fillGlyphVBO() {
    const VolumeList* vc = vcInport_.getData();
    tgt::svec3 dim = vc->at(0)->getDimensions();
    tgt::vec3 spacing = vc->at(0)->getSpacing()/2.f;
    tgt::mat4 v2w = vc->at(0)->getVoxelToWorldMatrix();

    GLint size = sizeof(GLfloat)*7;
    GLfloat* vertexElement = new GLfloat[7](); //dir(x,y,z), center(x,y,z), step
    GLint offset = 0;

    for (size_t i = 0; i < dim.x ; ++i) {
        for (size_t j = 0; j < dim.y; ++j) {
            for (size_t k = 0; k < dim.z; ++k){
                LGL_ERROR;
                setProgress(static_cast<float>(i*dim.y*dim.z+j*dim.z+k)/(dim.x*dim.y*dim.z-1));
                glBindBuffer(GL_ARRAY_BUFFER, glyphVBO_); //should not be necessary, but otherwise setProgress is buggy
                tgt::vec4 center = v2w*(tgt::vec4(static_cast<float>(i),static_cast<float>(j),static_cast<float>(k),1.f))+tgt::vec4(spacing,0.f);
                LGL_ERROR;
                vertexElement[0] = 0.f;
                vertexElement[1] = 0.f;
                vertexElement[2] = 0.f;
                vertexElement[3] = center.x;
                vertexElement[4] = center.y;
                vertexElement[5] = center.z;
                vertexElement[6] = 0.f;
                glBufferSubData(GL_ARRAY_BUFFER, offset, size, vertexElement);
                offset += size;
                LGL_ERROR;
                    for(size_t l = 0; l < vc->size(); l=l+volumeStepProp_.get()){
                        tgt::vec3 dis = vc->at(l)->getRepresentation<VolumeRAM_3xFloat>()->voxel(i,j,k);
                        vertexElement[0] = dis.x;
                        vertexElement[1] = dis.y;
                        vertexElement[2] = dis.z;
                        vertexElement[3] = center.x;
                        vertexElement[4] = center.y;
                        vertexElement[5] = center.z;
                        vertexElement[6] = static_cast<GLfloat>(l+1);
                        glBufferSubData(GL_ARRAY_BUFFER, offset, size, vertexElement);
                        offset += size;
                    }
            }
        }
    }
    delete[] vertexElement;
}

/*void RadarGlyphRenderer3D::fillIndexVBO() {
    const VolumeCollection* vc = vcInport_.getData();
    tgt::svec3 dim = vc->at(0)->getDimensions();
    tgt::vec3 halfSpacing = vc->at(0)->getSpacing()/2.f;
    tgt::mat4 v2w = vc->at(0)->getVoxelToWorldMatrix();

    GLint size = sizeof(GLfloat)*3;
    GLfloat* vertexElement = new GLfloat[3](); //dir(x,y,z), center(x,y,z), step
    GLint offset = 0;

    int l = activeIndexProp_.get();
    float factor = (l+1)*tgt::min(halfSpacing)/(static_cast<float>(vc->size()));

    for (size_t i = 0; i < dim.x ; ++i) {
        for (size_t j = 0; j < dim.y; ++j) {
            for (size_t k = 0; k < dim.z; ++k){
                tgt::vec4 center = v2w*(tgt::vec4(static_cast<float>(i),static_cast<float>(j),static_cast<float>(k),1.f))+tgt::vec4(halfSpacing,0.f);
                tgt::vec3 dis = vc->at(l)->getRepresentation<VolumeRAM_3xFloat>()->voxel(i,j,k);
                if(tgt::length(dis) > 0.f)
                    dis = tgt::normalize(dis);
                tgt::vec3 point = center.xyz() + (dis*factor);
                vertexElement[0] = point.x;
                vertexElement[1] = point.y;
                vertexElement[2] = point.z;
                glBufferSubData(GL_ARRAY_BUFFER, offset, size, vertexElement);
                offset += size;
            }
        }
    }
    delete[] vertexElement;
}*/

}   // namespace
