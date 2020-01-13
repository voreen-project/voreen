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

#include "tgt/glmath.h"
#include "tgt/shadermanager.h"
#include "tgt/gpucapabilities.h"
#include "voreen/core/utils/glsl.h"
#include "voreen/core/interaction/camerainteractionhandler.h"

#include <algorithm>

#include "bodyparts3drenderer.h"

namespace voreen{

    const std::string BodyParts3DRenderer::loggerCat_("voreen.touchtable.BodyParts3DRenderer");

    BodyParts3DRenderer::BodyParts3DRenderer()
        : RenderProcessor()
        , lightPosition_("lightPosition","Light Source Position", tgt::vec4(2.3f, 1.5f, 1.5f, 1.0f))
        , lightDiffuse_("diffuseLight", "diffuse Light", tgt::vec4(0.5f,0.5f,0.5f,1.0f))
        , lightSpecular_("specularLight", "specular Light", tgt::vec4(0.1176f,0.1176f,0.1176f,1.0f))
        , lightAmbient_("ambientLight", "ambient Light", tgt::vec4(0.3137f,0.3137f,0.3137f,1.0f))
        , action_("Action", "Action", false)
        , shadeMode_("shading","Shading",Processor::INVALID_PROGRAM)
        , shininess_("shininess","shininess", 1.0f, 1.0f, 100.0f)
        , restoreSignal_("restoreSignal","Restore Signal", -1,-INT_MAX ,INT_MAX)
        , shaderProgram_("BodyParts3Dshader", "BodyParts3D Shader", "bodyparts3d/bp3d_standard.frag", "bodyparts3d/bp3d_standard.vert","bodyparts3d/bp3d_standard.geom")
        , boundingBoxPort_(Port::OUTPORT,"BodyParts3DBoundingBoxPort","BodyParts3D BoundingBoxPort")
        , renderOutPort_(Port::OUTPORT, "image.entrypoints", "Entry-points Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
        , pickingPort_(Port::OUTPORT, "image.tmp" , "image.tmp", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
        , inPort_(Port::INPORT, "trianglemeshgeometry.geometry","Triangle Mesh Geometry Input")
        , textOutPort_(Port::OUTPORT, "textoutport", "Text Outport")
        , texturePort_(Port::OUTPORT,"bodypartstextureoutport","BodyParts Texture Outport")
        , privateRenderPort_(Port::OUTPORT,"bodypartsprivatertexrenderport", "BodyParts texture render port")
        , camera_("camera", "Camera", tgt::Camera(tgt::vec3(1.f,1.f,1.f), tgt::vec3(0.f,0.f,0.f), tgt::vec3(0.f, 1.f, 0.f)))
        , highlight_("highlightid","Hightlight ID", -1, -INT_MAX, INT_MAX,Processor::INVALID_RESULT, NumericProperty<int>::DYNAMIC)
        , revert_("revertcutbutton", "revert cut button", 0)
        , enableClipping_("EnableClipping", "Enable on-the-fly clipping", false)
        , planeNormal_("planeNormal", "Clipping plane normal", tgt::vec3(1.f,0.f,0.f),tgt::vec3(-1.f),tgt::vec3(1.f))
        , planeDistance_("planeDistance", "Clipping plane distance", 0.f, -10000.f, 10000.f)
        , planeManipulatorScale_("planeManipulatorScale", "PlaneManipulatorScale", 1.0f, 0.01f, 1.0f)
        , enablePicking_("enablepicking","enable Picking",true)
        , invertPlane_("invertPlane","invert",false)
    {
        GLSL::fillShadingModesProperty(shadeMode_);
        addProperty(shadeMode_);
        addProperty(camera_);
        addProperty(lightPosition_);
        addProperty(lightDiffuse_);
        addProperty(lightSpecular_);
        addProperty(lightAmbient_);
        addProperty(shaderProgram_);
        addProperty(enablePicking_);
        addProperty(shininess_);
        addProperty(highlight_);
        addProperty(action_);
        addProperty(revert_);
        addProperty(enableClipping_);
        addProperty(invertPlane_);
        addProperty(planeNormal_);
        addProperty(planeDistance_);
        addProperty(restoreSignal_);
        addProperty(planeManipulatorScale_);

        enableClipping_.setGroupID("clipping");
        invertPlane_.setGroupID("clipping");
        planeNormal_.setGroupID("clipping");
        planeDistance_.setGroupID("clipping");
        setPropertyGroupGuiName("clipping", "Clipping Parameters");

        revert_.setVisibleFlag(false);
        highlight_.setVisibleFlag(false);
        restoreSignal_.setVisibleFlag(false);
        planeDistance_.setVisibleFlag(false);
        planeNormal_.setVisibleFlag(false);
        planeManipulatorScale_.setVisibleFlag(false);

        cameraHandler_ = new CameraInteractionHandler("cameraHandler","Camera", &camera_);
        addInteractionHandler(cameraHandler_);

        camera_.setFrustum(tgt::Frustum(45.0f,1.0f,0.1f,20000.0f));

        revert_.onChange(MemberFunctionCallback<BodyParts3DRenderer>(this,&BodyParts3DRenderer::addGeometry));
        enableClipping_.onChange(MemberFunctionCallback<BodyParts3DRenderer>(this, &BodyParts3DRenderer::updatePropertyVisibilities));
        restoreSignal_.onChange(MemberFunctionCallback<BodyParts3DRenderer>(this, &BodyParts3DRenderer::restore));
        enablePicking_.onChange(MemberFunctionCallback<BodyParts3DRenderer>(this, &BodyParts3DRenderer::resetHighlight));


        addPort(inPort_);
        addPort(boundingBoxPort_);
        addPort(texturePort_);
        addPort(renderOutPort_);
        addPort(textOutPort_);
        addPrivateRenderPort(&pickingPort_);
        addPrivateRenderPort(&privateRenderPort_);
    }

    BodyParts3DRenderer::~BodyParts3DRenderer(){
        delete cameraHandler_;
    }

    bool BodyParts3DRenderer::isReady() const{
        return (renderOutPort_.isReady() && inPort_.isReady() && (shaderProgram_.hasValidShader() || shaderProgram_.requiresRebuild()));
    }

    void BodyParts3DRenderer::adjustRenderOutportSizes(){
        tgt::ivec2 size = tgt::ivec2(-1);

        RenderSizeReceiveProperty* entrySizeProp = renderOutPort_.getSizeReceiveProperty();

        if(renderOutPort_.isConnected()) {
            size = entrySizeProp->get();
        }

        if(size != tgt::ivec2(-1)){
            renderOutPort_.resize(size);
            pickingPort_.resize(size);
        }
    }

    std::string BodyParts3DRenderer::generateHeader(const tgt::GpuCapabilities::GlVersion* version){
        std::string headerSource = GLSL::generateStandardShaderHeader(version);
        headerSource += GLSL::getShaderDefine(shadeMode_.get(), "APPLY_SHADING");
        headerSource += "#define TRIANGLE_VEC3\n";
        return headerSource;
    }

    void BodyParts3DRenderer::compile() {
        shaderProgram_.setHeader(generateHeader());
        shaderProgram_.rebuild();
    }

    void BodyParts3DRenderer::initialize() {
        RenderProcessor::initialize();
        compile();
        idManager_.setRenderTarget(pickingPort_.getRenderTarget());
        idManager_.initializeTarget();

        privateRenderPort_.resize(tgt::ivec2(350,350));

        renderOutPort_.changeFormat(GL_RGBA16);
        pickingPort_.changeFormat(GL_RGBA16);
        privateRenderPort_.changeFormat(GL_RGBA16);

        pickingShaderProgram_ = ShdrMgr.loadSeparate("bodyparts3d/picking.vert", "bodyparts3d/bp3d_standard.geom","bodyparts3d/picking.frag", GLSL::generateStandardShaderHeader(), false);
        textureShaderProgram_ = ShdrMgr.loadSeparate("bodyparts3d/bp3d_touchtablewidget.vert", "bodyparts3d/bp3d_touchtablewidget.frag", generateHeader(), false);
        textureShaderProgram_->rebuild();
    }

    void BodyParts3DRenderer::beforeProcess(){
        RenderProcessor::beforeProcess();

        // compile program if needed
        if (getInvalidationLevel() >= Processor::INVALID_PROGRAM) {
            PROFILING_BLOCK("compile");
            compile();
        }
        LGL_ERROR;

        RenderPort& refPort = renderOutPort_;

        if(inPort_.hasChanged() && inPort_.hasData() && inPort_.getData()->size() > 0){
            const std::vector<std::pair<TriangleMeshGeometrySingleColor*,std::string> >* input = inPort_.getData();
            //reset bodyparts data stored
            portInput_.clear();
            removedParts_.clear();
            texturePort_.clear();
            for(std::vector<std::pair<tgt::Texture*, std::string> >::iterator texIter = textures_.begin(); texIter != textures_.end() ; ++texIter){
                texIter->first->downloadTexture();
                texIter->first->setCpuTextureData(0,false);
                delete texIter->first;
            }
            textures_.clear();
            portInput_.reserve(input->size());
            for(size_t i = 0; i < input->size() ; i++){
                portInput_.push_back(input->at(i));
            }
            texturePort_.clear();

            calculateBoundingBox();

            tgt::vec3 urb = boundingBox_.getURB();
            tgt::vec3 llf = boundingBox_.getLLF();

            //adjust planemanipulator size to boundingbox
            planeManipulatorScale_.set(1.6f*(std::min(std::min(urb.x-llf.x,urb.y-llf.y),urb.z-llf.z) + tgt::length(boundingBox_.diagonal()))/(2.0f*tgt::length(boundingBox_.diagonal())));

            //scale boundingbox down so that the longest side is bounded by [-1,1]
            float scale = 2.0f/(boundingBox_.getURB().z-boundingBox_.getLLF().z);
            transformationMatrix_ = tgt::mat4::createScale(tgt::vec3(scale,scale,scale))*tgt::mat4::createTranslation(-boundingBox_.center());

            //transform boundingbox accordingly
            boundingBox_ = boundingBox_.transform(transformationMatrix_);
            LINFO(boundingBox_);

            //Create an empty mesh to represent the boundingbox and set boundingbox data
            TriangleMeshGeometrySimple* t = new TriangleMeshGeometrySimple();
            t->addTriangle(Triangle<VertexBase>(VertexBase(boundingBox_.getLLF()),VertexBase(boundingBox_.getURB()),VertexBase(boundingBox_.center())));
            boundingBoxPort_.setData(t);

            //adjust camera to scene
            camera_.adaptInteractionToScene(boundingBox_);
            camera_.getTrackball().setCenter(tgt::vec3(0.0f));
            camera_.resetCameraFocusToTrackballCenter();
            float distance = std::max(std::max(boundingBox_.getURB().x-boundingBox_.getLLF().x,boundingBox_.getURB().y-boundingBox_.getLLF().y),boundingBox_.getURB().z-boundingBox_.getLLF().z);
            camera_.setPosition(tgt::vec3(0.0f,-1.5f * distance,0.0f));
            camera_.setUpVector(tgt::vec3(1.0f,0.0f,0.0f));

            resetHighlight();
        }
    }

    void BodyParts3DRenderer::process(){
        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.pushMatrix();
        MatStack.loadMatrix(camera_.get().getProjectionMatrix(renderOutPort_.getSize()));
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.pushMatrix();
        MatStack.loadMatrix(camera_.get().getViewMatrix());

        renderGeometry();
        LGL_ERROR;

        if(enablePicking_.get()){
            renderPickingPort();
            LGL_ERROR;
        }

        //set description for highlighted geometry
        if(highlight_.get() == -1 || highlight_.get() >= portInput_.size()){
            textOutPort_.setData("");
            highlight_.set(-1);
        }else{
            textOutPort_.setData(portInput_[highlight_.get()].second);
        }

        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.popMatrix();
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.popMatrix();
        LGL_ERROR;
    }

    void BodyParts3DRenderer::deinitialize() {
        //free resources
        ShdrMgr.dispose(pickingShaderProgram_);
        pickingShaderProgram_ = 0;

        ShdrMgr.dispose(textureShaderProgram_);
        textureShaderProgram_ = 0;

        RenderProcessor::deinitialize();
    }

    void BodyParts3DRenderer::renderGeometry(){
        tgt::Shader* shader = shaderProgram_.getShader();
        shader->activate();
        shader->setIgnoreUniformLocationError(true);
        shader->setUniform("transformationMatrix_",transformationMatrix_);
        shader->setUniform("lightSource_.position_", (camera_.get().getViewMatrix() * tgt::vec4(lightPosition_.get().xyz(), 0.0)).xyz());
        shader->setUniform("lightSource_.attenuation_", tgt::vec3(1.0f,0.0f,0.0f));
        shader->setUniform("lightSource_.ambientColor_", lightAmbient_.get().xyz());
        shader->setUniform("lightSource_.diffuseColor_", lightDiffuse_.get().xyz());
        shader->setUniform("lightSource_.specularColor_", lightSpecular_.get().xyz());
        shader->setUniform("shininess_", shininess_.get());
        shader->setUniform("enableClipping_",enableClipping_.get());
        if(enableClipping_.get()){
            if(invertPlane_.get()){
                shader->setUniform("plane_",tgt::vec4(-planeNormal_.get(),-planeDistance_.get()));
            }else{
                shader->setUniform("plane_",tgt::vec4(planeNormal_.get(),planeDistance_.get()));
            }
        }
        shader->setIgnoreUniformLocationError(false);
        LGL_ERROR;
        setGlobalShaderParameters(shader, &camera_.get());

        renderOutPort_.activateTarget();

        glClearDepth(1.0f);
        glDepthFunc(GL_LESS);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        //enable culling only if clipping is disabled
        //this is necessary because the clipped geometry will not be closed!
        if(!enableClipping_.get()){
            glEnable(GL_CULL_FACE);
            glCullFace(GL_BACK);
        }else{
            glDisable(GL_CULL_FACE);
        }

        for(int i = 0; i < portInput_.size() ; i++){
            shader->setUniform("vertexColor_", portInput_[i].first->getColor() );
            shader->setUniform("highlight_", i == highlight_.get() );
            portInput_[i].first->render();
            LGL_ERROR;
        }

        renderOutPort_.deactivateTarget();

        shader->deactivate();

        glDisable(GL_CULL_FACE);
        LGL_ERROR;
    }

    void BodyParts3DRenderer::renderPickingPort(){
        pickingShaderProgram_->activate();

        pickingShaderProgram_->setIgnoreUniformLocationError(true);
        pickingShaderProgram_->setUniform("transformationMatrix_",transformationMatrix_);
        pickingShaderProgram_->setUniform("enableClipping_",enableClipping_.get());
        if(enableClipping_.get()){
            if(invertPlane_.get()){
                pickingShaderProgram_->setUniform("plane_",tgt::vec4(-planeNormal_.get(),-planeDistance_.get()));
            }else{
                pickingShaderProgram_->setUniform("plane_",tgt::vec4(planeNormal_.get(),planeDistance_.get()));
            }
        }
        pickingShaderProgram_->setIgnoreUniformLocationError(false);

        tgt::Camera cam = camera_.get();
        setGlobalShaderParameters(pickingShaderProgram_, &cam);
        LGL_ERROR;

        pickingPort_.activateTarget();

        glClearDepth(1.0f);
        glDepthFunc(GL_LESS);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        //enable culling only if clipping is disabled
        //this is necessary because the clipped geometry will not be closed!
        if(!enableClipping_.get()){
            glEnable(GL_CULL_FACE);
            glCullFace(GL_BACK);
        }else{
            glDisable(GL_CULL_FACE);
        }

        for(int i = 0; i < portInput_.size() ; i++){
            idManager_.setGLColor(i);
            portInput_[i].first->render();
            LGL_ERROR;
        }

        pickingPort_.deactivateTarget();

        pickingShaderProgram_->deactivate();

        glDisable(GL_CULL_FACE);
        LGL_ERROR;
    }

    tgt::Texture* BodyParts3DRenderer::renderGeometryTexture(TriangleMeshGeometrySingleColor* geometry){
        tgt::Bounds bbox = geometry->getBoundingBox();
        float width = std::max(std::max(bbox.getURB().y-bbox.getLLF().y, bbox.getURB().z - bbox.getLLF().z)
            , bbox.getURB().x-bbox.getLLF().x)*1.5f;
        //adjust camera to scene
        tgt::Camera camera(tgt::vec3(0.0f,-width,0.0f),tgt::vec3(0.0f),tgt::vec3(0.0f,0.0f,1.0f));
        camera.setFrustum(tgt::Frustum(45.0f,1.0f,0.5f,5000.0f));

        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.pushMatrix();
        MatStack.loadMatrix(camera.getProjectionMatrix(privateRenderPort_.getSize()));
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.pushMatrix();
        MatStack.loadMatrix(camera.getViewMatrix());

        textureShaderProgram_->activate();

        setGlobalShaderParameters(textureShaderProgram_, &camera);
        textureShaderProgram_->setIgnoreUniformLocationError(true);

        textureShaderProgram_->setUniform("lightSource_.attenuation_", tgt::vec3(1.0f,0.0f,0.0f));
        textureShaderProgram_->setUniform("lightSource_.ambientColor_", lightAmbient_.get().xyz());
        textureShaderProgram_->setUniform("lightSource_.diffuseColor_", lightDiffuse_.get().xyz());
        textureShaderProgram_->setUniform("lightSource_.specularColor_", lightSpecular_.get().xyz());
        textureShaderProgram_->setUniform("shininess_", shininess_.get());
        textureShaderProgram_->setUniform("lightSource_.position_", camera.getPosition());
        textureShaderProgram_->setUniform("cameraPosition_", camera.getPosition());
        textureShaderProgram_->setUniform("center_",bbox.center());
        textureShaderProgram_->setUniform("vertexColor_", geometry->getColor());

        textureShaderProgram_->setIgnoreUniformLocationError(false);

        privateRenderPort_.activateTarget();
        privateRenderPort_.clearTarget();

        glClearDepth(1.0f);
        glDepthFunc(GL_LESS);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        geometry->render();
        LGL_ERROR;

        privateRenderPort_.deactivateTarget();

        textureShaderProgram_->deactivate();

        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.popMatrix();
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.popMatrix();
        LGL_ERROR;

        //create a copy of the texture
        int texwidth = privateRenderPort_.getSize().x;
        int texheight = privateRenderPort_.getSize().y;
        tgt::Texture* tex = privateRenderPort_.getColorTexture();
        tex->downloadTexture();
        tgt::Texture* newtex = new tgt::Texture(tgt::ivec3(texwidth,texheight,1),GL_RGBA,GL_RGBA,GL_FLOAT,tgt::Texture::LINEAR,tgt::Texture::REPEAT,tex->getCpuTextureData(),false);
        newtex->uploadTexture();
        return newtex;
    }

    void BodyParts3DRenderer::updatePropertyVisibilities(){
        //show clipping properties only if clipping is enabled
        invertPlane_.setVisibleFlag(enableClipping_.get());
        planeDistance_.setVisibleFlag(enableClipping_.get());
        planeNormal_.setVisibleFlag(enableClipping_.get());
    }

    void BodyParts3DRenderer::onEvent(tgt::Event* e){

        if(e->getEventType() == tgt::Event::MOUSEPRESSEVENT){
            if(action_.get())
                removeGeometry(dynamic_cast<tgt::MouseEvent*>(e));
            else
                highlightGeometry(dynamic_cast<tgt::MouseEvent*>(e));
        }
        else if(e->getEventType() == tgt::Event::KEYEVENT)
            addGeometry();
        Processor::onEvent(e);
    }

    void BodyParts3DRenderer::highlightGeometry(tgt::MouseEvent* e){
        tgt::ivec2 pos = tgt::ivec2(e->x(), e->viewport().y - e->y());
        tgt::col4 color = idManager_.getColorAtPos(pos);
        for(int i = 0; i < portInput_.size() ; i++){
            if(idManager_.getColorFromId(i) == color){
                highlight_.set(i);
                e->accept();
                break;
            }
        }
        invalidate();
    }

    void BodyParts3DRenderer::removeGeometry(tgt::MouseEvent* e){
        tgt::ivec2 pos = tgt::ivec2(e->x(), e->viewport().y - e->y());
        tgt::col4 color = idManager_.getColorAtPos(pos);
        for(int i = 0; i < portInput_.size() ; i++){
            if(idManager_.getColorFromId(i) == color){
                texturePort_.clear();
                removedParts_.push_back(portInput_.at(i));
                textures_.push_back(std::make_pair(renderGeometryTexture(portInput_.at(i).first),portInput_.at(i).second));
                portInput_.erase(portInput_.begin()+i);
                highlight_.set(-1);
                texturePort_.setData(&textures_,false);
                e->accept();
                break;
            }
        }
        invalidate();
    }

    void BodyParts3DRenderer::addGeometry(){
        if(removedParts_.size() != 0){
            texturePort_.clear();
            portInput_.push_back(removedParts_.at(removedParts_.size()-1));
            removedParts_.erase(removedParts_.end()-1);
            textures_.back().first->downloadTexture();
            delete textures_.back().first;
            textures_.erase(textures_.end()-1);
            if(removedParts_.size() > 0){ //dont set portdata if there is none!
                texturePort_.setData(&textures_,false);
            }
        }
        invalidate();
    }

    void BodyParts3DRenderer::restore(){
        if(restoreSignal_.get() >= 0 && restoreSignal_.get() < removedParts_.size() && removedParts_.size() != 0){
            texturePort_.clear();
            portInput_.push_back(removedParts_.at(restoreSignal_.get()));
            removedParts_.erase(removedParts_.begin() + restoreSignal_.get());

            textures_.at(restoreSignal_.get()).first->downloadTexture();
            delete textures_.at(restoreSignal_.get()).first;
            textures_.erase(textures_.begin() + restoreSignal_.get());

            if(textures_.size() > 0){ //dont set portdata if there is none!
                texturePort_.setData(&textures_,false);
            }
        }
        invalidate();
    }

    void BodyParts3DRenderer::resetHighlight(){
        highlight_.set(-1);
    }

    void BodyParts3DRenderer::calculateBoundingBox(){
        if(portInput_.size() > 0){
            float upper, lower, left, right, front, back;
            tgt::Bounds tmpBoundingBox = inPort_.getData()->at(0).first->getBoundingBox(false);
            upper = tmpBoundingBox.getURB().x;
            right = tmpBoundingBox.getURB().y;
            back = tmpBoundingBox.getURB().z;
            lower = tmpBoundingBox.getLLF().x;
            left = tmpBoundingBox.getLLF().y;
            front = tmpBoundingBox.getLLF().z;
            for(size_t i = 1; i < inPort_.getData()->size(); i++){
                //for the urb simply take the maximum in each direction and
                // for the llf the minimum
                tmpBoundingBox = inPort_.getData()->at(i).first->getBoundingBox(false);
                upper = std::max(upper,tmpBoundingBox.getURB().x);
                right = std::max(right,tmpBoundingBox.getURB().y);
                back = std::max(back,tmpBoundingBox.getURB().z);
                lower = std::min(lower,tmpBoundingBox.getLLF().x);
                left = std::min(left,tmpBoundingBox.getLLF().y);
                front = std::min(front,tmpBoundingBox.getLLF().z);
            }
            boundingBox_ = tgt::Bounds(tgt::vec3(upper,right,back),tgt::vec3(lower,left,front));
        }else{
            // a default boundingbox
            boundingBox_ = tgt::Bounds(tgt::vec3(1.0f),tgt::vec3(-1.0f));
        }
    }

} //namespace voreen
