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

#include "bodyparts3dmanager.h"
#include "voreen/core/utils/glsl.h"

namespace voreen{

    const std::string BodyParts3DManager::loggerCat_("voreen.touchtable.BodyParts3DRenderer");

    BodyParts3DManager::BodyParts3DManager()
        : RenderProcessor()
        , inPort_(Port::INPORT, "geometry.input", "Geometry Input")
        , pickingPort_(Port::OUTPORT, "picking.port", "Picking Port", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_STATIC)
        , privateRenderPort_(Port::OUTPORT, "private.render.port", "Private Renderport", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_STATIC)
        , textureOutPort_(Port::OUTPORT, "texture.outport", "texture outport")
        , outPort_(Port::OUTPORT, "geometry.output", "Geometry Output")
        , textOutPort_(Port::OUTPORT, "text.outport","Text Outport")
        , canvasSize_("canvasSize", "Canvas Size", tgt::ivec2(256), tgt::ivec2(32), tgt::ivec2(2 << 12))
        , camera_("camera","Camera",tgt::Camera(tgt::vec3(0.0f,-1.0f,0.0f),tgt::vec3(0.0f),tgt::vec3(0.0f,0.0f,1.0f)))
        , enableClipping_("enable.clipping","enable Clipping")
        , planeNormal_("plane.normal","plane normal", tgt::vec3(1.f,0.f,0.f),tgt::vec3(-1.f),tgt::vec3(1.f))
        , planeDistance_("plane.distance","plane distance",0.0f,-10000.0f,10000.0f)
        , invertPlane_("invert.plane","invert plane")
        , lightAmbient_("ambient.light", "ambient Light", tgt::vec4(0.3137f,0.3137f,0.3137f,1.0f))
        , lightDiffuse_("diffuse.light", "diffuse Light", tgt::vec4(0.5f,0.5f,0.5f,1.0f))
        , lightSpecular_("specular.light", "specular Light", tgt::vec4(0.1176f,0.1176f,0.1176f,1.0f))
        , shininess_("shininess", "shininess", 1.0f, 0.0f, 10000.0f)
        , action_("action", "action", false)
        , revert_("revert.action","undo")
    {
        addProperty(canvasSize_);

        addProperty(camera_);

        addProperty(lightAmbient_);
        addProperty(lightDiffuse_);
        addProperty(lightSpecular_);
        addProperty(shininess_);

        addProperty(enableClipping_);
        addProperty(planeNormal_);
        addProperty(planeDistance_);
        addProperty(invertPlane_);

        addProperty(action_);
        addProperty(revert_);

        enableClipping_.setGroupID("clipping");
        planeNormal_.setGroupID("clipping");
        planeDistance_.setGroupID("clipping");
        invertPlane_.setGroupID("clipping");
        setPropertyGroupGuiName("clipping","Clipping");

        lightAmbient_.setGroupID("lighting");
        lightDiffuse_.setGroupID("lighting");
        lightSpecular_.setGroupID("lighting");
        shininess_.setGroupID("lighting");
        setPropertyGroupGuiName("lighting", "Lighting");

        canvasSize_.setVisibleFlag(false);
        planeNormal_.setVisibleFlag(false);
        planeDistance_.setVisibleFlag(false);
        invertPlane_.setVisibleFlag(false);
        //action_.setVisibleFlag(false);

        enableClipping_.onChange(MemberFunctionCallback<BodyParts3DManager>(this, &BodyParts3DManager::updatePropertyVisibilities));
        canvasSize_.onChange(MemberFunctionCallback<BodyParts3DManager>(this, &BodyParts3DManager::adjustPickingRenderSize));
        revert_.onChange(MemberFunctionCallback<BodyParts3DManager>(this, &BodyParts3DManager::restore));

        cameraHandler_ = new CameraInteractionHandler("cameraHandler","Camera", &camera_);
        addInteractionHandler(cameraHandler_);

        highlight_ = -1;

        removedMeshes_ = new TriangleMeshGeometryCollectionBodyParts3D();

        addPort(inPort_);
        addPort(textureOutPort_);
        addPort(outPort_);
        addPort(textOutPort_);
        addPrivateRenderPort(pickingPort_);
        addPrivateRenderPort(privateRenderPort_);
    }

    BodyParts3DManager::~BodyParts3DManager(){
        delete cameraHandler_;
        delete removedMeshes_;
    }

    bool BodyParts3DManager::isReady() const{
        return (outPort_.isReady() && inPort_.isReady() && pickingShader_ && pickingShader_->isLinked());
    }

    void BodyParts3DManager::updatePropertyVisibilities(){
        planeNormal_.setVisibleFlag(enableClipping_.get());
        planeDistance_.setVisibleFlag(enableClipping_.get());
        invertPlane_.setVisibleFlag(enableClipping_.get());
    }

    void BodyParts3DManager::adjustPickingRenderSize(){
        pickingPort_.resize(canvasSize_.get().x, canvasSize_.get().y);
    }

    void BodyParts3DManager::initialize() {
        RenderProcessor::initialize();

        pickingPort_.resize(canvasSize_.get().x, canvasSize_.get().y);
        pickingPort_.changeFormat(GL_RGBA16);

        privateRenderPort_.resize(350, 350);
        privateRenderPort_.changeFormat(GL_RGBA16);

        idManager_.setRenderTarget(pickingPort_.getRenderTarget());
        idManager_.initializeTarget();

        pickingShader_ = ShdrMgr.loadSeparate("bodyparts3d/picking.vert", "bodyparts3d/bp3d_standard.geom", "bodyparts3d/picking.frag", GLSL::generateStandardShaderHeader(), false);
        pickingShader_->rebuild();

        textureShader_ = ShdrMgr.loadSeparate("bodyparts3d/bp3d_touchtablewidget.vert", "bodyparts3d/bp3d_touchtablewidget.frag", GLSL::generateStandardShaderHeader(), false);
        textureShader_->rebuild();
    }

    void BodyParts3DManager::beforeProcess(){
        RenderProcessor::beforeProcess();

        if(inPort_.hasChanged() && inPort_.hasData()){
            textureOutPort_.clear();
            for(std::vector<std::pair<tgt::Texture*, std::string> >::iterator iter = textures_.begin(); iter != textures_.end(); ++iter){
                iter->first->downloadTexture();
                iter->first->setCpuTextureData(0,false);
                delete iter->first;
            }
            textures_.clear();

            //const_cast is ok here because the underlying data is not const and we wont modify it anyway
            TriangleMeshGeometryCollectionBodyParts3D* portData = const_cast<TriangleMeshGeometryCollectionBodyParts3D*>(
                static_cast<const TriangleMeshGeometryCollectionBodyParts3D*>(inPort_.getData()));
            tgtAssert( portData , "Failed to recognize BodyParts Geometry");
            portData_ = portData;
            outPort_.setData(portData_, false);
        }
        if(!inPort_.hasData()){
            outPort_.clear();
        }
    }

    void BodyParts3DManager::process(){
        tgt::mat4 model =  portData_->getTransformationMatrix();
        tgt::mat4 view = camera_.get().getViewMatrix();
        tgt::mat4 projection = camera_.get().getProjectionMatrix(pickingPort_.getSize());

        pickingShader_->activate();

        tgt::plane p = tgt::plane(planeNormal_.get(),planeDistance_.get());
        tgt::mat4 invModelMatrix;
        model.invert(invModelMatrix);
        p = p.transform((invModelMatrix));

        pickingShader_->setUniform("MVP", projection*view*model);
        pickingShader_->setUniform("enableClipping_",enableClipping_.get());
        if(enableClipping_.get()){
            if(invertPlane_.get()){
                pickingShader_->setUniform("plane_", -p.toVec4());
            }else{
                pickingShader_->setUniform("plane_", p.toVec4());
            }
        }
        LGL_ERROR;

        idManager_.activateTarget();

        glClearDepth(1.0f);
        glDepthFunc(GL_LESS);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        if(!enableClipping_.get()){
            glEnable(GL_CULL_FACE);
            glCullFace(GL_BACK);
        }else{
            glDisable(GL_CULL_FACE);
        }

        for(size_t i = 0; i < portData_->getNumMeshes(); ++i){
            idManager_.setGLColor(static_cast<int>(i));
            portData_->getMesh(i)->renderUncolored();
            LGL_ERROR;
        }

        idManager_.deactivateTarget();

        pickingShader_->deactivate();

        glDisable(GL_CULL_FACE);
        LGL_ERROR;

        if(highlight_ >= 0 && highlight_ < portData_->getNumMeshes())
            textOutPort_.setData(portData_->getMesh(highlight_)->getTag());
        else
            textOutPort_.clear();

    }

    void BodyParts3DManager::deinitialize() {
        ShdrMgr.dispose(pickingShader_);
        pickingShader_ = 0;

        RenderProcessor::deinitialize();
    }

    tgt::Texture* BodyParts3DManager::renderGeometryTexture(TriangleMeshGeometryBodyParts3D* geometry){
        tgt::Bounds bbox = geometry->getBoundingBox();
        float width = std::max(std::max(bbox.getURB().y-bbox.getLLF().y, bbox.getURB().z - bbox.getLLF().z)
            , bbox.getURB().x-bbox.getLLF().x)*1.5f;
        //adjust camera to scene
        tgt::Camera camera(tgt::vec3(0.0f,-width,0.0f),tgt::vec3(0.0f),tgt::vec3(0.0f,0.0f,1.0f));
        camera.setFrustum(tgt::Frustum(45.0f,1.0f,0.5f,5000.0f));

        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.pushMatrix();
        MatStack.loadMatrix(tgt::mat4::createPerspective(tgt::PIf/4,1.0f,0.5,5000.0f));//camera.getProjectionMatrix(privateRenderPort_.getSize()));
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.pushMatrix();
        MatStack.loadMatrix(tgt::mat4::createLookAt(tgt::vec3(0.0f,-width,0.0f),tgt::vec3(0.0f),tgt::vec3(0.0f,0.0f,1.0f)));

        textureShader_->activate();

        setGlobalShaderParameters(textureShader_, &camera);
        textureShader_->setIgnoreUniformLocationError(true);

        textureShader_->setUniform("lightSource_.attenuation_", tgt::vec3(1.0f,0.0f,0.0f));
        textureShader_->setUniform("lightSource_.ambientColor_", lightAmbient_.get().xyz());
        textureShader_->setUniform("lightSource_.diffuseColor_", lightDiffuse_.get().xyz());
        textureShader_->setUniform("lightSource_.specularColor_", lightSpecular_.get().xyz());
        textureShader_->setUniform("shininess_", shininess_.get());
        textureShader_->setUniform("lightSource_.position_", camera.getPosition());
        textureShader_->setUniform("cameraPosition_", camera.getPosition());
        textureShader_->setUniform("center_",bbox.center());
        textureShader_->setUniform("vertexColor_", geometry->getColor());

        textureShader_->setIgnoreUniformLocationError(false);

        privateRenderPort_.activateTarget();
        privateRenderPort_.clearTarget();

        glClearDepth(1.0f);
        glDepthFunc(GL_LESS);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        geometry->render();
        LGL_ERROR;

        privateRenderPort_.deactivateTarget();

        textureShader_->deactivate();

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

    void BodyParts3DManager::onEvent(tgt::Event* e){
        if(e->getEventType() == tgt::Event::MOUSERELEASEEVENT){
            if(action_.get())
                removeGeometry(static_cast<tgt::MouseEvent*>(e));
            else
                highlightGeometry(static_cast<tgt::MouseEvent*>(e));
        }
        else if(e->getEventType() == tgt::Event::KEYEVENT &&
            static_cast<tgt::KeyEvent*>(e)->pressed() &&
                static_cast<tgt::KeyEvent*>(e)->keyCode() == tgt::KeyEvent::K_LCTRL)
            addGeometry();
        Processor::onEvent(e);
    }

    void BodyParts3DManager::highlightGeometry(tgt::MouseEvent* e){
        tgt::ivec2 pos = tgt::ivec2(e->x(), e->viewport().y - e->y());
        tgt::col4 color = idManager_.getColorAtPos(pos);
        for(int i = 0; i < portData_->getNumMeshes() ; ++i){
            if(idManager_.getColorFromId(i) == color){
                resetHighlight();
                portData_->getMesh(i)->setHighlight(true);
                highlight_ = i;
                e->accept();
                invalidate();
                break;
            }
        }
    }

    void BodyParts3DManager::removeGeometry(tgt::MouseEvent* e){
        tgt::ivec2 pos = tgt::ivec2(e->x(), e->viewport().y - e->y());
        tgt::col4 color = idManager_.getColorAtPos(pos);
        for(int i = 0; i < portData_->getNumMeshes() ; ++i){
            if(idManager_.getColorFromId(i) == color){
                removedMeshes_->addMesh(portData_->removeMesh(i));
                resetHighlight();
                highlight_ = -1;
                e->accept();
                break;
            }
        }
    }

    void BodyParts3DManager::addGeometry(){
        if(removedMeshes_->getNumMeshes() > 0){
            portData_->addMesh(removedMeshes_->removeMesh(removedMeshes_->getNumMeshes()-1));
            outPort_.setData(portData_, false);
            invalidate();
        }
    }

    void BodyParts3DManager::restore(){

    }

    void BodyParts3DManager::resetHighlight(){
        if(highlight_ >= 0 && highlight_ < portData_->getNumMeshes()){
            portData_->getMesh(highlight_)->setHighlight(false);
            highlight_ = -1;
            invalidate();
        }
    }

} //namespace voreen
