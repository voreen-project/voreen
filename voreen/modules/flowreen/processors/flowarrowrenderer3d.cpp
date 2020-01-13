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

//header
#include "flowarrowrenderer3d.h"
//volume
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"
//texture handling
#include "tgt/texturemanager.h"
#include "tgt/textureunit.h"
#include "tgt/camera.h"
//shader support test
#include "tgt/gpucapabilities.h"
#include "voreen/core/utils/glsl.h"

#include <iostream>
namespace voreen {

    const std::string FlowArrowRenderer3D::loggerCat_("voreen.bovenkamp.FlowArrowRenderer3D");

FlowArrowRenderer3D::FlowArrowRenderer3D()
    : RenderProcessor()
    , volumeInport_(Port::INPORT, "volumeInport","Velocity Inport")
    , renderOutport_(Port::OUTPORT, "renderOutport", "Arrow Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , shaderProp_("arrow.prg", "Arrow shader", "arrowshader.frag", "arrowshader.vert", "arrowshader.geom")
    , program_(0)
    , cameraProp_("cameraProp","Camera", tgt::Camera(tgt::vec3(0.0f, 0.0f, 50.f), tgt::vec3(0.0f, 0.0f, 0.0f),
                                         tgt::vec3(0.0f, 1.0f, 0.0f), 45.f,1.f,0.1f,500.f))
    , cameraHandler_(0)

    , arrowColorMode_("arrowColorMode","Derive arrow color...")
    , arrowComponentMode_("arrowComponentMode","Use from data...")
    , arrowSize_("arrowSize","Arrow Size",1.f,0.01f,10.f)
    , tfProp_("tfProp","Transfer Function")
    , xVoxelOffset_("xVoxelOffset", "X Shift",0.5f,0.f,1.f)
    , yVoxelOffset_("yVoxelOffset", "Y Shift",0.5f,0.f,1.f)
    , zVoxelOffset_("zVoxelOffset", "Z Shift",0.5f,0.f,1.f)

    , gpuGsSupport_(true), gpuErrorString_("")
{
    addPort(volumeInport_);
    addPort(renderOutport_);

    addProperty(arrowColorMode_);
        arrowColorMode_.addOption("fromtf", "from Arrow Length", FlowArrowRenderer3D::ARROW_LENGTH);
        arrowColorMode_.addOption("fromdir", "from Arrow Direction", FlowArrowRenderer3D::ARROW_DIRECTION);
        arrowColorMode_.setGroupID("arrow");
    addProperty(arrowSize_);
        arrowSize_.setGroupID("arrow");
    addProperty(tfProp_);
        tfProp_.setGroupID("arrow");
    addProperty(arrowComponentMode_);
        arrowComponentMode_.addOption("ac_all", "all Dimensions", FlowArrowRenderer3D::AC_ALL_DIM);
        arrowComponentMode_.addOption("ac_only_x", "only X Dimension", FlowArrowRenderer3D::AC_ONLY_X_DIM);
        arrowComponentMode_.addOption("ac_only_y", "only Y Dimension", FlowArrowRenderer3D::AC_ONLY_Y_DIM);
        arrowComponentMode_.addOption("ac_only_z", "only Z Dimension", FlowArrowRenderer3D::AC_ONLY_Z_DIM);
        arrowComponentMode_.setGroupID("arrow");
    setPropertyGroupGuiName("arrow","Arrow Color Settings");

    addProperty(xVoxelOffset_);
        xVoxelOffset_.setGroupID("pos");
    addProperty(yVoxelOffset_);
        yVoxelOffset_.setGroupID("pos");
    addProperty(zVoxelOffset_);
        zVoxelOffset_.setGroupID("pos");
    setPropertyGroupGuiName("pos","Arrow Position Settings");


    addProperty(cameraProp_);
    addProperty(shaderProp_);

    cameraHandler_ = new CameraInteractionHandler("cameraHandler", "Camera Handler", &cameraProp_);
    addInteractionHandler(cameraHandler_);

    ON_PROPERTY_CHANGE(arrowColorMode_,FlowArrowRenderer3D,changePropertyVisibility);
}

FlowArrowRenderer3D::~FlowArrowRenderer3D() {
    delete cameraHandler_;
}

Processor* FlowArrowRenderer3D::create() const {
    return new FlowArrowRenderer3D();
}


bool FlowArrowRenderer3D::isReady() const{
    if(!volumeInport_.isReady())
        return false;

    if(!program_ && gpuGsSupport_){
        LWARNING("No shader program compiled.");
        return false;
    }

    return true;
}


std::string FlowArrowRenderer3D::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string header = RenderProcessor::generateHeader();
    header += "#extension GL_EXT_geometry_shader4 : enable\n";
    header += "out vec4 FragData1;\n";
    return header;
}

void FlowArrowRenderer3D::initialize() {
    RenderProcessor::initialize();
    //test if GS is supportet
    std::stringstream strstr;
    if(GpuCaps.getShaderModel() <= tgt::GpuCapabilities::SHADER_MODEL_4){
        gpuGsSupport_ = false;
        strstr << "No GS support. GPU supports ShaderModel "<< GpuCaps.getShaderModel()
               << ".0, but 4.0 is needed. Please choose another render mode.";
        gpuErrorString_ = strstr.str();
    } else {
        if(GpuCaps.getMaxGeometryShaderVertices() < 30){
            gpuGsSupport_ = false;
            strstr << "The GS supports just " << GpuCaps.getMaxGeometryShaderVertices()
                   << " Vertices, but 30 are needed. Please choose another render mode.";
            gpuErrorString_ = strstr.str();
        }
    }

    if(gpuGsSupport_)
        compile();
}

void FlowArrowRenderer3D::deinitialize() {
    RenderProcessor::deinitialize();
}

void FlowArrowRenderer3D::compile() {
        shaderProp_.setHeader(generateHeader());
        shaderProp_.rebuild();
        if(shaderProp_.hasValidShader())
            program_ = shaderProp_.getShader();
        else
            program_ = 0;
    }

void FlowArrowRenderer3D::beforeProcess() {
        RenderProcessor::beforeProcess();

        // compile program if needed
        if (getInvalidationLevel() >= Processor::INVALID_PROGRAM && gpuGsSupport_) {
            PROFILING_BLOCK("compile");
            compile();
        }

        //adapt camera
        if(volumeInport_.hasChanged() && volumeInport_.hasData())
            cameraProp_.adaptInteractionToScene(volumeInport_.getData()->getBoundingBox().getBoundingBox());

        LGL_ERROR;
    }

void FlowArrowRenderer3D::changePropertyVisibility() {
    //tfProp_.setVisibleFlag(!arrowColorMode_.getValue());
}
//---------------------------------------------------------------------
//      process
//---------------------------------------------------------------------
void FlowArrowRenderer3D::process() {
    //load volume from inport
    const VolumeBase* volume = volumeInport_.getData();
    tgtAssert(volume, "no volume");
    const VolumeRAM_3xFloat* velocity = dynamic_cast<const VolumeRAM_3xFloat*>(volume->getRepresentation<VolumeRAM>());
    if(!velocity) {
        LERROR("Expected a 3xFloat volume!");
        return;
    }
    tgtAssert(velocity,"Nullpointer in volumeInport");
    //get Maximal Displacement
    float maxLength = velocity->maxNormalizedMagnitude();
    if(maxLength <= 0) return;

    //activate outport
    renderOutport_.activateTarget();
    renderOutport_.clearTarget();

    //set camera and lightning
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    cameraProp_.look(renderOutport_.getSize());
    glLightfv(GL_LIGHT0, GL_POSITION, cameraProp_.get().getPosition().elem);

    //bind tftexture
    tgt::TextureUnit tfUnit;
    tgt::Texture* tex = tfProp_.get()->getTexture();
    tfUnit.activate();
    tex->bind();

    //check, if GS should be used
    if(gpuGsSupport_){
        tgtAssert(program_,"No shader program compiled.");
        //activate shader
        glLightfv(GL_LIGHT0, GL_POSITION, tgt::vec4(cameraProp_.get().getPosition(),1.f).elem);
        program_->activate();
        program_->setUniform("voxelToWorldMatrix",volume->getVoxelToWorldMatrix());
        program_->setUniform("maxVelocity", maxLength);
        program_->setUniform("ColorTexture", tfUnit.getUnitNumber());
        program_->setUniform("hight_", tgt::min(volume->getSpacing())*arrowSize_.get());
        program_->setUniform("colorFromDir_", arrowColorMode_.getValue() == ARROW_LENGTH ? false : true);
        //position modified
        tgt::vec3 modPos(xVoxelOffset_.get(),yVoxelOffset_.get(),zVoxelOffset_.get());
        glBegin(GL_POINTS);
        //render arrows
        VRN_FOR_EACH_VOXEL(pos,tgt::svec3::zero,volume->getDimensions()){
            tgt::vec3 voxel = velocity->voxel(pos);
            switch(arrowComponentMode_.getValue()) {
            case AC_ONLY_X_DIM:
                voxel.y = 0.f; voxel.z = 0.f;
                break;
            case AC_ONLY_Y_DIM:
                voxel.x = 0.f; voxel.z = 0.f;
                break;
            case AC_ONLY_Z_DIM:
                voxel.x = 0.f; voxel.y = 0.f;
                break;
            default:
            break;
            }
            glNormal3fv(voxel.elem);
            glVertex3fv((tgt::vec3(pos)+modPos).elem);
        }
        glEnd();

        //deactivate shader
        program_->deactivate();
    } else {
        LERROR(gpuErrorString_);
    }

    //restore openGl context
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();
    //deactivate outport
    renderOutport_.deactivateTarget();
    LGL_ERROR;
}

} // namespace voreen
