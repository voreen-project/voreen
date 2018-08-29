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

#include "orientationoverlay.h"
#include "voreen/core/interaction/camerainteractionhandler.h"
#include "voreen/core/voreenapplication.h"

#include "tgt/texturemanager.h"
#include "tgt/gpucapabilities.h"
#include "tgt/glmath.h"
#include "tgt/textureunit.h"
#include "tgt/matrix.h"

#include "voreen/core/datastructures/geometry/trianglemeshgeometryindexed.h"

namespace voreen {

const std::string OrientationOverlay::loggerCat_("voreen.OrientationOverlay");

OrientationOverlay::OrientationOverlay()
    : ImageProcessor("image/orientationoverlay") //determines .frag shader
    //ports
    , inport_(Port::INPORT, "image.input", "Image Input")
    , outport_(Port::OUTPORT, "image.output", "Image Output")
    , privatePort_(Port::OUTPORT, "image.tmp", "image.tmp", false)
    //basic
    , enableProp_("enableProp","Enable",true)
    , orientationTypeProp_("orientationTypeProp","Overlay Type")
    , cameraProp_("camera", "Camera", tgt::Camera(tgt::vec3(0.f, 0.f, 3.5f), tgt::vec3(0.f, 0.f, 0.f), tgt::vec3(0.f, 1.f, 0.f)),
              false,true,500.f,Processor::INVALID_RESULT,Property::LOD_DEBUG)
    //position
    , shiftXProp_("shiftX", "Horizontal Position", 0.85f, 0.0f, 1.0f)
    , shiftYProp_("shiftY", "Vertical Position", 0.15f, 0.0f, 1.0f)
    //axis
    , axisSizeProp_("axisLength", "Axes Length", 1.f, 0.5f, 3.f,Processor::VALID)
    , axesAlignmentProp_("axesAlignmentProp", "Axes Alignment")
    , renderAxesLabelsProp_("renderAxesLabelsProp","Show Labels?",true)
    //cube
    , cubeSizeProp_("cubeSize", "Cube Size", 1.f, 0.5f, 3.f)
    //textures
    , filenameFrontProp_("filenameFront", "Front Texture", "Select texture",
                        VoreenApplication::app()->getUserDataPath("textures"), "*.jpg;*.png;*.bmp",
                        FileDialogProperty::OPEN_FILE)
    , filenameBackProp_("filenameBack", "Back Texture", "Select texture",
                        VoreenApplication::app()->getUserDataPath("textures"), "*.jpg;*.png;*.bmp",
                        FileDialogProperty::OPEN_FILE)
    , filenameTopProp_("filenameTop", "Top Texture", "Select texture",
                        VoreenApplication::app()->getUserDataPath("textures"), "*.jpg;*.png;*.bmp",
                        FileDialogProperty::OPEN_FILE)
    , filenameBottomProp_("filenameBottom", "Bottom Texture", "Select texture",
                        VoreenApplication::app()->getUserDataPath("textures"), "*.jpg;*.png;*.bmp",
                        FileDialogProperty::OPEN_FILE)
    , filenameLeftProp_("filenameLeft", "Left Texture", "Select texture",
                        VoreenApplication::app()->getUserDataPath("textures"), "*.jpg;*.png;*.bmp",
                        FileDialogProperty::OPEN_FILE)
    , filenameRightProp_("filenameRight", "Right Texture", "Select texture",
                        VoreenApplication::app()->getUserDataPath("textures"), "*.jpg;*.png;*.bmp",
                        FileDialogProperty::OPEN_FILE)
    , shaderProp_("geometry.prg", "Shader", "trianglemesh.frag", "trianglemesh.vert", "trianglemesh.geom",
                  Processor::INVALID_PROGRAM,Property::LOD_DEBUG)
    //helper
    , geometryMustBeRecreated_(true)
    , currentGeometry_(0), letterXGeometry_(0), letterYGeometry_(0), letterZGeometry_(0)
    , overlayBaseLength_(0.16f)
    , frontTex_(0), backTex_(0), topTex_(0), leftTex_(0), bottomTex_(0), rightTex_(0)
    , reloadTextures_(false), loadingTextures_(false)
{
    //ports
    addPort(inport_);
    addPort(outport_);
    addPrivateRenderPort(&privatePort_);
    //basic
    addProperty(enableProp_);
    addProperty(orientationTypeProp_);
        orientationTypeProp_.addOption("axis2d","Axis 2D (Slices)",OT_AXES_2D);
        orientationTypeProp_.addOption("axis3d","Axis 3D",OT_AXES_3D);
        orientationTypeProp_.addOption("col_cube","Colored Cube",OT_COLOR_CUBE);
        orientationTypeProp_.addOption("tex_cube","Textured Cube",OT_TEXTURE_CUBE);
        orientationTypeProp_.addOption("col_tex_cube","Colored/Textured Cube",OT_COLOR_TEXTURE_CUBE);
        orientationTypeProp_.onChange(MemberFunctionCallback<OrientationOverlay>(this, &OrientationOverlay::adjustPropertyVisibility));
    addProperty(cameraProp_);
    //position
    addProperty(shiftXProp_);
    shiftXProp_.setGroupID("position");
    addProperty(shiftYProp_);
    shiftYProp_.setGroupID("position");
    setPropertyGroupGuiName("position","Position");
    //axis
    addProperty(axisSizeProp_);
    axisSizeProp_.setGroupID("axis");
    axisSizeProp_.onChange(MemberFunctionCallback<OrientationOverlay>(this, &OrientationOverlay::invalidateGeometry));
    addProperty(axesAlignmentProp_);
    axesAlignmentProp_.addOption("xy-plane", "XY-Plane (axial)", XY_PLANE);
    axesAlignmentProp_.addOption("xz-plane", "XZ-Plane (coronal)", XZ_PLANE);
    axesAlignmentProp_.addOption("yz-plane", "YZ-Plane (sagittal)", YZ_PLANE);
    axesAlignmentProp_.onChange(MemberFunctionCallback<OrientationOverlay>(this, &OrientationOverlay::invalidateGeometry));
    axesAlignmentProp_.setGroupID("axis");
    addProperty(renderAxesLabelsProp_);
    renderAxesLabelsProp_.setGroupID("axis");
    setPropertyGroupGuiName("axis","Axis Settings");
    //cube
    addProperty(cubeSizeProp_);
    cubeSizeProp_.onChange(MemberFunctionCallback<OrientationOverlay>(this, &OrientationOverlay::invalidateGeometry));
    cubeSizeProp_.setGroupID("cube");
    setPropertyGroupGuiName("cube","Cube Settings");
    //textures
    addProperty(filenameFrontProp_);
    filenameFrontProp_.onChange(MemberFunctionCallback<OrientationOverlay>(this, &OrientationOverlay::reloadTextures));
    filenameFrontProp_.setGroupID("texture");
    addProperty(filenameBackProp_);
    filenameBackProp_.onChange(MemberFunctionCallback<OrientationOverlay>(this, &OrientationOverlay::reloadTextures));
    filenameBackProp_.setGroupID("texture");
    addProperty(filenameTopProp_);
    filenameTopProp_.onChange(MemberFunctionCallback<OrientationOverlay>(this, &OrientationOverlay::reloadTextures));
    filenameTopProp_.setGroupID("texture");
    addProperty(filenameBottomProp_);
    filenameBottomProp_.onChange(MemberFunctionCallback<OrientationOverlay>(this, &OrientationOverlay::reloadTextures));
    filenameBottomProp_.setGroupID("texture");
    addProperty(filenameLeftProp_);
    filenameLeftProp_.onChange(MemberFunctionCallback<OrientationOverlay>(this, &OrientationOverlay::reloadTextures));
    filenameLeftProp_.setGroupID("texture");
    addProperty(filenameRightProp_);
    filenameRightProp_.onChange(MemberFunctionCallback<OrientationOverlay>(this, &OrientationOverlay::reloadTextures));
    filenameRightProp_.setGroupID("texture");
    setPropertyGroupGuiName("texture","Texture Settings");
    addProperty(shaderProp_);

    // set initial texture names
    std::string texturePath = VoreenApplication::app()->getCoreResourcePath("textures");
    textureNames_[0] = texturePath + "/axial_t.png";
    textureNames_[1] = texturePath + "/axial_b.png";
    textureNames_[2] = texturePath + "/coronal_f.png";
    textureNames_[3] = texturePath + "/coronal_b.png";
    textureNames_[4] = texturePath + "/sagittal_l.png";
    textureNames_[5] = texturePath + "/sagittal_r.png";

    adjustPropertyVisibility();
}

OrientationOverlay::~OrientationOverlay() {

}

Processor* OrientationOverlay::create() const {
    return new OrientationOverlay();
}

void OrientationOverlay::initialize() {
    ImageProcessor::initialize();
    loadTextures();
}

void OrientationOverlay::deinitialize() {
    // dispose textures before sending this processor into nirvana
    if(frontTex_)
        TexMgr.dispose(frontTex_);
    if(backTex_)
        TexMgr.dispose(backTex_);
    if(bottomTex_)
        TexMgr.dispose(bottomTex_);
    if(leftTex_)
        TexMgr.dispose(leftTex_);
    if(topTex_)
        TexMgr.dispose(topTex_);
    if(rightTex_)
        TexMgr.dispose(rightTex_);

    //delete geom
    delete currentGeometry_;
    delete letterXGeometry_;
    delete letterYGeometry_;
    delete letterZGeometry_;

    ImageProcessor::deinitialize();
}

bool OrientationOverlay::isReady() const {
    return outport_.isReady();
}

//--------------------------------------------------------------------------------------------
//          process functions
//--------------------------------------------------------------------------------------------
void OrientationOverlay::beforeProcess() {
    ImageProcessor::beforeProcess();
    if (reloadTextures_)
        loadTextures();
    if(!shaderProp_.hasValidShader() || getInvalidationLevel() >= Processor::INVALID_PROGRAM) {
        shaderProp_.setHeader(generateHeader() + (currentGeometry_ ? currentGeometry_->getShaderDefines(): ""));
        shaderProp_.rebuild();
    }
}

void OrientationOverlay::process() {
    //shader must be valid
    if(!shaderProp_.hasValidShader()) {
        LERROR("Shader for geometry failed to compile");
        return;
    }

    if(geometryMustBeRecreated_)
        createOverlayGeometry();

    if (!inport_.isReady())
        outport_.activateTarget();
    else
        privatePort_.activateTarget();

    glClearDepth(1);

    glEnable(GL_DEPTH_TEST);

    // render cube and axis overlays
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if(enableProp_.get()) {
        renderGeometry();
    }

    if (!inport_.isReady())
        outport_.deactivateTarget();
    else
        privatePort_.deactivateTarget();

    // restore OpenGL state
    glCullFace(GL_BACK);
    glDisable(GL_CULL_FACE);
    LGL_ERROR;

    glEnable(GL_DEPTH_TEST);

    // now do the composition of the OrientationOverlay with the inport render data by shader
    if (inport_.isReady()) {
        outport_.activateTarget();
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // use the shader to draw cube and axis over inport render data
        // therefore bind rgba and depth values of both to textures
        tgt::TextureUnit colorUnit0, depthUnit0, colorUnit1, depthUnit1;
        privatePort_.bindTextures(colorUnit0.getEnum(), depthUnit0.getEnum());
        inport_.bindTextures(colorUnit1.getEnum(), depthUnit1.getEnum());

        // initialize shader
        program_->activate();
        setGlobalShaderParameters(program_);
        program_->setUniform("colorTexMe_", colorUnit0.getUnitNumber());
        program_->setUniform("depthTexMe_", depthUnit0.getUnitNumber());
        program_->setUniform("colorTexIn_", colorUnit1.getUnitNumber());
        program_->setUniform("depthTexIn_", depthUnit1.getUnitNumber());
        privatePort_.setTextureParameters(program_, "textureParametersMe_");
        inport_.setTextureParameters(program_, "textureParametersIn_");

        glDepthFunc(GL_ALWAYS);
        renderQuad(); // render quad primitive textured by fragment shader
        glDepthFunc(GL_LESS);
        outport_.deactivateTarget();

        program_->deactivate();
        glActiveTexture(GL_TEXTURE0);
    }

    LGL_ERROR;
}

//--------------------------------------------------------------------------------------------
//          render functions
//--------------------------------------------------------------------------------------------
void OrientationOverlay::renderGeometry() {
    tgt::Shader* prog = shaderProp_.getShader();
    prog->activate();
    //global settings
    setGlobalShaderParameters(prog, &(cameraProp_.get()));

    //set lightning
    prog->setUniform("enableLighting_", orientationTypeProp_.getValue() != OT_AXES_2D); // no light in 2d
    prog->setUniform("lightPositionEye_", tgt::vec3((shiftXProp_.get()*2.f-1.f),(shiftYProp_.get()*2.f-1.f),1.0f));
    prog->setUniform("shininess_", 5.f);
    prog->setUniform("lightSource_.ambientColor_", tgt::vec3(0.4f, 0.4f, 0.4f));
    prog->setUniform("lightSource_.diffuseColor_", tgt::vec3(0.4f, 0.4f, 0.4f));
    prog->setUniform("lightSource_.specularColor_", tgt::vec3(0.6f, 0.6f, 0.6f));

    // disable clipping
    prog->setUniform("enableClipping_", false);
    prog->setUniform("plane_", tgt::vec4(0.0f, 0.0f, 0.0f, 0.0f));

    // move/rotate axis right
    prog->setUniform("projectionMatrix_", tgt::mat4::createOrtho(-1,1,-1,1,-2,2));

    const tgt::Camera* currentCam = 0; //camera dependent on 2d or 3d
    tgt::Camera cam2D;
    // calculate viewMatrix
    tgt::mat4 viewMatrix = tgt::mat4::createTranslation(tgt::vec3(shiftXProp_.get()*2.0f-1.0f, shiftYProp_.get()*2.0f-1.0f, 0.f));
    viewMatrix *= tgt::mat4::createScale(tgt::vec3((float)outport_.getSize().y / (float)outport_.getSize().x, 1.f, 1.f));
    //is 2d?
    if(orientationTypeProp_.getValue() == OT_AXES_2D) {
        switch(axesAlignmentProp_.getValue()) {
            case XY_PLANE:
                cam2D.setPosition(tgt::vec3(0.f,0.f,-1.f)); cam2D.setFocus(tgt::vec3(0.f,0.f,0.f)); cam2D.setUpVector(tgt::vec3(0.f,-1.f,0.f));
                viewMatrix *= cam2D.getViewMatrix().getRotationalPart();
                break;
            case XZ_PLANE:
                cam2D.setPosition(tgt::vec3(0.f,-1.f,0.f)); cam2D.setFocus(tgt::vec3(0.f,0.f,0.f)); cam2D.setUpVector(tgt::vec3(0.f,0.f,1.f));
                viewMatrix *= cam2D.getViewMatrix().getRotationalPart();
                break;
            case YZ_PLANE:
                cam2D.setPosition(tgt::vec3(-1.f,0.f,0.f)); cam2D.setFocus(tgt::vec3(0.f,0.f,0.f)); cam2D.setUpVector(tgt::vec3(0.f,0.f,1.f));
                viewMatrix *= cam2D.getViewMatrix().getRotationalPart();
                break;
            case UNALIGNED_PLANE:
                LERROR("Unaligned Plane is not supported in 2D overlay mode");
                prog->deactivate();
                return;
        }
        currentCam = &(cam2D);
    } else {
        // use normal camera in 3d
        viewMatrix *= cameraProp_.get().getViewMatrix().getRotationalPart();
        currentCam = &(cameraProp_.get());
    }
    prog->setUniform("viewMatrix_", viewMatrix);

    // set modelMatrix
    prog->setUniform("modelMatrix_", currentGeometry_->getTransformationMatrix());

    // calculate normalMatrix out of view- and modelMatrix
    tgt::mat3 modelViewMat3 = (viewMatrix * currentGeometry_->getTransformationMatrix()).getRotationalPartMat3();
    tgt::mat3 normalMatrix;
    if(!modelViewMat3.invert(normalMatrix)) {
        LWARNING("Could not generate normal matrix out of current view / model matrix, using identity.");
        normalMatrix = tgt::mat3::identity;
    }
    normalMatrix = transpose(normalMatrix);
    prog->setUniform("normalMatrix_", normalMatrix);
    LGL_ERROR;

    tgt::TextureUnit texUnit;
    if(currentGeometry_->supportsTextureData()) {
        texUnit.activate();
        currentGeometry_->getTextureData()->bind();
        prog->setUniform("textures_", texUnit.getUnitNumber());
    }
    LGL_ERROR;

    currentGeometry_->render();

    //render labels if needed
    if(renderAxesLabelsProp_.get() && ((orientationTypeProp_.getValue() == OT_AXES_3D) || (orientationTypeProp_.getValue() == OT_AXES_2D))) {
        if(letterXGeometry_) {
            prog->setUniform("modelMatrix_", tgt::Matrix::rigidBodyTranformation(-currentCam->getLook(),currentCam->getUpVector(),
                                                                                 tgt::vec3(overlayBaseLength_*axisSizeProp_.get()*0.7f,0.f,0.f)));
            letterXGeometry_->render();
        }
        if(letterYGeometry_) {
            prog->setUniform("modelMatrix_", tgt::Matrix::rigidBodyTranformation(-currentCam->getLook(),currentCam->getUpVector(),
                                                                             tgt::vec3(0.f,overlayBaseLength_*axisSizeProp_.get()*0.7f,0.f)));
            letterYGeometry_->render();
        }
        if(letterZGeometry_) {
            prog->setUniform("modelMatrix_", tgt::Matrix::rigidBodyTranformation(-currentCam->getLook(),currentCam->getUpVector(),
                                                                             tgt::vec3(0.f,0.f,overlayBaseLength_*axisSizeProp_.get()*0.7f)));
            letterZGeometry_->render();
        }
    }


    prog->deactivate();
    LGL_ERROR;
}

void OrientationOverlay::createOverlayGeometry() {
    //create axes, if geometry is not present
    if(!geometryMustBeRecreated_)
        return;

    //clean up
    delete letterXGeometry_; letterXGeometry_ = 0;
    delete letterYGeometry_; letterYGeometry_ = 0;
    delete letterZGeometry_; letterZGeometry_ = 0;


    switch(orientationTypeProp_.getValue()) {
        case OT_AXES_2D:
            switch(axesAlignmentProp_.getValue()) {
            case XY_PLANE:
                createAxesGeometry(true,true,false);
                break;
            case XZ_PLANE:
                createAxesGeometry(true,false,true);
                break;
            case YZ_PLANE:
                createAxesGeometry(false,true,true);
                break;
            default:
                LERROR("Unaligned plane not supported in 2D mode");
                geometryMustBeRecreated_ = true;
                return;
            }
            break;
        case OT_AXES_3D:
            createAxesGeometry(true,true,true);
            break;
        case OT_COLOR_CUBE:
            createCubeGeometry(true,false);
            break;
        case OT_TEXTURE_CUBE:
            createCubeGeometry(false,true);
            break;
        case OT_COLOR_TEXTURE_CUBE:
            createCubeGeometry(true,true);
            break;
        default:
            tgtAssert(false,"should not get here");
    }

    //set helper to false
    geometryMustBeRecreated_ = false;
    std::string header = generateHeader();
    header += currentGeometry_->getShaderDefines();
    if(currentGeometry_->supportsTextureData()) {
        header += "#define USE_TEXTURE_ARRAY\n";
        if(currentGeometry_->supportsColors())
            header += "#define TEXTURE_MODE_MODULATE\n";
        else
            header += "#define TEXTURE_MODE_REPLACE\n";
    }
    shaderProp_.setHeader(header);
    shaderProp_.rebuild();
}


TriangleMeshGeometryUInt16IndexedColorNormal* OrientationOverlay::createArrowGeometry(tgt::vec4 col) {
    TriangleMeshGeometryUInt16IndexedColorNormal* geom =  new TriangleMeshGeometryUInt16IndexedColorNormal();

    float smallRadius = 0.005f*axisSizeProp_.get();
    float mainLength = overlayBaseLength_*axisSizeProp_.get();

    geom->addDiskGeometry(0.f,smallRadius, tgt::vec3(0.f,0.f,-0.4f*mainLength), col);
    geom->addCylinderGeometry(smallRadius, smallRadius, 0.75f*mainLength, tgt::vec3(0.f,0.f,-0.4f*mainLength), col);
    geom->addDiskGeometry(0.f,3*smallRadius,tgt::vec3(0.f,0.f,0.35f*mainLength), col);
    geom->addCylinderGeometry(3*smallRadius,0.f,0.25f*mainLength,tgt::vec3(0.f,0.f,0.35f*mainLength), col);

    return geom;
}

void OrientationOverlay::createAxesGeometry(bool createXAxis, bool createYAxis, bool createZAxis) {
    //clear old geom
    delete currentGeometry_;
    TriangleMeshGeometryUInt16IndexedColorNormal* tmp = new TriangleMeshGeometryUInt16IndexedColorNormal();

    //x axis
    if(createXAxis) {
        TriangleMeshGeometryUInt16IndexedColorNormal* xAxis = createArrowGeometry(tgt::vec4(1.f,0.f,0.f,1.f));
        xAxis->setTransformationMatrix(tgt::Matrix::createRotationY(1.570796327f));//pi/2
        tmp->addMesh(xAxis);
        delete xAxis;
        letterXGeometry_ =  new TriangleMeshGeometryUInt16IndexedColorNormal();
        static_cast<TriangleMeshGeometryUInt16IndexedColorNormal*>(letterXGeometry_)->addLetterX(0.03f*axisSizeProp_.get(),tgt::vec4(1.f,0.f,0.f,1.f));
    }
    //y axis
    if(createYAxis) {
        TriangleMeshGeometryUInt16IndexedColorNormal* yAxis = createArrowGeometry(tgt::vec4(0.f,1.f,0.f,1.f));
        yAxis->setTransformationMatrix(tgt::Matrix::createRotationX(-1.570796327f));//pi/2
        tmp->addMesh(yAxis);
        delete yAxis;
        letterYGeometry_ =  new TriangleMeshGeometryUInt16IndexedColorNormal();
        static_cast<TriangleMeshGeometryUInt16IndexedColorNormal*>(letterYGeometry_)->addLetterY(0.03f*axisSizeProp_.get(),tgt::vec4(0.f,1.f,0.f,1.f));
    }
    //z axis
    if(createZAxis) {
        TriangleMeshGeometryUInt16IndexedColorNormal* zAxis = createArrowGeometry(tgt::vec4(0.f,0.f,1.f,1.f));
        tmp->addMesh(zAxis);
        delete zAxis;
        letterZGeometry_ =  new TriangleMeshGeometryUInt16IndexedColorNormal();
        static_cast<TriangleMeshGeometryUInt16IndexedColorNormal*>(letterZGeometry_)->addLetterZ(0.03f*axisSizeProp_.get(),tgt::vec4(0.f,0.f,1.f,1.f));
    }

    currentGeometry_ =  tmp;
}

void OrientationOverlay::createCubeGeometry(bool colored, bool textured) {
      //clear old geom
    delete currentGeometry_;
    if(colored) {
        if(textured) {
            createCubeGeometry<TriangleMeshGeometryUInt16IndexedColorNormalTexCoord>();
        } else{
            createCubeGeometry<TriangleMeshGeometryUInt16IndexedColorNormal>();
        }
    } else {
        if(textured) {
            createCubeGeometry<TriangleMeshGeometryUInt16IndexedNormalTexCoord>();
        } else{
            tgtAssert(false,"not implemented yet!")
        }
    }
}

//--------------------------------------------------------------------------------------------
//          helper functions
//--------------------------------------------------------------------------------------------
void OrientationOverlay::adjustPropertyVisibility() {
    switch(orientationTypeProp_.getValue()) {
       case OT_AXES_2D:
            setPropertyGroupVisible("axis",true);
            setPropertyGroupVisible("cube",false);
            setPropertyGroupVisible("texture",false);
            break;
        case OT_AXES_3D:
            setPropertyGroupVisible("axis",true);
            setPropertyGroupVisible("cube",false);
            setPropertyGroupVisible("texture",false);
            axesAlignmentProp_.setVisibleFlag(false); // not this one in 3d
            break;
        case  OT_COLOR_CUBE:
            setPropertyGroupVisible("axis",false);
            setPropertyGroupVisible("cube",true);
            setPropertyGroupVisible("texture",false);
            break;
        case OT_TEXTURE_CUBE:
            setPropertyGroupVisible("axis",false);
            setPropertyGroupVisible("cube",true);
            setPropertyGroupVisible("texture",true);
            break;
        case OT_COLOR_TEXTURE_CUBE:
            setPropertyGroupVisible("axis",false);
            setPropertyGroupVisible("cube",true);
            setPropertyGroupVisible("texture",true);
            break;
        default:
            tgtAssert(false, "unknown orientation type");
    }
    //re-create geometry
    geometryMustBeRecreated_ = true;
}

void OrientationOverlay::reloadTextures() {
    reloadTextures_ = true;
    invalidate();
}

void OrientationOverlay::invalidateGeometry() {
    geometryMustBeRecreated_ = true;
    invalidate();
}

#ifdef VRN_MODULE_DEVIL
void OrientationOverlay::loadTextures() {

    if (loadingTextures_)
        return;

    if (tgt::TextureManager::isInited()) {

        loadingTextures_ = true;

        // first dispose textures
        TexMgr.dispose(topTex_);
        TexMgr.dispose(bottomTex_);
        TexMgr.dispose(frontTex_);
        TexMgr.dispose(backTex_);
        TexMgr.dispose(leftTex_);
        TexMgr.dispose(rightTex_);
        LGL_ERROR;

        // now try loading textures
        if (!filenameTopProp_.get().empty())
            topTex_ = TexMgr.load(filenameTopProp_.get());
        else if (textureNames_[0] != "") {
            topTex_ = TexMgr.load(textureNames_[0]);
            if (topTex_)
                filenameTopProp_.set(textureNames_[0]);
        }
        else
            topTex_ = 0;

        if (!filenameBottomProp_.get().empty())
            bottomTex_ = TexMgr.load(filenameBottomProp_.get());
        else if (textureNames_[1] != "") {
            bottomTex_ = TexMgr.load(textureNames_[1]);
            if (bottomTex_)
                filenameBottomProp_.set(textureNames_[1]);
        }
        else
            bottomTex_ = 0;

        if (!filenameFrontProp_.get().empty())
            frontTex_ = TexMgr.load(filenameFrontProp_.get());
        else if (textureNames_[2] != "") {
            frontTex_ = TexMgr.load(textureNames_[2]);
            if (frontTex_)
                filenameFrontProp_.set(textureNames_[2]);
        }
        else
            frontTex_ = 0;

        if (!filenameBackProp_.get().empty())
            backTex_ = TexMgr.load(filenameBackProp_.get());
        else if (textureNames_[3] != "") {
            backTex_ = TexMgr.load(textureNames_[3]);
            if (backTex_)
                filenameBackProp_.set(textureNames_[3]);
        }
        else
            backTex_ = 0;

        if (!filenameLeftProp_.get().empty())
            leftTex_ = TexMgr.load(filenameLeftProp_.get());
        else if (textureNames_[4] != "") {
            leftTex_ = TexMgr.load(textureNames_[4]);
            if (leftTex_)
                filenameLeftProp_.set(textureNames_[4]);
        }
        else
            leftTex_ = 0;

        if (!filenameRightProp_.get().empty())
            rightTex_ = TexMgr.load(filenameRightProp_.get());
        else if (textureNames_[5] != "") {
            rightTex_ = TexMgr.load(textureNames_[5]);
            if (rightTex_)
                filenameRightProp_.set(textureNames_[5]);
        }
        else
            rightTex_ = 0;

        LGL_ERROR;
        loadingTextures_ = false;
        reloadTextures_ = false;

        invalidate();
    }
    else {
        LWARNING("loadTextures(): TextureManager not initialized");
    }
}
#else
void OrientationOverlay::loadTextures() {
    LWARNING("OrientationOverlay needs module DeVil to load texture!");
}

#endif //VRN_MODULE_DEVIL

} // namespace
