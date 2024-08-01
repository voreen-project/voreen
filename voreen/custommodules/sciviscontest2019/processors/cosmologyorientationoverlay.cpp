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

#include "cosmologyorientationoverlay.h"
#include "voreen/core/interaction/camerainteractionhandler.h"
#include "voreen/core/voreenapplication.h"

#include "tgt/texturemanager.h"
#include "tgt/gpucapabilities.h"
#include "tgt/glmath.h"
#include "tgt/textureunit.h"
#include "tgt/matrix.h"

#include "voreen/core/datastructures/geometry/trianglemeshgeometryindexed.h"

#include "tgt/immediatemode/immediatemode.h"
#include "tgt/matrixstack.h"

using tgt::TextureUnit;
using tgt::vec3;
using tgt::vec4;

namespace voreen {

const std::string CosmologyOrientationOverlay::loggerCat_("voreen.CosmologyOrientationOverlay");

CosmologyOrientationOverlay::CosmologyOrientationOverlay()
    : ImageProcessor("image/orientationoverlay") //determines .frag shader
    //ports
    , inport_(Port::INPORT, "image.input", "Image Input")
    , outport_(Port::OUTPORT, "image.output", "Image Output")
    , privatePort_(Port::OUTPORT, "image.tmp", "image.tmp", false)
    //basic
    , enableProp_("enableProp","Enable",true)
    , enableProp2_("enableProp2","Enable2",true)
    , orientationTypeProp_("orientationTypeProp","Overlay Type")
    , cameraProp_("camera", "Camera", tgt::Camera(vec3(0.f, 0.f, 3.5f), vec3(0.f, 0.f, 0.f), vec3(0.f, 1.f, 0.f)),
              false,true,500.f,Processor::INVALID_RESULT,Property::LOD_DEBUG)
    //position
    , shiftXProp_("shiftX", "Horizontal Position", 0.85f, 0.0f, 1.0f)
    , shiftYProp_("shiftY", "Vertical Position", 0.15f, 0.0f, 1.0f)
    //axis
    , axisSizeProp_("axisLength", "Axes Length", 1.f, 0.5f, 3.f,Processor::VALID)
    //cube
    , cubeSizeProp_("cubeSize", "Cube Size", 1.f, 0.5f, 3.f)
    , colorCubeRotation_("colorcuberotation", "Color Cube Rotation", tgt::mat4::identity, tgt::mat4(-1.1f), 
            tgt::mat4(1.1f), Processor::INVALID_RESULT, NumericProperty<tgt::mat4>::STATIC, Property::LOD_DEBUG) 
    , shaderProp_("geometry.prg", "Shader", "trianglemesh.frag", "trianglemesh.vert", "trianglemesh.geom",
                  Processor::INVALID_PROGRAM,Property::LOD_DEBUG)
    //helper
    , geometryMustBeRecreated_(true)
    , axisBaseLength_(0.16f)
    , cubeBaseLength_(0.16f)
    , currentGeometry_(0)
{
    //ports
    addPort(inport_);
    addPort(outport_);
    addPrivateRenderPort(&privatePort_);
    //basic
    addProperty(enableProp_);
    addProperty(enableProp2_);
    addProperty(orientationTypeProp_);
        orientationTypeProp_.addOption("diagonals", "Diagonals", OT_DIAGONALS);
        orientationTypeProp_.addOption("diagonalsedges", "Diagonals and Edges", OT_DIAGONALS_EDGES);
        orientationTypeProp_.addOption("diagonalsaxes", "Diagonals and Axes", OT_DIAGONALS_AXES);
        orientationTypeProp_.addOption("axes","Axes",OT_AXES);
        orientationTypeProp_.addOption("col_cube","Colored Cube",OT_COLOR_CUBE);
        orientationTypeProp_.onChange(MemberFunctionCallback<CosmologyOrientationOverlay>(this, &CosmologyOrientationOverlay::adjustPropertyVisibility));
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
    axisSizeProp_.onChange(MemberFunctionCallback<CosmologyOrientationOverlay>(this, &CosmologyOrientationOverlay::invalidateGeometry));
    setPropertyGroupGuiName("axis","Axis Settings");
    //cube
    addProperty(cubeSizeProp_);
    cubeSizeProp_.onChange(MemberFunctionCallback<CosmologyOrientationOverlay>(this, &CosmologyOrientationOverlay::invalidateGeometry));
    cubeSizeProp_.setGroupID("cube");
    setPropertyGroupGuiName("cube","Cube Settings");

    colorCubeRotation_.setReadOnlyFlag(true);
    colorCubeRotation_.onChange(MemberFunctionCallback<CosmologyOrientationOverlay>(this, &CosmologyOrientationOverlay::rotationChanged));

    addProperty(colorCubeRotation_);

    addProperty(shaderProp_);


    adjustPropertyVisibility();
}

CosmologyOrientationOverlay::~CosmologyOrientationOverlay() {
    //TODO: delete geometry ?!
}

Processor* CosmologyOrientationOverlay::create() const {
    return new CosmologyOrientationOverlay();
}

void CosmologyOrientationOverlay::initialize() {
    ImageProcessor::initialize();
    //loadTextures();
}

void CosmologyOrientationOverlay::deinitialize() {
    //delete geom
    delete currentGeometry_;

    ImageProcessor::deinitialize();
}

bool CosmologyOrientationOverlay::isReady() const {
    return outport_.isReady();
}

//--------------------------------------------------------------------------------------------
//          process functions
//--------------------------------------------------------------------------------------------
void CosmologyOrientationOverlay::beforeProcess() {
    ImageProcessor::beforeProcess();
    //if (reloadTextures_)
    //    loadTextures();
    if(!shaderProp_.hasValidShader() || getInvalidationLevel() >= Processor::INVALID_PROGRAM) {
        shaderProp_.setHeader(generateHeader() + (currentGeometry_ ? currentGeometry_->getShaderDefines(): ""));
        shaderProp_.rebuild();
    }
}

void CosmologyOrientationOverlay::process() {
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

    if(enableProp_.get() || enableProp2_.get()) {
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

    // now do the composition of the CosmologyOrientationOverlay with the inport render data by shader
    if (inport_.isReady()) {
        outport_.activateTarget();
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // use the shader to draw cube and axis over inport render data
        // therefore bind rgba and depth values of both to textures
        TextureUnit colorUnit0, depthUnit0, colorUnit1, depthUnit1;
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
void CosmologyOrientationOverlay::renderGeometry() {

    tgt::Shader* prog = shaderProp_.getShader();
    prog->activate();
    //global settings
    setGlobalShaderParameters(prog, &(cameraProp_.get()));
    prog->setUniform("enableLighting_", true);
    prog->setUniform("lightPositionEye_", tgt::vec3((shiftXProp_.get()*2.f-1.f),(shiftYProp_.get()*2.f-1.f),1.0f));
    prog->setUniform("shininess_", 5.f);

    prog->setUniform("lightSource_.ambientColor_", tgt::vec3(0.4f, 0.4f, 0.4f));
    prog->setUniform("lightSource_.diffuseColor_", tgt::vec3(0.4f, 0.4f, 0.4f));
    prog->setUniform("lightSource_.specularColor_", tgt::vec3(0.6f, 0.6f, 0.6f));

    // draw the cube in the desired position

    // move/rotate axis right
    prog->setUniform("projectionMatrix_", tgt::mat4::createOrtho(-1,1,-1,1,-2,2));

    // calculate viewMatrix
    tgt::mat4 viewMatrix = tgt::mat4::createTranslation(tgt::vec3(shiftXProp_.get()*2.0f-1.0f, shiftYProp_.get()*2.0f-1.0f, 0.f));
    viewMatrix *= tgt::mat4::createScale(tgt::vec3((float)outport_.getSize().y / (float)outport_.getSize().x, 1.f, 1.f));
    viewMatrix *= cameraProp_.get().getViewMatrix().getRotationalPart();
    prog->setUniform("viewMatrix_", viewMatrix);

    // set modelMatrix
    prog->setUniform("modelMatrix_", currentGeometry_->getTransformationMatrix());

    // calculate normalMatrix out of view- and modelMatrix
    tgt::mat3 modelViewMat3 = (viewMatrix * currentGeometry_->getTransformationMatrix()).getRotationalPartMat3();
    //tgt::mat3 modelViewMat3 = (cameraProp_.get().getViewMatrix() * currentGeometry_->getTransformationMatrix()).getRotationalPartMat3();
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

    prog->deactivate();

    //TODO: use a different geometry
    if (orientationTypeProp_.getValue() == OT_DIAGONALS_EDGES) {
        //render the axes
        tgt::ImmediateMode& im = tgt::ImmediateMode::getRef();
        tgt::MatrixStack& ms = tgt::MatrixStack::getRef();

        ms.matrixMode(tgt::MatrixStack::MODELVIEW);
        ms.pushMatrix();
        ms.loadMatrix(viewMatrix * currentGeometry_->getTransformationMatrix());

        ms.matrixMode(tgt::MatrixStack::PROJECTION);
        ms.pushMatrix();
        ms.loadMatrix(tgt::mat4::createOrtho(-1,1,-1,1,-2,2));
        
        im.begin(tgt::ImmediateMode::LINES);
        
        for (int sign = -1; sign <= 1; sign += 2) {      
            for (int dim = -1; dim < 3; ++dim) {
                vec3 direction = normalize(static_cast<float>(sign) * vec3(1.f));
                vec3 startColor = static_cast<float>(sign) * vec3(0.5f);
                if (dim >= 0) {
                    size_t index = static_cast<size_t>(dim);
                    direction.elem[index] = -direction.elem[index];
                    startColor.elem[index] = -startColor.elem[index];
                }
                startColor += vec3(0.5f);

                vec3 startPos = 0.85f * direction * (axisBaseLength_ * axisSizeProp_.get());
            
                //find line end positions
                for (size_t i = 0; i < 3; ++i) {
                    if (direction.elem[i] < 0) {
                        vec3 endPos = direction;
                        endPos.elem[i] = -endPos.elem[i];
                        vec3 endColor = endPos * vec3(0.5f) + vec3(0.5f);

                        endPos = 0.85f * endPos * (axisBaseLength_ * axisSizeProp_.get());
                    
                        im.color(startColor);
                        im.vertex(startPos);

                        im.color(endColor);
                        im.vertex(endPos);
                    }
                } 
            }
        }

        im.end();

        ms.matrixMode(tgt::MatrixStack::PROJECTION);
        ms.popMatrix();
        ms.matrixMode(tgt::MatrixStack::MODELVIEW);
        ms.popMatrix();
    }

    LGL_ERROR;
}

void CosmologyOrientationOverlay::createOverlayGeometry() {
    //create axes, if geometry is not present
    if(!geometryMustBeRecreated_)
        return;

    switch(orientationTypeProp_.getValue()) {
        case OT_DIAGONALS_EDGES:
        case OT_DIAGONALS:
            createDiagonalsGeometry();
            break;
        case OT_AXES:
            createAxesGeometry();
            break;
        case OT_DIAGONALS_AXES:
            createDiagonalsAxesGeometry();
            break;
        case OT_COLOR_CUBE:
            createCubeGeometry();
            break;
        default:
            tgtAssert(false,"should not get here");
    }
    //set helper to false
    geometryMustBeRecreated_ = false;
    std::string header = generateHeader();
    header += currentGeometry_->getShaderDefines();

    shaderProp_.setHeader(header);
    shaderProp_.rebuild();
}


TriangleMeshGeometryUInt16IndexedColorNormal* CosmologyOrientationOverlay::createArrowGeometry(tgt::vec4 col) {
    TriangleMeshGeometryUInt16IndexedColorNormal* geom =  new TriangleMeshGeometryUInt16IndexedColorNormal();

    float smallRadius = 0.005f*axisSizeProp_.get();
    float mainLength = axisBaseLength_*axisSizeProp_.get();

    geom->addDiskGeometry(0.f,smallRadius, tgt::vec3(0.f,0.f,-0.4f*mainLength), col);
    geom->addCylinderGeometry(smallRadius, smallRadius, 0.75f*mainLength, tgt::vec3(0.f,0.f,-0.4f*mainLength), col);
    geom->addDiskGeometry(0.f,3*smallRadius,tgt::vec3(0.f,0.f,0.35f*mainLength), col);
    geom->addCylinderGeometry(3*smallRadius,0.f,0.25f*mainLength,tgt::vec3(0.f,0.f,0.35f*mainLength), col);

    return geom;
}

void CosmologyOrientationOverlay::createDiagonalsGeometry() {
    //clear old geom
    delete currentGeometry_;

    //create new geometry
    TriangleMeshGeometryUInt16IndexedColorNormal* mesh = new TriangleMeshGeometryUInt16IndexedColorNormal();
    for (int sign = -1; sign <= 1; sign += 2) {      
        for (int dim = -1; dim < 3; ++dim) {
            vec3 direction = normalize(static_cast<float>(sign) * vec3(1.f));
            vec3 color = static_cast<float>(sign) * vec3(0.5f);
            if (dim >= 0) {
                size_t index = static_cast<size_t>(dim);
                direction.elem[index] = -direction.elem[index];
                color.elem[index] = -color.elem[index];
            }
            color += vec3(0.5f);
            
            TriangleMeshGeometryUInt16IndexedColorNormal* axis = createArrowGeometry(vec4(color,1.f));
            vec3 rotationAxis = tgt::cross(vec3(0.f,0.f,1.f), direction);
            float rotationAngle = acos(tgt::dot(direction, vec3(0.f,0.f,1.f)));
            axis->setTransformationMatrix(tgt::Matrix::createRotation(rotationAngle, rotationAxis) * tgt::Matrix::createScale(vec3(0.85f)) 
                                          * tgt::Matrix::createTranslation(vec3(0.f, 0.f, axisBaseLength_ * axisSizeProp_.get() * 0.4f)));
            mesh->addMesh(axis);
            delete axis;
        }
    }

    // transform by the additional rotation matrix
    mesh->transform(tgt::mat4(colorCubeRotation_.get()));

    currentGeometry_ =  mesh;
}

void CosmologyOrientationOverlay::createAxesGeometry() {
    //clear old geom
    delete currentGeometry_;

    //create new geometry
    TriangleMeshGeometryUInt16IndexedColorNormal* mesh = new TriangleMeshGeometryUInt16IndexedColorNormal();
    for (int sign = -1; sign <= 1; sign += 2) {      
        for (int dim = 0; dim < 3; ++dim) {
            vec3 direction = vec3(0.f);
            direction.elem[dim] = static_cast<float>(sign) * 1.f;
            vec3 color = direction * vec3(0.5f);
            color += vec3(0.5f);
            
            TriangleMeshGeometryUInt16IndexedColorNormal* axis = createArrowGeometry(vec4(color,1.f));
            vec3 rotationAxis = tgt::cross(vec3(0.f,0.f,1.f), direction);
            if (tgt::length(rotationAxis) < 0.001f)
                rotationAxis = vec3(1.f, 0.f, 0.f);
            float rotationAngle = acos(tgt::dot(direction, vec3(0.f,0.f,1.f)));
            if (rotationAngle > 0.001f)
                axis->setTransformationMatrix(tgt::Matrix::createRotation(rotationAngle, rotationAxis) * tgt::Matrix::createScale(vec3(0.85f)) 
                                          * tgt::Matrix::createTranslation(vec3(0.f, 0.f, axisBaseLength_ * axisSizeProp_.get() * 0.4f)));
            else 
                axis->setTransformationMatrix(tgt::Matrix::createScale(vec3(0.85f)) * tgt::Matrix::createTranslation(vec3(0.f, 0.f, axisBaseLength_ * axisSizeProp_.get() * 0.4f)));

            mesh->addMesh(axis);
            delete axis;
        }
    }

    // transform by the additional rotation matrix
    mesh->transform(colorCubeRotation_.get());

    currentGeometry_ =  mesh;
}

void CosmologyOrientationOverlay::createDiagonalsAxesGeometry() {
    //clear old geom
    delete currentGeometry_;

    //create new geometry
    TriangleMeshGeometryUInt16IndexedColorNormal* mesh = new TriangleMeshGeometryUInt16IndexedColorNormal();
    for (int sign = -1; sign <= 1; sign += 2) {      
        for (int dim = 0; dim < 3; ++dim) {
            vec3 direction = vec3(0.f);
            direction.elem[dim] = static_cast<float>(sign) * 1.f;
            vec3 color = direction * vec3(0.5f);
            color += vec3(0.5f);
            
            TriangleMeshGeometryUInt16IndexedColorNormal* axis = createArrowGeometry(vec4(color,1.f));
            vec3 rotationAxis = tgt::cross(vec3(0.f,0.f,1.f), direction);
            if (tgt::length(rotationAxis) < 0.001f)
                rotationAxis = vec3(1.f, 0.f, 0.f);
            float rotationAngle = acos(tgt::dot(direction, vec3(0.f,0.f,1.f)));
            if (rotationAngle > 0.001f)
                axis->setTransformationMatrix(tgt::Matrix::createRotation(rotationAngle, rotationAxis) * tgt::Matrix::createScale(vec3(0.85f)) 
                                          * tgt::Matrix::createTranslation(vec3(0.f, 0.f, axisBaseLength_ * axisSizeProp_.get() * 0.4f)));
            else 
                axis->setTransformationMatrix(tgt::Matrix::createScale(vec3(0.85f)) * tgt::Matrix::createTranslation(vec3(0.f, 0.f, axisBaseLength_ * axisSizeProp_.get() * 0.4f)));

            mesh->addMesh(axis);
            delete axis;
        }
    }

    for (int sign = -1; sign <= 1; sign += 2) {      
        for (int dim = -1; dim < 3; ++dim) {
            vec3 direction = normalize(static_cast<float>(sign) * vec3(1.f));
            vec3 color = static_cast<float>(sign) * vec3(0.5f);
            if (dim >= 0) {
                size_t index = static_cast<size_t>(dim);
                direction.elem[index] = -direction.elem[index];
                color.elem[index] = -color.elem[index];
            }
            color += vec3(0.5f);
            
            TriangleMeshGeometryUInt16IndexedColorNormal* axis = createArrowGeometry(vec4(color,1.f));
            vec3 rotationAxis = tgt::cross(vec3(0.f,0.f,1.f), direction);
            float rotationAngle = acos(tgt::dot(direction, vec3(0.f,0.f,1.f)));
            axis->setTransformationMatrix(tgt::Matrix::createRotation(rotationAngle, rotationAxis) * tgt::Matrix::createScale(vec3(0.85f)) 
                                          * tgt::Matrix::createTranslation(vec3(0.f, 0.f, axisBaseLength_ * axisSizeProp_.get() * 0.4f)));
            mesh->addMesh(axis);
            delete axis;
        }
    }

    // transform by the additional rotation matrix
    mesh->transform(colorCubeRotation_.get());

    currentGeometry_ =  mesh;

}

void CosmologyOrientationOverlay::createCubeGeometry(/*bool colored, bool textured*/) {
    //clear old geom
    delete currentGeometry_;

    // create new geometry
    TriangleMeshGeometry<VertexColorNormal>* geom = new TriangleMeshGeometryColorNormal();
   
    for (int sign = -1; sign <= 1; sign += 2) {
        for (size_t dim = 0; dim < 3; ++dim) {
            vec3 normal = vec3(0.f);
            normal[dim] = 1.f * static_cast<float>(sign);

            VertexColorNormal vertices[4];
            
            //set fixed coordinate according to normal and sign
            for (size_t i = 0; i < 4; ++i) 
                vertices[i].pos_.elem[dim] = (sign < 0) ? -0.5 : 0.5;

            //set the other coordinates
            for (size_t incr = 1; incr <= 2; ++incr) {
                size_t curDim = (dim + incr) % 3;

                vertices[0].pos_[curDim] = -0.5f;
                vertices[1].pos_[curDim] = (incr == 1) ? 0.5f : -0.5f;
                vertices[2].pos_[curDim] = 0.5f;
                vertices[3].pos_[curDim] = (incr == 1) ? -0.5f : 0.5f;
            }

            //set normals and colors
            for (size_t i = 0; i < 4; ++i) {
                vertices[i].color_ = vec4(vertices[i].pos_ + vec3(0.5f), 1.f);
                vertices[i].normal_ = normal;
            }

            if (sign > 0) 
                geom->addQuad(vertices[0], vertices[1], vertices[2], vertices[3]);
            else 
                geom->addQuad(vertices[0], vertices[3], vertices[2], vertices[1]);
        }
    }

    geom->setTransformationMatrix(tgt::Matrix::createScale(vec3(cubeBaseLength_ * cubeSizeProp_.get())));  

    // transform by the additional rotation matrix
    geom->transform(colorCubeRotation_.get());

    currentGeometry_ = geom;
}

//--------------------------------------------------------------------------------------------
//          helper functions
//--------------------------------------------------------------------------------------------
void CosmologyOrientationOverlay::adjustPropertyVisibility() {
    switch(orientationTypeProp_.getValue()) {
        case OT_DIAGONALS_EDGES:
        case OT_DIAGONALS:
        case OT_DIAGONALS_AXES:    
        case OT_AXES:
            setPropertyGroupVisible("axis",true);
            setPropertyGroupVisible("cube",false);
            break;
        case  OT_COLOR_CUBE:
            setPropertyGroupVisible("axis",false);
            setPropertyGroupVisible("cube",true);
            break;
        default:
            tgtAssert(false, "unknown orientation type");
    }
    //re-create geometry
    geometryMustBeRecreated_ = true;
}

void CosmologyOrientationOverlay::invalidateGeometry() {
    geometryMustBeRecreated_ = true;
    invalidate();
}

void CosmologyOrientationOverlay::rotationChanged() {
    invalidateGeometry();
}


} // namespace
