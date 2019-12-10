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

#include "flowdirectionoverlay.h"

#include "tgt/glmath.h"
#include "tgt/matrix.h"

#include "voreen/core/datastructures/geometry/trianglemeshgeometryindexed.h"

#include "tgt/immediatemode/immediatemode.h"
#include "tgt/matrixstack.h"
#include "tgt/textureunit.h"

namespace voreen {

const std::string FlowDirectionOverlay::loggerCat_("voreen.FlowDirectionOverlay");

FlowDirectionOverlay::FlowDirectionOverlay()
    : ImageProcessor("image/orientationoverlay") //determines .frag shader
    //ports
    , inport_(Port::INPORT, "image.input", "Image Input")
    , outport_(Port::OUTPORT, "image.output", "Image Output")
    , privatePort_(Port::OUTPORT, "image.tmp", "image.tmp", false)
    //basic
    , enableProp_("enableProp","Enable",true)
    , directionTypeProp_("orientationTypeProp","Overlay Type")
    , cameraProp_("camera", "Camera", tgt::Camera(tgt::vec3(0.f, 0.f, 3.5f), tgt::vec3(0.f, 0.f, 0.f), tgt::vec3(0.f, 1.f, 0.f)),
              false,true,500.f,Processor::INVALID_RESULT,Property::LOD_DEBUG)
    , shaderProp_("geometry.prg", "Shader", "trianglemesh.frag", "trianglemesh.vert", "trianglemesh.geom",
                  Processor::INVALID_PROGRAM,Property::LOD_DEBUG)
    // additional visualization properties
    , sphereScaleFactor_("spherescalefactor", "Sphere Scale Factor", 0.5f, 0.1f, 1.f)
    , edgeLineWidth_("edgelinewidth", "Edge Line Width", 1.f, 1.f, 10.f)
    //position
    , shiftXProp_("shiftX", "Horizontal Position", 0.85f, 0.0f, 1.0f)
    , shiftYProp_("shiftY", "Vertical Position", 0.15f, 0.0f, 1.0f)
    , overlaySizeProp_("overlaySize", "Overlay Size", 1.f, 0.5f, 3.f,Processor::VALID)
    //rotation
    , colorRotationMatrix_("colorRotationMatrix", " To be linked Rotation", tgt::mat4::identity, tgt::mat4(-1.1f),
            tgt::mat4(1.1f), Processor::INVALID_RESULT, NumericProperty<tgt::mat4>::STATIC, Property::LOD_DEBUG)
    //helper
    , geometryMustBeRecreated_(true)
    , overlayBaseSize_(0.16f) //magic number
    , currentGeometry_(0)
    , sphereGeometry_(0)
{

    removeProperty(ImageProcessor::shaderProp_);

    //ports
    addPort(inport_);
    addPort(outport_);
    addPrivateRenderPort(&privatePort_);
    //basic
    addProperty(enableProp_);
    addProperty(directionTypeProp_);
        directionTypeProp_.addOption("diagonals", "Diagonals", FD_DIAGONALS);
        directionTypeProp_.addOption("diagonalsedges", "Diagonals and Edges", FD_DIAGONALS_EDGES);
        directionTypeProp_.addOption("diagonalsaxes", "Diagonals and Axes", FD_DIAGONALS_AXES);
        directionTypeProp_.addOption("axes","Axes",FD_AXES);
        directionTypeProp_.addOption("col_cube","Colored Cube",FD_COLOR_CUBE);
        directionTypeProp_.addOption("col_sphere", "Colored Sphere", FD_COLOR_SPHERE);
        directionTypeProp_.addOption("col_sphere_diag", "Colored Sphere and Diagonals", FD_COLOR_SPHERE_DIAGONALS);
        directionTypeProp_.addOption("col_sphere_diag_edges", "Colored Sphere, Diagonals, and Edges", FD_COLOR_SPHERE_DIAGONALS_EDGES);
        directionTypeProp_.onChange(MemberFunctionCallback<FlowDirectionOverlay>(this, &FlowDirectionOverlay::invalidateGeometry));
        directionTypeProp_.onChange(MemberFunctionCallback<FlowDirectionOverlay>(this, &FlowDirectionOverlay::adjustPropertyVisibility));
    //position
    addProperty(shiftXProp_);
        shiftXProp_.setGroupID("position");
    addProperty(shiftYProp_);
        shiftYProp_.setGroupID("position");
    addProperty(overlaySizeProp_);
        overlaySizeProp_.setGroupID("position");
        overlaySizeProp_.onChange(MemberFunctionCallback<FlowDirectionOverlay>(this, &FlowDirectionOverlay::invalidateGeometry));
    setPropertyGroupGuiName("position","Overlay Settings");
    addProperty(colorRotationMatrix_);
        colorRotationMatrix_.setReadOnlyFlag(true);
        colorRotationMatrix_.onChange(MemberFunctionCallback<FlowDirectionOverlay>(this, &FlowDirectionOverlay::invalidateGeometry));
        colorRotationMatrix_.setGuiName("rotation");
    setPropertyGroupGuiName("rotation","Rotation Settings");

    addProperty(sphereScaleFactor_);
    addProperty(edgeLineWidth_);
    //basic2
    addProperty(cameraProp_);
    addProperty(shaderProp_);
    addProperty(ImageProcessor::shaderProp_);

    adjustPropertyVisibility();
}

FlowDirectionOverlay::~FlowDirectionOverlay() {
    //TODO: delete geometry ?!
}

void FlowDirectionOverlay::initialize() {
    ImageProcessor::initialize();
}

void FlowDirectionOverlay::deinitialize() {
    //delete geom (GL Object)
    delete currentGeometry_;
    delete sphereGeometry_;
    ImageProcessor::deinitialize();
}

bool FlowDirectionOverlay::isReady() const {
    return outport_.isReady();
}

//--------------------------------------------------------------------------------------------
//          process functions
//--------------------------------------------------------------------------------------------
void FlowDirectionOverlay::beforeProcess() {
    ImageProcessor::beforeProcess();
    if(!shaderProp_.hasValidShader() || getInvalidationLevel() >= Processor::INVALID_PROGRAM) {
        // check the geometry type and set the shader defines
        std::string shaderDefines = "";
        if (directionTypeProp_.getValue() == FD_COLOR_SPHERE) {
            if (sphereGeometry_)
                shaderDefines = sphereGeometry_->getShaderDefines();
        }
        else if (currentGeometry_) {
            shaderDefines = currentGeometry_->getShaderDefines();
        }

        shaderProp_.setHeader(generateHeader() + shaderDefines);
        shaderProp_.rebuild();
    }
}

void FlowDirectionOverlay::process() {
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

    // now do the composition of the FlowDirectionOverlay with the inport render data by shader
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
//          Geometry Handling
//--------------------------------------------------------------------------------------------
void FlowDirectionOverlay::renderGeometry() {

    tgt::Shader* prog = shaderProp_.getShader();
    prog->activate();
    //global settings
    setGlobalShaderParameters(prog, &(cameraProp_.get()));
    // enable lighting for all geometries but colored cube and sphere
    bool enableShaderLighting = !(directionTypeProp_.getValue() == FD_COLOR_CUBE);
    prog->setUniform("enableLighting_", enableShaderLighting);
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

    if (directionTypeProp_.getValue() != FD_COLOR_SPHERE) {
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

        // not necessary since all of our data does not use textures
        /*tgt::TextureUnit texUnit;
        if(currentGeometry_->supportsTextureData()) {
            texUnit.activate();
            currentGeometry_->getTextureData()->bind();
            prog->setUniform("textures_", texUnit.getUnitNumber());
        }*/
        LGL_ERROR;

        currentGeometry_->render();
    }

    if ((directionTypeProp_.getValue() == FD_COLOR_SPHERE) || (directionTypeProp_.getValue() == FD_COLOR_SPHERE_DIAGONALS) || (directionTypeProp_.getValue() == FD_COLOR_SPHERE_DIAGONALS_EDGES)) {
        prog->setUniform("enableLighting_", false);
        // set modelMatrix
        tgt::mat4 geometryTransformation = sphereGeometry_->getTransformationMatrix();
        // do an additional scaling when rendering the sphere with diagonals
        if ((directionTypeProp_.getValue() == FD_COLOR_SPHERE_DIAGONALS) || (directionTypeProp_.getValue() == FD_COLOR_SPHERE_DIAGONALS_EDGES))
            geometryTransformation = geometryTransformation * tgt::mat4::createScale(tgt::vec3(sphereScaleFactor_.get()));
        prog->setUniform("modelMatrix_", geometryTransformation);

        // calculate normalMatrix out of view- and modelMatrix
        tgt::mat3 modelViewMat3 = (viewMatrix * geometryTransformation).getRotationalPartMat3();
        tgt::mat3 normalMatrix;
        if(!modelViewMat3.invert(normalMatrix)) {
            LWARNING("Could not generate normal matrix out of current view / model matrix, using identity.");
            normalMatrix = tgt::mat3::identity;
        }
        normalMatrix = transpose(normalMatrix);
        prog->setUniform("normalMatrix_", normalMatrix);
        LGL_ERROR;

        sphereGeometry_->render();
    }

    prog->deactivate();

    //TODO: use a different geometry
    if ((directionTypeProp_.getValue() == FD_DIAGONALS_EDGES) || (directionTypeProp_.getValue() == FD_COLOR_SPHERE_DIAGONALS_EDGES)) {
        //render the axes
        tgt::ImmediateMode& im = tgt::ImmediateMode::getRef();
        tgt::MatrixStack& ms = tgt::MatrixStack::getRef();

        ms.matrixMode(tgt::MatrixStack::MODELVIEW);
        ms.pushMatrix();
        ms.loadMatrix(viewMatrix * currentGeometry_->getTransformationMatrix());

        ms.matrixMode(tgt::MatrixStack::PROJECTION);
        ms.pushMatrix();
        ms.loadMatrix(tgt::mat4::createOrtho(-1,1,-1,1,-2,2));

        glLineWidth(edgeLineWidth_.get());

        im.begin(tgt::ImmediateMode::LINES);

        for (int sign = -1; sign <= 1; sign += 2) {
            for (int dim = -1; dim < 3; ++dim) {
                tgt::vec3 direction = tgt::normalize(static_cast<float>(sign) * tgt::vec3(1.f));
                tgt::vec3 startColor = static_cast<float>(sign) * tgt::vec3(0.5f);
                if (dim >= 0) {
                    size_t index = static_cast<size_t>(dim);
                    direction.elem[index] = -direction.elem[index];
                    startColor.elem[index] = -startColor.elem[index];
                }
                startColor += tgt::vec3(0.5f);

                tgt::vec3 startPos = 0.85f * direction * (overlayBaseSize_ * overlaySizeProp_.get());

                //find line end positions
                for (size_t i = 0; i < 3; ++i) {
                    if (direction.elem[i] < 0) {
                        tgt::vec3 endPos = direction;
                        endPos.elem[i] = -endPos.elem[i];
                        tgt::vec3 endColor = endPos * tgt::vec3(0.5f) + tgt::vec3(0.5f);

                        endPos = 0.85f * endPos * (overlayBaseSize_ * overlaySizeProp_.get());

                        im.color(startColor);
                        im.vertex(startPos);

                        im.color(endColor);
                        im.vertex(endPos);
                    }
                }
            }
        }

        im.end();

        im.color(tgt::vec4::one);

        ms.matrixMode(tgt::MatrixStack::PROJECTION);
        ms.popMatrix();
        ms.matrixMode(tgt::MatrixStack::MODELVIEW);
        ms.popMatrix();

        glLineWidth(1.f);
    }

    LGL_ERROR;
}

void FlowDirectionOverlay::adjustPropertyVisibility() {
    if (directionTypeProp_.getValue() == FD_DIAGONALS_EDGES ||
            directionTypeProp_.getValue() == FD_COLOR_SPHERE_DIAGONALS_EDGES) {
        edgeLineWidth_.setVisibleFlag(true);
    }
    else
        edgeLineWidth_.setVisibleFlag(false);

    if (directionTypeProp_.getValue() == FD_COLOR_SPHERE_DIAGONALS ||
            directionTypeProp_.getValue() == FD_COLOR_SPHERE_DIAGONALS_EDGES) {
        sphereScaleFactor_.setVisibleFlag(true);
    }
    else
        sphereScaleFactor_.setVisibleFlag(false);
}

void FlowDirectionOverlay::createOverlayGeometry() {
    //create axes, if geometry is not present
    if(!geometryMustBeRecreated_)
        return;

    switch(directionTypeProp_.getValue()) {
        case FD_DIAGONALS_EDGES:
        case FD_DIAGONALS:
            createDiagonalsGeometry();
            break;
        case FD_AXES:
            createAxesGeometry();
            break;
        case FD_DIAGONALS_AXES:
            createDiagonalsAxesGeometry();
            break;
        case FD_COLOR_CUBE:
            createCubeGeometry();
            break;
        case FD_COLOR_SPHERE:
            createSphereGeometry();
            break;
        case FD_COLOR_SPHERE_DIAGONALS:
        case FD_COLOR_SPHERE_DIAGONALS_EDGES:
            createSphereGeometry();
            createDiagonalsGeometry();
            break;
        default:
            tgtAssert(false,"should not get here");
    }
    //set helper to false
    geometryMustBeRecreated_ = false;
    std::string header = generateHeader();
    // check the geometry type and set the shader defines
    std::string shaderDefines = "";
    if (directionTypeProp_.getValue() == FD_COLOR_SPHERE) {
        if (sphereGeometry_)
            shaderDefines = sphereGeometry_->getShaderDefines();
    }
    else if (currentGeometry_) {
        shaderDefines = currentGeometry_->getShaderDefines();
    }
    header += shaderDefines;

    shaderProp_.setHeader(header);
    shaderProp_.rebuild();
}


TriangleMeshGeometryUInt16IndexedColorNormal* FlowDirectionOverlay::createArrowGeometry(tgt::vec4 col) {
    TriangleMeshGeometryUInt16IndexedColorNormal* geom =  new TriangleMeshGeometryUInt16IndexedColorNormal();

    float smallRadius = 0.005f*overlaySizeProp_.get();
    float mainLength = overlayBaseSize_* overlaySizeProp_.get();

    geom->addDiskGeometry(0.f,smallRadius, tgt::vec3(0.f,0.f,-0.4f*mainLength), col);
    geom->addCylinderGeometry(smallRadius, smallRadius, 0.75f*mainLength, tgt::vec3(0.f,0.f,-0.4f*mainLength), col);
    geom->addDiskGeometry(0.f,3*smallRadius,tgt::vec3(0.f,0.f,0.35f*mainLength), col);
    geom->addCylinderGeometry(3*smallRadius,0.f,0.25f*mainLength,tgt::vec3(0.f,0.f,0.35f*mainLength), col);

    return geom;
}

void FlowDirectionOverlay::createDiagonalsGeometry() {
    //clear old geom
    delete currentGeometry_;

    //create new geometry
    TriangleMeshGeometryUInt16IndexedColorNormal* mesh = new TriangleMeshGeometryUInt16IndexedColorNormal();
    for (int sign = -1; sign <= 1; sign += 2) {
        for (int dim = -1; dim < 3; ++dim) {
            tgt::vec3 direction = tgt::normalize(static_cast<float>(sign) * tgt::vec3(1.f));
            tgt::vec3 color = static_cast<float>(sign) * tgt::vec3(0.5f);
            if (dim >= 0) {
                size_t index = static_cast<size_t>(dim);
                direction.elem[index] = -direction.elem[index];
                color.elem[index] = -color.elem[index];
            }
            color += tgt::vec3(0.5f);

            TriangleMeshGeometryUInt16IndexedColorNormal* axis = createArrowGeometry(tgt::vec4(color,1.f));
            tgt::vec3 rotationAxis = tgt::cross(tgt::vec3(0.f,0.f,1.f), direction);
            float rotationAngle = acos(tgt::dot(direction, tgt::vec3(0.f,0.f,1.f)));
            axis->setTransformationMatrix(tgt::Matrix::createRotation(rotationAngle, rotationAxis) * tgt::Matrix::createScale(tgt::vec3(0.85f))
                                          * tgt::Matrix::createTranslation(tgt::vec3(0.f, 0.f, overlayBaseSize_ * overlaySizeProp_.get() * 0.4f)));
            mesh->addMesh(axis);
            delete axis;
        }
    }

    // transform by the additional rotation matrix
    mesh->transform(tgt::mat4(colorRotationMatrix_.get()));

    currentGeometry_ =  mesh;
}

void FlowDirectionOverlay::createAxesGeometry() {
    //clear old geom
    delete currentGeometry_;

    //create new geometry
    TriangleMeshGeometryUInt16IndexedColorNormal* mesh = new TriangleMeshGeometryUInt16IndexedColorNormal();
    for (int sign = -1; sign <= 1; sign += 2) {
        for (int dim = 0; dim < 3; ++dim) {
            tgt::vec3 direction = tgt::vec3(0.f);
            direction.elem[dim] = static_cast<float>(sign) * 1.f;
            tgt::vec3 color = direction * tgt::vec3(0.5f);
            color += tgt::vec3(0.5f);

            TriangleMeshGeometryUInt16IndexedColorNormal* axis = createArrowGeometry(tgt::vec4(color,1.f));
            tgt::vec3 rotationAxis = tgt::cross(tgt::vec3(0.f,0.f,1.f), direction);
            if (tgt::length(rotationAxis) < 0.001f)
                rotationAxis = tgt::vec3(1.f, 0.f, 0.f);
            float rotationAngle = acos(tgt::dot(direction, tgt::vec3(0.f,0.f,1.f)));
            if (rotationAngle > 0.001f)
                axis->setTransformationMatrix(tgt::Matrix::createRotation(rotationAngle, rotationAxis) * tgt::Matrix::createScale(tgt::vec3(0.85f))
                                          * tgt::Matrix::createTranslation(tgt::vec3(0.f, 0.f, overlayBaseSize_ * overlaySizeProp_.get() * 0.4f)));
            else
                axis->setTransformationMatrix(tgt::Matrix::createScale(tgt::vec3(0.85f)) *
                                            tgt::Matrix::createTranslation(tgt::vec3(0.f, 0.f, overlayBaseSize_ * overlaySizeProp_.get() * 0.4f)));

            mesh->addMesh(axis);
            delete axis;
        }
    }

    // transform by the additional rotation matrix
    mesh->transform(colorRotationMatrix_.get());

    currentGeometry_ =  mesh;
}

void FlowDirectionOverlay::createDiagonalsAxesGeometry() {
    //clear old geom
    delete currentGeometry_;

    //create new geometry
    TriangleMeshGeometryUInt16IndexedColorNormal* mesh = new TriangleMeshGeometryUInt16IndexedColorNormal();
    for (int sign = -1; sign <= 1; sign += 2) {
        for (int dim = 0; dim < 3; ++dim) {
            tgt::vec3 direction = tgt::vec3(0.f);
            direction.elem[dim] = static_cast<float>(sign) * 1.f;
            tgt::vec3 color = direction * tgt::vec3(0.5f);
            color += tgt::vec3(0.5f);

            TriangleMeshGeometryUInt16IndexedColorNormal* axis = createArrowGeometry(tgt::vec4(color,1.f));
            tgt::vec3 rotationAxis = tgt::cross(tgt::vec3(0.f,0.f,1.f), direction);
            if (tgt::length(rotationAxis) < 0.001f)
                rotationAxis = tgt::vec3(1.f, 0.f, 0.f);
            float rotationAngle = acos(tgt::dot(direction, tgt::vec3(0.f,0.f,1.f)));
            if (rotationAngle > 0.001f)
                axis->setTransformationMatrix(tgt::Matrix::createRotation(rotationAngle, rotationAxis) * tgt::Matrix::createScale(tgt::vec3(0.85f))
                                          * tgt::Matrix::createTranslation(tgt::vec3(0.f, 0.f, overlayBaseSize_ * overlaySizeProp_.get() * 0.4f)));
            else
                axis->setTransformationMatrix(tgt::Matrix::createScale(tgt::vec3(0.85f))
                                          * tgt::Matrix::createTranslation(tgt::vec3(0.f, 0.f, overlayBaseSize_ * overlaySizeProp_.get() * 0.4f)));

            mesh->addMesh(axis);
            delete axis;
        }
    }

    for (int sign = -1; sign <= 1; sign += 2) {
        for (int dim = -1; dim < 3; ++dim) {
            tgt::vec3 direction = tgt::normalize(static_cast<float>(sign) * tgt::vec3(1.f));
            tgt::vec3 color = static_cast<float>(sign) * tgt::vec3(0.5f);
            if (dim >= 0) {
                size_t index = static_cast<size_t>(dim);
                direction.elem[index] = -direction.elem[index];
                color.elem[index] = -color.elem[index];
            }
            color += tgt::vec3(0.5f);

            TriangleMeshGeometryUInt16IndexedColorNormal* axis = createArrowGeometry(tgt::vec4(color,1.f));
            tgt::vec3 rotationAxis = tgt::cross(tgt::vec3(0.f,0.f,1.f), direction);
            float rotationAngle = acos(tgt::dot(direction, tgt::vec3(0.f,0.f,1.f)));
            axis->setTransformationMatrix(tgt::Matrix::createRotation(rotationAngle, rotationAxis) * tgt::Matrix::createScale(tgt::vec3(0.85f))
                                          * tgt::Matrix::createTranslation(tgt::vec3(0.f, 0.f, overlayBaseSize_ * overlaySizeProp_.get() * 0.4f)));
            mesh->addMesh(axis);
            delete axis;
        }
    }

    // transform by the additional rotation matrix
    mesh->transform(colorRotationMatrix_.get());

    currentGeometry_ =  mesh;

}

void FlowDirectionOverlay::createCubeGeometry() {
    //clear old geom
    delete currentGeometry_;

    // create new geometry
    TriangleMeshGeometry<VertexColorNormal>* geom = new TriangleMeshGeometryColorNormal();

    for (int sign = -1; sign <= 1; sign += 2) {
        for (size_t dim = 0; dim < 3; ++dim) {
            tgt::vec3 normal = tgt::vec3(0.f);
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
                vertices[i].color_ = tgt::vec4(vertices[i].pos_ + tgt::vec3(0.5f), 1.f);
                vertices[i].normal_ = normal;
            }

            if (sign > 0)
                geom->addQuad(vertices[0], vertices[1], vertices[2], vertices[3]);
            else
                geom->addQuad(vertices[0], vertices[3], vertices[2], vertices[1]);
        }
    }

    geom->setTransformationMatrix(tgt::Matrix::createScale(tgt::vec3(overlayBaseSize_ * overlaySizeProp_.get())));

    // transform by the additional rotation matrix
    geom->transform(colorRotationMatrix_.get());

    currentGeometry_ = geom;
}

void FlowDirectionOverlay::createSphereGeometry() {
    //clear old geom
    delete sphereGeometry_;

    GlMeshGeometryUInt16ColorNormal* mesh = new GlMeshGeometryUInt16ColorNormal();

    // set sphere geometry manually because we also have to set the colors
    float radius = 0.5f;
    size_t numSlicesAndTiles = 200;

    // create the sphere vertices
    const float dfi = tgt::PIf / static_cast<float>(numSlicesAndTiles);
    const float dth = (2.0f * tgt::PIf) / static_cast<float>(numSlicesAndTiles);
    const float du  = 1.0f / static_cast<float>(numSlicesAndTiles);
    const float dv  = 1.0f / static_cast<float>(numSlicesAndTiles);
    std::vector<float>  cs_fi, cs_th, sn_fi, sn_th;

    float fi = 0.f;// PI / 2.0;
    float th = 0.f;
    float u  = 0.f;
    float v  = 0.f;

    for (size_t i = 0; i <= numSlicesAndTiles; ++i) {
        cs_fi.push_back(std::cos(fi));
        cs_th.push_back(std::cos(th));

        sn_fi.push_back(std::sin(fi));
        sn_th.push_back(std::sin(th));

        fi += dfi;
        th += dth;
    }

    for (size_t i = 0; i <= numSlicesAndTiles; ++i) {
        for (size_t j = 0; j <= numSlicesAndTiles; ++j) {

            size_t k = j % numSlicesAndTiles;
            tgt::vec3 normal = tgt::vec3(sn_fi[i] * cs_th[k], sn_fi[i] * sn_th[k], cs_fi[i]);
            tgt::vec3 position = radius * normal;
            tgt::vec2 uv = tgt::vec2(u, v);

            // compute the vertex color
            tgt::vec4 color = tgt::vec4(position + tgt::vec3(0.5f), 1.f);

            // transform the position
            position = (colorRotationMatrix_.get() * tgt::vec4(position, 1.f)).xyz();

            mesh->addVertex(position, normal, color, uv);

            u += du;
        }

        u = 0;
        v += dv;
    }

    // create the indices for the triangles
    uint16_t offset  = 0;
    uint16_t offsetDelta = static_cast<uint16_t>(numSlicesAndTiles) + 1;
    for (size_t j = 0; j < numSlicesAndTiles; ++j) {
        for (size_t i = 1; i <= numSlicesAndTiles; ++i) {
            uint16_t idx = offset + static_cast<uint16_t>(i);
            mesh->addIndex(idx - 1);
            mesh->addIndex(idx + numSlicesAndTiles);
            mesh->addIndex(idx);

            mesh->addIndex(idx);
            mesh->addIndex(idx + numSlicesAndTiles);
            mesh->addIndex(idx + numSlicesAndTiles + 1);
        }
        offset += offsetDelta;
    }

    mesh->setPrimitiveType(GL_TRIANGLE_STRIP);
    mesh->disablePrimitiveRestart();

    mesh->setTransformationMatrix(tgt::Matrix::createScale(tgt::vec3(overlayBaseSize_ * overlaySizeProp_.get())));

    sphereGeometry_ = mesh;
}

//--------------------------------------------------------------------------------------------
//          Callbacks
//--------------------------------------------------------------------------------------------
void FlowDirectionOverlay::invalidateGeometry() {
    geometryMustBeRecreated_ = true;
    invalidate();
}


} // namespace
