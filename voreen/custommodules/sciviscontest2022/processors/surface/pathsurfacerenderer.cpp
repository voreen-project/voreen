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


#include "pathsurfacerenderer.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

namespace voreen {

PathSurfaceRenderer::PathSurfaceRenderer():
    GeometryRenderer()
    ,pathSurfaceShaderProperty_("geometry_pathsurface.prg", "Pathsurface shader", "pathsurface.frag", "trianglemesh.vert", "trianglemesh.geom",
                Processor::INVALID_PROGRAM,Property::LOD_ADVANCED)
    ,enableTimelines_("enableTimelines","Enable Timelines",false)
    ,enableTransparency_("enableTransparency","Enable Transparency",false)
    ,texCoordYScale_("texCoordYScale","Scaling of Y Coordinate",1, 0, 5)
    ,timelineWidth_("timelineWidth","Width of Timeline",0.05f, 0.01f, 0.5f)
    ,timeInterval_("timeInterval", "Time Interval",{0,1}, 0, 100)
    ,gradient_("gradient", "Timeline gradient")
    ,useGradient_("useGradient", "Use Gradient:", false) {
    
    gradient_.setDomainFittingStrategy(TransFuncPropertyBase::DomainAutoFittingStrategy::FIT_DOMAIN_ALWAYS);

    addProperty(enableTimelines_);
    addProperty(enableTransparency_);
    addProperty(texCoordYScale_);
    addProperty(timelineWidth_);
    addProperty(timeInterval_);
    addProperty(useGradient_);
    addProperty(gradient_);

    addProperty(pathSurfaceShaderProperty_);

    ON_CHANGE_LAMBDA(enableTimelines_, [&](){
        shaderNeedsRebuild_ = true;
        renderTransparent();
    });
    ON_CHANGE_LAMBDA(enableTransparency_, [&](){
        shaderNeedsRebuild_ = true;
        renderTransparent();
    });
    ON_CHANGE_LAMBDA(texCoordYScale_, [&](){
        renderTransparent();
    });
    ON_CHANGE_LAMBDA(timelineWidth_, [&](){
        renderTransparent();
    });
    ON_CHANGE_LAMBDA(timeInterval_, [&](){
        // Update transfer function.
        float* data = new float[2];
        data[0] = timeInterval_.get().x;
        data[1] = timeInterval_.get().y;
        VolumeRAM_Float* representation = new VolumeRAM_Float(data, tgt::svec3(2,1,1));
        tfVolume_.reset(new Volume(representation, tgt::vec3::one, tgt::vec3::zero));
        tfVolume_->addDerivedData(new VolumeMinMax(data[0], data[1], data[0], data[1])); // To save time by not triggering a background thread.
        gradient_.setVolume(tfVolume_.get());
        gradient_.applyDomainFromData();
        renderTransparent();
    });
    ON_CHANGE_LAMBDA(gradient_, [&]() {
        renderTransparent();
    });
    ON_CHANGE_LAMBDA(useGradient_, [&]() {
        shaderNeedsRebuild_ = true;
        renderTransparent();
    });
}

void PathSurfaceRenderer::render() {
    renderTransparent();
}

void PathSurfaceRenderer::render(ShaderProperty& shaderProp) {
    renderTransparent();
}

void PathSurfaceRenderer::setupPathsurfaceShader() {
    std::string defines = generateShaderDefines();
    defines += "#define USE_TRANSPARENCY\n";
    if (useGradient_.get()) {
        defines += "#define USE_GRADIENT\n";
        defines += gradient_.get()->getShaderDefines();
    }

    pathSurfaceShaderProperty_.setHeader(defines);
    pathSurfaceShaderProperty_.rebuild();

    shaderNeedsRebuild_ = false;
    if(const TriangleMeshGeometryBase* tmgb = dynamic_cast<const TriangleMeshGeometryBase*>(inport_.getData())) {
        lastShaderConfiguration_ = static_cast<int>(tmgb->getVertexLayout());
    }
    else if (const GlMeshGeometryBase* glmgb = dynamic_cast<const GlMeshGeometryBase*>(inport_.getData())) {
        lastShaderConfiguration_ = static_cast<int>(glmgb->getVertexLayout());
    }
}

bool PathSurfaceRenderer::isReady() const {
    auto ready = inport_.isReady();
    ready &= outPort_.isConnected();
    return ready;
}

void PathSurfaceRenderer::renderTransparent() {
    auto& shaderProp = pathSurfaceShaderProperty_;

    if (!inport_.hasData())
        return;

    tgtAssert(inport_.hasData(), "No geometry");
    // check, if a shader is needed
    bool isTriangleMeshShaderNeeded = false;
    bool solidColorNeeded = isUsingSolidColor();
    const TriangleMeshGeometryBase* tmgb = nullptr;
    const GlMeshGeometryBase* glmgb = nullptr;
    if((tmgb = dynamic_cast<const TriangleMeshGeometryBase*>(inport_.getData()))) {
        isTriangleMeshShaderNeeded = true;
    }
    else if (glmgb = dynamic_cast<const GlMeshGeometryBase*>(inport_.getData())) {
        isTriangleMeshShaderNeeded =
                (glmgb->getPrimitiveType() == GL_TRIANGLE_STRIP) ||
                (glmgb->getPrimitiveType() == GL_TRIANGLE_FAN) ||
                (glmgb->getPrimitiveType() == GL_TRIANGLES) ||
                (glmgb->getPrimitiveType() == GL_TRIANGLE_STRIP_ADJACENCY) ||
                (glmgb->getPrimitiveType() == GL_TRIANGLES_ADJACENCY );
    }
    else {
        lastShaderConfiguration_ = -1;
    }

    bool shaderUpdateNeeded = shaderProp.requiresRebuild();

    if(isTriangleMeshShaderNeeded && texPort_.hasChanged()) {
        shaderUpdateNeeded = true;
    }

    if(shaderNeedsRebuild_ ||
                    (tmgb && (isTriangleMeshShaderNeeded &&
                    ((tmgb->getVertexLayout() != lastShaderConfiguration_)   ||
                    !shaderProp.hasValidShader()                          ||
                    getInvalidationLevel() >= Processor::INVALID_PROGRAM)))
                    || (glmgb && (isTriangleMeshShaderNeeded &&
                    ((glmgb->getVertexLayout() != lastShaderConfiguration_)   ||
                    !shaderProp.hasValidShader()                          ||
                    getInvalidationLevel() >= Processor::INVALID_PROGRAM)))
    )
    {
        shaderUpdateNeeded = true;
    }

    if(shaderUpdateNeeded)
        setupPathsurfaceShader();

    if(isTriangleMeshShaderNeeded && !shaderProp.hasValidShader()) {
        LERROR("Shader for geometry failed to compile");
        return;
    }

    glPolygonMode(GL_FRONT_AND_BACK, polygonMode_.getValue());
    GLenum cullFace = cullFace_.getValue();
    if(cullFace != GL_NONE) {
        glEnable(GL_CULL_FACE);
        glCullFace(cullFace);
    }
#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
    glShadeModel(GL_SMOOTH);
#endif

    if (polygonMode_.isSelected("point"))
        glPointSize(pointSize_.get());
    else if (polygonMode_.isSelected("line"))
        glLineWidth(lineWidth_.get());

#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
    if (!isTriangleMeshShaderNeeded && mapTexture_.get() && texPort_.isReady()) {
        texPort_.bindColorTexture(GL_TEXTURE0);
        glEnable(GL_TEXTURE_2D);
    }

    if (enableLighting_.get()) {
        glEnable(GL_LIGHTING);
        glEnable(GL_COLOR_MATERIAL);
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, tgt::vec4(0.f).elem);

        glEnable(GL_LIGHT0);
        glLightfv(GL_LIGHT0, GL_POSITION, lightPosition_.get().elem);
        glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient_.get().elem);
        glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse_.get().elem);
        glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular_.get().elem);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, materialShininess_.get());
        LGL_ERROR;
    }
#endif

    // for triangle meshes, we need to use a correctly layouted shader
    if(isTriangleMeshShaderNeeded) {
        tgt::Shader* prog = shaderProp.getShader();
        prog->activate();
        setGlobalShaderParameters(prog, &camera_);

        LGL_ERROR;
        prog->setUniform("timeInterval_", timeInterval_.get());
        prog->setUniform("timelineWidth_", timelineWidth_.get());
        prog->setUniform("texCoordYScale_", texCoordYScale_.get());
        prog->setUniform("enableTransparency_", enableTransparency_.get());
        prog->setUniform("enableTimelines_", enableTimelines_.get());
        LGL_ERROR;

        // texture data
        prog->setIgnoreUniformLocationError(true);
        tgt::TextureUnit texUnit;
        LGL_ERROR;
        if (useGradient_.get()) {
            texUnit.activate();
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_REPEAT);

            gradient_.get()->getTexture()->bind();
            gradient_.get()->setUniform(prog, "transFuncParam_", "transFuncTex_", texUnit.getUnitNumber());
        }
        else if (mapTexture_.get() && texPort_.isReady()) {
            texUnit.activate();
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_REPEAT);

            texPort_.bindColorTexture(texUnit.getEnum());
            prog->setUniform("texture_", texUnit.getUnitNumber());
        }
        LGL_ERROR;

        tgt::mat3 modelView = (camera_.getViewMatrix() * inport_.getData()->getTransformationMatrix()).getRotationalPartMat3();
        tgt::mat3 normalMatrix;
        if(!modelView.invert(normalMatrix)) {
            LWARNING("Could not generate normal matrix out of current view / model matrix, using identity.");
            normalMatrix = tgt::mat3::identity;
        }
        normalMatrix = transpose(normalMatrix);

        prog->setUniform("normalMatrix_", normalMatrix);
        prog->setUniform("modelMatrix_", inport_.getData()->getTransformationMatrix());
        prog->setUniform("viewMatrix_", camera_.getViewMatrix());
        prog->setUniform("projectionMatrix_", camera_.getProjectionMatrix(viewport_));
        prog->setUniform("enableLighting_", enableLighting_.get());
        LGL_ERROR;

        if (solidColorNeeded)
            prog->setUniform("solidColor_", solidColor_.get());

        // lighting
        if(enableLighting_.get()) {
            prog->setUniform("lightPositionEye_", (camera_.getViewMatrix() * tgt::vec4(lightPosition_.get().xyz(), 1.f)).xyz());
            prog->setUniform("lightSource_.ambientColor_", lightAmbient_.get().xyz());
            prog->setUniform("lightSource_.diffuseColor_", lightDiffuse_.get().xyz());
            prog->setUniform("lightSource_.specularColor_", lightSpecular_.get().xyz());
            //prog->setUniform("lightSource_.attenuation_", tgt::vec3(1.f, 0.f, 0.f));
            prog->setUniform("shininess_", materialShininess_.get());
        }
        prog->setIgnoreUnsetUniform("headPointerImage_"); // Disable warning for uniform that is already set in GeometryProcessor.
        prog->setIgnoreUniformLocationError(false);

        prog->setUniform("enableClipping_", enableClipping_.get());
        if (enableClipping_.get())
            prog->setUniform("plane_", tgt::vec4(normalize(planeNormal_.get()), planeDistance_.get()));
        else // set to zero to prevent undefined behavior
            prog->setUniform("plane_", tgt::vec4(0.0f, 0.0f, 0.0f, 0.0f));

    }

    LGL_ERROR;
    inport_.getData()->render();
    LGL_ERROR;

    if(isTriangleMeshShaderNeeded)
        shaderProp.getShader()->deactivate();

#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHT0);
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);
#endif
    glPointSize(1.f);
    glLineWidth(1.f);
    glColor4f(1.f, 1.f, 1.f, 1.f);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDisable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    tgt::TextureUnit::setZeroUnit();
}

void PathSurfaceRenderer::adjustPropertiesToInput() {
    GeometryRenderer::adjustPropertiesToInput();
    shaderNeedsRebuild_ = true;
}

}