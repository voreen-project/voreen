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

#include "geometryrenderer.h"
#include "voreen/core/datastructures/geometry/geometry.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"
#include "voreen/core/datastructures/meta/serializablevectormetadata.h"
#include "tgt/texturemanager.h"

namespace voreen {

GeometryRenderer::GeometryRenderer()
    : GeometryRendererBase()
    , inport_(Port::INPORT, "inport.geometry", "Geometry Input", false)
    , texPort_(Port::INPORT, "inport.texture", "Texture Input")
    , polygonMode_("polygonMode", "Polygon Mode")
    , lineWidth_("lineWidth", "Line Width", 1.f, 1.f, 20.f)
    , pointSize_("pointSize", "Point Size", 1.f, 1.f, 20.f)
    , solidColor_("solidColor", "Color", tgt::Color(1.f, 1.f, 1.f, 1.f))
    , mapTexture_("mapTexture", "Map Texture", false)
    , textureMode_("textureMode", "Texture Mode")
    , enableLighting_("enableLighting", "Enable Lighting", false)
    , lightPosition_("lightPosition", "Light Source Position", tgt::vec4(2.3f, 1.5f, 1.5f, 1.f), tgt::vec4(-10000), tgt::vec4(10000))
    , lightAmbient_("lightAmbient", "Ambient Light", tgt::Color(0.4f, 0.4f, 0.4f, 1.f))
    , lightDiffuse_("lightDiffuse", "Diffuse Light", tgt::Color(0.8f, 0.8f, 0.8f, 1.f))
    , lightSpecular_("lightSpecular", "Specular Light", tgt::Color(0.6f, 0.6f, 0.6f, 1.f))
    , materialShininess_("materialShininess", "Shininess", 60.f, 0.1f, 128.f)
    , shaderOpaqueProp_("geometry_opaque.prg", "Shader opaque", "trianglemesh.frag", "trianglemesh.vert", "trianglemesh.geom",
                    Processor::INVALID_PROGRAM,Property::LOD_DEBUG)
    , shaderTransparentProp_("geometry_transparent.prg", "Shader transparent", "trianglemesh.frag", "trianglemesh.vert", "trianglemesh.geom",
                    Processor::INVALID_PROGRAM,Property::LOD_DEBUG)
    , enableClipping_("enable.clipping", "Enable on-the-fly clipping", false)
    , planeNormal_("plane.normal", "Clipping plane normal", tgt::vec3(0.f, 0.f, 1.f), tgt::vec3(-1.f), tgt::vec3(1.f))
    , planeDistance_("plane.distance", "Clipping plane distance", 0.f, -10000.f, 10000.f)
    , lastShaderConfiguration_(-1)
    , shaderNeedsRebuild_(true)
    , fbo_(0)
{
    addPort(inport_);
    addPort(texPort_);

    polygonMode_.addOption("point", "Point", GL_POINT);
    polygonMode_.addOption("line",  "Line",  GL_LINE);
    polygonMode_.addOption("fill",  "Fill",  GL_FILL);
    polygonMode_.select("fill");
    polygonMode_.onChange(MemberFunctionCallback<GeometryRenderer>(this, &GeometryRenderer::updatePropertyVisibilities));
    addProperty(polygonMode_);

    addProperty(pointSize_);
    addProperty(lineWidth_);
    addProperty(solidColor_);

    mapTexture_.onChange(MemberFunctionCallback<GeometryRenderer>(this, &GeometryRenderer::invalidateShader));
    mapTexture_.onChange(MemberFunctionCallback<GeometryRenderer>(this, &GeometryRenderer::updatePropertyVisibilities));
    addProperty(mapTexture_);

    textureMode_.addOption("REPLACE", "replace");
    textureMode_.addOption("MODULATE", "modulate");
    textureMode_.addOption("DECAL", "decal");
    textureMode_.addOption("BLEND", "blend");
    textureMode_.selectByKey("REPLACE");
    textureMode_.onChange(MemberFunctionCallback<GeometryRenderer>(this, &GeometryRenderer::invalidateShader));
    addProperty(textureMode_);

    enableLighting_.onChange(MemberFunctionCallback<GeometryRenderer>(this, &GeometryRenderer::updatePropertyVisibilities));

    addProperty(enableLighting_);
    addProperty(lightPosition_);
    addProperty(lightAmbient_);
    addProperty(lightDiffuse_);
    addProperty(lightSpecular_);
    addProperty(materialShininess_);

    addProperty(shaderOpaqueProp_);
    addProperty(shaderTransparentProp_);
    addProperty(enableClipping_);
    addProperty(planeNormal_);
    addProperty(planeDistance_);

    enableClipping_.onChange(MemberFunctionCallback<GeometryRenderer>(this, &GeometryRenderer::updatePropertyVisibilities));

    // assign lighting properties to property group
    lightPosition_.setGroupID("lighting");
    lightAmbient_.setGroupID("lighting");
    lightDiffuse_.setGroupID("lighting");
    lightSpecular_.setGroupID("lighting");
    materialShininess_.setGroupID("lighting");
    setPropertyGroupGuiName("lighting", "Lighting Parameters");

    updatePropertyVisibilities();
}

Processor* GeometryRenderer::create() const {
    return new GeometryRenderer();
}

GeometryRenderer::~GeometryRenderer() {
}

void GeometryRenderer::setupShaders() {
    const std::string defines = generateShaderDefines();

    shaderOpaqueProp_.setHeader(defines);
    shaderOpaqueProp_.rebuild();

    shaderTransparentProp_.setHeader(defines + "#define USE_TRANSPARENCY\n");
    shaderTransparentProp_.rebuild();

    shaderNeedsRebuild_ = false;
    if(const TriangleMeshGeometryBase* tmgb = dynamic_cast<const TriangleMeshGeometryBase*>(inport_.getData())) {
        lastShaderConfiguration_ = static_cast<int>(tmgb->getVertexLayout());
    }
    else if (const GlMeshGeometryBase* glmgb = dynamic_cast<const GlMeshGeometryBase*>(inport_.getData())) {
        lastShaderConfiguration_ = static_cast<int>(glmgb->getVertexLayout());
    }
}

bool GeometryRenderer::isReady() const {
    return inport_.isReady();
}

void GeometryRenderer::invalidateShader() {
    shaderNeedsRebuild_ = true;
}

bool GeometryRenderer::usesTransparency() const {
    const Geometry* currentGeometry = inport_.getData();
    return currentGeometry && (currentGeometry->isTransparent() || (isUsingSolidColor() && solidColor_.get().a < 1));
}
void GeometryRenderer::render() {
    render(shaderOpaqueProp_);
}
void GeometryRenderer::renderTransparent() {
    render(shaderTransparentProp_);
}
tgt::Bounds GeometryRenderer::getBoundingBox() const {
    if (inport_.hasData())
        return inport_.getData()->getBoundingBox();

    return GeometryRendererBase::getBoundingBox();
}
void GeometryRenderer::render(ShaderProperty& shaderProp) {
    tgtAssert(inport_.hasData(), "No geometry");

    // check, if a shader is needed
    bool isTriangleMeshShaderNeeded = false;
    bool solidColorNeeded = isUsingSolidColor();
    const TriangleMeshGeometryBase* tmgb = 0;
    const GlMeshGeometryBase* glmgb = 0;
    if((tmgb = dynamic_cast<const TriangleMeshGeometryBase*>(inport_.getData()))) {
        isTriangleMeshShaderNeeded = true;
    }
    else if (glmgb = dynamic_cast<const GlMeshGeometryBase*>(inport_.getData())) {
        isTriangleMeshShaderNeeded = true;
    }
    else {
        lastShaderConfiguration_ = -1;
    }

    bool shaderUpdateNeeded = false;

    static bool texPortConnected = texPort_.hasData();
    if(isTriangleMeshShaderNeeded && texPortConnected != texPort_.hasData()) {
        texPortConnected = texPort_.hasData();
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
        setupShaders();

    if(isTriangleMeshShaderNeeded && !shaderProp.hasValidShader()) {
        LERROR("Shader for geometry failed to compile");
        return;
    }

    glPolygonMode(GL_FRONT_AND_BACK, polygonMode_.getValue());
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

        // texture data
        prog->setIgnoreUniformLocationError(true);
        tgt::TextureUnit texUnit;
        if(mapTexture_.get()) {
            if(texPort_.isReady()) {
                texUnit.activate();
                texPort_.bindColorTexture(texUnit.getEnum());
                prog->setUniform("texture_", texUnit.getUnitNumber());
            }
            else if(tmgb && tmgb->getTextureData()) {
                LGL_ERROR;
                texUnit.activate();
                tmgb->getTextureData()->bind();
                prog->setUniform("textures_", texUnit.getUnitNumber());
                LGL_ERROR;
            }
            else if(glmgb && glmgb->getTextureData()) {
                LGL_ERROR;
                texUnit.activate();
                glmgb->getTextureData()->bind();
                prog->setUniform("textures_", texUnit.getUnitNumber());
                LGL_ERROR;
            }
        }

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
    tgt::TextureUnit::setZeroUnit();
}

void GeometryRenderer::process() {
}

void GeometryRenderer::adjustPropertiesToInput() {
    const TriangleMeshGeometryBase* tmgb = dynamic_cast<const TriangleMeshGeometryBase*>(inport_.getData());
    const GlMeshGeometryBase* glmgb = dynamic_cast<const GlMeshGeometryBase*>(inport_.getData());
    if(tmgb || glmgb) {
        shaderOpaqueProp_.setVisibleFlag(true);
        shaderTransparentProp_.setVisibleFlag(true);
        enableClipping_.setVisibleFlag(true);
        planeNormal_.setVisibleFlag(enableClipping_.get());
        planeDistance_.setVisibleFlag(enableClipping_.get());
    }
    else {
        shaderOpaqueProp_.setVisibleFlag(false);
        shaderTransparentProp_.setVisibleFlag(false);
        enableClipping_.setVisibleFlag(false);
        planeNormal_.setVisibleFlag(false);
        planeDistance_.setVisibleFlag(false);
    }
    solidColor_.setVisibleFlag(isUsingSolidColor());
}

void GeometryRenderer::updatePropertyVisibilities() {
    pointSize_.setVisibleFlag(polygonMode_.isSelected("point"));
    lineWidth_.setVisibleFlag(polygonMode_.isSelected("line"));

    bool lighting = enableLighting_.get();
    setPropertyGroupVisible("lighting", lighting);

    textureMode_.setVisibleFlag(mapTexture_.get());
}

std::string GeometryRenderer::generateShaderDefines() {
    std::string header;
    if(const TriangleMeshGeometryBase* tmgb = dynamic_cast<const TriangleMeshGeometryBase*>(inport_.getData())) {
        header += tmgb->getShaderDefines();
        if(mapTexture_.get()) {
            if(texPort_.isReady())
                header += "#define USE_SINGLE_TEXTURE\n";
            else if(mapTexture_.get() && tmgb->supportsTextureData() && tmgb->hasTextureData())
                header += "#define USE_TEXTURE_ARRAY\n";

            header += "#define TEXTURE_MODE_" + textureMode_.getValue() + "\n";
        }
    }
    else if (const GlMeshGeometryBase* glmgb = dynamic_cast<const GlMeshGeometryBase*>(inport_.getData())) {
        header += glmgb->getShaderDefines();
        if(mapTexture_.get()) {
            if(texPort_.isReady())
                header += "#define USE_SINGLE_TEXTURE\n";
            else if(mapTexture_.get() && glmgb->supportsTextureData() && glmgb->hasTextureData())
                header += "#define USE_TEXTURE_ARRAY\n";

            header += "#define TEXTURE_MODE_" + textureMode_.getValue() + "\n";
        }
    }
    return header;
}
bool GeometryRenderer::isUsingSolidColor() const {
    bool solidColorNeeded = false;
    const TriangleMeshGeometryBase* tmgb = 0;
    const GlMeshGeometryBase* glmgb = 0;
    if((tmgb = dynamic_cast<const TriangleMeshGeometryBase*>(inport_.getData())))
        return !tmgb->supportsColors();
    else if (glmgb = dynamic_cast<const GlMeshGeometryBase*>(inport_.getData()))
        return !glmgb->supportsColors();
    else
        return false;
}

}

