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

#include "quadricrenderer.h"

namespace voreen {

using tgt::vec4;
using tgt::vec3;

QuadricRenderer::QuadricRenderer()
    : GeometryRendererBase()
    , enabled_("enabled", "Enabled", true)
    , quadricType_("quadricType", "Quadric Type")
    , position_("position", "Position", tgt::vec3(0.f), tgt::vec3(-99.0f), tgt::vec3(99.0f))
    , closeCylinder_("closecylinder", "Close Cylinder", false)
    , start_("start", "Start Position", tgt::vec3(0.f), tgt::vec3(-99.0f), tgt::vec3(99.0f))
    , end_("end", "End Position", tgt::vec3(1.f), tgt::vec3(-99.0f), tgt::vec3(99.0f))
    , radius_("radius", "Radius", 0.2f, 0.001f, 5.0f)
    , color_("color", "Color", tgt::Color(0.75f, 0.25f, 0.f, 1.f))
    , applyLighting_("applyLighting", "Apply Lighting", true)
    , lightPosition_("lightPosition", "Light Source Position", tgt::vec3(-2.f, 2.f, 2.f), tgt::vec3(-10), tgt::vec3(10))
    , lightAmbient_("lightAmbient", "Ambient Light", tgt::Color(0.4f, 0.4f, 0.4f, 1.f))
    , lightDiffuse_("lightDiffuse", "Diffuse Light", tgt::Color(0.6f, 0.6f, 0.6f, 1.f))
    , lightSpecular_("lightSpecular", "Specular Light", tgt::Color(0.4f, 0.4f, 0.4f, 1.f))
    , materialShininess_("materialShininess", "Shininess", 60.f, 0.1f, 128.f)
    , shader_(0)
{
    quadricType_.addOption("cylinder", "Cylinder");
    quadricType_.addOption("sphere",   "Sphere");

    addProperty(enabled_);
    addProperty(quadricType_);
    addProperty(closeCylinder_);
    addProperty(start_);
    addProperty(end_);
    addProperty(position_);
    addProperty(radius_);
    addProperty(color_);

    addProperty(applyLighting_);
    addProperty(lightPosition_);
    addProperty(lightAmbient_);
    addProperty(lightDiffuse_);
    addProperty(lightSpecular_);
    addProperty(materialShininess_);

    // assign lighting properties to property group
    lightPosition_.setGroupID("lighting");
    lightAmbient_.setGroupID("lighting");
    lightDiffuse_.setGroupID("lighting");
    lightSpecular_.setGroupID("lighting");
    materialShininess_.setGroupID("lighting");
    setPropertyGroupGuiName("lighting", "Lighting Parameters");

    quadricType_.onChange(MemberFunctionCallback<QuadricRenderer>(this, &QuadricRenderer::adjustPropertyVisibilities));

    quadricType_.onChange(MemberFunctionCallback<QuadricRenderer>(this, &QuadricRenderer::invalidateGeometry));
    closeCylinder_.onChange(MemberFunctionCallback<QuadricRenderer>(this, &QuadricRenderer::invalidateGeometry));
    start_.onChange(MemberFunctionCallback<QuadricRenderer>(this, &QuadricRenderer::invalidateGeometry));
    end_.onChange(MemberFunctionCallback<QuadricRenderer>(this, &QuadricRenderer::invalidateGeometry));
    radius_.onChange(MemberFunctionCallback<QuadricRenderer>(this, &QuadricRenderer::invalidateGeometry));
    position_.onChange(MemberFunctionCallback<QuadricRenderer>(this, &QuadricRenderer::invalidateGeometry));
    color_.onChange(MemberFunctionCallback<QuadricRenderer>(this, &QuadricRenderer::invalidateGeometry));

    applyLighting_.onChange(MemberFunctionCallback<QuadricRenderer>(this, &QuadricRenderer::adjustPropertyVisibilities));
}

Processor* QuadricRenderer::create() const {
    return new QuadricRenderer();
}

void QuadricRenderer::initialize() {
    shader_ = ShdrMgr.load("quadricrenderer", "", false);

    GeometryRendererBase::initialize();
    adjustPropertyVisibilities();
}

void QuadricRenderer::deinitialize() {
    ShdrMgr.dispose(shader_);
    GeometryRendererBase::deinitialize();
}

tgt::Bounds QuadricRenderer::getBoundingBox() const {
    if (!mesh_.isEmpty())
        return mesh_.getBoundingBox();

    return GeometryRendererBase::getBoundingBox();
}

void QuadricRenderer::render() {

    if (!enabled_.get())
        return;

    // first time: create initial mesh
    if (mesh_.isEmpty())
        invalidateGeometry();

    // activate shader
    shader_->activate();

    // set matrix stack uniforms
    shader_->setUniform("modelViewMatrixStack_", MatStack.getModelViewMatrix() * mesh_.getTransformationMatrix());
    shader_->setUniform("modelViewProjectionMatrixStack_", MatStack.getProjectionMatrix() * MatStack.getModelViewMatrix() * mesh_.getTransformationMatrix());
    tgt::mat4 viewInverse;
    tgt::mat4 mv = MatStack.getModelViewMatrix() * mesh_.getTransformationMatrix();
    mv.invert(viewInverse);
    shader_->setUniform("modelViewMatrixInverseStack_", viewInverse);

    // set lighting parameters
    shader_->setUniform("lightingEnabled_", applyLighting_.get());
    // set light source
    std::string prefix = "lightSource_.";
    shader_->setUniform(prefix+"position_", tgt::vec4(lightPosition_.get(), 1.f));
    shader_->setUniform(prefix+"ambientColor_", lightAmbient_.get().xyz());
    shader_->setUniform(prefix+"diffuseColor_", lightDiffuse_.get().xyz());
    shader_->setUniform(prefix+"specularColor_", lightSpecular_.get().xyz());

    // material
    prefix = "material_.";
    shader_->setUniform(prefix+"ambientColor_", color_.get().xyz());
    shader_->setUniform(prefix+"diffuseColor_", color_.get().xyz());
    shader_->setUniform(prefix+"specularColor_", color_.get().xyz());
    shader_->setUniform(prefix+"shininess_", materialShininess_.get());

    shader_->setUniform("color_", color_.get());

    mesh_.render();

    shader_->deactivate();

    LGL_ERROR;
}

void QuadricRenderer::adjustPropertyVisibilities() {
    bool cylinder = quadricType_.isSelected("cylinder");
    bool sphere = quadricType_.isSelected("sphere");
    bool lighting = applyLighting_.get();

    closeCylinder_.setVisibleFlag(cylinder);
    start_.setVisibleFlag(cylinder);
    end_.setVisibleFlag(cylinder);
    position_.setVisibleFlag(sphere);
    setPropertyGroupVisible("lighting", lighting);
}

void QuadricRenderer::invalidateGeometry() {

    if (!isInitialized())
        return;

    mesh_.clear();

    if (quadricType_.isSelected("cylinder")) {
        //calculate correct rotation matrix:
        vec3 rotz = normalize(end_.get() - start_.get());
        vec3 roty = normalize(vec3(rotz.y, -rotz.z, 0.0f));
        vec3 rotx = cross(roty, rotz);

        float m[16];

        m[0] = rotx.x;
        m[1] = rotx.y;
        m[2] = rotx.z;
        m[3] = 0.0f;

        m[4] = roty.x;
        m[5] = roty.y;
        m[6] = roty.z;
        m[7] = 0.0f;

        m[8] = rotz.x;
        m[9] = rotz.y;
        m[10] = rotz.z;
        m[11] = 0.0f;

        m[12] = 0.0f;
        m[13] = 0.0f;
        m[14] = 0.0f;
        m[15] = 1.0f;

        tgt::mat4 matrix(m);
        matrix = tgt::mat4::createTranslation(start_.get()) * tgt::transpose(matrix);

        float l = length(start_.get() - end_.get());

        mesh_.setCylinderGeometry(color_.get(), radius_.get(), radius_.get(), l, 200, 200, closeCylinder_.get(), closeCylinder_.get());
        mesh_.setTransformationMatrix(matrix);
    }
    else if (quadricType_.isSelected("sphere")) {
        mesh_.setSphereGeometry(radius_.get(), tgt::vec3::zero, color_.get(), 20);
        mesh_.setTransformationMatrix(tgt::mat4::createTranslation(position_.get()));
    }
    else {
        LERROR("Unknown quadric type: " << quadricType_.get());
    }
}

}
