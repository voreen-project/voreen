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

#include "pointlistrenderer.h"

#include "voreen/core/datastructures/geometry/pointlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

#include "tgt/immediatemode/immediatemode.h"

namespace voreen {

const std::string PointListRenderer::loggerCat_("voreen.PointListRenderer");

PointListRenderer::PointListRenderer()
    : GeometryRendererBase()
    , coordinateSystem_("coordinateSystem", "Coordinate System")
    , renderingPrimitiveProp_("renderingPrimitive", "Rendering Primitive")
    , color_("color", "Primitive Color", tgt::Color(0.75f, 0.25f, 0.f, 1.f))
    , depthTest_("depthTest", "Depth Test", true)
    , pointSize_("pointSize", "Point Size", 3.f, 1.f, 20.f)
    , pointSmooth_("pointSmooth", "Point Smooth", false)
    , sphereDiameter_("sphereDiameter", "Sphere Diameter", 1.f, 0.001f, 20.f)
    , sphereSlicesStacks_("sphereSlicesTiles", "Sphere Slices/Tiles", 20, 5, 100)
    , geometryInport_(Port::INPORT, "geometry.input", "Geometry Input")
    , mesh_(nullptr)
    , sphereShader_(0)
    , pointList_(nullptr)
    , transformation_(tgt::mat4::identity)
{
    // Coordinate systems
    coordinateSystem_.addOption("world", "World coordinates");
    coordinateSystem_.addOption("viewport", "Viewport coordinates");
    coordinateSystem_.addOption("untransformed", "No transformation");

    // Rendering primitives
    renderingPrimitiveProp_.addOption("points", "Points");
    renderingPrimitiveProp_.addOption("line-strip", "Line Strip");
    //renderingPrimitiveProp_.addOption("line-loop", "Line Loop");
    renderingPrimitiveProp_.addOption("spheres", "Spheres");
    renderingPrimitiveProp_.addOption("illuminated-spheres", "Illuminated Spheres");

    renderingPrimitiveProp_.onChange(MemberFunctionCallback<PointListRenderer>(this, &PointListRenderer::invalidateGeometry));
    sphereDiameter_.onChange(MemberFunctionCallback<PointListRenderer>(this, &PointListRenderer::invalidateGeometry));
    sphereSlicesStacks_.onChange(MemberFunctionCallback<PointListRenderer>(this, &PointListRenderer::invalidateGeometry));

    pointSize_.setStepping(0.5f);
    pointSize_.setNumDecimals(1);
    sphereDiameter_.setStepping(0.001f);
    sphereDiameter_.setNumDecimals(3);

    addProperty(coordinateSystem_);
    addProperty(renderingPrimitiveProp_);
    addProperty(color_);
    addProperty(depthTest_);
    addProperty(pointSize_);
    addProperty(pointSmooth_);
    addProperty(sphereDiameter_);
    addProperty(sphereSlicesStacks_);

    addPort(geometryInport_);
    ON_CHANGE(geometryInport_, PointListRenderer, invalidateGeometry);
}

tgt::Bounds PointListRenderer::getBoundingBox() const {
    if (geometryInport_.hasData())
        return geometryInport_.getData()->getBoundingBox();

    return GeometryRendererBase::getBoundingBox();
}

void PointListRenderer::process() {
    tgtAssert(geometryInport_.isReady(), "inport not ready");
}

void PointListRenderer::initialize() {
    GeometryRendererBase::initialize();
    sphereShader_ = ShdrMgr.load("pointlistrenderer_spheres", "", false);
    circleShader_ = ShdrMgr.load("pointlistrenderer_circles", "", false);
}

void PointListRenderer::deinitialize() {
    ShdrMgr.dispose(sphereShader_);
    ShdrMgr.dispose(circleShader_);
    mesh_.reset();
    GeometryRendererBase::deinitialize();
}

void PointListRenderer::render() {
    tgtAssert(geometryInport_.isReady(), "inport not ready");

    // abort on invalid geometry
    if (!pointList_)
        return;

    // set matrix stack and OpenGL state
    if (coordinateSystem_.get() != "world") {
        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.pushMatrix();
        MatStack.loadIdentity();
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.pushMatrix();
        MatStack.loadIdentity();
        if (coordinateSystem_.get() == "viewport") {
            MatStack.translate(-1.f, -1.f, 0.f);
            MatStack.scale(2.f/viewport_.x, 2.f/viewport_.y, 1.f);
        }
    }

    if (!depthTest_.get())
        glDisable(GL_DEPTH_TEST);

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.multMatrix(transformation_);
  
    // render according to mode
    if (renderingPrimitiveProp_.get() == "spheres" || renderingPrimitiveProp_.get() == "illuminated-spheres") { // render spheres

        // activate shader
        sphereShader_->activate();

        // set matrix stack uniforms
        sphereShader_->setUniform("modelViewMatrixStack_", MatStack.getModelViewMatrix());
        sphereShader_->setUniform("modelViewProjectionMatrixStack_", MatStack.getProjectionMatrix() * MatStack.getModelViewMatrix());
        tgt::mat4 viewInverse;
        MatStack.getModelViewMatrix().invert(viewInverse);
        sphereShader_->setUniform("modelViewMatrixInverseStack_", viewInverse);

        // set lighting parameters
        sphereShader_->setUniform("lightingEnabled_", renderingPrimitiveProp_.get() == "illuminated-spheres");
        // set light source
        std::string prefix = "lightSource_.";
        sphereShader_->setUniform(prefix + "position_", tgt::vec4(0.f, 1.f * sphereDiameter_.get(), 3.f * sphereDiameter_.get(), 1.f));
        sphereShader_->setUniform(prefix + "ambientColor_", tgt::vec3(0.3f, 0.3f, 0.3f));
        sphereShader_->setUniform(prefix + "diffuseColor_", tgt::vec3(0.6f, 0.6f, 0.6f));
        sphereShader_->setUniform(prefix + "specularColor_", tgt::vec3(0.6f, 0.6f, 0.6f));
        //sphereShader_->setUniform(prefix+"attenuation_", tgt::vec3(0.f, 0.f, 0.f));

        // material
        prefix = "material_.";
        sphereShader_->setUniform(prefix + "ambientColor_", color_.get().xyz());
        sphereShader_->setUniform(prefix + "diffuseColor_", color_.get().xyz());
        sphereShader_->setUniform(prefix + "specularColor_", color_.get().xyz());
        sphereShader_->setUniform(prefix + "shininess_", 75.f);

        sphereShader_->setUniform("color_", color_.get());

        // now render all of the points using the sphere mesh
        for (size_t i = 0; i<pointList_->size(); ++i) {
            // translate the sphere to its location
            sphereShader_->setUniform("translate_", (*pointList_)[i]);

            mesh_->render();
        }
        sphereShader_->deactivate();
    }
    else if (renderingPrimitiveProp_.get() == "points") {
        glPointSize(pointSize_.get());
        if (pointSmooth_.get()) {
#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
            glEnable(GL_POINT_SPRITE);
#endif
            circleShader_->activate();
            circleShader_->setUniform("modelViewProjectionMatrixStack_", MatStack.getProjectionMatrix() * MatStack.getModelViewMatrix());
            circleShader_->setUniform("color_", color_.get());
        }

        // Render Mesh.
        mesh_->render();

        if (pointSmooth_.get()) {
            circleShader_->deactivate();
#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
            glDisable(GL_POINT_SPRITE);
#endif
        }
        glPointSize(1.f);

    }
    else { // line-strip or line-loop
        glLineWidth(pointSize_.get());
        if (pointSmooth_.get())
            glEnable(GL_LINE_SMOOTH);

        // Render Mesh.
        mesh_->render();

        glLineWidth(1.f);
        glDisable(GL_LINE_SMOOTH);
    }
    

    MatStack.popMatrix();

    LGL_ERROR;

    // reset matrix stack and OpenGL state
    if (coordinateSystem_.get() != "world") {
        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.popMatrix();
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.popMatrix();
    }

    glEnable(GL_DEPTH_TEST);
}

void PointListRenderer::invalidateGeometry() {

    if (!isInitialized())
        return;

    // Delete old mesh.
    mesh_.reset();
    pointList_ = nullptr;

    const PointListGeometry<tgt::vec3>* pointListGeom = dynamic_cast<const PointListGeometry<tgt::vec3>* >(geometryInport_.getData());
    const PointSegmentListGeometry<tgt::vec3>* segmentListGeom = dynamic_cast<const PointSegmentListGeometry<tgt::vec3>* >(geometryInport_.getData());

    if (pointListGeom) {
        transformation_ = pointListGeom->getTransformationMatrix();
        pointList_ = &(pointListGeom->getData());
    }
    else if (segmentListGeom) {
        transformation_ = segmentListGeom->getTransformationMatrix();
        pointList_ = &(segmentListGeom->getPoints());
    }
    else {
        LWARNING("No valid point list geometry!");
        return;
    }

    if (renderingPrimitiveProp_.get() == "spheres" || renderingPrimitiveProp_.get() == "illuminated-spheres") {
        GlMeshGeometryUInt16ColorNormal* sphereMesh = new GlMeshGeometryUInt16ColorNormal;
        sphereMesh->setPrimitiveType(GL_TRIANGLE_STRIP);
        sphereMesh->setSphereGeometry(sphereDiameter_.get(), tgt::vec3::zero, color_.get(), static_cast<size_t>(sphereSlicesStacks_.get()));
        mesh_.reset(sphereMesh);
    }
    else {

        // We need position and color, only.
        GlMeshGeometryUInt16Color* primitiveMesh = new GlMeshGeometryUInt16Color;

        // Fill mesh.
        for (size_t i = 0; i < pointList_->size(); ++i)
            primitiveMesh->addVertex(VertexColor((*pointList_)[i], color_.get()));

        if (renderingPrimitiveProp_.get() == "line-strip")
            primitiveMesh->setPrimitiveType(GL_LINE_STRIP);
        else if (renderingPrimitiveProp_.get() == "line-loop")
            primitiveMesh->setPrimitiveType(GL_LINE_LOOP);
        else if (renderingPrimitiveProp_.get() == "points")
            primitiveMesh->setPrimitiveType(GL_POINTS);
        else
            tgtAssert(false, "Unknown primitive type");

        mesh_.reset(primitiveMesh);
    }
}

} // namespace voreen
