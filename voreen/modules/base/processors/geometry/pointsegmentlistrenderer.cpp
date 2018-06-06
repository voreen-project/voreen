/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "pointsegmentlistrenderer.h"

#include "voreen/core/datastructures/geometry/pointlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"

#include "tgt/immediatemode/immediatemode.h"

namespace voreen {

const std::string PointSegmentListRenderer::loggerCat_("voreen.PointSegmentListRenderer");

PointSegmentListRenderer::PointSegmentListRenderer()
    : GeometryRendererBase()
    , coordinateSystem_("coordinateSystem", "Coordinate System")
    , renderingPrimitiveProp_("renderingPrimitive", "Rendering Primitive")
    , applyUniformColor_("applyUniformColor", "Apply Uniform Color", false)
    , color_("color", "Primitive Color", tgt::Color(0.75f, 0.25f, 0.f, 1.f))
    , depthTest_("depthTest", "Depth Test", true)
    , pointSize_("pointSize", "Point Size", 3.f, 1.f, 20.f)
    , pointSmooth_("pointSmooth", "Point Smooth", false)
    , sphereDiameter_("sphereDiameter", "Sphere Diameter", 0.02f, 0.01f, 20.f)
    , sphereSlicesStacks_("sphereSlicesStacks", "Sphere Slices/Stacks", 20, 5, 100)
    , segmentList_("segmentList", "Segment list")
    , segmentColor_("segmentColor", "Segment color")
    , segmentEnabled_("segmentEnabled", "Segment Enabled")
    , geometryInport_(Port::INPORT, "geometry.input", "Geometry Input")
    , sphereShader_(0)
    , mesh_(nullptr)
    , segmentListGeom_(nullptr)
{
    // Coordinate systems
    coordinateSystem_.addOption("world", "World coordinates");
    coordinateSystem_.addOption("viewport", "Viewport coordinates");
    coordinateSystem_.addOption("untransformed", "No transformation");

    // Rendering primitives
    renderingPrimitiveProp_.addOption("points", "Points");
    renderingPrimitiveProp_.addOption("line-strip", "Line Strip");
    renderingPrimitiveProp_.addOption("line-loop", "Line Loop");
    renderingPrimitiveProp_.addOption("spheres", "Spheres");
    renderingPrimitiveProp_.addOption("illuminated-spheres", "Illuminated Spheres");

    // display list invalidation callback
    applyUniformColor_.onChange(MemberFunctionCallback<PointSegmentListRenderer>(this, &PointSegmentListRenderer::invalidateGeometry));
    renderingPrimitiveProp_.onChange(MemberFunctionCallback<PointSegmentListRenderer>(this, &PointSegmentListRenderer::invalidateGeometry));
    sphereDiameter_.onChange(MemberFunctionCallback<PointSegmentListRenderer>(this, &PointSegmentListRenderer::invalidateGeometry));
    sphereSlicesStacks_.onChange(MemberFunctionCallback<PointSegmentListRenderer>(this, &PointSegmentListRenderer::invalidateGeometry));

    pointSize_.setStepping(0.5f);
    pointSize_.setNumDecimals(1);
    sphereDiameter_.setStepping(0.001f);
    sphereDiameter_.setNumDecimals(3);

    addProperty(coordinateSystem_);
    addProperty(renderingPrimitiveProp_);
    addProperty(applyUniformColor_);
    addProperty(color_);
    addProperty(depthTest_);
    addProperty(pointSize_);
    addProperty(pointSmooth_);
    addProperty(sphereDiameter_);
    addProperty(sphereSlicesStacks_);

    addProperty(segmentList_);
        segmentList_.setGroupID("segment");
    addProperty(segmentEnabled_);
        segmentEnabled_.setGroupID("segment");
    addProperty(segmentColor_);
        segmentColor_.setGroupID("segment");
    setPropertyGroupGuiName("segment", "Segment configuration");

    addPort(geometryInport_);
    ON_CHANGE(geometryInport_, PointSegmentListRenderer, invalidateGeometry);

    ON_CHANGE_LAMBDA(segmentList_, [this]{
        updateSegmentUI();
    });

    ON_CHANGE_LAMBDA(segmentColor_, [this]{
        int selectedIndex = segmentList_.getSelectedIndex();
        if (selectedIndex != -1 && static_cast<size_t>(selectedIndex) < segments_.size()){
            Segment &seg = segments_[selectedIndex];
            if (seg.color_ != segmentColor_.get()){
                seg.color_ = segmentColor_.get();
                invalidateGeometry();
                updateSegmentListUI();
            }
        }
    });

    ON_CHANGE_LAMBDA(segmentEnabled_, [this]{
        int selectedIndex = segmentList_.getSelectedIndex();
        if (selectedIndex != -1 && static_cast<size_t>(selectedIndex) < segments_.size()){
            Segment &seg = segments_[selectedIndex];
            if (seg.enabled_ != segmentEnabled_.get()){
                seg.enabled_ = segmentEnabled_.get();
                invalidateGeometry();
                updateSegmentListUI();
            }
        }
    });

    updateSegmentUI();
}

PointSegmentListRenderer::~PointSegmentListRenderer() {
}

void PointSegmentListRenderer::process() {
    tgtAssert(geometryInport_.isReady(), "inport not ready");
}

tgt::Bounds PointSegmentListRenderer::getBoundingBox() const {
    if (geometryInport_.hasData())
        return geometryInport_.getData()->getBoundingBox();

    return GeometryRendererBase::getBoundingBox();
}

void PointSegmentListRenderer::initialize() {
    GeometryRendererBase::initialize();
    sphereShader_ = ShdrMgr.load("pointlistrenderer_spheres", "", false);
    circleShader_ = ShdrMgr.load("pointlistrenderer_circles", "", false);
}

void PointSegmentListRenderer::deinitialize() {
    ShdrMgr.dispose(sphereShader_);
    ShdrMgr.dispose(circleShader_);
    mesh_.reset();
    GeometryRendererBase::deinitialize();
}

void PointSegmentListRenderer::render() {

    tgtAssert(geometryInport_.isReady(), "inport not ready");

    // abort on invalid geometry
    if (!segmentListGeom_)
        return;

    // set geometry transformation matrix and get point list
    tgt::mat4 m = segmentListGeom_->getTransformationMatrix();
    const std::vector<std::vector<tgt::vec3> >& segmentList = segmentListGeom_->getData();

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
    MatStack.multMatrix(m);

    // render according to mode
    if(renderingPrimitiveProp_.get() == "spheres" || renderingPrimitiveProp_.get() == "illuminated-spheres") { // render spheres

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
        sphereShader_->setUniform(prefix+"position_", tgt::vec4(0.f, 1.f * sphereDiameter_.get(), 3.f * sphereDiameter_.get(), 1.f));
        sphereShader_->setUniform(prefix+"ambientColor_", tgt::vec3(0.3f, 0.3f, 0.3f));
        sphereShader_->setUniform(prefix+"diffuseColor_", tgt::vec3(0.6f, 0.6f, 0.6f));
        sphereShader_->setUniform(prefix+"specularColor_", tgt::vec3(0.6f, 0.6f, 0.6f));
        //sphereShader_->setUniform(prefix+"attenuation_", tgt::vec3(0.f, 0.f, 0.f));

        // material
        prefix = "material_.";
        sphereShader_->setUniform(prefix+"ambientColor_", color_.get().xyz());
        sphereShader_->setUniform(prefix+"diffuseColor_", color_.get().xyz());
        sphereShader_->setUniform(prefix+"specularColor_", color_.get().xyz());
        sphereShader_->setUniform(prefix+"shininess_", 75.f);

        sphereShader_->setUniform("color_", color_.get());

        // iterate over segments
        for (size_t i=0; i < segmentList.size(); ++i) {
            Segment seg = segments_.at(i);
            if (!seg.enabled_) continue;
            // apply segment color
            if (!applyUniformColor_.get()) {
                prefix = "material_.";
                sphereShader_->setUniform(prefix+"ambientColor_", seg.color_.xyz());
                sphereShader_->setUniform(prefix+"diffuseColor_", seg.color_.xyz());
                sphereShader_->setUniform(prefix+"specularColor_", seg.color_.xyz());
                sphereShader_->setUniform(prefix+"shininess_", 75.f);

                sphereShader_->setUniform("color_", seg.color_);
            }

            // now render all of the points using the sphere mesh
            for (size_t j=0; j<segmentList[i].size(); ++j) {
                // translate the sphere to its location
                sphereShader_->setUniform("translate_", segmentList[i][j]);

                mesh_->render();
            }
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

        // Render segments.
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

        // Render segments.
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

void PointSegmentListRenderer::invalidateGeometry() {

    if (!isInitialized())
        return;

    // Delete old mesh.
    mesh_.reset();

    segmentListGeom_ = dynamic_cast<const PointSegmentListGeometry<tgt::vec3>* >(geometryInport_.getData());
    if (!segmentListGeom_) {
        LWARNING("Invalid geometry. PointSegmentListGeometry<vec3> expected.");
        return;
    }

    // Update segment list.
    updateSegmentList();

    if (renderingPrimitiveProp_.get() == "spheres" || renderingPrimitiveProp_.get() == "illuminated-spheres") {
        GlMeshGeometryUInt16ColorNormal* sphereMesh = new GlMeshGeometryUInt16ColorNormal;
        sphereMesh->setPrimitiveType(GL_TRIANGLE_STRIP);
        sphereMesh->setSphereGeometry(sphereDiameter_.get(), tgt::vec3::zero, color_.get(), static_cast<size_t>(sphereSlicesStacks_.get()));
        mesh_.reset(sphereMesh);
    }
    else {

        // We need position and color, only.
        GlMeshGeometryUInt32Color* primitiveMesh = new GlMeshGeometryUInt32Color;
        primitiveMesh->enablePrimitiveRestart(); // By default, enable primitive restart.

        const std::vector<std::vector<tgt::vec3> >& segmentList = segmentListGeom_->getData();
        for (size_t i = 0; i < segmentList.size(); ++i) {
            Segment seg = segments_.at(i);
            if (!seg.enabled_) continue;
            
            // Extract color.
            tgt::vec4 color = applyUniformColor_.get() ? color_.get() : seg.color_;

            // Fill mesh.
            if (renderingPrimitiveProp_.get() != "points") { // line-strip, line-loop
                for (size_t j = 0; j < segmentList[i].size(); ++j) {
                    primitiveMesh->addIndex(static_cast<uint32_t>(primitiveMesh->getNumVertices()));
                    primitiveMesh->addVertex(VertexColor(segmentList[i][j], color));
                }
                primitiveMesh->addIndex(static_cast<uint32_t>(primitiveMesh->getPrimitiveRestartIndex()));
            }
            else { // for points, draw arrays to save memory
                for (size_t j = 0; j < segmentList[i].size(); ++j)
                    primitiveMesh->addVertex(VertexColor(segmentList[i][j], color));
            }
        }

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

void PointSegmentListRenderer::serialize(Serializer& s) const
{
    GeometryRendererBase::serialize(s);
    s.serialize("segments", segments_);
}

void PointSegmentListRenderer::deserialize(Deserializer& s)
{
    GeometryRendererBase::deserialize(s);
    s.optionalDeserialize("segments", segments_, std::vector<Segment>());
}

void PointSegmentListRenderer::updateSegmentList()
{
    std::vector<tgt::vec4> colorMap;
    colorMap.push_back(tgt::vec4(255,0,0,255) / 255.f);
    colorMap.push_back(tgt::vec4(0,255,0,255) / 255.f);
    colorMap.push_back(tgt::vec4(0,0,255,255) / 255.f);
    colorMap.push_back(tgt::vec4(255,0,255,255) / 255.f);
    colorMap.push_back(tgt::vec4(0,255,255,255) / 255.f);
    colorMap.push_back(tgt::vec4(255,255,0,255) / 255.f);
    colorMap.push_back(tgt::vec4(255,100,20,255) / 255.f);
    colorMap.push_back(tgt::vec4(250,200,150,255) / 255.f);
    colorMap.push_back(tgt::vec4(150,200,250,255) / 255.f);
    //colorMap.push_back(tgt::vec4(30,30,30,255) / 255.f);
    colorMap.push_back(tgt::vec4(200,200,127,255) / 255.f);    

    int oldSegmentCount = static_cast<int>(segments_.size());
    int newSegmentCount = static_cast<int>(segmentListGeom_->getNumSegments());

    bool needstoupdateprops = false;

    if (oldSegmentCount > newSegmentCount){
        // just throw away unused segments
        segments_.resize(newSegmentCount);
    }
    else if (oldSegmentCount < newSegmentCount){
        segments_.resize(newSegmentCount);
        for(int i = oldSegmentCount; i != newSegmentCount; i++){
            segments_[i].enabled_ = true;
            segments_[i].color_ = colorMap[i%colorMap.size()];
        }
    }
    if (segmentList_.getKeys().size() != segments_.size())
        updateSegmentListUI();
}

void PointSegmentListRenderer::updateSegmentListUI()
{
    std::vector<std::string> keys = segmentList_.getKeys();
    while (keys.size() > segments_.size()){
        segmentList_.removeOption(keys.back());
        keys.pop_back();
    }


    for(int i = static_cast<int>(keys.size()); i <= static_cast<int>(segments_.size()); i++){
        std::stringstream sk;
        sk << i;
        segmentList_.addOption(sk.str(), "", i);
    }

    for(size_t i = 0; i != segments_.size(); i++){
        Segment seg = segments_[i];
        std::stringstream ss;
        ss << i << " (color: " << seg.color_ << ", " << (seg.enabled_ ? "enabled" : "disabled") <<")";
        std::stringstream sk;
        sk << i;
        segmentList_.setOptionDescription(sk.str(), ss.str());
    }

    updateSegmentUI();
    segmentList_.updateWidgets();
}

void PointSegmentListRenderer::updateSegmentUI()
{
    int selectedIndex = segmentList_.getSelectedIndex();
    bool exists = selectedIndex != -1 && static_cast<size_t>(selectedIndex) < segments_.size();
    segmentColor_.setReadOnlyFlag(!exists);
    segmentEnabled_.setReadOnlyFlag(!exists);
    if (exists){
        Segment seg = segments_[selectedIndex];
        segmentEnabled_.set(seg.enabled_);
        segmentColor_.set(seg.color_);
    }
}


void PointSegmentListRenderer::Segment::serialize(Serializer& s) const
{
    s.serialize("color", color_);
    s.serialize("enabled", enabled_);
}

void PointSegmentListRenderer::Segment::deserialize(Deserializer& s)
{
    s.deserialize("color", color_);
    s.deserialize("enabled", enabled_);
}

} // namespace voreen
