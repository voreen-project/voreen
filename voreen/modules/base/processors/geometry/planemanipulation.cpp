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

#include "planemanipulation.h"
#include "voreen/core/utils/voreenqualitymode.h"
#include "tgt/immediatemode/immediatemode.h"

//------------------------------------------------------------------------//
//  Helper Funktions and Constants                                        //
//------------------------------------------------------------------------//
namespace {
    const float EPSILON = 0.000001f;                                ///< Float accuracy in this processor
    const tgt::vec4 LIGHT_POSITION = tgt::vec4(1.f,1.f,-1.f,1.f);   ///< The used light position in camera space
    const tgt::vec3 AMBIENT_COLOR  = tgt::vec3(0.3f, 0.3f, 0.3f);   ///< Used ambient color
    const tgt::vec3 DIFFUSE_COLOR  = tgt::vec3(0.6f, 0.6f, 0.6f);   ///< Used diffuse color
    const tgt::vec3 SPECULAR_COLOR = tgt::vec3(0.6f, 0.6f, 0.6f);   ///< Used speculat color

    /** Updates the camera/matrix uniforms. (except translate_)*/
    void setMatrixShaderUniforms(tgt::Shader* shader) {
        tgtAssert(shader->isActivated(),"Shader is not activated");
        // set matrix stack uniforms
        shader->setUniform("modelViewMatrixStack_", MatStack.getModelViewMatrix());
        shader->setUniform("modelViewProjectionMatrixStack_", MatStack.getProjectionMatrix() * MatStack.getModelViewMatrix());
        tgt::mat4 viewInverse;
        MatStack.getModelViewMatrix().invert(viewInverse);
        shader->setUniform("modelViewMatrixInverseStack_", viewInverse);
    }

    /**
     * Updates the color/material uniforms.
     * @param enableLighting (Dis)Enables lighting in the shader.
     * @param color The current color. If lighting is enabled, it is used as material color, else as solid color.
     */
    void setColorShaderUniforms(tgt::Shader* shader, bool enableLighting, const tgt::vec4& color) {
        // set lighting parameters
        shader->setUniform("lightingEnabled_", enableLighting);
        if(enableLighting) {
            // set light source
            std::string prefix = "lightSource_.";
            shader->setUniform(prefix+"position_", LIGHT_POSITION);
            shader->setUniform(prefix+"ambientColor_", AMBIENT_COLOR);
            shader->setUniform(prefix+"diffuseColor_", DIFFUSE_COLOR);
            shader->setUniform(prefix+"specularColor_", SPECULAR_COLOR);
            // material
            prefix = "material_.";
            shader->setUniform(prefix+"ambientColor_", color);
            shader->setUniform(prefix+"diffuseColor_", color);
            shader->setUniform(prefix+"specularColor_", color);
            shader->setUniform(prefix+"shininess_", 25.f);
            // set solid color (default values, prevent undefined behavior)
            shader->setUniform("color_", tgt::vec4::zero);
        } else {
            // set light source (default values, prevent undefined behavior)
            std::string prefix = "lightSource_.";
            shader->setUniform(prefix + "position_", tgt::vec4::zero);
            shader->setUniform(prefix + "ambientColor_", tgt::vec3::zero);
            shader->setUniform(prefix + "diffuseColor_", tgt::vec3::zero);
            shader->setUniform(prefix + "specularColor_", tgt::vec3::zero);
            // material (default values, prevent undefined behavior)
            prefix = "material_.";
            shader->setUniform(prefix + "ambientColor_", tgt::vec4::zero);
            shader->setUniform(prefix + "diffuseColor_", tgt::vec4::zero);
            shader->setUniform(prefix + "specularColor_", tgt::vec4::zero);
            shader->setUniform(prefix + "shininess_", 0.0f);
            //set solid color
            shader->setUniform("color_", color);
        }
    }

    /**
     * Sets all shader defines. camera + color (except translate_).
     * @param shader The used shader. It must be activated before this call.
     * @param setMatrices If true, the camera Matrices are set according to the current MatrixStack.
     * @param enableLighting (Dis)Enables lighting in the shader.
     * @param color The current color. If lighting is enabled, it is used as material color, else as solid color.
     */
    void setAllShaderUniforms(tgt::Shader* shader, bool enableLighting, const tgt::vec4& color) {
        tgtAssert(shader->isActivated(),"Shader is not activated");
        // set matrix stack uniforms
        setMatrixShaderUniforms(shader);
        // set color parameters
        setColorShaderUniforms(shader,enableLighting,color);
    }

    /**
     * Return the given line vector equation scalar for the projection of the given point.
     *
     * @param point the Point
     * @param linePositionVector position vector of line vector equation
     * @param lineDirectionVector direction vector of line vector equation
     *
     * @return distance between point and line
     */
    inline float getPointLineProjectionScalar(const tgt::vec2& point, const tgt::vec2& linePositionVector, const tgt::vec2& lineDirectionVector) {
        return tgt::dot(point - linePositionVector, lineDirectionVector) / tgt::lengthSq(lineDirectionVector);
    }
    /**
     * Determines the intersection point of the given edge and a plane equation.
     *
     * If the intersection point is outside of the edge (ie. not between start and end vertex),
     * the intersection point is omitted.
     *
     * @param intersectionPoint the variable in which the intersection point should be stored
     * @param startVertex start vertex of the edge
     * @param endVertex end vertex of the edge
     * @param plane the plane equation
     *
     * @return @c true if an intersection point within the edge is found, otherwise @c false.
     */
    bool getIntersectionPointInBB(tgt::vec3& intersectionPoint, const tgt::vec3& startVertex,
                                  const tgt::vec3& endVertex, const tgt::vec4& plane) {
        tgt::vec3 lineDirectionVector = tgt::normalize(endVertex - startVertex);
        float denominator = tgt::dot(lineDirectionVector, plane.xyz());
        if (std::abs(denominator) < EPSILON)
            return false;
        float numerator = tgt::dot((plane.xyz() * plane.w) - startVertex, plane.xyz());
        float t = numerator / denominator;
        tgt::vec3 intersection = startVertex + t * lineDirectionVector;
        // Is intersection point within AABB?
        float tEnd = tgt::length(endVertex - startVertex);
        if (t >= 0 && t <= tEnd) {
            intersectionPoint = intersection;
            return true;
        }
        return false;
    }
}

//------------------------------------------------------------------------//
//------------------------------------------------------------------------//
//  PlanemMnipulation Processor                                           //
//------------------------------------------------------------------------//
//------------------------------------------------------------------------//
namespace voreen {

    const std::string PlaneManipulation::loggerCat_("voreen.base.PlaneManipulation");

PlaneManipulation::PlaneManipulation()
    : GeometryRendererBase()
    //ports
    , volumeInport_(Port::INPORT, "volume", "Volume Input")
    , geometryInport_(Port::INPORT, "geometry", "Geometry Input")
    //properties
        //to be linked
    , enable_("enable", "Enable", true)
    , invert_("invert", "Invert", false)
    , planeNormal_("planeNormal", "Plane Normal", tgt::vec3(0, 1, 0), tgt::vec3(-1), tgt::vec3(1))
    , planePosition_("planePosition", "Plane Position", 0.0f, -1000000.f, 1000000.f, Processor::INVALID_RESULT, NumericProperty<float>::STATIC)
        //general
    , generalColor_("planeColor", "Used Color", tgt::vec4(0.0f, 0.0f, 1.0f, 1.f))
    , resetPlane_("resetplane", "Reset Plane")
        //plane
    , renderPlane_("renderGeometry", "Show Plane", true)
    , lineWidth_("lineWidth", "Line Width", 2.0f, 1.0f, 10.0f)
    , grabbedPlaneOpacity_("grabbedPlaneOpacity", "Grabbed Handle Opacity", 0.2f, 0.f, 1.f)
    , alwaysUseOpacity_("alwaysUseOpacity","Always use Opacity",false)
    , renderAnchors_("renderAnchors","Render Anchors",true)
        //manipulator
    , renderManipulator_("renderManipulator", "Show Handle", true)
    , manipulatorScale_("ManipulatorScale","Handle Scale",1.0f,0.01f,2.0f)
    , grabbedElementColor_("grabbedElementColor", "Grabbed Handle Element Color", tgt::vec4(0.9f, 0.9f, 0.9f, 1.0f))
    , blockedElementColor_("blockedElementColor", "Blocked Handle Element Color", tgt::vec4(0.75f, 0.05f, 0.05f, 1.0f))
    , resetManipulatorPos_("resetmanipulatorpos", "Reset Handle on Plane")
    , autoReset_("autoreset", "Auto Reset Handle", false, Processor::VALID,Property::LOD_ADVANCED)
    , manipulatorCenter_("manipulatorCenter", "Center in World Space",tgt::vec3::zero,tgt::vec3(-1000000.f),tgt::vec3(1000000.f),Processor::INVALID_RESULT,NumericProperty<tgt::vec3>::STATIC,Property::LOD_ADVANCED)
        //debug
    , positionRange_("positionrange", "Position Range", tgt::vec2(-1.f, 1.f), tgt::vec2(-1000000.f), tgt::vec2(1000000.f),
                     Processor::VALID, NumericProperty<tgt::vec2>::STATIC, Property::LOD_DEBUG)
    , glMeshShaderTransparentProp_("glMeshShaderTransparentProp","Mesh Shader (transparent)","planemanipulation_mesh.frag","planemanipulation_mesh.vert","",Processor::INVALID_PROGRAM,Property::LOD_DEBUG)
    , glMeshShaderOpaqueProp_("glMeshShaderOpaqueProp","Mesh Shader (opaque)","planemanipulation_mesh.frag","planemanipulation_mesh.vert","",Processor::INVALID_PROGRAM,Property::LOD_DEBUG)
    , iModeShaderTransparentProp_("iModeShaderTransparentProp_","Shader (transparent)","planemanipulation_imode.frag","planemanipulation_imode.vert","",Processor::INVALID_PROGRAM,Property::LOD_DEBUG)
    , iModeShaderOpaqueProp_("iModeShaderOpaqueProp","Shader (opaque)","planemanipulation_imode.frag","planemanipulation_imode.vert","",Processor::INVALID_PROGRAM,Property::LOD_DEBUG)
    , moveEventProp_(0)
        //geometries
    , anchorGeometry_(0)
    , tipGeometry_(0)
    , topCylinderGeometry_(0)
    , bottomCylinderGeometry_(0)
        //sizes
    , manipulatorLength_(1.f)
    , manipulatorDiameter_(0.04f)
    , manipulatorTipRadius_(0.05f)
        //event handling flags
    , manipulatorTopTipIsGrabbed_(false)
    , manipulatorBottomTipIsGrabbed_(false)
    , manipulatorCylinderIsGrabbed_(false)
    , wholeManipulatorIsGrabbed_(false)
    , invalidMovement_(false)
        //helper
    , actualPlane_(0.f) //stored as (normal.xyz,distPosition)
{
    //ports
    addPort(volumeInport_);
    addPort(geometryInport_);
    //properties
        //to be linked
    addProperty(enable_);
        enable_.setGroupID("link");
    addProperty(invert_);
        invert_.setGroupID("link");
    addProperty(planeNormal_);
        planeNormal_.onChange(MemberFunctionCallback<PlaneManipulation>(this, &PlaneManipulation::setManipulatorToStdPos));
        planeNormal_.setReadOnlyFlag(true);
        planeNormal_.setGroupID("link");
    addProperty(planePosition_);
        planePosition_.onChange(MemberFunctionCallback<PlaneManipulation>(this, &PlaneManipulation::setManipulatorToStdPos));
        planePosition_.setReadOnlyFlag(true);
        planePosition_.setGroupID("link");
    setPropertyGroupGuiName("link","Link with GeometryClippingProcessor");
        //general
    addProperty(generalColor_);
        generalColor_.setGroupID("general");
    addProperty(resetPlane_);
        resetPlane_.onClick(MemberFunctionCallback<PlaneManipulation>(this, &PlaneManipulation::resetPlane));
        resetPlane_.setGroupID("general");
    setPropertyGroupGuiName("general","General Settings");
        //plane
    addProperty(renderPlane_);
        renderPlane_.onChange(MemberFunctionCallback<PlaneManipulation>(this, &PlaneManipulation::updatePropertyVisibility));
        renderPlane_.setGroupID("plane");
    addProperty(lineWidth_);
        lineWidth_.setGroupID("plane");
    addProperty(grabbedPlaneOpacity_);
        grabbedPlaneOpacity_.setGroupID("plane");
    addProperty(alwaysUseOpacity_);
        alwaysUseOpacity_.setGroupID("plane");
    addProperty(renderAnchors_);
        renderAnchors_.setGroupID("plane");
    addProperty(manipulatorCenter_);
        manipulatorCenter_.setGroupID("plane");
        manipulatorCenter_.setNumDecimals(4);
        manipulatorCenter_.setReadOnlyFlag(true);
        manipulatorCenter_.onChange(MemberFunctionOneParameterCallback<PlaneManipulation,bool>(this, &PlaneManipulation::updateManipulatorElements,true));
    setPropertyGroupGuiName("plane","Plane Settings");
        //manipulator
    addProperty(renderManipulator_);
        renderManipulator_.onChange(MemberFunctionCallback<PlaneManipulation>(this, &PlaneManipulation::updatePropertyVisibility));
        renderManipulator_.setGroupID("mani");
    addProperty(manipulatorScale_);
        manipulatorScale_.onChange(MemberFunctionCallback<PlaneManipulation>(this, &PlaneManipulation::updateManipulatorScale));
        manipulatorScale_.setGroupID("mani");
    addProperty(grabbedElementColor_);
        grabbedElementColor_.setGroupID("mani");
    addProperty(blockedElementColor_);
        blockedElementColor_.setGroupID("mani");
    addProperty(resetManipulatorPos_);
        resetManipulatorPos_.onClick(MemberFunctionCallback<PlaneManipulation>(this, &PlaneManipulation::setManipulatorToStdPos));
        resetManipulatorPos_.setGroupID("mani");
    addProperty(autoReset_);
        autoReset_.setGroupID("mani");
    setPropertyGroupGuiName("mani","Manipulator Settings");
        //debug
    addProperty(positionRange_);
        positionRange_.setReadOnlyFlag(true);
    addProperty(glMeshShaderTransparentProp_);
    addProperty(glMeshShaderOpaqueProp_);
    addProperty(iModeShaderTransparentProp_);
    addProperty(iModeShaderOpaqueProp_);
    moveEventProp_ = new EventProperty<PlaneManipulation>(
        "mouseEvent.clipplaneManipulation", "Clipplane manipulation", this,
        &PlaneManipulation::planeManipulation,
        tgt::MouseEvent::MOUSE_BUTTON_ALL,
        tgt::MouseEvent::MOTION | tgt::MouseEvent::PRESSED | tgt::MouseEvent::RELEASED,
        tgt::Event::MODIFIER_NONE, false);
    addEventProperty(moveEventProp_);

    //update default visibility
    updatePropertyVisibility();
}

PlaneManipulation::~PlaneManipulation() {
    delete moveEventProp_;
}

Processor* PlaneManipulation::create() const {
    return new PlaneManipulation();
}

bool PlaneManipulation::isReady() const {
    return (volumeInport_.isReady() || geometryInport_.isReady());
}

void PlaneManipulation::initialize() {
    GeometryRendererBase::initialize();

    glMeshShaderTransparentProp_.setHeader("#define USE_TRANSPARENCY\n");
    glMeshShaderTransparentProp_.rebuild();

    glMeshShaderOpaqueProp_.rebuild();

    iModeShaderTransparentProp_.setHeader("#define USE_TRANSPARENCY\n");
    iModeShaderTransparentProp_.rebuild();

    iModeShaderOpaqueProp_.rebuild();
}

void PlaneManipulation::deinitialize() {
    //delete all geometries
    delete anchorGeometry_; delete tipGeometry_;
    delete topCylinderGeometry_; delete bottomCylinderGeometry_;
    GeometryRendererBase::deinitialize();
}

void PlaneManipulation::invalidate(int inv) {
    GeometryRendererBase::invalidate(inv);

    // check for initialization is necessary since otherwise the shader programs are rebuilt after deserialization
    // when the network is deserialized and the port connections are removed
    if (!isInitialized())
        return;

    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM) {
        if(!glMeshShaderTransparentProp_.getHeader().empty()) {
            glMeshShaderTransparentProp_.rebuild();
        }
        if(!glMeshShaderOpaqueProp_.getHeader().empty()) {
            glMeshShaderOpaqueProp_.rebuild();
        }
        if(!iModeShaderTransparentProp_.getHeader().empty()) {
            iModeShaderTransparentProp_.rebuild();
        }
        if(!iModeShaderOpaqueProp_.getHeader().empty()) {
            iModeShaderOpaqueProp_.rebuild();
        }
    }
}

tgt::Bounds PlaneManipulation::getBoundingBox() const {
    return worldSceneBounds_;
}

void PlaneManipulation::setIDManager(IDManager* idm) {
    if (idManager_ == idm)
        return;

    idManager_ = idm;
    if (idManager_) {
        //set objects: plane manipulator
        idManager_->registerObject(&manipulatorCenter_);
        idManager_->registerObject(&topManipulatorPos_);
        idManager_->registerObject(&bottomManipulatorPos_);
    }
}

void PlaneManipulation::renderPicking() {
    //Do nothing, if disabled
    tgt::Shader* meshShader = glMeshShaderOpaqueProp_.getShader();

    if (!idManager_ || !enable_.get() || !meshShader)
        return;

    meshShader->activate();

    //update camrea and color
    setAllShaderUniforms(meshShader,false,idManager_->getColorFromObjectFloat(&topManipulatorPos_));
    paintManipulatorTip(topManipulatorPos_, meshShader);

    //update only color (camera is already up to date)
    setColorShaderUniforms(meshShader,false,idManager_->getColorFromObjectFloat(&bottomManipulatorPos_));
    paintManipulatorTip(bottomManipulatorPos_, meshShader);

    //update only color (camera will be reset by paintManipulatorCylinderTop/BottomHalf())
    setColorShaderUniforms(meshShader,false,idManager_->getColorFromObjectFloat(&manipulatorCenter_));
    paintManipulatorCylinder(manipulatorCenter_.get(), manipulatorOrientation_, meshShader);

    meshShader->deactivate();
}

void PlaneManipulation::render() {
    render(false);
}
void PlaneManipulation::renderTransparent() {
    render(true);
}
bool PlaneManipulation::usesTransparency() const {
    return true;
}
void PlaneManipulation::render(bool useTransparency) {
    //Do nothing, if disabled
    if (!enable_.get())
        return;

    // check if at least one input is connected
    tgtAssert(volumeInport_.isReady() || geometryInport_.isReady(), "neither inport is ready");

    // check if bounding box changed
    checkInportsForBoundingBox();

    // render clipping plane
    if (renderPlane_.get()) {
        renderPlane(useTransparency);
    }

    // render manipulator
    if (renderManipulator_.get()) {
        paintManipulator(useTransparency);
    }
}

//------------------------------------------------------------------------//
//  Render Functions                                                      //
//------------------------------------------------------------------------//
void PlaneManipulation::renderPlane(bool useTransparency) {
    // Update anchors and lines, if plane changed...
    if (tgt::equal(actualPlane_, tgt::vec4(planeNormal_.get(), planePosition_.get())) != tgt::bvec4(true, true, true, true))
        updateAnchorsAndLines();

    tgt::Shader* iModeShader = useTransparency ? iModeShaderTransparentProp_.getShader() : iModeShaderOpaqueProp_.getShader();
    tgtAssert(iModeShader, "noShader");

    iModeShader->activate();
    iModeShader->setUniform("headPointerImage_", 0); // unued, but set default value to prevent undefined behavior
    iModeShader->setUniform("modelViewProjectionMatrixStack_", MatStack.getProjectionMatrix() * MatStack.getModelViewMatrix());
    // Render lines...
    iModeShader->setUniform("color_", generalColor_.get());
    glLineWidth(lineWidth_.get());
    IMode.begin(tgt::ImmediateMode::LINES);
        for (std::vector<Line>::iterator it = lines_.begin(); it != lines_.end(); ++it) {
            IMode.vertex(*it->first);
            IMode.vertex(*it->second);
        }
    IMode.end();
    glLineWidth(1.f);

    //check if plane should be rendered
    if (shouldRenderPlane()) {

        //prevent z fighting by using polygon offset
        glEnable(GL_POLYGON_OFFSET_FILL);

        float sign = invert_.get() ? -1.f : 1.f;
        float cos = sign * tgt::dot(manipulatorOrientation_, tgt::normalize(manipulatorCenter_.get() - camera_.getPositionWithOffsets()));
        if (cos < 0.f)
            glPolygonOffset(-2.f, -2.f);
        else
            glPolygonOffset(2.f, 2.f);

        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);// GL_DST_ALPHA);

        //render plane polygon
        iModeShader->setUniform("color_", tgt::vec4(generalColor_.get().xyz(), grabbedPlaneOpacity_.get()));
        IMode.begin(tgt::ImmediateMode::TRIANGLES);
            for (std::vector<Line>::iterator it = lines_.begin(); it != lines_.end(); ++it) {
                IMode.vertex(*lines_.begin()->first);
                IMode.vertex(*it->first);
                IMode.vertex(*it->second);
            }
        IMode.end();

        glBlendFunc(GL_ONE, GL_ZERO);
        glDisable(GL_BLEND);
        glPolygonOffset(0.f,0.f);
        glDisable(GL_POLYGON_OFFSET_FILL);
    }
    iModeShader->deactivate();

    tgt::Shader* meshShader = useTransparency ? glMeshShaderTransparentProp_.getShader() : glMeshShaderOpaqueProp_.getShader();

    //render anchors
    if (renderAnchors_.get() && meshShader) {

        meshShader->activate();
        //set all uniforms (camera + color)
        setAllShaderUniforms(meshShader, true, tgt::vec4(generalColor_.get().xyz(),1.f));
        meshShader->setUniform("headPointerImage_", 0); // unued, but set default value to prevent undefined behavior

        //update geometry lazy
        if(!anchorGeometry_) {
            anchorGeometry_ = new GlMeshGeometryUInt16Normal();
            anchorGeometry_->setSphereGeometry(manipulatorTipRadius_ / 2.5f,tgt::vec3::zero,generalColor_.get(),20); //color is not used (shader-material)
        }

        //render
        for (std::vector<Anchor>::iterator it = anchors_.begin(); it != anchors_.end(); ++it) {
            meshShader->setUniform("translate_", *it);
            anchorGeometry_->render();
        }

        meshShader->deactivate();
    }
}

void PlaneManipulation::paintManipulator(bool useTransparency) {
    //activate shader
    tgt::Shader* meshShader = useTransparency ? glMeshShaderTransparentProp_.getShader() : glMeshShaderOpaqueProp_.getShader();
    tgtAssert(meshShader, "no shader");
    meshShader->activate();

    tgt::vec4 renderColor;          //color that is determined below and is used for rendering

    //get cosines for determining opacity etc.
    float topCos = tgt::dot(manipulatorOrientation_, tgt::normalize(camera_.getPositionWithOffsets() - manipulatorCenter_.get()));
    float bottomCos = tgt::dot(-manipulatorOrientation_, tgt::normalize(camera_.getPositionWithOffsets() - manipulatorCenter_.get()));

    //render cylinder
    if (manipulatorCylinderIsGrabbed_ || wholeManipulatorIsGrabbed_) {
        if (invalidMovement_)
            renderColor = blockedElementColor_.get();
        else
            renderColor = grabbedElementColor_.get();
        renderColor.a = 1.f;

        //camera set by paintManipulatorTop/BottomHalf()
        setColorShaderUniforms(meshShader,true,renderColor);
        paintManipulatorCylinder(manipulatorCenter_.get(), manipulatorOrientation_, meshShader);
    }
    else {
        //if not grabbed: the half of the cylinder behind the plane should be rendered transparent, depending on the angle with the camera
        renderColor = generalColor_.get();
        renderColor.a = 1.f;

        if (topCos < -EPSILON * 10.f)
            renderColor.a = 1.f + std::max(-(topCos * topCos) - 0.075f, -0.96f);

        //camera set by paintManipulatorTopHalf()
        setColorShaderUniforms(meshShader,true,renderColor);
        paintManipulatorCylinderTopHalf(manipulatorCenter_.get(), manipulatorOrientation_, meshShader);

        renderColor.a = 1.f;

        if (bottomCos < -EPSILON * 10.f)
            renderColor.a = 1.f + std::max(-(bottomCos * bottomCos) - 0.075f, -0.96f);

        //camera set by paintManipulatorBottomHalf()
        setColorShaderUniforms(meshShader,true,renderColor);
        paintManipulatorCylinderBottomHalf(manipulatorCenter_.get(), manipulatorOrientation_, meshShader);
    }

    //render top tip
    if (manipulatorTopTipIsGrabbed_ || wholeManipulatorIsGrabbed_)
        if (invalidMovement_)
            renderColor = blockedElementColor_.get();
        else
            renderColor = grabbedElementColor_.get();
    else
        renderColor = generalColor_.get();

    //if not grabbed: render transparent if behind the plane (depending on angle with camera)
    renderColor.a = 1.f;
    if (topCos < -EPSILON * 10.f)
            renderColor.a = 1.f + std::max(-(topCos * topCos) - 0.075f, -0.96f);

    //Update camera + color
    setAllShaderUniforms(meshShader,true,renderColor);
    paintManipulatorTip(topManipulatorPos_, meshShader);

    //render bottom tip
    if (manipulatorBottomTipIsGrabbed_ || wholeManipulatorIsGrabbed_)
        if (invalidMovement_)
            renderColor = blockedElementColor_.get();
        else
            renderColor = grabbedElementColor_.get();
    else
        renderColor = generalColor_.get();

    //if not grabbed: render transparent if behind the plane (depending on angle with camera)
    renderColor.a = 1.f;
    if (bottomCos < -EPSILON * 10.f)
            renderColor.a = 1.f + std::max(-(bottomCos * bottomCos) - 0.075f, -0.96f);

    setColorShaderUniforms(meshShader,true,renderColor); //camera already set by top
    paintManipulatorTip(bottomManipulatorPos_, meshShader);

    meshShader->deactivate();
}

//------------------------------------------------------------------------//
//  Anchors and Lines Functions                                           //
//------------------------------------------------------------------------//
void PlaneManipulation::updateAnchorsAndLines() {
    tgtAssert((volumeInport_.getData() || geometryInport_.getData()), "no input volume or geometry");

    anchors_.clear();
    lines_.clear();

    //use physical bounding box and transform its vertices to world coordinates to compute the anchors and lines
    tgt::Bounds bBox = physicalSceneBounds_;

    tgt::vec3 physicalLLF = bBox.getLLF();
    tgt::vec3 physicalURB = bBox.getURB();

    tgt::vec3 worldLLF = (physicalToWorld_ * tgt::vec4(physicalLLF, 1.f)).xyz();
    tgt::vec3 worldURB = (physicalToWorld_ * tgt::vec4(physicalURB, 1.f)).xyz();
    tgt::vec3 worldLRF = (physicalToWorld_ * tgt::vec4(physicalURB.x, physicalLLF.y, physicalLLF.z, 1.f)).xyz();
    tgt::vec3 worldULF = (physicalToWorld_ * tgt::vec4(physicalLLF.x, physicalURB.y, physicalLLF.z, 1.f)).xyz();
    tgt::vec3 worldLLB = (physicalToWorld_ * tgt::vec4(physicalLLF.x, physicalLLF.y, physicalURB.z, 1.f)).xyz();
    tgt::vec3 worldURF = (physicalToWorld_ * tgt::vec4(physicalURB.x, physicalURB.y, physicalLLF.z, 1.f)).xyz();
    tgt::vec3 worldULB = (physicalToWorld_ * tgt::vec4(physicalLLF.x, physicalURB.y, physicalURB.z, 1.f)).xyz();
    tgt::vec3 worldLRB = (physicalToWorld_ * tgt::vec4(physicalURB.x, physicalLLF.y, physicalURB.z, 1.f)).xyz();

    //get plane
    tgt::vec4 plane = tgt::vec4(tgt::normalize(planeNormal_.get()), planePosition_.get());

    // Determine and add all anchors...
    addAnchor(worldLLF, worldLRF, plane);
    addAnchor(worldLLF, worldULF, plane);
    addAnchor(worldLLF, worldLLB, plane);
    addAnchor(worldURB, worldULB, plane);
    addAnchor(worldURB, worldLRB, plane);
    addAnchor(worldURB, worldURF, plane);
    addAnchor(worldLLB, worldLRB, plane);
    addAnchor(worldLLB, worldULB, plane);
    addAnchor(worldULF, worldURF, plane);
    addAnchor(worldULF, worldULB, plane);
    addAnchor(worldLRF, worldURF, plane);
    addAnchor(worldLRF, worldLRB, plane);

    // Less than 3 anchors (indicates that the clipping plane only touches the AABB)?
    if (anchors_.size() < 3)
        anchors_.clear();

    tgt::mat4 t;
    physicalToWorld_.getRotationalPart().invert(t);

    //determine the normal vectors for the planes
    //tgt::vec3 nx = tgt::normalize(worldLRF - worldLLF);
    //tgt::vec3 ny = tgt::normalize(worldULF - worldLLF);
    //tgt::vec3 nz = tgt::normalize(worldLLB - worldLLF);
    tgt::vec3 nx = tgt::normalize((tgt::transpose(t) * tgt::vec4(1.f, 0.f, 0.f, 1.f)).xyz());
    tgt::vec3 ny = tgt::normalize((tgt::transpose(t) * tgt::vec4(0.f, 1.f, 0.f, 1.f)).xyz());
    tgt::vec3 nz = tgt::normalize((tgt::transpose(t) * tgt::vec4(0.f, 0.f, 1.f, 1.f)).xyz());
    addLine(tgt::vec4( nx ,  tgt::dot(nx, worldURB)));
    addLine(tgt::vec4( nx ,  tgt::dot(nx, worldLLF)));
    addLine(tgt::vec4( ny ,  tgt::dot(ny, worldURB)));
    addLine(tgt::vec4( ny ,  tgt::dot(ny, worldLLF)));
    addLine(tgt::vec4( nz ,  tgt::dot(nz, worldURB)));
    addLine(tgt::vec4( nz ,  tgt::dot(nz, worldLLF)));


    // Update actual plane equation to ensure that this update function is only called when necessary...
    actualPlane_ = tgt::vec4(planeNormal_.get(), planePosition_.get());

    // No lines found?
    if (lines_.size() == 0)
        return;

    //
    // Order lines and their anchors to CCW polygon vertex order which is expected by OpenGL...
    //

    // Sort lines to produce contiguous vertex order...
    bool reverseLastLine = false;
    for (std::vector<Line>::size_type i = 0; i < lines_.size() - 1; ++i) {
        for (std::vector<Line>::size_type j = i + 1; j < lines_.size(); ++j) {
            Anchor* connectionAnchor;
            if (reverseLastLine)
                connectionAnchor = lines_.at(i).first;
            else
                connectionAnchor = lines_.at(i).second;

            if (std::abs(tgt::length(*connectionAnchor - *(lines_.at(j).first))) < EPSILON) {
                std::swap(lines_.at(i + 1), lines_.at(j));
                reverseLastLine = false;
                break;
            }
            else if (std::abs(tgt::length(*connectionAnchor - *(lines_.at(j).second))) < EPSILON) {
                std::swap(lines_.at(i + 1), lines_.at(j));
                reverseLastLine = true;
                break;
            }
        }
    }

    std::vector<Anchor*> vertices;
    // Convert sorted line list to sorted vertex list...
    for (std::vector<Line>::size_type i = 0; i < lines_.size(); ++i) {
        bool reverseLine = i != 0 && std::abs(tgt::length(*vertices.at(vertices.size() - 1) - *(lines_.at(i).first))) >= EPSILON;

        Anchor* first = (reverseLine ? lines_.at(i).second : lines_.at(i).first);
        Anchor* second = (reverseLine ? lines_.at(i).first : lines_.at(i).second);

        if (i == 0)
            vertices.push_back(first);

        if (i < (vertices.size() - 1))
            vertices.push_back(second);
    }

    // Convert vertex order to counter clockwise if necessary...
    tgt::vec3 normal(0, 0, 0);
    for (std::vector<Anchor*>::size_type i = 0; i < vertices.size(); ++i)
        normal += tgt::cross(*vertices[i], *vertices[(i + 1) % vertices.size()]);
    normal = tgt::normalize(normal);

    if (tgt::dot(plane.xyz(), normal) < 0) {
        std::reverse(lines_.begin(), lines_.end());
        for (std::vector<Line>::iterator it = lines_.begin(); it != lines_.end(); ++it)
            std::swap(it->first, it->second);
    }
}

void PlaneManipulation::addAnchor(const tgt::vec3& startVertex, const tgt::vec3& endVertex, const tgt::vec4& plane) {
    tgt::vec3 intersectionPoint;
    // No intersection point withing the bounding box?
    if (!getIntersectionPointInBB(intersectionPoint, startVertex, endVertex, plane))
        return;

    // Omit anchor, if it already exist...
    for (std::vector<Anchor>::iterator it = anchors_.begin(); it != anchors_.end(); ++it)
        if (std::abs(tgt::length(*it - intersectionPoint)) < EPSILON)
            return;

    anchors_.push_back(intersectionPoint);
}

void PlaneManipulation::addLine(const tgt::vec4& plane) {

    Anchor* begin = 0;
    Anchor* end = 0;

    for (std::vector<Anchor>::iterator it = anchors_.begin(); it != anchors_.end(); ++it) {
        // Is anchor on given plane? TODO: epsilon should not be multiplied, find a better solution for this to account for numerical problems
        if (std::abs(tgt::dot(plane.xyz(), *it) - plane.w) < EPSILON * 100.f) {
            if (!begin)
                begin = &(*it);
            else {
                end = &(*it);
                break;
            }
        }
    }
    // Was line start and end point found?
    if (begin && end)
        lines_.push_back(std::make_pair(begin, end));
}

//------------------------------------------------------------------------//
//  Manipulator Functions                                                 //
//------------------------------------------------------------------------//
void PlaneManipulation::paintManipulatorCylinder(const tgt::vec3& center, const tgt::vec3& direction, tgt::Shader* shader) {
    paintManipulatorCylinderTopHalf(center, direction, shader);
    paintManipulatorCylinderBottomHalf(center, direction, shader);
}

void PlaneManipulation::paintManipulatorCylinderTopHalf(const tgt::vec3& center, const tgt::vec3& direction, tgt::Shader* shader) {
    MatStack.pushMatrix();

    // 3. translate to new center
    MatStack.translate(center.x, center.y, center.z);

    // 2. rotate to new direction
    float angle = static_cast<float>(std::acos(tgt::dot(tgt::normalize(direction), tgt::vec3(0.f, 0.f, 1.f))));
    tgt::vec3 rotationAxis = tgt::cross(tgt::normalize(direction), tgt::vec3(0.f, 0.f, 1.f));

    if (tgt::equal(rotationAxis, tgt::vec3(0.f)) != tgt::bvec3(true, true, true)) {
        tgt::mat4 matrix = tgt::mat4::createRotation(-angle, rotationAxis);
        MatStack.multMatrix(matrix);
    }
    else {
        //if there can no rotation axis be found the direction lies in positive or negative z direction
        //for negative direction: rotate around x axis, else: do nothing
        if (direction.z < 0.f)
            MatStack.rotate(180.f, 1.f, 0.f, 0.f);
    }

    //update camera settings
    setMatrixShaderUniforms(shader);
    shader->setUniform("translate_", tgt::vec3(0.f));

    //create geometry lazy
    if (!topCylinderGeometry_) {
        topCylinderGeometry_ = new GlMeshGeometryUInt16Normal();
        topCylinderGeometry_->setCylinderGeometry(generalColor_.get(),manipulatorDiameter_, manipulatorDiameter_ / 2.f,
                                                 manipulatorLength_ / 2.f,20,20,false,false); //color not used (shader -> material)
    }

    //render
    topCylinderGeometry_->render();

    MatStack.popMatrix();
}

void PlaneManipulation::paintManipulatorCylinderBottomHalf(const tgt::vec3& center, const tgt::vec3& direction, tgt::Shader* shader) {
    MatStack.pushMatrix();

    // 3. translate to new center
    MatStack.translate(center.x, center.y, center.z);

    // 2. rotate to new direction
    float angle = static_cast<float>(std::acos(tgt::dot(tgt::normalize(direction), tgt::vec3(0.f, 0.f, 1.f))));
    tgt::vec3 rotationAxis = tgt::cross(tgt::normalize(direction), tgt::vec3(0.f, 0.f, 1.f));

    if (tgt::equal(rotationAxis, tgt::vec3(0.f)) != tgt::bvec3(true, true, true)) {
        tgt::mat4 matrix = tgt::mat4::createRotation(-angle, rotationAxis);
        MatStack.multMatrix(matrix);
    }
    else {
        //if there can no rotation axis be found the direction lies in positive or negative z direction
        //for negative direction: rotate around x axis, else: do nothing
        if (direction.z < 0.f)
            MatStack.rotate(180.f, 1.f, 0.f, 0.f);
    }

    // 1. cylinder: translate center to origin
    MatStack.translate(0.f, 0.f, -manipulatorLength_ / 2.f);

    //update camera settings
    setMatrixShaderUniforms(shader);
    shader->setUniform("translate_", tgt::vec3(0.f));

    //create geometry lazy
    if (!bottomCylinderGeometry_) {
        bottomCylinderGeometry_ = new GlMeshGeometryUInt16Normal();
        bottomCylinderGeometry_->setCylinderGeometry(generalColor_.get(),manipulatorDiameter_ / 2.f, manipulatorDiameter_,
                                                     manipulatorLength_ / 2.f,20,20,false,false); //color not used (shader -> material)
    }

    //render
    bottomCylinderGeometry_->render();

    MatStack.popMatrix();
}

void PlaneManipulation::paintManipulatorTip(const tgt::vec3& center, tgt::Shader* shader) {
    if(!tipGeometry_) {
        tipGeometry_ = new GlMeshGeometryUInt16Normal();
        tipGeometry_->setSphereGeometry(manipulatorTipRadius_,tgt::vec3::zero,generalColor_.get(),32); //color not used (shader->color)
    }
    shader->setUniform("translate_", center);
    tipGeometry_->render();
}

//------------------------------------------------------------------------//
//  Mouse Events and Callbacks                                            //
//------------------------------------------------------------------------//
void PlaneManipulation::planeManipulation(tgt::MouseEvent* e) {

    e->ignore();

    if (!idManager_)
        return;


    //no plane manipulation if rendering of manipulator is not enabled
    if (!enable_.get() || !renderManipulator_.get())
        return;

    //no plane manipulation / error message if no inport is connected
    if (!(volumeInport_.getData() || geometryInport_.getData())) {
        //LWARNING("volume or geometry inport has to be connected, ignoring events");
        return;
    }

    if (e->button() == tgt::MouseEvent::MOUSE_BUTTON_LEFT) {
        if (e->action() & tgt::MouseEvent::PRESSED) {
            const void* obj = idManager_->getObjectAtPos(tgt::ivec2(e->coord().x, e->viewport().y - e->coord().y));

            //check selection
            if (&topManipulatorPos_ == obj) {
                manipulatorTopTipIsGrabbed_ = true;
                QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
                e->accept();

                //check if manipulating this tip is allowed by comparing direction of view vector to manipulator orientation via dot product
                float check = tgt::dot(manipulatorOrientation_, camera_.getPositionWithOffsets() - manipulatorCenter_.get());

                if (check < -EPSILON*1000.f)
                    invalidMovement_ = true;
                else
                    invalidMovement_ = false;

                invalidate();
            }
            else if (&bottomManipulatorPos_ == obj) {
                manipulatorBottomTipIsGrabbed_ = true;
                QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
                e->accept();

                //check if manipulating this tip is allowed by comparing direction of view vector to manipulator orientation via dot product
                float check = tgt::dot(manipulatorOrientation_, camera_.getPositionWithOffsets() - manipulatorCenter_.get());

                if (check > EPSILON*1000.f)
                    invalidMovement_ = true;
                else
                    invalidMovement_ = false;

                invalidate();
            }
            else if (&manipulatorCenter_ == obj) {
                manipulatorCylinderIsGrabbed_ = true;
                QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
                e->accept();

                //compute scale factor to avoid numerical problems with small spacing
                scale_ = manipulatorLength_ / 200.f;

                //compute offset
                tgt::vec2 mousePos = tgt::vec2(static_cast<float>(e->coord().x), e->viewport().y-static_cast<float>(e->coord().y));

                tgt::vec3 begin = getWindowPos(manipulatorCenter_.get());
                tgt::vec3 end = getWindowPos(manipulatorCenter_.get() + manipulatorOrientation_ * scale_);

                tgt::vec3 point = begin;
                tgt::vec3 direction = end - begin;

                shiftOffset_ = getPointLineProjectionScalar(mousePos, point.xy(), direction.xy());

                invalidate();
            }
        }

        if (e->action() & tgt::MouseEvent::MOTION) {
            if (manipulatorTopTipIsGrabbed_ || manipulatorBottomTipIsGrabbed_ || manipulatorCylinderIsGrabbed_)
                e->accept();

            if (manipulatorCylinderIsGrabbed_) {
                // Determine mouse move offset...
                tgt::vec2 mousePos = tgt::vec2(static_cast<float>(e->coord().x), e->viewport().y-static_cast<float>(e->coord().y));

                tgt::vec3 begin = getWindowPos(manipulatorCenter_.get());
                tgt::vec3 end = getWindowPos(manipulatorCenter_.get() + manipulatorOrientation_ * scale_);

                tgt::vec3 point = begin;
                tgt::vec3 direction = end - begin;

                float t = getPointLineProjectionScalar(mousePos, point.xy(), direction.xy()) - shiftOffset_;

                //check if new manipulator position is valid (depending on the plane position)
                tgt::vec3 tmpCenter =  manipulatorCenter_.get() + t * manipulatorOrientation_ * scale_;
                float tmpPosition = tgt::dot(tmpCenter, manipulatorOrientation_);

                if ((positionRange_.get().x <= tmpPosition) && (tmpPosition <= positionRange_.get().y)) {
                    planePosition_.set(tmpPosition);
                    manipulatorCenter_.set(tmpCenter);
                    invalidMovement_ = false;
                }
                else {
                    invalidMovement_ = true;
                }

                invalidate();
            }
            else if (manipulatorTopTipIsGrabbed_ || manipulatorBottomTipIsGrabbed_) {

                //check if manipulating this tip is allowed by comparing direction of view vector to manipulator orientation via dot product
                float check = tgt::dot(manipulatorOrientation_, camera_.getPositionWithOffsets() - manipulatorCenter_.get());

                if ((manipulatorTopTipIsGrabbed_ && !(check >= -EPSILON*1000.f)) || (manipulatorBottomTipIsGrabbed_ && !(check <= EPSILON*1000.f))) {
                    invalidMovement_ = true;
                    invalidate();
                    return;
                }
                else
                    invalidMovement_ = false;

                //get inverse view and projection matrices
                tgt::mat4 projectionMatrixInverse;
                camera_.getProjectionMatrix(e->viewport()).invert(projectionMatrixInverse);
                tgt::mat4 viewMatrixInverse = camera_.getViewMatrixInverse();

                //normalize mouse position
                tgt::vec2 pos = tgt::vec2(static_cast<float>(e->coord().x), e->viewport().y-static_cast<float>(e->coord().y)) / tgt::vec2(e->viewport());
                pos *= tgt::vec2(2.f);
                pos -= tgt::vec2(1.f);
                //get line of mouse position in world coordinates
                tgt::vec4 startWorld = viewMatrixInverse *  (projectionMatrixInverse *  tgt::vec4(pos, -1.f, 1.f));
                tgt::vec4 endWorld = viewMatrixInverse *  (projectionMatrixInverse *  tgt::vec4(pos, 1.f, 1.f));

                startWorld *= 1.f/startWorld.w;
                endWorld *= 1.f/endWorld.w;

                tgt::vec3 lineNormal = tgt::normalize(endWorld.xyz() - startWorld.xyz());

                //translate line so that sphere around manipulator is at origin
                tgt::vec3 startWorldTranslated = startWorld.xyz() - manipulatorCenter_.get();
                //tgt::vec3 endWorldTranslated = endWorld.xyz() - manipulatorCenter_;

                float radius = manipulatorLength_ / 2.f;

                //check if line intersects sphere
                float dis = static_cast<float>(std::pow(tgt::dot(lineNormal, startWorldTranslated), 2.f))
                        - tgt::dot(startWorldTranslated, startWorldTranslated) + static_cast<float>(std::pow(radius, 2.f));

                tgt::vec3 tmpOrientation; // computed orientation

                if (dis >= 0.f) {
                    //hits sphere -> compute position on sphere
                    float t = -tgt::dot(lineNormal, startWorldTranslated) - std::sqrt(dis);
                    tgt::vec3 pos = (startWorld.xyz()) + t * lineNormal;
                    //compute orientation using the position
                    if (manipulatorTopTipIsGrabbed_)
                        tmpOrientation = tgt::normalize(pos - manipulatorCenter_.get());
                    else
                        tmpOrientation = tgt::normalize(manipulatorCenter_.get() - pos);
                }
                else {
                    //sphere has not been hit -> compute intersection with plane aligned with viewing plane through sphere center
                    tgt::vec3 planeNormal = tgt::normalize(camera_.getPositionWithOffsets() - manipulatorCenter_.get());
                    float t = tgt::dot(manipulatorCenter_.get() - startWorld.xyz(), planeNormal) / tgt::dot(lineNormal, planeNormal);
                    tgt::vec3 intersectionPoint = startWorld.xyz() + lineNormal * t;
                    //project to the intersection circle of sphere and plane
                    tgt::vec3 projectionPoint = tgt::normalize(intersectionPoint - manipulatorCenter_.get()) * radius + manipulatorCenter_.get();
                    //compute orientation
                    if (manipulatorTopTipIsGrabbed_)
                        tmpOrientation = tgt::normalize(projectionPoint - manipulatorCenter_.get());
                    else
                        tmpOrientation = tgt::normalize(manipulatorCenter_.get() - projectionPoint);
                }

                //check if the orientation is valid for the current handle position (to ensure that the plane position is within the valid range)
                float tmpPosition = tgt::dot(manipulatorCenter_.get(), tmpOrientation);

                if ((positionRange_.get().x <= tmpPosition) && (tmpPosition <= positionRange_.get().y)) {
                    manipulatorOrientation_ = tmpOrientation;
                    updateManipulatorElements();
                    planeNormal_.set(manipulatorOrientation_);
                    planePosition_.set(tmpPosition);
                }
                else {
                    invalidMovement_ = true;
                }
                invalidate();
            }
        }

        if (e->action() & tgt::MouseEvent::RELEASED) {
            if(manipulatorTopTipIsGrabbed_ || manipulatorBottomTipIsGrabbed_ || manipulatorCylinderIsGrabbed_) {
                // Accept mouse event, turn off interaction mode, ungrab anchor and lines, and repaint...
                e->accept();
                QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);
                manipulatorTopTipIsGrabbed_ = false;
                manipulatorBottomTipIsGrabbed_ = false;
                manipulatorCylinderIsGrabbed_ = false;

                if (autoReset_.get() && !anchors_.empty())
                    setManipulatorToStdPos();

                invalidMovement_ = false;

                invalidate();
            }
        }
    }

    if (e->button() == tgt::MouseEvent::MOUSE_BUTTON_RIGHT) {

        if (e->action() & tgt::MouseEvent::PRESSED) {
            const void* obj = idManager_->getObjectAtPos(tgt::ivec2(e->coord().x, e->viewport().y - e->coord().y));

            //check selection
            if ((obj == &topManipulatorPos_) || (obj == &bottomManipulatorPos_) || (obj == &manipulatorCenter_)) {
                wholeManipulatorIsGrabbed_ = true;
                e->accept();

                //compute scale factor to avoid numerical problems with small spacing
                scale_ = manipulatorLength_ / 200.f;

                //compute offset
                tgt::vec2 mousePos = tgt::vec2(static_cast<float>(e->coord().x), e->viewport().y-static_cast<float>(e->coord().y));

                tgt::vec3 begin = getWindowPos(manipulatorCenter_.get());
                tgt::vec3 end = getWindowPos(manipulatorCenter_.get() + manipulatorOrientation_ * scale_);

                tgt::vec3 point = begin;
                tgt::vec3 direction = end - begin;

                shiftOffset_ = getPointLineProjectionScalar(mousePos, point.xy(), direction.xy());

                //check if the manipulator may be moved on the plane
                if (std::abs(tgt::dot(manipulatorOrientation_, tgt::normalize(camera_.getLook()))) < 0.2f)
                   invalidMovement_ = true;

                invalidate();
            }
        }

        if (e->action() & tgt::MouseEvent::MOTION) {

            if (wholeManipulatorIsGrabbed_) {

                //check if the manipulator may be moved on the plane
                float check = tgt::dot(manipulatorOrientation_, tgt::normalize(camera_.getLook()));

                if (std::abs(check) < 0.2f) {
                    invalidMovement_ = true;
                    invalidate();
                    return;
                }
                else
                    invalidMovement_ = false;

                e->accept();

                if (anchors_.empty())
                    return;

                //get line through mouse position:
                //1. get matrices
                tgt::mat4 projectionMatrixInverse;
                camera_.getProjectionMatrix(e->viewport()).invert(projectionMatrixInverse);
                tgt::mat4 viewMatrixInverse = camera_.getViewMatrixInverse();

                //2. normalize mouse position
                tgt::vec2 pos = tgt::vec2(static_cast<float>(e->coord().x), e->viewport().y-static_cast<float>(e->coord().y)) / tgt::vec2(e->viewport());
                pos *= tgt::vec2(2.f);
                pos -= tgt::vec2(1.f);
                //3. get line of mouse position in world coordinates
                tgt::vec4 startWorld = viewMatrixInverse *  (projectionMatrixInverse *  tgt::vec4(pos, -1.f, 1.f));
                tgt::vec4 endWorld = viewMatrixInverse *  (projectionMatrixInverse *  tgt::vec4(pos, 1.f, 1.f));

                startWorld *= 1.f/startWorld.w;
                endWorld *= 1.f/endWorld.w;

                tgt::vec3 lineNormal = tgt::normalize(endWorld.xyz() - startWorld.xyz());

                //get the plane and compute the intersection point
                tgt::vec3 planeNormal = tgt::normalize(planeNormal_.get());
                float t = tgt::dot(manipulatorCenter_.get() + manipulatorOrientation_ * scale_ * shiftOffset_ - startWorld.xyz(), planeNormal)
                                        / tgt::dot(lineNormal, planeNormal);

                tgt::vec3 intersectionPoint = startWorld.xyz() + lineNormal * t;

                //check if new position is valid
                tgt::vec3 tmpCenter = intersectionPoint - shiftOffset_ * scale_ * manipulatorOrientation_;

                if ((tgt::lessThanEqual(tmpCenter, manipulatorBounds_.getURB()) == tgt::bvec3(true))
                        && tgt::greaterThanEqual(tmpCenter, manipulatorBounds_.getLLF()) == tgt::bvec3(true))
                {
                    manipulatorCenter_.set(tmpCenter);
                }
                else
                    invalidMovement_ = true;

                invalidate();
            }
        }

        if (e->action() & tgt::MouseEvent::RELEASED) {

            if (wholeManipulatorIsGrabbed_) {
                // Accept mouse event, turn off interaction mode, ungrab anchor and lines, and repaint...
                e->accept();
                wholeManipulatorIsGrabbed_ = false;

                invalidMovement_ = false;

                invalidate();
            }
        }
    }
}

void PlaneManipulation::checkInportsForBoundingBox() {

    if (!volumeInport_.getData() && !geometryInport_.getData())
        return;

    //get bounding box, first try to get from volume, else get from geometry
    tgt::Bounds bBoxWorld;
    tgt::Bounds bBoxPhysical;
    if (volumeInport_.getData()) {
        bBoxWorld = volumeInport_.getData()->getBoundingBox().getBoundingBox();
        bBoxPhysical = volumeInport_.getData()->getBoundingBox(false).getBoundingBox(false);
    }
    else {
        bBoxWorld = geometryInport_.getData()->getBoundingBox();
        bBoxPhysical = geometryInport_.getData()->getBoundingBox(false);
    }

    //check if scene bounds have changed -> update range for plane position and size of manipulator according to current scene
    if ((!worldSceneBounds_.isDefined() ||
                (bBoxWorld.getLLF() != worldSceneBounds_.getLLF())
                || (bBoxWorld.getURB() != worldSceneBounds_.getURB()))
        || (!physicalSceneBounds_.isDefined() ||
                (bBoxPhysical.getLLF() != physicalSceneBounds_.getLLF())
                || (bBoxPhysical.getURB() != physicalSceneBounds_.getURB())))
    {

        //set new bounds
        worldSceneBounds_ = bBoxWorld;
        manipulatorBounds_ = tgt::Bounds(bBoxWorld.center() + (bBoxWorld.getLLF() - bBoxWorld.center()) * 6.f,
                bBoxWorld.center() + (bBoxWorld.getURB() - bBoxWorld.center()) * 6.f);

        physicalSceneBounds_ = bBoxPhysical;

        //get transformations
        if (volumeInport_.getData()) {
            physicalToWorld_ = volumeInport_.getData()->getPhysicalToWorldMatrix();
        }
        else {
            physicalToWorld_ = geometryInport_.getData()->getTransformationMatrix();
        }

        //check every corner of world scene bounds and compute position range
        tgt::vec3 llf = bBoxWorld.getLLF();
        tgt::vec3 urb = bBoxWorld.getURB();
        float length = std::max(tgt::length(llf), tgt::length(urb));

        tgt::vec3 llb = tgt::vec3(llf.x, llf.y, urb.z);
        tgt::vec3 urf = tgt::vec3(urb.x, urb.y, llf.z);
        length = std::max(length, std::max(tgt::length(llb), tgt::length(urf)));

        tgt::vec3 lrf = tgt::vec3(llf.x, urb.y, llf.z);
        tgt::vec3 ulb = tgt::vec3(urb.x, llf.y, urb.z);
        length = std::max(length, std::max(tgt::length(lrf), tgt::length(ulb)));

        tgt::vec3 lrb = tgt::vec3(llf.x, urb.y, urb.z);
        tgt::vec3 ulf = tgt::vec3(urb.x, llf.y, llf.z);
        length = std::max(length, std::max(tgt::length(lrb), tgt::length(ulf)));

        if (length == 0.f)
            length = 1.f;

        positionRange_.set(tgt::vec2(-length, length));

        //check if plane position is valid and adjust if necessary
        if (planePosition_.get() < -length)
            planePosition_.set(-length);
        else if (planePosition_.get() > length)
            planePosition_.set(length);

        //if the bounding box changed: reset the manipulator
        setManipulatorToStdPos();
        updateManipulatorScale();
        updateAnchorsAndLines();
    }
}

void PlaneManipulation::updateManipulatorScale(){
    //compute length and diameter of manipulator handle according to scene size
    float diagonal = 1.f;
    if (volumeInport_.getData()) {
        diagonal = tgt::length(physicalSceneBounds_.diagonal());
    } else {
        diagonal = tgt::length(worldSceneBounds_.diagonal());
    }

    manipulatorLength_ = 0.25f * diagonal * manipulatorScale_.get();
    manipulatorDiameter_ = 0.05f * manipulatorLength_;
    manipulatorTipRadius_ = 0.05f * manipulatorLength_;
    updateManipulatorElements();
    //clear old geometry
    delete anchorGeometry_; anchorGeometry_ = 0;
    delete tipGeometry_; tipGeometry_ = 0;
    delete topCylinderGeometry_; topCylinderGeometry_ = 0;
    delete bottomCylinderGeometry_; bottomCylinderGeometry_ = 0;
}

void PlaneManipulation::setManipulatorToStdPos() {

    //do not change position if the manipulator is currently in use
    if (manipulatorTopTipIsGrabbed_ || manipulatorBottomTipIsGrabbed_ || manipulatorCylinderIsGrabbed_)
        return;

    tgt::vec3 normal = tgt::normalize(planeNormal_.get());
    tgt::vec3 center = worldSceneBounds_.center();
    tgt::vec3 pointOnPlane = normal * planePosition_.get();
    float c = tgt::dot(normal, center - pointOnPlane);

    manipulatorOrientation_ = normal;
    manipulatorCenter_.set(center - normal * c);

    invalidate();
}

void PlaneManipulation::updateManipulatorElements(bool updatePosition) {
    topManipulatorPos_ = manipulatorCenter_.get() + tgt::normalize(manipulatorOrientation_) * (manipulatorLength_ / 2.f);
    bottomManipulatorPos_ = manipulatorCenter_.get()  - tgt::normalize(manipulatorOrientation_) * (manipulatorLength_ / 2.f);

    if(updatePosition && isInitialized()) { //hack to prevent deserialization warnings
        //callbacks are blocked, since the planePosition resets the center otherwise
        planePosition_.blockCallbacks(true);
        planePosition_.set(tgt::dot(manipulatorCenter_.get(), manipulatorOrientation_));
        planePosition_.blockCallbacks(false);
    }
}

void PlaneManipulation::resetPlane() {
    //set normal to -z in physical coordinates
    planeNormal_.set(tgt::normalize((physicalToWorld_.getRotationalPart() * tgt::vec4(0.f,0.f,-1.f, 1.f)).xyz()));

    //set position to center of bounding box in world coordinates
    float position = tgt::dot(worldSceneBounds_.center(), planeNormal_.get());
    planePosition_.set(position);
}

void PlaneManipulation::updatePropertyVisibility() {
    //plane settings
    lineWidth_.setVisibleFlag(renderPlane_.get());
    grabbedPlaneOpacity_.setVisibleFlag(renderPlane_.get());
    alwaysUseOpacity_.setVisibleFlag(renderPlane_.get());
    renderAnchors_.setVisibleFlag(renderPlane_.get());
    //manipulator settings
    manipulatorScale_.setVisibleFlag(renderManipulator_.get());
    grabbedElementColor_.setVisibleFlag(renderManipulator_.get());
    blockedElementColor_.setVisibleFlag(renderManipulator_.get());
    resetManipulatorPos_.setVisibleFlag(renderManipulator_.get());
    autoReset_.setVisibleFlag(renderManipulator_.get());
}

bool PlaneManipulation::shouldRenderPlane() const {
    return (alwaysUseOpacity_.get() && grabbedPlaneOpacity_.get() > 0) ||
        (manipulatorTopTipIsGrabbed_ || manipulatorBottomTipIsGrabbed_ || manipulatorCylinderIsGrabbed_ || wholeManipulatorIsGrabbed_)
        && !invalidMovement_;
}

} //namespace voreen
