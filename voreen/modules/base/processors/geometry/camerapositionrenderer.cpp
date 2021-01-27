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

#include "camerapositionrenderer.h"

namespace voreen {

CameraPositionRenderer::CameraPositionRenderer()
    : GeometryRendererBase()
    , sizeReferencePort_(Port::INPORT,"sizeReferencePort","Reference Port")
    , displayCameraProp_("displayCameraProp", "Camera to be linked", tgt::Camera(tgt::vec3(0.f, 0.f, 3.5f), tgt::vec3(0.f, 0.f, 0.f), tgt::vec3(0.f, 1.f, 0.f)))
    , myShaderProp_("myShaderProp", "Shader", "camposition.frag", "camposition.vert")
    , renderFrustumProp_("enableFrustumProp", "Render frustum", true)
    , frustLineColorProp_("frustLineColorProp", "Color of frustum lines", tgt::vec4(0.f, 0.5f, 1.f, 1.f), Processor::INVALID_RESULT, Property::LOD_DEFAULT, false)
    , frustLineWidthProp_("frustLineWidthProp", "Width of frustum lines", 1.0f, 1.0f, 5.0f)
    , renderPlanesProp_("enablePlanesProp", "Render near and far plane", true)
    , opacPlanesProp_("opacPlanesProp", "Plane opacity")
    , optionModelProp_("optionModelProp", "Display camera model", true)
    , opacModelProp_("opacModelProp", "Camera model opacity")
    , colorModelProp_("colorModelProp", "Color of camera model", tgt::vec4(0.0f, 0.0f, 1.0f, 1.f), Processor::INVALID_RESULT, Property::LOD_DEFAULT, false)
    , lineModelProp_("lineModelProp", "Width of model lines", 1.0f, 1.0f, 5.0f)
    , scalingModelProp_("scalingModelProp", "Scale of camera model", 1.0f, 0.1f, 50.0f)
    , optionArrowsProp_("optionArrowsProp", "Display orientation of linked camera", true)
    , scalingArrowsProp_("scalingArrowsProp",  "Scale of direction arrows", 1.0f, 0.1f, 50.0f)
    , xLengthProp_("xLengthProp", "Length of x axis", 1.f, -50.f, 50.f)
    , yLengthProp_("yLengthProp", "Length of y axis", 1.f, -50.f, 50.f)
    , zLengthProp_("zLengthProp", "Length of z axis", 1.f, -50.f, 50.f)
{
    //this processor should block events by default
    outPort_.setBlockEvents(true);

    addPort(sizeReferencePort_);

    // build the GUI
    addProperty(displayCameraProp_);
    addProperty(myShaderProp_);

    // frustum-related settings
    addProperty(renderFrustumProp_);
    renderFrustumProp_.setGroupID("frustum");
    renderFrustumProp_.onChange(MemberFunctionCallback<CameraPositionRenderer>(this, &CameraPositionRenderer::onRenderFrustumChanged));
    onRenderFrustumChanged();
    addProperty(frustLineColorProp_);
    frustLineColorProp_.setGroupID("frustum");
    addProperty(frustLineWidthProp_);
    frustLineWidthProp_.setGroupID("frustum");
    addProperty(renderPlanesProp_);
    renderPlanesProp_.setGroupID("frustum");
    renderPlanesProp_.onChange(MemberFunctionCallback<CameraPositionRenderer>(this, &CameraPositionRenderer::onRenderPlanesChanged));
    onRenderPlanesChanged();
    addProperty(opacPlanesProp_);
    opacPlanesProp_.setGroupID("frustum");
    setPropertyGroupGuiName("frustum", "Frustum and Planes Settings");

    // camera model-related settings
    optionModelProp_.addOption("off", "Off");
    optionModelProp_.addOption("dot", "Dot");
    optionModelProp_.addOption("lines", "Lines");
    optionModelProp_.addOption("model", "Model");
    addProperty(optionModelProp_);
    optionModelProp_.setGroupID("model");
    optionModelProp_.onChange(MemberFunctionCallback<CameraPositionRenderer>(this, &CameraPositionRenderer::onRenderModelChanged));
    onRenderModelChanged();
    addProperty(opacModelProp_);
    opacModelProp_.setGroupID("model");
    addProperty(colorModelProp_);
    colorModelProp_.setGroupID("model");
    addProperty(lineModelProp_);
    lineModelProp_.setGroupID("model");
    addProperty(scalingModelProp_);
    scalingModelProp_.setGroupID("model");
    setPropertyGroupGuiName("model", "Camera Model Settings");

    // arrows/axis-related settings
    optionArrowsProp_.addOption("off", "Off");
    optionArrowsProp_.addOption("arrowsxyz", "Direction Arrows (xyz)");
    optionArrowsProp_.addOption("arrowsxy", "Direction Arrows (xy)");
    optionArrowsProp_.addOption("coords", "Coordinate System");
    optionArrowsProp_.select("arrowsxyz");
    addProperty(optionArrowsProp_);
    optionArrowsProp_.setGroupID("arrows");
    optionArrowsProp_.onChange(MemberFunctionCallback<CameraPositionRenderer>(this, &CameraPositionRenderer::onRenderArrowsChanged));
    onRenderArrowsChanged();
    addProperty(scalingArrowsProp_);
    scalingArrowsProp_.setGroupID("arrows");
    addProperty(xLengthProp_);
    xLengthProp_.setGroupID("arrows");
    addProperty(yLengthProp_);
    yLengthProp_.setGroupID("arrows");
    addProperty(zLengthProp_);
    zLengthProp_.setGroupID("arrows");
    setPropertyGroupGuiName("arrows", "Arrows and Coordinate System Settings");
}

CameraPositionRenderer::~CameraPositionRenderer() {
}

Processor* CameraPositionRenderer::create() const {
    return new CameraPositionRenderer();
}

std::string CameraPositionRenderer::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string header = GeometryRendererBase::generateHeader(version);
    return header;
}

void CameraPositionRenderer::initialize() {

    GeometryRendererBase::initialize();

    myShaderProp_.setHeader(generateHeader());
    myShaderProp_.rebuild();

    // generate OpenGL buffers and arrays
    glGenBuffers(3,vbo_);
    glGenVertexArrays(3,vao_);

    // build vertex arrays
    createFrustumBuffer();
    createCameraModelBuffer();
    createArrowsBuffer();
}


void CameraPositionRenderer::deinitialize() {
    // cleanup
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glDeleteVertexArrays(3,vao_);
    glDeleteBuffers(3,vbo_);

    GeometryRendererBase::deinitialize();
}

void CameraPositionRenderer::createFrustumBuffer(){
    // frustum corners
    tgt::vec4 llf = tgt::vec4(-1.f, -1.f, -1.f, 1.f);
    tgt::vec4 ulf = tgt::vec4(-1.f, 1.f, -1.f, 1.f);
    tgt::vec4 llb = tgt::vec4(-1.f, -1.f, 1.f, 1.f);
    tgt::vec4 ulb = tgt::vec4(-1.f, 1.f, 1.f, 1.f);
    tgt::vec4 lrf = tgt::vec4(1.f, -1.f, -1.f, 1.f);
    tgt::vec4 urf = tgt::vec4(1.f, 1.f, -1.f, 1.f);
    tgt::vec4 lrb = tgt::vec4(1.f, -1.f, 1.f, 1.f);
    tgt::vec4 urb = tgt::vec4(1.f, 1.f, 1.f, 1.f);

    // define array
    tgt::vec4 frustumcorners[24] = { llf, ulf, lrf, urf, llb, ulb, lrb, urb, llf, lrf, ulf, urf, llb, lrb, ulb, urb, llf, llb, ulf, ulb, lrf, lrb, urf, urb };

    // setup VBO
    int SIZE_FPOS = sizeof(frustumcorners);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_[0]);
    glBufferData(GL_ARRAY_BUFFER, SIZE_FPOS, (GLvoid*)frustumcorners, GL_STATIC_DRAW);

    // setup VAO
    glBindVertexArray(vao_[0]);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);

    //clean-up
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void CameraPositionRenderer::createCameraModelBuffer(){
    // camera model, rendered as truncated pyramid without top and bottom
    tgt::vec4 zero = tgt::vec4(0.f, 0.f, 0.f, 1.f);
    tgt::vec4 llf = tgt::vec4(-0.1f, -0.1f, 0.5f, 1.f);
    tgt::vec4 ulf = tgt::vec4(-0.1f, 0.1f, 0.5f, 1.f);
    tgt::vec4 llb = tgt::vec4(-0.3f, -0.3f, 0.f, 1.f);
    tgt::vec4 ulb = tgt::vec4(-0.3f, 0.3f, 0.f, 1.f);
    tgt::vec4 lrf = tgt::vec4(0.1f, -0.1f, 0.5f, 1.f);
    tgt::vec4 urf = tgt::vec4(0.1f, 0.1f, 0.5f, 1.f);
    tgt::vec4 lrb = tgt::vec4(0.3f, -0.3f, 0.f, 1.f);
    tgt::vec4 urb = tgt::vec4(0.3f, 0.3f, 0.f, 1.f);

    // define array
    tgt::vec4 cameramodel[11] = { llf, llb, ulf, ulb, urf, urb, lrf, lrb, llf, llb, zero };

    // setup VBO
    int SIZE_FPOS = sizeof(cameramodel);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_[1]);
    glBufferData(GL_ARRAY_BUFFER, SIZE_FPOS, (GLvoid*)cameramodel, GL_STATIC_DRAW);

    // setup VAO
    glBindVertexArray(vao_[1]);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);

    //clean-up
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void CameraPositionRenderer::createArrowsBuffer(){
    //lines and arrowheads in 3 directions
    tgt::vec4 zero = tgt::vec4(0.f, 0.f, 0.f, 1.f);
    tgt::vec4 x1 = tgt::vec4(-1.f, 0.f, 0.f, 1.f);
    tgt::vec4 x2 = tgt::vec4( 1.f, 0.f, 0.f, 1.f);
    tgt::vec4 x3 = tgt::vec4(0.8f, 0.1f, 0.f, 1.f);
    tgt::vec4 x4 = tgt::vec4(0.8f, 0.f, 0.1f, 1.f);
    tgt::vec4 x5 = tgt::vec4(0.8f,-0.1f, 0.f, 1.f);
    tgt::vec4 x6 = tgt::vec4(0.8f, 0.f,-0.1f, 1.f);

    tgt::vec4 y1 = tgt::vec4( 0.f,-1.f, 0.f, 1.f);
    tgt::vec4 y2 = tgt::vec4( 0.f, 1.f, 0.f, 1.f);
    tgt::vec4 y3 = tgt::vec4(0.1f, 0.8f, 0.f, 1.f);
    tgt::vec4 y4 = tgt::vec4(0.f, 0.8f, 0.1f, 1.f);
    tgt::vec4 y5 = tgt::vec4(-0.1f, 0.8f, 0.f, 1.f);
    tgt::vec4 y6 = tgt::vec4(0.f, 0.8f,-0.1f, 1.f);

    tgt::vec4 z1 = tgt::vec4( 0.f, 0.f,-1.f, 1.f);
    tgt::vec4 z2 = tgt::vec4( 0.f, 0.f, 1.f, 1.f);
    tgt::vec4 z3 = tgt::vec4(0.1f, 0.f, 0.8f, 1.f);
    tgt::vec4 z4 = tgt::vec4(0.f, 0.1f, 0.8f, 1.f);
    tgt::vec4 z5 = tgt::vec4(-0.1f, 0.f, 0.8f, 1.f);
    tgt::vec4 z6 = tgt::vec4(0.f,-0.1f, 0.8f, 1.f);

    // define array
    tgt::vec4 arrowlines[24] = { x1, zero, x2, x3, x4, x5, x6, x3, y1, zero, y2, y3, y4, y5, y6, y3, z1, zero, z2, z3, z4, z5, z6, z3 };

    // setup VBO
    int SIZE_FPOS = sizeof(arrowlines);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_[2]);
    glBufferData(GL_ARRAY_BUFFER, SIZE_FPOS, (GLvoid*)arrowlines, GL_STATIC_DRAW);

    // setup VAO
    glBindVertexArray(vao_[2]);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);

    //clean-up
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void CameraPositionRenderer::renderFrustum(){
    // bind VAO
    glBindVertexArray(vao_[0]);

    // set uniforms
    myShaderProp_.getShader()->setUniform("color_", frustLineColorProp_.get());
    myShaderProp_.getShader()->setUniform("invTransform_", pvInv_);
    myShaderProp_.getShader()->setUniform("scale_", 1.0f);

    // draw frustum
    glDrawArrays(GL_LINES, 0, 24);

    // draw planes if necessary
    if (renderPlanesProp_.get()){
        tgt::vec4 col = frustLineColorProp_.get();
        col[3] = opacPlanesProp_.get();
        myShaderProp_.getShader()->setUniform("color_", col);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        glDrawArrays(GL_TRIANGLE_STRIP, 4, 4);
    }

    //clean-up
    glBindVertexArray(0);
}

void CameraPositionRenderer::renderCameraModel(){
    // bind VAO
    glBindVertexArray(vao_[1]);

    // set uniforms
    tgt::vec4 col = colorModelProp_.get();
    col[3] = opacModelProp_.get();
    myShaderProp_.getShader()->setUniform("color_", col);
    myShaderProp_.getShader()->setUniform("invTransform_", vInv_);
    

    // draw model...
    if (optionModelProp_.isSelected("model")){
        myShaderProp_.getShader()->setUniform("scale_", scalingModelProp_.get());
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 10);
    }
    else if (optionModelProp_.isSelected("lines")){
        myShaderProp_.getShader()->setUniform("scale_", scalingModelProp_.get());
        glDrawArrays(GL_LINES, 0, 9);
    }
    else {
        myShaderProp_.getShader()->setUniform("scale_", 1.0f);
        glDrawArrays(GL_POINTS, 10, 1);
    }

    //clean-up
    glBindVertexArray(0);
}

void CameraPositionRenderer::renderArrows(){

    // bind VAO
    glBindVertexArray(vao_[2]);
    
    // set transformation matrix
    myShaderProp_.getShader()->setUniform("invTransform_", vInv_);

    // draw direction arrows...
    if (optionArrowsProp_.isSelected("arrowsxyz") || optionArrowsProp_.isSelected("arrowsxy")){
        myShaderProp_.getShader()->setUniform("scale_", scalingArrowsProp_.get());
        myShaderProp_.getShader()->setUniform("color_", tgt::vec4(1.f, 0.f, 0.f, 1.f));
        glDrawArrays(GL_LINE_STRIP, 0, 3);
        glDrawArrays(GL_TRIANGLE_FAN, 2, 6);
        myShaderProp_.getShader()->setUniform("color_", tgt::vec4(0.f, 1.f, 0.f, 1.f));
        glDrawArrays(GL_LINE_STRIP, 8, 3);
        glDrawArrays(GL_TRIANGLE_FAN, 10, 6);
        if (optionArrowsProp_.isSelected("arrowsxyz")){
            myShaderProp_.getShader()->setUniform("color_", tgt::vec4(0.f, 0.f, 1.f, 1.f));
            glDrawArrays(GL_LINE_STRIP, 16, 3);
            glDrawArrays(GL_TRIANGLE_FAN, 18, 6);
        }
    }
    // ...or draw coordinate system
    else if (optionArrowsProp_.isSelected("coords")){
        myShaderProp_.getShader()->setUniform("color_", tgt::vec4(1.f, 0.f, 0.f, 1.f));
        myShaderProp_.getShader()->setUniform("scale_", xLengthProp_.get());
        glDrawArrays(GL_LINES, 1, 2);
        myShaderProp_.getShader()->setUniform("color_", tgt::vec4(0.f, 1.f, 0.f, 1.f));
        myShaderProp_.getShader()->setUniform("scale_", yLengthProp_.get());
        glDrawArrays(GL_LINES, 9, 2);
        myShaderProp_.getShader()->setUniform("color_", tgt::vec4(0.f, 0.f, 1.f, 1.f));
        myShaderProp_.getShader()->setUniform("scale_", zLengthProp_.get());
        glDrawArrays(GL_LINES, 17, 2);
    }

    //clean-up
    glBindVertexArray(0);
}

// main rendering function
void CameraPositionRenderer::render(){

        // activate shader(s)
        myShaderProp_.getShader()->activate();

        // get secondary camera
        tgt::Camera cam = displayCameraProp_.get();
        setGlobalShaderParameters(myShaderProp_.getShader(), &cam);

        // calculate and pass transformation matricess
        tgt::mat4 projectionView = displayCameraProp_.get().getProjectionMatrix(sizeReferencePort_.getSize())*displayCameraProp_.get().getViewMatrix();
        tgt::mat4 viewMatrix = displayCameraProp_.get().getViewMatrix();
        projectionView.invert(pvInv_);
        viewMatrix.invert(vInv_);
        tgt::mat4 cameraMatrix = camera_.getProjectionMatrix(viewport_)*camera_.getViewMatrix();
        myShaderProp_.getShader()->setUniform("cMat_", cameraMatrix);

        //render
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        if (renderFrustumProp_.get()) {
            glLineWidth(frustLineWidthProp_.get());
            renderFrustum();
            glLineWidth(1.0f);
        }

        if (!optionArrowsProp_.isSelected("off")){
            renderArrows();
        }

        if (!optionModelProp_.isSelected("off")) {
            glLineWidth(lineModelProp_.get());
            renderCameraModel();
            glLineWidth(1.0f);
        }

        //reset blending settings
        glBlendFunc(GL_ONE, GL_ZERO);
        glDisable(GL_BLEND);

        // deactivate shader(s)
        myShaderProp_.getShader()->deactivate();
        LGL_ERROR;
}

void CameraPositionRenderer::onRenderFrustumChanged(){
    if (renderFrustumProp_.get()){
        frustLineColorProp_.setVisibleFlag(true);
        frustLineWidthProp_.setVisibleFlag(true);
        renderPlanesProp_.setVisibleFlag(true);
        opacPlanesProp_.setVisibleFlag(true);
    }
    else{
        frustLineColorProp_.setVisibleFlag(false);
        frustLineWidthProp_.setVisibleFlag(false);
        renderPlanesProp_.setVisibleFlag(false);
        opacPlanesProp_.setVisibleFlag(false);
    }
}

void CameraPositionRenderer::onRenderPlanesChanged(){
    if (renderPlanesProp_.get() && renderFrustumProp_.get()){
        opacPlanesProp_.setVisibleFlag(true);
    }
    else{
        opacPlanesProp_.setVisibleFlag(false);
    }
}

void CameraPositionRenderer::onRenderModelChanged(){
    if (optionModelProp_.isSelected("model")){
        opacModelProp_.setVisibleFlag(true);
        colorModelProp_.setVisibleFlag(true);
        lineModelProp_.setVisibleFlag(false);
        scalingModelProp_.setVisibleFlag(true);
    }
    else if (optionModelProp_.isSelected("lines")){
        opacModelProp_.setVisibleFlag(false);
        colorModelProp_.setVisibleFlag(true);
        lineModelProp_.setVisibleFlag(true);
        scalingModelProp_.setVisibleFlag(true);
    }
    else if (optionModelProp_.isSelected("dot")){
        opacModelProp_.setVisibleFlag(true);
        colorModelProp_.setVisibleFlag(true);
        lineModelProp_.setVisibleFlag(false);
        scalingModelProp_.setVisibleFlag(false);
    }
    else if (optionModelProp_.isSelected("off")){
        opacModelProp_.setVisibleFlag(false);
        colorModelProp_.setVisibleFlag(false);
        lineModelProp_.setVisibleFlag(false);
        scalingModelProp_.setVisibleFlag(false);
    }
}

void CameraPositionRenderer::onRenderArrowsChanged(){
    if (optionArrowsProp_.isSelected("off")){
        scalingArrowsProp_.setVisibleFlag(false);
        xLengthProp_.setVisibleFlag(false);
        yLengthProp_.setVisibleFlag(false);
        zLengthProp_.setVisibleFlag(false);
    }
    else if (optionArrowsProp_.isSelected("coords")){
        scalingArrowsProp_.setVisibleFlag(false);
        xLengthProp_.setVisibleFlag(true);
        yLengthProp_.setVisibleFlag(true);
        zLengthProp_.setVisibleFlag(true);
    }
    else{
        scalingArrowsProp_.setVisibleFlag(true);
        xLengthProp_.setVisibleFlag(false);
        yLengthProp_.setVisibleFlag(false);
        zLengthProp_.setVisibleFlag(false);
    }
}

} //namespace
