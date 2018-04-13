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

#include "slicepositionrenderer.h"

#include "tgt/glmath.h"
#include "tgt/quadric.h"

namespace voreen {

using tgt::vec4;
using tgt::vec3;

SlicePositionRenderer::SlicePositionRenderer()
    : GeometryRendererBase()
    //port
    , inport_(Port::INPORT, "volume", "Volume Input")
    //properties
    , enable_("enable", "Enable", true)
        //x slice
    , renderXSlice_("renderXSlice", "Render X Slice", true)
    , xSliceIndexProp_("xSliceIndex", "X Slice Number", 0, 0, 10000)
    , xColor_("xColor", "X Color", tgt::vec4(1.0f, 0.0f, 0.0f, 1.0f))
    , lineXWidth_("lineXWidth","Line X Width",1.f,1.f,10.f)
    , alphaFactorXPlane_("alphaFactorXPlane"," X Plane Alpha", 0.f, 0.f, 1.f)
        //y slice
    , renderYSlice_("renderYSlice", "Render Y Slice", true)
    , ySliceIndexProp_("ySliceIndex", "Y Slice Number", 0, 0, 10000)
    , yColor_("yColor", "Y Color", tgt::vec4(0.0f, 1.0f, 0.0f, 1.0f))
    , lineYWidth_("lineYWidth","Line Y Width",1.f,1.f,10.f)
    , alphaFactorYPlane_("alphaFactorYPlane"," Y Plane Alpha", 0.f, 0.f, 1.f)
        //z slice
    , renderZSlice_("renderZSlice", "Render Z Slice", true)
    , zSliceIndexProp_("zSliceIndex", "Z Slice Number", 0, 0, 10000)
    , lineZWidth_("lineZWidth","Line Z Width",1.f,1.f,10.f)
    , zColor_("zColor", "Z Color", tgt::vec4(0.0f, 0.0f, 1.0f, 1.0f))
    , alphaFactorZPlane_("alphaFactorZPlane"," Z Plane Alpha", 0.f, 0.f, 1.f)
        //Shaders (For both the same files, but we will set different headers)
    , planeShaderOpaque_("planeShaderOpaque", "Plane shader opaque", "slicepositionrenderer.frag", "slicepositionrenderer.vert")
    , planeShaderTransparent_("planeShaderTransparent", "Plane shader transparent", "slicepositionrenderer.frag", "slicepositionrenderer.vert")
    , vao_(0)
    , vbo_(0)
{
    //ports
    addPort(inport_);
    //properties
    addProperty(enable_);
        // x slice
    addProperty(renderXSlice_);
    renderXSlice_.setGroupID("x");
    addProperty(xSliceIndexProp_);
    xSliceIndexProp_.setGroupID("x");
    addProperty(xColor_);
    xColor_.setGroupID("x");
    addProperty(lineXWidth_);
    lineXWidth_.setGroupID("x");
    addProperty(alphaFactorXPlane_);
    alphaFactorXPlane_.setGroupID("x");
    setPropertyGroupGuiName("x","X Slice");
        // y slice
    addProperty(renderYSlice_);
    renderYSlice_.setGroupID("y");
    addProperty(ySliceIndexProp_);
    ySliceIndexProp_.setGroupID("y");
    addProperty(yColor_);
    yColor_.setGroupID("y");
    addProperty(lineYWidth_);
    lineYWidth_.setGroupID("y");
    addProperty(alphaFactorYPlane_);
    alphaFactorYPlane_.setGroupID("y");
    setPropertyGroupGuiName("y","Y Slice");
        // z slice
    addProperty(renderZSlice_);
    renderZSlice_.setGroupID("z");
    addProperty(zSliceIndexProp_);
    zSliceIndexProp_.setGroupID("z");
    addProperty(zColor_);
    zColor_.setGroupID("z");
    addProperty(lineZWidth_);
    lineZWidth_.setGroupID("z");
    addProperty(alphaFactorZPlane_);
    alphaFactorZPlane_.setGroupID("z");
    setPropertyGroupGuiName("z","Z Slice");

    addProperty(planeShaderOpaque_);
    addProperty(planeShaderTransparent_);

    ON_PROPERTY_CHANGE(renderXSlice_,SlicePositionRenderer,togglePropertyVisibility);
    ON_PROPERTY_CHANGE(renderYSlice_,SlicePositionRenderer,togglePropertyVisibility);
    ON_PROPERTY_CHANGE(renderZSlice_,SlicePositionRenderer,togglePropertyVisibility);
    togglePropertyVisibility();
}

void SlicePositionRenderer::initialize() {
    GeometryRendererBase::initialize();

    planeShaderTransparent_.setHeader("#define USE_TRANSPARENCY \n");
    planeShaderTransparent_.rebuild();
    planeShaderOpaque_.rebuild();

    glGenBuffers(1, &vbo_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(tgt::vec3)*4, NULL, GL_DYNAMIC_DRAW);

    glGenVertexArrays(1, &vao_);
    glBindVertexArray(vao_);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glBindVertexArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}
void SlicePositionRenderer::deinitialize() {
    glDeleteBuffers(1, &vbo_);
    glDeleteVertexArrays(1, &vao_);

    GeometryRendererBase::deinitialize();
}
void SlicePositionRenderer::process() {}

tgt::Bounds SlicePositionRenderer::getBoundingBox() const {
    if (inport_.hasData())
        return inport_.getData()->getBoundingBox().getBoundingBox();

    return GeometryRendererBase::getBoundingBox();
}

void SlicePositionRenderer::setVBOData(const std::array<tgt::vec3, 4>& vertices) {
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(tgt::vec3)*4, vertices.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

bool SlicePositionRenderer::usesTransparency() const {
    return (renderXSlice_.get() && ((xColor_.get().a > 0 && xColor_.get().a < 1) || (alphaFactorXPlane_.get() > 0 && alphaFactorXPlane_.get() < 1)))
        || (renderYSlice_.get() && ((yColor_.get().a > 0 && yColor_.get().a < 1) || (alphaFactorYPlane_.get() > 0 && alphaFactorYPlane_.get() < 1)))
        || (renderZSlice_.get() && ((zColor_.get().a > 0 && zColor_.get().a < 1) || (alphaFactorZPlane_.get() > 0 && alphaFactorZPlane_.get() < 1)));
}
void SlicePositionRenderer::renderTransparent() {
    render(planeShaderTransparent_.getShader());
}
void SlicePositionRenderer::render() {
    render(planeShaderOpaque_.getShader());
}
void SlicePositionRenderer::render(tgt::Shader* shader) {
    if(!shader) {
        LWARNING("No Shader.");
    }
    if (!inport_.isReady() || !enable_.get())
        return;

    MatStack.pushMatrix();
    MatStack.multMatrix(inport_.getData()->getPhysicalToWorldMatrix());

    tgt::vec3 geomLlf = inport_.getData()->getLLF();
    tgt::vec3 geomUrb = inport_.getData()->getURB();

    // Convert voxel coordinates to physical coordinates for slice positions
    tgt::vec3 slicePositionVector(
                static_cast<float>(xSliceIndexProp_.get()),
                static_cast<float>(ySliceIndexProp_.get()),
                static_cast<float>(zSliceIndexProp_.get())
            );
    slicePositionVector = inport_.getData()->getVoxelToPhysicalMatrix() * slicePositionVector;
    float xSlice = slicePositionVector.x;
    float ySlice = slicePositionVector.y;
    float zSlice = slicePositionVector.z;

    if(lineXWidth_.get() != 1.f || lineYWidth_.get() != 1.f || lineZWidth_.get() != 1.f)
        glEnable(GL_LINE_SMOOTH);

    shader->activate();
    glBindVertexArray(vao_);
    shader->setIgnoreUniformLocationError(true);
    shader->setIgnoreUnsetUniform("headPointerImage_"); // Only transparent shader uses this uniform.
    shader->setIgnoreUniformLocationError(false);
    shader->setUniform("transform_", MatStack.getProjectionMatrix()*MatStack.getModelViewMatrix());
    // X Slice
    if (xColor_.get().a > 0.f && renderXSlice_.get()) {
        tgt::vec4 color = xColor_.get();
        shader->setUniform("color_", color);
        glLineWidth(lineXWidth_.get());

        std::array<tgt::vec3, 4> tmp = {
                tgt::vec3(xSlice, geomUrb.y, geomUrb.z),
                tgt::vec3(xSlice, geomLlf.y, geomUrb.z),
                tgt::vec3(xSlice, geomLlf.y, geomLlf.z),
                tgt::vec3(xSlice, geomUrb.y, geomLlf.z)
                };
        setVBOData(tmp);
        glDrawArrays(GL_LINE_LOOP, 0, 4);
        if (alphaFactorXPlane_.get() > 0.f) {
            color.w *= alphaFactorXPlane_.get();
            shader->setUniform("color_", color);
            glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
        }
    }
    // Y Slice
    if (yColor_.get().a > 0.f && renderYSlice_.get()) {
        tgt::vec4 color = yColor_.get();
        shader->setUniform("color_", color);
        glLineWidth(lineYWidth_.get());

        std::array<tgt::vec3, 4> tmp = {
                tgt::vec3(geomUrb.x, ySlice, geomUrb.z),
                tgt::vec3(geomLlf.x, ySlice, geomUrb.z),
                tgt::vec3(geomLlf.x, ySlice, geomLlf.z),
                tgt::vec3(geomUrb.x, ySlice, geomLlf.z)
                };
        setVBOData(tmp);
        glDrawArrays(GL_LINE_LOOP, 0, 4);
        if (alphaFactorYPlane_.get() > 0.f) {
            color.w *= alphaFactorYPlane_.get();
            shader->setUniform("color_", color);
            glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
        }
    }
    // Z Slice
    if (zColor_.get().a > 0.f && renderZSlice_.get()) {
        tgt::vec4 color = zColor_.get();
        shader->setUniform("color_", color);
        glLineWidth(lineZWidth_.get());

        std::array<tgt::vec3, 4> tmp = {
                tgt::vec3(geomUrb.x, geomUrb.y, zSlice),
                tgt::vec3(geomLlf.x, geomUrb.y, zSlice),
                tgt::vec3(geomLlf.x, geomLlf.y, zSlice),
                tgt::vec3(geomUrb.x, geomLlf.y, zSlice)
                };
        setVBOData(tmp);
        glDrawArrays(GL_LINE_LOOP, 0, 4);
        if (alphaFactorZPlane_.get() > 0.f) {
            color.w *= alphaFactorZPlane_.get();
            shader->setUniform("color_", color);
            glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
        }
    }
    glBindVertexArray(0);
    shader->deactivate();

    glDisable(GL_LINE_SMOOTH);

    MatStack.popMatrix();

    glLineWidth(1.f);
    glColor4f(1.f, 1.f, 1.f, 1.f);
}

void SlicePositionRenderer::invalidate(int inv) {
    GeometryRendererBase::invalidate(inv);

    if (inport_.hasChanged() && inport_.hasData()) {
        tgt::ivec3 numSlices = inport_.getData()->getDimensions();

        xSliceIndexProp_.setMaxValue(numSlices.x-1);
        ySliceIndexProp_.setMaxValue(numSlices.y-1);
        zSliceIndexProp_.setMaxValue(numSlices.z-1);
    }
}

void SlicePositionRenderer::togglePropertyVisibility() {
    /*xColor_.setVisibleFlag(renderXSlice_.get());
    xSliceIndexProp_.setVisibleFlag(renderXSlice_.get());
    lineXWidth_.setVisibleFlag(renderXSlice_.get());
    alphaFactorXPlane_.setVisibleFlag(renderXSlice_.get());

    yColor_.setVisibleFlag(renderYSlice_.get());
    ySliceIndexProp_.setVisibleFlag(renderYSlice_.get());
    lineYWidth_.setVisibleFlag(renderYSlice_.get());
    alphaFactorYPlane_.setVisibleFlag(renderYSlice_.get());

    zColor_.setVisibleFlag(renderZSlice_.get());
    zSliceIndexProp_.setVisibleFlag(renderZSlice_.get());
    lineZWidth_.setVisibleFlag(renderZSlice_.get());
    alphaFactorZPlane_.setVisibleFlag(renderZSlice_.get());*/
}

void SlicePositionRenderer::adjustPropertiesToInput() {
    if(inport_.getData()) {
        tgt::svec3 dim = inport_.getData()->getDimensions();
        xSliceIndexProp_.setMaxValue(dim.x - 1);
        ySliceIndexProp_.setMaxValue(dim.y - 1);
        zSliceIndexProp_.setMaxValue(dim.z - 1);
    }
}

} // namespace

