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

#include "flowslicerenderer3d.h"
#include "voreen/core/interaction/camerainteractionhandler.h"

namespace voreen {

FlowSliceRenderer3D::FlowSliceRenderer3D()
    : FlowSliceRenderer(),
    useXYSliceProp_("useXYSliceProp", "render xy-slice:", true),
    useXZSliceProp_("useXZSliceProp", "render xz-slice:", false),
    useZYSliceProp_("useZYSliceProp", "render zy slice:", false),
    sliceNoXYProp_("sliceNoXYProp", "xy-slice number:", 1, 1, 100),
    sliceNoXZProp_("sliceNoXZProp", "xz-slice number:", 1, 1, 100),
    sliceNoZYProp_("sliceNoZYProp", "zy-slice number:", 1, 1, 100),
    camProp_("camera", "Camera", tgt::Camera(tgt::vec3(0.0f, 0.0f, 3.5f), tgt::vec3(0.0f, 0.0f, 0.0f), tgt::vec3(0.0f, 1.0f, 0.0f))),
    cameraHandler_(0),
    textureXY_(0),
    textureXZ_(0),
    textureZY_(0),
    rebuildTextureXY_(true),
    rebuildTextureXZ_(true),
    rebuildTextureZY_(true),
    alignment_(FlowSliceRenderer::PLANE_XY)
{
    zoom_ = FlowSliceRenderer::ZOOM_1X;

    MemberFunctionCallback<FlowSliceRenderer3D> invalidateXY(this,
        &FlowSliceRenderer3D::invalidateXYTexture);
    MemberFunctionCallback<FlowSliceRenderer3D> invalidateXZ(this,
        &FlowSliceRenderer3D::invalidateXZTexture);
    MemberFunctionCallback<FlowSliceRenderer3D> invalidateZY(this,
        &FlowSliceRenderer3D::invalidateZYTexture);

    useXYSliceProp_.onChange(invalidateXY);
    sliceNoXYProp_.onChange(invalidateXY);
    useXZSliceProp_.onChange(invalidateXZ);
    sliceNoXZProp_.onChange(invalidateXZ);
    useZYSliceProp_.onChange(invalidateZY);
    sliceNoZYProp_.onChange(invalidateZY);

    addProperty(useBoundingBoxProp_);
    addProperty(boundingBoxColorProp_);
    addProperty(useCoordinateAxisProp_);
    addProperty(useXYSliceProp_);
    addProperty(sliceNoXYProp_);
    addProperty(useXZSliceProp_);
    addProperty(sliceNoXZProp_);
    addProperty(useZYSliceProp_);
    addProperty(sliceNoZYProp_);
    addProperty(camProp_);

    cameraHandler_ = new CameraInteractionHandler("cameraHandler", "Camera Handler", &camProp_);
    addInteractionHandler(cameraHandler_);
}

FlowSliceRenderer3D::~FlowSliceRenderer3D() {
    delete textureXY_;
    delete textureXZ_;
    delete textureZY_;
    delete cameraHandler_;
}

void FlowSliceRenderer3D::process() {
    bool handleChanged = false;
    const VolumeFlow3D* volFlow = checkVolumeHandleForFlowData<VolumeFlow3D>(volInport_.getData(), handleChanged);
    if (volFlow == 0) {
        LERROR("process(): supplied volume does not contain compatible flow data!\n");
        return;
    }

    tgt::vec2 viewportSize = imgOutport_.getSize();
    imgOutport_.activateTarget("FlowSliceRenderer3D::process()");

    const Flow3D& flow = volFlow->getFlow3D();
    flowDimensions_ = flow.dimensions_;

    if (handleChanged == true) {
        updateNumSlices();  // validate the currently set values and adjust them if necessary
        rebuildTexture_ = true;
    }

    // important: save current camera state before using the processor's camera or
    // successive processors will use those settings!
    //
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();

    camProp_.look(imgOutport_.getSize());

    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

    if (useBoundingBoxProp_.get() == true)
        renderBoundingBox(boundingBoxColorProp_.get());

    if (useCoordinateAxisProp_.get() == true)
        renderCoordinateAxis();

    if (useXYSliceProp_.get() == true) {
        const float sliceNoXY = static_cast<float>(sliceNoXYProp_.get() - 1);
        renderSlice(flow, sliceNoXY, viewportSize, FlowSliceRenderer::PLANE_XY);
        rebuildTextureXY_ = false;
    }

    if (useXZSliceProp_.get() == true) {
        const float sliceNoXZ = static_cast<float>(sliceNoXZProp_.get() - 1);
        renderSlice(flow, sliceNoXZ, viewportSize, FlowSliceRenderer::PLANE_XZ);
        rebuildTextureXZ_ = false;
    }

    if (useZYSliceProp_.get() == true) {
        const float sliceNoZY = static_cast<float>(sliceNoZYProp_.get() - 1);
        renderSlice(flow, sliceNoZY, viewportSize, FlowSliceRenderer::PLANE_ZY);
        rebuildTextureZY_ = false;
    }

    rebuildTexture_ = false;

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();

    imgOutport_.deactivateTarget();
}

// protected methods
//

std::string FlowSliceRenderer3D::generateShaderHeader() {
    std::ostringstream oss;
    oss << FlowSliceRenderer::generateShaderHeader();
    oss << "#define PLANE " << static_cast<size_t>(alignment_) << std::endl;
    return oss.str();
}

// private methods
//

void FlowSliceRenderer3D::renderSlice(const Flow3D& flow3D, const float sliceNo,
                                      const tgt::vec2& viewportSize,
                                      const SliceAlignment& sliceAlignment)
{
    const tgt::vec2 thresholds = (flow3D.maxMagnitude_ * thresholdProp_.get() / 100.0f);

    // get permutation of vector components for that slice alignment, permute the
    // components of vertices for the textured quad and map them into the 3D flow
    // volume's boundaries
    //
    permutation_ = getCoordinatePermutation(sliceAlignment);
    tgt::vec2 textureSize = (flow3D.getFlowSliceDimensions(permutation_) - tgt::ivec2(1));
    tgt::vec3 ll = mapToFlowBoundingBox(
        permuteComponents(tgt::vec3(0.0f, 0.0f, sliceNo), permutation_));
    tgt::vec3 lr = mapToFlowBoundingBox(
        permuteComponents(tgt::vec3(textureSize.x, 0.0f, sliceNo), permutation_));
    tgt::vec3 ul = mapToFlowBoundingBox(
        permuteComponents(tgt::vec3(0.0f, textureSize.y, sliceNo), permutation_));
    tgt::vec3 ur = mapToFlowBoundingBox(
        permuteComponents(tgt::vec3(textureSize, sliceNo), permutation_));

    tgt::Texture** textureAddr = 0;
    bool prevRebuildState = rebuildTexture_;
    switch (sliceAlignment) {
        default:
        case FlowSliceRenderer::PLANE_XY:
            textureAddr = &textureXY_;
            rebuildTexture_ = rebuildTextureXY_ || rebuildTexture_;
            break;
        case FlowSliceRenderer::PLANE_XZ:
            textureAddr = &textureXZ_;
            rebuildTexture_ = rebuildTextureXZ_ || rebuildTexture_;
            break;
        case FlowSliceRenderer::PLANE_ZY:
            textureAddr = &textureZY_;
            rebuildTexture_ = rebuildTextureZY_ || rebuildTexture_;
            break;
    }

    if ((rebuildTexture_ == true) && (textureAddr != 0)) {
        if (sliceAlignment != alignment_) {
            alignment_ = sliceAlignment;    // used for shader header generation...
            rebuildShader();
        }

        delete *textureAddr;
        *textureAddr = 0;
        switch (techniqueProp_->getValue()) {
            default:
            case TECHNIQUE_COLOR_CODING:
            case TECHNIQUE_COLOR_CODING_PROJECTED:
                *textureAddr = renderColorCodingTexture(flow3D, static_cast<size_t>(sliceNo),
                    textureSize, privatePort1_, thresholds,
                    (techniqueProp_->getValue() == TECHNIQUE_COLOR_CODING_PROJECTED));
                break;

            case TECHNIQUE_ARROW_PLOT_RAND:
            case TECHNIQUE_ARROW_PLOT_GRID:
                *textureAddr = renderArrowPlotTexture(flow3D.extractSlice(permutation_,
                    static_cast<size_t>(sliceNo)), textureSize, privatePort1_, thresholds);
                break;

            case TECHNIQUE_SPOTNOISE:
            case TECHNIQUE_SPOTNOISE_PROJECTED:
                {
                    std::vector<RenderPort*> tempPorts;
                    tempPorts.push_back(&privatePort1_);
                    tempPorts.push_back(&privatePort2_);

                    *textureAddr = renderSpotNoiseTexture(flow3D, static_cast<size_t>(sliceNo),
                        textureSize, viewportSize, tempPorts,
                        (techniqueProp_->getValue() == TECHNIQUE_SPOTNOISE_PROJECTED));
                }
                break;

            case TECHNIQUE_INTEGRATE_DRAW:
            case TECHNIQUE_FAST_LIC:
            case TECHNIQUE_INTEGRATE_DRAW_PROJECTED:
            case TECHNIQUE_FAST_LIC_PROJECTED:
                *textureAddr = renderFlowTexture(flow3D, static_cast<size_t>(sliceNo), 1,
                    techniqueProp_->getValue(), thresholds);
                break;
        }   // switch

        if (*textureAddr != 0)
            (*textureAddr)->uploadTexture();
        else {
            std::cout << "Error: FlowSliceRenderer3D::renderSlice(): slice creation failed\n";
            return;
        }
    } else if (textureAddr == 0) {
        std::cout << "Fatal Error: FlowSliceRenderer3D::renderSlice(): textureAddr == 0! Cannot proceed.\n";
        return;
    }

    imgOutport_.activateTarget("FlowSliceRenderer3D::process()");

    glActiveTexture(GL_TEXTURE0);
    glEnable(GL_TEXTURE_2D);
    if ((*textureAddr != 0))
        glBindTexture(GL_TEXTURE_2D, (*textureAddr)->getId());
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

    glBegin(GL_QUADS);
        glMultiTexCoord2f(GL_TEXTURE0, 0.0f, 0.0f);
        glVertex3fv(ll.elem);

        glMultiTexCoord2f(GL_TEXTURE0, 1.0f, 0.0f);
        glVertex3fv(lr.elem);

        glMultiTexCoord2f(GL_TEXTURE0, 1.0f, 1.0f);
        glVertex3fv(ur.elem);

        glMultiTexCoord2f(GL_TEXTURE0, 0.0f, 1.0f);
        glVertex3fv(ul.elem);
    glEnd();

    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glDisable(GL_TEXTURE_2D);
    rebuildTexture_ = prevRebuildState;
    imgOutport_.deactivateTarget();
}

void FlowSliceRenderer3D::updateNumSlices() {
    // set number of slice for axial plane (along z-axis)
    //
    size_t numSlicesXY = flowDimensions_.z;
    if (numSlicesXY > static_cast<size_t>(sliceNoXYProp_.getMinValue())) {
        sliceNoXYProp_.setMaxValue(static_cast<int>(numSlicesXY));
        if (static_cast<size_t>(sliceNoXYProp_.get()) > numSlicesXY)
            sliceNoXYProp_.set(static_cast<int>(numSlicesXY) / 2);
        sliceNoXYProp_.updateWidgets();
    }

    // set number of slices for coronal plane (along y-axis)
    //
    size_t numSlicesXZ = flowDimensions_.y;
    if (numSlicesXZ > static_cast<size_t>(sliceNoXZProp_.getMinValue())) {
        sliceNoXZProp_.setMaxValue(static_cast<int>(numSlicesXZ));
        if (static_cast<size_t>(sliceNoXZProp_.get()) > numSlicesXZ)
            sliceNoXZProp_.set(static_cast<int>(numSlicesXZ) / 2);
        sliceNoXZProp_.updateWidgets();
    }

    // set number of slices for sagittal plane (along x-axis)
    //
    size_t numSlicesZY = flowDimensions_.x;
    if (numSlicesZY > static_cast<size_t>(sliceNoZYProp_.getMinValue())) {
        sliceNoZYProp_.setMaxValue(static_cast<int>(numSlicesZY));
        if (static_cast<size_t>(sliceNoZYProp_.get()) > numSlicesZY)
            sliceNoZYProp_.set(static_cast<int>(numSlicesZY) / 2);
        sliceNoZYProp_.updateWidgets();
    }
}

inline tgt::vec3 FlowSliceRenderer3D::permuteComponents(const tgt::vec3& input, const tgt::ivec3& permutation) {
    return tgt::vec3(input[permutation.x], input[permutation.y], input[permutation.z]);
}

}   // namespace
