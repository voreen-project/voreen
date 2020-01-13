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

#include "flowslicerenderer2d.h"

namespace voreen {

FlowSliceRenderer2D::FlowSliceRenderer2D()
    : FlowSliceRenderer(),
    alignmentProp_(0),
    sliceIndexProp_("sliceIndex", "slice number: ", 1, 1, 100),
    textureZoomProp_(0),
    sliceTexture_(0),
    alignment_(PLANE_XY),
    numSlices_(0)
{

    alignmentProp_ = new OptionProperty<SliceAlignment>("sliceAlignment", "slice alignment: ");
    alignmentProp_->addOption("xy-plane", "xy-plane", PLANE_XY);
    alignmentProp_->addOption("xz-plane", "xz-plane", PLANE_XZ);
    alignmentProp_->addOption("zy-plane", "zy-plane", PLANE_ZY);
    alignmentProp_->onChange(
        MemberFunctionCallback<FlowSliceRenderer2D>(this, &FlowSliceRenderer2D::onSliceAlignmentChange));

    textureZoomProp_ = new OptionProperty<TextureZoom>("textureZoom", "texture zoom: ");
    textureZoomProp_->addOption("full size", "full size", ZOOM_FULL_SIZE);
    textureZoomProp_->addOption("1x", "1x", ZOOM_1X);
    textureZoomProp_->addOption("2x", "2x", ZOOM_2X);
    textureZoomProp_->addOption("4x", "4x", ZOOM_4X);
    textureZoomProp_->addOption("8x", "8x", ZOOM_8X);
    textureZoomProp_->onChange(
        MemberFunctionCallback<FlowSliceRenderer2D>(this, &FlowSliceRenderer2D::onTextureZoomChange));

    MemberFunctionCallback<FlowSliceRenderer2D> invalidateAction(this, &FlowSliceRenderer2D::invalidateTexture);
    sliceIndexProp_.onChange(invalidateAction);

    addProperty(textureZoomProp_);
    addProperty(alignmentProp_);
    addProperty(sliceIndexProp_);

    onTextureZoomChange();
}

FlowSliceRenderer2D::~FlowSliceRenderer2D() {
    delete alignmentProp_;
    delete textureZoomProp_;
    delete sliceTexture_;
}

void FlowSliceRenderer2D::process() {
    bool handleChanged = false;
    const VolumeFlow3D* volFlow = checkVolumeHandleForFlowData<VolumeFlow3D>(volInport_.getData(), handleChanged);
    if (volFlow == 0) {
        LERROR("process(): supplied volume does not contain compatible flow data!\n");
        return;
    }

    tgt::vec2 viewportSize = imgOutport_.getSize();
    const Flow3D& flow = volFlow->getFlow3D();
    flowDimensions_ = flow.dimensions_;

    if (handleChanged == true) {
        updateNumSlices();  // validate the currently set values and adjust them if necessary
        rebuildTexture_ = true;
    }

    const GLint sliceNo = static_cast<GLint>(sliceIndexProp_.get() - 1);
    tgt::vec2 thresholds(flow.maxMagnitude_ * thresholdProp_.get() / 100.0f);

    // set slice's width and height according to currently slice alignment
    //
    tgt::vec2 texDim2D = static_cast<tgt::vec2>(flow.getFlowSliceDimensions(permutation_));
    float scaling = getTextureScalingFactor(viewportSize, texDim2D);
    tgt::vec2 outputSize(texDim2D * scaling);  // size of the quad on which the texture will be rendered
    tgt::vec2 slicePos = (viewportSize - outputSize) / 2.0f;    // offset for the quad on which the texture will be rendered

    if (rebuildTexture_ == true) {
        delete sliceTexture_;
        sliceTexture_ = 0;

        switch (techniqueProp_->getValue()) {
            default:
            case TECHNIQUE_COLOR_CODING:
            case TECHNIQUE_COLOR_CODING_PROJECTED:
                sliceTexture_ = renderColorCodingTexture(flow, sliceNo, outputSize, privatePort1_,
                    thresholds, (techniqueProp_->getValue() == TECHNIQUE_COLOR_CODING_PROJECTED));
                break;

            case TECHNIQUE_ARROW_PLOT_RAND:
            case TECHNIQUE_ARROW_PLOT_GRID:
                    sliceTexture_ = renderArrowPlotTexture(flow.extractSlice(permutation_, sliceNo),
                        outputSize, privatePort1_, thresholds);
                break;

            case TECHNIQUE_SPOTNOISE:
            case TECHNIQUE_SPOTNOISE_PROJECTED:
                {
                    std::vector<RenderPort*> tempPorts;
                    tempPorts.push_back(&privatePort1_);
                    tempPorts.push_back(&privatePort2_);

                    sliceTexture_ = renderSpotNoiseTexture(flow, static_cast<size_t>(sliceNo),
                        outputSize, viewportSize, tempPorts,
                        (techniqueProp_->getValue() == TECHNIQUE_SPOTNOISE_PROJECTED));
                }
                break;

            case TECHNIQUE_INTEGRATE_DRAW:
            case TECHNIQUE_FAST_LIC:
            case TECHNIQUE_INTEGRATE_DRAW_PROJECTED:
            case TECHNIQUE_FAST_LIC_PROJECTED:
                sliceTexture_ = renderFlowTexture(flow, sliceNo, static_cast<int>(tgt::round(scaling)),
                    techniqueProp_->getValue(), thresholds);
                break;
        }   // switch
    }   // if (rebuildTexture_)

    GLuint texID = 0;
    if (sliceTexture_ != 0) {
        sliceTexture_->uploadTexture();
        glBindTexture(GL_TEXTURE_2D, 0);
        rebuildTexture_ = false;
        texID = sliceTexture_->getId();
    }

    imgOutport_.activateTarget("FlowSliceRenderer2D::process()");

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.loadIdentity();
    glOrtho(0.0f, viewportSize.x, 0.0f, viewportSize.y, -1.0f, 1.0f);
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.loadIdentity();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glActiveTexture(GL_TEXTURE0);
    if (texID != 0) {
        // Bind the texture containing the flow image, if one has been generated
        //
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, texID);
        glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    }

    glPushAttrib(GL_TEXTURE_BIT);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    MatStack.translate(slicePos.x, slicePos.y, 0.0f);
    glBegin(GL_QUADS);
        glMultiTexCoord2f(GL_TEXTURE0, 0.0f, 0.0f);
        glVertex2f(0.0f, 0.0f);

        glMultiTexCoord2f(GL_TEXTURE0, 1.0f, 0.0f);
        glVertex2f(outputSize.x, 0.0f);

        glMultiTexCoord2f(GL_TEXTURE0, 1.0f, 1.0f);
        glVertex2f(outputSize.x, outputSize.y);

        glMultiTexCoord2f(GL_TEXTURE0, 0.0f, 1.0f);
        glVertex2f(0.0f, outputSize.y);
    glEnd();

    glPopAttrib();
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glDisable(GL_TEXTURE_2D);

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();

    imgOutport_.deactivateTarget();
}

// protected methods
//

std::string FlowSliceRenderer2D::generateShaderHeader() {
    std::ostringstream oss;
    oss << FlowSliceRenderer::generateShaderHeader();
    oss << "#define PLANE " << static_cast<size_t>(alignment_) << std::endl;
    return oss.str();
}

// private methods
//

float FlowSliceRenderer2D::getTextureScalingFactor(const tgt::vec2& viewportSize,
                                                   const tgt::vec2& textureSize)
{
    float scaleWidth = viewportSize.x / textureSize.x;
    float scaleHeight = viewportSize.y / textureSize.y;
    float scaling = (scaleWidth <= scaleHeight) ? scaleWidth : scaleHeight;
    if (zoom_ != ZOOM_FULL_SIZE) {
        // find minimal scaling factor
        //
        float zoom = static_cast<float>(zoom_);
        if (zoom <= scaling)
            scaling = zoom;
        else {
            zoom_ = ZOOM_FULL_SIZE;
            textureZoomProp_->set("");
            textureZoomProp_->updateWidgets();
        }
    }
    return scaling;
}

void FlowSliceRenderer2D::onSliceAlignmentChange() {
    alignment_ = alignmentProp_->getValue();
    permutation_ = getCoordinatePermutation(alignment_);
    rebuildTexture_ = true;
    updateNumSlices();
    rebuildShader();
}

void FlowSliceRenderer2D::onTextureZoomChange() {
    zoom_ = textureZoomProp_->getValue();
    rebuildTexture_ = true;
}

void FlowSliceRenderer2D::updateNumSlices() {
    numSlices_ = static_cast<int>(flowDimensions_[alignment_]);
    if (numSlices_ > sliceIndexProp_.getMinValue()) {
        sliceIndexProp_.setMaxValue(numSlices_);
        if (sliceIndexProp_.get() > numSlices_)
            sliceIndexProp_.set(numSlices_ / 2);
        sliceIndexProp_.updateWidgets();
    }
}

}   // namespace
