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

#include "flowslicerenderer.h"

#include "tgt/glmath.h"
#include "tgt/texture.h"

#include "../utils/flowmath.h"
#include "../datastructures/deprecated/streamlinetexture.h"

#include <sstream>
#include <string>
#include <vector>

namespace voreen {

FlowSliceRenderer::FlowSliceRenderer()
    : RenderProcessor(),
    techniqueProp_(0),
    numSeedingPositionsProp_("numSeedingPositions", "# seeding positions: ", 10, 1, 10000),
    spotTesselationProp_("spotTesselation", "spot tesselation: ", 16, 2, 100),
    spotScaleProp_("spotScale", "spot size (%): ", 5, 1, 100),
    useSpotFilteringProp_("useSpotFilteringProp", "filter spots: ", true),
    useSpotBendingProp_("useSpotBendingProp", "use spotBending: ", true),
    useAdaptiveKernelSizeProp_("useAdaptiveKernelSizeProp", "use adaptive kernel size: ", false),
    kernelSizeProp_("kernelSizeProp", "max. convolution kernel size: ", 7, 1, 100),
    pixelSamplingProp_("pixelSamplingProp", "pixel sampling: ", 50, 1, 10000),
    arrowSizeProp_("arrowSizeProp", "arrow size (pixel / zoom): ", 10, 1, 2000),
    arrowSpacingProp_("arrowSpacingProp", "arrow spacing (x arrow size): ", 2, 1, 100),
    zoom_(ZOOM_1X),
    permutation_(0, 1, 2),
    shader_(0),
    randomIntensities_(0),
    randomPositions_(0),
    noiseTexture_(0),
    flow2DTexture_(0),
    flow3DTexture_(0),
    rebuildTexture_(true),
    arrowList_(0),
    volInport_(Port::INPORT, "volumehandle.volumehandle", "Volume Input"),
    imgOutport_(Port::OUTPORT, "image.outport", "Image Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER),
    privatePort1_(Port::OUTPORT, "image.temp1", "image.temp1", false),
    privatePort2_(Port::OUTPORT, "image.temp2", "image.temp2", false),
    portGroup_(true)
{
    techniqueProp_ = new OptionProperty<RenderingTechnique>("visualizationMode", "technique:");
    techniqueProp_->addOption("color coding", "color coding",
        TECHNIQUE_COLOR_CODING);
    techniqueProp_->addOption("Integrate & Draw", "Integrate & Draw",
        TECHNIQUE_INTEGRATE_DRAW);
    techniqueProp_->addOption("fast LIC", "fast LIC",
        TECHNIQUE_FAST_LIC);
    techniqueProp_->addOption("arrow plot (random positions)", "arrow plot (random positions)",
        TECHNIQUE_ARROW_PLOT_RAND);
    techniqueProp_->addOption("arrow plot (regular grid)", "arrow plot (regular grid)",
        TECHNIQUE_ARROW_PLOT_GRID);
    techniqueProp_->addOption("spot noise", "spot noise",
        TECHNIQUE_SPOTNOISE);
    techniqueProp_->addOption("color coding (projected)", "color coding (projected)",
        TECHNIQUE_COLOR_CODING_PROJECTED);
    techniqueProp_->addOption("Integrate & Draw (projected)", "Integrate & Draw (projected)",
        TECHNIQUE_INTEGRATE_DRAW_PROJECTED);
    techniqueProp_->addOption("fast LIC (projected)", "fast LIC (projected)",
        TECHNIQUE_FAST_LIC_PROJECTED);
    techniqueProp_->addOption("spot noise (projected)", "spot noise (projected)",
        TECHNIQUE_SPOTNOISE_PROJECTED);
    techniqueProp_->onChange(
        MemberFunctionCallback<FlowSliceRenderer>(this, &FlowSliceRenderer::onTechniqueChange));

    MemberFunctionCallback<FlowSliceRenderer> invalidateAction(this, &FlowSliceRenderer::invalidateTexture);

    OptionProperty<ColorCodingAbility::ColorCodingMode>& colorCodingModeProp =
        colorCoding_.getColorCodingModeProp();
    colorCodingModeProp.onChange(
        MemberFunctionCallback<FlowSliceRenderer>(this, &FlowSliceRenderer::onColorCodingChange));

    IntOptionProperty& colorTableProp = colorCoding_.getColorTableProp();
    colorTableProp.onChange(
        MemberFunctionCallback<FlowSliceRenderer>(this, &FlowSliceRenderer::onColorCodingChange));

    FloatVec4Property& colorProp = colorCoding_.getColorProp();
    colorProp.onChange(MemberFunctionCallback<FlowSliceRenderer>(this, &FlowSliceRenderer::onColorCodingChange));

    thresholdProp_.onChange(invalidateAction);
    maxStreamlineLengthProp_.onChange(invalidateAction);

    numSeedingPositionsProp_.onChange(
        MemberFunctionCallback<FlowSliceRenderer>(this, &FlowSliceRenderer::onNumberSeedingPositionsChange));
    spotScaleProp_.onChange(invalidateAction);

    spotTesselationProp_.onChange(invalidateAction);
    useSpotFilteringProp_.onChange(invalidateAction);
    useSpotBendingProp_.onChange(invalidateAction);
    useAdaptiveKernelSizeProp_.onChange(invalidateAction);
    kernelSizeProp_.onChange(invalidateAction);
    pixelSamplingProp_.onChange(invalidateAction);
    arrowSizeProp_.onChange(invalidateAction);
    arrowSpacingProp_.onChange(invalidateAction);

    addProperty(techniqueProp_);
    addProperty(maxStreamlineLengthProp_);
    addProperty(thresholdProp_);
    addProperty(colorCodingModeProp);
    addProperty(colorTableProp);
    addProperty(colorProp);
    addProperty(numSeedingPositionsProp_);
    addProperty(spotTesselationProp_);
    addProperty(spotScaleProp_);
    addProperty(useSpotFilteringProp_);
    addProperty(useSpotBendingProp_);
    addProperty(useAdaptiveKernelSizeProp_);
    addProperty(kernelSizeProp_);
    addProperty(pixelSamplingProp_);
    addProperty(arrowSizeProp_);
    addProperty(arrowSpacingProp_);

    addPort(volInport_);
    addPort(imgOutport_);

    addPrivateRenderPort(&privatePort1_);
    addPrivateRenderPort(&privatePort2_);

    onTechniqueChange();
    toggleProperties();
    colorCoding_.onColorCodingModeChange();
}

FlowSliceRenderer::~FlowSliceRenderer() {
    delete techniqueProp_;

    delete [] randomIntensities_;
    delete [] randomPositions_;
    delete noiseTexture_;
    delete flow2DTexture_;
}

void FlowSliceRenderer::deinitialize() {
    if (shader_ != 0) {
        shader_->deactivate();
        ShdrMgr.dispose(shader_);
        shader_ = 0;
    }
    portGroup_.deinitialize();
    RenderProcessor::deinitialize();
}

void FlowSliceRenderer::initialize() {
    FlowreenProcessor::init();

    RenderProcessor::initialize();

    portGroup_.initialize();
    portGroup_.addPort(privatePort1_);
    portGroup_.addPort(privatePort2_);
}

// protected methods
//

tgt::ivec3 FlowSliceRenderer::getCoordinatePermutation(const SliceAlignment& alignment) {
    switch (alignment) {
        case PLANE_XY:
            return tgt::ivec3(0, 1, 2);
        case PLANE_XZ:
            return tgt::ivec3(0, 2, 1);
        case PLANE_ZY:
            return tgt::ivec3(2, 1, 0);
    }
    return tgt::ivec3(0, 0, 0);
}

void FlowSliceRenderer::beginOrthographicRendering(const tgt::ivec2& viewport) const {
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    // Set a viewport of the size of the generated texture: this is ugly but still better
    // than resizing the entire texture container twice. Shaders have to take the scaling
    // of the viewport into account when calculating texture coordinates (see Spot Noise)
    //
    glViewport(0, 0, viewport.x, viewport.y);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.loadIdentity();
    glOrtho(0.0, viewport.x, 0.0, viewport.y, -1.0, 1.0);
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.loadIdentity();

    glDisable(GL_DEPTH_TEST);
}

void FlowSliceRenderer::endOrthographicRendering() const {

    glEnable(GL_DEPTH_TEST);

    glPopAttrib();
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();
}

void FlowSliceRenderer::prepareFlowTexture2D(const Flow2D& flow2D) {
    delete flow2DTexture_;
    flow2DTexture_ = 0;
    tgt::vec3* pixels = flow2D.getNormalizedFlow();
    flow2DTexture_ = new tgt::Texture(tgt::ivec3(flow2D.dimensions_, 1), GL_RGB, GL_RGB16, GL_FLOAT, tgt::Texture::NEAREST,
                                      tgt::Texture::REPEAT, reinterpret_cast<GLubyte*>(pixels), true);
    flow2DTexture_->uploadTexture();
    glBindTexture(GL_TEXTURE_2D, 0);
}

void FlowSliceRenderer::prepareFlowTexture3D() {
    if ((currentVolume_ != 0) && (currentVolume_->getRepresentation<VolumeGL>() != 0))
        flow3DTexture_ = currentVolume_->getRepresentation<VolumeGL>()->getTexture();
    else
        flow3DTexture_ = 0;
}

tgt::Texture* FlowSliceRenderer::copyTargetToTexture(RenderPort& renderPort,
                                                     const tgt::vec2& textureSize) const
{
    renderPort.activateTarget();

    //glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    tgt::ivec2 size = static_cast<tgt::ivec2>(textureSize);
    tgt::Texture* texture = new tgt::Texture(tgt::vec3(textureSize, 1), GL_RGBA, GL_RGBA16,
        GL_FLOAT, tgt::Texture::NEAREST);
    texture->alloc(true);
    glReadPixels(0, 0, size.x, size.y, GL_RGBA, GL_FLOAT, texture->getCpuTextureData());

    //renderPort.deactivateTarget();
    return texture;
}

tgt::Texture* FlowSliceRenderer::renderArrowPlotTexture(const Flow2D& flow,
                                                        const tgt::vec2& textureSize,
                                                        RenderPort& tempPort,
                                                        const tgt::vec2& thresholds)
{
    prepareColorCoding2DShader(flow, thresholds);

    tempPort.activateTarget();
    beginOrthographicRendering(static_cast<tgt::ivec2>(textureSize));

    if ((rebuildTexture_ == false) && (arrowList_ != 0))
        glCallList(arrowList_);
    else {
        if (techniqueProp_->getValue() == TECHNIQUE_ARROW_PLOT_GRID)
            buildArrowGridDisplayList(flow, textureSize, thresholds);
        else
            buildArrowDisplayList(flow, textureSize, thresholds);
    }

    if ((shader_ != 0) && (shader_->isActivated() == true))
        shader_->deactivate();

    endOrthographicRendering();

    return copyTargetToTexture(tempPort, static_cast<tgt::ivec2>(textureSize));
}

tgt::Texture* FlowSliceRenderer::renderColorCodingTexture(const Flow3D& flow, const size_t sliceNo,
                                                          const tgt::vec2& textureSize,
                                                          RenderPort& tempPort,
                                                          const tgt::vec2& thresholds,
                                                          const bool projected)
{
    tempPort.resize(static_cast<tgt::ivec2>(textureSize));
    if (projected == true) {
        prepareColorCoding3DShader(flow, sliceNo, thresholds);
        glPushAttrib(GL_TEXTURE_BIT);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    } else
        prepareColorCoding2DShader(flow.extractSlice(permutation_, sliceNo), thresholds);

    tempPort.activateTarget();
    beginOrthographicRendering(static_cast<tgt::ivec2>(textureSize));

    glActiveTexture(GL_TEXTURE0);
    glBegin(GL_QUADS);
        glMultiTexCoord2f(GL_TEXTURE0, 0.0f, 0.0f);
        glVertex2f(0.0f, 0.0f);

        glMultiTexCoord2f(GL_TEXTURE0, 1.0f, 0.0f);
        glVertex2f(textureSize.x, 0.0f);

        glMultiTexCoord2f(GL_TEXTURE0, 1.0f, 1.0f);
        glVertex2fv(textureSize.elem);

        glMultiTexCoord2f(GL_TEXTURE0, 0.0f, 1.0f);
        glVertex2f(0.0f, textureSize.y);
    glEnd();

    if ((shader_ != 0) && (shader_->isActivated() == true))
        shader_->deactivate();

    endOrthographicRendering();

    if (projected == true)
        glPopAttrib();

    return copyTargetToTexture(tempPort, static_cast<tgt::ivec2>(textureSize));
}

tgt::Texture* FlowSliceRenderer::renderSpotNoiseTexture(const Flow3D& flow, const size_t sliceNo,
                                                        const tgt::vec2& textureSize,
                                                        const tgt::vec2& viewportSize,
                                                        const std::vector<RenderPort*>& tempPorts,
                                                        const bool projected)
{

    std::string vertexShader = (projected == true) ? "spotnoise3d.vert" : "spotnoise2d.vert";

    if ((loadShader(vertexShader, "spotnoise.frag") == false) || (rebuildTexture_ == false))
        return 0;

    // Scaling factor for texture coordinate calculation by the shader:
    // texture coordinates have to take the "viewport" into account as the
    // render target might be much larger than the generated texture.
    //
    const tgt::vec2 texScale = textureSize / viewportSize;
    const size_t numSpots = static_cast<size_t>(numSeedingPositionsProp_.get());
    const size_t spotTesselation = static_cast<size_t>(spotTesselationProp_.get());

    shader_->activate();
    if (projected == true) {
        prepareFlowTexture3D();

        shader_->setUniform("flowDimensions_", static_cast<tgt::vec3>(flow.dimensions_));
        shader_->setUniform("sliceNo_", static_cast<GLint>(sliceNo));
        shader_->setUniform("vectorField3D_", 0);
    } else {
        Flow2D flow2D = flow.extractSlice(permutation_, sliceNo);
        prepareFlowTexture2D(flow2D);

        shader_->setUniform("flowDimensions_", static_cast<tgt::vec2>(flow2D.dimensions_));
        shader_->setUniform("vectorField2D_", 0);
    }

    shader_->setUniform("textureDimensions_", textureSize);
    shader_->setUniform("textureScaling_", texScale);
    shader_->setUniform("useFilter_", useSpotFilteringProp_.get());
    shader_->setUniform("spotScale_", (spotScaleProp_.get() / 100.0f));
    shader_->setUniform("minValue_", flow.minValue_);
    shader_->setUniform("maxValue_", flow.maxValue_);
    shader_->setUniform("maxMagnitude_", flow.maxMagnitude_);
    shader_->setUniform("spotColor_", tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f));
    shader_->setUniform("spotNoiseTex_", 1);
    shader_->setUniform("useSpotBending_", useSpotBendingProp_.get());

    initRandomValues(numSpots);

    portGroup_.reattachTargets();   // This is necessary to attach any texture to the FBO at all.
    portGroup_.activateTargets();   // This ensures the activation and the binding of the FBO.

    glClearColor(1.0f, 0.5f, 0.0f, 1.0f);
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

    beginOrthographicRendering(static_cast<tgt::ivec2>(textureSize));

    glActiveTexture(GL_TEXTURE0);
    if ((projected == true) && (flow3DTexture_ != 0))
        glBindTexture(GL_TEXTURE_3D, flow3DTexture_->getId());
    else if (flow2DTexture_ != 0)
        glBindTexture(GL_TEXTURE_2D, flow2DTexture_->getId());

    glActiveTexture(GL_TEXTURE1);
    for (unsigned int i = 0; i < numSpots; ++i) {
        // bind alternatingly one of the two texture as the input texture from which the
        // shader can read.
        //
        tempPorts[(i % 2)]->bindColorTexture();
        shader_->setUniform("randPos_", randomPositions_[i]);
        shader_->setUniform("randIntensity_", randomIntensities_[i]);

        renderSpotMesh(spotTesselation);
    }

    if ((shader_ != 0) && (shader_->isActivated() == true))
        shader_->deactivate();

    endOrthographicRendering();

    return copyTargetToTexture(*(tempPorts[numSpots % 2]), textureSize);
}

tgt::Texture* FlowSliceRenderer::renderFlowTexture(const Flow3D& flow, const size_t sliceNo,
                                                   const int textureScaling,
                                                   const RenderingTechnique& technique,
                                                   const tgt::vec2& thresholds)
{
    const size_t sampling = static_cast<size_t>(pixelSamplingProp_.get());

    // get the size of the slice and calculate the output texture's size, based
    // on the nearest integral scaling factor.
    //
    tgt::ivec2 sliceSize = flow.getFlowSliceDimensions(permutation_);

    float* pixels = 0;
    switch (technique) {
        case TECHNIQUE_INTEGRATE_DRAW:
            pixels = StreamlineTexture<float>::integrateDraw(flow.extractSlice(permutation_, sliceNo),
                textureScaling, sampling, thresholds);
            break;

        case TECHNIQUE_FAST_LIC:
            {
                Flow2D flow2D = flow.extractSlice(permutation_, sliceNo);
                const int kernelSize = kernelSizeProp_.get();
                SimpleTexture<float> noiseTexture(flow2D.dimensions_, true);
                noiseTexture.createWhiteNoise();
                pixels = StreamlineTexture<float>::fastLIC(flow2D, noiseTexture,
                    textureScaling, sampling, kernelSize, thresholds, useAdaptiveKernelSizeProp_.get());
            }
            break;

        case TECHNIQUE_INTEGRATE_DRAW_PROJECTED:
            pixels = integrateDraw(flow, sliceNo, textureScaling);
            break;

        case TECHNIQUE_FAST_LIC_PROJECTED:
            pixels = fastLIC(flow, sliceNo, textureScaling);
            break;

        default:
            LERROR("renderFlowSlice(): unsupport rendering technique #'"
                << static_cast<int>(technique) << "'!");
            return 0;
    }

    if (pixels == 0)
        return 0;

    return new tgt::Texture(tgt::ivec3(sliceSize * textureScaling, 1), GL_LUMINANCE, GL_LUMINANCE, GL_FLOAT,
                            tgt::Texture::NEAREST, tgt::Texture::REPEAT, reinterpret_cast<GLubyte*>(pixels), true);
}

void FlowSliceRenderer::initNoiseTexture(const tgt::ivec2& size) {
    if ((rebuildTexture_ == true) || (noiseTexture_ == 0)) {
        delete noiseTexture_;
        noiseTexture_ = new SimpleTexture<unsigned char>(size, true);
        noiseTexture_->createWhiteNoise();
    }
}

void FlowSliceRenderer::initRandomValues(const size_t count) {
    if (randomIntensities_ == 0) {
        randomIntensities_ = new float[count];
        for (size_t i = 0; i < count; ++i)
            randomIntensities_[i] = 0.2f * (FlowMath::uniformRandom() - 0.5f); // [-0.1, 0.1]
    }

    if (randomPositions_ == 0) {
        randomPositions_ = new tgt::vec2[count];
        for (size_t i = 0; i < count; ++i)
            randomPositions_[i] = FlowMath::uniformRandomVec2();
    }
}

tgt::vec2 FlowSliceRenderer::reseedRandomPosition(const size_t validPositions) {
    tgt::vec2 randVec = FlowMath::uniformRandomVec2();

    // Use a "die" to determine wether a completely new random position
    // will be taken or wether an exisiting one will be used.
    // When the probability for a new position is low, the seeding positions
    // seem be prone to cluster at single location.
    //
    int die = rand() % 6;
    if ((die < 3) || (validPositions <= 1) || (randomPositions_ == 0))
        return randVec;

    // if there are already random positions which lead to a vector-field value not being
    // zero or which falls within the limits defined by thresholds, take this
    // position to generate another.
    // Use the position and add some random offset to it.
    //
    size_t index = rand() % validPositions;
    randVec *= ((rand() % 10) + 1) / 100.0f;

    return tgt::clamp((randomPositions_[index] + randVec), tgt::vec2::zero, tgt::vec2(1.0f));
}

std::string FlowSliceRenderer::generateShaderHeader() {
    std::string header = RenderProcessor::generateHeader();
    header += colorCoding_.getShaderDefines();
    return header;
}

bool FlowSliceRenderer::loadShader(const std::string& vertexShaderName,
                                   const std::string& fragmentShaderName)
{
    if (shader_ == 0) {
        try {
            shader_ = ShdrMgr.loadSeparate(vertexShaderName.c_str(), fragmentShaderName.c_str(),
                generateShaderHeader(), false);
        } catch(tgt::Exception) {
            LERROR("Failed to load shaders '" <<  vertexShaderName << "' & '"
                << fragmentShaderName << "'!");
            return false;
        }
        rebuildShader();
    }
    return true;
}

bool FlowSliceRenderer::rebuildShader() {
    if (shader_ == 0)
        return false;

    shader_->deactivate();
    shader_->setHeaders(generateShaderHeader());
    return shader_->rebuild();
}

// private methods
//

void FlowSliceRenderer::buildArrowDisplayList(const Flow2D& flow, const tgt::vec2& textureSize,
                                              const tgt::vec2& thresholds)
{
    const size_t numStreamlines = static_cast<size_t>(numSeedingPositionsProp_.get());

    // integration length
    //
    const float length = maxStreamlineLengthProp_.get();

    // stepwidth for integrator
    //
    const float stepwidth = 0.5f;

    // spacing between the arrows in multiples of arrow length (in pixels)
    //
    const size_t spacing = static_cast<size_t>(tgt::round(arrowSpacingProp_.get() / stepwidth));

    // size of an arrow (in pixels / zoomfactor)
    //
    const int arrowSize = arrowSizeProp_.get();

    const tgt::vec2 dim = static_cast<tgt::vec2>(flow.dimensions_);

    // scaling factor for a single arrow relative to texture size
    //
    const float scaling = tgt::min(textureSize / dim);

    initRandomValues(numStreamlines);

    if (arrowList_ != 0)
        glDeleteLists(arrowList_, 1);
    arrowList_ = glGenLists(1);

    size_t i = 0;
    size_t numTries = 0;    // number of total tries to get a valid seeding position
    const size_t maxTries = numStreamlines * 5; // maxiaml number of tries per streamline
    glNewList(arrowList_, GL_COMPILE);
    for (; ((i < numStreamlines) && (numTries < maxTries)); ++numTries) {
        const tgt::vec2 r0(randomPositions_[i] * dim);

        std::deque<tgt::vec2> streamlineD =
            FlowMath::computeStreamlineRungeKutta(flow, r0, length, stepwidth, 0, thresholds);

        // If the streamline contains only one element, the seeding position lead
        // to a result being a null-vector (zero velocity) or a result which did
        // not fall into the range defined by the thresholds.
        // Therefore, the last seeding position is taken and a new one at the
        // end of the array of random positions is generated.
        //
        if (streamlineD.size() <= 1) {
            randomPositions_[i] = reseedRandomPosition(i);
            continue;
        }
        std::vector<tgt::vec2> streamline = FlowMath::dequeToVector(streamlineD);
        const size_t step = spacing * arrowSize;
        for (size_t n = 0; n < streamline.size(); n += step) {
            const tgt::vec2& r = streamline[n];

            tgt::vec2 texCoord = r / dim;
            glMultiTexCoord2fv(GL_TEXTURE0, texCoord.elem);

            tgt::mat3 trafo = flow.getTransformationMatrix(r, static_cast<float>(arrowSize));
            renderArrow(trafo * scaling);
        }
        ++i;
    }   // for (i
    glEndList();

    if (numTries >= maxTries) {
        LINFO("Only " << i << " streamlines could be created from valid random seeding positions. \
Giving up after " << numTries << " tries.\n");
    }
    if (arrowList_ != 0)
        glCallList(arrowList_);
}


void FlowSliceRenderer::buildArrowGridDisplayList(const Flow2D& flow, const tgt::vec2& textureSize,
                                              const tgt::vec2& thresholds)
{
    // size of an arrow (in pixels / zoomfactor)
    //
    const float arrowSize = static_cast<float>(arrowSizeProp_.get());

    // scaling factor for a single arrow relative to texture size
    //
    const tgt::vec2 dim = static_cast<tgt::vec2>(flow.dimensions_);
    const float scaling = tgt::min(textureSize / dim);

    // spacing between the arrows in multiples of arrow length (in pixels)
    //
    int spacing = arrowSpacingProp_.get() * arrowSizeProp_.get();
    if (spacing <= 0)
        spacing = 1;
    else if (spacing > tgt::min(flow.dimensions_))
        spacing = tgt::min(flow.dimensions_);

    tgt::ivec2 grid = flow.dimensions_ / spacing;

    if (arrowList_ != 0)
        glDeleteLists(arrowList_, 1);
    arrowList_ = glGenLists(1);

    glNewList(arrowList_, GL_COMPILE);
    for (int y = 0; y <= grid.y; ++y) {
        for (int x = 0; x <= grid.x; ++x) {
            tgt::ivec2 r(x * spacing, y * spacing);
            if (flow.isInsideBoundings(r) == false)
                continue;

            tgt::vec2 v = flow.lookupFlow(r);
            if (v == tgt::vec2::zero)
                continue;

            if (thresholds != tgt::vec2::zero) {
                float magnitude = tgt::length(v);
                if ((magnitude < thresholds.x) || (magnitude > thresholds.y))
                    continue;
            }

            tgt::vec2 texCoord(r.x / dim.x, r.y / dim.y);
            glMultiTexCoord2fv(GL_TEXTURE0, texCoord.elem);

            tgt::mat3 trafo = flow.getTransformationMatrix(r, arrowSize);
            renderArrow(trafo * scaling);
        }   // for (x
    }   // for (y
    glEndList();

    if (arrowList_ != 0)
        glCallList(arrowList_);
}

float* FlowSliceRenderer::fastLIC(const Flow3D& flow, const size_t sliceNo, const int textureScaling)
{
    const tgt::ivec2 inputTexSize = flow.getFlowSliceDimensions(permutation_);
    initNoiseTexture(inputTexSize);

    const tgt::ivec2 outputTexSize(inputTexSize * textureScaling);

    // initialize output texture and counter for hits per pixel on output texture
    //
    const size_t numOutputPixels = outputTexSize.x * outputTexSize.y;
    int* numHits = new int[numOutputPixels];
    float* output = new float[numOutputPixels];
    memset(numHits, 0, numOutputPixels * sizeof(int));
    memset(output, 0, numOutputPixels * sizeof(float));

    const float stepSize = 0.5f / textureScaling;   // stepwidth for streamline integration
    const int kernelSize = kernelSizeProp_.get();
    const float k = 1.0f / (((2.0f * kernelSize) + 1.0f) * textureScaling);
    const int delta = pixelSamplingProp_.get(); // 1/rate of pixels to be sampled on output texture

    size_t numStreamlines = 0;
    for (int y = 0; y < outputTexSize.y; y += delta) {
        for (int x = 0; x < outputTexSize.x; x += delta) {
            int index = y * outputTexSize.x + x;
            if (numHits[index] > 0)
                continue;

            // get the coordinates of the pixel in the input texture which corresponds
            // to this position in the output texture and calculate its position within
            // the flow.
            //
            tgt::ivec2 r0Input = tgt::ivec2(x, y) / textureScaling;
            tgt::ivec2 errorInput(0, 0);
            tgt::vec3 r0 =
                flow.slicePosToFlowPos(r0Input, inputTexSize, permutation_, sliceNo, &errorInput);

            if (flow.lookupFlow(r0) == tgt::vec3::zero)
                continue;

            // also determine the round-off error which occurs if the flow positions was
            // converted back directly to the coordinates of the output textures.
            //
            tgt::ivec2 errorOutput(0, 0);
            flow.slicePosToFlowPos(tgt::ivec2(x, y), outputTexSize, permutation_, sliceNo, &errorOutput);

            // start streamline computation
            //
            ++numStreamlines;
            int indexR0 = 0;
            std::deque<tgt::vec3> deque = FlowMath::computeStreamlineRungeKutta(flow, r0, 150.0f, stepSize, &indexR0);
            std::vector<tgt::vec3> streamline = FlowMath::dequeToVector(deque);

            // copy the streamline for second coordinate conversion
            //
            std::vector<tgt::vec3> streamlineCopy(streamline);

            // convert the streamline into dimensions of the input texture
            //
            std::vector<tgt::ivec2> streamlineInput =
                flow.flowPosToSlicePos(streamline, inputTexSize, permutation_, errorInput);

            // also convert the streamline into dimensions of the output texture
            //
            std::vector<tgt::ivec2> streamlineOutput =
                flow.flowPosToSlicePos(streamlineCopy, outputTexSize, permutation_, errorOutput);

            // calculate initial intensity for the starting pixel
            //
            float intensity0 = k * initialFastLICIntensity(indexR0, kernelSize, streamlineInput);

            // determine the affected pixel in the output texture and add the
            // initial intensity
            //
            tgt::ivec2& outputTexCoord = streamlineOutput[indexR0];
            size_t pixel = outputTexCoord.y * outputTexSize.x + outputTexCoord.x;
            output[pixel] += intensity0;
            ++numHits[pixel];

            // trace streamline in forward direction and update intensity
            //
            float intensity = intensity0;
            int left = indexR0 + kernelSize + 1;
            int right = indexR0 - kernelSize;
            const int numPoints = static_cast<int>(streamlineInput.size());

            for (int i = (indexR0 + 1); i < numPoints; ++i, ++left, ++right) {
                int l = (left >= numPoints) ? (numPoints - 1) : left;
                const tgt::ivec2& a = streamlineInput[l];

                int r = (right <= 0) ? 0 : right;
                const tgt::ivec2& b = streamlineInput[r];

                intensity += (((*noiseTexture_)[a] / 255.0f) - ((*noiseTexture_)[b] / 255.0f)) * k;

                outputTexCoord = streamlineOutput[i];
                pixel = outputTexCoord.y * outputTexSize.x + outputTexCoord.x;
                ++numHits[pixel];
                output[pixel] += intensity;
            }

            // trace streamline in backward direction and update intensity
            //
            intensity = intensity0;
            left = indexR0 - kernelSize - 1;
            right = indexR0 + kernelSize;
            for (int i = (indexR0 - 1); i >= 0; --i, --left, --right) {
                int l = (left <= 0) ? 0 : left;
                const tgt::ivec2& a = streamlineInput[l];

                int r = (right >= numPoints) ? (numPoints - 1) : right;
                const tgt::ivec2& b = streamlineInput[r];

                intensity += (((*noiseTexture_)[a] / 255.0f) - ((*noiseTexture_)[b] / 255.0f)) * k;

                outputTexCoord = streamlineOutput[i];
                pixel = outputTexCoord.y * outputTexSize.x + outputTexCoord.x;
                ++numHits[pixel];
                output[pixel] += intensity;
            }
        }
    }

    size_t unhitPixels = 0;
    for (size_t i = 0; i < numOutputPixels; ++i) {
        if (numHits[i] > 1)
            output[i] /= numHits[i];
        else
            if (numHits[i] <= 0)
                ++unhitPixels;
    }
    std::cout << "# streamlines = " << numStreamlines << ", # unhit pixels = " << unhitPixels
        << " (" << static_cast<float>(100 * unhitPixels) / static_cast<float>(numOutputPixels) << " %)\n";
    delete [] numHits;
    return output;
}

float FlowSliceRenderer::initialFastLICIntensity(const int indexR0, const size_t kernelSize,
                                                   const std::vector<tgt::ivec2>& streamline) const
{
    const tgt::ivec2& texCoord = streamline[indexR0];
    int intensity = (*noiseTexture_)[texCoord];

    for (int n = 1; n <= static_cast<int>(kernelSize); ++n) {
        if (static_cast<size_t>(indexR0 + n) < streamline.size()) {
            const tgt::ivec2& texCoord = streamline[indexR0 + n];
            intensity += (*noiseTexture_)[texCoord];
        }
        if ((indexR0 - n) >= 0) {
            const tgt::ivec2& texCoord = streamline[indexR0 - n];
            intensity += (*noiseTexture_)[texCoord];
        }
    }
    return (intensity / 255.0f);
}

float* FlowSliceRenderer::integrateDraw(const Flow3D& flow, const size_t sliceNo, const int textureScaling)
{
    const tgt::ivec2 sliceSize = flow.getFlowSliceDimensions(permutation_);
    const tgt::ivec2 outputTexSize(sliceSize * textureScaling);

    const int texSize = outputTexSize.x * outputTexSize.y;
    int* numHits = new int[texSize];
    float* output = new float[texSize];
    memset(numHits, 0, texSize * sizeof(int));
    memset(output, 0, texSize * sizeof(float));

    const int delta = pixelSamplingProp_.get();
    size_t numStreamlines = 0;
    const float stepSize = (0.5f / textureScaling);

    for (int y = 0; y < outputTexSize.y; y += delta) {
        for (int x = 0; x < outputTexSize.x; x += delta) {
            int index = y * outputTexSize.x + x;
            if (numHits[index] > 0)
                continue;

            // By transforming r0 = [x, y] from viewport coordinates to voxel coordinates
            // for streamline calculation round-off errors will occur, so that the
            // transformation of the points on the streamline from voxel coordinates
            // backwards to viewport coordinates will result in a different r0' = [x', y']
            // for the original coordinate r0.
            // This happens, for example, if the viewport is of different size than the
            // slice through the flow so that one voxel corresponds to n x n pixels.
            // But as it is necessary to have exactly the same pixel hit by the computed
            // streamline, the error between r0 and r0' error = (r0 - r0') is calculated
            // and added to each transfomed point on the streamline. This will ensure,
            // that during the next iteration, if pixel r0 has already been hit, it will
            // not result in the re-computation of the streamline.
            //
            tgt::ivec2 error(0, 0);
            tgt::vec3 r0 =
                flow.slicePosToFlowPos(tgt::ivec2(x, y), outputTexSize, permutation_, sliceNo, &error);

            if (flow.lookupFlow(r0) == tgt::vec3::zero)
                continue;

            ++numStreamlines;
            float gray = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
            std::deque<tgt::vec3> streamline =
                FlowMath::computeStreamlineRungeKutta(flow, r0, 150.0f, stepSize);
            std::deque<tgt::ivec2> streamlineVP =
                flow.flowPosToSlicePos(streamline, outputTexSize, permutation_, error);

            for (std::deque<tgt::ivec2>::const_iterator it = streamlineVP.begin();
                it != streamlineVP.end(); ++it)
            {
                const tgt::ivec2& p = *it;
                size_t pixel = (p.y * outputTexSize.x) + p.x;
                output[pixel] += (++numHits[pixel] * gray);
            }   // for (it
        }   // for (x
    }   // for (y

    size_t unhitPixels = 0;
    for (int i = 0; i < texSize; ++i) {
        if (numHits[i] > 1)
            output[i] /= (numHits[i] * (numHits[i] + 1) / 2.0f);
        else
            if (numHits[i] <= 0)
                ++unhitPixels;
    }
    std::cout << "# streamlines = " << numStreamlines << ", #unhit pixels = " << unhitPixels
        << " (" << static_cast<float>(100 * unhitPixels) / static_cast<float>(texSize) << " %)\n";
    delete [] numHits;
    return output;
}

void FlowSliceRenderer::invalidateTexture() {
    rebuildTexture_ = true;
}

void FlowSliceRenderer::onColorCodingChange() {
    rebuildTexture_ = true;
    rebuildShader();
}

void FlowSliceRenderer::onNumberSeedingPositionsChange() {
    delete [] randomIntensities_;
    randomIntensities_ = 0;
    delete [] randomPositions_;
    randomPositions_ = 0;
    rebuildTexture_ = true;
}

void FlowSliceRenderer::onTechniqueChange() {
    // Dispose the shader as another shader or no shader
    // needs to be used for the selected technique.
    //
    if (shader_ != 0) {
        shader_->deactivate();
        ShdrMgr.dispose(shader_);
        shader_ = 0;
    }

    // Also clear all kinds of random values. If a technique
    // has been selected requiring one them, it will re-
    // initialize them.
    //
    delete [] randomIntensities_;
    randomIntensities_ = 0;
    delete [] randomPositions_;
    randomPositions_ = 0;
    delete noiseTexture_;
    noiseTexture_ = 0;

    toggleProperties();

    // indicate that the texture has to be rebuild (used by Integrate & Draw
    // and FastLIC)
    //
    rebuildTexture_ = true;
}

void FlowSliceRenderer::prepareColorCoding2DShader(const Flow2D& flow2D,
                                                      const tgt::vec2& thresholds)
{
    if ((rebuildTexture_ == true) || (flow2DTexture_ == 0))
        prepareFlowTexture2D(flow2D);

    if ((loadShader("colorcoding.vert", "colorcoding2d.frag") == true)
        && (flow2DTexture_ != 0))
    {
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, flow2DTexture_->getId());

        shader_->activate();
        if (colorCoding_.getColorCodingMode() != ColorCodingAbility::COLOR_CODING_MODE_MONOCHROME) {
            shader_->setUniform("vectorField2D_", 1);
            shader_->setUniform("minValue_", flow2D.minValue_);
            shader_->setUniform("maxValue_", flow2D.maxValue_);
            if (colorCoding_.getColorCodingMode() != ColorCodingAbility::COLOR_CODING_MODE_DIRECTION) {
                shader_->setUniform("maxMagnitude_", flow2D.maxMagnitude_);
                shader_->setUniform("thresholds_", thresholds);
            }
        } else {
            const tgt::vec4& color = colorCoding_.getColorProp().get();
            shader_->setUniform("color_", color);
        }
    }
}

void FlowSliceRenderer::prepareColorCoding3DShader(const Flow3D& flow3D, const size_t sliceNo,
                                                   const tgt::vec2& thresholds)
{
    if ((rebuildTexture_ == true) || (flow3DTexture_ == 0))
        prepareFlowTexture3D();

    if ((loadShader("colorcoding.vert", "colorcoding3d.frag") == true)
        && (flow3DTexture_ != 0))
    {
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_3D, flow3DTexture_->getId());

        shader_->activate();
        shader_->setUniform("vectorField3D_", 1);
        shader_->setUniform("flowDimensions_", static_cast<tgt::vec3>(flow3D.dimensions_));
        shader_->setUniform("sliceNo_", static_cast<GLint>(sliceNo));
        shader_->setUniform("minValue_", flow3D.minValue_);
        shader_->setUniform("maxValue_", flow3D.maxValue_);
        if (colorCoding_.getColorCodingMode() != ColorCodingAbility::COLOR_CODING_MODE_DIRECTION) {
            shader_->setUniform("maxMagnitude_", flow3D.maxMagnitude_);
            shader_->setUniform("thresholds_", thresholds);
            shader_->setUniform("useBilinearInterpolation_", false);
        }
    }
}

void FlowSliceRenderer::renderArrow(const tgt::mat3& transformation) const {
    // basic geometry for the arrow (triangle + quad with tip pointing along
    // y-axis and bounding sizes of [0.5, 1.0], centered at [0.0, 0.0])
    //
    tgt::vec3 tri[7] = {
        tgt::vec3(-0.25f, 0.0f, 1.0f),  // left on base of triangle
        tgt::vec3(0.25f, 0.0f, 1.0f),   // right on base of triangle
        tgt::vec3(0.0f, 0.5f, 1.0f),    // top of triangle

        tgt::vec3(-0.1f, -0.5f, 1.0f),  // lower left corner of quad
        tgt::vec3(0.1f, -0.5f, 1.0f),   // lower right corner of quad
        tgt::vec3(0.1f, 0.0f, 1.0f),    // upper right corner of quad
        tgt::vec3(-0.1f, 0.0f, 1.0f)    // upper left corner of quad
    };

    // apply transformation to every vertex
    //
    for (size_t i = 0; i < 7; ++i)
        tri[i] = transformation * tri[i];

    glBegin(GL_TRIANGLES);
        glVertex2fv(tri[0]._xy.elem);
        glVertex2fv(tri[1]._xy.elem);
        glVertex2fv(tri[2]._xy.elem);
    glEnd();

    glBegin(GL_QUADS);
        glVertex2fv(tri[3]._xy.elem);
        glVertex2fv(tri[4]._xy.elem);
        glVertex2fv(tri[5]._xy.elem);
        glVertex2fv(tri[6]._xy.elem);
    glEnd();
}

void FlowSliceRenderer::renderSpotMesh(const size_t tesselation) {
    if (tesselation <= 1) {
        glBegin(GL_QUADS);
            glVertex2f(0.0f, 0.0f);
            glVertex2f(1.0f, 0.0f);
            glVertex2f(1.0f, 1.0f);
            glVertex2f(0.0f, 1.0f);
        glEnd();
        return;
    }

    const float f = 1.0f / static_cast<float>(tesselation);
    for (size_t iy = 0; iy < tesselation; ++iy) {
        glBegin(GL_QUAD_STRIP);
        for (size_t ix = 0; ix <= tesselation; ++ix) {
            glVertex2f(ix * f, (iy + 1) * f);
            glVertex2f(ix * f, iy * f);
        }
        glEnd();
    }
}

void FlowSliceRenderer::toggleProperties() {
    const size_t numProps = 15;
    OptionProperty<ColorCodingAbility::ColorCodingMode>& colorCodingModeProp =
        colorCoding_.getColorCodingModeProp();
    IntOptionProperty& colorTableProp = colorCoding_.getColorTableProp();
    FloatVec4Property& colorProp = colorCoding_.getColorProp();
    Property* arrProps[numProps] = {
        &colorCodingModeProp,       // 0
        &colorTableProp,            // 1
        &numSeedingPositionsProp_,  // 2
        &spotTesselationProp_,      // 3
        &spotScaleProp_,            // 4
        &useSpotFilteringProp_,     // 5
        &useSpotBendingProp_,       // 6
        &useAdaptiveKernelSizeProp_, // 7
        &kernelSizeProp_,           // 8
        &pixelSamplingProp_,        // 9
        &arrowSizeProp_,            // 10
        &arrowSpacingProp_,         // 11
        &maxStreamlineLengthProp_,  // 12
        &thresholdProp_,            // 13
        &colorProp                  // 14
    };

    bool isVisible[numProps] = {false};
    switch (techniqueProp_->getValue()) {
        case TECHNIQUE_COLOR_CODING:
        case TECHNIQUE_COLOR_CODING_PROJECTED:
            isVisible[0] = true;
            isVisible[1] = true;
            isVisible[14] = true;
            isVisible[13] = true;
            break;
        case TECHNIQUE_SPOTNOISE:
        case TECHNIQUE_SPOTNOISE_PROJECTED:
            isVisible[2] = true;
            isVisible[3] = true;
            isVisible[4] = true;
            isVisible[5] = true;
            isVisible[6] = true;
            break;
        case TECHNIQUE_FAST_LIC:
        case TECHNIQUE_FAST_LIC_PROJECTED:
            isVisible[7] = true;
            isVisible[8] = true;    // no break here!
        case TECHNIQUE_INTEGRATE_DRAW:
        case TECHNIQUE_INTEGRATE_DRAW_PROJECTED:
            isVisible[9] = true;
            isVisible[12] = true;
            isVisible[13] = true;
            break;
        case TECHNIQUE_ARROW_PLOT_RAND:
            isVisible[2] = true;
            isVisible[12] = true;
            // no break here!
        case TECHNIQUE_ARROW_PLOT_GRID:
            isVisible[0] = true;
            isVisible[1] = true;
            isVisible[14] = true;
            isVisible[10] = true;
            isVisible[11] = true;
            isVisible[13] = true;
            break;
        default:
            break;
    }

    for (size_t i = 0; i < numProps; ++i)
        arrProps[i]->setVisibleFlag(isVisible[i]);
}

}   // namespace
