/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "cpuraycaster.h"

#include "voreen/core/datastructures/transfunc/1d/1dkeys/transfunc1dkeys.h"
#include "voreen/core/datastructures/transfunc/2d/2dprimitives/transfunc2dprimitives.h"

#ifdef VRN_MODULE_OPENMP
#include "omp.h"
#endif

namespace voreen {

using tgt::vec2;
using tgt::vec3;
using tgt::vec4;
using tgt::svec3;

const float SAMPLING_BASE_INTERVAL_RCP = 200.0;

CPURaycaster::CPURaycaster()
  : VolumeRaycaster()
  , gradientVolumePort_(Port::INPORT, "volumehandle.gradientvolumehandle", "Gradient Volume Input")
  , tFuncType_("tFuncType", "TF Type")
  , tFunc1DKeys_("tFunc1DKeys", "1D Color Map")
  , tFunc2DPrimitives_("tFunc2DPrimitives", "2D Color Map")
  , texFilterMode_("textureFilterMode_", "Texture Filtering")
  , preIntegrationTableSize_("preIntegrationTableSize", "Width of pre-integration table")
{
    addPort(gradientVolumePort_);

    addProperty(tFuncType_);
    addProperty(tFunc1DKeys_);
        tFuncType_.addTransFuncProperty(&tFunc1DKeys_);
    addProperty(tFunc2DPrimitives_);
        tFuncType_.addTransFuncProperty(&tFunc2DPrimitives_);

    // volume texture filtering
    texFilterMode_.addOption("nearest", "Nearest",  GL_NEAREST);
    texFilterMode_.addOption("linear",  "Linear",   GL_LINEAR);
    texFilterMode_.selectByKey("linear");
    addProperty(texFilterMode_);

    addProperty(classificationMode_);

    preIntegrationTableSize_.addOption("deriveFromBitDepth", "Derive from bit depth", 0);
    preIntegrationTableSize_.addOption("256", "256", 256);
    preIntegrationTableSize_.addOption("512", "512", 512);
    preIntegrationTableSize_.addOption("1024", "1024", 1024);
    preIntegrationTableSize_.select("deriveFromBitDepth");
    addProperty(preIntegrationTableSize_);

}

Processor* CPURaycaster::create() const {
    return new CPURaycaster();
}

bool CPURaycaster::isReady() const {
    return (volumeInport_.hasData() && entryPort_.isConnected() && exitPort_.isConnected());
}

void CPURaycaster::initialize() {
    VolumeRaycaster::initialize();
}

void CPURaycaster::process() {

    tgtAssert(volumeInport_.getData()->getRepresentation<VolumeRAM>(), "no input volume");

    tFunc1DKeys_.setVolume(volumeInport_.getData());
    tFunc2DPrimitives_.setVolume(volumeInport_.getData());
    LGL_ERROR;

    // retrieve intensity volume
    const VolumeRAM* volume = volumeInport_.getData()->getRepresentation<VolumeRAM>();
    tgtAssert(volume, "no input volume");
    tgt::svec3 volDim = volume->getDimensions();

    // use dimension with the highest resolution for calculating the sampling step size
    float samplingStepSize = 1.f / (tgt::max(volDim) * samplingRate_.get());

    // retrieve tf texture
    tgt::Texture* tfTexture = 0;
    const PreIntegrationTable* table = 0;

    CpuRaycasterClassificationMode mode = TF;

    switch(tFuncType_.getValue()) {
    case TFT_1DKEYS:
    {
        if (classificationMode_.getKey() == "pre-integrated-cpu") {
            table = tFunc1DKeys_.get()->getPreIntegrationTable(samplingStepSize, static_cast<size_t>(preIntegrationTableSize_.getValue()), false);
            mode = PREINTEGRATED;
        }
        else if(classificationMode_.getKey() == "pre-integrated-approximation") {
            table = tFunc1DKeys_.get()->getPreIntegrationTable(samplingStepSize, static_cast<size_t>(preIntegrationTableSize_.getValue()), true);
            mode = PREINTEGRATED;
        }
        else if(classificationMode_.getKey() == "pre-integrated-gpu") {
            //table = tfi->getPreIntegrationTable(samplingStepSize, static_cast<size_t>(preIntegrationTableSize_.getValue()), false);
            LERROR("CPU raycaster does not support computing the pre-integration table on the GPU. Falling back to CPU.");
            classificationMode_.select("pre-integrated-cpu");
            invalidate();
            return;
        }
        tfTexture = tFunc1DKeys_.get()->getTexture();
        break;
    }
    case TFT_2DPRIMITIVES:
    {
        if (!gradientVolumePort_.hasData() || gradientVolumePort_.getData()->getRepresentation<VolumeRAM>()->getNumChannels() < 3) {
            LERROR("To use 2D tfs a RGB or RGBA gradient volume is needed");
            return;
        }
        if (gradientVolumePort_.getData()->getRepresentation<VolumeRAM>()->getDimensions() != volumeInport_.getData()->getRepresentation<VolumeRAM>()->getDimensions()) {
            LERROR("Gradient volume dimensions differ from intensity volume dimensions");
            return;
        }
        if (startsWith(classificationMode_.getKey(), "pre-integrated")) {
            LERROR("Pre-integration cannot be used with 2D tfs.");
            return;
        }
        tfTexture = tFunc2DPrimitives_.get()->getTexture();
        break;
    }
    default:
        LWARNING("Unsupported texture type!");
        return;
    }
    tfTexture->downloadTexture();

        // activate outport
    outport1_.activateTarget();
    outport1_.clearTarget();
    LGL_ERROR;


    // create output buffer
    tgt::vec4* output = new tgt::vec4[entryPort_.getSize().x * entryPort_.getSize().y];

    // download entry/exit point textures
    tgt::vec4* entryBuffer = reinterpret_cast<tgt::vec4*>(
        entryPort_.getColorTexture()->downloadTextureToBuffer(GL_RGBA, GL_FLOAT));
    tgt::vec4* exitBuffer = reinterpret_cast<tgt::vec4*>(
        exitPort_.getColorTexture()->downloadTextureToBuffer(GL_RGBA, GL_FLOAT));
    LGL_ERROR;

    // iterate over viewport and perform ray casting for each fragment
    for (int y=0; y < entryPort_.getSize().y; ++y) {
#ifdef VRN_MODULE_OPENMP
        #pragma omp parallel for schedule(dynamic) shared(volume,samplingStepSize,table,output,entryBuffer,exitBuffer,tfTexture,y,mode)
#endif
        for (int x=0; x < entryPort_.getSize().x; ++x) {
            vec4 gl_FragColor = vec4(0.f);
            int p = (y * entryPort_.getSize().x + x);

            vec4 frontPos = entryBuffer[p];
            vec4 backPos = exitBuffer[p];

            if ((frontPos == vec4(0.0)) && (backPos == vec4(0.0))) {
                //background needs no raycasting
            }
            else {
                //fragCoords are lying inside the boundingbox
                gl_FragColor = directRendering(frontPos.xyz(), backPos.xyz(), tfTexture, volume, samplingStepSize, table,mode, volumeInport_.getData()->getTextureToVoxelMatrix());
            }
            output[p] = gl_FragColor;
        }
    }
    delete[] entryBuffer;
    delete[] exitBuffer;

    tgt::Texture outputTex(tgt::ivec3(entryPort_.getSize().x, entryPort_.getSize().y, 1), GL_RGBA, GL_RGBA, GL_FLOAT,
                           tgt::Texture::LINEAR,tgt::Texture::REPEAT,(GLubyte*)output,true);
    outputTex.uploadTexture();

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    renderQuad();

    LGL_ERROR;

    outport1_.deactivateTarget();
    LGL_ERROR;
}

vec4 CPURaycaster::directRendering(const vec3& first, const vec3& last, tgt::Texture* tfTexture, const VolumeRAM* volume, float samplingStepSize, const PreIntegrationTable* table, CpuRaycasterClassificationMode mode, tgt::mat4 textureToVoxelMatrix) {

    tgt::svec3 volDim = volume->getDimensions();

    //get real world mapping
    RealWorldMapping rwm = volumeInport_.getData()->getRealWorldMapping();

    // calculate ray parameters
    float tend;
    float t = 0.0f;
    vec3 direction = last - first;
    // if direction is a nullvector the entry- and exitparams are the same
    // so special handling for tend is needed, otherwise we divide by zero
    // furthermore both for-loops will cause only 1 pass overall.
    // The test whether last and first are nullvectors is already done in main-function
    // but however the framerates are higher with this test.
    if (direction == vec3(0.0f) && last != vec3(0.0f) && first != vec3(0.0f))
        tend = samplingStepSize / 2.0f;
    else {
        tend = length(direction);
        direction = normalize(direction);
    }

    float lastIntensity = 0.f; //for pre-integration

    // ray-casting loop
    vec4 result = vec4(0.0f);
    float depthT = -1.0f;
    bool finished = false;
    for (int loop=0; !finished && loop<255*255; ++loop) {
        vec3 sample = first + t * direction;

        // convert sample pos to voxel position (including clamping) for accessing the volume
        vec3 sampleVoxelPos = tgt::clamp((textureToVoxelMatrix * tgt::vec4(sample, 1.f)).xyz(), tgt::vec3::zero, tgt::vec3(volDim) - tgt::vec3::one);

        float intensity = 0.f;
        if (texFilterMode_.getValue() == GL_NEAREST)
            intensity = volume->getVoxelNormalized(tgt::iround(sampleVoxelPos));
        else if (texFilterMode_.getValue() == GL_LINEAR)
            intensity = volume->getVoxelNormalizedLinear(sampleVoxelPos);
        else
            LERROR("Unknown texture filter mode");
        vec4 color = vec4(intensity);
        //apply realworld mapping and TF domain
        intensity = rwm.normalizedToRealWorld(intensity);

        //pre-integration
        if (mode == PREINTEGRATED) {
            if (!table)
                return vec4(0.f);
            intensity = tFunc1DKeys_.get()->realWorldToNormalized(intensity);
            color = table->classify(lastIntensity, intensity);
            lastIntensity = intensity;
        }
        else if (mode == TF) {
            if (tFuncType_.getValue() == TFT_1DKEYS) {
                intensity = tFunc1DKeys_.get()->realWorldToNormalized(intensity);
                //no shading is applied
                color = apply1DTF(tfTexture, intensity);
            }
            else {
                const VolumeRAM* volumeGradient = gradientVolumePort_.getData()->getRepresentation<VolumeRAM>();
                tgt::vec3 grad = tgt::vec3::zero;
                if (texFilterMode_.getValue() == GL_NEAREST) {
                    tgt::ivec3 iSample = tgt::iround(sampleVoxelPos);
                    grad.x = volumeGradient->getVoxelNormalized(iSample, 0);
                    grad.y = volumeGradient->getVoxelNormalized(iSample, 1);
                    grad.z = volumeGradient->getVoxelNormalized(iSample, 2);
                }
                else if (texFilterMode_.getValue() == GL_LINEAR) {
                    grad.x = volumeGradient->getVoxelNormalizedLinear(sampleVoxelPos, 0);
                    grad.y = volumeGradient->getVoxelNormalizedLinear(sampleVoxelPos, 1);
                    grad.z = volumeGradient->getVoxelNormalizedLinear(sampleVoxelPos, 2);
                }
                else
                    LERROR("Unknown texture filter mode");
                intensity = tFunc2DPrimitives_.get()->realWorldToNormalized(intensity,0);
                float gradMag = tFunc2DPrimitives_.get()->realWorldToNormalized(tgt::length(grad),1);
                color = apply2DTF(tfTexture, intensity, gradMag);
            }
            // apply opacity correction to accomodate for variable sampling intervals
            color.a = 1.f - pow(1.f - color.a, samplingStepSize * SAMPLING_BASE_INTERVAL_RCP);
        }
        else {
            //no classification
        }

        // perform compositing
        if (color.a > 0.0f) {

            //actual compositing
            result.xyz() = result.xyz() + (1.0f - result.a) * color.a * color.xyz();
            result.a = result.a + (1.0f - result.a) * color.a;
        }

        // save first hit ray parameter for depth value calculation
        if (depthT < 0.0f && result.a > 0.0f)
            depthT = t;

        // early ray termination
        if (result.a >= 0.95f) {
             result.a = 1.0f;
             finished = true;
        }

        t += samplingStepSize;
        finished = finished || (t > tend);

    } // ray-casting loop

    // calculate depth value from ray parameter (todo)
    // gl_FragDepth = 1.0;
    // if (depthT >= 0.0)
    //  gl_FragDepth = calculateDepthValue(depthT / tend);

    return result;
}

vec4 CPURaycaster::apply1DTF(tgt::Texture* tfTexture, float intensity) {
    int widthMinusOne = tfTexture->getWidth()-1;
    tgt::vec4 value = tgt::vec4(tfTexture->texel<tgt::col4>(static_cast<size_t>(tgt::clamp(tgt::iround(intensity * widthMinusOne), 0, widthMinusOne)))) / 255.f;
    return value;
}

vec4 CPURaycaster::apply2DTF(tgt::Texture* tfTexture, float intensity, float gradientMagnitude) {
    vec4 value = vec4(tfTexture->texel<tgt::vec4>(size_t(intensity * (tfTexture->getWidth()-1)),
        size_t(gradientMagnitude * (tfTexture->getHeight()-1))));
    return value;
}

} // namespace voreen
