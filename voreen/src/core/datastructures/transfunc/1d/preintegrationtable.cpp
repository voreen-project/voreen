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

#include "voreen/core/datastructures/transfunc/1d/preintegrationtable.h"

#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"

#include "math.h"
#include "tgt/texture.h"
#include "tgt/immediatemode/immediatemode.h"

namespace voreen {

PreIntegrationTable::PreIntegrationTable(TransFunc1D* transFunc, size_t resolution, float d, bool useIntegral, bool computeOnGPU, tgt::Shader* program)
    : transFunc_(transFunc), resolution_(resolution), samplingStepSize_(d), useIntegral_(useIntegral), computeOnGPU_(computeOnGPU), table_(0), tex_(0), program_(program)
{

    //check values
    if (resolution <= 1)
        resolution_ = 256;
    if (d <= 0.f)
        samplingStepSize_ = 1.f;

    //initialize render target if necessary
    if (computeOnGPU_) {
        renderTarget_.initialize(GL_RGBA16, GL_DEPTH_COMPONENT24);
        renderTarget_.setDebugLabel("Pre-Integration Table");
        renderTarget_.resize(tgt::vec2(static_cast<float>(resolution_)));
    }
}

PreIntegrationTable::~PreIntegrationTable() {

    // delete table if no texture has been created (otherwise the texture deletes the data)
    if (!tex_ || computeOnGPU_)
        delete[] table_;

    //deinitialize render target or delete texture
    if (computeOnGPU_)
        renderTarget_.deinitialize();
    else
        delete tex_;
}

void PreIntegrationTable::computeTable() const {
    if (!transFunc_)
        return;

    //buffer for TF values
    tgt::vec4* tfBuffer = new tgt::vec4[resolution_];

    int front_end = tgt::iround(transFunc_->getThreshold().x * resolution_);
    int back_start = tgt::iround(transFunc_->getThreshold().y * resolution_);
    //all values before front_end and after back_start are set to zero
    //all other values remain the same
    for (int i = 0; i < front_end; ++i)
        tfBuffer[i] = tgt::vec4(0.f);

    for (int i = front_end; i < back_start; ++i) {
        //fetch current value from TF
        float intensity = static_cast<float>(i) / static_cast<float>(resolution_ - 1);
        tgt::vec4 value;
        if(transFunc_->getDataType() == TransFuncBase::TF_FLOAT) {
            tgt::Vector4<GLfloat> floatValue = transFunc_->getMappingForValueFloat(transFunc_->applyGammaToIntensity(intensity));
            value = tgt::vec4(floatValue);
        } else {
            tgt::Vector4<GLubyte> ubyteValue = transFunc_->getMappingForValueUByte(transFunc_->applyGammaToIntensity(intensity));
            value = tgt::vec4(ubyteValue) / 255.f;
        }
        if(transFunc_->getAlphaMode() == TransFuncBase::TF_ONE_ALPHA && value.a != 0.f)
            value.a = 1.f;
        else if (transFunc_->getAlphaMode() == TransFuncBase::TF_ZERO_ALPHA)
            value.a = 0.f;
        tfBuffer[i] = value;

    }

    for (int i = back_start; i < static_cast<int>(resolution_); ++i)
        tfBuffer[i] = tgt::vec4(0.f);

    if(!useIntegral_) { // Correct (but slow) calculation of PI-table:
        int lookupindex = 0;
        for (int sb = 0; sb < static_cast<int>(resolution_); ++sb) {
            for (int sf = 0; sf < static_cast<int>(resolution_); ++sf) {
                if (sb != sf) {
                    float scale = 1.0f / (fabs(static_cast<float>(sb - sf)) + 1);

                    int incr = 1;
                    if(sb < sf)
                        incr = -1;

                    tgt::vec4 result = tgt::vec4(0.0f);
                    for(int s = sf; (incr == 1 ? s<=sb : s>=sb) && (result.a < 0.95); s += incr) {

                        tgt::vec4 curCol = tfBuffer[s];

                        if (curCol.a > 0.0f) {
                            // apply opacity correction to accomodate for variable sampling intervals
                            curCol.a = 1.f - pow(1.f - curCol.a, samplingStepSize_ * 200.0f * scale);

                            //actual compositing
                            tgt::vec3 result_rgb = result.xyz() + (1.0f - result.a) * curCol.a * curCol.xyz();
                            result.a = result.a + (1.0f - result.a) * curCol.a;

                            result.xyz() = result_rgb;
                        }
                    }
                    result.xyz() /= std::max(result.a, 0.001f);

                    //result.a = 1.f - pow(1.f - result.a, 1.0f / (samplingStepSize_ * 200.0f));
                    table_[lookupindex] = tgt::clamp(result, 0.f, 1.f);
                } else {
                    tgt::vec4 result = tgt::clamp(tfBuffer[sf], 0.f, 1.f);

                    // apply opacity correction to accomodate for variable sampling intervals
                    result.a = 1.f - pow(1.f - result.a, samplingStepSize_ * 200.0f);

                    table_[lookupindex] = result;
                }
                lookupindex++;
            }
        }
    }
    else { //faster version using integral functions, see Real-Time Volume Graphics, p96
        //compute integral functions
        tgt::vec4* intFunc = new tgt::vec4[resolution_];

        tgt::vec4 accumResult(0.f);
        tgt::vec4 curColor;

        for (int i = 0; i < static_cast<int>(resolution_); ++i) {
            //fetch current value from TF
            //float nIndex = static_cast<float>(i) / static_cast<float>(resolution_ - 1);
            //vec4 curCol = apply1DTF(nIndex);
            tgt::vec4 curCol = tfBuffer[i];

            //calculate new integral function
            if (curCol.a > 0.0f) {
                //actual compositing
                accumResult.xyz() += curCol.xyz() * curCol.a;
                accumResult.a += curCol.a;
            }
            intFunc[i] = accumResult;
        }
        float factor;

        int lookupindex = 0;
        int endIndex = static_cast<int>(resolution_);

        // compute look-up table from integral functions
        for (int sb = 0; sb < endIndex; ++sb)
            for (int sf = 0; sf < endIndex; ++sf) {

                int smin = std::min(sb, sf);
                int smax = std::max(sb, sf);

                tgt::vec4 col;
                if (smax != smin) {
                    factor = samplingStepSize_ * 200.f / static_cast<float>(smax - smin);

                    col.xyz() = (intFunc[smax].xyz() - intFunc[smin].xyz()) * factor;
                    col.a = 1.f - exp(-(intFunc[smax].a - intFunc[smin].a) * factor);

                    col.xyz() /= std::max(col.a, 0.001f);

                    //col.a = 1.f - pow(1.f - col.a, 1.0f / (samplingStepSize_ * 200.0f));
                } else {
                    //float nIndex = static_cast<float>(smin) / static_cast<float>(resolution_ - 1);
                    //col = apply1DTF(nIndex);
                    col = tfBuffer[smin];
                    // apply opacity correction
                    col.a = 1.f - pow(1.f - col.a, samplingStepSize_ * 200.0f);

                }
                table_[lookupindex] = tgt::clamp(col, 0.f, 1.f);
                lookupindex++;
            }
        delete[] intFunc;
    }

    delete[] tfBuffer;
}

void PreIntegrationTable::computeTableGPU() const {

    // Get currently active shader.
    tgt::Shader* previousShader = ShdrMgr.getActiveShader();

    //render pre-integration texture into render target
    renderTarget_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //get texture of transfer function
    tgt::TextureUnit transferUnit;
    transferUnit.activate();

    transFunc_->getTexture()->bind();
    transFunc_->getTexture()->enable();

    // activate shader program
    program_->activate();
    IMode.setMatstackUniforms(program_);

    program_->setUniform("samplingStepSize_", samplingStepSize_);

    //set transfer function texture
    program_->setUniform("tfTex_", transferUnit.getUnitNumber());

    bool oldIgnoreError = program_->getIgnoreUniformLocationError();
    program_->setIgnoreUniformLocationError(true);
    program_->setUniform("texParams_.dimensions_", tgt::vec2((float)resolution_));
    program_->setUniform("texParams_.dimensionsRCP_", tgt::vec2(1.f) / tgt::vec2((float)resolution_));
    program_->setUniform("texParams_.matrix_", tgt::mat4::identity);
    program_->setIgnoreUniformLocationError(oldIgnoreError);

    RenderProcessor::renderQuad();

    //clean up
    transFunc_->getTexture()->disable();
    renderTarget_.deactivateTarget();
    program_->deactivate();

    // Restore shader.
    if (previousShader)
        previousShader->activate();

    tgt::TextureUnit::setZeroUnit();
    LGL_ERROR;

    //set output texture
    tex_ = renderTarget_.getColorTexture();
}

tgt::vec4 PreIntegrationTable::classify(float fs, float fe) const {

    //lazy computation
    if (!table_) {
        //create the table
        table_ = new tgt::vec4[resolution_ * resolution_];
        computeTable();
    }

    if (!table_)
        return tgt::vec4(0.f,0.f,0.f,0.f);
    else {
        //compute indices to table
        int i1 = tgt::iround(fs * (resolution_ - 1));
        int i2 = tgt::iround(fe * (resolution_ - 1));

        //check indices
        i1 = tgt::clamp(i1, 0, static_cast<int>(resolution_) - 1);
        i2 = tgt::clamp(i2, 0, static_cast<int>(resolution_) - 1);

        //table lookup
        return table_[i2 * resolution_ + i1];
    }
}

const tgt::vec4 * const PreIntegrationTable::getTable() const {

    if (computeOnGPU_)
        return 0;
    else if (!table_) {     //lazy computation
        //create the table
        table_ = new tgt::vec4[resolution_ * resolution_];
        computeTable();
    }

    return table_;
}

void PreIntegrationTable::createTexFromTable() const {
    if (tex_) {
        delete tex_;
        tex_ = 0;
        table_ = 0;     // texture also deletes table...
    }

    //lazy computation
    if (!table_) {
        //create the table
        table_ = new tgt::vec4[resolution_ * resolution_];
        computeTable();
    }

    tex_ = new tgt::Texture(tgt::ivec3(static_cast<int>(resolution_),static_cast<int>(resolution_),1), GL_RGBA, GL_RGBA32F, GL_FLOAT,
                            tgt::Texture::LINEAR, tgt::Texture::CLAMP, reinterpret_cast<GLubyte*>(table_), true);
    LGL_ERROR;
    tex_->uploadTexture();
}

const tgt::Texture* PreIntegrationTable::getTexture() const {

    //check if texture has to be created and if computation should be made on the gpu
    if (!tex_) {
        if (computeOnGPU_ /*|| useIntegral_*/)
            computeTableGPU();
        else
            createTexFromTable();
    }

    tgtAssert(tex_, "No texture");

    return tex_;
}

float PreIntegrationTable::getSamplingStepSize() const {
    return samplingStepSize_;
}

bool PreIntegrationTable::usesIntegral() const {
    return useIntegral_;
}

size_t PreIntegrationTable::getDimension() const {
    return resolution_;
}

} //namespace
