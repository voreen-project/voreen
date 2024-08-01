/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"

#include "voreen/core/datastructures/transfunc/1d/preintegrationtable.h"

#include "voreen/core/utils/glsl.h"

namespace voreen {

const std::string TransFunc1D::loggerCat_("voreen.TransFunc1D");

TransFunc1D::TransFunc1D(int width, DataType dataType, tgt::Texture::Filter filter)
    : TransFuncBase(width, 1, 1, dataType, filter)
    , gammaValue_(1.f)
    , domain_(0.f, 1.f)
    , threshold_(0.f, 1.f)
    , preIntegrationTableMap_(NULL, 10)
    , preIntegrationProgram_(0)
{
    preIntegrationTableMap_.setTransFunc(this);
}

TransFunc1D::~TransFunc1D() {
    preIntegrationTableMap_.clear();

    if (preIntegrationProgram_)
        ShdrMgr.dispose(preIntegrationProgram_);
}

void TransFunc1D::setMemberValuesFrom(const TransFuncBase* transfunc) {
    tgtAssert(transfunc, "null pointer passed");
    TransFuncBase::setMemberValuesFrom(transfunc);

    if(const TransFunc1D* tmpFunc = dynamic_cast<const TransFunc1D*>(transfunc)) {
        gammaValue_ = tmpFunc->gammaValue_;
        domain_ = tmpFunc->domain_;
        threshold_ = tmpFunc->threshold_;
    } else {
        LWARNING("setMemberValuesFrom(): passed parameter is not of type TransFunc1D");
    }
}

bool TransFunc1D::compareTo(const TransFuncBase& tf) const {
    if(const TransFunc1D* tf1d = dynamic_cast<const TransFunc1D*>(&tf)) {
        return ((TransFuncBase::compareTo(tf))     &&
                (gammaValue_ == tf1d->gammaValue_)  &&
                (domain_ == tf1d->domain_)         &&
                (threshold_ == tf1d->threshold_));
    }
    return false;
}

    //--------------------------------------
    //  texture handling
    //--------------------------------------
void TransFunc1D::invalidateTexture() {
    TransFuncBase::invalidateTexture();
    preIntegrationTableMap_.clear();
}

GLvoid* TransFunc1D::createTextureData() {
    int numValues = tgt::hmul(dimensions_);

    int threshold_front_end  = tgt::iround(threshold_.x*dimensions_.x);
    int threshold_back_start = tgt::iround(threshold_.y*dimensions_.x);

    GLvoid* returnArray = 0;
    if (dataType_ == TF_UBYTE) {
        tgt::Vector4<GLubyte>* ubyteData = new tgt::Vector4<GLubyte>[numValues]; //always GL_RGBA
        for(int x = 0; x < dimensions_.x; x++) {
            if(x < threshold_front_end || x > threshold_back_start) {
                // set zero
                ubyteData[x] = tgt::Vector4<GLubyte>(0,0,0,0);
            } else {
                ubyteData[x] = getMappingForValueUByte(applyGammaToIntensity(static_cast<float>(x) / (dimensions_.x-1)));
            }
        }
        returnArray = reinterpret_cast<GLvoid*>(ubyteData);
    } else { //TF_FLOAT
        tgt::Vector4<GLfloat>* floatData = new tgt::Vector4<GLfloat>[numValues]; //always GL_RGBA
        for(int x = 0; x < dimensions_.x; x++) {
            if(x < threshold_front_end || x > threshold_back_start) {
                // set zero
                floatData[x] = tgt::Vector4<GLfloat>(0,0,0,0);
            } else {
                floatData[x] = getMappingForValueFloat(applyGammaToIntensity(static_cast<float>(x) / (dimensions_.x-1)));
            }
        }
        returnArray = reinterpret_cast<GLvoid*>(floatData);
    }
    return returnArray;
}

    //--------------------------------------
    //  pre-integration handling
    //--------------------------------------
const PreIntegrationTable* TransFunc1D::getPreIntegrationTable(float samplingStepSize, size_t dimension, bool useIntegral, bool computeOnGPU) {
    if (computeOnGPU && !preIntegrationProgram_) {
        preIntegrationProgram_ = ShdrMgr.loadSeparate("passthrough.vert", "preintegration/preintegration.frag", GLSL::generateStandardShaderHeader(0), false);
        tgtAssert(preIntegrationProgram_, "no preintegration shader program");
    }

    return preIntegrationTableMap_.getPreIntegrationTable(samplingStepSize,dimension,useIntegral,computeOnGPU, preIntegrationProgram_);
}

    //--------------------------------------
    //  gamma handling
    //--------------------------------------
void TransFunc1D::setGammaValue(float gamma) {
    if(gamma != gammaValue_) {
        gammaValue_ = gamma;
        invalidateTexture();
    }
}

float TransFunc1D::getGammaValue() const {
    return gammaValue_;
}

float TransFunc1D::applyGammaToIntensity(float value) const {
    if(gammaValue_ != 1.f) {
        return powf(value,gammaValue_);
    }
    return value;
}
    //--------------------------------------
    //  domain handling
    //--------------------------------------
void TransFunc1D::setDomain(float lower, float upper) {
    domain_ = tgt::vec2(lower,upper);
}

void TransFunc1D::setDomain(const tgt::vec2& domain) {
    setDomain(domain.x,domain.y);
}

tgt::vec2 TransFunc1D::getDomain() const {
    return domain_;
}

    //--------------------------------------
    //  threshold handling
    //--------------------------------------
void TransFunc1D::setThreshold(float lower, float upper) {
    threshold_ = tgt::vec2(lower,upper);
    invalidateTexture();
}

void TransFunc1D::setThreshold(const tgt::vec2& threshold) {
    setThreshold(threshold.x,threshold.y);
}

tgt::vec2 TransFunc1D::getThreshold() const {
    return threshold_;
}

    //--------------------------------------
    //  shader defines
    //--------------------------------------
void TransFunc1D::setUniform(tgt::Shader* shader, const std::string& uniform, const std::string& uniformTex, const GLint texUnit) {
    tgtAssert(shader, "Null pointer passed");
    bool oldIgnoreError = shader->getIgnoreUniformLocationError();
    shader->setIgnoreUniformLocationError(true);

    shader->setUniform(uniformTex, texUnit);
    shader->setUniform(uniform + ".domainLower_", tgt::vec3(domain_.x, 0.0f, 0.0f));
    shader->setUniform(uniform + ".domainUpper_", tgt::vec3(domain_.y, 0.0f, 0.0f));

    shader->setIgnoreUniformLocationError(oldIgnoreError);
}

    //--------------------------------------
    //  load and save
    //--------------------------------------
void TransFunc1D::serialize(Serializer& s) const {
    TransFuncBase::serialize(s);
    s.serialize("gammaValue", gammaValue_);
    s.serialize("domain", domain_);
    s.serialize("threshold", threshold_);
}

void TransFunc1D::deserialize(Deserializer& s) {
    TransFuncBase::deserialize(s);
    s.optionalDeserialize("gammaValue", gammaValue_, 1.f);
    tgt::vec2 tmp(0.f,1.f);
    s.optionalDeserialize("domain", domain_, tmp);
    s.optionalDeserialize("threshold", threshold_, tmp);
}

    //--------------------------------------
    //  handle real world mapping
    //--------------------------------------
float TransFunc1D::realWorldToNormalized(float rw) const {
    return TransFuncBase::realWorldToNormalized(rw, domain_);
}

float TransFunc1D::normalizedToRealWorld(float n) const {
    return TransFuncBase::normalizedToRealWorld(n, domain_);
}

} // namespace voreen
