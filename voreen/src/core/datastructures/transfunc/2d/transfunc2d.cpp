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

#include "voreen/core/datastructures/transfunc/2d/transfunc2d.h"

#include "tgt/shadermanager.h"

namespace voreen {

const std::string TransFunc2D::loggerCat_("voreen.TransFunc2D");

TransFunc2D::TransFunc2D(int width, int height, DataType dataType, tgt::Texture::Filter filter)
    : TransFuncBase(width, height, 1, dataType, filter)
{
    gammaValues_ = tgt::vec2(1.f,1.f);
    domains_[0] = tgt::vec2(0.f, 1.f); domains_[1] = tgt::vec2(0.f, 1.f);
    thresholds_[0] = tgt::vec2(0.f, 1.f); thresholds_[1] = tgt::vec2(0.f, 1.f);
}

TransFunc2D::~TransFunc2D() {
}

void TransFunc2D::setMemberValuesFrom(const TransFuncBase* transfunc) {
    tgtAssert(transfunc, "null pointer passed");
    TransFuncBase::setMemberValuesFrom(transfunc);

    if(const TransFunc2D* tmpFunc = dynamic_cast<const TransFunc2D*>(transfunc)) {
        gammaValues_ = tmpFunc->gammaValues_;
        domains_[0] = tmpFunc->domains_[0];
        domains_[1] = tmpFunc->domains_[1];
        thresholds_[0] = tmpFunc->thresholds_[0];
        thresholds_[1] = tmpFunc->thresholds_[1];
    } else {
        LWARNING("setMemberValuesFrom(): passed parameter is not of type TransFunc2D");
    }
}

bool TransFunc2D::compareTo(const TransFuncBase& tf) const {
    if(const TransFunc2D* tf2d = dynamic_cast<const TransFunc2D*>(&tf)) {
        return ((TransFuncBase::compareTo(tf))             &&
                (gammaValues_ == tf2d->gammaValues_)       &&
                (domains_[0] == tf2d->domains_[0])         &&
                (domains_[1] == tf2d->domains_[1])         &&
                (thresholds_[0] == tf2d->thresholds_[0])   &&
                (thresholds_[1] == tf2d->thresholds_[1]));
    }
    return false;
}

    //--------------------------------------
    //  texture handling
    //--------------------------------------
GLvoid* TransFunc2D::createTextureData() {
   int numValues = tgt::hmul(dimensions_);

    int x_threshold_front_end  = tgt::iround(thresholds_[0].x*dimensions_.x);
    int x_threshold_back_start = tgt::iround(thresholds_[0].y*dimensions_.x);
    int y_threshold_front_end  = tgt::iround(thresholds_[1].x*dimensions_.y);
    int y_threshold_back_start = tgt::iround(thresholds_[1].y*dimensions_.y);

    GLvoid* returnArray = 0;
    if (dataType_ == TF_UBYTE) {
        tgt::Vector4<GLubyte>* ubyteData = new tgt::Vector4<GLubyte>[numValues]; //always GL_RGBA
        for(int y = 0; y < dimensions_.y; y++) {
            if(y < y_threshold_front_end || y > y_threshold_back_start) {
                // set zero
                std::fill_n(&ubyteData[y*dimensions_.x],dimensions_.x,tgt::Vector4<GLubyte>(0,0,0,0));
            } else {
                for(int x = 0; x < dimensions_.x; x++) {
                    if(x < x_threshold_front_end || x > x_threshold_back_start) {
                        // set zero
                        ubyteData[y*dimensions_.x+x] = tgt::Vector4<GLubyte>(0,0,0,0);
                    } else {
                        tgt::vec2 values = applyGammaToValues(static_cast<float>(x) / (dimensions_.x-1), static_cast<float>(y) / (dimensions_.y-1));
                        ubyteData[y*dimensions_.x+x] = getMappingForValueUByte(values.x,values.y);
                    }
                }
            }
        }
        returnArray = reinterpret_cast<GLvoid*>(ubyteData);
    } else { //TF_FLOAT
        tgt::Vector4<GLfloat>* floatData = new tgt::Vector4<GLfloat>[numValues]; //always GL_RGBA
        for(int y = 0; y < dimensions_.y; y++) {
            if(y < y_threshold_front_end || y > y_threshold_back_start) {
                // set zero
                std::fill_n(&floatData[y*dimensions_.x],dimensions_.x,tgt::Vector4<GLfloat>(0.f,0.f,0.f,0.f));
            } else {
                for(int x = 0; x < dimensions_.x; x++) {
                    if(x < x_threshold_front_end || x > x_threshold_back_start) {
                        // set zero
                        floatData[y*dimensions_.x+x] = tgt::Vector4<GLfloat>(0.f,0.f,0.f,0.f);
                    } else {
                        tgt::vec2 values = applyGammaToValues(static_cast<float>(x) / (dimensions_.x-1), static_cast<float>(y) / (dimensions_.y-1));
                        floatData[y*dimensions_.x+x] = getMappingForValueFloat(values.x, values.y);
                    }
                }
            }
        }
        returnArray = reinterpret_cast<GLvoid*>(floatData);
    }
    return returnArray;
}

    //--------------------------------------
    //  gamma handling
    //--------------------------------------
void TransFunc2D::setGammaValue(tgt::vec2 gamma) {
    if(gamma != gammaValues_) {
        gammaValues_ = gamma;
        invalidateTexture();
    }
}

tgt::vec2 TransFunc2D::getGammaValue() const {
    return gammaValues_;
}

tgt::vec2 TransFunc2D::applyGammaToValues(float x, float y) const {
    tgt::vec2 result(x,y);
    if(gammaValues_.x != 1.f) {
        result.x = powf(x,gammaValues_.x);
    }
    if(gammaValues_.y != 1.f) {
        result.y = powf(y,gammaValues_.y);
    }
    return result;
}

    //--------------------------------------
    //  domain handling
    //--------------------------------------
void TransFunc2D::setDomain(float lower, float upper, size_t dimension) {
    tgtAssert(dimension < 2, "dimension must be 0 or 1");
    domains_[dimension] = tgt::vec2(lower,upper);


}

void TransFunc2D::setDomain(const tgt::vec2& domain, size_t dimension) {
    setDomain(domain.x,domain.y,dimension);
}

tgt::vec2 TransFunc2D::getDomain(size_t dimension) const {
    tgtAssert(dimension < 2, "dimension must be 0 or 1");
    return domains_[dimension];
}

    //--------------------------------------
    //  threshold handling
    //--------------------------------------
void TransFunc2D::setThreshold(float lower, float upper, size_t dimension) {
    tgtAssert(dimension < 2, "dimension must be 0 or 1");
    thresholds_[dimension] = tgt::vec2(lower,upper);
    invalidateTexture();
}

void TransFunc2D::setThreshold(const tgt::vec2& thresholds, size_t dimension) {
    setThreshold(thresholds.x,thresholds.y,dimension);
}

tgt::vec2 TransFunc2D::getThreshold(size_t dimension) const {
    tgtAssert(dimension < 2, "dimension must be 0 or 1");
    return thresholds_[dimension];
}

    //--------------------------------------
    //  shader defines
    //--------------------------------------
void TransFunc2D::setUniform(tgt::Shader* shader, const std::string& uniform, const std::string& uniformTex, const GLint texUnit) {
    tgtAssert(shader, "Null pointer passed");
    bool oldIgnoreError = shader->getIgnoreUniformLocationError();
    shader->setIgnoreUniformLocationError(true);
    shader->setUniform(uniformTex, texUnit);

    shader->setUniform(uniform + ".domainLower_", tgt::vec3(domains_[0].x, domains_[1].x, 0.0f));
    shader->setUniform(uniform + ".domainUpper_", tgt::vec3(domains_[0].y, domains_[1].y, 0.0f));

    shader->setIgnoreUniformLocationError(oldIgnoreError);
}

    //--------------------------------------
    //  load and save
    //--------------------------------------
void TransFunc2D::serialize(Serializer& s) const {
    TransFuncBase::serialize(s);
    s.serialize("gammaValues", gammaValues_);
    tgt::vec4 tmp(domains_[0].x,domains_[0].y,domains_[1].x,domains_[1].y);
    s.serialize("domains", tmp);
    tmp = tgt::vec4(thresholds_[0].x,thresholds_[0].y,thresholds_[1].x,thresholds_[1].y);
    s.serialize("thresholds", tmp);
}

void TransFunc2D::deserialize(Deserializer& s) {
    TransFuncBase::deserialize(s);
    s.optionalDeserialize("gammaValues", gammaValues_, tgt::vec2(1.f,1.f));
    tgt::vec4 defaultValue(0.f, 1.f, 0.f ,1.f);
    tgt::vec4 tmp;
    s.optionalDeserialize("domains", tmp, defaultValue);
    domains_[0] = tgt::vec2(tmp.x,tmp.y); domains_[1] = tgt::vec2(tmp.z,tmp.a);
    s.optionalDeserialize("thresholds", tmp, defaultValue);
    thresholds_[0] = tgt::vec2(tmp.x,tmp.y); thresholds_[1] = tgt::vec2(tmp.z,tmp.a);
}

    //--------------------------------------
    //  handle real world mapping
    //--------------------------------------
float TransFunc2D::realWorldToNormalized(float rw, int dimension) const {
    return TransFuncBase::realWorldToNormalized(rw, domains_[dimension]);
}

tgt::vec2 TransFunc2D::realWorldToNormalized(tgt::vec2 rw) const {
    return tgt::vec2(realWorldToNormalized(rw.x, 0), realWorldToNormalized(rw.y, 1));
}

float TransFunc2D::normalizedToRealWorld(float n, int dimension) const {
    return TransFuncBase::normalizedToRealWorld(n, domains_[dimension]);
}

tgt::vec2 TransFunc2D::normalizedToRealWorld(tgt::vec2 n) const {
    return tgt::vec2(normalizedToRealWorld(n.x, 0), normalizedToRealWorld(n.y, 1));
}


} // namespace voreen
