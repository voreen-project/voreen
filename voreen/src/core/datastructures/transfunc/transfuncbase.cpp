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

#include "voreen/core/datastructures/transfunc/transfuncbase.h"

#include "tgt/glcontextmanager.h"

namespace voreen {

const std::string TransFuncBase::loggerCat_("voreen.TransFuncBase");

TransFuncBase::TransFuncBase(int width, int height, int depth, DataType dataType, tgt::Texture::Filter filter)
    : tex_(0)
    , dimensions_(width, height, depth)
    , dataType_(dataType)
    , filter_(filter)
    , textureInvalid_(true)
    , alphaMode_(TF_USE_ALPHA)
{
}

TransFuncBase::~TransFuncBase() {
    if (tex_) {
        delete tex_;
        tex_ = 0;
        LGL_ERROR;
    }
}

tgt::ivec3 TransFuncBase::getDimensions() const {
    return dimensions_;
}

void TransFuncBase::setMemberValuesFrom(const TransFuncBase* transfunc) {
    tgtAssert(transfunc, "null pointer passed");

    dimensions_ = transfunc->dimensions_;
    dataType_ = transfunc->dataType_;
    filter_ = transfunc->filter_;
    alphaMode_ = transfunc->alphaMode_;

    //Crash in using Voreen without OpenGL support
    /*tgt::Texture* tex = const_cast<TransFuncBase*>(transfunc)->getTexture();
    GLubyte* pixels = tex->getPixelData();
    int factor = tex->getBpp();
    //update texture
    createTex();
    GLubyte* newpixels = new GLubyte[tgt::hmul(dimensions_) * factor];
    for (int i=0;i<dimensions_.x*dimensions_.y*dimensions_.z*factor;i++) {
        newpixels[i] = pixels[i];
    }
    tex_->setPixelData(newpixels);*/
}

bool TransFuncBase::compareTo(const TransFuncBase& tf) const {
    return ((dimensions_ == tf.dimensions_)         &&
            (filter_     == tf.filter_)             &&
            (dataType_   == tf.dataType_)   &&
            (alphaMode_  == tf.alphaMode_));
}

    //--------------------------------------
    //  handle tf texture
    //--------------------------------------
TransFuncBase::DataType TransFuncBase::getDataType() const {
    return dataType_;
}

void TransFuncBase::invalidateTexture() {
    textureInvalid_ = true;
}

bool TransFuncBase::isTextureValid() const {
    return !textureInvalid_;
}

void TransFuncBase::resize(int width, int height, int depth) {
    int maxTexSize = 1024;
    if (tgt::Singleton<tgt::GpuCapabilities>::isInited()) {
        if (depth == 1)
            maxTexSize = GpuCaps.getMaxTextureSize();
        else
            maxTexSize = GpuCaps.getMax3DTextureSize();
    }
    if (maxTexSize < width)
        width = maxTexSize;

    if (maxTexSize < height)
        height = maxTexSize;

    if (maxTexSize < depth)
        depth = maxTexSize;

    if (width != dimensions_.x) {
        dimensions_.x = width;
        invalidateTexture();
    }
    if (height != dimensions_.y) {
        dimensions_.y = height;
        invalidateTexture();
    }
    if (depth != dimensions_.z) {
        dimensions_.z = depth;
        invalidateTexture();
    }
}


void TransFuncBase::createTex() {
    delete tex_;

    tex_ = new tgt::Texture(dimensions_, (GLint)GL_RGBA, (GLint)GL_RGBA, (GLenum)(dataType_ == TF_FLOAT ? GL_FLOAT : GL_UNSIGNED_BYTE), filter_,tgt::Texture::CLAMP);
    LGL_ERROR;
}

void TransFuncBase::updateTexture() {
    if (!tex_ || (tex_->getDimensions() != dimensions_))
        createTex();

    if (!tex_) {
        LERROR("Failed to create texture");
        return;
    }

    if(textureInvalid_) {
        GLvoid* tmp = createTextureData();
        applyAlpha(tmp);
        tex_->setCpuTextureData(reinterpret_cast<GLubyte*>(tmp),true);
    }

    tex_->uploadTexture();
    textureInvalid_ = false;
}

tgt::Texture* TransFuncBase::getTexture() const {
    if (textureInvalid_)
        const_cast<TransFuncBase*>(this)->updateTexture(); //HACK
    return tex_;
}

    //--------------------------------------
    //  handle alpha
    //--------------------------------------
void TransFuncBase::setAlphaMode(AlphaMode mode) {
    if(mode != alphaMode_) {
        alphaMode_ = mode;
        invalidateTexture();
    }
}

TransFuncBase::AlphaMode TransFuncBase::getAlphaMode() const {
    return alphaMode_;
}

void TransFuncBase::applyAlpha(GLvoid* data) const {
    if(alphaMode_ == TF_USE_ALPHA) return;

    int texelCount = tgt::hmul(dimensions_);

    GLubyte* ubyteData = 0;
    GLfloat* floatData = 0;

    if(dataType_ == TF_UBYTE) {
        ubyteData = reinterpret_cast<GLubyte*>(data);
        //we always have GL_RGBA && GL_UNSIGNED_BYTE
        if(alphaMode_ == TF_ONE_ALPHA) {
            for(int i = 0; i < texelCount; i++) {
                if(*(ubyteData + 3 + i*4)) //replace only, if unequal zero
                *(ubyteData + 3 + i*4) = static_cast<GLubyte>(255);
            }
        } else { // TF_ZERO_ALPHA
            for(int i = 0; i < texelCount; i++)
                *(ubyteData + 3 + i*4) = static_cast<GLubyte>(0);
        }
    } else {
        floatData = reinterpret_cast<GLfloat*>(data);
        //we always have GL_RGBA && GL_FLOAT
        if(alphaMode_ == TF_ONE_ALPHA) {
            for(int i = 0; i < texelCount; i++) {
                if(*(floatData + 3 + i*4))
                    *(floatData + 3 + i*4) = static_cast<GLfloat>(1.f);
            }
        } else { // TF_ZERO_ALPHA
            for(int i = 0; i < texelCount; i++)
                *(floatData + 3 + i*4) = static_cast<GLfloat>(0.f);
        }
    }
}

    //--------------------------------------
    //  load and save
    //--------------------------------------
void TransFuncBase::serialize(Serializer& s) const {
    s.serialize("alphaMode", static_cast<int>(alphaMode_));
    s.serialize("texDimensions", dimensions_);
    s.serialize("filter", filter_);
    s.serialize("dataType", dataType_);
}

void TransFuncBase::deserialize(Deserializer& s) {
    int tmp;
    s.optionalDeserialize("alphaMode", tmp, static_cast<int>(TF_USE_ALPHA));
    alphaMode_ = static_cast<TransFuncBase::AlphaMode>(tmp);
    s.optionalDeserialize("texDimensions", dimensions_, dimensions_);
    s.optionalDeserialize("filter", tmp, static_cast<int>(tgt::Texture::NEAREST));
    filter_ = static_cast<tgt::Texture::Filter>(tmp);
    s.optionalDeserialize("dataType", tmp, static_cast<int>(TF_UBYTE));
    dataType_ = static_cast<DataType>(tmp);
    invalidateTexture();
}

    //--------------------------------------
    //  shader defines
    //--------------------------------------
std::string TransFuncBase::getShaderDefines(const std::string& defineName) const {
    return ("#define " + defineName + " " + getSamplerType() + "\n");
}

    //--------------------------------------
    //  real world mapping helper
    //--------------------------------------
float TransFuncBase::realWorldToNormalized(float rw, const tgt::vec2& domain) {
    tgtAssert(domain.x <= domain.y, "invalid transfer function domain");
    if (rw < domain.x)
        return 0.f;
    else if (rw > domain.y)
        return 1.f;
    else if(domain.x == domain.y)
        return 0.f;
    else {
        if (domain.y < domain.x) { //< handle invalid domain gracefully in release mode
            LERROR("Invalid transfer function domain:" << domain);
            return 1.f;
        }
        return (rw - domain.x) / (domain.y - domain.x);
    }
}

float TransFuncBase::normalizedToRealWorld(float n, const tgt::vec2& domain) {
    tgtAssert(domain.x <= domain.y, "invalid domain");
    return domain.x + (domain.y - domain.x) * n;
}

} // namespace voreen
