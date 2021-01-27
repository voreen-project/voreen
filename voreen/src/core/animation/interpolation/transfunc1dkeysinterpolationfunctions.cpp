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

#include "voreen/core/animation/interpolation/transfunc1dkeysinterpolationfunctions.h"
#include "voreen/core/datastructures/transfunc/1d/1dkeys/utils/transfuncmappingkey.h"
#include "voreen/core/animation/interpolation/basicfloatinterpolation.h"
#include "voreen/core/animation/interpolation/basicintinterpolation.h"
#include <math.h>

namespace voreen {

GLubyte* TransFunc1DKeysInterpolationFunctionBase::convertTextureToRGBA(tgt::ivec3 dim, GLubyte* texture, GLuint inputformat) {
    GLubyte* data = new GLubyte[4 * dim.x * dim.y * dim.z];
    for (int x = 0; x < dim.x; ++x) {
        for (int y = 0; y < dim.y; ++y) {
            for (int z = 0; z < dim.z; ++z) {
                switch (inputformat) {
                case GL_RED:
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+0] = texture[x*dim.y*dim.z+y*dim.z+z];
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+1] = 0;
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+2] = 0;
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+3] = 1;
                    break;
                case GL_GREEN:
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+0] = 0;
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+1] = texture[x*dim.y*dim.z+y*dim.z+z];
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+2] = 0;
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+3] = 1;
                    break;
                case GL_BLUE:
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+0] = 0;
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+1] = 0;
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+2] = texture[x*dim.y*dim.z+y*dim.z+z];
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+3] = 1;
                    break;
                case GL_ALPHA:
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+0] = 0;
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+1] = 0;
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+2] = 0;
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+3] = texture[x*dim.y*dim.z+y*dim.z+z];
                    break;
                case GL_RGB:
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+0] = texture[3*(x*dim.y*dim.z+y*dim.z+z)+0];
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+1] = texture[3*(x*dim.y*dim.z+y*dim.z+z)+1];
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+2] = texture[3*(x*dim.y*dim.z+y*dim.z+z)+2];
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+3] = 1;
                    break;
                case GL_RGBA:
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+0] = texture[4*(x*dim.y*dim.z+y*dim.z+z)+0];
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+1] = texture[4*(x*dim.y*dim.z+y*dim.z+z)+1];
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+2] = texture[4*(x*dim.y*dim.z+y*dim.z+z)+2];
                    data[(x*dim.y*dim.z+y*dim.z+z)*4+3] = texture[4*(x*dim.y*dim.z+y*dim.z+z)+3];
                    break;
                }
            }
        }
    }
    return data;
}

//Resizes a given 3D Textur in RGBA to another dimension by trilinear interpolation
GLubyte* TransFunc1DKeysInterpolationFunctionBase::changeTextureDimension(tgt::ivec3 in_dim, tgt::ivec3 out_dim, GLubyte* indata) {
    GLubyte* outdata;
    if (in_dim.x != out_dim.x) {
        outdata = new GLubyte[4 * out_dim.x * in_dim.y * in_dim.z];
        for (int x = 0; x < out_dim.x; ++x) {
            float x_position = static_cast<float>(in_dim.x-1) / static_cast<float>(out_dim.x-1) * static_cast<float>(x);
            int x1 = static_cast<int>(std::floor(x_position));
            float a1 = 1-(x_position-x1);
            int x2 = static_cast<int>(std::ceil(x_position));
            float a2 = 1-a1;
            int x_value = -1;
            if (x1 == x2)
                x_value = x1;
            if (std::abs(a1) < 0.001f)
                x_value = x2;
            if (std::abs(a2) < 0.001f)
                x_value = x1;
            if (x_value >= 0) {
                for (int y = 0; y < in_dim.y; ++y) {
                    for (int z = 0; z < in_dim.z; ++z) {
                        outdata[4*(x*in_dim.y*in_dim.z+y*in_dim.z+z)+0] = indata[4*(x_value*in_dim.y*in_dim.z+y*in_dim.z+z)+0];
                        outdata[4*(x*in_dim.y*in_dim.z+y*in_dim.z+z)+1] = indata[4*(x_value*in_dim.y*in_dim.z+y*in_dim.z+z)+1];
                        outdata[4*(x*in_dim.y*in_dim.z+y*in_dim.z+z)+2] = indata[4*(x_value*in_dim.y*in_dim.z+y*in_dim.z+z)+2];
                        outdata[4*(x*in_dim.y*in_dim.z+y*in_dim.z+z)+3] = indata[4*(x_value*in_dim.y*in_dim.z+y*in_dim.z+z)+3];
                    }
                }
            }
            else {
                for (int y = 0; y < in_dim.y; ++y) {
                    for (int z = 0; z < in_dim.z; ++z) {
                        outdata[4*(x*in_dim.y*in_dim.z+y*in_dim.z+z)+0] = static_cast<GLubyte>(a1*indata[4*(x1*in_dim.y*in_dim.z+y*in_dim.z+z)+0]+a2*indata[4*(x2*in_dim.y*in_dim.z+y*in_dim.z+z)+0]);
                        outdata[4*(x*in_dim.y*in_dim.z+y*in_dim.z+z)+1] = static_cast<GLubyte>(a1*indata[4*(x1*in_dim.y*in_dim.z+y*in_dim.z+z)+1]+a2*indata[4*(x2*in_dim.y*in_dim.z+y*in_dim.z+z)+1]);
                        outdata[4*(x*in_dim.y*in_dim.z+y*in_dim.z+z)+2] = static_cast<GLubyte>(a1*indata[4*(x1*in_dim.y*in_dim.z+y*in_dim.z+z)+2]+a2*indata[4*(x2*in_dim.y*in_dim.z+y*in_dim.z+z)+2]);
                        outdata[4*(x*in_dim.y*in_dim.z+y*in_dim.z+z)+3] = static_cast<GLubyte>(a1*indata[4*(x1*in_dim.y*in_dim.z+y*in_dim.z+z)+3]+a2*indata[4*(x2*in_dim.y*in_dim.z+y*in_dim.z+z)+3]);
                    }
                }
            }
        }
        in_dim.x = out_dim.x;
        delete indata;
        indata = outdata;
    }
    if (in_dim.y != out_dim.y) {
        outdata = new GLubyte[4 * out_dim.x * out_dim.y * in_dim.z];
        for (int y = 0; y < out_dim.y; ++y) {
            float y_position = static_cast<float>(in_dim.y-1) / static_cast<float>(out_dim.y-1) * static_cast<float>(y);
            int y1 = static_cast<int>(std::floor(y_position));
            float a1 = 1-(y_position-y1);
            int y2 = static_cast<int>(std::ceil(y_position));
            float a2 = 1 - a1;
            int y_value = -1;
            if (y1 == y2)
                y_value = y1;
            if (std::abs(a1) < 0.001f)
                y_value = y2;
            if (std::abs(a2) < 0.001f)
                y_value = y1;
            if (y_value >= 0) {
                for (int x = 0; x < out_dim.x; ++x) {
                    for (int z = 0; z < in_dim.z; ++z) {
                        outdata[4*(x*out_dim.y*in_dim.z+y*in_dim.z+z)+0] = indata[4*(x*in_dim.y*in_dim.z+y_value*in_dim.z+z)+0];
                        outdata[4*(x*out_dim.y*in_dim.z+y*in_dim.z+z)+1] = indata[4*(x*in_dim.y*in_dim.z+y_value*in_dim.z+z)+1];
                        outdata[4*(x*out_dim.y*in_dim.z+y*in_dim.z+z)+2] = indata[4*(x*in_dim.y*in_dim.z+y_value*in_dim.z+z)+2];
                        outdata[4*(x*out_dim.y*in_dim.z+y*in_dim.z+z)+3] = indata[4*(x*in_dim.y*in_dim.z+y_value*in_dim.z+z)+3];
                    }
                }
            }
            else {
                for (int x = 0; x < out_dim.x; ++x) {
                    for (int z = 0; z < in_dim.z; ++z) {
                        outdata[4*(x*out_dim.y*in_dim.z+y*in_dim.z+z)+0] = static_cast<GLubyte>(a1*indata[4*(x*in_dim.y*in_dim.z+y1*in_dim.z+z)+0]+a2*indata[4*(x*in_dim.y*in_dim.z+y2*in_dim.z+z)+0]);
                        outdata[4*(x*out_dim.y*in_dim.z+y*in_dim.z+z)+1] = static_cast<GLubyte>(a1*indata[4*(x*in_dim.y*in_dim.z+y1*in_dim.z+z)+1]+a2*indata[4*(x*in_dim.y*in_dim.z+y2*in_dim.z+z)+1]);
                        outdata[4*(x*out_dim.y*in_dim.z+y*in_dim.z+z)+2] = static_cast<GLubyte>(a1*indata[4*(x*in_dim.y*in_dim.z+y1*in_dim.z+z)+2]+a2*indata[4*(x*in_dim.y*in_dim.z+y2*in_dim.z+z)+2]);
                        outdata[4*(x*out_dim.y*in_dim.z+y*in_dim.z+z)+3] = static_cast<GLubyte>(a1*indata[4*(x*in_dim.y*in_dim.z+y1*in_dim.z+z)+3]+a2*indata[4*(x*in_dim.y*in_dim.z+y2*in_dim.z+z)+3]);
                    }
                }
            }
        }
        in_dim.y = out_dim.y;
        delete indata;
        indata = outdata;
    }
    if (in_dim.z != out_dim.z) {
        outdata = new GLubyte[4 * out_dim.x * out_dim.y * out_dim.z];
        for (int z = 0; z < out_dim.z; ++z) {
            float z_position = static_cast<float>(in_dim.z-1) / static_cast<float>(out_dim.z-1) * static_cast<float>(z);
            int z1 = static_cast<int>(std::floor(z_position));
            float a1 = 1-(z_position-z1);
            int z2 = static_cast<int>(std::ceil(z_position));
            float a2 = 1-a1;
            int z_value = -1;
            if (z1 == z2)
                z_value = z1;
            if (std::abs(a1) < 0.001f)
                z_value = z2;
            if (std::abs(a2) < 0.001f)
                z_value = z1;
            if (z_value >= 0) {
                for (int x = 0; x < out_dim.x; ++x) {
                    for (int y = 0; y < out_dim.y; ++y) {
                        outdata[4*(x*out_dim.y*out_dim.z+y*out_dim.z+z)+0] = indata[4*(x*in_dim.y*in_dim.z+y*in_dim.z+z_value)+0];
                        outdata[4*(x*out_dim.y*out_dim.z+y*out_dim.z+z)+1] = indata[4*(x*in_dim.y*in_dim.z+y*in_dim.z+z_value)+1];
                        outdata[4*(x*out_dim.y*out_dim.z+y*out_dim.z+z)+2] = indata[4*(x*in_dim.y*in_dim.z+y*in_dim.z+z_value)+2];
                        outdata[4*(x*out_dim.y*out_dim.z+y*out_dim.z+z)+3] = indata[4*(x*in_dim.y*in_dim.z+y*in_dim.z+z_value)+3];
                    }
                }
            }
            else {
                for (int x = 0; x < out_dim.x; ++x) {
                    for (int y = 0; y < out_dim.y; ++y) {
                        outdata[4*(x*out_dim.y*out_dim.z+y*out_dim.z+z)+0] = static_cast<GLubyte>(a1*indata[4*(x*in_dim.y*in_dim.z+y*in_dim.z+z1)+0]+a2*indata[4*(x*in_dim.y*in_dim.z+y*in_dim.z+z2)+0]);
                        outdata[4*(x*out_dim.y*out_dim.z+y*out_dim.z+z)+1] = static_cast<GLubyte>(a1*indata[4*(x*in_dim.y*in_dim.z+y*in_dim.z+z1)+1]+a2*indata[4*(x*in_dim.y*in_dim.z+y*in_dim.z+z2)+1]);
                        outdata[4*(x*out_dim.y*out_dim.z+y*out_dim.z+z)+2] = static_cast<GLubyte>(a1*indata[4*(x*in_dim.y*in_dim.z+y*in_dim.z+z1)+2]+a2*indata[4*(x*in_dim.y*in_dim.z+y*in_dim.z+z2)+2]);
                        outdata[4*(x*out_dim.y*out_dim.z+y*out_dim.z+z)+3] = static_cast<GLubyte>(a1*indata[4*(x*in_dim.y*in_dim.z+y*in_dim.z+z1)+3]+a2*indata[4*(x*in_dim.y*in_dim.z+y*in_dim.z+z2)+3]);
                    }
                }
            }
        }
        in_dim.y = out_dim.y;
        delete indata;
        indata = outdata;
    }
    return indata;
}

//-------------------------------------------------------------------------------------------------

TransFunc1DKeysInterpolationFunction::TransFunc1DKeysInterpolationFunction() {}

std::string TransFunc1DKeysInterpolationFunction::getGuiName() const {
    return "linear interpolation";
}

std::string TransFunc1DKeysInterpolationFunction::getCategory() const {
    return "default linear (keywise if possible)";
}

TransFunc1DKeys* TransFunc1DKeysInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    if (!startvalue || !endvalue) {
        LERROR("Null pointer passed");
        return 0;
    }

    tgt::ivec3 dimensions_start = startvalue->getDimensions();
    tgt::ivec3 dimensions_end = endvalue->getDimensions();
    tgt::ivec3 dim = tgt::ivec3();
    dim.x = std::max(dimensions_start.x,dimensions_end.x);
    dim.y = std::max(dimensions_start.y,dimensions_end.y);
    dim.z = std::max(dimensions_start.z,dimensions_end.z);

    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys(dim.x);
        func->setAlphaMode(startvalue->getAlphaMode());

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::linearInterpolation(t1.x, t2.x, time),
                                BasicFloatInterpolation::linearInterpolation(t1.y, t2.y, time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
           LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::linearInterpolation(d1.x, d2.x, time),
                        BasicFloatInterpolation::linearInterpolation(d1.y, d2.y, time));

        float gamma1 = startvalue->getGammaValue();
        float gamma2 = endvalue->getGammaValue();
        func->setGammaValue(BasicFloatInterpolation::linearInterpolation(gamma1, gamma2, time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>((1-time)*(*it1)->getColorL().r + time*(*it2)->getColorL().r);
            col.g = static_cast<uint8_t>((1-time)*(*it1)->getColorL().g + time*(*it2)->getColorL().g);
            col.b = static_cast<uint8_t>((1-time)*(*it1)->getColorL().b + time*(*it2)->getColorL().b);
            col.a = static_cast<uint8_t>((1-time)*(*it1)->getColorL().a + time*(*it2)->getColorL().a);
            key->setColorL(col);

            col.r = static_cast<uint8_t>((1-time)*(*it1)->getColorR().r + time*(*it2)->getColorR().r);
            col.g = static_cast<uint8_t>((1-time)*(*it1)->getColorR().g + time*(*it2)->getColorR().g);
            col.b = static_cast<uint8_t>((1-time)*(*it1)->getColorR().b + time*(*it2)->getColorR().b);
            col.a = static_cast<uint8_t>((1-time)*(*it1)->getColorR().a + time*(*it2)->getColorR().a);
            key->setColorR(col);

            key->setIntensity((1-time)*(*it1)->getIntensity() + time*(*it2)->getIntensity());

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    LERROR("No equal key size!");
    return 0;
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysInterpolationFunction::create() const {
    return new TransFunc1DKeysInterpolationFunction();
}

//-------------------------------------------------------------------------------------------------

TransFunc1DKeysStartInterpolationFunction::TransFunc1DKeysStartInterpolationFunction() {}

std::string TransFunc1DKeysStartInterpolationFunction::getGuiName() const {
    return "focus on startvalue";
}

std::string TransFunc1DKeysStartInterpolationFunction::getCategory() const {
    return "boolean";
}

TransFunc1DKeys* TransFunc1DKeysStartInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->clone());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->clone());
}
InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysStartInterpolationFunction::create() const {
    return new TransFunc1DKeysStartInterpolationFunction();
}

TransFunc1DKeysEndInterpolationFunction::TransFunc1DKeysEndInterpolationFunction() {}

std::string TransFunc1DKeysEndInterpolationFunction::getGuiName() const {
    return "focus on endvalue";
}

std::string TransFunc1DKeysEndInterpolationFunction::getCategory() const {
    return "boolean";
}

TransFunc1DKeys* TransFunc1DKeysEndInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    if (time > 0.f)
        return static_cast<TransFunc1DKeys*>(endvalue->clone());
    else
        return static_cast<TransFunc1DKeys*>(startvalue->clone());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysEndInterpolationFunction::create() const {
    return new TransFunc1DKeysEndInterpolationFunction();
}

TransFunc1DKeysStartEndInterpolationFunction::TransFunc1DKeysStartEndInterpolationFunction() {}

std::string TransFunc1DKeysStartEndInterpolationFunction::getGuiName() const {
    return "bisection";
}

std::string TransFunc1DKeysStartEndInterpolationFunction::getCategory() const {
    return "boolean";
}

TransFunc1DKeys* TransFunc1DKeysStartEndInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    if (time < 0.5f)
        return static_cast<TransFunc1DKeys*>(startvalue->clone());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->clone());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysStartEndInterpolationFunction::create() const {
    return new TransFunc1DKeysStartEndInterpolationFunction();
}

TransFunc1DKeysKeyWiseInterpolationFunction::TransFunc1DKeysKeyWiseInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseInterpolationFunction::getGuiName() const {
    return "keywise linear";
}

std::string TransFunc1DKeysKeyWiseInterpolationFunction::getCategory() const {
    return "keywise linear";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold((1-time)*t1.x+time*t2.x,
                           (1-time)*t1.y+time*t2.y);

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain((1-time)*d1.x+time*d2.x,
                        (1-time)*d1.y+time*d2.y);

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>((1-time)*(*it1)->getColorL().r + time*(*it2)->getColorL().r);
            col.g = static_cast<uint8_t>((1-time)*(*it1)->getColorL().g + time*(*it2)->getColorL().g);
            col.b = static_cast<uint8_t>((1-time)*(*it1)->getColorL().b + time*(*it2)->getColorL().b);
            col.a = static_cast<uint8_t>((1-time)*(*it1)->getColorL().a + time*(*it2)->getColorL().a);
            key->setColorL(col);

            col.r = static_cast<uint8_t>((1-time)*(*it1)->getColorR().r + time*(*it2)->getColorR().r);
            col.g = static_cast<uint8_t>((1-time)*(*it1)->getColorR().g + time*(*it2)->getColorR().g);
            col.b = static_cast<uint8_t>((1-time)*(*it1)->getColorR().b + time*(*it2)->getColorR().b);
            col.a = static_cast<uint8_t>((1-time)*(*it1)->getColorR().a + time*(*it2)->getColorR().a);
            key->setColorR(col);

            key->setIntensity((1-time)*(*it1)->getIntensity()+time*(*it2)->getIntensity());

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseInterpolationFunction();
}

TransFunc1DKeysKeyWiseQuadInInterpolationFunction::TransFunc1DKeysKeyWiseQuadInInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseQuadInInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string TransFunc1DKeysKeyWiseQuadInInterpolationFunction::getCategory() const {
    return "keywise quadratic";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseQuadInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::inQuadInterpolation(t1.x,t2.x,time),
                           BasicFloatInterpolation::inQuadInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::inQuadInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::inQuadInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inQuadInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inQuadInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inQuadInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inQuadInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inQuadInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inQuadInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inQuadInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inQuadInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::inQuadInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseQuadInInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseQuadInInterpolationFunction();
}

TransFunc1DKeysKeyWiseQuadOutInterpolationFunction::TransFunc1DKeysKeyWiseQuadOutInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseQuadOutInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string TransFunc1DKeysKeyWiseQuadOutInterpolationFunction::getCategory() const {
    return "keywise quadratic";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseQuadOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::outQuadInterpolation(t1.x,t2.x,time),
                            BasicFloatInterpolation::outQuadInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
                            BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outQuadInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outQuadInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outQuadInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outQuadInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outQuadInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outQuadInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outQuadInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outQuadInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::outQuadInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseQuadOutInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseQuadOutInterpolationFunction();
}

TransFunc1DKeysKeyWiseQuadInOutInterpolationFunction::TransFunc1DKeysKeyWiseQuadInOutInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseQuadInOutInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string TransFunc1DKeysKeyWiseQuadInOutInterpolationFunction::getCategory() const {
    return "keywise quadratic";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseQuadInOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::inOutQuadInterpolation(t1.x,t2.x,time),
                            BasicFloatInterpolation::inOutQuadInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::inOutQuadInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::inOutQuadInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuadInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuadInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuadInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuadInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuadInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuadInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuadInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuadInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::inOutQuadInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseQuadInOutInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseQuadInOutInterpolationFunction();
}

TransFunc1DKeysKeyWiseQuadOutInInterpolationFunction::TransFunc1DKeysKeyWiseQuadOutInInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseQuadOutInInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string TransFunc1DKeysKeyWiseQuadOutInInterpolationFunction::getCategory() const {
    return "keywise quadratic";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseQuadOutInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::outInQuadInterpolation(t1.x,t2.x,time),
                           BasicFloatInterpolation::outInQuadInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::outInQuadInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::outInQuadInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outInQuadInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outInQuadInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outInQuadInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outInQuadInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outInQuadInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outInQuadInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outInQuadInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outInQuadInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::outInQuadInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseQuadOutInInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseQuadOutInInterpolationFunction();
}

TransFunc1DKeysKeyWiseCubicInInterpolationFunction::TransFunc1DKeysKeyWiseCubicInInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseCubicInInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string TransFunc1DKeysKeyWiseCubicInInterpolationFunction::getCategory() const {
    return "keywise cubic";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseCubicInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::inQuadInterpolation(t1.x,t2.x,time),
                            BasicFloatInterpolation::inQuadInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::inQuadInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::inQuadInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inCubicInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inCubicInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inCubicInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inCubicInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inCubicInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inCubicInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inCubicInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inCubicInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::inCubicInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseCubicInInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseCubicInInterpolationFunction();
}

TransFunc1DKeysKeyWiseCubicOutInterpolationFunction::TransFunc1DKeysKeyWiseCubicOutInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseCubicOutInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string TransFunc1DKeysKeyWiseCubicOutInterpolationFunction::getCategory() const {
    return "keywise cubic";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseCubicOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::outCubicInterpolation(t1.x,t2.x,time),
                           BasicFloatInterpolation::outCubicInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::outCubicInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::outCubicInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outCubicInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outCubicInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outCubicInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outCubicInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outCubicInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outCubicInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outCubicInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outCubicInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::outCubicInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseCubicOutInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseCubicOutInterpolationFunction();
}

TransFunc1DKeysKeyWiseCubicInOutInterpolationFunction::TransFunc1DKeysKeyWiseCubicInOutInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseCubicInOutInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string TransFunc1DKeysKeyWiseCubicInOutInterpolationFunction::getCategory() const {
    return "keywise cubic";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseCubicInOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::inOutCubicInterpolation(t1.x,t2.x,time),
                           BasicFloatInterpolation::inOutCubicInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::inOutCubicInterpolation(d1.x,d2.x,time),
                            BasicFloatInterpolation::inOutCubicInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inOutCubicInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inOutCubicInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inOutCubicInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inOutCubicInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inOutCubicInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inOutCubicInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inOutCubicInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inOutCubicInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::inOutCubicInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseCubicInOutInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseCubicInOutInterpolationFunction();
}

TransFunc1DKeysKeyWiseCubicOutInInterpolationFunction::TransFunc1DKeysKeyWiseCubicOutInInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseCubicOutInInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string TransFunc1DKeysKeyWiseCubicOutInInterpolationFunction::getCategory() const {
    return "keywise cubic";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseCubicOutInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::outInCubicInterpolation(t1.x,t2.x,time),
                           BasicFloatInterpolation::outInCubicInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::outInCubicInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::outInCubicInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outInCubicInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outInCubicInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outInCubicInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outInCubicInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outInCubicInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outInCubicInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outInCubicInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outInCubicInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::outInCubicInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseCubicOutInInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseCubicOutInInterpolationFunction();
}

TransFunc1DKeysKeyWiseQuartInInterpolationFunction::TransFunc1DKeysKeyWiseQuartInInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseQuartInInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string TransFunc1DKeysKeyWiseQuartInInterpolationFunction::getCategory() const {
    return "keywise quartetic";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseQuartInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::inQuartInterpolation(t1.x,t2.x,time),
                           BasicFloatInterpolation::inQuartInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::inQuartInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::inQuartInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inQuartInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inQuartInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inQuartInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inQuartInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inQuartInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inQuartInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inQuartInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inQuartInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::inQuartInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseQuartInInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseQuartInInterpolationFunction();
}

TransFunc1DKeysKeyWiseQuartOutInterpolationFunction::TransFunc1DKeysKeyWiseQuartOutInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseQuartOutInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string TransFunc1DKeysKeyWiseQuartOutInterpolationFunction::getCategory() const {
    return "keywise quartetic";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseQuartOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {

    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::outQuartInterpolation(t1.x,t2.x,time),
                            BasicFloatInterpolation::outQuartInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::outQuartInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::outQuartInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outQuartInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outQuartInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outQuartInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outQuartInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outQuartInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outQuartInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outQuartInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outQuartInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::outQuartInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseQuartOutInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseQuartOutInterpolationFunction();
}

TransFunc1DKeysKeyWiseQuartInOutInterpolationFunction::TransFunc1DKeysKeyWiseQuartInOutInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseQuartInOutInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string TransFunc1DKeysKeyWiseQuartInOutInterpolationFunction::getCategory() const {
    return "keywise quartetic";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseQuartInOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::inOutQuartInterpolation(t1.x,t2.x,time),
                            BasicFloatInterpolation::inOutQuartInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::inOutQuartInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::inOutQuartInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuartInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuartInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuartInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuartInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuartInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuartInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuartInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuartInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::inOutQuartInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseQuartInOutInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseQuartInOutInterpolationFunction();
}

TransFunc1DKeysKeyWiseQuartOutInInterpolationFunction::TransFunc1DKeysKeyWiseQuartOutInInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseQuartOutInInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string TransFunc1DKeysKeyWiseQuartOutInInterpolationFunction::getCategory() const {
    return "keywise quartetic";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseQuartOutInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::outInQuartInterpolation(t1.x,t2.x,time),
                            BasicFloatInterpolation::outInQuartInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::outInQuartInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::outInQuartInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outInQuartInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outInQuartInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outInQuartInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outInQuartInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outInQuartInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outInQuartInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outInQuartInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outInQuartInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::outInQuartInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseQuartOutInInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseQuartOutInInterpolationFunction();
}

TransFunc1DKeysKeyWiseQuintInInterpolationFunction::TransFunc1DKeysKeyWiseQuintInInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseQuintInInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string TransFunc1DKeysKeyWiseQuintInInterpolationFunction::getCategory() const {
    return "keywise quintic";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseQuintInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {

    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::inQuintInterpolation(t1.x,t2.x,time),
                            BasicFloatInterpolation::inQuintInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::inQuintInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::inQuintInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inQuintInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inQuintInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inQuintInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inQuintInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inQuintInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inQuintInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inQuintInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inQuintInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::inQuintInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseQuintInInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseQuintInInterpolationFunction();
}

TransFunc1DKeysKeyWiseQuintOutInterpolationFunction::TransFunc1DKeysKeyWiseQuintOutInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseQuintOutInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string TransFunc1DKeysKeyWiseQuintOutInterpolationFunction::getCategory() const {
    return "keywise quintic";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseQuintOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {

    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::outQuintInterpolation(t1.x,t2.x,time),
                            BasicFloatInterpolation::outQuintInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::outQuintInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::outQuintInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outQuintInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outQuintInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outQuintInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outQuintInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outQuintInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outQuintInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outQuintInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outQuintInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::outQuintInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseQuintOutInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseQuintOutInterpolationFunction();
}

TransFunc1DKeysKeyWiseQuintInOutInterpolationFunction::TransFunc1DKeysKeyWiseQuintInOutInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseQuintInOutInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string TransFunc1DKeysKeyWiseQuintInOutInterpolationFunction::getCategory() const {
    return "keywise quintic";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseQuintInOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {

    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::inOutQuintInterpolation(t1.x,t2.x,time),
                            BasicFloatInterpolation::inOutQuintInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::inOutQuintInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::inOutQuintInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuintInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuintInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuintInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuintInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuintInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuintInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuintInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inOutQuintInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::inOutQuintInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseQuintInOutInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseQuintInOutInterpolationFunction();
}

TransFunc1DKeysKeyWiseQuintOutInInterpolationFunction::TransFunc1DKeysKeyWiseQuintOutInInterpolationFunction() {
}

std::string TransFunc1DKeysKeyWiseQuintOutInInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string TransFunc1DKeysKeyWiseQuintOutInInterpolationFunction::getCategory() const {
    return "keywise quintic";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseQuintOutInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {

    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::outInQuintInterpolation(t1.x,t2.x,time),
                            BasicFloatInterpolation::outInQuintInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::outInQuintInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::outInQuintInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outInQuintInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outInQuintInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outInQuintInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outInQuintInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outInQuintInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outInQuintInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outInQuintInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outInQuintInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::outInQuintInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseQuintOutInInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseQuintOutInInterpolationFunction();
}

TransFunc1DKeysKeyWiseSineInInterpolationFunction::TransFunc1DKeysKeyWiseSineInInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseSineInInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string TransFunc1DKeysKeyWiseSineInInterpolationFunction::getCategory() const {
    return "keywise sineousidal";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseSineInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {

    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::inSineInterpolation(t1.x,t2.x,time),
                            BasicFloatInterpolation::inSineInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::inSineInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::inSineInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inSineInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inSineInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inSineInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inSineInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inSineInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inSineInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inSineInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inSineInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::inSineInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseSineInInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseSineInInterpolationFunction();
}

TransFunc1DKeysKeyWiseSineOutInterpolationFunction::TransFunc1DKeysKeyWiseSineOutInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseSineOutInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string TransFunc1DKeysKeyWiseSineOutInterpolationFunction::getCategory() const {
    return "keywise sineousidal";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseSineOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::outSineInterpolation(t1.x,t2.x,time),
                            BasicFloatInterpolation::outSineInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::outSineInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::outSineInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outSineInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outSineInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outSineInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outSineInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outSineInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outSineInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outSineInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outSineInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::outSineInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseSineOutInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseSineOutInterpolationFunction();
}

TransFunc1DKeysKeyWiseSineInOutInterpolationFunction::TransFunc1DKeysKeyWiseSineInOutInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseSineInOutInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string TransFunc1DKeysKeyWiseSineInOutInterpolationFunction::getCategory() const {
    return "keywise sineousidal";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseSineInOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::inOutSineInterpolation(t1.x,t2.x,time),
                            BasicFloatInterpolation::inOutSineInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::inOutSineInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::inOutSineInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inOutSineInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inOutSineInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inOutSineInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inOutSineInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inOutSineInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inOutSineInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inOutSineInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inOutSineInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::inOutSineInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseSineInOutInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseSineInOutInterpolationFunction();
}

TransFunc1DKeysKeyWiseSineOutInInterpolationFunction::TransFunc1DKeysKeyWiseSineOutInInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseSineOutInInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string TransFunc1DKeysKeyWiseSineOutInInterpolationFunction::getCategory() const {
    return "keywise sineousidal";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseSineOutInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::outInSineInterpolation(t1.x,t2.x,time),
                            BasicFloatInterpolation::outInSineInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");


        func->setDomain(BasicFloatInterpolation::outInSineInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::outInSineInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outInSineInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outInSineInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outInSineInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outInSineInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outInSineInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outInSineInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outInSineInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outInSineInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::outInSineInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseSineOutInInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseSineOutInInterpolationFunction();
}

TransFunc1DKeysKeyWiseExponentInInterpolationFunction::TransFunc1DKeysKeyWiseExponentInInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseExponentInInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string TransFunc1DKeysKeyWiseExponentInInterpolationFunction::getCategory() const {
    return "keywise exponential";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseExponentInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::inExponentInterpolation(t1.x,t2.x,time),
                            BasicFloatInterpolation::inExponentInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if (d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::inExponentInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::inExponentInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inExponentInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inExponentInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inExponentInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inExponentInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inExponentInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inExponentInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inExponentInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inExponentInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::inExponentInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseExponentInInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseExponentInInterpolationFunction();
}

TransFunc1DKeysKeyWiseExponentOutInterpolationFunction::TransFunc1DKeysKeyWiseExponentOutInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseExponentOutInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string TransFunc1DKeysKeyWiseExponentOutInterpolationFunction::getCategory() const {
    return "keywise exponential";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseExponentOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::outExponentInterpolation(t1.x,t2.x,time),
                            BasicFloatInterpolation::outExponentInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if(d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::outExponentInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::outExponentInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outExponentInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outExponentInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outExponentInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outExponentInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outExponentInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outExponentInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outExponentInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outExponentInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::outExponentInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseExponentOutInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseExponentOutInterpolationFunction();
}

TransFunc1DKeysKeyWiseExponentInOutInterpolationFunction::TransFunc1DKeysKeyWiseExponentInOutInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseExponentInOutInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string TransFunc1DKeysKeyWiseExponentInOutInterpolationFunction::getCategory() const {
    return "keywise exponential";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseExponentInOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::inOutExponentInterpolation(t1.x,t2.x,time),
                            BasicFloatInterpolation::inOutExponentInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if(d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::inOutExponentInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::inOutExponentInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inOutExponentInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inOutExponentInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inOutExponentInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inOutExponentInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inOutExponentInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inOutExponentInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inOutExponentInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inOutExponentInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::inOutExponentInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }

    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseExponentInOutInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseExponentInOutInterpolationFunction();
}

TransFunc1DKeysKeyWiseExponentOutInInterpolationFunction::TransFunc1DKeysKeyWiseExponentOutInInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseExponentOutInInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string TransFunc1DKeysKeyWiseExponentOutInInterpolationFunction::getCategory() const {
    return "keywise exponential";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseExponentOutInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {

    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::outInExponentInterpolation(t1.x,t2.x,time),
                            BasicFloatInterpolation::outInExponentInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if(d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::outInExponentInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::outInExponentInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outInExponentInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outInExponentInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outInExponentInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outInExponentInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outInExponentInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outInExponentInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outInExponentInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outInExponentInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::outInExponentInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseExponentOutInInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseExponentOutInInterpolationFunction();
}

TransFunc1DKeysKeyWiseCircInInterpolationFunction::TransFunc1DKeysKeyWiseCircInInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseCircInInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string TransFunc1DKeysKeyWiseCircInInterpolationFunction::getCategory() const {
    return "keywise circular";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseCircInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {

    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::inCircInterpolation(t1.x,t2.x,time),
                            BasicFloatInterpolation::inCircInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if(d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::inCircInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::inCircInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inCircInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inCircInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inCircInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inCircInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inCircInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inCircInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inCircInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inCircInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::inCircInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseCircInInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseCircInInterpolationFunction();
}

TransFunc1DKeysKeyWiseCircOutInterpolationFunction::TransFunc1DKeysKeyWiseCircOutInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseCircOutInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string TransFunc1DKeysKeyWiseCircOutInterpolationFunction::getCategory() const {
    return "keywise circular";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseCircOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::outCircInterpolation(t1.x,t2.x,time),
                            BasicFloatInterpolation::outCircInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if(d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::outCircInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::outCircInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outCircInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outCircInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outCircInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outCircInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outCircInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outCircInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outCircInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outCircInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::outCircInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseCircOutInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseCircOutInterpolationFunction();
}

TransFunc1DKeysKeyWiseCircInOutInterpolationFunction::TransFunc1DKeysKeyWiseCircInOutInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseCircInOutInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string TransFunc1DKeysKeyWiseCircInOutInterpolationFunction::getCategory() const {
    return "keywise circular";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseCircInOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {

    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::inOutCircInterpolation(t1.x,t2.x,time),
                            BasicFloatInterpolation::inOutCircInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if(d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::inOutCircInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::inOutCircInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inOutCircInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inOutCircInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inOutCircInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inOutCircInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::inOutCircInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::inOutCircInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::inOutCircInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::inOutCircInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::inOutCircInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }
    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseCircInOutInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseCircInOutInterpolationFunction();
}

TransFunc1DKeysKeyWiseCircOutInInterpolationFunction::TransFunc1DKeysKeyWiseCircOutInInterpolationFunction() {}

std::string TransFunc1DKeysKeyWiseCircOutInInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string TransFunc1DKeysKeyWiseCircOutInInterpolationFunction::getCategory() const {
    return "keywise circular";
}

TransFunc1DKeys* TransFunc1DKeysKeyWiseCircOutInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
    std::vector<TransFuncMappingKey*> keys1 = startvalue->getKeys();
    std::vector<TransFuncMappingKey*> keys2 = endvalue->getKeys();
    if (keys1.size() == keys2.size()) {
        TransFunc1DKeys* func = new TransFunc1DKeys();

        tgt::vec2 t1 = startvalue->getThreshold();
        tgt::vec2 t2 = endvalue->getThreshold();

        func->setThreshold(BasicFloatInterpolation::outInCircInterpolation(t1.x,t2.x,time),
                            BasicFloatInterpolation::outInCircInterpolation(t1.y,t2.y,time));

        tgt::vec2 d1 = startvalue->getDomain();
        tgt::vec2 d2 = endvalue->getDomain();
        if(d1 != d2)
            LDEBUG("Transfer functions have different domains...interpolation is (probably) incorrect.");

        func->setDomain(BasicFloatInterpolation::outInCircInterpolation(d1.x,d2.x,time),
                        BasicFloatInterpolation::outInCircInterpolation(d1.y,d2.y,time));

        func->clearKeys();
        std::vector<TransFuncMappingKey*>::iterator it1 = keys1.begin();
        std::vector<TransFuncMappingKey*>::iterator it2 = keys2.begin();
        while ((it1 != keys1.end()) && (it2 != keys2.end())) {
            tgt::col4 col = tgt::col4();
            TransFuncMappingKey* key = new TransFuncMappingKey(0, col);
            key->setSplit((*it1)->isSplit()||(*it2)->isSplit(), true);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outInCircInterpolation((*it1)->getColorL().r,(*it2)->getColorL().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outInCircInterpolation((*it1)->getColorL().g,(*it2)->getColorL().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outInCircInterpolation((*it1)->getColorL().b,(*it2)->getColorL().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outInCircInterpolation((*it1)->getColorL().a,(*it2)->getColorL().a,time));
            key->setColorL(col);

            col.r = static_cast<uint8_t>(BasicFloatInterpolation::outInCircInterpolation((*it1)->getColorR().r,(*it2)->getColorR().r,time));
            col.g = static_cast<uint8_t>(BasicFloatInterpolation::outInCircInterpolation((*it1)->getColorR().g,(*it2)->getColorR().g,time));
            col.b = static_cast<uint8_t>(BasicFloatInterpolation::outInCircInterpolation((*it1)->getColorR().b,(*it2)->getColorR().b,time));
            col.a = static_cast<uint8_t>(BasicFloatInterpolation::outInCircInterpolation((*it1)->getColorR().a,(*it2)->getColorR().a,time));
            key->setColorR(col);

            key->setIntensity(BasicFloatInterpolation::outInCircInterpolation((*it1)->getIntensity(),(*it2)->getIntensity(),time));

            func->addKey(key);

            it1++;
            it2++;
        }
        return func;
    }

    if (time < 1.f)
        return static_cast<TransFunc1DKeys*>(startvalue->create());
    else
        return static_cast<TransFunc1DKeys*>(endvalue->create());
}

InterpolationFunction<TransFunc1DKeys*>* TransFunc1DKeysKeyWiseCircOutInInterpolationFunction::create() const {
    return new TransFunc1DKeysKeyWiseCircOutInInterpolationFunction();
}

//TransFuncTextureLinearInterpolationFunction::TransFuncTextureLinearInterpolationFunction() {}
//
//std::string TransFuncTextureLinearInterpolationFunction::getGuiName() const {
//    return "texturebased linear";
//}
//
//std::string TransFuncTextureLinearInterpolationFunction::getCategory() const {
//    return "texture linear";
//}
//
//TransFunc1DKeys* TransFuncTextureLinearInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::linearInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::linearInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::linearInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureLinearInterpolationFunction::create() const {
//    return new TransFuncTextureLinearInterpolationFunction();
//}
//
//TransFuncTextureQuadInInterpolationFunction::TransFuncTextureQuadInInterpolationFunction() {}
//
//std::string TransFuncTextureQuadInInterpolationFunction::getGuiName() const {
//    return "easing in";
//}
//
//std::string TransFuncTextureQuadInInterpolationFunction::getCategory() const {
//    return "texture quadratic";
//}
//
//TransFunc1DKeys* TransFuncTextureQuadInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::inQuadInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureQuadInInterpolationFunction::create() const {
//    return new TransFuncTextureQuadInInterpolationFunction();
//}
//
//TransFuncTextureQuadOutInterpolationFunction::TransFuncTextureQuadOutInterpolationFunction() {}
//
//std::string TransFuncTextureQuadOutInterpolationFunction::getGuiName() const {
//    return "easing out";
//}
//
//std::string TransFuncTextureQuadOutInterpolationFunction::getCategory() const {
//    return "texture quadratic";
//}
//
//TransFunc1DKeys* TransFuncTextureQuadOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::outQuadInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureQuadOutInterpolationFunction::create() const {
//    return new TransFuncTextureQuadOutInterpolationFunction();
//}
//
//TransFuncTextureQuadInOutInterpolationFunction::TransFuncTextureQuadInOutInterpolationFunction() {}
//
//std::string TransFuncTextureQuadInOutInterpolationFunction::getGuiName() const {
//    return "first easing in, then easing out";
//}
//
//std::string TransFuncTextureQuadInOutInterpolationFunction::getCategory() const {
//    return "texture quadratic";
//}
//
//TransFunc1DKeys* TransFuncTextureQuadInOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::inOutQuadInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureQuadInOutInterpolationFunction::create() const {
//    return new TransFuncTextureQuadInOutInterpolationFunction();
//}
//
//TransFuncTextureQuadOutInInterpolationFunction::TransFuncTextureQuadOutInInterpolationFunction() {}
//
//std::string TransFuncTextureQuadOutInInterpolationFunction::getGuiName() const {
//    return "first easing out, then easing in";
//}
//
//std::string TransFuncTextureQuadOutInInterpolationFunction::getCategory() const {
//    return "texture quadratic";
//}
//
//TransFunc1DKeys* TransFuncTextureQuadOutInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::outInQuadInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureQuadOutInInterpolationFunction::create() const {
//    return new TransFuncTextureQuadOutInInterpolationFunction();
//}
//
//TransFuncTextureCubicInInterpolationFunction::TransFuncTextureCubicInInterpolationFunction() {}
//
//std::string TransFuncTextureCubicInInterpolationFunction::getGuiName() const {
//    return "easing in";
//}
//
//std::string TransFuncTextureCubicInInterpolationFunction::getCategory() const {
//    return "texture cubic";
//}
//
//TransFunc1DKeys* TransFuncTextureCubicInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::inCubicInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureCubicInInterpolationFunction::create() const {
//    return new TransFuncTextureCubicInInterpolationFunction();
//}
//
//TransFuncTextureCubicOutInterpolationFunction::TransFuncTextureCubicOutInterpolationFunction() {}
//
//std::string TransFuncTextureCubicOutInterpolationFunction::getGuiName() const {
//    return "easing out";
//}
//
//std::string TransFuncTextureCubicOutInterpolationFunction::getCategory() const {
//    return "texture cubic";
//}
//
//TransFunc1DKeys* TransFuncTextureCubicOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::outCubicInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureCubicOutInterpolationFunction::create() const {
//    return new TransFuncTextureCubicOutInterpolationFunction();
//}
//
//TransFuncTextureCubicInOutInterpolationFunction::TransFuncTextureCubicInOutInterpolationFunction() {}
//
//std::string TransFuncTextureCubicInOutInterpolationFunction::getGuiName() const {
//    return "first easing in, then easing out";
//}
//
//std::string TransFuncTextureCubicInOutInterpolationFunction::getCategory() const {
//    return "texture cubic";
//}
//
//TransFunc1DKeys* TransFuncTextureCubicInOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::inOutCubicInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureCubicInOutInterpolationFunction::create() const {
//    return new TransFuncTextureCubicInOutInterpolationFunction();
//}
//
//TransFuncTextureCubicOutInInterpolationFunction::TransFuncTextureCubicOutInInterpolationFunction() {}
//
//std::string TransFuncTextureCubicOutInInterpolationFunction::getGuiName() const {
//    return "first easing out, then easing in";
//}
//
//std::string TransFuncTextureCubicOutInInterpolationFunction::getCategory() const {
//    return "texture cubic";
//}
//
//TransFunc1DKeys* TransFuncTextureCubicOutInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::outInCubicInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureCubicOutInInterpolationFunction::create() const {
//    return new TransFuncTextureCubicOutInInterpolationFunction();
//}
//
//TransFuncTextureQuartInInterpolationFunction::TransFuncTextureQuartInInterpolationFunction() {}
//
//std::string TransFuncTextureQuartInInterpolationFunction::getGuiName() const {
//    return "easing in";
//}
//
//std::string TransFuncTextureQuartInInterpolationFunction::getCategory() const {
//    return "texture quartetic";
//}
//
//TransFunc1DKeys* TransFuncTextureQuartInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::inQuartInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureQuartInInterpolationFunction::create() const {
//    return new TransFuncTextureQuartInInterpolationFunction();
//}
//
//TransFuncTextureQuartOutInterpolationFunction::TransFuncTextureQuartOutInterpolationFunction() {}
//
//std::string TransFuncTextureQuartOutInterpolationFunction::getGuiName() const {
//    return "easing out";
//}
//
//std::string TransFuncTextureQuartOutInterpolationFunction::getCategory() const {
//    return "texture quartetic";
//}
//
//TransFunc1DKeys* TransFuncTextureQuartOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::outQuartInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureQuartOutInterpolationFunction::create() const {
//    return new TransFuncTextureQuartOutInterpolationFunction();
//}
//
//TransFuncTextureQuartInOutInterpolationFunction::TransFuncTextureQuartInOutInterpolationFunction() {}
//
//std::string TransFuncTextureQuartInOutInterpolationFunction::getGuiName() const {
//    return "first easing in, then easing out";
//}
//
//std::string TransFuncTextureQuartInOutInterpolationFunction::getCategory() const {
//    return "texture quartetic";
//}
//
//TransFunc1DKeys* TransFuncTextureQuartInOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::inOutQuartInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureQuartInOutInterpolationFunction::create() const {
//    return new TransFuncTextureQuartInOutInterpolationFunction();
//}
//
//TransFuncTextureQuartOutInInterpolationFunction::TransFuncTextureQuartOutInInterpolationFunction() {}
//
//std::string TransFuncTextureQuartOutInInterpolationFunction::getGuiName() const {
//    return "first easing out, then easing in";
//}
//
//std::string TransFuncTextureQuartOutInInterpolationFunction::getCategory() const {
//    return "texture quartetic";
//}
//
//TransFunc1DKeys* TransFuncTextureQuartOutInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::outInQuartInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureQuartOutInInterpolationFunction::create() const {
//    return new TransFuncTextureQuartOutInInterpolationFunction();
//}
//
//TransFuncTextureQuintInInterpolationFunction::TransFuncTextureQuintInInterpolationFunction() {}
//
//std::string TransFuncTextureQuintInInterpolationFunction::getGuiName() const {
//    return "easing in";
//}
//
//std::string TransFuncTextureQuintInInterpolationFunction::getCategory() const {
//    return "texture quintic";
//}
//
//TransFunc1DKeys* TransFuncTextureQuintInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::inQuintInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureQuintInInterpolationFunction::create() const {
//    return new TransFuncTextureQuintInInterpolationFunction();
//}
//
//TransFuncTextureQuintOutInterpolationFunction::TransFuncTextureQuintOutInterpolationFunction() {}
//
//std::string TransFuncTextureQuintOutInterpolationFunction::getGuiName() const {
//    return "easing out";
//}
//
//std::string TransFuncTextureQuintOutInterpolationFunction::getCategory() const {
//    return "texture quintic";
//}
//
//TransFunc1DKeys* TransFuncTextureQuintOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::outQuintInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureQuintOutInterpolationFunction::create() const {
//    return new TransFuncTextureQuintOutInterpolationFunction();
//}
//
//TransFuncTextureQuintInOutInterpolationFunction::TransFuncTextureQuintInOutInterpolationFunction() {}
//
//std::string TransFuncTextureQuintInOutInterpolationFunction::getGuiName() const {
//    return "first easing in, then easing out";
//}
//
//std::string TransFuncTextureQuintInOutInterpolationFunction::getCategory() const {
//    return "texture quintic";
//}
//
//TransFunc1DKeys* TransFuncTextureQuintInOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::inOutQuintInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureQuintInOutInterpolationFunction::create() const {
//    return new TransFuncTextureQuintInOutInterpolationFunction();
//}
//
//TransFuncTextureQuintOutInInterpolationFunction::TransFuncTextureQuintOutInInterpolationFunction() {}
//
//std::string TransFuncTextureQuintOutInInterpolationFunction::getGuiName() const {
//    return "first easing out, then easing in";
//}
//
//std::string TransFuncTextureQuintOutInInterpolationFunction::getCategory() const {
//    return "texture quintic";
//}
//
//TransFunc1DKeys* TransFuncTextureQuintOutInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::outInQuintInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureQuintOutInInterpolationFunction::create() const {
//    return new TransFuncTextureQuintOutInInterpolationFunction();
//}
//
//TransFuncTextureSineInInterpolationFunction::TransFuncTextureSineInInterpolationFunction() {}
//
//std::string TransFuncTextureSineInInterpolationFunction::getGuiName() const {
//    return "easing in";
//}
//
//std::string TransFuncTextureSineInInterpolationFunction::getCategory() const {
//    return "texture sineousidal";
//}
//
//TransFunc1DKeys* TransFuncTextureSineInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::inSineInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureSineInInterpolationFunction::create() const {
//    return new TransFuncTextureSineInInterpolationFunction();
//}
//
//TransFuncTextureSineOutInterpolationFunction::TransFuncTextureSineOutInterpolationFunction() {}
//
//std::string TransFuncTextureSineOutInterpolationFunction::getGuiName() const {
//    return "easing out";
//}
//
//std::string TransFuncTextureSineOutInterpolationFunction::getCategory() const {
//    return "texture sineousidal";
//}
//
//TransFunc1DKeys* TransFuncTextureSineOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::outSineInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureSineOutInterpolationFunction::create() const {
//    return new TransFuncTextureSineOutInterpolationFunction();
//}
//
//TransFuncTextureSineInOutInterpolationFunction::TransFuncTextureSineInOutInterpolationFunction() {}
//
//std::string TransFuncTextureSineInOutInterpolationFunction::getGuiName() const {
//    return "first easing in, then easing out";
//}
//
//std::string TransFuncTextureSineInOutInterpolationFunction::getCategory() const {
//    return "texture sineousidal";
//}
//TransFunc1DKeys* TransFuncTextureSineInOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::inOutSineInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureSineInOutInterpolationFunction::create() const {
//    return new TransFuncTextureSineInOutInterpolationFunction();
//}
//
//TransFuncTextureSineOutInInterpolationFunction::TransFuncTextureSineOutInInterpolationFunction() {}
//
//std::string TransFuncTextureSineOutInInterpolationFunction::getGuiName() const {
//    return "first easing out, then easing in";
//}
//
//std::string TransFuncTextureSineOutInInterpolationFunction::getCategory() const {
//    return "texture sineousidal";
//}
//
//TransFunc1DKeys* TransFuncTextureSineOutInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::outInSineInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureSineOutInInterpolationFunction::create() const {
//    return new TransFuncTextureSineOutInInterpolationFunction();
//}
//
//TransFuncTextureExponentInInterpolationFunction::TransFuncTextureExponentInInterpolationFunction() {}
//
//std::string TransFuncTextureExponentInInterpolationFunction::getGuiName() const {
//    return "easing in";
//}
//
//std::string TransFuncTextureExponentInInterpolationFunction::getCategory() const {
//    return "texture exponential";
//}
//
//TransFunc1DKeys* TransFuncTextureExponentInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::inExponentInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureExponentInInterpolationFunction::create() const {
//    return new TransFuncTextureExponentInInterpolationFunction();
//}
//
//TransFuncTextureExponentOutInterpolationFunction::TransFuncTextureExponentOutInterpolationFunction() {}
//
//std::string TransFuncTextureExponentOutInterpolationFunction::getGuiName() const {
//    return "easing out";
//}
//
//std::string TransFuncTextureExponentOutInterpolationFunction::getCategory() const {
//    return "texture exponential";
//}
//
//TransFunc1DKeys* TransFuncTextureExponentOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::outExponentInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureExponentOutInterpolationFunction::create() const {
//    return new TransFuncTextureExponentOutInterpolationFunction();
//}
//
//TransFuncTextureExponentInOutInterpolationFunction::TransFuncTextureExponentInOutInterpolationFunction() {}
//
//std::string TransFuncTextureExponentInOutInterpolationFunction::getGuiName() const {
//    return "first easing in, then easing out";
//}
//
//std::string TransFuncTextureExponentInOutInterpolationFunction::getCategory() const {
//    return "texture exponential";
//}
//
//TransFunc1DKeys* TransFuncTextureExponentInOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::inOutExponentInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureExponentInOutInterpolationFunction::create() const {
//    return new TransFuncTextureExponentInOutInterpolationFunction();
//}
//
//TransFuncTextureExponentOutInInterpolationFunction::TransFuncTextureExponentOutInInterpolationFunction() {}
//
//std::string TransFuncTextureExponentOutInInterpolationFunction::getGuiName() const {
//    return "first easing out, then easing in";
//}
//
//std::string TransFuncTextureExponentOutInInterpolationFunction::getCategory() const {
//    return "texture exponential";
//}
//
//TransFunc1DKeys* TransFuncTextureExponentOutInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::outInExponentInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureExponentOutInInterpolationFunction::create() const {
//    return new TransFuncTextureExponentOutInInterpolationFunction();
//}
//
//TransFuncTextureCircInInterpolationFunction::TransFuncTextureCircInInterpolationFunction() {}
//
//std::string TransFuncTextureCircInInterpolationFunction::getGuiName() const {
//    return "easing in";
//}
//
//std::string TransFuncTextureCircInInterpolationFunction::getCategory() const {
//    return "texture circular";
//}
//
//TransFunc1DKeys* TransFuncTextureCircInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::inCircInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureCircInInterpolationFunction::create() const {
//    return new TransFuncTextureCircInInterpolationFunction();
//}
//
//TransFuncTextureCircOutInterpolationFunction::TransFuncTextureCircOutInterpolationFunction() {}
//
//std::string TransFuncTextureCircOutInterpolationFunction::getGuiName() const {
//    return "easing out";
//}
//
//std::string TransFuncTextureCircOutInterpolationFunction::getCategory() const {
//    return "texture circular";
//}
//
//TransFunc1DKeys* TransFuncTextureCircOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::outCircInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureCircOutInterpolationFunction::create() const {
//    return new TransFuncTextureCircOutInterpolationFunction();
//}
//
//TransFuncTextureCircInOutInterpolationFunction::TransFuncTextureCircInOutInterpolationFunction() {}
//
//std::string TransFuncTextureCircInOutInterpolationFunction::getGuiName() const {
//    return "first easing in, then easing out";
//}
//
//std::string TransFuncTextureCircInOutInterpolationFunction::getCategory() const {
//    return "texture circular";
//}
//
//TransFunc1DKeys* TransFuncTextureCircInOutInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::inOutCircInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureCircInOutInterpolationFunction::create() const {
//    return new TransFuncTextureCircInOutInterpolationFunction();
//}
//
//TransFuncTextureCircOutInInterpolationFunction::TransFuncTextureCircOutInInterpolationFunction() {}
//
//std::string TransFuncTextureCircOutInInterpolationFunction::getGuiName() const {
//    return "first easing out, then easing in";
//}
//
//std::string TransFuncTextureCircOutInInterpolationFunction::getCategory() const {
//    return "texture circular";
//}
//
//TransFunc1DKeys* TransFuncTextureCircOutInInterpolationFunction::interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const {
//    float a2 = BasicFloatInterpolation::outInCircInterpolation(0, 1, time);
//    float a1 = 1-a2;
//
//    // new dimensions: maxima of each dimension (x,y,z)
//    tgt::ivec3 dimensions_start = startvalue->getDimensions();
//    tgt::ivec3 dimensions_end = endvalue->getDimensions();
//    tgt::ivec3 dim = tgt::ivec3();
//    dim.x = std::max(dimensions_start.x,dimensions_end.x);
//    dim.y = std::max(dimensions_start.y,dimensions_end.y);
//    dim.z = std::max(dimensions_start.z,dimensions_end.z);
//
//    GLubyte* texture1 = startvalue->getPixelData();
//    GLubyte* texture2 = endvalue->getPixelData();
//    texture1 = convertTextureToRGBA(startvalue->getDimensions(), texture1 , startvalue->getFormat());
//    texture2 = convertTextureToRGBA(endvalue->getDimensions(), texture2 , endvalue->getFormat());
//    texture1 = changeTextureDimension(startvalue->getDimensions(), dim, texture1);
//    texture2 = changeTextureDimension(startvalue->getDimensions(), dim, texture2);
//
//    TransFunc1DKeys* func = new TransFunc1DKeys(dim.x,dim.y,dim.z);
//
//    tgt::vec2 d1 = startvalue->getDomain();
//    tgt::vec2 d2 = endvalue->getDomain();
//    func->setDomain(BasicFloatInterpolation::outQuadInterpolation(d1.x,d2.x,time),
//                    BasicFloatInterpolation::outQuadInterpolation(d1.y,d2.y,time));
//
//    GLubyte* texture = new GLubyte[4*dim.x*dim.y*dim.z];
//    for (int x = 0; x < dim.x; ++x) {
//        for (int y = 0; y < dim.y; ++y) {
//            for (int z = 0; z < dim.z; ++z) {
//                for (int i = 0; i < 4; ++i) {
//                    float f = (a1*texture1[4*(x*dim.y*dim.z+y*dim.z+z)+i]+a2*texture2[4*(x*dim.y*dim.z+y*dim.z+z)+i]);
//                    GLubyte b = static_cast<GLubyte>(f);
//                    texture[4*(x*dim.y*dim.z+y*dim.z+z)+i] = b;
//                }
//            }
//        }
//    }
//    func->setPixelData(texture);
//    return func;
//}
//
//InterpolationFunction<TransFunc1DKeys*>* TransFuncTextureCircOutInInterpolationFunction::create() const {
//    return new TransFuncTextureCircOutInInterpolationFunction();
//}

} // namespace voreen
