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

#include "voreen/core/animation/interpolation/transfunc2dprimitivesinterpolationfunctions.h"

namespace voreen {

GLubyte* TransFunc2DPrimitivesInterpolationFunctionBase::convertTextureToRGBA(tgt::ivec3 dim, GLubyte* texture, GLuint inputformat) {
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
GLubyte* TransFunc2DPrimitivesInterpolationFunctionBase::changeTextureDimension(tgt::ivec3 in_dim, tgt::ivec3 out_dim, GLubyte* indata) {
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

TransFunc2DPrimitivesStartInterpolationFunction::TransFunc2DPrimitivesStartInterpolationFunction() {}

std::string TransFunc2DPrimitivesStartInterpolationFunction::getGuiName() const {
    return "focus on startvalue";
}

std::string TransFunc2DPrimitivesStartInterpolationFunction::getCategory() const {
    return "boolean";
}

TransFunc2DPrimitives* TransFunc2DPrimitivesStartInterpolationFunction::interpolate(TransFunc2DPrimitives* startvalue, TransFunc2DPrimitives* endvalue, float time) const {
    if (time < 1.f)
        return static_cast<TransFunc2DPrimitives*>(startvalue->clone());
    else
        return static_cast<TransFunc2DPrimitives*>(endvalue->clone());
}
InterpolationFunction<TransFunc2DPrimitives*>* TransFunc2DPrimitivesStartInterpolationFunction::create() const {
    return new TransFunc2DPrimitivesStartInterpolationFunction();
}

TransFunc2DPrimitivesEndInterpolationFunction::TransFunc2DPrimitivesEndInterpolationFunction() {}

std::string TransFunc2DPrimitivesEndInterpolationFunction::getGuiName() const {
    return "focus on endvalue";
}

std::string TransFunc2DPrimitivesEndInterpolationFunction::getCategory() const {
    return "boolean";
}

TransFunc2DPrimitives* TransFunc2DPrimitivesEndInterpolationFunction::interpolate(TransFunc2DPrimitives* startvalue, TransFunc2DPrimitives* endvalue, float time) const {
    if (time > 0.f)
        return static_cast<TransFunc2DPrimitives*>(endvalue->clone());
    else
        return static_cast<TransFunc2DPrimitives*>(startvalue->clone());
}

InterpolationFunction<TransFunc2DPrimitives*>* TransFunc2DPrimitivesEndInterpolationFunction::create() const {
    return new TransFunc2DPrimitivesEndInterpolationFunction();
}

TransFunc2DPrimitivesStartEndInterpolationFunction::TransFunc2DPrimitivesStartEndInterpolationFunction() {}

std::string TransFunc2DPrimitivesStartEndInterpolationFunction::getGuiName() const {
    return "bisection";
}

std::string TransFunc2DPrimitivesStartEndInterpolationFunction::getCategory() const {
    return "boolean";
}

TransFunc2DPrimitives* TransFunc2DPrimitivesStartEndInterpolationFunction::interpolate(TransFunc2DPrimitives* startvalue, TransFunc2DPrimitives* endvalue, float time) const {
    if (time < 0.5f)
        return static_cast<TransFunc2DPrimitives*>(startvalue->clone());
    else
        return static_cast<TransFunc2DPrimitives*>(endvalue->clone());
}

InterpolationFunction<TransFunc2DPrimitives*>* TransFunc2DPrimitivesStartEndInterpolationFunction::create() const {
    return new TransFunc2DPrimitivesStartEndInterpolationFunction();
}

} // namespace voreen
