/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2018 University of Muenster, Germany,           *
 * Department of Computer Science.                                    *
 *                                                                    *
 * This file is part of the tgt library. This library is free         *
 * software; you can redistribute it and/or modify it under the terms *
 * of the GNU Lesser General Public License version 2.1 as published  *
 * by the Free Software Foundation.                                   *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       *
 * GNU Lesser General Public License for more details.                *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License in the file "LICENSE.txt" along with this library.         *
 * If not, see <http://www.gnu.org/licenses/>.                        *
 *                                                                    *
 **********************************************************************/

#include "tgt/tgt_gl.h"

#include "tgt/logmanager.h"
#include "tgt/glcontextmanager.h"

#include <math.h>

namespace tgt {

static void printGlDebugLog(int line, const char* file) {
#ifdef TGT_DEBUG
    if(GpuCaps.isExtensionSupported("KHR_debug")) {
        std::ostringstream loggerCat, locationInfo;
        loggerCat << "gl-error:" << file << ':' << line;
        locationInfo << " File: " << file << "@" << line;

        while(true) {
            GLint msgLength;
            glGetIntegerv(GL_DEBUG_NEXT_LOGGED_MESSAGE_LENGTH, &msgLength);
            if(msgLength == 0) {
                break;
            }

            GLsizei bufSize = msgLength;
            GLenum source;
            GLenum type;
            GLuint id;
            GLenum severity;

            GLsizei length;
            std::string messageLog(bufSize, 0);

            glGetDebugMessageLog(1 /*one message*/, bufSize, &source, &type, &id, &severity, &length, &messageLog[0]);

            LogLevel loglevel;
            switch(severity) {
                case GL_DEBUG_SEVERITY_NOTIFICATION:
                    loglevel = tgt::Debug;
                    break;
                case GL_DEBUG_SEVERITY_LOW:
                    loglevel = tgt::Debug;
                    break;
                case GL_DEBUG_SEVERITY_MEDIUM:
                    loglevel = tgt::Debug;
                    break;
                case GL_DEBUG_SEVERITY_HIGH:
                    loglevel = tgt::Error;
                    break;
                default:
                    tgtAssert(false, "Invalid OpenGL debug severity level");
                    loglevel = tgt::Error;
            }

            LogMgr.log(loggerCat.str(), loglevel, messageLog, locationInfo.str());
        }
    }
#endif
}

GLenum _lGLError(int line, const char* file) {

    GLenum error = glGetError();
    if (error != GL_NO_ERROR) {
        const GLubyte* exp = gluErrorString(error);

        std::ostringstream msg, info, loggerCat;
        info << " File: " << file << "@" << line;
        msg << (exp ? (const char*)exp : "unknown");
        loggerCat << "gl-error:" << file << ':' << line;

        LogMgr.log(loggerCat.str(), tgt::Error, msg.str(), info.str());
    }

    printGlDebugLog(line, file);

    return error;
}

// GLU code taken from Mesa3D http://www.mesa3d.org
// Licensed under MIT License
//
// Copyright (C) 1999-2007  Brian Paul   All Rights Reserved.
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:

static int __gluInvertMatrixd(const GLdouble m[16], GLdouble invOut[16])
{
    double inv[16], det;
    int i;

    inv[0] =   m[5]*m[10]*m[15] - m[5]*m[11]*m[14] - m[9]*m[6]*m[15]
    + m[9]*m[7]*m[14] + m[13]*m[6]*m[11] - m[13]*m[7]*m[10];
    inv[4] =  -m[4]*m[10]*m[15] + m[4]*m[11]*m[14] + m[8]*m[6]*m[15]
    - m[8]*m[7]*m[14] - m[12]*m[6]*m[11] + m[12]*m[7]*m[10];
    inv[8] =   m[4]*m[9]*m[15] - m[4]*m[11]*m[13] - m[8]*m[5]*m[15]
    + m[8]*m[7]*m[13] + m[12]*m[5]*m[11] - m[12]*m[7]*m[9];
    inv[12] = -m[4]*m[9]*m[14] + m[4]*m[10]*m[13] + m[8]*m[5]*m[14]
    - m[8]*m[6]*m[13] - m[12]*m[5]*m[10] + m[12]*m[6]*m[9];
    inv[1] =  -m[1]*m[10]*m[15] + m[1]*m[11]*m[14] + m[9]*m[2]*m[15]
    - m[9]*m[3]*m[14] - m[13]*m[2]*m[11] + m[13]*m[3]*m[10];
    inv[5] =   m[0]*m[10]*m[15] - m[0]*m[11]*m[14] - m[8]*m[2]*m[15]
    + m[8]*m[3]*m[14] + m[12]*m[2]*m[11] - m[12]*m[3]*m[10];
    inv[9] =  -m[0]*m[9]*m[15] + m[0]*m[11]*m[13] + m[8]*m[1]*m[15]
    - m[8]*m[3]*m[13] - m[12]*m[1]*m[11] + m[12]*m[3]*m[9];
    inv[13] =  m[0]*m[9]*m[14] - m[0]*m[10]*m[13] - m[8]*m[1]*m[14]
    + m[8]*m[2]*m[13] + m[12]*m[1]*m[10] - m[12]*m[2]*m[9];
    inv[2] =   m[1]*m[6]*m[15] - m[1]*m[7]*m[14] - m[5]*m[2]*m[15]
    + m[5]*m[3]*m[14] + m[13]*m[2]*m[7] - m[13]*m[3]*m[6];
    inv[6] =  -m[0]*m[6]*m[15] + m[0]*m[7]*m[14] + m[4]*m[2]*m[15]
    - m[4]*m[3]*m[14] - m[12]*m[2]*m[7] + m[12]*m[3]*m[6];
    inv[10] =  m[0]*m[5]*m[15] - m[0]*m[7]*m[13] - m[4]*m[1]*m[15]
    + m[4]*m[3]*m[13] + m[12]*m[1]*m[7] - m[12]*m[3]*m[5];
    inv[14] = -m[0]*m[5]*m[14] + m[0]*m[6]*m[13] + m[4]*m[1]*m[14]
    - m[4]*m[2]*m[13] - m[12]*m[1]*m[6] + m[12]*m[2]*m[5];
    inv[3] =  -m[1]*m[6]*m[11] + m[1]*m[7]*m[10] + m[5]*m[2]*m[11]
    - m[5]*m[3]*m[10] - m[9]*m[2]*m[7] + m[9]*m[3]*m[6];
    inv[7] =   m[0]*m[6]*m[11] - m[0]*m[7]*m[10] - m[4]*m[2]*m[11]
    + m[4]*m[3]*m[10] + m[8]*m[2]*m[7] - m[8]*m[3]*m[6];
    inv[11] = -m[0]*m[5]*m[11] + m[0]*m[7]*m[9] + m[4]*m[1]*m[11]
    - m[4]*m[3]*m[9] - m[8]*m[1]*m[7] + m[8]*m[3]*m[5];
    inv[15] =  m[0]*m[5]*m[10] - m[0]*m[6]*m[9] - m[4]*m[1]*m[10]
    + m[4]*m[2]*m[9] + m[8]*m[1]*m[6] - m[8]*m[2]*m[5];

    det = m[0]*inv[0] + m[1]*inv[4] + m[2]*inv[8] + m[3]*inv[12];
    if (det == 0)
        return GL_FALSE;

    det = 1.0 / det;

    for (i = 0; i < 16; i++)
        invOut[i] = inv[i] * det;

    return GL_TRUE;
}

static void __gluMultMatricesd(const GLdouble a[16], const GLdouble b[16],
                               GLdouble r[16])
{
    int i, j;

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            r[i*4+j] =
            a[i*4+0]*b[0*4+j] +
            a[i*4+1]*b[1*4+j] +
            a[i*4+2]*b[2*4+j] +
            a[i*4+3]*b[3*4+j];
        }
    }
}

static void __gluMultMatrixVecd(const GLdouble matrix[16], const GLdouble in[4],
                                GLdouble out[4])
{
    int i;

    for (i=0; i<4; i++) {
        out[i] =
        in[0] * matrix[0*4+i] +
        in[1] * matrix[1*4+i] +
        in[2] * matrix[2*4+i] +
        in[3] * matrix[3*4+i];
    }
}

GLint
gluProject(GLdouble objx, GLdouble objy, GLdouble objz,
           const GLdouble modelMatrix[16],
           const GLdouble projMatrix[16],
           const GLint viewport[4],
           GLdouble *winx, GLdouble *winy, GLdouble *winz)
{
    double in[4];
    double out[4];

    in[0]=objx;
    in[1]=objy;
    in[2]=objz;
    in[3]=1.0;
    __gluMultMatrixVecd(modelMatrix, in, out);
    __gluMultMatrixVecd(projMatrix, out, in);
    if (in[3] == 0.0) return(GL_FALSE);
    in[0] /= in[3];
    in[1] /= in[3];
    in[2] /= in[3];
    /* Map x, y and z to range 0-1 */
    in[0] = in[0] * 0.5 + 0.5;
    in[1] = in[1] * 0.5 + 0.5;
    in[2] = in[2] * 0.5 + 0.5;

    /* Map x,y to viewport */
    in[0] = in[0] * viewport[2] + viewport[0];
    in[1] = in[1] * viewport[3] + viewport[1];

    *winx=in[0];
    *winy=in[1];
    *winz=in[2];
    return(GL_TRUE);
}

GLint
gluUnProject(GLdouble winx, GLdouble winy, GLdouble winz,
             const GLdouble modelMatrix[16],
             const GLdouble projMatrix[16],
             const GLint viewport[4],
             GLdouble *objx, GLdouble *objy, GLdouble *objz)
{
    double finalMatrix[16];
    double in[4];
    double out[4];

    __gluMultMatricesd(modelMatrix, projMatrix, finalMatrix);
    if (!__gluInvertMatrixd(finalMatrix, finalMatrix)) return(GL_FALSE);

    in[0]=winx;
    in[1]=winy;
    in[2]=winz;
    in[3]=1.0;

    /* Map x and y from window coordinates */
    in[0] = (in[0] - viewport[0]) / viewport[2];
    in[1] = (in[1] - viewport[1]) / viewport[3];

    /* Map to range -1 to 1 */
    in[0] = in[0] * 2 - 1;
    in[1] = in[1] * 2 - 1;
    in[2] = in[2] * 2 - 1;

    __gluMultMatrixVecd(finalMatrix, in, out);
    if (out[3] == 0.0) return(GL_FALSE);
    out[0] /= out[3];
    out[1] /= out[3];
    out[2] /= out[3];
    *objx = out[0];
    *objy = out[1];
    *objz = out[2];
    return(GL_TRUE);
}

struct token_string
{
    GLuint Token;
    const char *String;
};
static const struct token_string Errors[] = {
    { GL_NO_ERROR, "no error" },
    { GL_INVALID_ENUM, "invalid enumerant" },
    { GL_INVALID_VALUE, "invalid value" },
    { GL_INVALID_OPERATION, "invalid operation" },
    { GL_INVALID_FRAMEBUFFER_OPERATION, "invalid framebuffer operation" },
    { GL_OUT_OF_MEMORY, "out of memory" },
    { static_cast<GLuint>(~0), NULL } /* end of list indicator */
};
const GLubyte*
gluErrorString(GLenum errorCode)
{
    int i;
    for (i = 0; Errors[i].String; i++) {
        if (Errors[i].Token == errorCode)
            return (const GLubyte *) Errors[i].String;
    }
    return (const GLubyte *) "unknown error";
}

} // namespace tgt
