/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2024 University of Muenster, Germany,           *
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

#ifndef VRN_IMGL_H
#define VRN_IMGL_H

#include "tgt/types.h"

/**
 *  This is the compatibility header to map
 *  legacy immediate mode calls to tgt::ImmediateMode.
 */


enum {
    GL_POLYGON
};

TGT_API void  glBegin(GLenum  mode);
TGT_API void  glEnd(void);

TGT_API void  glVertex2s(GLshort x, GLshort y);
TGT_API void  glVertex2i(GLint x, GLint y);
TGT_API void  glVertex2f(GLfloat x, GLfloat y);
TGT_API void  glVertex2d(GLdouble x, GLdouble y);
TGT_API void  glVertex3s(GLshort x, GLshort y, GLshort z);
TGT_API void  glVertex3i(GLint x, GLint y, GLint z);
TGT_API void  glVertex3f(GLfloat x, GLfloat y, GLfloat z);
TGT_API void  glVertex3d(GLdouble x, GLdouble y, GLdouble z);
TGT_API void  glVertex4s(GLshort x, GLshort y, GLshort z, GLshort w);
TGT_API void  glVertex4i(GLint x, GLint y, GLint z, GLint w);
TGT_API void  glVertex4f(GLfloat x, GLfloat y, GLfloat z, GLfloat w);
TGT_API void  glVertex4d(GLdouble x, GLdouble y, GLdouble z, GLdouble w);

TGT_API void  glVertex2sv(const GLshort *v);
TGT_API void  glVertex2iv(const GLint *v);
TGT_API void  glVertex2fv(const GLfloat *v);
TGT_API void  glVertex2dv(const GLdouble *v);
TGT_API void  glVertex3sv(const GLshort *v);
TGT_API void  glVertex3iv(const GLint *v);
TGT_API void  glVertex3fv(const GLfloat *v);
TGT_API void  glVertex3dv(const GLdouble *v);
TGT_API void  glVertex4sv(const GLshort *v);
TGT_API void  glVertex4iv(const GLint *v);
TGT_API void  glVertex4fv(const GLfloat *v);
TGT_API void  glVertex4dv(const GLdouble *v);

TGT_API void  glTexCoord1s(GLshort s);
TGT_API void  glTexCoord1i(GLint s);
TGT_API void  glTexCoord1f(GLfloat s);
TGT_API void  glTexCoord1d(GLdouble s);
TGT_API void  glTexCoord2s(GLshort s, GLshort t);
TGT_API void  glTexCoord2i(GLint s, GLint t);
TGT_API void  glTexCoord2f(GLfloat s, GLfloat t);
TGT_API void  glTexCoord2d(GLdouble s, GLdouble t);
TGT_API void  glTexCoord3s(GLshort s, GLshort t, GLshort r);
TGT_API void  glTexCoord3i(GLint s, GLint t, GLint r);
TGT_API void  glTexCoord3f(GLfloat s, GLfloat t, GLfloat r);
TGT_API void  glTexCoord3d(GLdouble s, GLdouble t, GLdouble r);
TGT_API void  glTexCoord4s(GLshort s, GLshort t, GLshort r, GLshort q);
TGT_API void  glTexCoord4i(GLint s, GLint t, GLint r, GLint q);
TGT_API void  glTexCoord4f(GLfloat s, GLfloat t, GLfloat r, GLfloat q);
TGT_API void  glTexCoord4d(GLdouble s, GLdouble t, GLdouble r, GLdouble q);

TGT_API void  glTexCoord1sv(const GLshort *v);
TGT_API void  glTexCoord1iv(const GLint *v);
TGT_API void  glTexCoord1fv(const GLfloat *v);
TGT_API void  glTexCoord1dv(const GLdouble *v);
TGT_API void  glTexCoord2sv(const GLshort *v);
TGT_API void  glTexCoord2iv(const GLint *v);
TGT_API void  glTexCoord2fv(const GLfloat *v);
TGT_API void  glTexCoord2dv(const GLdouble *v);
TGT_API void  glTexCoord3sv(const GLshort *v);
TGT_API void  glTexCoord3iv(const GLint *v);
TGT_API void  glTexCoord3fv(const GLfloat *v);
TGT_API void  glTexCoord3dv(const GLdouble *v);
TGT_API void  glTexCoord4sv(const GLshort *v);
TGT_API void  glTexCoord4iv(const GLint *v);
TGT_API void  glTexCoord4fv(const GLfloat *v);
TGT_API void  glTexCoord4dv(const GLdouble *v);

TGT_API void  glNormal3b(GLbyte nx, GLbyte ny, GLbyte nz);
TGT_API void  glNormal3d(GLdouble nx, GLdouble ny, GLdouble nz);
TGT_API void  glNormal3f(GLfloat nx, GLfloat ny, GLfloat nz);
TGT_API void  glNormal3i(GLint nx, GLint ny, GLint nz);
TGT_API void  glNormal3s(GLshort nx, GLshort ny, GLshort nz);

TGT_API void  glNormal3bv(const GLbyte *v);
TGT_API void  glNormal3dv(const GLdouble *v);
TGT_API void  glNormal3fv(const GLfloat *v);
TGT_API void  glNormal3iv(const GLint *v);
TGT_API void  glNormal3sv(const GLshort *v);

TGT_API void  glColor3b(GLbyte red, GLbyte green, GLbyte blue);
TGT_API void  glColor3s(GLshort red, GLshort green, GLshort blue);
TGT_API void  glColor3i(GLint red, GLint green, GLint blue);
TGT_API void  glColor3f(GLfloat red, GLfloat green, GLfloat blue);
TGT_API void  glColor3d(GLdouble red, GLdouble green, GLdouble blue);
TGT_API void  glColor3ub(GLubyte red, GLubyte green, GLubyte blue);
TGT_API void  glColor3us(GLushort red, GLushort green, GLushort blue);
TGT_API void  glColor3ui(GLuint red, GLuint green, GLuint blue);
TGT_API void  glColor4b(GLbyte red, GLbyte green, GLbyte blue, GLbyte alpha);
TGT_API void  glColor4s(GLshort red, GLshort green, GLshort blue, GLshort alpha);
TGT_API void  glColor4i(GLint red, GLint green, GLint blue, GLint alpha);
TGT_API void  glColor4f(GLfloat red, GLfloat green, GLfloat blue, GLfloat alpha);
TGT_API void  glColor4d(GLdouble red, GLdouble green, GLdouble blue, GLdouble alpha);
TGT_API void  glColor4ub(GLubyte red, GLubyte green, GLubyte blue, GLubyte alpha);
TGT_API void  glColor4us(GLushort red, GLushort green, GLushort blue, GLushort alpha);
TGT_API void  glColor4ui(GLuint red, GLuint green, GLuint blue, GLuint alpha);

TGT_API void  glColor3bv(const GLbyte *v);
TGT_API void  glColor3sv(const GLshort *v);
TGT_API void  glColor3iv(const GLint *v);
TGT_API void  glColor3fv(const GLfloat *v);
TGT_API void  glColor3dv(const GLdouble *v);
TGT_API void  glColor3ubv(const GLubyte *v);
TGT_API void  glColor3usv(const GLushort *v);
TGT_API void  glColor3uiv(const GLuint *v);
TGT_API void  glColor4bv(const GLbyte *v);
TGT_API void  glColor4sv(const GLshort *v);
TGT_API void  glColor4iv(const GLint *v);
TGT_API void  glColor4fv(const GLfloat *v);
TGT_API void  glColor4dv(const GLdouble *v);
TGT_API void  glColor4ubv(const GLubyte *v);
TGT_API void  glColor4usv(const GLushort *v);
TGT_API void  glColor4uiv(const GLuint *v);

#endif // VRN_VRN_IMGL_H_H
