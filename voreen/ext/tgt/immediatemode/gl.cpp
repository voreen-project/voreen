/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2021 University of Muenster, Germany,           *
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
#include "tgt/immediatemode/immediatemode.h"

using namespace tgt;

// Shortcut for the normalize method
#define s(x) IMode.normalize(x)

void glBegin(GLenum mode) {
    IMode.begin(static_cast<ImmediateMode::PolygonMode>(mode));
}

void glEnd(void) {
    IMode.end();
}

void glVertex2s(GLshort x, GLshort y) {
    IMode.vertex(vec2(x,y));
}
void glVertex2i(GLint x, GLint y) {
    IMode.vertex(vec2(x,y));
}
void glVertex2f(GLfloat x, GLfloat y) {
    IMode.vertex(vec2(x,y));
}
void glVertex2d(GLdouble x, GLdouble y) {
    IMode.vertex(vec2(x,y));
}
void glVertex3s(GLshort x, GLshort y, GLshort z) {
    IMode.vertex(vec2(x,y));
}
void glVertex3i(GLint x, GLint y, GLint z) {
    IMode.vertex(vec3(x, y, z));
}
void glVertex3f(GLfloat x, GLfloat y, GLfloat z) {
    IMode.vertex(vec3(x, y, z));
}
void glVertex3d(GLdouble x, GLdouble y, GLdouble z) {
    IMode.vertex(vec3(x, y, z));
}
void glVertex4s(GLshort x, GLshort y, GLshort z, GLshort w) {
    IMode.vertex(vec4(x, y, z, w));
}
void glVertex4i(GLint x, GLint y, GLint z, GLint w) {
    IMode.vertex(vec4(x, y, z, w));
}
void glVertex4f(GLfloat x, GLfloat y, GLfloat z, GLfloat w) {
    IMode.vertex(vec4(x, y, z, w));
}
void glVertex4d(GLdouble x, GLdouble y, GLdouble z, GLdouble w) {
    IMode.vertex(vec4(x, y, z, w));
}
void glVertex2sv(const GLshort *v) {
    IMode.vertex(vec2(*v, *(v+1)));
}
void glVertex2iv(const GLint *v) {
    IMode.vertex(vec2(*v, *(v+1)));
}
void glVertex2fv(const GLfloat *v) {
    IMode.vertex(vec2(*v, *(v+1)));
}
void glVertex2dv(const GLdouble *v) {
    IMode.vertex(vec2(*v, *(v+1)));
}
void glVertex3sv(const GLshort *v) {
    IMode.vertex(vec3(*v, *(v+1), *(v+2)));
}
void glVertex3iv(const GLint *v) {
    IMode.vertex(vec3(*v, *(v+1), *(v+2)));
}
void glVertex3fv(const GLfloat *v) {
    IMode.vertex(vec3(*v, *(v+1), *(v+2)));
}
void glVertex3dv(const GLdouble *v) {
    IMode.vertex(vec3(*v, *(v+1), *(v+2)));
}
void glVertex4sv(const GLshort *v) {
    IMode.vertex(vec4(*v, *(v+1), *(v+2), *(v+3)));
}
void glVertex4iv(const GLint *v) {
    IMode.vertex(vec4(*v, *(v+1), *(v+2), *(v+3)));
}
void glVertex4fv(const GLfloat *v) {
    IMode.vertex(vec4(*v, *(v+1), *(v+2), *(v+3)));
}
void glVertex4dv(const GLdouble *v) {
    IMode.vertex(vec4(*v, *(v+1), *(v+2), *(v+3)));
}

void glTexCoord1s(GLshort s) {
    IMode.texcoord(s);
}
void glTexCoord1i(GLint s) {}
void glTexCoord1f(GLfloat s) {
    IMode.texcoord(s);
}
void glTexCoord1d(GLdouble s) {
    IMode.texcoord(s);
}
void glTexCoord2s(GLshort s, GLshort t) {
    IMode.texcoord(vec2(s, t));
}
void glTexCoord2i(GLint s, GLint t) {
    IMode.texcoord(vec2(s, t));
}
void glTexCoord2f(GLfloat s, GLfloat t) {
    IMode.texcoord(vec2(s, t));
}
void glTexCoord2d(GLdouble s, GLdouble t) {
    IMode.texcoord(vec2(s, t));
}
void glTexCoord3s(GLshort s, GLshort t, GLshort r) {
    IMode.texcoord(vec3(s, t, r));
}
void glTexCoord3i(GLint s, GLint t, GLint r) {
    IMode.texcoord(vec3(s, t, r));
}
void glTexCoord3f(GLfloat s, GLfloat t, GLfloat r) {
    IMode.texcoord(vec3(s, t, r));
}
void glTexCoord3d(GLdouble s, GLdouble t, GLdouble r) {
    IMode.texcoord(vec3(s, t, r));
}
void glTexCoord4s(GLshort s, GLshort t, GLshort r, GLshort q) {
    IMode.texcoord(vec4(s, t, r, q));
}
void glTexCoord4i(GLint s, GLint t, GLint r, GLint q) {
    IMode.texcoord(vec4(s, t, r, q));
}
void glTexCoord4f(GLfloat s, GLfloat t, GLfloat r, GLfloat q) {
    IMode.texcoord(vec4(s, t, r, q));
}
void glTexCoord4d(GLdouble s, GLdouble t, GLdouble r, GLdouble q) {
    IMode.texcoord(vec4(s, t, r, q));
}
void glTexCoord1sv(const GLshort *v) {
     IMode.texcoord(*v);
}
void glTexCoord1iv(const GLint *v) {
    IMode.texcoord(*v);
}
void glTexCoord1fv(const GLfloat *v) {
    IMode.texcoord(*v);
}
void glTexCoord1dv(const GLdouble *v) {
    IMode.texcoord(*v);
}
void glTexCoord2sv(const GLshort *v) {
    IMode.texcoord(vec2(*v, *(v+1)));
}
void glTexCoord2iv(const GLint *v) {
    IMode.texcoord(vec2(*v, *(v+1)));
}
void glTexCoord2fv(const GLfloat *v) {
    IMode.texcoord(vec2(*v, *(v+1)));
}
void glTexCoord2dv(const GLdouble *v) {
    IMode.texcoord(vec2(*v, *(v+1)));
}
void glTexCoord3sv(const GLshort *v) {
    IMode.texcoord(vec3(*v, *(v+1), *(v+2)));
}
void glTexCoord3iv(const GLint *v) {
    IMode.texcoord(vec3(*v, *(v+1), *(v+2)));
}
void glTexCoord3fv(const GLfloat *v) {
    IMode.texcoord(vec3(*v, *(v+1), *(v+2)));
}
void glTexCoord3dv(const GLdouble *v) {
    IMode.texcoord(vec3(*v, *(v+1), *(v+2)));
}
void glTexCoord4sv(const GLshort *v) {
    IMode.texcoord(vec4(*v, *(v+1), *(v+2), *(v+3)));
}
void glTexCoord4iv(const GLint *v) {
    IMode.texcoord(vec4(*v, *(v+1), *(v+2), *(v+3)));
}
void glTexCoord4fv(const GLfloat *v) {
    IMode.texcoord(vec4(*v, *(v+1), *(v+2), *(v+3)));
}
void glTexCoord4dv(const GLdouble *v) {
    IMode.texcoord(vec4(*v, *(v+1), *(v+2), *(v+3)));
}

void glNormal3b(GLbyte nx, GLbyte ny, GLbyte nz) {
    IMode.normal(vec3(s(nx), s(ny), s(nz)));
}
void glNormal3d(GLdouble nx, GLdouble ny, GLdouble nz) {
    IMode.normal(vec3(nx, ny, nz));
}
void glNormal3f(GLfloat nx, GLfloat ny, GLfloat nz) {
    IMode.normal(vec3(nx, ny, nz));
}
void glNormal3i(GLint nx, GLint ny, GLint nz) {
    IMode.normal(vec3(s(nx), s(ny), s(nz)));
}
void glNormal3s(GLshort nx, GLshort ny, GLshort nz) {
    IMode.normal(vec3(s(nx), s(ny), s(nz)));
}

void glNormal3bv(const GLbyte *v) {
    IMode.normal(vec3(s(*v), s(*(v+1)), s(*(v+2))));
}
void glNormal3dv(const GLdouble *v) {
    IMode.normal(vec3(*v, *(v+1), *(v+2)));
}
void glNormal3fv(const GLfloat *v) {
    IMode.normal(vec3(*v, *(v+1), *(v+2)));
}
void glNormal3iv(const GLint *v) {
    IMode.normal(vec3(s(*v), s(*(v+1)), s(*(v+2))));
}
void glNormal3sv(const GLshort *v) {
    IMode.normal(vec3(s(*v), s(*(v+1)), s(*(v+2))));
}

void glColor3b(GLbyte red, GLbyte green, GLbyte blue) {
    IMode.color(vec3(s(red), s(green), s(blue)));
}
void glColor3s(GLshort red, GLshort green, GLshort blue) {
    IMode.color(vec3(s(red), s(green), s(blue)));
}
void glColor3i(GLint red, GLint green, GLint blue) {
    IMode.color(vec3(s(red), s(green), s(blue)));
}
void glColor3f(GLfloat red, GLfloat green, GLfloat blue) {
    IMode.color(vec3(red, green, blue));
}
void glColor3d(GLdouble red, GLdouble green, GLdouble blue) {
    IMode.color(vec3(red, green, blue));
}
void glColor3ub(GLubyte red, GLubyte green, GLubyte blue) {
    IMode.color(vec3(s(red), s(green), s(blue)));
}
void glColor3us(GLushort red, GLushort green, GLushort blue) {
    IMode.color(vec3(s(red), s(green), s(blue)));
}
void glColor3ui(GLuint red, GLuint green, GLuint blue) {
    IMode.color(vec3(s(red), s(green), s(blue)));
}
void glColor4b(GLbyte red, GLbyte green, GLbyte blue, GLbyte alpha) {
    IMode.color(vec3(s(red), s(green), s(blue)));
}
void glColor4s(GLshort red, GLshort green, GLshort blue, GLshort alpha) {
    IMode.color(vec3(s(red), s(green), s(blue)));
}
void glColor4i(GLint red, GLint green, GLint blue, GLint alpha) {
    IMode.color(vec3(s(red), s(green), s(blue)));
}
void glColor4f(GLfloat red, GLfloat green, GLfloat blue, GLfloat alpha) {
    IMode.color(vec4(red, green, blue, alpha));
}
void glColor4d(GLdouble red, GLdouble green, GLdouble blue, GLdouble alpha) {
    IMode.color(vec4(red, green, blue, alpha));
}
void glColor4ub(GLubyte red, GLubyte green, GLubyte blue, GLubyte alpha) {
    IMode.color(vec3(s(red), s(green), s(blue)));
}
void glColor4us(GLushort red, GLushort green, GLushort blue, GLushort alpha) {
    IMode.color(vec3(s(red), s(green), s(blue)));
}
void glColor4ui(GLuint red, GLuint green, GLuint blue, GLuint alpha) {
    IMode.color(vec3(s(red), s(green), s(blue)));
}

void glColor3bv(const GLbyte *v) {
    IMode.color(vec3(s(*v), s(*(v+1)), s(*(v+2))));
}
void glColor3sv(const GLshort *v) {
    IMode.color(vec3(s(*v), s(*(v+1)), s(*(v+2))));
}
void glColor3iv(const GLint *v) {
    IMode.color(vec3(s(*v), s(*(v+1)), s(*(v+2))));
}
void glColor3fv(const GLfloat *v) {
    IMode.color(vec3(*v, *(v+1), *(v+2)));
}
void glColor3dv(const GLdouble *v) {
    IMode.color(vec3(*v, *(v+1), *(v+2)));
}
void glColor3ubv(const GLubyte *v) {
    IMode.color(vec3(s(*v), s(*(v+1)), s(*(v+2))));
}
void glColor3usv(const GLushort *v) {
    IMode.color(vec3(s(*v), s(*(v+1)), s(*(v+2))));
}
void glColor3uiv(const GLuint *v) {
    IMode.color(vec3(s(*v), s(*(v+1)), s(*(v+2))));
}
void glColor4bv(const GLbyte *v) {
    IMode.color(vec4(s(*v), s(*(v+1)), s(*(v+2)), s(*(v+3))));
}
void glColor4sv(const GLshort *v) {
    IMode.color(vec4(s(*v), s(*(v+1)), s(*(v+2)), s(*(v+3))));
}
void glColor4iv(const GLint *v) {
    IMode.color(vec4(s(*v), s(*(v+1)), s(*(v+2)), s(*(v+3))));
}
void glColor4fv(const GLfloat *v) {
    IMode.color(vec4(*v, *(v+1), *(v+2), *(v+3)));
}
void glColor4dv(const GLdouble *v) {
    IMode.color(vec4(*v, *(v+1), *(v+2), *(v+3)));
}
void glColor4ubv(const GLubyte *v) {
    IMode.color(vec4(s(*v), s(*(v+1)), s(*(v+2)), s(*(v+3))));
}
void glColor4usv(const GLushort *v) {
    IMode.color(vec4(s(*v), s(*(v+1)), s(*(v+2)), s(*(v+3))));
}
void glColor4uiv(const GLuint *v) {
    IMode.color(vec4(s(*v), s(*(v+1)), s(*(v+2)), s(*(v+3))));
}

