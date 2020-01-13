/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2020 University of Muenster, Germany,           *
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

#ifndef TGT_GLMATH_H
#define TGT_GLMATH_H

#include "tgt/tgt_gl.h"

#include "tgt/matrix.h"
#include "tgt/matrixstack.h"
#include "tgt/plane.h"


/*
    The purpose of this header is to divide math and OpenGL functionality.
    Files like vector.h should not contain any OpenGL dependencies.
*/

namespace tgt {

#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE

/**
 * Transposes m and loads this matrix to the current matrix stack.
 * Thus a conversion from C++ to OpenGL has been accomplished.
 */
inline void loadMatrix(const mat4& m) {
#if defined(__APPLE__) || defined(TGT_WORKAROUND_GLLOADTRANSPOSEMATRIXF)
    glLoadMatrixf(transpose(m).elem);
#else
    //glLoadTransposeMatrixf(m.elem);
    glLoadMatrixf(transpose(m).elem);
#endif
}

/**
 * Does _NOT_ transpose m and loads this matrix to the current matrix stack.
 * Thus a transposed conversion from C++ to OpenGL has been accomplished.
 */
inline void loadTransposeMatrix(const mat4& m) {
    glLoadMatrixf(m.elem);
}

/**
 * Transposes m and multiplies this matrix to the current matrix stack.
 * Thus a conversion from C++ to OpenGL has been accomplished.
 */
inline void multMatrix(const mat4& m) {
    glMultTransposeMatrixf(m.elem);
}

/**
 * Does _NOT_ transposes m and multiplies this matrix to the current matrix stack.
 * Thus a transposed conversion from C++ to OpenGL has been accomplished.
 */
inline void multTransposeMatrix(const mat4& m) {
    glMultMatrixf(m.elem);
}

/**
 * Gets the current projection matrix and transposes it.
 * Thus a conversion from C++ to OpenGL has been accomplished.
 */
inline mat4 getProjectionMatrix() {
    mat4 m;
    glGetFloatv(GL_TRANSPOSE_PROJECTION_MATRIX, m.elem);
    return m;
}

/**
 * Gets the current projection matrix and does _NOT_ transpose it.
 * Thus a transposed conversion from C++ to OpenGL has been accomplished.
 */
inline mat4 getTransposeProjectionMatrix() {
    mat4 m;
    glGetFloatv(GL_PROJECTION_MATRIX, m.elem);
    return m;
}

/**
 * Gets the current model view matrix and transposes it.
 * Thus a conversion from C++ to OpenGL has been accomplished.
 */
inline mat4 getModelViewMatrix() {
    mat4 m;
    glGetFloatv(GL_TRANSPOSE_MODELVIEW_MATRIX, m.elem);
    return m;
}

/**
 * Gets the current model view matrix and does _NOT_ transpose it.
 * Thus a transposed conversion from C++ to OpenGL has been accomplished.
 */
inline mat4 getTransposeModelViewMatrix() {
    mat4 m;
    glGetFloatv(GL_MODELVIEW_MATRIX, m.elem);
    return m;
}

/**
 * Gets the current texture matrix and transposes it.
 * Thus a conversion from C++ to OpenGL has been accomplished.
 */
inline mat4 getTextureMatrix() {
    mat4 m;
    glGetFloatv(GL_TRANSPOSE_TEXTURE_MATRIX, m.elem);
    return m;
}

/**
 * Gets the current texture matrix and does _NOT_ transpose it.
 * Thus a transposed conversion from C++ to OpenGL has been accomplished.
 */
inline mat4 getTransposeTextureMatrix() {
    mat4 m;
    glGetFloatv(GL_TEXTURE_MATRIX, m.elem);
    return m;
}

/// sets the current clip plane in OpenGL
inline void setClipPlane(GLenum planeName, const plane& p) {
    dplane dp(p);
    glClipPlane(planeName, dp.n.elem);
}

/// sets the current clip plane in OpenGL
inline void setClipPlane(GLenum planeName, const dplane& p) {
    glClipPlane(planeName, p.n.elem);
}

/// gets the current viewport from OpenGL
inline ivec4 getViewport() {
    ivec4 v;
    glGetIntegerv(GL_VIEWPORT, (GLint *) v.elem);
    return v;
}

#endif

/// gets the current viewport from OpenGL
inline void getViewport(ivec2& pos, ivec2& size) {
    ivec4 v;
    glGetIntegerv(GL_VIEWPORT, (GLint *) v.elem);
    pos = v.xy();
    size = v.zw();
}

/*
    glVertex* wrapper
*/

// float wrapper

/// calls the appropriate version of glVertex*
inline void vertex(const Vector2<float>& v) {
    glVertex2fv(v.elem);
}

/// calls the appropriate version of glVertex*
inline void vertex(const Vector3<float>& v) {
    glVertex3fv(v.elem);
}

/// calls the appropriate version of glVertex*
inline void vertex(const Vector4<float>& v) {
    glVertex4fv(v.elem);
}

// double wrapper

/// calls the appropriate version of glVertex*
inline void vertex(const Vector2<double>& v) {
    glVertex2dv(v.elem);
}

/// calls the appropriate version of glVertex*
inline void vertex(const Vector3<double>& v) {
    glVertex3dv(v.elem);
}

/// calls the appropriate version of glVertex*
inline void vertex(const Vector4<double>& v) {
    glVertex4dv(v.elem);
}

// int wrapper

/// calls the appropriate version of glVertex*
inline void vertex(const Vector2<GLint>& v) {
    glVertex2iv(v.elem);
}

/// calls the appropriate version of glVertex*
inline void vertex(const Vector3<GLint>& v) {
    glVertex3iv(v.elem);
}

/// calls the appropriate version of glVertex*
inline void vertex(const Vector4<GLint>& v) {
    glVertex4iv(v.elem);
}

// short wrapper

/// calls the appropriate version of glVertex*
inline void vertex(const Vector2<GLshort>& v) {
    glVertex2sv(v.elem);
}

/// calls the appropriate version of glVertex*
inline void vertex(const Vector3<GLshort>& v) {
    glVertex3sv(v.elem);
}

/// calls the appropriate version of glVertex*
inline void vertex(const Vector4<GLshort>& v) {
    glVertex4sv(v.elem);
}

/*
    glTexCoord* wrapper
*/

// float wrapper

/// calls the appropriate version of glTexCoord*
inline void texCoord(const Vector2<float>& v) {
    glTexCoord2fv(v.elem);
}

/// calls the appropriate version of glTexCoord*
inline void texCoord(const Vector3<float>& v) {
    glTexCoord3fv(v.elem);
}

/// calls the appropriate version of glTexCoord*
inline void texCoord(const Vector4<float>& v) {
    glTexCoord4fv(v.elem);
}

// double wrapper

/// calls the appropriate version of glTexCoord*
inline void texCoord(const Vector2<double>& v) {
    glTexCoord2dv(v.elem);
}

/// calls the appropriate version of glTexCoord*
inline void texCoord(const Vector3<double>& v) {
    glTexCoord3dv(v.elem);
}

/// calls the appropriate version of glTexCoord*
inline void texCoord(const Vector4<double>& v) {
    glTexCoord4dv(v.elem);
}

// int wrapper

/// calls the appropriate version of glTexCoord*
inline void texCoord(const Vector2<GLint>& v) {
    glTexCoord2iv(v.elem);
}

/// calls the appropriate version of glTexCoord*
inline void texCoord(const Vector3<GLint>& v) {
    glTexCoord3iv(v.elem);
}

/// calls the appropriate version of glTexCoord*
inline void texCoord(const Vector4<GLint>& v) {
    glTexCoord4iv(v.elem);
}

// short wrapper

/// calls the appropriate version of glTexCoord*
inline void texCoord(const Vector2<GLshort>& v) {
    glTexCoord2sv(v.elem);
}

/// calls the appropriate version of glTexCoord*
inline void texCoord(const Vector3<GLshort>& v) {
    glTexCoord3sv(v.elem);
}

/// calls the appropriate version of glTexCoord*
inline void texCoord(const Vector4<GLshort>& v) {
    glTexCoord4sv(v.elem);
}

#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE

/*
    glMultiTexCoord* wrapper
*/

// float wrapper

/// calls the appropriate version of glMultiTexCoord*
inline void multiTexCoord(GLenum target, const Vector2<float>& v) {
    glMultiTexCoord2fv(target, v.elem);
}

/// calls the appropriate version of glMultiTexCoord*
inline void multiTexCoord(GLenum target, const Vector3<float>& v) {
    glMultiTexCoord3fv(target, v.elem);
}

/// calls the appropriate version of glMultiTexCoord*
inline void multiTexCoord(GLenum target, const Vector4<float>& v) {
    glMultiTexCoord4fv(target, v.elem);
}

// double wrapper

/// calls the appropriate version of glMultiTexCoord*
inline void multiTexCoord(GLenum target, const Vector2<double>& v) {
    glMultiTexCoord2dv(target, v.elem);
}

/// calls the appropriate version of glMultiTexCoord*
inline void multiTexCoord(GLenum target, const Vector3<double>& v) {
    glMultiTexCoord3dv(target, v.elem);
}

/// calls the appropriate version of glMultiTexCoord*
inline void multiTexCoord(GLenum target, const Vector4<double>& v) {
    glMultiTexCoord4dv(target, v.elem);
}

// int wrapper

/// calls the appropriate version of glMultiTexCoord*
inline void multiTexCoord(GLenum target, const Vector2<GLint>& v) {
    glMultiTexCoord2iv(target, v.elem);
}

/// calls the appropriate version of glMultiTexCoord*
inline void multiTexCoord(GLenum target, const Vector3<GLint>& v) {
    glMultiTexCoord3iv(target, v.elem);
}

/// calls the appropriate version of glMultiTexCoord*
inline void multiTexCoord(GLenum target, const Vector4<GLint>& v) {
    glMultiTexCoord4iv(target, v.elem);
}

// short wrapper

/// calls the appropriate version of glMultiTexCoord*
inline void multiTexCoord(GLenum target, const Vector2<GLshort>& v) {
    glMultiTexCoord2sv(target, v.elem);
}

/// calls the appropriate version of glMultiTexCoord*
inline void multiTexCoord(GLenum target, const Vector3<GLshort>& v) {
    glMultiTexCoord3sv(target, v.elem);
}

/// calls the appropriate version of glMultiTexCoord*
inline void multiTexCoord(GLenum target, const Vector4<GLshort>& v) {
    glMultiTexCoord4sv(target, v.elem);
}

#endif

/*
    glNormal* wrapper
*/

/// calls the appropriate version of glNormal*
inline void normal(const Vector3<GLbyte>& v) {
    glNormal3bv(v.elem);
}

/// calls the appropriate version of glNormal*
inline void normal(const Vector3<GLshort>& v) {
    glNormal3sv(v.elem);
}

/// calls the appropriate version of glNormal*
inline void normal(const Vector3<GLint>& v) {
    glNormal3iv(v.elem);
}

/// calls the appropriate version of glNormal*
inline void normal(const Vector3<float>& v) {
    glNormal3fv(v.elem);
}

/// calls the appropriate version of glNormal*
inline void normal(const Vector3<double>& v) {
    glNormal3dv(v.elem);
}


#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE

/*
    wrapper for glTranslate* and glScale*
*/

/// calls the appropriate version of glTranslate*
inline void translate(const Vector3<float>& v) {
    glTranslatef(v.x, v.y, v.z);
}

/// calls the appropriate version of glTranslate*
inline void translate(const Vector3<double>& v) {
    glTranslated(v.x, v.y, v.z);
}

/// calls the appropriate version of glTranslate*.
/// z-translation is set to 0.f
inline void translate(const Vector2<float>& v) {
    glTranslatef(v.x, v.y, 0.f);
}

/// calls the appropriate version of glTranslate*
/// z-translation is set to 0.0
inline void translate(const Vector2<double>& v) {
    glTranslated(v.x, v.y, 0.0);
}

/// calls the appropriate version of glScale*
inline void scale(const Vector3<float>& v) {
    glScalef(v.x, v.y, v.z);
}

/// calls the appropriate version of glScale*
inline void scale(const Vector3<double>& v) {
    glScaled(v.x, v.y, v.z);
}

/// calls the appropriate version of glScale*.
/// z-scale factor is set to 1.f
inline void scale(const Vector2<float>& v) {
    glScalef(v.x, v.y, 1.f);
}

/// calls the appropriate version of glScale*.
/// z-scale factor is set to 1.0
inline void scale(const Vector2<double>& v) {
    glScaled(v.x, v.y, 1.0);
}

#endif

} // namespace tgt

#endif //TGT_GLMATH_H
