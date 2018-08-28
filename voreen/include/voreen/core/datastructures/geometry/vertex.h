/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_VERTEX_H
#define VRN_VERTEX_H

#include "tgt/vector.h"
#include "tgt/plane.h"
#include "tgt/tgt_gl.h"

namespace voreen {

/*
 * Base class for vertices stored in meshes.
 * Most of the methods have to be reimplemented for subclasses.
 *
 * __WARNING__
 * The methods are intentionally not made virtual. Adding any virtual methods to this struct
 * increases the its size by sizeof(pointer to virtual table), leading to a failure when passing
 * an array of these structs to OpenGL within a vertex buffer object as the data ist not packed
 * tightly anymore.
 * __WARNING__
 */
struct VertexBase {
    enum VertexLayout {
        SIMPLE    = 0,
        COLOR     = 1,
        NORMAL    = 2,
        TEXCOORD  = 4,
        COLOR_NORMAL    = COLOR | NORMAL,
        NORMAL_TEXCOORD = NORMAL | TEXCOORD,
        COLOR_TEXCOORD  = COLOR | TEXCOORD,
        COLOR_NORMAL_TEXCOORD = COLOR | NORMAL | TEXCOORD
    };

    /// The vertex layout of this vertex class
    static const VertexLayout layout = SIMPLE;

    tgt::vec3 pos_;

    VertexBase() {}
    VertexBase(tgt::vec3 pos) : pos_(pos) {}
    VertexBase(tgt::vec3 pos, tgt::vec4 color, tgt::vec3 normal, tgt::vec2 texCoord, tgt::ivec2 texIndex) : pos_(pos) {}

    /// Compares this vertex to another vertex. Reimplement in subclass.
    bool equals(const VertexBase& other, double epsilon = 1e-5) const;

    /// Sets up vertex attributs pointers for rendering. Reimplement in subclass.
    static void setupVertexAttributePointers(size_t stride = 0);

    /// De-initializes vertex attributes pointers. Reimplement in subclass.
    static void disableVertexAttributePointers();

    /**
     * Returns the distance between the vertex geometry and the given plane.
     *
     * @note Use the @c epsilon parameter to change the accuracy at which
     *       the vertex geometry lies on the given plane.
     *
     * @param plane the plane
     * @param epsilon the accuracy at which the vertex geometry lies on the given plane
     *
     * @returns distance between vertex geometry and given plane
     */
    double getDistanceToPlane(const tgt::plane& plane, double epsilon = 1e-5) const;

    /// Interpolates two vertices of this type. Reimplement in subclass.
    static VertexBase interpolate(const VertexBase& v1, const VertexBase& v2, float t);

    /// Sets the normal of this vertex, used in clipping code. Default implementation does nothing.
    void setNormal(tgt::vec3 n) {}
    tgt::vec3 getNormal() {return tgt::vec3::zero;}
    void setColor(tgt::vec4 col) {};
    tgt::vec4 getColor() {return tgt::vec4::zero;}
    void setTexCoord(tgt::vec2 texCoord) {};
    tgt::vec2 getTexCoord() {return tgt::vec2::zero;}
    void setTexIndex(tgt::ivec2 texIndex) {};
    tgt::ivec2 getTexIndex() {return tgt::ivec2::zero;}
};

//-------------------------------------------------------------------------------------------------

/// A vertex with an additional vec3 attribute (e.g., texture coordinate).
struct VertexNormal : public VertexBase {
    /// The vertex layout of this vertex class
    static const VertexLayout layout = NORMAL;

    tgt::vec3 normal_;

    VertexNormal() {}
    VertexNormal(tgt::vec3 pos, tgt::vec3 normal);
    VertexNormal(tgt::vec3 pos, tgt::vec4 color, tgt::vec3 normal, tgt::vec2 texCoord, tgt::ivec2 texIndex);

    bool equals(const VertexNormal& other, double epsilon = 1e-5) const;

    static void setupVertexAttributePointers();
    static void disableVertexAttributePointers();

    static VertexNormal interpolate(const VertexNormal& v1, const VertexNormal& v2, float t);

    void setNormal(tgt::vec3 n) {normal_ = n;}
    tgt::vec3 getNormal() {return normal_;}
};

//-------------------------------------------------------------------------------------------------

/// A vertex with an additional vec3 attribute (e.g., texture coordinate).
struct VertexColor : public VertexBase {
    /// The vertex layout of this vertex class
    static const VertexLayout layout = COLOR;

    tgt::vec4 color_;

    VertexColor() {}
    VertexColor(tgt::vec3 pos, tgt::vec4 color);
    VertexColor(tgt::vec3 pos, tgt::vec4 color, tgt::vec3 normal, tgt::vec2 texCoord, tgt::ivec2 texIndex);

    bool equals(const VertexColor& other, double epsilon = 1e-5) const;

    static void setupVertexAttributePointers();
    static void disableVertexAttributePointers();

    static VertexColor interpolate(const VertexColor& v1, const VertexColor& v2, float t);

    void setColor(tgt::vec4 color) {color_ = color;}
    tgt::vec4 getColor() {return color_;}
};

//-------------------------------------------------------------------------------------------------

/// A vertex with an additional vec3 attribute (e.g., texture coordinate).
struct VertexTexCoord : public VertexBase {
    /// The vertex layout of this vertex class
    static const VertexLayout layout = TEXCOORD;

    tgt::vec2  texCoord_;
    tgt::ivec2 texIndex_;

    VertexTexCoord() {}
    VertexTexCoord(tgt::vec3 pos, tgt::vec2 texCoord, tgt::ivec2 texIndex);
    VertexTexCoord(tgt::vec3 pos, tgt::vec4 color, tgt::vec3 normal, tgt::vec2 texCoord, tgt::ivec2 texIndex);

    bool equals(const VertexTexCoord& other, double epsilon = 1e-5) const;

    static void setupVertexAttributePointers();
    static void disableVertexAttributePointers();

    static VertexTexCoord interpolate(const VertexTexCoord& v1, const VertexTexCoord& v2, float t);

    void setTexCoord(tgt::vec2 texCoord) {texCoord_ = texCoord;}
    tgt::vec2 getTexCoord() {return texCoord_;}
    void setTexIndex(tgt::ivec2 texIndex) {texIndex_ = texIndex;}
    tgt::ivec2 getTexIndex() {return texIndex_;}
};

//-------------------------------------------------------------------------------------------------

/// A vertex with an additional vec4 and vec3 attribute (e.g., RGBA color and texture coordinate).
struct VertexColorNormal : public VertexBase {
    /// The vertex layout of this vertex class
    static const VertexLayout layout = COLOR_NORMAL;

    tgt::vec4 color_;
    tgt::vec3 normal_;

    VertexColorNormal() {}
    VertexColorNormal(tgt::vec3 pos, tgt::vec4 color, tgt::vec3 normal);
    VertexColorNormal(tgt::vec3 pos, tgt::vec4 color, tgt::vec3 normal, tgt::vec2 texCoord, tgt::ivec2 texIndex);

    bool equals(const VertexColorNormal& other, double epsilon = 1e-5) const;

    static void setupVertexAttributePointers();

    static void disableVertexAttributePointers();

    static VertexColorNormal interpolate(const VertexColorNormal& v1, const VertexColorNormal& v2, float t);

    void setColor(tgt::vec4 color) {color_ = color;}
    void setNormal(tgt::vec3 n) {normal_ = n;}
    tgt::vec4 getColor() {return color_;}
    tgt::vec3 getNormal() {return normal_;}
};

//-------------------------------------------------------------------------------------------------
//
struct VertexNormalTexCoord : public VertexBase {
    /// The vertex layout of this vertex class
    static const VertexLayout layout = NORMAL_TEXCOORD;

    tgt::vec3 normal_;
    tgt::vec2  texCoord_;
    tgt::ivec2 texIndex_;

    VertexNormalTexCoord() {}
    VertexNormalTexCoord(tgt::vec3 pos, tgt::vec3 normal, tgt::vec2 texCoord, tgt::ivec2 texIndex);
    VertexNormalTexCoord(tgt::vec3 pos, tgt::vec4 color, tgt::vec3 normal, tgt::vec2 texCoord, tgt::ivec2 texIndex);

    bool equals(const VertexNormalTexCoord& other, double epsilon = 1e-5) const;

    static void setupVertexAttributePointers();

    static void disableVertexAttributePointers();

    static VertexNormalTexCoord interpolate(const VertexNormalTexCoord& v1, const VertexNormalTexCoord& v2, float t);

    void setNormal(tgt::vec3 n) {normal_ = n;}
    void setTexCoord(tgt::vec2 texCoord) {texCoord_ = texCoord;}
    void setTexIndex(tgt::ivec2 texIndex) {texIndex_ = texIndex;}
    tgt::vec3 getNormal() {return normal_;}
    tgt::vec2 getTexCoord() {return texCoord_;}
    tgt::ivec2 getTexIndex() {return texIndex_;}
};

//-------------------------------------------------------------------------------------------------

/// A vertex with an additional vec4 and vec3 attribute (e.g., RGBA color and texture coordinate).
struct VertexColorTexCoord : public VertexBase {
    /// The vertex layout of this vertex class
    static const VertexLayout layout = COLOR_TEXCOORD;

    tgt::vec4 color_;
    tgt::vec2 texCoord_;
    tgt::ivec2 texIndex_;

    VertexColorTexCoord() {}
    VertexColorTexCoord(tgt::vec3 pos, tgt::vec4 color, tgt::vec2 texCoord, tgt::ivec2 texIndex);
    VertexColorTexCoord(tgt::vec3 pos, tgt::vec4 color, tgt::vec3 normal, tgt::vec2 texCoord, tgt::ivec2 texIndex);


    bool equals(const VertexColorTexCoord& other, double epsilon = 1e-5) const;

    static void setupVertexAttributePointers();

    static void disableVertexAttributePointers();

    static VertexColorTexCoord interpolate(const VertexColorTexCoord& v1, const VertexColorTexCoord& v2, float t);

    void setColor(tgt::vec4 color) {color_ = color;}
    void setTexCoord(tgt::vec2 texCoord) {texCoord_ = texCoord;}
    void setTexIndex(tgt::ivec2 texIndex) {texIndex_ = texIndex;}
    tgt::vec4 getColor() {return color_;}
    tgt::vec2 getTexCoord() {return texCoord_;}
    tgt::ivec2 getTexIndex() {return texIndex_;}
};

//-------------------------------------------------------------------------------------------------

/// A vertex with an additional vec4 and vec3 attribute (e.g., RGBA color and texture coordinate).
struct VertexColorNormalTexCoord : public VertexBase {
    /// The vertex layout of this vertex class
    static const VertexLayout layout = COLOR_NORMAL_TEXCOORD;

    tgt::vec4 color_;
    tgt::vec3 normal_;
    tgt::vec2 texCoord_;
    tgt::ivec2 texIndex_;

    VertexColorNormalTexCoord() {}
    VertexColorNormalTexCoord(tgt::vec3 pos, tgt::vec4 color, tgt::vec3 normal, tgt::vec2 texCoord, tgt::ivec2 texIndex);

    bool equals(const VertexColorNormalTexCoord& other, double epsilon = 1e-5) const;

    static void setupVertexAttributePointers();

    static void disableVertexAttributePointers();

    static VertexColorNormalTexCoord interpolate(const VertexColorNormalTexCoord& v1, const VertexColorNormalTexCoord& v2, float t);

    void setColor(tgt::vec4 color) {color_ = color;}
    void setNormal(tgt::vec3 n) {normal_ = n;}
    void setTexCoord(tgt::vec2 texCoord) {texCoord_ = texCoord;}
    void setTexIndex(tgt::ivec2 texIndex) {texIndex_ = texIndex;}
    tgt::vec3 getNormal() {return normal_;}
    tgt::vec4 getColor() {return color_;}
    tgt::vec2 getTexCoord() {return texCoord_;}
    tgt::ivec2 getTexIndex() {return texIndex_;}
};

} // namespace

#endif  //VRN_VERTEX_H
