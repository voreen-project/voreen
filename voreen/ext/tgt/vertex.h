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

#ifndef TGT_VERTEX_H
#define TGT_VERTEX_H


#include <vector>

#include "tgt/tgt_gl.h"
#include "tgt/vector.h"

//TODO inherit from tgt::Renderable

namespace tgt {

/*
    Basic Vertices
*/

struct VertexV2F {
    vec2 point_;

    VertexV2F() {}

    VertexV2F(const vec2& point)
      : point_(point)
    {}

    VertexV2F operator * (float factor) const {
        return VertexV2F( point_ * factor );
    }

    VertexV2F operator + (const VertexV2F& v) const {
        return VertexV2F( point_ + v.point_);
    }
};

struct VertexV3F {
    vec3 point_;

    VertexV3F(const vec3& point)
      : point_(point)
    {}

    VertexV3F() {};

    VertexV3F operator * (float factor) const {
        return VertexV3F( point_ * factor );
    }

    VertexV3F operator + (const VertexV3F& v) const {
        return VertexV3F( point_ + v.point_);
    }
};

struct VertexC4UB_V2F {
    col4 color_;
    vec2 point_;

    VertexC4UB_V2F(const col4& color, const vec2& point)
      : color_(color),
        point_(point)
    {}
};

struct VertexC4UB_V3F {
    col4 color_;
    vec3 point_;

    VertexC4UB_V3F(const col4& color, const vec3& point)
      : color_(color),
        point_(point)
    {}
};

struct VertexC3F_V3F {
    vec3 color_;
    vec3 point_;

    VertexC3F_V3F(const vec3& color, const vec3& point)
      : color_(color),
        point_(point)
    {}
};

struct VertexN3F_V3F {
    vec3 normal_;
    vec3 point_;

    VertexN3F_V3F(const vec3& normal, const vec3& point)
      : normal_(normal),
        point_(point)
    {}
};

struct VertexC4F_N3F_V3F {
    col4 color_;
    vec3 normal_;
    vec3 point_;

    VertexC4F_N3F_V3F(const col4& color, const vec3& normal, const vec3& point)
      : color_(color),
        normal_(normal),
        point_(point)
    {}
};

struct VertexT2F_V3F {
    vec2 texCoord_;
    vec3 point_;
};

struct VertexT4F_V4F {
    vec4 texCoord_;
    vec4 point_;
};

struct VertexT2F_C4UB_V3F {
    vec2 texCoords_;
    col4 color_;
    vec3 point_;
};

struct VertexT2F_C3F_V3F {
    vec2 texCoord_;
    vec3 color_;
    vec3 point_;
};

struct VertexT2F_N3F_V3F {
    vec2 texCoord_;
    vec3 normal_;
    vec3 point_;

    VertexT2F_N3F_V3F() {}

    VertexT2F_N3F_V3F(const vec2& texCoord, const vec3& normal, const vec3& point)
      : texCoord_(texCoord),
        normal_(normal),
        point_(point)
    {}

    VertexT2F_N3F_V3F operator * (float factor) const {
        return VertexT2F_N3F_V3F( texCoord_ * factor, normal_ * factor, point_ * factor );
    }

    VertexT2F_N3F_V3F operator + (const VertexT2F_N3F_V3F& v) const {
        return VertexT2F_N3F_V3F( texCoord_ + v.texCoord_, normal_ + v.normal_, point_ + v.point_);
    }
};

struct VertexT2F_C4F_N3F_V3F {
    vec2 texCoord_;
    vec4 color_;
    vec3 normal_;
    vec3 point_;

    VertexT2F_C4F_N3F_V3F(const vec2& texCoord, const vec4& color, const vec3& normal, const vec3& point)
      : texCoord_(texCoord),
        color_(color),
        normal_(normal),
        point_(point)
    {}


    VertexT2F_C4F_N3F_V3F operator * (float factor) const {
        return VertexT2F_C4F_N3F_V3F( texCoord_ * factor, color_, normal_ * factor, point_ * factor );
    }

    VertexT2F_C4F_N3F_V3F operator + (const VertexT2F_C4F_N3F_V3F& v) const {
        return VertexT2F_C4F_N3F_V3F( texCoord_ + v.texCoord_, color_, normal_ + v.normal_, point_ + v.point_);
    }

};

struct VertexT4F_C4F_N3F_V4F {
    vec4 texCoord_;
    vec4 color_;
    vec3 normal_;
    vec4 point_;

    VertexT4F_C4F_N3F_V4F(const vec4& texCoord, const vec4& color, const vec3& normal, const vec4& point)
      : texCoord_(texCoord),
        color_(color),
        normal_(normal),
        point_(point)
    {}
};

/*
    Interleaved
*/

struct InterleavedV2F : public std::vector<VertexV2F> {
    GLenum primitiveMode_;

    enum {
        FORMAT = GL_V2F
    };

    typedef VertexV2F VertexFormat;

    InterleavedV2F(GLenum primitiveMode = GL_TRIANGLES)
      : std::vector<VertexV2F>(),
        primitiveMode_(primitiveMode)
    {}

    void render() {
        glInterleavedArrays(FORMAT, 0, &operator[](0));
        glDrawArrays(primitiveMode_, 0, static_cast<GLsizei>(size()));
        glDisableClientState(GL_VERTEX_ARRAY);
    }
};

struct InterleavedV3F : public std::vector<VertexV3F> {
    GLenum primitiveMode_;

    enum {
        FORMAT = GL_V3F
    };

    typedef VertexV3F VertexFormat;

    InterleavedV3F(GLenum primitiveMode = GL_TRIANGLES)
      : std::vector<VertexV3F>(),
        primitiveMode_(primitiveMode)
    {}

    void render() {
        glInterleavedArrays(FORMAT, 0, &operator[](0));
        glDrawArrays(primitiveMode_, 0, static_cast<GLsizei>(size()));
        glDisableClientState(GL_VERTEX_ARRAY);
    }
};

struct InterleavedC4UB_V2F : public std::vector<VertexC4UB_V2F> {
    GLenum primitiveMode_;

    enum {
        FORMAT = GL_C4UB_V2F
    };

    typedef VertexC4UB_V2F VertexFormat;

    InterleavedC4UB_V2F(GLenum primitiveMode = GL_TRIANGLES)
      : std::vector<VertexC4UB_V2F>(),
        primitiveMode_(primitiveMode)
    {}

    void render() {
        glInterleavedArrays(FORMAT, 0, &operator[](0));
        glDrawArrays(primitiveMode_, 0, static_cast<GLsizei>(size()));
        glDisableClientState(GL_COLOR_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
};


struct InterleavedC4UB_V3F : public std::vector<VertexC4UB_V3F> {
    GLenum primitiveMode_;

    enum {
        FORMAT = GL_C4UB_V3F
    };

    typedef VertexC4UB_V3F VertexFormat;

    InterleavedC4UB_V3F(GLenum primitiveMode = GL_TRIANGLES)
      : std::vector<VertexC4UB_V3F>(),
        primitiveMode_(primitiveMode)
    {}

    void render() {
        glInterleavedArrays(FORMAT, 0, &operator[](0));
        glDrawArrays(primitiveMode_, 0, static_cast<GLsizei>(size()));
        glDisableClientState(GL_COLOR_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
};

struct InterleavedC3F_V3F : public std::vector<VertexC3F_V3F> {
    GLenum primitiveMode_;

    enum {
        FORMAT = GL_C3F_V3F
    };

    typedef VertexC3F_V3F VertexFormat;

    InterleavedC3F_V3F(GLenum primitiveMode = GL_TRIANGLES)
      : std::vector<VertexC3F_V3F>(),
        primitiveMode_(primitiveMode)
    {}

    void render() {
        glInterleavedArrays(FORMAT, 0, &operator[](0));
        glDrawArrays(primitiveMode_, 0, static_cast<GLsizei>(size()));
        glDisableClientState(GL_COLOR_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
};

struct InterleavedN3F_V3F : public std::vector<VertexN3F_V3F> {
    GLenum primitiveMode_;

    enum {
        FORMAT = GL_N3F_V3F
    };

    typedef VertexN3F_V3F VertexFormat;

    InterleavedN3F_V3F(GLenum primitiveMode = GL_TRIANGLES)
      : std::vector<VertexN3F_V3F>(),
        primitiveMode_(primitiveMode)
    {}

    void render() {
        glInterleavedArrays(FORMAT, 0, &operator[](0));
        glDrawArrays(primitiveMode_, 0, static_cast<GLsizei>(size()));
        glDisableClientState(GL_NORMAL_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
};

struct InterleavedC4F_N3F_V3F : public std::vector<VertexC4F_N3F_V3F> {
    GLenum primitiveMode_;

    enum {
        FORMAT = GL_C4F_N3F_V3F
    };

    typedef VertexC4F_N3F_V3F VertexFormat;

    InterleavedC4F_N3F_V3F(GLenum primitiveMode = GL_TRIANGLES)
      : std::vector<VertexC4F_N3F_V3F>(),
        primitiveMode_(primitiveMode)
    {}

    void render() {
        glInterleavedArrays(FORMAT, 0, &operator[](0));
        glDrawArrays(primitiveMode_, 0, static_cast<GLsizei>(size()));
        glDisableClientState(GL_COLOR_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
};

struct InterleavedT2F_V3F : public std::vector<VertexT2F_V3F> {
    GLenum primitiveMode_;

    enum {
        FORMAT = GL_T2F_V3F
    };

    typedef VertexT2F_V3F VertexFormat;

    InterleavedT2F_V3F(GLenum primitiveMode = GL_TRIANGLES)
      : std::vector<VertexT2F_V3F>(),
        primitiveMode_(primitiveMode)
    {}

    void render() {
        glInterleavedArrays(FORMAT, 0, &operator[](0));
        glDrawArrays(primitiveMode_, 0, static_cast<GLsizei>(size()));
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
};

struct InterleavedT4F_V4F : public std::vector<VertexT4F_V4F> {
    GLenum primitiveMode_;

    enum {
        FORMAT = GL_T4F_V4F
    };

    typedef VertexT4F_V4F VertexFormat;

    InterleavedT4F_V4F(GLenum primitiveMode = GL_TRIANGLES)
      : std::vector<VertexT4F_V4F>(),
        primitiveMode_(primitiveMode)
    {}

    void render() {
        glInterleavedArrays(FORMAT, 0, &operator[](0));
        glDrawArrays(primitiveMode_, 0, static_cast<GLsizei>(size()));
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
};

struct InterleavedT2F_C4UB_V3F : public std::vector<VertexT2F_C4UB_V3F> {
    GLenum primitiveMode_;

    enum {
        FORMAT = GL_T2F_C4UB_V3F
    };

    typedef VertexT2F_C4UB_V3F VertexFormat;

    InterleavedT2F_C4UB_V3F(GLenum primitiveMode = GL_TRIANGLES)
      : std::vector<VertexT2F_C4UB_V3F>(),
        primitiveMode_(primitiveMode)
    {}

    void render() {
        glInterleavedArrays(FORMAT, 0, &operator[](0));
        glDrawArrays(primitiveMode_, 0, static_cast<GLsizei>(size()));
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
        glDisableClientState(GL_COLOR_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
};

struct InterleavedT2F_C3F_V3F : public std::vector<VertexT2F_C3F_V3F> {
    GLenum primitiveMode_;

    enum {
        FORMAT = GL_T2F_C3F_V3F
    };

    typedef VertexT2F_C3F_V3F VertexFormat;

    InterleavedT2F_C3F_V3F(GLenum primitiveMode = GL_TRIANGLES)
      : std::vector<VertexT2F_C3F_V3F>(),
        primitiveMode_(primitiveMode)
    {}

    void render() {
        glInterleavedArrays(FORMAT, 0, &operator[](0));
        glDrawArrays(primitiveMode_, 0, static_cast<GLsizei>(size()));
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
        glDisableClientState(GL_COLOR_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
};

struct InterleavedT2F_N3F_V3F : public std::vector<VertexT2F_N3F_V3F> {
    GLenum primitiveMode_;

    enum {
        FORMAT = GL_T2F_N3F_V3F
    };

    typedef VertexT2F_N3F_V3F VertexFormat;


    InterleavedT2F_N3F_V3F(GLenum primitiveMode = GL_TRIANGLES)
      : std::vector<VertexT2F_N3F_V3F>(),
        primitiveMode_(primitiveMode)
    {}

    void render() {
        glInterleavedArrays(FORMAT, 0, &operator[](0));
        glDrawArrays(primitiveMode_, 0, static_cast<GLsizei>(size()));
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
};

struct InterleavedT2F_C4F_N3F_V3F : public std::vector<VertexT2F_C4F_N3F_V3F> {
    GLenum primitiveMode_;

    enum {
        FORMAT = GL_T2F_C4F_N3F_V3F
    };

    typedef VertexT2F_C4F_N3F_V3F VertexFormat;


    InterleavedT2F_C4F_N3F_V3F(GLenum primitiveMode = GL_TRIANGLES)
      : std::vector<VertexT2F_C4F_N3F_V3F>(),
        primitiveMode_(primitiveMode)
    {}

    void render() {
        glInterleavedArrays(FORMAT, 0, &operator[](0));
        glDrawArrays(primitiveMode_, 0, static_cast<GLsizei>(size()));
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
        glDisableClientState(GL_COLOR);
        glDisableClientState(GL_NORMAL_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
};

struct InterleavedT4F_C4F_N3F_V4F : public std::vector<VertexT4F_C4F_N3F_V4F> {
    GLenum primitiveMode_;

    enum {
        FORMAT = GL_T4F_C4F_N3F_V4F
    };

    typedef VertexT4F_C4F_N3F_V4F VertexFormat;


    InterleavedT4F_C4F_N3F_V4F(GLenum primitiveMode = GL_TRIANGLES)
      : std::vector<VertexT4F_C4F_N3F_V4F>(),
        primitiveMode_(primitiveMode)
    {}

    void render() {
        glInterleavedArrays(FORMAT, 0, &operator[](0));
        glDrawArrays(primitiveMode_, 0, static_cast<GLsizei>(size()));
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
        glDisableClientState(GL_COLOR);
        glDisableClientState(GL_NORMAL_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
};

/*
    VertexData
*/
/*

//Do we need this?

struct VertexDataV2F {
    GLenum              primitiveMode_;
    std::vector<vec2>   points_;
};

struct VertexDataV3F {
    GLenum              primitiveMode_;
    std::vector<vec3>   points_;
};

struct VertexDataC4UB_V2F {
    GLenum              primitiveMode_;
    std::vector<col4>   colors_;
    std::vector<vec2>   points_;
};

struct VertexDataC4UB_V3F {
    GLenum              primitiveMode_;
    std::vector<col4>   colors_;
    std::vector<vec3>   points_;
};

struct VertexDataC3F_V3F {
    GLenum              primitiveMode_;
    std::vector<vec3>   colors_;
    std::vector<vec3>   points_;
};

struct VertexDataN3F_V3F {
    GLenum              primitiveMode_;
    std::vector<vec3>   normals_;
    std::vector<vec3>   points_;
};

struct VertexDataC4F_N3F_V3F {
    GLenum              primitiveMode_;
    std::vector<col4>   colors_;
    std::vector<vec3>   normals_;
    std::vector<vec3>   points_;
};

struct VertexDataT2F_V3F {
    GLenum              primitiveMode_;
    std::vector<vec2>   texCoords_;
    std::vector<vec3>   points_;
};

struct VertexDataT4F_V4F {
    GLenum              primitiveMode_;
    std::vector<vec4>   texCoords_;
    std::vector<vec4>   points_;
};

struct VertexDataT2F_C4UB_V3F {
    GLenum              primitiveMode_;
    std::vector<vec2>   texCoords_;
    std::vector<col4>   colors_;
    std::vector<vec3>   points_;
};

struct VertexDataT2F_C3F_V3F {
    GLenum     primitiveMode_;
    std::vector<vec2> texCoords_;
    std::vector<vec3> colors_;
    std::vector<vec3> points_;
};

struct VertexDataT2F_N3F_V3F {
    GLenum              primitiveMode_;
    std::vector<vec2>   texCoords_;
    std::vector<vec3>   normals_;
    std::vector<vec3>   points_;
};

struct VertexDataT2F_C4F_N3F_V3F {
    GLenum              primitiveMode_;
    std::vector<vec2>   texCoords_;
    std::vector<vec4>   colors_;
    std::vector<vec3>   normals_;
    std::vector<vec3>   points_;
};

struct VertexDataT4F_C4F_N3F_V4F {
    GLenum              primitiveMode_;
    std::vector<vec4>   texCoords_;
    std::vector<vec4>   colors_;
    std::vector<vec3>   normals_;
    std::vector<vec4>   points_;
};
*/
}

#endif //tgtVERTEX_H
