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

#ifndef TGT_TESSELATOR_H
#define TGT_TESSELATOR_H

#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE

#include <vector>

#include "tgt/assert.h"
#include "tgt/singleton.h"
#include "tgt/vector.h"
#include "tgt/vertex.h"

/*
    TODO Edgeflag support

    TODO Perhaps there are problems in extractTriangleMesh with the winding direction

    TODO Perhaps problems with vertexTypes where point_ is only a vec2 or even a vec4
*/

namespace tgt {

class Tesselator;
#ifdef DLL_TEMPLATE_INST
template class TGT_API Singleton<Tesselator>;
#endif

class TGT_API Tesselator : public Singleton<Tesselator> {
private:

#ifdef WIN32
    #define TGT_TESSELATORAPI __stdcall
#else //WIN32
    #define TGT_TESSELATORAPI
#endif //WIN32

    typedef void (TGT_TESSELATORAPI *fnp_gluCallback)();

    GLUtesselator* tesselator_;

    /*
        Callbacks for GLU
    */
    template<class InterleavedType>
    struct Callbacks {
        static void TGT_TESSELATORAPI beginData(GLenum type, void* polygonData );
        static void TGT_TESSELATORAPI endData(void* polygonData);
    //  static void TGT_TESSELATORAPI edgeFlagData(GLboolean flag, void* polygonData );
        static void TGT_TESSELATORAPI vertexData(void* vertexData, void* polygonData );
        static void TGT_TESSELATORAPI combine    (GLdouble coords[3], void* vertexData[4],
                                                GLfloat weight[4], void** outData);
    };

#ifdef TGT_DEBUG
    static void errorData(GLenum error, void* polygonData);
#endif //TGT_DEBUG

// Prevent namespace clashes with X11
#ifdef WindingRule
#undef WindingRule
#endif

public:

    enum WindingRule {
        ODD         = GLU_TESS_WINDING_ODD,
        NONZERO     = GLU_TESS_WINDING_NONZERO,
        NEGATIVE    = GLU_TESS_WINDING_NEGATIVE,
        POSITIVE    = GLU_TESS_WINDING_POSITIVE,
        ABS_GEQ_TWO = GLU_TESS_WINDING_ABS_GEQ_TWO
    };

    Tesselator()
      : tesselator_(gluNewTess()) {
    }

    virtual ~Tesselator() {
        gluDeleteTess(tesselator_);
    }

    /*
        getters and setters
    */
    void setWindingRule(WindingRule windingRule) {
        gluTessProperty(tesselator_, GLU_TESS_WINDING_RULE, windingRule);
    }

    WindingRule getWindingRule() const;

    void setBoundaryOnly(bool boundaryOnly) {
        double value = boundaryOnly ? GL_TRUE : GL_FALSE;
        gluTessProperty(tesselator_, GLU_TESS_BOUNDARY_ONLY, value);
    }

    bool getBoundaryOnly() const {
        GLdouble result;
        gluGetTessProperty(tesselator_, GLU_TESS_BOUNDARY_ONLY, &result);
        return result == GL_TRUE ? true : false;
    }

    void setTolerance(double value) {
        gluTessProperty(tesselator_, GLU_TESS_TOLERANCE, value);
    }

    GLdouble getTolerance() const {
        GLdouble value;
        gluGetTessProperty(tesselator_, GLU_TESS_TOLERANCE, &value);
        return value;
    }

    /*
        tesselation methods
    */

    /**
        tessalte a simple non-concave polygon
    */
    template<class InterleavedType>
    std::vector<InterleavedType> tesselate(
        const InterleavedType& vertices,
        const vec3& polygonNormal = vec3(0.f, 0.f, 0.f),
        WindingRule windingRule = ODD) const;

    /**
        tessalte more complex polygons
    */
    template<class InterleavedType>
    std::vector<InterleavedType> tesselate(
        const std::vector<InterleavedType>& contours,
        const vec3& polygonNormal = vec3(0.f, 0.f, 0.f),
        WindingRule = ODD) const;

    template<class InterleavedType>
    static InterleavedType extractTriangleMesh(const InterleavedType& polygon);

    template<class InterleavedType>
    static InterleavedType extractTriangleMesh(const std::vector<InterleavedType>& polys);
};

/*
    gluCallbacks
*/

template<class InterleavedType>
void TGT_TESSELATORAPI Tesselator::Callbacks<InterleavedType>::beginData(GLenum type, void* polygonData ) {
    //Don't try this at home
    ((std::vector<InterleavedType>*) polygonData)->push_back( InterleavedType(type) );
}

// template<class InterleavedType> TGT_TESSELATORAPI
// void Tesselator::edgeFlagData(tgtUNUSED GLboolean flag, tgtUNUSED void* polygonData ) {
//     std::cout << "hallo" << std::endl;
// }

template<class InterleavedType>
void TGT_TESSELATORAPI Tesselator::Callbacks<InterleavedType>::endData(void* polygonData) {
    std::vector<InterleavedType> contours
        = *((std::vector<InterleavedType>*) polygonData);

    for (size_t contourIdx = 0; contourIdx < contours.size(); ++contourIdx) {
        std::cout << std::endl;
        for (size_t i = 0; i < contours[contourIdx].size(); ++i)
            std::cout << contours[contourIdx][i].point_ << std::endl;
    }
}

template<class InterleavedType>
void TGT_TESSELATORAPI Tesselator::Callbacks<InterleavedType>::vertexData(void* vertexData, void* polygonData) {
    std::vector<InterleavedType>* polys = ((std::vector<InterleavedType>*) polygonData);

    polys->back().push_back(
        *((typename InterleavedType::VertexFormat*) vertexData));
}

template<class InterleavedType>
void TGT_TESSELATORAPI Tesselator::Callbacks<InterleavedType>::combine(GLdouble coords[3], void* vertexData[4],
        GLfloat weight[4], void** outData) {

    typename InterleavedType::VertexFormat** vertex    =    (typename InterleavedType::VertexFormat**) vertexData;
    typename InterleavedType::VertexFormat*  newVertex
        = (typename InterleavedType::VertexFormat*) malloc( sizeof(typename InterleavedType::VertexFormat) );

    *newVertex =  *(vertex[0]) * weight[0];
    if (weight[1] != 0)
        *newVertex = *newVertex + *(vertex[1]) * weight[1];
    if (weight[2] != 0)
        *newVertex = *newVertex + *(vertex[2]) * weight[2];
    if (weight[3] != 0)
        *newVertex = *newVertex + *(vertex[3]) * weight[3];

    for (size_t i = 0; i < newVertex->point_.size; ++i)
        newVertex->point_[i] = static_cast<float>(coords[i]);

    *outData = newVertex;
}

/*
    Tesselation methods
*/

template<class InterleavedType>
std::vector<InterleavedType> Tesselator::tesselate(
        const InterleavedType& vertices,
        const vec3& polygonNormal,
        WindingRule windingRule) const {

    std::vector<InterleavedType> polygons;

    //set polygon normal
    gluTessNormal(tesselator_, polygonNormal.x, polygonNormal.y, polygonNormal.z);
    //set winding rule
    gluTessProperty(tesselator_, GLU_TESS_WINDING_RULE, windingRule);

    //register callbacks
    gluTessCallback(tesselator_, GLU_TESS_BEGIN_DATA,     (fnp_gluCallback) Callbacks<InterleavedType>::beginData);
    gluTessCallback(tesselator_, GLU_TESS_VERTEX_DATA,    (fnp_gluCallback) Callbacks<InterleavedType>::vertexData);
    gluTessCallback(tesselator_, GLU_TESS_COMBINE,        (fnp_gluCallback) Callbacks<InterleavedType>::combine);
    gluTessCallback(tesselator_, GLU_TESS_END_DATA,       NULL);
    // gluTessCallback(tesselator_, GLU_TESS_END_DATA,       (fnp_gluCallback) Callbacks<InterleavedType>::endData);
    gluTessCallback(tesselator_, GLU_TESS_EDGE_FLAG_DATA, NULL);

#ifdef TGT_DEBUG
    gluTessCallback(tesselator_, GLU_TESS_ERROR_DATA,     (fnp_gluCallback) &errorData);
#endif //TGT_DEBUG

    const size_t size = vertices[0].point_.size;

    gluTessBeginPolygon(tesselator_, &polygons);
        gluTessBeginContour(tesselator_);

            for (size_t i = 0; i < vertices.size(); ++i) {
                //build a vec3d vector
                Vector3d vec3d(0.0, 0.0, 0.0);
                for (size_t j = 0; j < size; ++j)
                    vec3d.elem[j] = vertices[i].point_.elem[j];

                gluTessVertex(tesselator_,
                    vec3d.elem,
                    // I promise you not to change the data
                    const_cast<typename InterleavedType::VertexFormat*>(&vertices[i]));
            }

        gluTessEndContour(tesselator_);
    gluTessEndPolygon(tesselator_);

    return polygons;
}

template<class InterleavedType>
std::vector<InterleavedType> Tesselator::tesselate(
        const std::vector<InterleavedType>& contours,
        const vec3& polygonNormal,
        WindingRule windingRule) const {

    std::vector<InterleavedType> polygons;

    //set polygon normal
    gluTessNormal(tesselator_, polygonNormal.x, polygonNormal.y, polygonNormal.z);
    //set winding rule
    gluTessProperty(tesselator_, GLU_TESS_WINDING_RULE, windingRule);

    //register callbacks
    gluTessCallback(tesselator_, GLU_TESS_BEGIN_DATA,     (fnp_gluCallback) Callbacks<InterleavedType>::beginData);
/*    gluTessCallback(tesselator_, GLU_TESS_END_DATA,       (fnp_gluCallback) Callbacks<InterleavedType>::endData);*/
    gluTessCallback(tesselator_, GLU_TESS_END_DATA,       NULL);
    gluTessCallback(tesselator_, GLU_TESS_VERTEX_DATA,    (fnp_gluCallback) Callbacks<InterleavedType>::vertexData);
    gluTessCallback(tesselator_, GLU_TESS_COMBINE,        (fnp_gluCallback) Callbacks<InterleavedType>::combine);
    gluTessCallback(tesselator_, GLU_TESS_EDGE_FLAG_DATA, NULL);

#ifdef TGT_DEBUG
    gluTessCallback(tesselator_, GLU_TESS_ERROR_DATA,     (fnp_gluCallback) &errorData);
#endif //TGT_DEBUG

    gluTessBeginPolygon(tesselator_, &polygons);

    const size_t size = contours[0][0].point_.size;

        for (size_t contourIdx = 0; contourIdx < contours.size(); ++contourIdx) {
            gluTessBeginContour(tesselator_);
                for (size_t i = 0; i < contours[contourIdx].size(); ++i) {
                    //build a vec3d vector
                    Vector3d vec3d(0.0, 0.0, 0.0);
                    for (size_t j = 0; j < size; ++j)
                        vec3d.elem[j] = contours[contourIdx][i].point_.elem[j];

                    gluTessVertex(
                        tesselator_,
                        vec3d.elem,
                        // I promise you not to change the data
                        const_cast<typename InterleavedType::VertexFormat*>(&contours[contourIdx][i]));
                }
            gluTessEndContour(tesselator_);
        }
    gluTessEndPolygon(tesselator_);

    return polygons;
}

template<class InterleavedType>
InterleavedType Tesselator::extractTriangleMesh(const InterleavedType& polygon) {
    if (polygon.primitiveMode_ == GL_TRIANGLES)
        return polygon;

    InterleavedType result(GL_TRIANGLES);

    result.push_back(polygon[0]);
    result.push_back(polygon[1]);
    result.push_back(polygon[2]);

    if (polygon.primitiveMode_ == GL_TRIANGLE_STRIP) {
        for (size_t i = 3; i < polygon.size(); ++i) {
            result.push_back(polygon[i-2]);
            result.push_back(polygon[i-1]);
            result.push_back(polygon[i]);
        }
        return result;
    }

    if (polygon.primitiveMode_ == GL_TRIANGLE_FAN) {
        for (size_t i = 3; i < polygon.size(); ++i) {
            result.push_back(polygon[0]);
            result.push_back(polygon[i-1]);
            result.push_back(polygon[i]);
        }
        return result;
    }

    tgtAssert(false, "This primitive mode is not supported.");
    return polygon;
}

template<class InterleavedType>
InterleavedType Tesselator::extractTriangleMesh(const std::vector<InterleavedType>& polys) {
    InterleavedType result(GL_TRIANGLES);

    for (size_t i = 0; i < polys.size(); ++i) {
        InterleavedType temp = extractTriangleMesh(polys[i]);
        result.insert(result.end(), temp.begin(), temp.end());
    }

    return result;
}

}

#define TessMgr tgt::Singleton<tgt::Tesselator>::getRef()

#endif //VRN_OPENGL_COMPATIBILITY_PROFILE

#endif //TGT_TESSELATOR_H
