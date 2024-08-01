/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_TRIANGLEMESHGEOMETRY_H
#define VRN_TRIANGLEMESHGEOMETRY_H

#include <vector>

#include "tgt/tgt_gl.h"
#include "tgt/glmath.h"
#include "tgt/matrixstack.h"
#include "tgt/texture.h"

#include "geometry.h"
#include "vertex.h"

#include "voreen/core/io/serialization/xmlserializer.h"
#include "voreen/core/io/serialization/xmldeserializer.h"

namespace voreen {

/// Base class for all TriangleMeshGeometries.
class VRN_CORE_API TriangleMeshGeometryBase : public Geometry {
public:

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

    TriangleMeshGeometryBase();
    ~TriangleMeshGeometryBase();

    virtual size_t getNumTriangles() const = 0;
    virtual bool isEmpty() const = 0;
    virtual void clear() = 0;
    virtual VertexLayout getVertexLayout() const = 0;

    bool supportsNormals() const;
    bool supportsColors() const;
    bool supportsTextureData() const;

    virtual void calculateNormals() = 0;

    void loadTextureData();
    void loadTextureData(const std::string& filename);
    void loadTextureData(const std::vector<std::string>& filenames);
    void setTextureData(tgt::Texture* tex, bool owner = false);
    void clearTextureData();

    tgt::Texture* getTextureData() const;
    bool hasTextureData() const;

    /// Additionally serialize texture filenames
    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

    std::string getShaderDefines() const;

    /**
     * Indicates whether this geometry or parts of it are transparent.
     */
    virtual bool isTransparent() const;

    /**
     * Mark this geometry as (not) using transparent components.
     * A user can specify this per TriangleMeshGeometryBase object
     * e.g. when he/she adds a semitransparent texture.
     */
    void setHasTransparentComponents(bool hasTransparentComponents);

protected:
    void setTextureDataAndFilenames(tgt::Texture* tex, const std::vector<std::string>& filenames, bool owner = false);
    tgt::Texture* textureData_;
    std::vector<std::string> textureFilenames_;
    bool textureOwner_;
    /// Used to indicate whether this geometry contains transparent components (e.g. transparent vertex colours, transparent regions in textures etc.)
    bool hasTransparentComponents_;
};

/// Generic triangle, parameterized by the vertex-type. (@see VertexBase)
template <class V>
struct Triangle {
    V v_[3];

    typedef V VertexType;

    Triangle() {}
    Triangle(V v1, V v2, V v3)
    {
      v_[0] = v1;
      v_[1] = v2;
      v_[2] = v3;
    }
};

/*
 * Generic triangle-mesh, storing a vector of triangles.
 * Template argument is the vertex (not triangle!) type.
 */
template <class V>
class TriangleMeshGeometry : public TriangleMeshGeometryBase {
public:
    typedef Triangle<V> TriangleType;
    typedef V VertexType;

    TriangleMeshGeometry();

    /// Clears the vector and deletes the OpenGL buffer if necessary.
    virtual ~TriangleMeshGeometry();

    virtual size_t getNumTriangles() const;

    virtual bool isEmpty() const;

    virtual void clear();

    virtual void calculateNormals();

    /// Adds a triangle to the mesh. Flags the bounding box and OpenGL buffer as invalid.
    void addTriangle(const TriangleType& t);

    const TriangleType& getTriangle(size_t i) const;

    const TriangleType& operator[] (size_t i) const;

    /// Modifies a triangle of the mesh. Flags the bounding box and OpenGL buffer as invalid.
    void setTriangle(const TriangleType& t, size_t i);

    /// Returns the bounding box in model or transformed coordinates. The BB is cached internally.
    virtual tgt::Bounds getBoundingBox(bool transformed = true) const;

    /// Serializes the triangles as binary blob.
    virtual void serialize(Serializer& s) const;

    virtual void deserialize(Deserializer& s);

    /// Flags the bounding box and OpenGL buffer as invalid.
    void invalidate();

    /*
     * Renders the mesh using OpenGL by binding the buffer, setting appropriate vertex attribute pointers and calling glDrawArrays.
     *
     * A shader with appropriate vertexattribute bindings has to be activated before calling render().
     */
    virtual void render() const;

    /**
     * Returns true, if the passed Geometry is a TriangleMeshGeometry of the same type and all its vertices are equal to this one's.
     */
    virtual bool equals(const Geometry* geometry, double epsilon = 1e-5) const;

    /// Triangulates a convex polygon
    void triangulate(const std::vector<VertexType>& poly);

    /// Triangulates a quad given by four vertices and add the two triangles to the mesh.
    void addQuad(const VertexType& v1, const VertexType& v2, const VertexType& v3, const VertexType& v4);

    /// Add the triangles of another mesh to this mesh. Vertices are transformed if necessary.
    void addMesh(const TriangleMeshGeometry<V>* mesh);

    /// Clips all triangles against the plane and closes the mesh. Works only for convex meshes.
    virtual void clip(const tgt::plane& clipPlane, double epsilon = 1e-5);

    virtual std::unique_ptr<Geometry> clone() const;

protected:
    void initialize() const;

    virtual void updateBoundingBox() const;

    /// Updates the OpenGL buffer.
    virtual void updateBuffer() const;

    void clipTriangle(const TriangleType& in, std::vector<TriangleType>& out, std::vector< std::pair<VertexType, VertexType> >& edgeList, const tgt::plane& clipPlane, double epsilon = 1e-5);

    mutable tgt::Bounds boundingBox_;
    mutable bool boundingBoxValid_;

    mutable GLuint vertexArrayID_;
    mutable GLuint bufferObjectID_;
    mutable bool bufferObjectValid_;

    mutable bool isInitialized_;

    std::vector<TriangleType> triangles_;
};

template <class V>
TriangleMeshGeometry<V>::TriangleMeshGeometry()
    : boundingBox_()
    , boundingBoxValid_(false)
    , vertexArrayID_(0)
    , bufferObjectID_(0)
    , bufferObjectValid_(false)
    , isInitialized_(false)
{}

template <class V>
TriangleMeshGeometry<V>::~TriangleMeshGeometry() {
    clear();

    //clean up buffers
    if(isInitialized_) {
        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER,0);
        glDeleteVertexArrays(1, &vertexArrayID_);
        glDeleteBuffers(1, &bufferObjectID_);
    }

    vertexArrayID_ = 0;
    bufferObjectID_ = 0;
    bufferObjectValid_ = false;
}

template <class V>
void TriangleMeshGeometry<V>::initialize() const {
    tgtAssert(!isInitialized_, "initialize is called twice");
    glGenVertexArrays(1, &vertexArrayID_);
    glGenBuffers(1, &bufferObjectID_);
    isInitialized_ = true;
}

template <class V>
size_t TriangleMeshGeometry<V>::getNumTriangles() const {
    return triangles_.size();
}

template <class V>
bool TriangleMeshGeometry<V>::isEmpty() const {
    return (triangles_.size() == 0);
}

template <class V>
void TriangleMeshGeometry<V>::clear() {
    triangles_.clear();

    invalidate();
}

template <class V>
void TriangleMeshGeometry<V>::calculateNormals() {
    if(!supportsNormals())
        return;
    //TODO:
    LERROR("Not implemented yet!");
}

template <class V>
void TriangleMeshGeometry<V>::addTriangle(const TriangleType& t) {
    triangles_.push_back(t);
    invalidate();
}

template <class V>
const typename TriangleMeshGeometry<V>::TriangleType& TriangleMeshGeometry<V>::getTriangle(size_t i) const {
    tgtAssert(i < triangles_.size(), "Invalid triangle index");
    return triangles_[i];
}

template <class V>
const typename TriangleMeshGeometry<V>::TriangleType& TriangleMeshGeometry<V>::operator[] (size_t i) const {
    tgtAssert(i < triangles_.size(), "Invalid triangle index");
    return triangles_[i];
}

template <class V>
void TriangleMeshGeometry<V>::setTriangle(const TriangleType& t, size_t i) {
    tgtAssert(i < triangles_.size(), "Invalid triangle index");
    triangles_[i] = t;

    invalidate();
}

template <class V>
tgt::Bounds TriangleMeshGeometry<V>::getBoundingBox(bool transformed) const {
    if(!boundingBoxValid_)
        updateBoundingBox();

    if(transformed)
        return boundingBox_.transform(getTransformationMatrix());
    else
        return boundingBox_;
}

template <class V>
void TriangleMeshGeometry<V>::serialize(Serializer& s) const {
    TriangleMeshGeometryBase::serialize(s);
    s.serializeBinaryBlob("triangles", triangles_);
}

template <class V>
void TriangleMeshGeometry<V>::deserialize(Deserializer& s) {
    TriangleMeshGeometryBase::deserialize(s);
    s.deserializeBinaryBlob("triangles", triangles_);
}

template <class V>
void TriangleMeshGeometry<V>::invalidate() {
    boundingBoxValid_ = false;
    bufferObjectValid_ = false;
}

template <class V>
void TriangleMeshGeometry<V>::render() const {
    updateBuffer();

    if(isEmpty())
        return;

   // MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
   // MatStack.pushMatrix();
   // MatStack.multMatrix(getTransformationMatrix());
    if(!isInitialized_)
        initialize();
    glBindVertexArray(vertexArrayID_);
    glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(triangles_.size() * 3));

   // MatStack.popMatrix();
}

template <class V>
bool TriangleMeshGeometry<V>::equals(const Geometry* geometry, double epsilon) const {
    const TriangleMeshGeometry<V>* triMesh = dynamic_cast<const TriangleMeshGeometry<V>*>(geometry);

    if(triMesh) {
        if(getNumTriangles() != triMesh->getNumTriangles())
            return false;

        for(size_t i=0; i<triangles_.size(); i++) {
            if (!triangles_[i].v_[0].equals(triMesh->triangles_[i].v_[0], epsilon) ||
                !triangles_[i].v_[1].equals(triMesh->triangles_[i].v_[1], epsilon) ||
                !triangles_[i].v_[2].equals(triMesh->triangles_[i].v_[2], epsilon)  )
            {
                return false;
            }
        }
        return true;
    }
    else
        return false;
}

template <class V>
void TriangleMeshGeometry<V>::triangulate(const std::vector<VertexType>& poly) {
    if(poly.size() < 3)
        return;

    for(size_t i=2; i<poly.size(); i++)
        addTriangle(TriangleType(poly[0], poly[i-1], poly[i]));
}

template <class V>
void TriangleMeshGeometry<V>::addQuad(const VertexType& v1, const VertexType& v2, const VertexType& v3, const VertexType& v4) {
    addTriangle(TriangleType(v1, v2, v3));
    addTriangle(TriangleType(v1, v3, v4));
}

template <class V>
void TriangleMeshGeometry<V>::addMesh(const TriangleMeshGeometry<V>* mesh) {
    if(getTransformationMatrix() == mesh->getTransformationMatrix()) {
        //Just copy the triangles:
        for(size_t i=0; i<mesh->getNumTriangles(); i++)
            addTriangle(mesh->getTriangle(i));
    }
    else {
        tgt::mat4 m;
        getTransformationMatrix().invert(m);
        m = m * mesh->getTransformationMatrix();

        tgt::mat3 normalTransform;
        bool transformNormals = getVertexLayout() & NORMAL;
        if(transformNormals) {
            m.getRotationalPartMat3().invert(normalTransform);
            normalTransform = transpose(normalTransform);
        }

        for(size_t i=0; i<mesh->getNumTriangles(); i++) {
            TriangleType t = mesh->getTriangle(i);
            t.v_[0].pos_ = m * t.v_[0].pos_;
            t.v_[1].pos_ = m * t.v_[1].pos_;
            t.v_[2].pos_ = m * t.v_[2].pos_;

            if(transformNormals) {
                t.v_[0].setNormal(normalTransform * t.v_[0].getNormal());
                t.v_[1].setNormal(normalTransform * t.v_[1].getNormal());
                t.v_[2].setNormal(normalTransform * t.v_[2].getNormal());
            }
            addTriangle(t);
        }
    }
}

template <class V>
void TriangleMeshGeometry<V>::clip(const tgt::plane& clipPlane, double epsilon) {
    tgtAssert(epsilon >= 0.0, "negative epsilon");

    tgt::plane pl = clipPlane.transform(getInvertedTransformationMatrix());

    tgt::Bounds b = getBoundingBox(false);
    if(!b.intersects(pl)) {
        // The plane isn't intersecting the boundingbox, so its either left as-is or removed:
        if(pl.distance(b.center()) < 0.0) {
            return;
        }
        else {
            clear();
            return;
        }
    }

    std::vector<TriangleType> clippedTriangles;
    std::vector< std::pair<VertexType, VertexType> > edgeList;

    // Clip all faces...
    for (size_t i=0; i<getNumTriangles(); i++) {
        clipTriangle(getTriangle(i), clippedTriangles, edgeList, pl, epsilon);
    }

    std::swap(triangles_, clippedTriangles);
    invalidate();

    // Is closing necessary?
    if (edgeList.size() > 1) {
        // Sort edges to produce contiguous vertex order...
        for (size_t i = 0; i < edgeList.size() - 1; ++i) {
            VertexType connectionVertex = edgeList.at(i).second;
            for (size_t j = i + 1; j < edgeList.size(); ++j) {
                if (distance(edgeList.at(j).first.pos_, connectionVertex.pos_) < epsilon) {
                    std::swap(edgeList.at(i + 1), edgeList.at(j));
                    break;
                }
                else if (distance(edgeList.at(j).second.pos_, connectionVertex.pos_) < epsilon) {
                    VertexType tmp = edgeList.at(j).first;
                    edgeList.at(j).first = edgeList.at(j).second;
                    edgeList.at(j).second = tmp;
                    std::swap(edgeList.at(i + 1), edgeList.at(j));
                    break;
                }
            }
        }
        // Set normal of all vertices to plane normal:
        for (size_t i = 0; i < edgeList.size(); ++i) {
            edgeList.at(i).first.setNormal(pl.n);
            edgeList.at(i).second.setNormal(pl.n);
        }

        // for concavity test: find best suited projection plane
        int index = (tgt::maxElem(tgt::abs(pl.n)) + 1) % 3;
        float signum = 0.f;
        bool foundConcavity = false;

        // Convert sorted edge list to sorted vertex list. Test for concavities in the resulting polygon.
        std::vector<VertexType> closingFaceVertices;
        for (size_t i = 0; i < edgeList.size(); ++i) {
            // find vec2 perpendicular to projection of edge onto the selected plane
            tgt::vec2 perp = normalize(tgt::vec2(edgeList.at(i).first.pos_[(index + 1) % 3] - edgeList.at(i).second.pos_[(index + 1) % 3], edgeList.at(i).second.pos_[index] - edgeList.at(i).first.pos_[index]));
            tgt::vec3 controlVec = tgt::vec3(edgeList.at((i+1) % edgeList.size()).first.pos_ - edgeList.at(i).first.pos_);
            if(signum == 0.f && i != edgeList.size() - 1)
                signum = tgt::sign(tgt::dot(perp, tgt::vec2(controlVec[index], controlVec[(index + 1) % 3])));

            if(signum != 0.f) {
                float res = tgt::dot(perp, tgt::vec2(controlVec[index], controlVec[(index + 1) % 3]));
                if(signum * res < 0.f && std::abs(res) > epsilon)
                    foundConcavity = true;
            }

            if (i == 0)
                closingFaceVertices.push_back(VertexType::interpolate(edgeList.at(edgeList.size() - 1).second, edgeList.at(i).first, 0.5f));
            else
                closingFaceVertices.push_back(VertexType::interpolate(edgeList.at(i - 1).second, edgeList.at(i).first, 0.5f));
        }

        if(foundConcavity)
            LWARNINGC("trianglemeshgeometry.clip", "Trying to clip concave geometry, may result in undefined behaviour. Possible reason: optimized proxy geometry with optimization level >= visible bricks.");

        // Convert vertex order to counter clockwise if necessary...
        tgt::vec3 closingFaceNormal(0, 0, 0);
        for (size_t i = 0; i < closingFaceVertices.size(); ++i)
            closingFaceNormal += tgt::cross(closingFaceVertices[i].pos_, closingFaceVertices[(i + 1) % closingFaceVertices.size()].pos_);

        closingFaceNormal = tgt::normalize(closingFaceNormal);

        if (tgt::dot(pl.n, closingFaceNormal) < 0)
            std::reverse(closingFaceVertices.begin(), closingFaceVertices.end());

        if(closingFaceVertices.size() > 2) {
            triangulate(closingFaceVertices);
        }
    }
}

template <class V>
void TriangleMeshGeometry<V>::updateBoundingBox() const {
    tgt::Bounds bb;
    for(size_t i=0; i<triangles_.size(); i++) {
        bb.addPoint(triangles_[i].v_[0].pos_);
        bb.addPoint(triangles_[i].v_[1].pos_);
        bb.addPoint(triangles_[i].v_[2].pos_);
    }
    boundingBox_ = bb;
    boundingBoxValid_ = true;
}

template <class V>
void TriangleMeshGeometry<V>::updateBuffer() const {
    if(bufferObjectValid_)
        return;

    if(isEmpty())
        return;

    if(!isInitialized_)
        initialize();

    glBindVertexArray(vertexArrayID_);
    glBindBuffer(GL_ARRAY_BUFFER, bufferObjectID_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(TriangleType) * triangles_.size(), &triangles_[0], GL_STATIC_DRAW);
    VertexType::setupVertexAttributePointers();
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    bufferObjectValid_ = true;
}

template <class V>
void TriangleMeshGeometry<V>::clipTriangle(const TriangleType& in, std::vector<TriangleType>& out, std::vector< std::pair<VertexType, VertexType> >& edgeList, const tgt::plane& clipPlane, double epsilon) {
    tgtAssert(epsilon >= 0.0, "negative epsilon");

    double lastDistance = in.v_[0].getDistanceToPlane(clipPlane, epsilon);

    std::vector<VertexType> outVertices;
    std::vector<VertexType> onPlaneVertices;

    // Process face edges...
    for (size_t i = 0; i < 3; ++i) {
        double distance = in.v_[(i + 1) % 3].getDistanceToPlane(clipPlane, epsilon);

        // Keep both vertices?
        if (lastDistance <= 0 && distance <= 0) {
            // If processing the first edge, insert first vertex...
            if (i == 0)
                outVertices.push_back(in.v_[i]);

            // If NOT processing the last edge, insert second vertex...
            if (i < 2)
                outVertices.push_back(in.v_[i + 1]);
        }
        // Discard first vertex, but keep second vertex?
        else if (lastDistance > 0 && distance <= 0) {
            // If NOT clipplane intersection vertex and second vertex are equal, insert clipplane intersection vertex...
            VertexType intersectionVertex = VertexType::interpolate(in.v_[i], in.v_[(i + 1) % 3], static_cast<float>(lastDistance / (lastDistance - distance)));
            if (!in.v_[(i + 1) % 3].equals(intersectionVertex, epsilon)) {
                outVertices.push_back(intersectionVertex);
                onPlaneVertices.push_back(intersectionVertex);
            }
            else
                onPlaneVertices.push_back(in.v_[(i + 1) % 3]);


            // If NOT processing the last edge, insert second vertex...
            if (i < 2)
                outVertices.push_back(in.v_[i + 1]);
        }
        // Keep first vertex, but discard second vertex?
        else if (lastDistance <= 0 && distance > 0) {
            // If processing the first edge, insert first vertex...
            if (i == 0)
                outVertices.push_back(in.v_[i]);

            //// If NOT clipplane intersection vertex and first vertex are equal, insert clipplane intersection vertex...
            VertexType intersectionVertex = VertexType::interpolate(in.v_[i], in.v_[(i + 1) % 3], static_cast<float>(lastDistance / (lastDistance - distance)));
            if (!in.v_[i].equals(intersectionVertex, epsilon)) {
                outVertices.push_back(intersectionVertex);
                onPlaneVertices.push_back(intersectionVertex);
            }
            else
                onPlaneVertices.push_back(in.v_[i]);
        }

        lastDistance = distance;
    }

    // Create triangles from output vertices:
    if(outVertices.size() == 3) {
        out.push_back(TriangleType(outVertices[0], outVertices[1], outVertices[2]));
    }
    else if(outVertices.size() == 4) {
        out.push_back(TriangleType(outVertices[0], outVertices[1], outVertices[2]));
        out.push_back(TriangleType(outVertices[0], outVertices[2], outVertices[3]));
    }

    if(onPlaneVertices.size() == 2) {
        if(distance(onPlaneVertices[0].pos_, onPlaneVertices[1].pos_) > epsilon)
            edgeList.push_back(std::pair<VertexType, VertexType>(onPlaneVertices[0], onPlaneVertices[1]));
    }
    else if(onPlaneVertices.size() == 0) {
    }
    else {
        tgtAssert(false, "Should not come here.");
    }
}

template<class V>
std::unique_ptr<Geometry> TriangleMeshGeometry<V>::clone() const {
    TriangleMeshGeometry<V>* newGeom = static_cast<TriangleMeshGeometry<V>*>(create());
    newGeom->triangles_ = triangles_;
    newGeom->setTransformationMatrix(getTransformationMatrix());
    newGeom->getMetaDataContainer() = MetaDataContainer(getMetaDataContainer());
    newGeom->setTextureDataAndFilenames(textureData_, textureFilenames_, false);
    return std::unique_ptr<Geometry>(newGeom);
}

//-------------------------------------------------------------------------------------------------

/// A triangle mesh with vertex type VertexNormal.
class VRN_CORE_API TriangleMeshGeometryNormal : public TriangleMeshGeometry<VertexNormal> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryNormal(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryNormal"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return NORMAL; }

    /// Creates a mesh containing a cube specified by two vertices.
    static TriangleMeshGeometryNormal* createCube(VertexType llfVertex, VertexType urbVertex);

    /// Adds a cube to this mesh.
    void addCube(VertexType llfVertex, VertexType urbVertex);
protected:
};

//-------------------------------------------------------------------------------------------------

/// A triangle mesh with vertex type VertexColorNormal
class VRN_CORE_API TriangleMeshGeometryColorNormal : public TriangleMeshGeometry<VertexColorNormal> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryColorNormal(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryColorNormal"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return COLOR_NORMAL; }

    static TriangleMeshGeometryColorNormal* createCube(tgt::vec3 coordLlf, tgt::vec3 coordUrb, tgt::vec3 colorLlf, tgt::vec3 colorUrb, float alpha, tgt::vec3 texLlf, tgt::vec3 texUrb);

    /// Adds a cube to this mesh.
    void addCube(VertexNormal llfVertex, VertexNormal urbVertex);

    void addQuad(const VertexNormal& v1, const VertexNormal& v2, const VertexNormal& v3, const VertexNormal& v4);

    /// Creates a cube with color and normals:
    static TriangleMeshGeometryColorNormal* createCube(tgt::vec3 coordLlf, tgt::vec3 coordUrb, tgt::vec3 colorLlf, tgt::vec3 colorUrb, float alpha);
};

//-------------------------------------------------------------------------------------------------

/// A triangle mesh with vertex type VertexBase (i.e., only position)
class VRN_CORE_API TriangleMeshGeometrySimple : public TriangleMeshGeometry<VertexBase> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometrySimple(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometrySimple"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return SIMPLE; }

    static TriangleMeshGeometrySimple* createCube(tgt::vec3 coordLlf, tgt::vec3 coordUrb);
};

//-------------------------------------------------------------------------------------------------

/// A triangle mesh with vertex type VertexColor
class VRN_CORE_API TriangleMeshGeometryColor : public TriangleMeshGeometry<VertexColor> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryColor(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryColor"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return COLOR; }

    /// Creates a mesh containing a cube specified by two vertices.
    static TriangleMeshGeometryColor* createCube(tgt::vec3 coordLlf, tgt::vec3 coordUrb, tgt::vec3 color, float alpha);
};

//-------------------------------------------------------------------------------------------------

/// A triangle mesh with vertex type VertexTexCoord
class VRN_CORE_API TriangleMeshGeometryTexCoord : public TriangleMeshGeometry<VertexTexCoord> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryTexCoord(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryTexCoord"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return TEXCOORD; }
};

//-------------------------------------------------------------------------------------------------

/// A triangle mesh with vertex type VertexColorTexCoord
class VRN_CORE_API TriangleMeshGeometryColorTexCoord : public TriangleMeshGeometry<VertexColorTexCoord> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryColorTexCoord(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryColorTexCoord"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return COLOR_TEXCOORD; }
};

//-------------------------------------------------------------------------------------------------

/// A triangle mesh with vertex type VertexNormalTexCoord
class VRN_CORE_API TriangleMeshGeometryNormalTexCoord : public TriangleMeshGeometry<VertexNormalTexCoord> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryNormalTexCoord(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryNormalTexCoord"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return NORMAL_TEXCOORD; }
};

//-------------------------------------------------------------------------------------------------

/// A triangle mesh with vertex type VertexColorNormalTexCoord
class VRN_CORE_API TriangleMeshGeometryColorNormalTexCoord : public TriangleMeshGeometry<VertexColorNormalTexCoord> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryColorNormalTexCoord(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryColorNormalTexCoord"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return COLOR_NORMAL_TEXCOORD; }
};

} // namespace

#endif  //VRN_TRIANGLEMESHGEOMETRY_H
