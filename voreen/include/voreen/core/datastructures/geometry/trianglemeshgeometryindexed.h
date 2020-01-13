/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_TRIANGLEMESHGEOMETRYINDEXED_H
#define VRN_TRIANGLEMESHGEOMETRYINDEXED_H

#include <vector>

#include "tgt/tgt_gl.h"
#include "tgt/tgt_math.h"
#include "tgt/glmath.h"

#include "trianglemeshgeometry.h"
#include "vertex.h"

#include "voreen/core/io/serialization/xmlserializer.h"
#include "voreen/core/io/serialization/xmldeserializer.h"

namespace voreen {

    /**
     * Base class of indexed geometry
     */
class VRN_CORE_API TriangleMeshGeometryIndexedBase : public TriangleMeshGeometryBase {
public:
    enum IndexType {
        INDEX_UINT16,
        INDEX_UINT32
    };

    virtual IndexType getIndexType() const = 0;
};

/*
 * Use this class to render indexed vertex meshes. It is faster to render large amounts of trianlges sharing
 * vertices with this class, but clipping is not implemented for now.
 */
template <class I, class V>
class VRN_CORE_API TriangleMeshGeometryIndexed : public TriangleMeshGeometryIndexedBase {
public:
    typedef Triangle<I> TriangleType;
    typedef V VertexType;

    /** SimpleGeometryQuality */
    enum SimpleGeometryQuality {
        QUALITY_WORST = 0,
        QUALITY_VERY_POOR = 1,
        QUALITY_POOR = 2,
        QUALITY_NORMAL = 3,
        QUALITY_GOOD = 4,
        QUALITY_VERY_GOOD = 5,
        QUALITY_PEREFCT = 6
    };

    TriangleMeshGeometryIndexed();

    /// Clears the vector and deletes the OpenGL buffers if necessary.
    virtual ~TriangleMeshGeometryIndexed();

    virtual size_t getNumTriangles() const;
    virtual bool isEmpty() const;
    virtual void clear();

    virtual void calculateNormals();

    /// Adds a triangle to the mesh by passing three indices. Flags the bounding box and OpenGL buffer as invalid.
    void addTriangle(const tgt::Vector3<I>& t);
    void addTriangle(const TriangleType& t);

    const TriangleType& getTriangle(size_t i) const;

    const TriangleType& operator[] (size_t i) const;

    /// Modifies a triangle of the mesh. Flags the bounding box and OpenGL buffer as invalid.
    void setTriangle(const TriangleType& t, size_t i);

    /// Returns the bounding box in model or transformed coordinates. The BB is cached internally.
    virtual tgt::Bounds getBoundingBox(bool transformed = true) const;

    /// Serializes the triangles and indices as binary blob.
    virtual void serialize(Serializer& s) const;

    virtual void deserialize(Deserializer& s);

    /// Flags the bounding box and OpenGL buffer as invalid.
    void invalidate();

    /*
     * Renders the mesh using OpenGL by binding the buffer, setting appropriate vertex attribute pointers and calling glDrawElements.
     *
     * A shader with appropriate vertexattribute bindings has to be activated before calling render().
     */
    virtual void render() const;

    /**
     * Returns true, if the passed Geometry is a TriangleMeshGeometry of the same type and all its vertices are equal to this one's.
     */
    virtual bool equals(const Geometry* geometry, double epsilon = 1e-5) const;

    /// Triangulates a convex polygon
    void triangulate(const std::vector<I>& poly);

    /// Triangulates a quad given by four vertices and add the two triangles to the mesh.
    void addQuad(const I& v1, const I& v2, const I& v3, const I& v4);

    /// Add the triangles of another mesh to this mesh. Vertices are transformed if necessary.
    void addMesh(const TriangleMeshGeometryIndexed<I, V>* mesh);

    /** Simple geometries */
    void addQuadGeometry(float xWidth = 1.f, float yHeight = 1.f, tgt::vec3 center = tgt::vec3::zero, tgt::vec4 color = tgt::vec4::one);
    void addCuboidGeometry(float xWidth = 1.f, float yHeight = 1.f, float zDepth = 1.f, tgt::vec3 center = tgt::vec3::zero, tgt::vec4 color = tgt::vec4::one);
    void addDiskGeometry(float innerRadius = 0.f, float outerRadius = 1.f, tgt::vec3 center = tgt::vec3::zero, tgt::vec4 color = tgt::vec4::one, SimpleGeometryQuality quality = QUALITY_NORMAL);
    void addCylinderGeometry(float bottomRadius = 1.f, float topRadius = 1.f, float length = 1.f, tgt::vec3 center = tgt::vec3::zero, tgt::vec4 color = tgt::vec4::one, SimpleGeometryQuality quality = QUALITY_NORMAL);

    // TODO: implement actual sphere gometry, this does only generate an octahedron
    void addSphereGeometry(float radius = 1.f, tgt::vec3 center = tgt::vec3::zero, tgt::vec4 color = tgt::vec4::one, SimpleGeometryQuality quality = QUALITY_NORMAL);
    //void addSphereGeometry(float radius, tgt::vec3 center, tgt::vec4 color, size_t numSlicesAndTiles);

    void addLetterX(float height, tgt::vec4 color = tgt::vec4::one);
    void addLetterY(float height, tgt::vec4 color = tgt::vec4::one);
    void addLetterZ(float height, tgt::vec4 color = tgt::vec4::one);
private:
    void addVertex(size_t index, tgt::vec3 pos, tgt::vec3 normal, tgt::vec4 color = tgt::vec4::one, tgt::vec2 texCoord = tgt::vec2::zero, tgt::ivec2 texIndex = tgt::ivec2::zero);
public:
    virtual std::unique_ptr<Geometry> clone() const;

    /// Used to pass the vertex set for this geometry. For now, this method _must_ be called before any triangles are added.
    void setVertices(const std::vector<VertexType>& vertices) {
        vertices_ = vertices;
    }

    const std::vector<VertexType>& getVertices() const {
        return vertices_;
    }

protected:

    void initialize() const;
    virtual void updateBoundingBox() const;

    /// Updates the OpenGL buffer.
    virtual void updateBuffer() const;

    mutable tgt::Bounds boundingBox_;
    mutable bool boundingBoxValid_;

    mutable GLuint vertexArrayID_;
    mutable GLuint bufferObjectID_;
    mutable bool bufferObjectValid_;

    mutable GLuint indexObject_;
    mutable bool indexObjectValid_;

    mutable bool isInitialized_;

    std::vector<TriangleType> triangles_;
    std::vector<VertexType> vertices_;
    GLuint indexDataType_;

    static const std::string loggerCat_;
};


template<class I, class V>
const std::string TriangleMeshGeometryIndexed<I,V>::loggerCat_ = "voreen.TriangleMeshGeometryIndexed";

template <class I, class V>
TriangleMeshGeometryIndexed<I, V>::TriangleMeshGeometryIndexed()
    : boundingBox_()
    , boundingBoxValid_(false)
    , bufferObjectID_(0)
    , vertexArrayID_(0)
    , bufferObjectValid_(false)
    , indexObject_(0)
    , indexObjectValid_(false)
    , isInitialized_(false)
{
    switch(sizeof(I)) {
        case 1:
            indexDataType_ = GL_UNSIGNED_BYTE;
            break;
        case 2:
            indexDataType_ = GL_UNSIGNED_SHORT;
            break;
        case 4:
            indexDataType_ = GL_UNSIGNED_INT;
            break;
        default:
            tgtAssert(false, "Unsupported index datatype in TriangleMeshGeometryIndexed!");
    }
}

template <class I, class V>
void TriangleMeshGeometryIndexed<I, V>::initialize() const {
    tgtAssert(!isInitialized_, "initialize is called twice");
    glGenVertexArrays(1, &vertexArrayID_);
    glGenBuffers(1, &bufferObjectID_);
    glGenBuffers(1, &indexObject_);
    isInitialized_ = true;
}

template <class I, class V>
TriangleMeshGeometryIndexed<I, V>::~TriangleMeshGeometryIndexed() {
    clear();

    if(!isInitialized_)
        return;

    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER,0);
    glDeleteBuffers(1, &bufferObjectID_);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glDeleteBuffers(1, &indexObject_);
    glDeleteVertexArrays(1, &vertexArrayID_);
}

template <class I, class V>
size_t TriangleMeshGeometryIndexed<I, V>::getNumTriangles() const {
    return triangles_.size();
}

template <class I, class V>
bool TriangleMeshGeometryIndexed<I, V>::isEmpty() const {
    return (triangles_.size() == 0);
}

template <class I, class V>
void TriangleMeshGeometryIndexed<I, V>::clear() {
    triangles_.clear();
    vertices_.clear();

    invalidate();
}

template <class I, class V>
void TriangleMeshGeometryIndexed<I, V>::calculateNormals() {
    if(!supportsNormals())
        return;
    //TODO:
    LERROR("Not implemented yet!");
}

template <class I, class V>
void TriangleMeshGeometryIndexed<I, V>::addTriangle(const tgt::Vector3<I>& t) {
    if (vertices_.empty()) {
        LERROR("No vertices have been set; cannot add triangle.");
        return;
    }
    tgtAssert(t.x < vertices_.size() && t.y < vertices_.size() && t.z < vertices_.size(), "Invalid triangle indices");
    triangles_.push_back(TriangleType(t.x, t.y, t.z));
    invalidate();
}

template <class I, class V>
void voreen::TriangleMeshGeometryIndexed<I, V>::addTriangle(const TriangleType& t) {
    if (vertices_.empty()) {
        LERROR("No vertices have been set; cannot add triangle.");
        return;
    }
    tgtAssert(t.v_[0] < vertices_.size() && t.v_[1] < vertices_.size() && t.v_[2] < vertices_.size(), "Invalid triangle indices");
    triangles_.push_back(t);
    invalidate();
}

template <class I, class V>
const typename TriangleMeshGeometryIndexed<I, V>::TriangleType& TriangleMeshGeometryIndexed<I, V>::getTriangle(size_t i) const {
    tgtAssert(i < triangles_.size(), "Invalid triangle index");
    return triangles_[i];
}

template <class I, class V>
const typename TriangleMeshGeometryIndexed<I, V>::TriangleType& TriangleMeshGeometryIndexed<I, V>::operator[] (size_t i) const {
    tgtAssert(i < triangles_.size(), "Invalid triangle index");
    return triangles_[i];
}

template <class I, class V>
void TriangleMeshGeometryIndexed<I, V>::setTriangle(const TriangleType& t, size_t i) {
    tgtAssert(i < triangles_.size(), "Invalid triangle index");
    triangles_[i] = t;

    invalidate();
}

template <class I, class V>
tgt::Bounds TriangleMeshGeometryIndexed<I, V>::getBoundingBox(bool transformed) const {
    if(!boundingBoxValid_)
        updateBoundingBox();

    if(transformed)
        return boundingBox_.transform(getTransformationMatrix());
    else
        return boundingBox_;
}

template <class I, class V>
void TriangleMeshGeometryIndexed<I, V>::serialize(Serializer& s) const {
    TriangleMeshGeometryBase::serialize(s);
    s.serializeBinaryBlob("triangles", triangles_);
    s.serializeBinaryBlob("vertices", vertices_);
}

template <class I, class V>
void TriangleMeshGeometryIndexed<I, V>::deserialize(Deserializer& s) {
    TriangleMeshGeometryBase::deserialize(s);
    s.deserializeBinaryBlob("triangles", triangles_);
    s.deserializeBinaryBlob("vertices", vertices_);
}

template <class I, class V>
void TriangleMeshGeometryIndexed<I, V>::invalidate() {
    boundingBoxValid_ = false;
    bufferObjectValid_ = false;
    indexObjectValid_ = false;
}

template <class I, class V>
void TriangleMeshGeometryIndexed<I, V>::render() const {
    updateBuffer();

    if(isEmpty())
        return;

   // MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
   // MatStack.pushMatrix();
   // MatStack.multMatrix(getTransformationMatrix());

    glBindVertexArray(vertexArrayID_);
    glBindBuffer(GL_ARRAY_BUFFER, bufferObjectID_);
    VertexType::setupVertexAttributePointers();

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexObject_);
    glDrawElements(GL_TRIANGLES, static_cast<GLsizei>(triangles_.size() * 3), indexDataType_, NULL);

    LGL_ERROR;

    VertexType::disableVertexAttributePointers();

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

   // MatStack.popMatrix();
}

template <class I, class V>
bool TriangleMeshGeometryIndexed<I, V>::equals(const Geometry* geometry, double epsilon) const {
    const TriangleMeshGeometryIndexed<I, V>* triMesh = dynamic_cast<const TriangleMeshGeometryIndexed<I, V>*>(geometry);

    if(triMesh) {
        if(getNumTriangles() != triMesh->getNumTriangles())
            return false;

        for(size_t i=0; i<triangles_.size(); i++) {
            if (triangles_[i].v_[0] != triMesh->triangles_[i].v_[0] ||
                triangles_[i].v_[1] != triMesh->triangles_[i].v_[1] ||
                triangles_[i].v_[2] != triMesh->triangles_[i].v_[2]  )
            {
                return false;
            }
        }

        for(size_t i=0; i<vertices_.size(); i++) {
            if (!vertices_[i].equals(triMesh->vertices_[i]))
                return false;
        }
        return true;
    }
    else
        return false;
}

template <class I, class V>
void TriangleMeshGeometryIndexed<I, V>::triangulate(const std::vector<I>& poly) {
    if(poly.size() < 3)
        return;

    for(size_t i=2; i<poly.size(); i++)
        addTriangle(TriangleType(poly[0], poly[i-1], poly[i]));
}

template <class I, class V>
void TriangleMeshGeometryIndexed<I, V>::addQuad(const I& v1, const I& v2, const I& v3, const I& v4) {
    addTriangle(TriangleType(v1, v2, v3));
    addTriangle(TriangleType(v1, v3, v4));
}

template <class I, class V>
void TriangleMeshGeometryIndexed<I, V>::addMesh(const TriangleMeshGeometryIndexed<I, V>* mesh) {
    //get old vertices
    size_t startVertex = vertices_.size();
    //add new vertices
    vertices_.insert(vertices_.end(),mesh->vertices_.begin(),mesh->vertices_.end());
    if(getTransformationMatrix() != mesh->getTransformationMatrix()) {
        tgt::mat4 m;
        getTransformationMatrix().invert(m);
        m = m * mesh->getTransformationMatrix();

        tgt::mat3 normalTransform;
        bool transformNormals = getVertexLayout() & NORMAL;
        if(transformNormals) {
            m.getRotationalPartMat3().invert(normalTransform);
            normalTransform = transpose(normalTransform);
        }

        //transform vertices
        for(size_t i = startVertex; i < vertices_.size(); i++) {
            vertices_[i].pos_ = m * vertices_[i].pos_;
            if(transformNormals)
                vertices_[i].setNormal(normalize(normalTransform * vertices_[i].getNormal()));
        }
    }

    //copy the triangles and shift them
    for(size_t i=0; i<mesh->getNumTriangles(); i++) {
        TriangleType t = mesh->getTriangle(i);
        t.v_[0] += (I)startVertex; t.v_[1] += (I)startVertex; t.v_[2] += (I)startVertex;
        addTriangle(t);
    }

    //handle textures
    if(supportsTextureData()) {
        //fix vertices
        for(size_t i = startVertex; i < vertices_.size(); i++) {
            tgt::ivec2 texIndex = vertices_[i].getTexIndex();
            texIndex.x = texIndex.x << textureFilenames_.size();
            vertices_[i].setTexIndex(texIndex);
        }
        //LDEBUG("last vertex : " << vertices_.back().getTexIndex().x);
        textureFilenames_.insert(textureFilenames_.end(),mesh->textureFilenames_.begin(),mesh->textureFilenames_.end());
        //update textures
        loadTextureData();
    }
}

template <class I, class V>
void TriangleMeshGeometryIndexed<I, V>::updateBoundingBox() const {
    tgt::Bounds bb;
    for(size_t i=0; i<triangles_.size(); i++) {
        bb.addPoint(vertices_.at(triangles_[i].v_[0]).pos_);
        bb.addPoint(vertices_.at(triangles_[i].v_[1]).pos_);
        bb.addPoint(vertices_.at(triangles_[i].v_[2]).pos_);
    }
    boundingBox_ = bb;
    boundingBoxValid_ = true;
}

template <class I, class V>
void TriangleMeshGeometryIndexed<I, V>::updateBuffer() const {
    if(bufferObjectValid_ && indexObjectValid_)
        return;

    if(isEmpty())
        return;

    if(!isInitialized_)
        initialize();

    glBindVertexArray(vertexArrayID_);
    glBindBuffer(GL_ARRAY_BUFFER, bufferObjectID_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(VertexType) * vertices_.size(), &vertices_[0], GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexObject_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(TriangleType) * triangles_.size(), &triangles_[0], GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    bufferObjectValid_ = true;
    indexObjectValid_ = true;
}

template<class I, class V>
std::unique_ptr<Geometry> TriangleMeshGeometryIndexed<I, V>::clone() const {
    TriangleMeshGeometryIndexed<I, V>* newGeom = static_cast<TriangleMeshGeometryIndexed<I, V>*>(create());
    newGeom->vertices_ = vertices_;
    newGeom->triangles_ = triangles_;
    newGeom->setTransformationMatrix(getTransformationMatrix());
    newGeom->getMetaDataContainer() = MetaDataContainer(getMetaDataContainer());
    newGeom->setTextureDataAndFilenames(textureData_, textureFilenames_, false);
    return std::unique_ptr<Geometry>(newGeom);
}

template<class I, class V>
void TriangleMeshGeometryIndexed<I, V>::addVertex(size_t index, tgt::vec3 pos, tgt::vec3 normal, tgt::vec4 color, tgt::vec2 texCoord, tgt::ivec2 texIndex) {
    tgtAssert(index < vertices_.size(), "index out of bounds");
    vertices_[index] = V(pos,color,normal,texCoord,texIndex);
}

template<class I, class V>
void TriangleMeshGeometryIndexed<I, V>::addQuadGeometry(float xWidth, float yHeight, tgt::vec3 center, tgt::vec4 color) {
    //check parameter
    tgtAssert(xWidth > 0, "xWidth must be greater zero!");
    tgtAssert(yHeight > 0, "yHeight must be greater zero!");
    if(xWidth <= 0.f) {
        xWidth = 1.f;
        LWARNING("xWidth must be greater zero. Set to 1.f.");
    }
    if(yHeight <= 0.f) {
        yHeight = 1.f;
        LWARNING("yHeight must be greater zero. Set to 1.f.");
    }

    //resize vectors
    size_t startVertex = vertices_.size();
    vertices_.resize(startVertex+4);

    float xHalf = xWidth/2.f;
    float yHalf = yHeight/2.f;

    //create base verticies
    // + z side
    tgt::vec3 normal = tgt::vec3(0.f,0.f,1.f);
    addVertex(startVertex,tgt::vec3(xHalf,yHalf,0),normal,color,tgt::vec2(1.f,1.f), tgt::ivec2(1 << textureFilenames_.size(),1));
    addVertex(startVertex+1,tgt::vec3(-xHalf,yHalf,0),normal,color,tgt::vec2(0.f,1.f), tgt::ivec2(1 << textureFilenames_.size(),1));
    addVertex(startVertex+2,tgt::vec3(-xHalf,-yHalf,0),normal,color,tgt::vec2(0.f,0.f), tgt::ivec2(1 << textureFilenames_.size(),1));
    addVertex(startVertex+3,tgt::vec3(xHalf,-yHalf,0),normal,color,tgt::vec2(1.f,0.f), tgt::ivec2(1 << textureFilenames_.size(), 1));

    //shift by center
    for(size_t i = startVertex; i < vertices_.size(); i++)
        vertices_[i].pos_ += center;

    //create indices
    for(I a = I(startVertex); a <= I(startVertex+3) ; a +=4) {
        addTriangle(TriangleType(a,a+1,a+2)); addTriangle(TriangleType(a+2,a+3,a));
    }
    //update buffers
    invalidate();
}

template<class I, class V>
void TriangleMeshGeometryIndexed<I, V>::addCuboidGeometry(float xWidth, float yHeight, float zDepth, tgt::vec3 center, tgt::vec4 color) {
    //check parameter
    tgtAssert(xWidth > 0, "xWidth must be greater zero!");
    tgtAssert(yHeight > 0, "yHeight must be greater zero!");
    tgtAssert(zDepth > 0, "zDepth must be greater zero!");
    if(xWidth <= 0.f) {
        xWidth = 1.f;
        LWARNING("xWidth must be greater zero. Set to 1.f.");
    }
    if(yHeight <= 0.f) {
        yHeight = 1.f;
        LWARNING("yHeight must be greater zero. Set to 1.f.");
    }
    if(zDepth <= 0.f) {
        zDepth = 1.f;
        LWARNING("zDepth must be greater zero. Set to 1.f.");
    }

    //magic numbers
    size_t ind = 0;         ///< index for vertex construction
    size_t totalVert = 24;  ///< number of all vertices

    //resize vectors
    size_t startVertex = vertices_.size();
    vertices_.resize(startVertex+totalVert);

    float xHalf = xWidth/2.f;
    float yHalf = yHeight/2.f;
    float zHalf = zDepth/2.f;

    //create base verticies
    // + x side
    tgt::vec3 normal = tgt::vec3(1.f,0.f,0.f);
    addVertex(startVertex+(ind++),tgt::vec3(xHalf,yHalf,zHalf),normal,color);
    addVertex(startVertex+(ind++),tgt::vec3(xHalf,-yHalf,zHalf),normal,color);
    addVertex(startVertex+(ind++),tgt::vec3(xHalf,-yHalf,-zHalf),normal,color);
    addVertex(startVertex+(ind++),tgt::vec3(xHalf,yHalf,-zHalf),normal,color);
    // - x side
    normal = tgt::vec3(-1.f,0.f,0.f);
    addVertex(startVertex+(ind++),tgt::vec3(-xHalf,yHalf,zHalf),normal,color);
    addVertex(startVertex+(ind++),tgt::vec3(-xHalf,yHalf,-zHalf),normal,color);
    addVertex(startVertex+(ind++),tgt::vec3(-xHalf,-yHalf,-zHalf),normal,color);
    addVertex(startVertex+(ind++),tgt::vec3(-xHalf,-yHalf,zHalf),normal,color);
    // + y side
    normal = tgt::vec3(0.f,1.f,0.f);
    addVertex(startVertex+(ind++),tgt::vec3(xHalf,yHalf,zHalf),normal,color);
    addVertex(startVertex+(ind++),tgt::vec3(xHalf,yHalf,-zHalf),normal,color);
    addVertex(startVertex+(ind++),tgt::vec3(-xHalf,yHalf,-zHalf),normal,color);
    addVertex(startVertex+(ind++),tgt::vec3(-xHalf,yHalf,zHalf),normal,color);
    // - y side
    normal = tgt::vec3(0.f,-1.f,0.f);
    addVertex(startVertex+(ind++),tgt::vec3(xHalf,-yHalf,zHalf),normal,color);
    addVertex(startVertex+(ind++),tgt::vec3(-xHalf,-yHalf,zHalf),normal,color);
    addVertex(startVertex+(ind++),tgt::vec3(-xHalf,-yHalf,-zHalf),normal,color);
    addVertex(startVertex+(ind++),tgt::vec3(xHalf,-yHalf,-zHalf),normal,color);
    // + z side
    normal = tgt::vec3(0.f,0.f,1.f);
    addVertex(startVertex+(ind++),tgt::vec3(xHalf,yHalf,zHalf),normal,color);
    addVertex(startVertex+(ind++),tgt::vec3(-xHalf,yHalf,zHalf),normal,color);
    addVertex(startVertex+(ind++),tgt::vec3(-xHalf,-yHalf,zHalf),normal,color);
    addVertex(startVertex+(ind++),tgt::vec3(xHalf,-yHalf,zHalf),normal,color);
    // - z side
    normal = tgt::vec3(0.f,0.f,-1.f);
    addVertex(startVertex+(ind++),tgt::vec3(xHalf,yHalf,-zHalf),normal,color);
    addVertex(startVertex+(ind++),tgt::vec3(xHalf,-yHalf,-zHalf),normal,color);
    addVertex(startVertex+(ind++),tgt::vec3(-xHalf,-yHalf,-zHalf),normal,color);
    addVertex(startVertex+(ind++),tgt::vec3(-xHalf,yHalf,-zHalf),normal,color);

    //shift by center
    for(size_t i = startVertex; i < vertices_.size(); i++)
        vertices_[i].pos_ += center;

    //create indices
    for(I a = startVertex; a < startVertex+totalVert ; a +=4) {
        addTriangle(TriangleType(a,a+1,a+2)); addTriangle(TriangleType(a+2,a+3,a));
    }
    //update buffers
    invalidate();
}

template<class I, class V>
void TriangleMeshGeometryIndexed<I, V>::addDiskGeometry(float innerRadius, float outerRadius, tgt::vec3 center, tgt::vec4 color, SimpleGeometryQuality quality) {
    //check parameters
    tgtAssert(innerRadius >= 0, "innerRadius must be greater or equal zero!");
    tgtAssert(outerRadius > 0, "outerRadius must be greater or equal zero!");
    tgtAssert(innerRadius < outerRadius, "innerRadius must be smaller than outerRadius");
    if(innerRadius < 0.f) {
        innerRadius = 0.f;
        LWARNING("innerRadius must be greater zero. Set to 0.f.");
    }
    if(outerRadius < 0.f) {
        outerRadius = 1.f;
        LWARNING("outerRadius must be greater zero. Set to 1.f.");
    }
    if(innerRadius >= outerRadius) {
        innerRadius = 0.f;
        outerRadius = 1.f;
        LWARNING("innerRadius must be smaller than outerRadius. Set innerRadius = 0.f, outerRadius = 1.f.");
    }
    //resize vectors
    size_t startVertex = vertices_.size();
    size_t numVertices = static_cast<size_t>(8*std::pow(2.f,static_cast<float>(quality)));
    vertices_.resize(startVertex+numVertices);
    //set default normal
    tgt::vec3 defaultNormal(0.f,0.f,-1.f);
    //help variables
    size_t quadStep = numVertices/8;
    size_t beginOfInner = numVertices/2;

    //create base verticies
    addVertex(startVertex+0,          tgt::vec3(outerRadius,0.f,0.f),defaultNormal,color);
    addVertex(startVertex+quadStep,   tgt::vec3(0.f,-outerRadius,0.f),defaultNormal,color);
    addVertex(startVertex+2*quadStep, tgt::vec3(-outerRadius,0.f,0.f),defaultNormal,color);
    addVertex(startVertex+3*quadStep, tgt::vec3(0.f,outerRadius,0.f),defaultNormal,color);

    addVertex(startVertex+4*quadStep, tgt::vec3(innerRadius,0.f,0.f),defaultNormal,color);
    addVertex(startVertex+5*quadStep, tgt::vec3(0.f,-innerRadius,0.f),defaultNormal,color);
    addVertex(startVertex+6*quadStep, tgt::vec3(-innerRadius,0.f,0.f),defaultNormal,color);
    addVertex(startVertex+7*quadStep, tgt::vec3(0.f,innerRadius,0.f),defaultNormal,color);

    //add quality vertices
    for(int q = 1; q <= quality; q++) {
        for(int i = 0; i < 4*std::pow(2.f,q-1.f); i++) {
            size_t currentIndex = startVertex+static_cast<size_t>((i+0.5f)*quadStep);
            //handle outer
            addVertex(currentIndex, tgt::normalize(vertices_[startVertex+(i*quadStep)].pos_ + vertices_[startVertex+((i+1)*quadStep == beginOfInner ? 0 : (i+1)*quadStep)].pos_)*outerRadius,defaultNormal, color);
            //handle inner (special handling of innerRadius_ = 0)
            addVertex(beginOfInner+currentIndex,(innerRadius == 0.f ? tgt::vec3::zero : tgt::normalize(vertices_[startVertex+beginOfInner+(i*quadStep)].pos_ + vertices_[startVertex+beginOfInner+((i+1)*quadStep == beginOfInner ? 0 : (i+1)*quadStep)].pos_)*innerRadius), defaultNormal, color);
        }
        //in next iteration, use half step size
        quadStep /= 2;
    }

    //shift by center
    for(size_t i = startVertex; i < vertices_.size(); i++)
        vertices_[i].pos_ += center;

    //create indices
    for(I i = 0; i < beginOfInner; i++) {
        addTriangle(TriangleType(I(startVertex + i), I(startVertex + i + beginOfInner), I(startVertex + (static_cast<size_t>(i + 1) == beginOfInner ? 0 : i + 1))));
        addTriangle(TriangleType(I(startVertex + i + beginOfInner), I(startVertex + (static_cast<size_t>(i + 1) == beginOfInner ? beginOfInner : i + beginOfInner + 1)),
                                 I(startVertex + (static_cast<size_t>(i + 1) == beginOfInner ? 0 : i + 1))));
    }
    //invalidate geometry
    invalidate();
}

template<class I, class V>
void TriangleMeshGeometryIndexed<I, V>::addCylinderGeometry(float bottomRadius, float topRadius, float length, tgt::vec3 center, tgt::vec4 color, SimpleGeometryQuality quality) {
    tgtAssert(bottomRadius >= 0.f, "bottomRadius must be greater or equal zero!");
    tgtAssert(topRadius >= 0.f, "topRadius must be greater or equal zero!");
    tgtAssert(bottomRadius != 0.f || topRadius != 0.f, "Add least one radii must be greater zero!");
    tgtAssert(length > 0.f, "length must be greater zero!");
    if(bottomRadius < 0.f) {
        bottomRadius = 1.f;
        LWARNING("bottomRadius must be greater zero. Set to 1.f.");
    }
    if(topRadius < 0.f) {
        topRadius = 1.f;
        LWARNING("topRadius must be greater zero. Set to 1.f.");
    }
    if(bottomRadius == 0.f && topRadius == 0.f) {
        bottomRadius = 1.f;
        topRadius = 1.f;
        LWARNING("bottomRadius or topRadius must be unequal zero. Set bottomRadius = 1.f, topRadius = 1.f.");
    }
    if(length <= 0.f) {
        length = 1.f;
        LWARNING("length must be greater zero. Set to 1.f.");
    }
        //resize vectors
    size_t startVertex = vertices_.size();
    size_t numVertices = static_cast<size_t>(8*std::pow(2.f,static_cast<float>(quality)));
    vertices_.resize(startVertex+numVertices);
    //help variables
    size_t quadStep = numVertices/8;
    size_t beginOfTop = numVertices/2;

    //create normals
    tgt::vec3 yPosNormal = tgt::vec3(0,topRadius,length)-tgt::vec3(0,bottomRadius,0);
    tgt::vec3 yNegNormal; tgt::vec3 xPosNormal; tgt::vec3 xNegNormal;
    float tmp = yPosNormal.y;
    yPosNormal.y = yPosNormal.z;
    yPosNormal.z = -tmp;
    yPosNormal = tgt::normalize(yPosNormal);
    yNegNormal = tgt::vec3(0,-yPosNormal.y,yPosNormal.z);
    xPosNormal = tgt::vec3(yPosNormal.y,0,yPosNormal.z);
    xNegNormal = tgt::vec3(-yPosNormal.y,0,yPosNormal.z);

    //create base verticies
    addVertex(startVertex+0,          tgt::vec3(bottomRadius,0.f,0.f),xPosNormal, color);
    addVertex(startVertex+1*quadStep, tgt::vec3(0.f,-bottomRadius,0.f),yNegNormal, color);
    addVertex(startVertex+2*quadStep, tgt::vec3(-bottomRadius,0.f,0.f),xNegNormal, color);
    addVertex(startVertex+3*quadStep, tgt::vec3(0.f,bottomRadius,0.f),yPosNormal, color);

    addVertex(startVertex+4*quadStep, tgt::vec3(topRadius,0.f,length),xPosNormal, color);
    addVertex(startVertex+5*quadStep, tgt::vec3(0.f,-topRadius,length),yNegNormal, color);
    addVertex(startVertex+6*quadStep, tgt::vec3(-topRadius,0.f,length),xNegNormal, color);
    addVertex(startVertex+7*quadStep, tgt::vec3(0.f,topRadius,length),yPosNormal, color);

    //add quality vertices
    for(int q = 1; q <= quality; q++) {
        for(int i = 0; i < 4*std::pow(2.f,q-1.f); i++) {
            size_t currentIndex = startVertex+static_cast<size_t>((i+0.5f)*quadStep);
            //calculate normal
            tgt::vec3 normal = ((getVertexLayout() & NORMAL) == 0 ? tgt::vec3::zero : tgt::normalize(vertices_[startVertex+(i*quadStep)].getNormal() + vertices_[startVertex+((i+1)*quadStep == beginOfTop ? 0 : (i+1)*quadStep)].getNormal()));
            //handle bottom (special handling of bottomRadius_ = 0)
            addVertex(currentIndex,(bottomRadius == 0.f ? tgt::vec3::zero : tgt::normalize(vertices_[startVertex+i*quadStep].pos_ + vertices_[startVertex+((i+1)*quadStep == beginOfTop ? 0 : (i+1)*quadStep)].pos_)*bottomRadius),normal, color);
            //handle top (special handling of topRadius_ = 0)
            addVertex(beginOfTop+currentIndex, (topRadius == 0.f ? tgt::vec3(0.f,0.f,length) : tgt::vec3( tgt::normalize(vertices_[startVertex+beginOfTop+(i*quadStep)].pos_.xy() + vertices_[startVertex+beginOfTop+((i+1)*quadStep == beginOfTop ? 0 : (i+1)*quadStep)].pos_.xy())*topRadius,length)),normal, color);
        }
        //in next iteration, use half step size
        quadStep /= 2;
    }

    //shift by center
    for(size_t i = startVertex; i < vertices_.size(); i++)
        vertices_[i].pos_ += center;

    //create indices
    for(I i = 0; i < beginOfTop; i++) {
        addTriangle(TriangleType(I(startVertex + i),
                                 I(startVertex + i + beginOfTop),
                                 I(startVertex + (static_cast<size_t>(i + 1) == beginOfTop ? 0 : i + 1))));
        addTriangle(TriangleType(I(startVertex + i + beginOfTop),
                                 I(startVertex + (static_cast<size_t>(i + 1) == beginOfTop ? beginOfTop : i + beginOfTop + 1)),
                                 I(startVertex + (static_cast<size_t>(i + 1) == beginOfTop ? 0 : i + 1))));
    }
    //invalidate geometry
    invalidate();
}

template<class I, class V>
void TriangleMeshGeometryIndexed<I, V>::addSphereGeometry(float radius, tgt::vec3 center, tgt::vec4 color, SimpleGeometryQuality quality) {

    //addSphereGeometry(radius, center, color, static_cast<size_t>(quality + 3));

    // TODO: workaround which does not actually generate spheres...
    //resize vectors
    size_t startVertex = vertices_.size();
    vertices_.resize(startVertex+6);

    //create vertices
    addVertex(startVertex+0, tgt::vec3(radius,0.f,0.f),tgt::vec3(1.f,0.f,0.f),color);
    addVertex(startVertex+1, tgt::vec3(-radius,0.f,0.f),tgt::vec3(-1.f,0.f,0.f),color);
    addVertex(startVertex+2, tgt::vec3(0.f,radius,0.f),tgt::vec3(0.f,1.f,0.f),color);
    addVertex(startVertex+3, tgt::vec3(0.f,-radius,0.f),tgt::vec3(0.f,-1.f,0.f),color);
    addVertex(startVertex+4, tgt::vec3(0.f,0.f,radius),tgt::vec3(0.f,0.f,1.f),color);
    addVertex(startVertex+5, tgt::vec3(0.f,0.f,-radius),tgt::vec3(0.f,0.f,-1.f),color);

    //shift by center
    for(size_t i = startVertex; i < vertices_.size(); i++)
        vertices_[i].pos_ += center;

    //create indices
    addTriangle(TriangleType(I(startVertex + 2), I(startVertex + 0), I(startVertex + 4)));
    addTriangle(TriangleType(I(startVertex + 2), I(startVertex + 0), I(startVertex + 5)));
    addTriangle(TriangleType(I(startVertex + 3), I(startVertex + 0), I(startVertex + 4)));
    addTriangle(TriangleType(I(startVertex + 3), I(startVertex + 0), I(startVertex + 5)));
    addTriangle(TriangleType(I(startVertex + 2), I(startVertex + 1), I(startVertex + 4)));
    addTriangle(TriangleType(I(startVertex + 2), I(startVertex + 1), I(startVertex + 5)));
    addTriangle(TriangleType(I(startVertex + 3), I(startVertex + 1), I(startVertex + 4)));
    addTriangle(TriangleType(I(startVertex + 3), I(startVertex + 1), I(startVertex + 5)));

    invalidate();
}

/*
 * broken, therefore commented, FIXME: this uses Triangle Strips...
 *
template<class I, class V>
void TriangleMeshGeometryIndexed<I,V>::addSphereGeometry(float radius, tgt::vec3 center, tgt::vec4 color, size_t numSlicesAndTiles) {
    //check parameter
    tgtAssert(radius > 0, "Radius must be greater zero!");
    if(radius < 0.f) {
        radius = 1.f;
        LWARNING("Radius must be greater zero. Set to 1.f.");
    }

    tgtAssert(numSlicesAndTiles > 2, "Number of slices and tiles must be > 2");
    if (numSlicesAndTiles < 3) {
        numSlicesAndTiles = 3;
        LWARNING("Number of slices and tiles must be > 2");
    }

    // create the sphere vertices
    const float dfi = tgt::PIf / static_cast<float>(numSlicesAndTiles);
    const float dth = (2.0f * tgt::PIf) / static_cast<float>(numSlicesAndTiles);
    const float du  = 1.0f / static_cast<float>(numSlicesAndTiles);
    const float dv  = 1.0f / static_cast<float>(numSlicesAndTiles);
    std::vector<float>  cs_fi, cs_th, sn_fi, sn_th;

    float fi = 0.f;// PI / 2.0;
    float th = 0.f;
    float u  = 0.f;
    float v  = 0.f;

    for (size_t i = 0; i <= numSlicesAndTiles; ++i) {
        cs_fi.push_back(std::cos(fi));
        cs_th.push_back(std::cos(th));

        sn_fi.push_back(std::sin(fi));
        sn_th.push_back(std::sin(th));

        fi += dfi;
        th += dth;
    }

    //resize vectors
    size_t startVertex = vertices_.size();
    vertices_.resize(startVertex + (numSlicesAndTiles + 1) * (numSlicesAndTiles + 1));

    size_t currentVertexIndex = startVertex;
    for (size_t i = 0; i <= numSlicesAndTiles; ++i) {
        for (size_t j = 0; j <= numSlicesAndTiles; ++j) {
            //Vertex vx;
            size_t k = j % numSlicesAndTiles;
            tgt::vec3 normal = tgt::vec3(sn_fi[i] * cs_th[k], sn_fi[i] * sn_th[k], cs_fi[i]);
            tgt::vec3 position = radius * normal;
            tgt::vec2 uv = tgt::vec2(u, v);

            addVertex(currentVertexIndex, position, normal, color, uv);
            currentVertexIndex++;

            u += du;
        }

        u = 0;
        v += dv;
    }

    //shift vertices by center
    for(size_t i = startVertex; i < vertices_.size(); i++)
        vertices_[i].pos_ += center;

    // create the indices for the triangles
    I offset  = static_cast<I>(startVertex);
    for (size_t j = 0; j < numSlicesAndTiles; ++j) {
        for (size_t i = 1; i <= numSlicesAndTiles; ++i) {
            I idx = offset + static_cast<I>(i);

            tgt::Vector3<I> t1(idx - 1, idx + static_cast<I>(numSlicesAndTiles), idx);
            addTriangle(t1);

            tgt::Vector3<I> t2(idx, idx + static_cast<I>(numSlicesAndTiles), idx + static_cast<I>(numSlicesAndTiles) + 1);
            addTriangle(t2);
        }
        offset += numSlicesAndTiles + 1;
    }

    invalidate();
}*/


template<class I, class V>
void TriangleMeshGeometryIndexed<I, V>::addLetterX(float height, tgt::vec4 color) {
    //check parameter
    tgtAssert(height > 0, "Height must be greater zero!");
    if(height < 0.f) {
        height = 1.f;
        LWARNING("Height must be greater zero. Set to 1.0f.");
    }

    //resize vectors
    size_t startVertex = vertices_.size();
    vertices_.resize(startVertex+48);

    //create verticies
        //front
    addVertex(startVertex+ 0, tgt::vec3(-0.46f, 0.5f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 1, tgt::vec3(-0.22f, 0.5f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 2, tgt::vec3( 0.22f, 0.5f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 3, tgt::vec3( 0.46f, 0.5f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 4, tgt::vec3(-0.46f,-0.5f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 5, tgt::vec3(-0.22f,-0.5f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 6, tgt::vec3( 0.22f,-0.5f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 7, tgt::vec3( 0.46f,-0.5f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
        //back
    addVertex(startVertex+ 8, tgt::vec3(-0.46f, 0.5f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+ 9, tgt::vec3(-0.22f, 0.5f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+10, tgt::vec3( 0.22f, 0.5f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+11, tgt::vec3( 0.46f, 0.5f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+12, tgt::vec3(-0.46f,-0.5f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+13, tgt::vec3(-0.22f,-0.5f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+14, tgt::vec3( 0.22f,-0.5f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+15, tgt::vec3( 0.46f,-0.5f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
        //top
    addVertex(startVertex+16,  tgt::vec3(-0.46f, 0.5f,-0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
    addVertex(startVertex+17,  tgt::vec3(-0.22f, 0.5f,-0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
    addVertex(startVertex+18,  tgt::vec3(-0.22f, 0.5f, 0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
    addVertex(startVertex+19,  tgt::vec3(-0.46f, 0.5f, 0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
     addVertex(startVertex+20, tgt::vec3( 0.22f, 0.5f,-0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
    addVertex(startVertex+21,  tgt::vec3( 0.46f, 0.5f,-0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
    addVertex(startVertex+22,  tgt::vec3( 0.46f, 0.5f, 0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
    addVertex(startVertex+23,  tgt::vec3( 0.22f, 0.5f, 0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
        //bottom
    addVertex(startVertex+24,  tgt::vec3(-0.46f,-0.5f,-0.125f), tgt::vec3( 0.f,-1.f, 0.f),color);
    addVertex(startVertex+25,  tgt::vec3(-0.22f,-0.5f,-0.125f), tgt::vec3( 0.f,-1.f, 0.f),color);
    addVertex(startVertex+26,  tgt::vec3(-0.22f,-0.5f, 0.125f), tgt::vec3( 0.f,-1.f, 0.f),color);
    addVertex(startVertex+27,  tgt::vec3(-0.46f,-0.5f, 0.125f), tgt::vec3( 0.f,-1.f, 0.f),color);
     addVertex(startVertex+28, tgt::vec3( 0.22f,-0.5f,-0.125f), tgt::vec3( 0.f,-1.f, 0.f),color);
    addVertex(startVertex+29,  tgt::vec3( 0.46f,-0.5f,-0.125f), tgt::vec3( 0.f,-1.f, 0.f),color);
    addVertex(startVertex+30,  tgt::vec3( 0.46f,-0.5f, 0.125f), tgt::vec3( 0.f,-1.f, 0.f),color);
    addVertex(startVertex+31,  tgt::vec3( 0.22f,-0.5f, 0.125f), tgt::vec3( 0.f,-1.f, 0.f),color);
        //left
    tgt::vec3 norm = tgt::normalize(tgt::vec3(-1.f,-0.68f,0.f));
    addVertex(startVertex+32, tgt::vec3(-0.46f, 0.5f,-0.125f), norm,color);
    addVertex(startVertex+33, tgt::vec3( 0.22f,-0.5f,-0.125f), norm,color);
    addVertex(startVertex+34, tgt::vec3( 0.22f,-0.5f, 0.125f), norm,color);
    addVertex(startVertex+35, tgt::vec3(-0.46f, 0.5f, 0.125f), norm,color);
              norm = tgt::normalize(tgt::vec3(-1.f, 0.68f,0.f));
    addVertex(startVertex+36, tgt::vec3(-0.46f,-0.5f,-0.125f), norm,color);
    addVertex(startVertex+37, tgt::vec3( 0.22f, 0.5f,-0.125f), norm,color);
    addVertex(startVertex+38, tgt::vec3( 0.22f, 0.5f, 0.125f), norm,color);
    addVertex(startVertex+39, tgt::vec3(-0.46f,-0.5f, 0.125f), norm,color);
        //right
              norm = tgt::normalize(tgt::vec3( 1.f, 0.68f,0.f));
    addVertex(startVertex+40, tgt::vec3(-0.22f, 0.5f,-0.125f), norm,color);
    addVertex(startVertex+41, tgt::vec3( 0.46f,-0.5f,-0.125f), norm,color);
    addVertex(startVertex+42, tgt::vec3( 0.46f,-0.5f, 0.125f), norm,color);
    addVertex(startVertex+43, tgt::vec3(-0.22f, 0.5f, 0.125f), norm,color);
             norm = tgt::normalize(tgt::vec3( 1.f,-0.68f,0.f));
    addVertex(startVertex+44, tgt::vec3(-0.22f,-0.5f,-0.125f), norm,color);
    addVertex(startVertex+45, tgt::vec3( 0.46f, 0.5f,-0.125f), norm,color);
    addVertex(startVertex+46, tgt::vec3( 0.46f, 0.5f, 0.125f), norm,color);
    addVertex(startVertex+47, tgt::vec3(-0.22f,-0.5f, 0.125f), norm,color);

    //resize by height
    for(size_t i = startVertex; i < vertices_.size(); i++)
        vertices_[i].pos_ *= height;

    //create indices
        //front
    addTriangle(TriangleType(I(startVertex + 0), I(startVertex + 1), I(startVertex + 6)));
    addTriangle(TriangleType(I(startVertex + 6), I(startVertex + 1), I(startVertex + 7)));
    addTriangle(TriangleType(I(startVertex + 2), I(startVertex + 3), I(startVertex + 5)));
    addTriangle(TriangleType(I(startVertex + 5), I(startVertex + 4), I(startVertex + 2)));
        //back
    addTriangle(TriangleType(I(startVertex + 11), I(startVertex + 10), I(startVertex + 13)));
    addTriangle(TriangleType(I(startVertex + 10), I(startVertex + 12), I(startVertex + 13)));
    addTriangle(TriangleType(I(startVertex +  9), I(startVertex +  8), I(startVertex + 14)));
    addTriangle(TriangleType(I(startVertex + 14), I(startVertex + 15), I(startVertex +  9)));
        //top
    addTriangle(TriangleType(I(startVertex + 16), I(startVertex + 17), I(startVertex + 18)));
    addTriangle(TriangleType(I(startVertex + 18), I(startVertex + 19), I(startVertex + 16)));
    addTriangle(TriangleType(I(startVertex + 20), I(startVertex + 21), I(startVertex + 22)));
    addTriangle(TriangleType(I(startVertex + 22), I(startVertex + 23), I(startVertex + 20)));
        //bottom
    addTriangle(TriangleType(I(startVertex + 24), I(startVertex + 25), I(startVertex + 26)));
    addTriangle(TriangleType(I(startVertex + 26), I(startVertex + 27), I(startVertex + 24)));
    addTriangle(TriangleType(I(startVertex + 28), I(startVertex + 29), I(startVertex + 30)));
    addTriangle(TriangleType(I(startVertex + 30), I(startVertex + 31), I(startVertex + 28)));
        //left
    addTriangle(TriangleType(I(startVertex + 32), I(startVertex + 33), I(startVertex + 34)));
    addTriangle(TriangleType(I(startVertex + 34), I(startVertex + 35), I(startVertex + 32)));
    addTriangle(TriangleType(I(startVertex + 36), I(startVertex + 37), I(startVertex + 38)));
    addTriangle(TriangleType(I(startVertex + 38), I(startVertex + 39), I(startVertex + 36)));
        //right
    addTriangle(TriangleType(I(startVertex + 40), I(startVertex + 41), I(startVertex + 42)));
    addTriangle(TriangleType(I(startVertex + 42), I(startVertex + 43), I(startVertex + 40)));
    addTriangle(TriangleType(I(startVertex + 44), I(startVertex + 45), I(startVertex + 46)));
    addTriangle(TriangleType(I(startVertex + 46), I(startVertex + 47), I(startVertex + 44)));

    invalidate();
}

template<class I, class V>
void TriangleMeshGeometryIndexed<I, V>::addLetterY(float height, tgt::vec4 color) {
    //check parameter
    tgtAssert(height > 0, "Height must be greater zero!");
    if(height < 0.f) {
        height = 1.f;
        LWARNING("Height must be greater zero. Set to 1.0f.");
    }

    //resize vectors
    size_t startVertex = vertices_.size();
    vertices_.resize(startVertex+52);

    //create verticies
        //front
    addVertex(startVertex+ 0, tgt::vec3(-0.46f, 0.5f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 1, tgt::vec3(-0.22f, 0.5f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 2, tgt::vec3( 0.22f, 0.5f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 3, tgt::vec3( 0.46f, 0.5f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 4, tgt::vec3(-0.10f,-0.1f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 5, tgt::vec3( 0.10f,-0.1f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 6, tgt::vec3(-0.10f,-0.5f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 7, tgt::vec3( 0.10f,-0.5f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
        //back
    addVertex(startVertex+ 8, tgt::vec3(-0.46f, 0.5f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+ 9, tgt::vec3(-0.22f, 0.5f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+10, tgt::vec3( 0.22f, 0.5f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+11, tgt::vec3( 0.46f, 0.5f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+12, tgt::vec3(-0.10f,-0.1f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+13, tgt::vec3( 0.10f,-0.1f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+14, tgt::vec3(-0.10f,-0.5f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+15, tgt::vec3( 0.10f,-0.5f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
        //top
    addVertex(startVertex+16, tgt::vec3(-0.46f, 0.5f,-0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
    addVertex(startVertex+17, tgt::vec3(-0.22f, 0.5f,-0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
    addVertex(startVertex+18, tgt::vec3(-0.22f, 0.5f, 0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
    addVertex(startVertex+19, tgt::vec3(-0.46f, 0.5f, 0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
    addVertex(startVertex+20, tgt::vec3( 0.22f, 0.5f,-0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
    addVertex(startVertex+21, tgt::vec3( 0.46f, 0.5f,-0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
    addVertex(startVertex+22, tgt::vec3( 0.46f, 0.5f, 0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
    addVertex(startVertex+23, tgt::vec3( 0.22f, 0.5f, 0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
        //bottom
    addVertex(startVertex+24, tgt::vec3(-0.10f,-0.5f,-0.125f), tgt::vec3( 0.f,-1.f, 0.f),color);
    addVertex(startVertex+25, tgt::vec3( 0.10f,-0.5f,-0.125f), tgt::vec3( 0.f,-1.f, 0.f),color);
    addVertex(startVertex+26, tgt::vec3( 0.10f,-0.5f, 0.125f), tgt::vec3( 0.f,-1.f, 0.f),color);
    addVertex(startVertex+27, tgt::vec3(-0.10f,-0.5f, 0.125f), tgt::vec3( 0.f,-1.f, 0.f),color);
        //right arm (normal approx)
    tgt::vec3 norm = tgt::normalize(tgt::vec3(-0.6f, 0.34f,0.f));
    addVertex(startVertex+28, tgt::vec3(-0.10f,-0.1f,-0.125f), norm,color);
    addVertex(startVertex+29, tgt::vec3( 0.22f, 0.5f,-0.125f), norm,color);
    addVertex(startVertex+30, tgt::vec3( 0.22f, 0.5f, 0.125f), norm,color);
    addVertex(startVertex+31, tgt::vec3(-0.10f,-0.1f, 0.125f), norm,color);
              norm = tgt::normalize(tgt::vec3(0.6f, -0.34f,0.f));
    addVertex(startVertex+32, tgt::vec3( 0.10f,-0.1f,-0.125f), norm,color);
    addVertex(startVertex+33, tgt::vec3( 0.46f, 0.5f,-0.125f), norm,color);
    addVertex(startVertex+34, tgt::vec3( 0.46f, 0.5f, 0.125f), norm,color);
    addVertex(startVertex+35, tgt::vec3( 0.10f,-0.1f, 0.125f), norm,color);
        //left arm (approx normal)
              norm = tgt::normalize(tgt::vec3(-0.6f,-0.34f,0.f));
    addVertex(startVertex+36, tgt::vec3(-0.46f, 0.5f,-0.125f), norm,color);
    addVertex(startVertex+37, tgt::vec3(-0.10f,-0.1f,-0.125f), norm,color);
    addVertex(startVertex+38, tgt::vec3(-0.10f,-0.1f, 0.125f), norm,color);
    addVertex(startVertex+39, tgt::vec3(-0.46f, 0.5f, 0.125f), norm,color);
              norm = tgt::normalize(tgt::vec3(0.6f, 0.34f,0.f));
    addVertex(startVertex+40, tgt::vec3(-0.22f, 0.5f,-0.125f), norm,color);
    addVertex(startVertex+41, tgt::vec3( 0.10f,-0.1f,-0.125f), norm,color);
    addVertex(startVertex+42, tgt::vec3( 0.10f,-0.1f, 0.125f), norm,color);
    addVertex(startVertex+43, tgt::vec3(-0.22f, 0.5f, 0.125f), norm,color);
        //left and right
    addVertex(startVertex+44, tgt::vec3(-0.10f,-0.5f,-0.125f), tgt::vec3(-1.f, 0.f, 0.f),color);
    addVertex(startVertex+45, tgt::vec3(-0.10f,-0.1f,-0.125f), tgt::vec3(-1.f, 0.f, 0.f),color);
    addVertex(startVertex+46, tgt::vec3(-0.10f,-0.1f, 0.125f), tgt::vec3(-1.f, 0.f, 0.f),color);
    addVertex(startVertex+47, tgt::vec3(-0.10f,-0.5f, 0.125f), tgt::vec3(-1.f, 0.f, 0.f),color);
    addVertex(startVertex+48, tgt::vec3( 0.10f,-0.5f,-0.125f), tgt::vec3( 1.f, 0.f, 0.f),color);
    addVertex(startVertex+49, tgt::vec3( 0.10f,-0.1f,-0.125f), tgt::vec3( 1.f, 0.f, 0.f),color);
    addVertex(startVertex+50, tgt::vec3( 0.10f,-0.1f, 0.125f), tgt::vec3( 1.f, 0.f, 0.f),color);
    addVertex(startVertex+51, tgt::vec3( 0.10f,-0.5f, 0.125f), tgt::vec3( 1.f, 0.f, 0.f),color);


    //resize by height
    for(size_t i = startVertex; i < vertices_.size(); i++)
        vertices_[i].pos_ *= height;

    //create indices
    addTriangle(TriangleType(I(startVertex +  0), I(startVertex +  1), I(startVertex +  4)));
    addTriangle(TriangleType(I(startVertex +  1), I(startVertex +  5), I(startVertex +  4)));
    addTriangle(TriangleType(I(startVertex +  2), I(startVertex +  3), I(startVertex +  5)));
    addTriangle(TriangleType(I(startVertex +  5), I(startVertex +  4), I(startVertex +  2)));
    addTriangle(TriangleType(I(startVertex +  4), I(startVertex +  5), I(startVertex +  6)));
    addTriangle(TriangleType(I(startVertex +  5), I(startVertex +  7), I(startVertex +  6)));
    addTriangle(TriangleType(I(startVertex + 11), I(startVertex + 10), I(startVertex + 13)));
    addTriangle(TriangleType(I(startVertex + 10), I(startVertex + 12), I(startVertex + 13)));
    addTriangle(TriangleType(I(startVertex +  9), I(startVertex +  8), I(startVertex + 12)));
    addTriangle(TriangleType(I(startVertex + 12), I(startVertex + 13), I(startVertex +  9)));
    addTriangle(TriangleType(I(startVertex + 13), I(startVertex + 12), I(startVertex + 14)));
    addTriangle(TriangleType(I(startVertex + 14), I(startVertex + 15), I(startVertex + 13)));
        //top
    addTriangle(TriangleType(I(startVertex + 16), I(startVertex + 17), I(startVertex + 18)));
    addTriangle(TriangleType(I(startVertex + 18), I(startVertex + 19), I(startVertex + 16)));
    addTriangle(TriangleType(I(startVertex + 20), I(startVertex + 21), I(startVertex + 22)));
    addTriangle(TriangleType(I(startVertex + 22), I(startVertex + 23), I(startVertex + 20)));
        //bottom
    addTriangle(TriangleType(I(startVertex + 24), I(startVertex + 25), I(startVertex + 26)));
    addTriangle(TriangleType(I(startVertex + 26), I(startVertex + 27), I(startVertex + 24)));
        //right
    addTriangle(TriangleType(I(startVertex + 28), I(startVertex + 29), I(startVertex + 30)));
    addTriangle(TriangleType(I(startVertex + 30), I(startVertex + 31), I(startVertex + 28)));
    addTriangle(TriangleType(I(startVertex + 32), I(startVertex + 33), I(startVertex + 34)));
    addTriangle(TriangleType(I(startVertex + 34), I(startVertex + 35), I(startVertex + 32)));
        //left
    addTriangle(TriangleType(I(startVertex + 36), I(startVertex + 37), I(startVertex + 38)));
    addTriangle(TriangleType(I(startVertex + 38), I(startVertex + 39), I(startVertex + 36)));
    addTriangle(TriangleType(I(startVertex + 40), I(startVertex + 41), I(startVertex + 42)));
    addTriangle(TriangleType(I(startVertex + 42), I(startVertex + 43), I(startVertex + 40)));
        //left and right
    addTriangle(TriangleType(I(startVertex + 44), I(startVertex + 45), I(startVertex + 46)));
    addTriangle(TriangleType(I(startVertex + 46), I(startVertex + 47), I(startVertex + 44)));
    addTriangle(TriangleType(I(startVertex + 48), I(startVertex + 49), I(startVertex + 50)));
    addTriangle(TriangleType(I(startVertex + 50), I(startVertex + 51), I(startVertex + 48)));

    invalidate();
}

template<class I, class V>
void TriangleMeshGeometryIndexed<I, V>::addLetterZ(float height, tgt::vec4 color) {
    //check parameter
    tgtAssert(height > 0, "Height must be greater zero!");
    if(height < 0.f) {
        height = 1.f;
        LWARNING("Height must be greater zero. Set to 1.0f.");
    }

    //resize vectors
    size_t startVertex = vertices_.size();
    vertices_.resize(startVertex+60);

    //create verticies
        //front
    addVertex(startVertex+ 0, tgt::vec3(-0.46f, 0.5f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 1, tgt::vec3( 0.46f, 0.5f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 2, tgt::vec3(-0.46f, 0.3f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 3, tgt::vec3( 0.22f, 0.3f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 4, tgt::vec3( 0.46f, 0.3f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 5, tgt::vec3(-0.46f,-0.3f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 6, tgt::vec3(-0.22f,-0.3f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 7, tgt::vec3( 0.46f,-0.3f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 8, tgt::vec3(-0.46f,-0.5f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
    addVertex(startVertex+ 9, tgt::vec3( 0.46f,-0.5f, 0.125f), tgt::vec3( 0.f, 0.f, 1.f),color);
        //back
    addVertex(startVertex+10, tgt::vec3(-0.46f, 0.5f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+11, tgt::vec3( 0.46f, 0.5f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+12, tgt::vec3(-0.46f, 0.3f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+13, tgt::vec3( 0.22f, 0.3f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+14, tgt::vec3( 0.46f, 0.3f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+15, tgt::vec3(-0.46f,-0.3f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+16, tgt::vec3(-0.22f,-0.3f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+17, tgt::vec3( 0.46f,-0.3f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+18, tgt::vec3(-0.46f,-0.5f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
    addVertex(startVertex+19, tgt::vec3( 0.46f,-0.5f,-0.125f), tgt::vec3( 0.f, 0.f,-1.f),color);
        //top
    addVertex(startVertex+20, tgt::vec3(-0.46f, 0.5f,-0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
    addVertex(startVertex+21, tgt::vec3( 0.46f, 0.5f,-0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
    addVertex(startVertex+22, tgt::vec3( 0.46f, 0.5f, 0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
    addVertex(startVertex+23, tgt::vec3(-0.46f, 0.5f, 0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
    addVertex(startVertex+24, tgt::vec3(-0.22f,-0.3f,-0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
    addVertex(startVertex+25, tgt::vec3( 0.46f,-0.3f,-0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
    addVertex(startVertex+26, tgt::vec3( 0.46f,-0.3f, 0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
    addVertex(startVertex+27, tgt::vec3(-0.22f,-0.3f, 0.125f), tgt::vec3( 0.f, 1.f, 0.f),color);
        //bottom
    addVertex(startVertex+28, tgt::vec3(-0.46f, 0.3f,-0.125f), tgt::vec3( 0.f,-1.f, 0.f),color);
    addVertex(startVertex+29, tgt::vec3( 0.22f, 0.3f,-0.125f), tgt::vec3( 0.f,-1.f, 0.f),color);
    addVertex(startVertex+30, tgt::vec3( 0.22f, 0.3f, 0.125f), tgt::vec3( 0.f,-1.f, 0.f),color);
    addVertex(startVertex+31, tgt::vec3(-0.46f, 0.3f, 0.125f), tgt::vec3( 0.f,-1.f, 0.f),color);
    addVertex(startVertex+32, tgt::vec3(-0.46f,-0.5f,-0.125f), tgt::vec3( 0.f,-1.f, 0.f),color);
    addVertex(startVertex+33, tgt::vec3( 0.46f,-0.5f,-0.125f), tgt::vec3( 0.f,-1.f, 0.f),color);
    addVertex(startVertex+34, tgt::vec3( 0.46f,-0.5f, 0.125f), tgt::vec3( 0.f,-1.f, 0.f),color);
    addVertex(startVertex+35, tgt::vec3(-0.46f,-0.5f, 0.125f), tgt::vec3( 0.f,-1.f, 0.f),color);
        //left
    addVertex(startVertex+36, tgt::vec3(-0.46f, 0.3f,-0.125f), tgt::vec3(-1.f, 0.f, 0.f),color);
    addVertex(startVertex+37, tgt::vec3(-0.46f, 0.5f,-0.125f), tgt::vec3(-1.f, 0.f, 0.f),color);
    addVertex(startVertex+38, tgt::vec3(-0.46f, 0.5f, 0.125f), tgt::vec3(-1.f, 0.f, 0.f),color);
    addVertex(startVertex+39, tgt::vec3(-0.46f, 0.3f, 0.125f), tgt::vec3(-1.f, 0.f, 0.f),color);
    addVertex(startVertex+40, tgt::vec3(-0.46f,-0.5f,-0.125f), tgt::vec3(-1.f, 0.f, 0.f),color);
    addVertex(startVertex+41, tgt::vec3(-0.46f,-0.3f,-0.125f), tgt::vec3(-1.f, 0.f, 0.f),color);
    addVertex(startVertex+42, tgt::vec3(-0.46f,-0.3f, 0.125f), tgt::vec3(-1.f, 0.f, 0.f),color);
    addVertex(startVertex+43, tgt::vec3(-0.46f,-0.5f, 0.125f), tgt::vec3(-1.f, 0.f, 0.f),color);
        //right
    addVertex(startVertex+44, tgt::vec3( 0.46f, 0.5f,-0.125f), tgt::vec3( 1.f, 0.f, 0.f),color);
    addVertex(startVertex+45, tgt::vec3( 0.46f, 0.3f,-0.125f), tgt::vec3( 1.f, 0.f, 0.f),color);
    addVertex(startVertex+46, tgt::vec3( 0.46f, 0.3f, 0.125f), tgt::vec3( 1.f, 0.f, 0.f),color);
    addVertex(startVertex+47, tgt::vec3( 0.46f, 0.5f, 0.125f), tgt::vec3( 1.f, 0.f, 0.f),color);
    addVertex(startVertex+48, tgt::vec3( 0.46f,-0.3f,-0.125f), tgt::vec3( 1.f, 0.f, 0.f),color);
    addVertex(startVertex+49, tgt::vec3( 0.46f,-0.5f,-0.125f), tgt::vec3( 1.f, 0.f, 0.f),color);
    addVertex(startVertex+50, tgt::vec3( 0.46f,-0.5f, 0.125f), tgt::vec3( 1.f, 0.f, 0.f),color);
    addVertex(startVertex+51, tgt::vec3( 0.46f,-0.3f, 0.125f), tgt::vec3( 1.f, 0.f, 0.f),color);
        //diagonal
    tgt::vec3 norm = tgt::normalize(tgt::vec3(-0.6f, 0.68f, 0.f));
    addVertex(startVertex+52, tgt::vec3(-0.46f,-0.3f,-0.125f), norm,color);
    addVertex(startVertex+53, tgt::vec3( 0.22f, 0.3f,-0.125f), norm,color);
    addVertex(startVertex+54, tgt::vec3( 0.22f, 0.3f, 0.125f), norm,color);
    addVertex(startVertex+55, tgt::vec3(-0.46f,-0.3f, 0.125f), norm,color);
              norm = tgt::normalize(tgt::vec3( 0.6f,-0.68f,0.f));
    addVertex(startVertex+56, tgt::vec3(-0.22f,-0.3f,-0.125f), norm,color);
    addVertex(startVertex+57, tgt::vec3( 0.46f, 0.3f,-0.125f), norm,color);
    addVertex(startVertex+58, tgt::vec3( 0.46f, 0.3f, 0.125f), norm,color);
    addVertex(startVertex+59, tgt::vec3(-0.22f,-0.3f, 0.125f), norm,color);


    //resize by height
    for(size_t i = startVertex; i < vertices_.size(); i++)
        vertices_[i].pos_ *= height;

    //create indices
        //front
    addTriangle(TriangleType(I(startVertex + 0), I(startVertex + 1), I(startVertex + 2)));
    addTriangle(TriangleType(I(startVertex + 1), I(startVertex + 4), I(startVertex + 2)));
    addTriangle(TriangleType(I(startVertex + 3), I(startVertex + 4), I(startVertex + 6)));
    addTriangle(TriangleType(I(startVertex + 6), I(startVertex + 5), I(startVertex + 3)));
    addTriangle(TriangleType(I(startVertex + 5), I(startVertex + 7), I(startVertex + 9)));
    addTriangle(TriangleType(I(startVertex + 9), I(startVertex + 8), I(startVertex + 5)));
        //back
    addTriangle(TriangleType(I(startVertex + 11), I(startVertex + 10), I(startVertex + 12)));
    addTriangle(TriangleType(I(startVertex + 12), I(startVertex + 14), I(startVertex + 11)));
    addTriangle(TriangleType(I(startVertex + 14), I(startVertex + 13), I(startVertex + 16)));
    addTriangle(TriangleType(I(startVertex + 13), I(startVertex + 15), I(startVertex + 16)));
    addTriangle(TriangleType(I(startVertex + 17), I(startVertex + 15), I(startVertex + 18)));
    addTriangle(TriangleType(I(startVertex + 18), I(startVertex + 19), I(startVertex + 17)));
        //top
    addTriangle(TriangleType(I(startVertex + 20), I(startVertex + 21), I(startVertex + 22)));
    addTriangle(TriangleType(I(startVertex + 22), I(startVertex + 23), I(startVertex + 20)));
    addTriangle(TriangleType(I(startVertex + 24), I(startVertex + 25), I(startVertex + 26)));
    addTriangle(TriangleType(I(startVertex + 26), I(startVertex + 27), I(startVertex + 24)));
        //bottom
    addTriangle(TriangleType(I(startVertex + 28), I(startVertex + 29), I(startVertex + 30)));
    addTriangle(TriangleType(I(startVertex + 30), I(startVertex + 31), I(startVertex + 28)));
    addTriangle(TriangleType(I(startVertex + 32), I(startVertex + 33), I(startVertex + 34)));
    addTriangle(TriangleType(I(startVertex + 34), I(startVertex + 35), I(startVertex + 32)));
        //left
    addTriangle(TriangleType(I(startVertex + 36), I(startVertex + 37), I(startVertex + 38)));
    addTriangle(TriangleType(I(startVertex + 38), I(startVertex + 39), I(startVertex + 36)));
    addTriangle(TriangleType(I(startVertex + 40), I(startVertex + 41), I(startVertex + 42)));
    addTriangle(TriangleType(I(startVertex + 42), I(startVertex + 43), I(startVertex + 40)));
        //right
    addTriangle(TriangleType(I(startVertex + 44), I(startVertex + 45), I(startVertex + 46)));
    addTriangle(TriangleType(I(startVertex + 46), I(startVertex + 47), I(startVertex + 44)));
    addTriangle(TriangleType(I(startVertex + 48), I(startVertex + 49), I(startVertex + 50)));
    addTriangle(TriangleType(I(startVertex + 50), I(startVertex + 51), I(startVertex + 48)));
        //diagonal
    addTriangle(TriangleType(I(startVertex + 52), I(startVertex + 53), I(startVertex + 54)));
    addTriangle(TriangleType(I(startVertex + 54), I(startVertex + 55), I(startVertex + 52)));
    addTriangle(TriangleType(I(startVertex + 56), I(startVertex + 57), I(startVertex + 58)));
    addTriangle(TriangleType(I(startVertex + 58), I(startVertex + 59), I(startVertex + 56)));

    invalidate();
}
//-------------------------------------------------------------------------------------------------
/// A triangle mesh with vertex type VertexBase with indices of type uint16_t <=> up to 2^16 triangles.
class VRN_CORE_API TriangleMeshGeometryUInt16IndexedSimple : public TriangleMeshGeometryIndexed<uint16_t, VertexBase> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryUInt16IndexedSimple(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryUInt16IndexedSimple"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return SIMPLE;}
    virtual TriangleMeshGeometryIndexedBase::IndexType getIndexType() const {return INDEX_UINT16;}
protected:
};

/// A triangle mesh with vertex type VertexColor with indices of type uint16_t <=> up to 2^16 triangles.
class VRN_CORE_API TriangleMeshGeometryUInt16IndexedColor : public TriangleMeshGeometryIndexed<uint16_t, VertexColor> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryUInt16IndexedColor(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryUInt16IndexedColor"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return COLOR; }
    virtual TriangleMeshGeometryIndexedBase::IndexType getIndexType() const {return INDEX_UINT16;}
protected:
};

/// A triangle mesh with vertex type VertexNormal with indices of type uint16_t <=> up to 2^16 triangles.
class VRN_CORE_API TriangleMeshGeometryUInt16IndexedNormal : public TriangleMeshGeometryIndexed<uint16_t, VertexNormal> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryUInt16IndexedNormal(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryUInt16IndexedNormal"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return NORMAL; }
    virtual TriangleMeshGeometryIndexedBase::IndexType getIndexType() const {return INDEX_UINT16;}
protected:
};

/// A triangle mesh with vertex type VertexTexCoord with indices of type uint16_t <=> up to 2^16 triangles.
class VRN_CORE_API TriangleMeshGeometryUInt16IndexedTexCoord : public TriangleMeshGeometryIndexed<uint16_t, VertexTexCoord> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryUInt16IndexedTexCoord(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryUInt16IndexedTexCoord"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return TEXCOORD; }
    virtual TriangleMeshGeometryIndexedBase::IndexType getIndexType() const {return INDEX_UINT16;}
protected:
};

/// A triangle mesh with vertex type VertexColorNormal with indices of type uint16_t <=> up to 2^16 triangles.
class VRN_CORE_API TriangleMeshGeometryUInt16IndexedColorNormal : public TriangleMeshGeometryIndexed<uint16_t, VertexColorNormal> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryUInt16IndexedColorNormal(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryUInt16IndexedColorNormal"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return COLOR_NORMAL; }
    virtual TriangleMeshGeometryIndexedBase::IndexType getIndexType() const {return INDEX_UINT16;}
protected:
};

/// A triangle mesh with vertex type VertexTextureNormal with indices of type uint16_t <=> up to 2^16 triangles.
class VRN_CORE_API TriangleMeshGeometryUInt16IndexedNormalTexCoord : public TriangleMeshGeometryIndexed<uint16_t, VertexNormalTexCoord> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryUInt16IndexedNormalTexCoord(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryUInt16IndexedNormaltexCoord"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return NORMAL_TEXCOORD; }
    virtual TriangleMeshGeometryIndexedBase::IndexType getIndexType() const {return INDEX_UINT16;}
protected:
};

/// A triangle mesh with vertex type VertexNormal with indices of type uint16_t <=> up to 2^16 triangles.
class VRN_CORE_API TriangleMeshGeometryUInt16IndexedColorTexCoord : public TriangleMeshGeometryIndexed<uint16_t, VertexColorTexCoord> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryUInt16IndexedColorTexCoord(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryUInt16IndexedColorTexCoord"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return COLOR_TEXCOORD; }
    virtual TriangleMeshGeometryIndexedBase::IndexType getIndexType() const {return INDEX_UINT16;}
protected:
};

/// A triangle mesh with vertex type VertexColorNormalTexCoord with indices of type uint16_t <=> up to 2^16 triangles.
class VRN_CORE_API TriangleMeshGeometryUInt16IndexedColorNormalTexCoord : public TriangleMeshGeometryIndexed<uint16_t, VertexColorNormalTexCoord> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryUInt16IndexedColorNormalTexCoord(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryUInt16IndexedColorNormalTexCoord"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return COLOR_NORMAL_TEXCOORD; }
    virtual TriangleMeshGeometryIndexedBase::IndexType getIndexType() const {return INDEX_UINT16;}
protected:
};

//-------------------------------------------------------------------------------------------------
/// A triangle mesh with vertex type VertexBase with indices of type uint32_t <=> up to 2^32 triangles.
class VRN_CORE_API TriangleMeshGeometryUInt32IndexedSimple : public TriangleMeshGeometryIndexed<uint32_t, VertexBase> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryUInt32IndexedSimple(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryUInt32IndexedSimple"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return SIMPLE; }
    virtual TriangleMeshGeometryIndexedBase::IndexType getIndexType() const {return INDEX_UINT32;}
protected:
};

/// A triangle mesh with vertex type VertexColor with indices of type uint32_t <=> up to 2^32 triangles.
class VRN_CORE_API TriangleMeshGeometryUInt32IndexedColor : public TriangleMeshGeometryIndexed<uint32_t, VertexColor> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryUInt32IndexedColor(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryUInt32IndexedColor"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return COLOR; }
    virtual TriangleMeshGeometryIndexedBase::IndexType getIndexType() const {return INDEX_UINT32;}
protected:
};

/// A triangle mesh with vertex type VertexNormal with indices of type uint32_t <=> up to 2^32 triangles.
class VRN_CORE_API TriangleMeshGeometryUInt32IndexedNormal : public TriangleMeshGeometryIndexed<uint32_t, VertexNormal> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryUInt32IndexedNormal(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryUInt32IndexedNormal"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return NORMAL; }
    virtual TriangleMeshGeometryIndexedBase::IndexType getIndexType() const {return INDEX_UINT32;}
protected:
};

/// A triangle mesh with vertex type VertexTexCoord with indices of type uint32_t <=> up to 2^32 triangles.
class VRN_CORE_API TriangleMeshGeometryUInt32IndexedTexCoord : public TriangleMeshGeometryIndexed<uint32_t, VertexTexCoord> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryUInt32IndexedTexCoord(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryUInt32IndexedTexCoord"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return TEXCOORD; }
    virtual TriangleMeshGeometryIndexedBase::IndexType getIndexType() const {return INDEX_UINT32;}
protected:
};


/// A triangle mesh with vertex type VertexColorNormal with indices of type uint32_t <=> up to 2^32 triangles.
class VRN_CORE_API TriangleMeshGeometryUInt32IndexedColorNormal : public TriangleMeshGeometryIndexed<uint32_t, VertexColorNormal> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryUInt32IndexedColorNormal(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryUInt32IndexedColorNormal"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return COLOR_NORMAL; }
    virtual TriangleMeshGeometryIndexedBase::IndexType getIndexType() const {return INDEX_UINT32;}
protected:
};

/// A triangle mesh with vertex type VertexTextureNormal with indices of type uint32_t <=> up to 2^32 triangles.
class VRN_CORE_API TriangleMeshGeometryUInt32IndexedNormalTexCoord : public TriangleMeshGeometryIndexed<uint32_t, VertexNormalTexCoord> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryUInt32IndexedNormalTexCoord(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryUInt32IndexedNormaltexCoord"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return NORMAL_TEXCOORD; }
    virtual TriangleMeshGeometryIndexedBase::IndexType getIndexType() const {return INDEX_UINT32;}
protected:
};

/// A triangle mesh with vertex type VertexNormal with indices of type uint32_t <=> up to 2^32 triangles.
class VRN_CORE_API TriangleMeshGeometryUInt32IndexedColorTexCoord : public TriangleMeshGeometryIndexed<uint32_t, VertexColorTexCoord> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryUInt32IndexedColorTexCoord(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryUInt32IndexedColorTexCoord"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return COLOR_TEXCOORD; }
    virtual TriangleMeshGeometryIndexedBase::IndexType getIndexType() const {return INDEX_UINT32;}
protected:
};

/// A triangle mesh with vertex type VertexColorNormalTexCoord with indices of type uint32_t <=> up to 2^32 triangles.
class VRN_CORE_API TriangleMeshGeometryUInt32IndexedColorNormalTexCoord : public TriangleMeshGeometryIndexed<uint32_t, VertexColorNormalTexCoord> {
public:
    virtual Geometry* create() const { return new TriangleMeshGeometryUInt32IndexedColorNormalTexCoord(); }
    virtual std::string getClassName() const { return "TriangleMeshGeometryUInt32IndexedColorNormalTexCoord"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return COLOR_NORMAL_TEXCOORD; }
    virtual TriangleMeshGeometryIndexedBase::IndexType getIndexType() const {return INDEX_UINT32;}
protected:
};




} // namespace

#endif  //VRN_TRIANGLEMESHGEOMETRY_H
