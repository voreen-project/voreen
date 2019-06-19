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

#ifndef VRN_GLMESHGEOMETRY_H
#define VRN_GLMESHGEOMETRY_H

#include "geometry.h"
#include "vertex.h"

#include "tgt/tgt_gl.h"
#include "tgt/tgt_math.h"
#include "tgt/glmath.h"
#include "tgt/texturemanager.h"
#include "tgt/texture2darray.h"
#include "tgt/immediatemode/immediatemode.h"

#include "voreen/core/utils/glsl.h"

#include "voreen/core/io/serialization/xmlserializer.h"
#include "voreen/core/io/serialization/xmldeserializer.h"

namespace voreen {

/**
 * Geometry class for OpenGL rendering. Supports both indexed and non-indexed rendering and all core profile primitive types.
 */
class VRN_CORE_API GlMeshGeometryBase : public Geometry {

public:

    /// if indexing is used, the indices can be of different types
    enum IndexType {
        INDEX_UINT16,
        INDEX_UINT32
    };

    GlMeshGeometryBase();
    virtual ~GlMeshGeometryBase();

    //----------------- base class methods
    bool supportsNormals() const;
    bool supportsColors() const;
    bool supportsTextureData() const;

    /**
     * Returns a boolean type which is true if the render()-method uses indexed drawing.
     *
     * By default, indexed drawing is disabled. It will be enabled when adding the first index to the data, and disabled either manually or by clearing the index data.
     */
    bool usesIndexedDrawing() const;

    /**
     * If previously manually disabled, this method will re-enable indexed drawing (only works if getNumIndices() > 0).
     */
    void enableIndexedDrawing();

    /**
     * This will disable indexed drawing until either enableIndexedDrawing() is called or a new index is added / new index data is set.
     * This can be used to (temporarily) ignore the index data, e.g. for rendering the vertices as points.
     */
    void disableIndexedDrawing();

    // texture data stuff
    void loadTextureData();
    void loadTextureData(const std::string& filename);
    void loadTextureData(const std::vector<std::string>& filenames);
    void setTextureData(tgt::Texture* tex, bool owner = false);
    void clearTextureData();

    tgt::Texture* getTextureData() const;
    bool hasTextureData() const;

    /// Serializes the vertices and indices as binary blob, additionally serializes texture filenames (serialization of vertices and indices is handled in the subclass!)
    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& d);

    /**
     * Generate Shader defines for the VertexLayout of this object.
     * Returns a string that sets the shader defines for using colors, normals, and/or texture coordinates.
     *
     * @param useTransparency Determines whether defines should prepare transparency information
     */
    std::string getShaderDefines(bool useTransparency = false) const;

    /// returns the primitive type set as standard for rendering (e.g., GL_TRIANGLE_STRIP, etc.)
    virtual GLenum getPrimitiveType() const;

    /// sets the primitive type that is used as default for rendering
    virtual void setPrimitiveType(GLenum primitiveType);

    // TODO: implement clipping (currently outputs an error message)
    //virtual void clip(const tgt::plane& clipPlane, double epsilon = 1e-5);


    //--------------------  methods to be implemented in subclass

    /// Returns the bounding box in model or transformed coordinates. The BB is cached internally.
    virtual tgt::Bounds getBoundingBox(bool transformed = true) const = 0;

    virtual IndexType getIndexType() const = 0;

    virtual size_t getNumVertices() const = 0;
    virtual size_t getNumIndices() const = 0;

    /// this returns true if the vertices are empty, even if there are elements in the index list
    virtual bool isEmpty() const = 0;

    /**
     * Clears both the vertices and the indices. Will also disable indexed drawing (until new indices are added).
     */
    virtual void clear() = 0;

    /**
     * Different vertex attribute layouts are supported. However, the order and vertex attribute locations in the shader are fixed as follows:
     *
     * layout(location = 0) in vec4 vertexPosition;
     * // when using colors:
     * layout(location = 1) in vec4 vertexColor;
     * // when using normals:
     * layout(location = 2) in vec3 vertexNormal;
     * // when using texture coordinates:
     * layout(location = 3) in vec2 vertexTexCoord;
     */
    virtual VertexBase::VertexLayout getVertexLayout() const = 0;

    virtual void addVertex(tgt::vec3 pos, tgt::vec3 normal = tgt::vec3(0.f, 0.f, 1.f), tgt::vec4 color = tgt::vec4::one, tgt::vec2 texCoord = tgt::vec2::zero, tgt::ivec2 texIndex = tgt::ivec2::zero) = 0;

    /**
     * Renders the mesh using OpenGL by binding the buffer, settings appropriate vertex attribute pointers and calling glDrawArrays or glDrawElements, depending on if an index buffer is used.
     *
     * A shader with appropriate vertex attribute bindings has to be activated before calling render().
     */
    virtual void render() const = 0;

    /*
     * Renders the mesh using OpenGL by binding the buffer, setting appropriate vertex attribute pointers and calling lDrawArrays or glDrawElements, depending on if an index buffer is used.
     *
     * A shader with appropriate vertexattribute bindings has to be activated before calling render().
     *
     * @param primitiveType allows to set a specific primitive type for the rendering (e.g., GL_TRIANGLES, GL_TRIANGLE_STRIP)
     */
    virtual void render(GLenum primitiveType) const = 0;

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

    /**
     * Generate Shader defines for the given VertexLayout.
     * Returns a string that sets the shader defines for using colors, normals, and/or texture coordinates.
     */
    static std::string getShaderDefines(VertexBase::VertexLayout layout, bool transparent);

    static bool supportsNormals(VertexBase::VertexLayout layout);
    static bool supportsColors(VertexBase::VertexLayout layout);
    static bool supportsTextureData(VertexBase::VertexLayout layout);

    //-------------- texture stuff

    void setTextureDataAndFilenames(tgt::Texture* tex, const std::vector<std::string>& filenames, bool owner = false);
    tgt::Texture* textureData_;
    std::vector<std::string> textureFilenames_;
    bool textureOwner_;

    //-------------- OpenGL type stuff

    GLenum primitiveType_; ///< primitive type which is used when calling render(), default is GL_TRIANGLES

    bool usesIndexedDrawing_;

    //-------------- transparency

    bool hasTransparentComponents_; ///< Used to indicate whether this geometry contains transparent components (e.g. transparent vertex colours, transparent regions in textures etc.)
};

// implementation of base class: see .cpp file


///////////// TEMPLATE SUBCLASSES ////////////////////////////////////////

/**
 * Template class for implementing the GlMeshGeometry with various vertex types.
 * The first template argument ist the type of the index buffer (uint16_t or uint32_t), the second one is the VertexType (see vertex.h).
 *
 * Although the template provides most of the implementation, its specialized concrete subclasses should be used.
 */
template <class I, class V>
class VRN_CORE_API GlMeshGeometry : public GlMeshGeometryBase {

public:
    typedef V VertexType;

    GlMeshGeometry();

    /// Clears the vector and deletes the OpenGL buffers if necessary.
    virtual ~GlMeshGeometry();

    virtual std::unique_ptr<Geometry> clone() const;

    /**
     * Returns true, if the passed Geometry is a GlMeshGeometry of the same type and all its vertices and indices are equal to this one's.
     */
    virtual bool equals(const Geometry* geometry, double epsilon = 1e-5) const;

    //---------------------- primitive restart for indexed geometries

    /// returns the primitive restart index (which is fixed in this class currently to prevent problems with adding new geometry)
    virtual I getPrimitiveRestartIndex() const;

    /// enables the primitive restart using the fixed index for this mesh (default: disabled)
    virtual void enablePrimitiveRestart();

    /// disables the primitive restart for this mesh (default behavior)
    virtual void disablePrimitiveRestart();

    /// returns a bool if primitive restart is enabled for this mesh
    virtual bool isPrimitiveRestartEnabled() const;

    //------------------- vertex and index buffer methods

    /// adds a vertex at the end of the vertex list
    virtual void addVertex(tgt::vec3 pos, tgt::vec3 normal = tgt::vec3(0.f, 0.f, 1.f), tgt::vec4 color = tgt::vec4::one, tgt::vec2 texCoord = tgt::vec2::zero, tgt::ivec2 texIndex = tgt::ivec2::zero);

    /// adds a vertex with the correct template type at the end of the vertex list
    virtual void addVertex(const VertexType& v);

    /**
     * Get the vertex at index i.
     *
     * Caution: i must be < numVertices()
     */
    const VertexType& getVertex(size_t i) const;

    /**
     * Overwrite the vertex at index i.
     *
     * Caution: i must be < numVertices()
     */
    void setVertex(size_t i, const VertexType& vertex) const;

    /// Used to pass the vertex set for this geometry. This will not (!) touch the index buffer.
    void setVertices(const std::vector<VertexType>& vertices);

    const std::vector<VertexType>& getVertices() const;

    /// adds an index at the end of the index list, also automatically activates indexed rendering (can be deactivated manually)
    void addIndex(I index);

    /// Used to pass a list of indices to this geometry. Does not change the vertex data itself. Automatically activates indexed rendering if the index list is not empty, and deactivates it otherwise.
    void setIndices(const std::vector<I>& indices);

    const std::vector<I>& getIndices() const;

    /// Returns the OpenGL data type of the indices (has to be implemented by special subclasses)
    virtual GLenum getIndexDataType() const = 0;

    //------------------- methods implemented from base class

    /// Returns the bounding box in model or transformed coordinates. The BB is cached internally.
    virtual tgt::Bounds getBoundingBox(bool transformed = true) const;

    virtual size_t getNumVertices() const;
    virtual size_t getNumIndices() const;

    /// this returns true if the vertices are empty, even if there are elements in the index list
    virtual bool isEmpty() const;

    /**
     * Clears both the vertices and the indices. Will also disable indexed drawing (until new indices are added).
     */
    virtual void clear();

    /**
     * Renders the mesh using OpenGL by binding the buffer, settings appropriate vertex attribute pointers and calling glDrawArrays or glDrawElements, depending on if an index buffer is used.
     *
     * A shader with appropriate vertex attribute bindings can be activated before calling render().
     * Otherwise the mesh will be rendered with a default shader.
     */
    virtual void render() const;

    /**
     * Renders the mesh using OpenGL by binding the buffer, setting appropriate vertex attribute pointers and calling glDrawArrays or glDrawElements, depending on if an index buffer is used.
     *
     * A shader with appropriate vertexattribute bindings has to be activated before calling render().
     * Otherwise the mesh will be rendered with a default shader.
     *
     * @param primitiveType allows to set a specific primitive type for the rendering (e.g., GL_TRIANGLES, GL_TRIANGLE_STRIP)
     */
    virtual void render(GLenum primitiveType) const;

    /**
     * Renders the mesh by binding a default shader for _this_ class and calling render() on the mesh.
     *
     * The VertexLayout of the mesh has to offer a superset of features that this class' VertexLayout offers.
     *
     * @param mesh The mesh to be rendered.
     * @param useTransparency if enabled, an OIT shader is used for rendering
     */
    template<class MeshV>
    static void renderDefault(const GlMeshGeometry<I, MeshV>& mesh, bool useTransparency = false);

    /**
     * Renders the mesh by binding a default shader for _this_ class and calling render() on the mesh.
     *
     * The VertexLayout of the mesh has to offer a superset of features that this class' VertexLayout offers.
     *
     * @param mesh The mesh to be rendered.
     * @param primitiveType allows to set a specific primitive type for the rendering (e.g., GL_TRIANGLES, GL_TRIANGLE_STRIP)
     * @param useTransparency if enabled, an OIT shader is used for rendering
     */
    template<class MeshVertexLayout>
    static void renderDefault(const GlMeshGeometry<I, MeshVertexLayout>& mesh, GLenum primitiveType, bool useTransparency = false);

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& d);

    virtual VertexBase::VertexLayout getVertexLayout() const;

    void writeObj(std::ostream& s) const;

    // Method is not implemented here to force use of subclasses
    //virtual IndexType getIndexType() const = 0;

    //------------------- member functions for mesh creation

    /// Clears the mesh and creates a sphere mesh with the specified parameters (using index buffers and GL_TRIANGLE_STRIP)
    void setSphereGeometry(float radius, tgt::vec3 center, tgt::vec4 color, size_t numSlicesAndTiles);

    /// Clears the mesh and creates a cylinder mesh with the specified parameters (using index buffers and GL_TRIANGLES)
    void setCylinderGeometry(tgt::vec4 color, float lowerRadius, float upperRadius, float height, size_t slices, size_t tiles, bool lowerCap, bool upperCap);

    // The following mesh creation functions are
    // analogous to those in tgt/quadrics.h
    // All meshes are centered at (0,0,0).

    // 2D meshes:

    /// Clears the mesh and creates a disk in the xy-plane with the specified parameters.
    void setDiskGeometry(float innerRadius, float outerRadius, size_t numSlices);

    /// Clears the mesh and creates an equilateral triangle in the xy-plane with the specified edge length.
    void setTriangleGeometry(float edgeLength);

    /// Clears the mesh and creates a rectangle in the xy-plane with the specified parameters.
    void setRectangleGeometry(float width, float height);

    // 3D meshes:

    /// Clears the mesh and creates a tetrahedron with the specified edge length.
    void setTetraHedronGeometry(float edgeLength);

    /// Clears the mesh and creates a cuboid in with the specified parameters.
    void setCuboidGeometry(float width, float height, float depth);

protected:

    /// helper method for setCylinderGeometry(...)
    void createCylinderIndices(size_t tiles, size_t slices, size_t offset, std::vector<I>& indices);

    void initialize() const;
    virtual void updateBoundingBox() const;

    /// Updates the OpenGL buffer.
    virtual void updateBuffer() const;

    /// Flags the bounding box and OpenGL buffer as invalid.
    void invalidate();

    mutable tgt::Bounds boundingBox_;
    mutable bool boundingBoxValid_;

    mutable GLuint vertexArrayID_;
    mutable GLuint bufferObjectID_;
    mutable bool bufferObjectValid_;

    mutable GLuint indexObject_;
    mutable bool indexObjectValid_;

    mutable bool isInitialized_;

    I primitiveRestartIndex_;
    bool primitiveRestartEnabled_;  ///< primitive restart index is disabled by default

    std::vector<I> indices_;
    std::vector<VertexType> vertices_;

    static const std::string loggerCat_;

    /// Get the default shader for the VertexLayout of this class
    static tgt::Shader& getDefaultShader();
    static tgt::Shader& getDefaultShaderTransparent();

private:

    //Static shaders shared for objects of this particular GlMeshGeometry types
    static tgt::Shader* defaultShader_;
    static tgt::Shader* defaultShaderTransparent_;
};


//////////////////// IMPLEMENTATION /////////////////////////////////////////////////////////

template <class I, class V>
const std::string GlMeshGeometry<I,V>::loggerCat_ = "voreen.GlMeshGeometry";

template <class I, class V>
tgt::Shader* GlMeshGeometry<I,V>::defaultShader_ = nullptr;

template <class I, class V>
tgt::Shader* GlMeshGeometry<I,V>::defaultShaderTransparent_ = nullptr;

template <class I, class V>
GlMeshGeometry<I, V>::GlMeshGeometry()
    : GlMeshGeometryBase()
    , boundingBox_()
    , boundingBoxValid_(false)
    , bufferObjectID_(0)
    , vertexArrayID_(0)
    , bufferObjectValid_(false)
    , indexObject_(0)
    , indexObjectValid_(false)
    , isInitialized_(false)
    , primitiveRestartIndex_(std::numeric_limits<I>::max())
    , primitiveRestartEnabled_(false)
{

}

template <class I, class V>
void GlMeshGeometry<I, V>::initialize() const {
    tgtAssert(!isInitialized_, "initialize is called twice");
    glGenVertexArrays(1, &vertexArrayID_);
    glGenBuffers(1, &bufferObjectID_);
    glGenBuffers(1, &indexObject_);
    isInitialized_ = true;
}

template <class I, class V>
GlMeshGeometry<I, V>::~GlMeshGeometry() {
    clear();

    if (isInitialized_) {
        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER,0);
        glDeleteBuffers(1, &bufferObjectID_);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glDeleteBuffers(1, &indexObject_);
        glDeleteVertexArrays(1, &vertexArrayID_);
    }
}

template <class I, class V>
bool GlMeshGeometry<I, V>::equals(const Geometry* geometry, double epsilon) const {
    const GlMeshGeometry<I, V>* triMesh = dynamic_cast<const GlMeshGeometry<I, V>*>(geometry);

    if(triMesh) {
        // check number of vertices, indices, and primitive type
        if (getNumVertices() != triMesh->getNumVertices())
            return false;

        if (getNumIndices() != triMesh->getNumIndices())
            return false;

        if (getPrimitiveType() != triMesh->getPrimitiveType())
            return false;

        // check individual indices and vertices
        for (size_t i = 0; i < getNumIndices(); ++i) {
            if (getIndices().at(i) != triMesh->getIndices().at(i))
                return false;
        }

        for (size_t i = 0; i < getNumVertices(); ++i) {
            if (!getVertex(i).equals(triMesh->getVertex(i), epsilon))
                return false;
        }

        return true;
    }
    else
        return false;
}

template <class I, class V>
void GlMeshGeometry<I, V>::updateBoundingBox() const {
    tgt::Bounds bb;
    for(auto& it : vertices_)
        bb.addPoint(it.pos_);
    boundingBox_ = bb;
    boundingBoxValid_ = true;
}

template <class I, class V>
void GlMeshGeometry<I, V>::updateBuffer() const {
    if(bufferObjectValid_ && indexObjectValid_)
        return;

    if(isEmpty())
        return;

    if(!isInitialized_)
        initialize();

    glBindVertexArray(vertexArrayID_);
    glBindBuffer(GL_ARRAY_BUFFER, bufferObjectID_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(VertexType) * vertices_.size(), vertices_.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexObject_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(I) * indices_.size(), indices_.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    bufferObjectValid_ = true;
    indexObjectValid_ = true;
}

template<class I, class V>
std::unique_ptr<Geometry> GlMeshGeometry<I, V>::clone() const {
    GlMeshGeometry<I, V>* newGeom = static_cast<GlMeshGeometry<I, V>*>(create());
    newGeom->vertices_ = vertices_;
    newGeom->indices_ = indices_;
    newGeom->setTransformationMatrix(getTransformationMatrix());
    newGeom->getMetaDataContainer() = MetaDataContainer(getMetaDataContainer());
    newGeom->setTextureDataAndFilenames(textureData_, textureFilenames_, false);

    newGeom->setPrimitiveType(primitiveType_);
    newGeom->usesIndexedDrawing_ = usesIndexedDrawing_;

    return std::unique_ptr<Geometry>(newGeom);
}

template <class I, class V>
I GlMeshGeometry<I, V>::getPrimitiveRestartIndex() const {
    return primitiveRestartIndex_;
}

template <class I, class V>
void GlMeshGeometry<I, V>::enablePrimitiveRestart() {
    primitiveRestartEnabled_ = true;
}

template <class I, class V>
void GlMeshGeometry<I, V>::disablePrimitiveRestart() {
    primitiveRestartEnabled_ = false;
}

template <class I, class V>
bool GlMeshGeometry<I, V>::isPrimitiveRestartEnabled() const {
    return primitiveRestartEnabled_;
}

template <class I, class V>
void GlMeshGeometry<I, V>::invalidate() {
    boundingBoxValid_ = false;
    bufferObjectValid_ = false;
    indexObjectValid_ = false;
}

template <class I, class V>
void GlMeshGeometry<I,V>::addVertex(tgt::vec3 pos, tgt::vec3 normal, tgt::vec4 color, tgt::vec2 texCoord, tgt::ivec2 texIndex) {
    tgtAssert(vertices_.size() + 1 < static_cast<size_t>(getPrimitiveRestartIndex()), "Too many vertices for numeric range of index type");

    vertices_.push_back(V(pos,color,normal,texCoord,texIndex));
}

template <class I, class V>
void GlMeshGeometry<I,V>::addVertex(const V& v) {
    tgtAssert(vertices_.size() + 1 < static_cast<size_t>(getPrimitiveRestartIndex()), "Too many vertices for numeric range of index type");

    vertices_.push_back(v);
}

template <class I, class V>
const V& GlMeshGeometry<I,V>::getVertex(size_t i) const {
    tgtAssert(i < getNumVertices(), "invalid vertex index");

    return vertices_.at(i);
}

template <class I, class V>
void GlMeshGeometry<I,V>::setVertex(size_t i, const VertexType& vertex) const {
    tgtAssert(i < getNumVertices(), "invalid vertex index");

    vertices_.at(i) = vertex;
    invalidate();
}

template <class I, class V>
void GlMeshGeometry<I,V>::setVertices(const std::vector<VertexType>& vertices) {
    vertices_ = vertices;
    invalidate();
}

template <class I, class V>
const std::vector<V>& GlMeshGeometry<I,V>::getVertices() const {
    return vertices_;
}

template <class I, class V>
void GlMeshGeometry<I,V>::addIndex(I index) {
    indices_.push_back(index);
    usesIndexedDrawing_ = true;
    invalidate();
}

template <class I, class V>
void GlMeshGeometry<I,V>::setIndices(const std::vector<I>& indices) {
    indices_ = indices;
    usesIndexedDrawing_ = !(indices_.empty());
    invalidate();
}

template <class I, class V>
const std::vector<I>& GlMeshGeometry<I, V>::getIndices() const {
    return indices_;
}

template <class I, class V>
tgt::Shader& GlMeshGeometry<I, V>::getDefaultShader() {
    if(!defaultShader_) {
        tgtAssert(ShdrMgr.isInited(), "Shadermanager not initialized")
        tgt::Shader* newShader = ShdrMgr.loadSeparate("meshgeometrydefault.vert", "meshgeometrydefault.frag", getShaderDefines(V::layout, false), false);
        ShdrMgr.registerStaticShader(newShader);
        defaultShader_ = newShader;
    }
    tgtAssert(defaultShader_, "Loading mesh shader failed");
    return *defaultShader_;
}

template <class I, class V>
tgt::Shader& GlMeshGeometry<I, V>::getDefaultShaderTransparent() {
    if(!defaultShaderTransparent_) {
        tgtAssert(ShdrMgr.isInited(), "Shadermanager not initialized")
        tgt::Shader* newShader = ShdrMgr.loadSeparate("meshgeometrydefault.vert", "meshgeometrydefault.frag", getShaderDefines(V::layout, true), false);
        ShdrMgr.registerStaticShader(newShader);
        defaultShaderTransparent_ = newShader;
    }
    tgtAssert(defaultShaderTransparent_, "Loading mesh shader failed");
    return *defaultShaderTransparent_;
}

namespace {
    void writeVec4AsObj(std::ostream& s, const std::string& prefix, const tgt::vec4& v) {
        s << prefix << " " << v.x << " " << v.y << " " << v.z << " " << v.w << "\n";
    }
    void writeVec3AsObj(std::ostream& s, const std::string& prefix, const tgt::vec3& v) {
        s << prefix << " " << v.x << " " << v.y << " " << v.z << " " << "\n";
    }
    void writeVec2AsObj(std::ostream& s, const std::string& prefix, const tgt::vec2& v) {
        s << prefix << " " << v.x << " " << v.y << "\n";
    }
    template <class I, class V>
    void writeIndex(std::ostream& s, I i) {
        switch(V::layout & (VertexBase::TEXCOORD | VertexBase::NORMAL)) {
            case VertexBase::SIMPLE: {
                s << i;
                break;
            }
            case VertexBase::TEXCOORD: {
                s << i << "/" << i;
                break;
            }
            case VertexBase::NORMAL: {
                s << i << "//" << i;
                break;
            }
            case VertexBase::NORMAL_TEXCOORD: {
                s << i << "/" << i << "/" << i;
                break;
            }
            default:
                {
                    tgtAssert(false, "Invalid layout");
                }
        }
    }
}

template <class I, class V>
void GlMeshGeometry<I,V>::writeObj(std::ostream& s) const {
    // Define all positions
    for(auto& v : vertices_) {
        writeVec4AsObj(s, "v", v.pos_);
    }

    // Define texcoords
    if((V::layout & VertexBase::TEXCOORD) == VertexBase::TEXCOORD) {
        for(auto& v : vertices_) {
            writeVec2AsObj(s, "v", v.getTexCoord());
        }
    }

    // Define normals
    if((V::layout & VertexBase::NORMAL) == VertexBase::NORMAL) {
        for(auto& v : vertices_) {
            writeVec3AsObj(s, "v", v.getNormal());
        }
    }

    // Colors are not supported in obj files :-(

    //Write faces via indices
    switch(getPrimitiveType()) {
        case GL_TRIANGLES:
            {
                std::vector<I> collected;
                for(auto& i : indices_) {
                    if(i == primitiveRestartIndex_) {
                        collected.clear();
                    } else {
                        collected.push_back(i);
                    }
                    tgtAssert(collected.size() <= 3, "Invalid number of collected indices");
                    if(collected.size() == 3) {
                        s << "f";
                        writeIndex(s, collected[0]);
                        writeIndex(s, collected[1]);
                        writeIndex(s, collected[2]);

                        collected.clear();
                    }
                }
                tgtAssert(collected.empty(), "Remaining indices");
                break;
            }
        case GL_TRIANGLE_STRIP:
            {
                std::deque<I> collected;
                for(auto& i : indices_) {
                    if(i == primitiveRestartIndex_) {
                        collected.clear();
                    } else {
                        collected.push_back(i);
                    }
                    tgtAssert(collected.size() <= 3, "Invalid number of collected indices");
                    if(collected.size() == 3) {
                        s << "f";
                        writeIndex(s, collected[0]);
                        writeIndex(s, collected[1]);
                        writeIndex(s, collected[2]);

                        collected.pop_front();
                    }
                }
                break;
            }
        case GL_TRIANGLE_FAN:
            {
                std::deque<I> collected;
                for(auto& i : indices_) {
                    if(i == primitiveRestartIndex_) {
                        collected.clear();
                    } else {
                        collected.push_back(i);
                    }
                    tgtAssert(collected.size() <= 3, "Invalid number of collected indices");
                    if(collected.size() == 3) {
                        s << "f";
                        writeIndex(s, collected[0]);
                        writeIndex(s, collected[1]);
                        writeIndex(s, collected[2]);

                        // Not the most efficient, but the simplest for now...
                        I first = collected.front();
                        collected.pop_front();
                        collected.pop_front();
                        collected.push_front(first);
                    }
                }
                break;
            }
        default:
            {
                throw tgt::Exception("Unsupported geometry type for .obj export");
            }
    }
}

//------------------------ implemented from base class


template <class I, class V>
VertexBase::VertexLayout GlMeshGeometry<I,V>::getVertexLayout() const {
    return V::layout;
}

//virtual IndexType getIndexType() const = 0; // implemented in specific subclasses

template <class I, class V>
size_t GlMeshGeometry<I, V>::getNumVertices() const {
    return vertices_.size();
}

template <class I, class V>
size_t GlMeshGeometry<I, V>::getNumIndices() const {
    return indices_.size();
}

template <class I, class V>
bool GlMeshGeometry<I, V>::isEmpty() const {
    return (vertices_.size() == 0);
}

template <class I, class V>
void GlMeshGeometry<I, V>::clear() {
    vertices_.clear();
    indices_.clear();

    usesIndexedDrawing_ = false;
    primitiveRestartEnabled_ = false;

    invalidate();
}

template <class I, class V>
void GlMeshGeometry<I, V>::render() const {
    render(primitiveType_);
}

template <class I, class V>
void GlMeshGeometry<I, V>::render(GLenum primitiveType) const {

    // If there is no shader bound, we use the default one
    if(tgt::Shader::getCurrentProgram() == 0) {
        renderDefault(*this, primitiveType);
        return;
    }

    updateBuffer();

    if(isEmpty())
        return;

    glBindVertexArray(vertexArrayID_);
    glBindBuffer(GL_ARRAY_BUFFER, bufferObjectID_);
    VertexType::setupVertexAttributePointers();

    if (usesIndexedDrawing_) {

        if (primitiveRestartEnabled_) {
            glEnable(GL_PRIMITIVE_RESTART);
            glPrimitiveRestartIndex(static_cast<GLuint>(getPrimitiveRestartIndex()));
        }

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexObject_);
        glDrawElements(primitiveType, static_cast<GLsizei>(indices_.size()), getIndexDataType(), NULL);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glDisable(GL_PRIMITIVE_RESTART);
    }
    else {
        glDrawArrays(primitiveType, 0, static_cast<GLsizei>(vertices_.size()));
    }

    LGL_ERROR;

    VertexType::disableVertexAttributePointers();

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    LGL_ERROR;
}

template <class I, class V>
template <class MeshV>
void GlMeshGeometry<I,V>::renderDefault(const GlMeshGeometry<I, MeshV>& mesh, bool useTransparency) {
    renderDefault(mesh, mesh.getPrimitiveType(), useTransparency);
}

template <class I, class V>
template <class MeshV>
void GlMeshGeometry<I,V>::renderDefault(const GlMeshGeometry<I, MeshV>& mesh, GLenum primitiveType, bool useTransparency) {
    static_assert((V::layout | MeshV::layout) == MeshV::layout, "VertexLayout of mesh has to be a superset of the VertexLayout of this class.");

    // The layout considered for rendering
    VertexBase::VertexLayout layout = V::layout;

    tgt::Shader& shader = useTransparency ? getDefaultShaderTransparent() : getDefaultShader();
    shader.activate();

    tgt::mat4 modelView(MatStack.getModelViewMatrix());
    shader.setUniform("modelViewMatrix_", modelView);
    shader.setUniform("projectionMatrix_", MatStack.getProjectionMatrix());


    if(!supportsColors(layout) && !supportsTextureData(layout)) {
        shader.setUniform("solidColor_", IMode.getCurrentColor());
    }

    if(supportsNormals(layout)) {
        tgt::mat3 normalMatrix;
        if(!modelView.getRotationalPartMat3().invert(normalMatrix)) {
            LWARNING("Could not generate normal matrix out of current view / model matrix, using identity.");
            normalMatrix = tgt::mat3::identity;
        }
        normalMatrix = transpose(normalMatrix);
        shader.setUniform("normalMatrix_", normalMatrix);

        shader.setUniform("enableLighting_", true);

        tgt::ImmediateMode::LightSource lightSource = IMode.getLightSource();
        shader.setUniform("lightPositionEye_", lightSource.position.xyz());
        shader.setUniform("lightSource_.ambientColor_", lightSource.ambientColor);
        shader.setUniform("lightSource_.diffuseColor_", lightSource.diffuseColor);
        shader.setUniform("lightSource_.specularColor_", lightSource.specularColor);
        //shader.setUniform("lightSource_.attenuation_", tgt::vec3(1.f, 0.f, 0.f));

        tgt::ImmediateMode::Material material = IMode.getMaterial();
        shader.setUniform("shininess_", material.shininess);
    }

    tgt::TextureUnit texUnit;
    if(supportsTextureData(layout)) {
        if(mesh.hasTextureData()) {
            //If this mesh has a texture, use it.
            texUnit.activate();
            tgt::Texture* texture = mesh.getTextureData();
            tgtAssert(texture, "No texture");
            texture->enable();
            texture->bind();

            // Will be used for inactive samplers.
            tgt::TextureUnit unusedTexUnit;

            // Check if texture is a simple texture or a texture array.
            if(tgt::Texture2DArray* arrayTexture = dynamic_cast<tgt::Texture2DArray*>(texture)) {
                // Set and configure the shader to use a 2D texture array:
                shader.setUniform("useTextureArray_", true);
                shader.setUniform("texture_", unusedTexUnit.getUnitNumber());
                shader.setUniform("textures_", texUnit.getUnitNumber());
            } else {
                // Set and configure the shader to use a single texture:
                shader.setUniform("useTextureArray_", false);
                shader.setUniform("texture_", texUnit.getUnitNumber());
                shader.setUniform("textures_", unusedTexUnit.getUnitNumber());
            }

        } else if(IMode.textureMode() == tgt::ImmediateMode::TEX2D) {
            // Otherwise behave like IMode:
            // Use a single texture currently bound to texture unit zero.
            shader.setUniform("useTextureArray_", false);
            shader.setUniform("texture_", 0);
            shader.setUniform("textures_", 1 /* Unused, but valid */);
        } else {
            LERROR("Mesh does not have a texture and TEX2D mode is not selected.");
            return;
        }
    }
    LGL_ERROR;

    mesh.render(primitiveType);
    LGL_ERROR;

    shader.deactivate();

    if(supportsTextureData(layout)) {
        if(mesh.hasTextureData()) {
            tgt::Texture* texture = mesh.getTextureData();
            tgtAssert(texture, "No texture");
            texture->disable();
        }
    }
    LGL_ERROR;
}

template <class I, class V>
tgt::Bounds GlMeshGeometry<I, V>::getBoundingBox(bool transformed) const {
    if(!boundingBoxValid_)
        updateBoundingBox();

    if(transformed)
        return boundingBox_.transform(getTransformationMatrix());
    else
        return boundingBox_;
}

template <class I, class V>
void GlMeshGeometry<I, V>::serialize(Serializer& s) const {
    GlMeshGeometryBase::serialize(s);
    s.serializeBinaryBlob("vertices", vertices_);
    s.serializeBinaryBlob("indices", indices_);

    s.serialize("primitiveRestartEnabled", primitiveRestartEnabled_);
    s.serialize("primitiveRestartIndex", primitiveRestartIndex_);
}

template <class I, class V>
void GlMeshGeometry<I, V>::deserialize(Deserializer& d) {
    GlMeshGeometryBase::deserialize(d);
    d.deserializeBinaryBlob("vertices", vertices_);
    d.deserializeBinaryBlob("indices", indices_);

    d.deserialize("primitiveRestartEnabled", primitiveRestartEnabled_);
    d.deserialize("primitiveRestartIndex", primitiveRestartIndex_);
}

template<class I, class V>
void GlMeshGeometry<I,V>::setSphereGeometry(float radius, tgt::vec3 center, tgt::vec4 color, size_t numSlicesAndTiles) {

    clear();

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

    for (size_t i = 0; i <= numSlicesAndTiles; ++i) {
        for (size_t j = 0; j <= numSlicesAndTiles; ++j) {

            size_t k = j % numSlicesAndTiles;
            tgt::vec3 normal = tgt::vec3(sn_fi[i] * cs_th[k], sn_fi[i] * sn_th[k], cs_fi[i]);
            tgt::vec3 position = radius * normal;
            tgt::vec2 uv = tgt::vec2(u, v);

            addVertex(position, normal, color, uv);

            u += du;
        }

        u = 0;
        v += dv;
    }

    //shift vertices by center (TODO: maybe better do this using a transformation?)
    for(auto& vertex : vertices_)
        vertex.pos_ += center;

    // create the indices for the triangles
    I offset  = 0;
    I offsetDelta = static_cast<I>(numSlicesAndTiles) + 1;
    for (size_t j = 0; j < numSlicesAndTiles; ++j) {
        for (size_t i = 1; i <= numSlicesAndTiles; ++i) {
            I idx = offset + static_cast<I>(i);
            addIndex(idx - 1);
            addIndex(idx + static_cast<I>(numSlicesAndTiles));
            addIndex(idx);

            addIndex(idx);
            addIndex(idx + static_cast<I>(numSlicesAndTiles));
            addIndex(idx + static_cast<I>(numSlicesAndTiles) + 1);
        }
        offset += offsetDelta;
    }

    setPrimitiveType(GL_TRIANGLE_STRIP);
    disablePrimitiveRestart();

    invalidate();
}


template<class I, class V>
void GlMeshGeometry<I,V>::setCylinderGeometry(tgt::vec4 color, float lowerRadius, float upperRadius, float height, size_t slices, size_t tiles, bool lowerCap, bool upperCap) {

    clear();

    //check parameter
    tgtAssert(height > 0, "Height must be greater zero!");
    if(height <= 0.f) {
        height = 1.f;
        LWARNING("Height must be greater zero. Set to 1.f.");
    }

    tgtAssert(slices > 2, "Number of slices must be > 2");
    if (slices < 3) {
        slices = 3;
        LWARNING("Number of slices must be > 2, setting to 3");
    }

    tgtAssert(tiles > 0, "Number of tiles must be > 0");
    if (tiles == 0) {
        tiles = 1;
        LWARNING("Number of tiles must be > 0, setting to 1");
    }

    // create the vertices
    const float dth = (2.0f * tgt::PIf) / slices;
    const float dh = height / tiles;

    const float du = 1.0f / slices;
    const float dv = 1.0f / tiles;
    const float dt = 1.0f;
    const float dr = dv * (lowerRadius - upperRadius);

    std::vector<float> cs_th, sn_th;
    float th = 0.0f;
    float u  = 0.f;
    float v  = 0.f;
    float h  = 0.f;

    for (size_t i = 0; i <= slices; ++i) {
        cs_th.push_back(std::cos(th));
        sn_th.push_back(std::sin(th));
        th += dth;
    }

    if (upperCap) {
        // center
        for (size_t j = 0; j < slices + 1; ++j) {
            tgt::vec3 position = tgt::vec3(0.f, 0.f, height);
            tgt::vec3 normal = tgt::vec3(0.f, 0.f, 1.f);
            // TODO: use texture indices for the three areas
            tgt::vec2 uv= tgt::vec2(0.5f, 0.5f/*, CYLINDER_TOP*/);

            addVertex(position, normal, color, uv);
        }

        // outer ring
        for (size_t j = 0; j <= slices; ++j) {
            size_t k = j % slices;
            tgt::vec2 n = tgt::vec2(cs_th[k], sn_th[k]);
            tgt::vec3 normal = tgt::vec3(0.f, 0.f, 1.f);
            tgt::vec3 position = tgt::vec3(upperRadius * n, height);
            // TODO: use texture indices for the three areas
            tgt::vec2 uv = tgt::vec2(0.5f + 0.5f * n.x, 0.5f - 0.5f * n.y /*, CYLINDER_TOP*/);

            addVertex(position, normal, color, uv);
            u += du;
        }
        u = 0;
    }

    for (size_t i = 0; i <= tiles; ++i) {
        float r = (1.0f - v) * upperRadius + v * lowerRadius;

        for (size_t j = 0; j <= slices; ++j) {
            size_t k = j % slices;
            tgt::vec2 n = tgt::vec2(cs_th[k], sn_th[k]);
            tgt::vec3 t1 = tgt::normalize(tgt::vec3(-dr*n, dh));
            tgt::vec3 t2 = tgt::normalize(tgt::cross(t1, tgt::vec3(n,0.f)));

            tgt::vec3 normal = tgt::normalize(tgt::cross(t2, t1));
            tgt::vec3 position = tgt::vec3(r * n, height - h);
            // TODO: use texture indices for the three areas
            tgt::vec2 uv = tgt::vec2(u, v /*, CYLINDER_MANTLE*/);

            addVertex(position, normal, color, uv);

            u += du;
        }

        u = 0;
        h += dh;
        v += dv;
    }

    if (lowerCap) {

        // outer ring
        for (size_t j = 0; j <= slices; ++j) {
            size_t k = j % slices;
            tgt::vec2 n = tgt::vec2(cs_th[k], sn_th[k]);
            tgt::vec3 normal = tgt::vec3(0.f, 0.f, -1.f);
            tgt::vec3 position = tgt::vec3(lowerRadius * n, 0.f);
            // TODO: use texture indices for the three areas
            tgt::vec2 uv = tgt::vec2(0.5f) + 0.5f * n; // /*, CYLINDER_BOTTOM)*/;

            addVertex(position, normal, color, uv);
            u += du;
        }

        // center
        for (size_t j = 0; j < slices + 1; ++j) {
            tgt::vec3 position = tgt::vec3::zero;
            tgt::vec3 normal = tgt::vec3(0.f, 0.f, -1.f);
            // TODO: use texture indices for the three areas
            tgt::vec2 uv= tgt::vec2(0.5f, 0.5f/*, CYLINDER_BOTTOM*/);

            addVertex(position, normal, color, uv);
        }
    }

    // create the index buffer
    size_t offset = 0;

    if (upperCap) {
        createCylinderIndices(2, slices + 1, offset, indices_);
        offset += 2 * (slices + 1);
    }

    createCylinderIndices(tiles, slices + 1, offset, indices_);
    offset += tiles * (slices + 1);

    if (lowerCap) {
        createCylinderIndices(2, slices + 1, offset, indices_);
        offset += 2 * (slices + 1);
    }

    setPrimitiveType(GL_TRIANGLES);
    disablePrimitiveRestart();
    invalidate();
}

template<class I, class V>
void GlMeshGeometry<I,V>::setDiskGeometry(float innerRadius, float outerRadius, size_t numSlices) {
    tgtAssert(innerRadius >= 0 && outerRadius > 0 && innerRadius < outerRadius, "Invalid radius or radii");

    tgt::vec3 normal(0, 0, 1);
    tgt::vec4 color = tgt::vec4::one;

    float relativeInnerRadius = innerRadius/outerRadius;

    if(innerRadius == 0) {
        setPrimitiveType(GL_TRIANGLE_FAN);
    } else {
        setPrimitiveType(GL_TRIANGLE_STRIP);
    }

    clear();
    for(size_t i=0; i <= numSlices; ++i) {
        float s = sin(tgt::PIf*2*i/numSlices);
        float c = cos(tgt::PIf*2*i/numSlices);
        tgt::vec3 vertPos(s, c, 0);
        tgt::vec2 texPos = tgt::vec2(s, c)*0.5f + tgt::vec2(0.5f);
        addVertex(outerRadius*vertPos, normal, color, texPos);
        if(innerRadius > 0) {
            addVertex(innerRadius*vertPos, normal, color, relativeInnerRadius*texPos);
        }
    }
}

template<class I, class V>
void GlMeshGeometry<I,V>::setRectangleGeometry(float width, float height) {
    tgt::vec3 normal(0, 0, 1);
    tgt::vec4 color = tgt::vec4::one;

    float half_width = width/2;
    float half_height = height/2;

    setPrimitiveType(GL_TRIANGLE_FAN);

    clear();
    addVertex(tgt::vec3(-half_width,-half_height,0), normal, color, tgt::vec2(0.f, 0.f));
    addVertex(tgt::vec3(-half_width, half_height,0), normal, color, tgt::vec2(0.f, 1.f));
    addVertex(tgt::vec3( half_width, half_height,0), normal, color, tgt::vec2(1.f, 1.f));
    addVertex(tgt::vec3( half_width,-half_height,0), normal, color, tgt::vec2(1.f, 0.f));
}

template<class I, class V>
void GlMeshGeometry<I,V>::setTriangleGeometry(float edgeLength) {
    tgt::vec3 normal(0, 0, 1);
    tgt::vec4 color = tgt::vec4::one;

    float y_offset = -edgeLength*sqrt(3.0)/6;

    setPrimitiveType(GL_TRIANGLES);

    clear();
    addVertex(tgt::vec3(edgeLength/2,y_offset,0), normal, color, tgt::vec2(1.f, 0.f));
    addVertex(tgt::vec3(0,edgeLength*sqrt(3.0)/2 + y_offset,0), normal, color, tgt::vec2(0.5f, 1.f));
    addVertex(tgt::vec3(-edgeLength/2,y_offset,0), normal, color, tgt::vec2(0.f, 0.f));
}

template<class I, class V>
void GlMeshGeometry<I,V>::setTetraHedronGeometry(float edgeLength) {
    tgt::vec3 normal(0, 0, 1);
    tgt::vec4 color = tgt::vec4::one;

    float z_offset = -edgeLength/sqrt(24.0);
    float y_offset = -edgeLength*sqrt(3.0)/6;

    setPrimitiveType(GL_TRIANGLES);

    // Prepare vertices
    tgt::vec3 v0( edgeLength/2,                          y_offset,   z_offset); // /=======================
    tgt::vec3 v1(            0, edgeLength*sqrt(3.0)/2 + y_offset,   z_offset); // bottom triangle vertices
    tgt::vec3 v2(-edgeLength/2,                          y_offset,   z_offset); // \=======================
    tgt::vec3 v3(            0,                                 0, 1+z_offset); // top vertex

    // Prepare normals of all 4 sides
    tgt::vec3 n012 = tgt::normalize(tgt::cross(v2-v0, v1-v0));
    tgt::vec3 n031 = tgt::normalize(tgt::cross(v1-v0, v3-v0));
    tgt::vec3 n023 = tgt::normalize(tgt::cross(v3-v0, v2-v0));
    tgt::vec3 n132 = tgt::normalize(tgt::cross(v2-v1, v3-v1));

    // Prepare texture coordinates
    tgt::vec2 t0(1.0f, 0.0f);
    tgt::vec2 t1(0.5f, 1.0f);
    tgt::vec2 t2(0.0f, 0.0f);

    clear();

    //Assuming culling of anti-clockwise triangles

    addVertex(v0, n012, color, t0);
    addVertex(v1, n012, color, t1);
    addVertex(v2, n012, color, t2);

    addVertex(v0, n031, color, t0);
    addVertex(v3, n031, color, t1);
    addVertex(v1, n031, color, t2);

    addVertex(v0, n023, color, t0);
    addVertex(v2, n023, color, t1);
    addVertex(v3, n023, color, t2);

    addVertex(v1, n132, color, t0);
    addVertex(v3, n132, color, t1);
    addVertex(v2, n132, color, t2);
}

template<class I, class V>
void GlMeshGeometry<I,V>::setCuboidGeometry(float width, float height, float depth) {
    tgt::vec4 color = tgt::vec4::one;

    setPrimitiveType(GL_TRIANGLES);

    float half_width = width/2;
    float half_height = height/2;
    float half_depth = depth/2;

    clear();

    //Assuming culling of anti-clockwise triangles

    //Bottom
    addVertex(tgt::vec3( half_width,-half_height,-half_depth), tgt::vec3( 0, 0,-1), color, tgt::vec2(1.f, 0.f));
    addVertex(tgt::vec3( half_width, half_height,-half_depth), tgt::vec3( 0, 0,-1), color, tgt::vec2(1.f, 1.f));
    addVertex(tgt::vec3(-half_width, half_height,-half_depth), tgt::vec3( 0, 0,-1), color, tgt::vec2(0.f, 1.f));
    addVertex(tgt::vec3(-half_width,-half_height,-half_depth), tgt::vec3( 0, 0,-1), color, tgt::vec2(0.f, 0.f));

    addIndex(0);
    addIndex(1);
    addIndex(2);

    addIndex(0);
    addIndex(2);
    addIndex(3);

    //Top
    addVertex(tgt::vec3(-half_width,-half_height, half_depth), tgt::vec3( 0, 0, 1), color, tgt::vec2(0.f, 0.f));
    addVertex(tgt::vec3(-half_width, half_height, half_depth), tgt::vec3( 0, 0, 1), color, tgt::vec2(0.f, 1.f));
    addVertex(tgt::vec3( half_width, half_height, half_depth), tgt::vec3( 0, 0, 1), color, tgt::vec2(1.f, 1.f));
    addVertex(tgt::vec3( half_width,-half_height, half_depth), tgt::vec3( 0, 0, 1), color, tgt::vec2(1.f, 0.f));

    addIndex(4);
    addIndex(5);
    addIndex(6);

    addIndex(4);
    addIndex(6);
    addIndex(7);

    //Front
    addVertex(tgt::vec3(-half_width,-half_height,-half_depth), tgt::vec3( 0,-1, 0), color, tgt::vec2(0.f, 0.f));
    addVertex(tgt::vec3(-half_width,-half_height, half_depth), tgt::vec3( 0,-1, 0), color, tgt::vec2(0.f, 1.f));
    addVertex(tgt::vec3( half_width,-half_height, half_depth), tgt::vec3( 0,-1, 0), color, tgt::vec2(1.f, 1.f));
    addVertex(tgt::vec3( half_width,-half_height,-half_depth), tgt::vec3( 0,-1, 0), color, tgt::vec2(1.f, 0.f));

    addIndex(8);
    addIndex(9);
    addIndex(10);

    addIndex(8);
    addIndex(10);
    addIndex(11);

    //Back
    addVertex(tgt::vec3( half_width, half_height,-half_depth), tgt::vec3( 0, 1, 0), color, tgt::vec2(1.f, 0.f));
    addVertex(tgt::vec3( half_width, half_height, half_depth), tgt::vec3( 0, 1, 0), color, tgt::vec2(1.f, 1.f));
    addVertex(tgt::vec3(-half_width, half_height, half_depth), tgt::vec3( 0, 1, 0), color, tgt::vec2(0.f, 1.f));
    addVertex(tgt::vec3(-half_width, half_height,-half_depth), tgt::vec3( 0, 1, 0), color, tgt::vec2(0.f, 0.f));

    addIndex(12);
    addIndex(13);
    addIndex(14);

    addIndex(12);
    addIndex(14);
    addIndex(15);

    //Left
    addVertex(tgt::vec3(-half_width,-half_height,-half_depth), tgt::vec3(-1, 0, 0), color, tgt::vec2(0.f, 0.f));
    addVertex(tgt::vec3(-half_width, half_height,-half_depth), tgt::vec3(-1, 0, 0), color, tgt::vec2(0.f, 1.f));
    addVertex(tgt::vec3(-half_width, half_height, half_depth), tgt::vec3(-1, 0, 0), color, tgt::vec2(1.f, 1.f));
    addVertex(tgt::vec3(-half_width,-half_height, half_depth), tgt::vec3(-1, 0, 0), color, tgt::vec2(1.f, 0.f));

    addIndex(16);
    addIndex(17);
    addIndex(18);

    addIndex(16);
    addIndex(18);
    addIndex(19);

    //Right
    addVertex(tgt::vec3( half_width,-half_height, half_depth), tgt::vec3( 1, 0, 0), color, tgt::vec2(1.f, 0.f));
    addVertex(tgt::vec3( half_width, half_height, half_depth), tgt::vec3( 1, 0, 0), color, tgt::vec2(1.f, 1.f));
    addVertex(tgt::vec3( half_width, half_height,-half_depth), tgt::vec3( 1, 0, 0), color, tgt::vec2(0.f, 1.f));
    addVertex(tgt::vec3( half_width,-half_height,-half_depth), tgt::vec3( 1, 0, 0), color, tgt::vec2(0.f, 0.f));

    addIndex(20);
    addIndex(21);
    addIndex(22);

    addIndex(20);
    addIndex(22);
    addIndex(23);
}


template<class I, class V>
void GlMeshGeometry<I,V>::createCylinderIndices(size_t tiles, size_t slices, size_t offset, std::vector<I>& indices) {
    for (size_t j = 0; j < tiles; ++j) {
        for (size_t i = 1; i < slices; ++i) {
            I idx = static_cast<I>(offset + i);
            addIndex(static_cast<I>(idx - 1));
            addIndex(static_cast<I>(idx + slices - 1));
            addIndex(static_cast<I>(idx));

            addIndex(static_cast<I>(idx));
            addIndex(static_cast<I>(idx + slices - 1));
            addIndex(static_cast<I>(idx + slices));
        }
        offset += slices;
    }
}


/////////////////// OLD CODE ////////////////////////////////////////////////////////////////

// old mesh methods
/*
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
        addTriangle(TriangleType(I(startVertex + i), I(startVertex + i + beginOfInner), I(startVertex + (i + 1 == beginOfInner ? 0 : i + 1))));
        addTriangle(TriangleType(I(startVertex + i + beginOfInner), I(startVertex + (i + 1 == beginOfInner ? beginOfInner : i + beginOfInner + 1)),
                                 I(startVertex + (i + 1 == beginOfInner ? 0 : i + 1))));
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
                                 I(startVertex + (i + 1 == beginOfTop ? 0 : i + 1))));
        addTriangle(TriangleType(I(startVertex + i + beginOfTop),
                                 I(startVertex + (i + 1 == beginOfTop ? beginOfTop : i + beginOfTop + 1)),
                                 I(startVertex + (i + 1 == beginOfTop ? 0 : i + 1))));
    }
    //invalidate geometry
    invalidate();
}

template<class I, class V>
void TriangleMeshGeometryIndexed<I, V>::addSphereGeometry(float radius, tgt::vec3 center, tgt::vec4 color, SimpleGeometryQuality quality) {
    addSphereGeometry(radius, center, color, static_cast<size_t>(quality + 3));
}



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
}*/


/////////////////// ACTUAL CLASSES TO BE USED ///////////////////////////////////////////////
//-------------------------------------------------------------------------------------------

/// A mesh with vertex type VertexBase with indices of type uint16_t.
class VRN_CORE_API GlMeshGeometryUInt16Simple : public GlMeshGeometry<uint16_t, VertexBase> {
public:
    virtual Geometry* create() const { return new GlMeshGeometryUInt16Simple(); }
    virtual std::string getClassName() const { return "GlMeshGeometryUInt16Simple"; }
    virtual IndexType getIndexType() const {return INDEX_UINT16;}
    virtual GLenum getIndexDataType() const {return GL_UNSIGNED_SHORT;}
};

/// A mesh with vertex type VertexColor with indices of type uint16_t.
class VRN_CORE_API GlMeshGeometryUInt16Color : public GlMeshGeometry<uint16_t, VertexColor> {
public:
    virtual Geometry* create() const { return new GlMeshGeometryUInt16Color(); }
    virtual std::string getClassName() const { return "GlMeshGeometryUInt16Color"; }
    virtual IndexType getIndexType() const {return INDEX_UINT16;}
    virtual GLenum getIndexDataType() const {return GL_UNSIGNED_SHORT;}
};

/// A mesh with vertex type VertexNormal with indices of type uint16_t.
class VRN_CORE_API GlMeshGeometryUInt16Normal : public GlMeshGeometry<uint16_t, VertexNormal> {
public:
    virtual Geometry* create() const { return new GlMeshGeometryUInt16Normal(); }
    virtual std::string getClassName() const { return "GlMeshGeometryUInt16Normal"; }
    virtual IndexType getIndexType() const {return INDEX_UINT16;}
    virtual GLenum getIndexDataType() const {return GL_UNSIGNED_SHORT;}
};

/// A mesh with vertex type VertexTexCoord with indices of type uint16_t.
class VRN_CORE_API GlMeshGeometryUInt16TexCoord : public GlMeshGeometry<uint16_t, VertexTexCoord> {
public:
    virtual Geometry* create() const { return new GlMeshGeometryUInt16TexCoord(); }
    virtual std::string getClassName() const { return "GlMeshGeometryUInt16TexCoord"; }
    virtual IndexType getIndexType() const {return INDEX_UINT16;}
    virtual GLenum getIndexDataType() const {return GL_UNSIGNED_SHORT;}
};

/// A mesh with vertex type VertexColorNormal with indices of type uint16_t.
class VRN_CORE_API GlMeshGeometryUInt16ColorNormal : public GlMeshGeometry<uint16_t, VertexColorNormal> {
public:
    virtual Geometry* create() const { return new GlMeshGeometryUInt16ColorNormal(); }
    virtual std::string getClassName() const { return "GlMeshGeometryUInt16ColorNormal"; }
    virtual IndexType getIndexType() const {return INDEX_UINT16;}
    virtual GLenum getIndexDataType() const {return GL_UNSIGNED_SHORT;}
};

/// A mesh with vertex type VertexTextureNormal with indices of type uint16_t.
class VRN_CORE_API GlMeshGeometryUInt16NormalTexCoord : public GlMeshGeometry<uint16_t, VertexNormalTexCoord> {
public:
    virtual Geometry* create() const { return new GlMeshGeometryUInt16NormalTexCoord(); }
    virtual std::string getClassName() const { return "GlMeshGeometryUInt16NormalTexCoord"; }
    virtual IndexType getIndexType() const {return INDEX_UINT16;}
    virtual GLenum getIndexDataType() const {return GL_UNSIGNED_SHORT;}
};

/// A mesh with vertex type VertexNormal with indices of type uint16_t.
class VRN_CORE_API GlMeshGeometryUInt16ColorTexCoord : public GlMeshGeometry<uint16_t, VertexColorTexCoord> {
public:
    virtual Geometry* create() const { return new GlMeshGeometryUInt16ColorTexCoord(); }
    virtual std::string getClassName() const { return "GlMeshGeometryUInt16ColorTexCoord"; }
    virtual IndexType getIndexType() const {return INDEX_UINT16;}
    virtual GLenum getIndexDataType() const {return GL_UNSIGNED_SHORT;}
};

/// A mesh with vertex type VertexColorNormalTexCoord with indices of type uint16_t.
class VRN_CORE_API GlMeshGeometryUInt16ColorNormalTexCoord : public GlMeshGeometry<uint16_t, VertexColorNormalTexCoord> {
public:
    virtual Geometry* create() const { return new GlMeshGeometryUInt16ColorNormalTexCoord(); }
    virtual std::string getClassName() const { return "GlMeshGeometryUInt16ColorNormalTexCoord"; }
    virtual IndexType getIndexType() const {return INDEX_UINT16;}
    virtual GLenum getIndexDataType() const {return GL_UNSIGNED_SHORT;}
};

//-------------------------------------------------------------------------------------------

/// A mesh with vertex type VertexBase with indices of type uint32_t.
class VRN_CORE_API GlMeshGeometryUInt32Simple : public GlMeshGeometry<uint32_t, VertexBase> {
public:
    virtual Geometry* create() const { return new GlMeshGeometryUInt32Simple(); }
    virtual std::string getClassName() const { return "GlMeshGeometryUInt32Simple"; }
    virtual IndexType getIndexType() const {return INDEX_UINT32;}
    virtual GLenum getIndexDataType() const {return GL_UNSIGNED_INT;}
};

/// A mesh with vertex type VertexColor with indices of type uint32_t.
class VRN_CORE_API GlMeshGeometryUInt32Color : public GlMeshGeometry<uint32_t, VertexColor> {
public:
    virtual Geometry* create() const { return new GlMeshGeometryUInt32Color(); }
    virtual std::string getClassName() const { return "GlMeshGeometryUInt32Color"; }
    virtual IndexType getIndexType() const {return INDEX_UINT32;}
    virtual GLenum getIndexDataType() const {return GL_UNSIGNED_INT;}
};

/// A mesh with vertex type VertexNormal with indices of type uint32_t.
class VRN_CORE_API GlMeshGeometryUInt32Normal : public GlMeshGeometry<uint32_t, VertexNormal> {
public:
    virtual Geometry* create() const { return new GlMeshGeometryUInt32Normal(); }
    virtual std::string getClassName() const { return "GlMeshGeometryUInt32Normal"; }
    virtual IndexType getIndexType() const {return INDEX_UINT32;}
    virtual GLenum getIndexDataType() const {return GL_UNSIGNED_INT;}
};

/// A mesh with vertex type VertexTexCoord with indices of type uint32_t.
class VRN_CORE_API GlMeshGeometryUInt32TexCoord : public GlMeshGeometry<uint32_t, VertexTexCoord> {
public:
    virtual Geometry* create() const { return new GlMeshGeometryUInt32TexCoord(); }
    virtual std::string getClassName() const { return "GlMeshGeometryUInt32TexCoord"; }
    virtual IndexType getIndexType() const {return INDEX_UINT32;}
    virtual GLenum getIndexDataType() const {return GL_UNSIGNED_INT;}
};

/// A mesh with vertex type VertexColorNormal with indices of type uint32_t.
class VRN_CORE_API GlMeshGeometryUInt32ColorNormal : public GlMeshGeometry<uint32_t, VertexColorNormal> {
public:
    virtual Geometry* create() const { return new GlMeshGeometryUInt32ColorNormal(); }
    virtual std::string getClassName() const { return "GlMeshGeometryUInt32ColorNormal"; }
    virtual IndexType getIndexType() const {return INDEX_UINT32;}
    virtual GLenum getIndexDataType() const {return GL_UNSIGNED_INT;}
};

/// A mesh with vertex type VertexTextureNormal with indices of type uint32_t.
class VRN_CORE_API GlMeshGeometryUInt32NormalTexCoord : public GlMeshGeometry<uint32_t, VertexNormalTexCoord> {
public:
    virtual Geometry* create() const { return new GlMeshGeometryUInt32NormalTexCoord(); }
    virtual std::string getClassName() const { return "GlMeshGeometryUInt32NormalTexCoord"; }
    virtual IndexType getIndexType() const {return INDEX_UINT32;}
    virtual GLenum getIndexDataType() const {return GL_UNSIGNED_INT;}
};

/// A mesh with vertex type VertexNormal with indices of type uint32_t.
class VRN_CORE_API GlMeshGeometryUInt32ColorTexCoord : public GlMeshGeometry<uint32_t, VertexColorTexCoord> {
public:
    virtual Geometry* create() const { return new GlMeshGeometryUInt32ColorTexCoord(); }
    virtual std::string getClassName() const { return "GlMeshGeometryUInt32ColorTexCoord"; }
    virtual IndexType getIndexType() const {return INDEX_UINT32;}
    virtual GLenum getIndexDataType() const {return GL_UNSIGNED_INT;}
};

/// A mesh with vertex type VertexColorNormalTexCoord with indices of type uint32_t.
class VRN_CORE_API GlMeshGeometryUInt32ColorNormalTexCoord : public GlMeshGeometry<uint32_t, VertexColorNormalTexCoord> {
public:
    virtual Geometry* create() const { return new GlMeshGeometryUInt32ColorNormalTexCoord(); }
    virtual std::string getClassName() const { return "GlMeshGeometryUInt32ColorNormalTexCoord"; }
    virtual IndexType getIndexType() const {return INDEX_UINT32;}
    virtual GLenum getIndexDataType() const {return GL_UNSIGNED_INT;}
};


} // namespace

#endif  //VRN_GLMESHGEOMETRY_H
