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

#ifndef VRN_MESHLISTGEOMETRY_H
#define VRN_MESHLISTGEOMETRY_H

#include <vector>

#include "voreen/core/datastructures/geometry/geometry.h"
#include "voreen/core/datastructures/geometry/meshgeometry.h"

namespace voreen {

/**
 * Represents a mesh geometry containing a list of mesh geometries.
 *
 * @par
 * Here is a short example of using the @c FaceGeometry:
 * @code
 * MeshListGeometry meshes;
 * meshes.addMesh(MeshGeometry::createCube());
 * meshes.addMesh(MeshGeometry::createCube(tgt::vec3(1, -1, 1), tgt::vec3(3, 1, -1)));
 *
 * meshes.clip(tgt::vec4(0, 1, 0, 0));
 * meshes.render();
 * @endcode
 *
 * @see Geometry
 * @see MeshGeometry
 * @see FaceGeometry
 * @see VertexGeometry
 */
class VRN_CORE_API MeshListGeometry : public Geometry {
public:
    /**
     * Type of the mesh geometry list.
     */
    typedef std::vector<MeshGeometry> MeshListType;

    /**
     * Type of the mesh geometry list iterator.
     */
    typedef MeshListType::iterator iterator;

    /**
     * Type of the const mesh geometry list iterator.
     */
    typedef MeshListType::const_iterator const_iterator;

    /**
     * Default constructor.
     */
    MeshListGeometry();

    virtual ~MeshListGeometry();

    virtual Geometry* create() const;

    virtual std::string getClassName() const { return "MeshListGeometry"; }

    /**
     * Returns the number of mesh geometries contained by this mesh list geometry.
     *
     * @returns mesh geometry count
     */
    size_t getMeshCount() const;

    /**
     * Returns the number of mesh faces contained by this mesh list geometry,
     * i.e. the sum of its meshes' face counts.
     *
     * @see MeshGeometry::getFaceCount
     */
    size_t getFaceCount() const;

    /**
     * Returns the number of mesh vertices contained by this mesh list geometry,
     * i.e. the sum of its meshes' vertex counts.
     *
     * @see MeshGeometry::getVertexCount
     */
    size_t getVertexCount() const;

    /**
     * Returns true, if the geometry's vertex count is zero.
     *
     * @see getVertexCount
     */
    bool empty() const;

    /**
     * Adds the given mesh geometry to this mesh list geometry.
     *
     * @param mesh the mesh geometry
     */
    void addMesh(const MeshGeometry& mesh);

    /**
     * Returns a const reference to the mesh geometry at the given @c index.
     *
     * @param index the mesh geometry index
     *
     * @returns the mesh geometry at the given @c index
     */
    const MeshGeometry& getMesh(size_t index) const;

    /**
     * Returns a the mesh geometry at the given @c index.
     *
     * @param index the mesh geometry index
     *
     * @returns the mesh geometry at the given @c index
     */
    MeshGeometry& getMesh(size_t index);

    /**
     * Removes all mesh geometries from this mesh list geometry.
     */
    void clear();

    /**
     * Returns an iterator pointing to the first mesh geometry of this mesh list geometry.
     *
     * @returns an iterator pointing to the first mesh geometry
     */
    iterator begin();

    /**
     * Returns a const iterator pointing to the first mesh geometry of this mesh list geometry.
     *
     * @returns a const iterator pointing to the first mesh geometry
     */
    const_iterator begin() const;

    /**
     * Returns an iterator pointing behind the last mesh geometry of this mesh list geometry.
     *
     * @returns an iterator pointing behind the last mesh geometry
     */
    iterator end();

    /**
     * Returns a const iterator pointing behind the last mesh geometry of this mesh list geometry.
     *
     * @returns a const iterator pointing behind the last mesh geometry
     */
    const_iterator end() const;

    /**
     * Returns the mesh geometry at the given @c index.
     *
     * @param index the mesh geometry index
     *
     * @returns the mesh geometry at the given @c index
     */
    MeshGeometry& operator[](size_t index);
    const MeshGeometry& operator[](size_t index) const;

    /**
     * Returns the axis-aligned bounding box of the union
     * of the bounding boxes of the contained meshes.
     */
    virtual tgt::Bounds getBoundingBox(bool transformed = true) const;

    /**
     * @see Geometry::render
     */
    virtual void render() const;

    /**
     * Calls the @c clip function of each contained mesh geometry.
     *
     * @see MeshGeometry::clip
     *
     * @param clipPlane an arbitrary clipping plane
     * @param closingFaces the clipping faces generated for closing the clipped meshes
     * @param epsilon the accuracy for vertex geometry comparison
     */
    void clip(const tgt::plane& clipPlane, MeshListGeometry& closingFaces, double epsilon = 1e-5);

    /// @overload
    virtual void clip(const tgt::plane& clipPlane, double epsilon = 1e-5);

    /**
     * Returns true, if all meshes of the passed MeshListGeometry are equal to this one's.
     *
     * @param meshList the mesh list to compare
     * @param epsilon maximum distance at which two vertices are to be considered equal
     */
    bool equals(const MeshListGeometry& meshList, double epsilon = 1e-5) const;

    /**
     * Returns true, if the passed Geometry is a MeshListGeometry
     * and all its meshes are equal to this one's.
     *
     * @see Geometry::equals
     */
    virtual bool equals(const Geometry* geometry, double epsilon = 1e-5) const;

    virtual void serialize(Serializer& s) const;

    virtual void deserialize(Deserializer& s);

private:
    /**
     * Mesh geometry list.
     */
    MeshListType meshes_;
};

}    // namespace

#endif  //VRN_MESHLISTGEOMETRY_H
