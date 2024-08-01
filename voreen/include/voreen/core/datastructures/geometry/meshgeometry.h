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

#ifndef VRN_MESHGEOMETRY_H
#define VRN_MESHGEOMETRY_H

#include <vector>

#include "tgt/plane.h"

#include "voreen/core/datastructures/geometry/geometry.h"
#include "voreen/core/datastructures/geometry/facegeometry.h"

namespace voreen {

/**
 * Represents a mesh geometry containing a list of face geometries.
 *
 * @note The following assumption has to be fulfilled:
 *       (1) the face list represents a closed convex polyhedron.
 *
 * @par
 * Here is a short example of using the @c MeshGeometry:
 * @code
 * MeshGeometry cube = MeshGeometry::createCube();
 * cube.clip(tgt::vec4(0, 1, 0, 0));
 * cube.render();
 * @endcode
 *
 * @see FaceGeometry
 * @see VertexGeometry
 * @see MeshListGeometry
 */
class VRN_CORE_API MeshGeometry : public Geometry {
public:
    /**
     * Type of the face geometry list.
     */
    typedef std::vector<FaceGeometry> FaceListType;

    /**
     * Type of the face geometry list iterator.
     */
    typedef FaceListType::iterator iterator;

    /**
     * Type of the const face geometry list iterator.
     */
    typedef FaceListType::const_iterator const_iterator;

    /**
     * Default constructor.
     */
    MeshGeometry();

    virtual ~MeshGeometry();

    virtual Geometry* create() const;

    virtual std::string getClassName() const { return "MeshGeometry"; }

    /**
     * Creates a cube mesh with given dimension and properties.
     *
     * @param coordLlf the vertex coordinates of the lower left front vertex
     * @param coordUrb the vertex coordinates of the upper right back vertex
     * @param texLlf the texture coordinates of the lower left front vertex
     * @param texUrb the texture coordinates of the upper right back vertex
     * @param colorLlf the vertex color of the lower left front vertex
     * @param colorUrb the vertex color of the upper right back vertex
     * @param alpha the alpha value used in the color
     *
     * @returns the cube mesh
     */
    static MeshGeometry createCube(
        tgt::vec3 coordLlf = tgt::vec3(0.f, 0.f, 0.f),
        tgt::vec3 coordUrb = tgt::vec3(1.f, 1.f, 1.f),
        tgt::vec3 texLlf = tgt::vec3(0.f, 0.f, 0.f),
        tgt::vec3 texUrb = tgt::vec3(1.f, 1.f, 1.f),
        tgt::vec3 colorLlf = tgt::vec3(0.f, 0.f, 0.f),
        tgt::vec3 colorUrb = tgt::vec3(0.f, 0.f, 0.f),
        float alpha = 1.f
        );

    /**
     * Creates a cube mesh with given dimension and properties.
     *
     * @param coordLlf the vertex coordinates of the lower left front vertex
     * @param coordUrb the vertex coordinates of the upper right back vertex
     * @param texLlf the texture coordinates of the lower left front vertex
     * @param texUrb the texture coordinates of the upper right back vertex
     * @param colorLlf the vertex color of the lower left front vertex
     * @param colorUrb the vertex color of the upper right back vertex
     * @param normalTop the normal direction of the top face
     * @param normalFront the normal direction of the front face
     * @param normalLeft the normal direction of the left face
     * @param normalBack the normal direction of the back face
     * @param normalRight the normal direction of the right face
     * @param normalBottom the normal direction of the bottom face
     * @param alpha the alpha value used in the color
     *
     * @returns the cube mesh
     */
    static MeshGeometry createCube(
        tgt::vec3 coordLlf,
        tgt::vec3 coordUrb,
        tgt::vec3 texLlf,
        tgt::vec3 texUrb,
        tgt::vec3 colorLlf,
        tgt::vec3 colorUrb,
        tgt::vec3 normalTop,
        tgt::vec3 normalFront,
        tgt::vec3 normalLeft,
        tgt::vec3 normalBack,
        tgt::vec3 normalRight,
        tgt::vec3 normalBottom,
        float alpha = 1.f
        );
    /**
     * Returns the number of face geometries contained by this mesh geometry.
     *
     * @returns face geometry count
     */
    size_t getFaceCount() const;

    /**
     * Returns the number of vertices contained by this mesh geometry,
     * i.e. the sum of its faces' vertex counts.
     *
     * @returns face vertex count
     */
    size_t getVertexCount() const;

    /**
     * Returns true, if the mesh geometry's vertex count is zero.
     *
     * @see getVertexCount
     */
    bool empty() const;

    /**
     * Adds the given face geometry to this mesh geometry.
     *
     * @param face the face geometry
     */
    void addFace(const FaceGeometry& face);

    /**
     * Returns a const reference to the face geometry at the given @c index.
     *
     * @param index the face geometry index
     *
     * @returns the face geometry at the given @c index
     */
    const FaceGeometry& getFace(size_t index) const;

    /**
     * Returns the face geometry at the given @c index.
     *
     * @param index the face geometry index
     *
     * @returns the face geometry at the given @c index
     */
    FaceGeometry& getFace(size_t index);

    /**
     * Removes all face geometries from this mesh geometry.
     */
    void clear();

    /**
     * Returns an iterator pointing to the first face geometry of this mesh geometry.
     *
     * @returns an iterator pointing to the first face geometry
     */
    iterator begin();

    /**
     * Returns a const iterator pointing to the first face geometry of this mesh geometry.
     *
     * @returns a const iterator pointing to the first face geometry
     */
    const_iterator begin() const;

    /**
     * Returns an iterator pointing behind the last face geometry of this face geometry.
     *
     * @returns an iterator pointing behind the last face geometry
     */
    iterator end();

    /**
     * Returns a const iterator pointing behind the last face geometry of this face geometry.
     *
     * @returns a const iterator pointing behind the last face geometry
     */
    const_iterator end() const;

    /**
     * Returns the face geometry at the given @c index.
     *
     * @param index the face geometry index
     *
     * @returns the face geometry at the given @c index
     */
    FaceGeometry& operator[](size_t index);

    /**
     * @see Geometry::render
     */
    virtual void render() const;

    /**
     * Clips the mesh geometry by the given arbitrary clipping plane.
     *
     * @par
     * The arbitrary clipping plane is represented by the equation
     * NX - c = 0 where N is the normal vector of the clipping plane
     * and corresponds to the x-, y-, and z-coordinate of the given vector,
     * and c is a scalar corresponding to the w-coordinate of the given
     * vector and X is any point on the plane.
     *
     * @note It is not necessary that the clipping plane normal vector
     *       is normalized, but it may improve computational accuracy.
     *
     * @note Use the @c epsilon parameter to change the accuracy of
     *       vertex- / vertex geometry and vertex geometry / clipping plane comparison.
     *
     * @note The clipping algorithm is based on the paper "Clipping a Mesh Against a Plane"
     *       by David Eberly. For further information see:
     *       http://www.geometrictools.com/Documentation/ClipMesh.pdf
     *
     * @param clipPlane the arbitrary clipping plane
     * @param closingMesh the clipping face generated for closing the clipped mesh
     * @param epsilon the accuracy for vertex geometry comparison
     */
    void clip(const tgt::plane& clipPlane, MeshGeometry& closingMesh, double epsilon = 1e-5);

    /// @overload
    virtual void clip(const tgt::plane& clipPlane, double epsilon = 1e-5);

    /**
     * Returns true, if all faces of the passed MeshGeometry are equal to this one's.
     *
     * @param mesh the mesh to compare
     * @param epsilon maximum distance at which two vertices are to be considered equal
     */
    bool equals(const MeshGeometry& mesh, double epsilon = 1e-5) const;

    /**
     * Returns true, if the passed Geometry is a MeshGeometry
     * and all its faces are equal to this one's.
     *
     * @see Geometry::equals
     */
    virtual bool equals(const Geometry* geometry, double epsilon = 1e-5) const;

    /**
     * Returns the axis-aligned bounding box of the union
     * of the bounding boxes of the mesh's faces.
     */
    virtual tgt::Bounds getBoundingBox(bool transformed = true) const;

    virtual void serialize(Serializer& s) const;

    virtual void deserialize(Deserializer& s);

private:
    static void createCubeFaces(FaceGeometry& topFace, FaceGeometry& frontFace, FaceGeometry& leftFace,
                                FaceGeometry& backFace,FaceGeometry& rightFace, FaceGeometry& bottomFace,
                                tgt::vec3 coordLlf = tgt::vec3(0.f, 0.f, 0.f),
                                tgt::vec3 coordUrb = tgt::vec3(1.f, 1.f, 1.f),
                                tgt::vec3 texLlf = tgt::vec3(0.f, 0.f, 0.f),
                                tgt::vec3 texUrb = tgt::vec3(1.f, 1.f, 1.f),
                                tgt::vec3 colorLlf = tgt::vec3(0.f, 0.f, 0.f),
                                tgt::vec3 colorUrb = tgt::vec3(0.f, 0.f, 0.f),
                                float alpha = 1.f);
    /**
     * Face geometry list.
     */
    FaceListType faces_;
};

} // namespace

#endif  //VRN_MESHGEOMETRY_H
