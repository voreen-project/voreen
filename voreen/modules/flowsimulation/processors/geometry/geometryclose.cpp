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

#include "geometryclose.h"

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include "../../ext/octree/Octree.hpp"
#include "../../ext/halfedge/trimesh.h"

#define ACCESS_VEC_IDX(V, I) \
    template <> \
    struct access<V, I> { \
        static float get(const V& p) { \
            return p.pos_[I]; \
        } \
    };

#define ACCESS_VEC(V) \
    ACCESS_VEC_IDX(V, 0) \
    ACCESS_VEC_IDX(V, 1) \
    ACCESS_VEC_IDX(V, 2)


namespace unibn {
namespace traits {

ACCESS_VEC(voreen::GlMeshGeometryUInt32Simple::VertexType)
ACCESS_VEC(voreen::GlMeshGeometryUInt32Color::VertexType)
ACCESS_VEC(voreen::GlMeshGeometryUInt32Normal::VertexType)
ACCESS_VEC(voreen::GlMeshGeometryUInt32TexCoord::VertexType)
ACCESS_VEC(voreen::GlMeshGeometryUInt32ColorTexCoord::VertexType)
ACCESS_VEC(voreen::GlMeshGeometryUInt32NormalTexCoord::VertexType)
ACCESS_VEC(voreen::GlMeshGeometryUInt32ColorNormal::VertexType)
ACCESS_VEC(voreen::GlMeshGeometryUInt32ColorNormalTexCoord::VertexType)

}
}

#undef ACCESS_VEC
#undef ACCESS_VEC_IDX

namespace voreen {

/**
 * Creates indices for this mesh, if not already available.
 * This will enable indexed drawing.
 * The vertices stay untouched, such that non-indexed drawing yields the same result.
 *
 * The indices can be optimized such that duplicated vertices will have
 * the very same index pointing to them.
 *
 * @param geometry Geometry for that indices shall be created
 * @param epsilon threshold for neighborhood determination
 * @param optimize determines, whether indices should be optimized
 */
template<typename T>
void createIndices(T* geometry, float epsilon, bool optimize) {

    typedef uint32_t I;
    typedef typename T::VertexType V;

    std::vector<I> indices = geometry->getIndices();
    const std::vector<V>& vertices = geometry->getVertices();

    // Add trivial indices.
    if(indices.empty()) {
        indices.resize(vertices.size());
        std::iota(indices.begin(), indices.end(), 0);
    }

    if(optimize) {

        // initializing the Octree with points from point cloud.
        unibn::Octree<V> octree;
        unibn::OctreeParams params;
        octree.initialize(vertices, indices, params);

        std::unordered_set<I> seenAlready;
        for (size_t i = 0; i < indices.size(); i++) {

            if(seenAlready.find(i) != seenAlready.end()) {
                continue;
            }

            I& index = indices[i];
            const V& vertex = vertices[index];

            std::vector<I> neighbors;
            octree.template radiusNeighbors<unibn::L2Distance<V>>(vertex, epsilon, neighbors);
            seenAlready.insert(index);
            seenAlready.insert(neighbors.begin(), neighbors.end());

            // TODO: we assume that index positions are equal to their value here!
            for(I neighborIdx : neighbors) {
                indices[neighborIdx] = index;
            }
        }
    }

    geometry->setIndices(indices);
}

template<typename T>
size_t fill(T* geometry, float epsilon) {

    typedef uint32_t I;
    typedef typename T::VertexType V;

    // Ensure the indices are set correctly.
    bool usesIndices = geometry->getNumIndices() > 0;
    if(!usesIndices) {
        createIndices(geometry, epsilon, true);
    }

    // Build the mesh.
    trimesh::trimesh_t mesh;
    {
        size_t numIndices = geometry->getNumIndices();
        std::vector<trimesh::triangle_t> triangles(numIndices / 3);
        for (size_t i = 0; i < geometry->getNumIndices(); i += 3) {
            trimesh::triangle_t triangle;
            triangle.v[0] = geometry->getIndices()[i + 0];
            triangle.v[1] = geometry->getIndices()[i + 1];
            triangle.v[2] = geometry->getIndices()[i + 2];
            triangles[i / 3] = triangle;
        }

        // We no longer need the indices.
        if(!usesIndices) {
            geometry->setIndices(std::vector<I>());
        }

        std::vector<trimesh::edge_t> edges;
        trimesh::unordered_edges_from_triangles(triangles.size(), &triangles[0], edges);

        mesh.build(numIndices, triangles.size(), &triangles[0], edges.size(), &edges[0]);
    }

    // Close holes by performing an edge loop around the null-pointing faces.
    std::set<trimesh::index_t> boundary;
    {
        std::vector<trimesh::index_t> boundary_tmp = mesh.boundary_vertices();
        boundary.insert(boundary_tmp.begin(), boundary_tmp.end());
    }

    // Outer loop iterates holes.
    size_t numHoles = 0;
    auto start = boundary.begin();
    while (start != boundary.end()) {

        // Inner loop iterates vertices for each hole.
        std::vector<trimesh::index_t> ring;
        trimesh::index_t current = *start;

        while(true) {
            tgtAssert(boundary.find(current) != boundary.end(), "Not part of boundary");

            // Update ring.
            ring.push_back(current);
            boundary.erase(current);

            // Find the next appropriate halfedge.
            trimesh::index_t halfEdgeIndex = -1;
            std::vector<trimesh::index_t> neighbours = mesh.vertex_vertex_neighbors(current);
            for(trimesh::index_t neighbour : neighbours) {
                if(boundary.find(neighbour) != boundary.end()) {
                    trimesh::index_t candidateHalfEdgeIndex = mesh.directed_edge2he_index(current, neighbour);
                    if(mesh.halfedge(candidateHalfEdgeIndex).face == -1) {
                        tgtAssert(halfEdgeIndex == -1, "more than one candidate found");
                        halfEdgeIndex = candidateHalfEdgeIndex;
#ifndef VRN_DEBUG
                        break;
#endif
                    }
                }
            }

            // TODO: Having found no halfedge could also be caused by an incorrect data structure
            //tgtAssert(halfEdgeIndex >= 0, "No halfedge found");
            if(halfEdgeIndex == -1)
                break;

            // Update current.
            current = mesh.halfedge(halfEdgeIndex).to_vertex;
        }

        size_t centerIndex = ring[0]; // Arbitrary choice.
        V center = geometry->getVertices()[centerIndex];
        if(!usesIndices) {
            for (size_t i = 1; i < ring.size() - 1; i++) {
                geometry->addVertex(geometry->getVertex(ring[i]));
                geometry->addVertex(geometry->getVertex(ring[(i + 1) % ring.size()]));
                geometry->addVertex(center);
            }
        }
        else {
            for (size_t i = 1; i < ring.size() - 1; i++) {
                geometry->addIndex(ring[i]);
                geometry->addIndex(ring[(i + 1) % ring.size()]);
                geometry->addIndex(centerIndex);
            }
        }

        start = boundary.begin();
        numHoles++;
    }
    tgtAssert(boundary.empty(), "Some boundary vertices unprocessed");

    return numHoles;
}


const std::string GeometryClose::loggerCat_("voreen.flowsimulation.GeometryClose");

GeometryClose::GeometryClose()
    : Processor()
    , inport_(Port::INPORT, "geometry.input", "Geometry Input")
    , outport_(Port::OUTPORT, "geometry.output", "Geometry Output")
    , enabled_("enabled", "Enable", true)
    , epsilon_("epsilon", "Epsilon", 1e-4f, 1e-7f, 1e-3f, Processor::INVALID_RESULT,FloatProperty::STATIC, Property::LOD_DEBUG)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(enabled_);
    addProperty(epsilon_);
    epsilon_.setTracking(false);
    epsilon_.adaptDecimalsToRange(7);
}

GeometryClose::~GeometryClose()
{}

Processor* GeometryClose::create() const {
    return new GeometryClose();
}

void GeometryClose::process() {
    const Geometry* inputGeometry = inport_.getData();

    if (!enabled_.get()) {
        outport_.setData(inputGeometry, false);
        return;
    }

    std::unique_ptr<Geometry> outputGeometry = inputGeometry->clone();
    GlMeshGeometryBase* geometry = dynamic_cast<GlMeshGeometryBase*>(outputGeometry.get());
    if(geometry->getNumVertices() == 0) {
        LERROR("Geometry is empty!");
        outport_.setData(nullptr);
        return;
    }

    size_t numHoles = 0;
    if(auto geom = dynamic_cast<GlMeshGeometryUInt32Simple*>(geometry)) {
        numHoles = fill<GlMeshGeometryUInt32Simple>(geom, epsilon_.get());
    }
    else if(auto geom = dynamic_cast<GlMeshGeometryUInt32Color*>(geometry)) {
        numHoles = fill<GlMeshGeometryUInt32Color>(geom, epsilon_.get());
    }
    else if(auto geom = dynamic_cast<GlMeshGeometryUInt32Normal*>(geometry)) {
        numHoles = fill<GlMeshGeometryUInt32Normal>(geom, epsilon_.get());
    }
    else if(auto geom = dynamic_cast<GlMeshGeometryUInt32TexCoord*>(geometry)) {
        numHoles = fill<GlMeshGeometryUInt32TexCoord>(geom, epsilon_.get());
    }
    else if(auto geom = dynamic_cast<GlMeshGeometryUInt32NormalTexCoord*>(geometry)) {
        numHoles = fill<GlMeshGeometryUInt32NormalTexCoord>(geom, epsilon_.get());
    }
    else if(auto geom = dynamic_cast<GlMeshGeometryUInt32ColorNormal*>(geometry)) {
        numHoles = fill<GlMeshGeometryUInt32ColorNormal>(geom, epsilon_.get());
    }
    else if(auto geom = dynamic_cast<GlMeshGeometryUInt32ColorTexCoord*>(geometry)) {
        numHoles = fill<GlMeshGeometryUInt32ColorTexCoord>(geom, epsilon_.get());
    }
    else if(auto geom = dynamic_cast<GlMeshGeometryUInt32ColorNormalTexCoord*>(geometry)) {
        numHoles = fill<GlMeshGeometryUInt32ColorNormalTexCoord>(geom, epsilon_.get());
    }
    else {
        LERROR("Currently only GlMeshGeometryUint32* types supported");
        outport_.setData(nullptr);
        return;
    }

    if(numHoles > 0) {
        LINFO("Automatically closed " << numHoles << ((numHoles > 1) ? " holes." : " hole."));
        outputGeometry->setTransformationMatrix(inputGeometry->getTransformationMatrix());
        outport_.setData(outputGeometry.release());
    }
    else {
        LINFO("No holes found.");
        outport_.setData(inputGeometry, false);
        // outputGeometry will be removed automatically.
    }
}

}  //namespace
