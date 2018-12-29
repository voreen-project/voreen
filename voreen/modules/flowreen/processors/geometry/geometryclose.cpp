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

#include "geometryclose.h"

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include "../../ext/halfedge/trimesh.h"

namespace voreen {

const std::string GeometryClose::loggerCat_("voreen.base.GeometryClose");

    GeometryClose::GeometryClose()
        : Processor()
        , inport_(Port::INPORT, "geometry.geometry", "Geometry Input")
        , outport_(Port::OUTPORT, "geometry.clippedgeometry", "Clipped Geometry Output")
        , enabled_("enabled", "Enable", true)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(enabled_);
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
    GlMeshGeometryUInt32Normal* geometry = dynamic_cast<GlMeshGeometryUInt32Normal*>(outputGeometry.get());
    if(!geometry) {
        LERROR("Currently only GlMeshGeometryUInt32Normal supported!");
        outport_.setData(nullptr);
        return;
    }

    // Ensure the indices are set correctly.
    bool usesIndices = geometry->getNumIndices() > 0;
    geometry->createIndices(true);

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
            geometry->setIndices(std::vector<uint32_t>());
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
        VertexNormal center;

        while(true) {
            tgtAssert(boundary.find(current) != boundary.end(), "Not part of boundary");

            // Update ring.
            ring.push_back(current);
            boundary.erase(current);

            // Update center.
            center.pos_ += geometry->getVertex(current).pos_;
            center.normal_ += geometry->getVertex(current).normal_;

            // Find the next appropriate halfedge.
            trimesh::index_t halfEdgeIndex = -1;
            std::vector<trimesh::index_t> neighbours = mesh.vertex_vertex_neighbors(current);
            for(trimesh::index_t neighbour : neighbours) {
                if(boundary.find(neighbour) != boundary.end()) {
                    trimesh::index_t candidateHalfEdgeIndex = mesh.directed_edge2he_index(current, neighbour);
                    if(mesh.halfedge(candidateHalfEdgeIndex).face == -1) {
                        tgtAssert(halfEdgeIndex == -1, "more than one candidate found");
                        halfEdgeIndex = candidateHalfEdgeIndex;
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

        center.pos_ /= static_cast<float>(ring.size());
        center.normal_ = tgt::normalize(center.normal_);
        if(!usesIndices) {
            for (size_t i = 0; i < ring.size(); i++) {
                geometry->addVertex(geometry->getVertex(ring[i]));
                geometry->addVertex(geometry->getVertex(ring[(i + 1) % ring.size()]));
                geometry->addVertex(center);
            }
        }
        else {
            geometry->addVertex(center);
            size_t centerIndex = geometry->getVertices().size();
            for (size_t i = 0; i < ring.size(); i++) {
                geometry->addIndex(ring[i]);
                geometry->addIndex(ring[(i + 1) % ring.size()]);
                geometry->addIndex(centerIndex);
            }
        }

        start = boundary.begin();
        numHoles++;
    }
    tgtAssert(boundary.empty(), "Some boundary vertices unprocessed");

    LINFO("Automatically closed " << numHoles << " holes.");

    outport_.setData(outputGeometry.release());
}

}  //namespace
