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

#include "geometrysmoothnormals.h"

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include "../../ext/octree/Octree.hpp"

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

ACCESS_VEC(voreen::GlMeshGeometryUInt32Normal::VertexType)
ACCESS_VEC(voreen::GlMeshGeometryUInt32NormalTexCoord::VertexType)
ACCESS_VEC(voreen::GlMeshGeometryUInt32ColorNormal::VertexType)
ACCESS_VEC(voreen::GlMeshGeometryUInt32ColorNormalTexCoord::VertexType)

}
}

#undef ACCESS_VEC
#undef ACCESS_VEC_IDX

namespace voreen {


template<typename T>
void smooth(T* geometry, float epsilon) {

    typedef uint32_t I;
    typedef typename T::VertexType V;

    const std::vector<V>& vertices = geometry->getVertices();

    // Initializing the Octree with points from point cloud.
    unibn::Octree<V> octree;
    unibn::OctreeParams params;
    octree.initialize(vertices, params);

    std::unordered_set<I> seenAlready;
    for (I index = 0; index < vertices.size(); index++) {

        if(seenAlready.find(index) != seenAlready.end()) {
            continue;
        }

        const V& vertex = vertices[index];

        std::vector<I> neighbors;
        octree.template radiusNeighbors<unibn::L2Distance<V>>(vertex, epsilon, neighbors);
        neighbors.push_back(index);
        seenAlready.insert(neighbors.begin(), neighbors.end());

        // Calculate average normal.
        tgt::vec3 smoothNormal = tgt::vec3::zero;
        for(I idx : neighbors) {
            smoothNormal  += vertices[idx].normal_;
        }
        smoothNormal /= static_cast<float>(neighbors.size());

        // Set average normal.
        for(I idx : neighbors) {
            V vertex = vertices[idx];
            vertex.normal_ = smoothNormal;
            geometry->setVertex(idx, vertex);
        }
    }
}


const std::string GeometrySmoothNormals::loggerCat_("voreen.flowsimulation.GeometrySmoothNormals");

GeometrySmoothNormals::GeometrySmoothNormals()
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

GeometrySmoothNormals::~GeometrySmoothNormals()
{}

Processor* GeometrySmoothNormals::create() const {
    return new GeometrySmoothNormals();
}

void GeometrySmoothNormals::process() {
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

    if(auto geom = dynamic_cast<GlMeshGeometryUInt32Normal*>(geometry)) {
        smooth<GlMeshGeometryUInt32Normal>(geom, epsilon_.get());
    }
    else if(auto geom = dynamic_cast<GlMeshGeometryUInt32NormalTexCoord*>(geometry)) {
        smooth<GlMeshGeometryUInt32NormalTexCoord>(geom, epsilon_.get());
    }
    else if(auto geom = dynamic_cast<GlMeshGeometryUInt32ColorNormal*>(geometry)) {
        smooth<GlMeshGeometryUInt32ColorNormal>(geom, epsilon_.get());
    }
    else if(auto geom = dynamic_cast<GlMeshGeometryUInt32ColorNormalTexCoord*>(geometry)) {
        smooth<GlMeshGeometryUInt32ColorNormalTexCoord>(geom, epsilon_.get());
    }
    else {
        LERROR("Currently only GlMeshGeometryUint32*Normal* types supported");
        outport_.setData(nullptr);
        return;
    }

    outport_.setData(outputGeometry.release());
}

}  //namespace
