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

#include "geometrymerge.h"

#include "voreen/core/datastructures/geometry/geometrysequence.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

namespace {

using namespace voreen;

template<typename T>
struct TypedGeometryPair {

    TypedGeometryPair(const Geometry* lhs, const Geometry* rhs)
            : lhs_(dynamic_cast<const T*>(lhs)), rhs_(dynamic_cast<const T*>(rhs)) {}

    operator bool() {
        return lhs_ != nullptr && rhs_ != nullptr;
    }

    const T* lhs_;
    const T* rhs_;
};

template<typename T>
std::unique_ptr<GlMeshGeometryBase> mergeGeometriesTyped(const T* lhs, const T* rhs) {
    tgtAssert(lhs, "lhs null");
    tgtAssert(rhs, "rhs null");

    auto merged = lhs->clone().release();
    auto mergedTyped = static_cast<T*>(merged);

    // Transform all vertex position to world space.
    auto vertices = mergedTyped->getVertices();
    for (auto& vertex: vertices) {
        vertex.pos_ = mergedTyped->getTransformationMatrix() * vertex.pos_;
    }
    mergedTyped->setTransformationMatrix(tgt::mat4::identity);

    // Copy over and transform vertices.
    for (auto vertex: rhs->getVertices()) {
        vertex.pos_ = rhs->getTransformationMatrix() * vertex.pos_;
        vertices.emplace_back(vertex);
    }
    mergedTyped->setVertices(vertices);

    // Add and adjust indices, if required.
    if (lhs->usesIndexedDrawing()) {
        auto indexOffset = lhs->getNumVertices();
        for (auto index: rhs->getIndices()) {
            mergedTyped->addIndex(index + indexOffset);
        }
    }

    return std::unique_ptr<GlMeshGeometryBase>(mergedTyped);
}


std::unique_ptr<GlMeshGeometryBase> mergeGeometries(const GlMeshGeometryBase* lhs, const GlMeshGeometryBase* rhs) {
    if (lhs->getPrimitiveType() != rhs->getPrimitiveType()) {
        return nullptr;
    }

    if (lhs->getVertexLayout() != rhs->getVertexLayout()) {
        return nullptr;
    }

    if (lhs->getIndexType() != rhs->getIndexType()) {
        return nullptr;
    }

    if (lhs->usesIndexedDrawing() != rhs->usesIndexedDrawing()) {
        return nullptr;
    }

    if (auto typedGeometry = TypedGeometryPair<GlMeshGeometryUInt32Normal>(lhs, rhs)) {
        return mergeGeometriesTyped(typedGeometry.lhs_, typedGeometry.rhs_);
    }

    return nullptr;
}

}

namespace voreen {

const std::string GeometryMerge::loggerCat_("voreen.flowsimulation.GeometryMerge");

GeometryMerge::GeometryMerge()
    : Processor()
    , inport_(Port::INPORT, "geometry.input", "Geometry Input")
    , outport_(Port::OUTPORT, "geometry.output", "Geometry Output")
    , enabled_("enabled", "Enable", true)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(enabled_);
}

GeometryMerge::~GeometryMerge()
{}

Processor* GeometryMerge::create() const {
    return new GeometryMerge();
}

void GeometryMerge::process() {
    const Geometry* inputGeometry = inport_.getData();

    if (!enabled_.get()) {
        outport_.setData(inputGeometry, false);
        return;
    }

    if(auto geometrySequence = dynamic_cast<const GeometrySequence*>(inport_.getData())) {

        std::vector<const GlMeshGeometryBase*> geometries;

        for(size_t i=0; i<geometrySequence->getNumGeometries(); i++) {
            auto geometry = dynamic_cast<const GlMeshGeometryBase*>(geometrySequence->getGeometry(i));
            if(geometry) {
                geometries.push_back(geometry);
            }
            else {
                LERROR("GeometrySequence contains non-GlMeshGeometry");
                outport_.setData(nullptr);
                return;
            }
        }

        if(geometries.empty()) {
            LERROR("No GlMeshGeometry input");
            outport_.setData(nullptr);
        }
        else if(geometries.size() == 1) {
            outport_.setData(geometries.front(), false);
        }
        else {
            auto geometry = geometries.front();
            for (size_t i=1; i<geometries.size(); i++) {
                auto combined = mergeGeometries(geometry, geometries.at(i));

                if(!combined) {
                    LERROR("Encountered Mesh of unexpected Type");
                    return;
                }

                geometry = combined.release();
            }
            outport_.setData(geometry, true);
        }
    }
    else if(auto geometryData = dynamic_cast<const GlMeshGeometryBase*>(inport_.getData())) {
        outport_.setData(geometryData, false);
    }

    outport_.setData(nullptr);
}

}  //namespace
