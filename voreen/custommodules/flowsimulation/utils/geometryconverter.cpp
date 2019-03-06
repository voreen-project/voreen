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

#include "geometryconverter.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include <cstdio>

namespace voreen {

bool exportGeometryToSTL(const Geometry* geometry, const std::string& path) {

    if(const GlMeshGeometryBase* data = dynamic_cast<const GlMeshGeometryBase*>(geometry)) {

        if(data->getPrimitiveType() != GL_TRIANGLES) {
            std::cout << "Currently only triangular meshes allowed" << std::endl;
            return false;
        }

        std::vector<VertexBase> vertices(data->getNumVertices());
        std::vector<uint32_t> indices;
        if(const GlMeshGeometryUInt32Simple* geom = dynamic_cast<const GlMeshGeometryUInt32Simple*>(geometry)) {
            vertices = geom->getVertices();
            indices = geom->getIndices();
        } else if(const GlMeshGeometryUInt32Normal* geom = dynamic_cast<const GlMeshGeometryUInt32Normal*>(geometry)) {
            std::copy(geom->getVertices().begin(), geom->getVertices().end(), vertices.begin());
            indices = geom->getIndices();
        } else if(const GlMeshGeometryUInt32NormalTexCoord* geom = dynamic_cast<const GlMeshGeometryUInt32NormalTexCoord*>(geometry)) {
            std::copy(geom->getVertices().begin(), geom->getVertices().end(), vertices.begin());
            indices = geom->getIndices();
        } else if(const GlMeshGeometryUInt32ColorNormal* geom = dynamic_cast<const GlMeshGeometryUInt32ColorNormal*>(geometry)) {
            std::copy(geom->getVertices().begin(), geom->getVertices().end(), vertices.begin());
            indices = geom->getIndices();
        } else if(const GlMeshGeometryUInt32TexCoord* geom = dynamic_cast<const GlMeshGeometryUInt32TexCoord*>(geometry)) {
            std::copy(geom->getVertices().begin(), geom->getVertices().end(), vertices.begin());
            indices = geom->getIndices();
        } else if(const GlMeshGeometryUInt32ColorNormalTexCoord* geom = dynamic_cast<const GlMeshGeometryUInt32ColorNormalTexCoord*>(geometry)) {
            std::copy(geom->getVertices().begin(), geom->getVertices().end(), vertices.begin());
            indices = geom->getIndices();
        } else {
            std::cout << "Unsupported geometry" << std::endl;
            return false;
        }

        std::ofstream f(path.c_str());
        f << "solid ascii " << path << "\n";

        tgt::mat4 m = geometry->getTransformationMatrix();
        if(data->usesIndexedDrawing()) {
            tgtAssert(data->getNumIndices() % 3 == 0, "No triangle mesh");
            for (size_t i = 0; i < data->getNumIndices(); i+=3) {

                tgt::vec3 v0 = m * vertices[indices[i+0]].pos_;
                tgt::vec3 v1 = m * vertices[indices[i+1]].pos_;
                tgt::vec3 v2 = m * vertices[indices[i+2]].pos_;

                tgt::vec3 normal = tgt::normalize(tgt::cross(v0-v1, v0-v2));

                f << "facet normal " << normal.x << " " << normal.y << " " << normal.z << "\n";
                f << "    outer loop\n";
                f << "        vertex " << v0.x << " " << v0.y << " " << v0.z << "\n";
                f << "        vertex " << v1.x << " " << v1.y << " " << v1.z << "\n";
                f << "        vertex " << v2.x << " " << v2.y << " " << v2.z << "\n";
                f << "    endloop\n";
                f << "endfacet\n";
            }
        }
        else {
            tgtAssert(data->getNumVertices() % 3 == 0, "No triangle mesh");
            for (size_t i = 0; i < data->getNumVertices(); i+=3) {

                tgt::vec3 v0 = m * vertices[i+0].pos_;
                tgt::vec3 v1 = m * vertices[i+1].pos_;
                tgt::vec3 v2 = m * vertices[i+2].pos_;

                tgt::vec3 normal = tgt::normalize(tgt::cross(v0-v1, v0-v2));

                f << "facet normal " << normal.x << " " << normal.y << " " << normal.z << "\n";
                f << "    outer loop\n";
                f << "        vertex " << v0.x << " " << v0.y << " " << v0.z << "\n";
                f << "        vertex " << v1.x << " " << v1.y << " " << v1.z << "\n";
                f << "        vertex " << v2.x << " " << v2.y << " " << v2.z << "\n";
                f << "    endloop\n";
                f << "endfacet\n";
            }
        }

        f.close();
        return true;
    }
    else {
        std::cout << "Geometry not supported!" << std::endl;
        return false;
    }
}

}
