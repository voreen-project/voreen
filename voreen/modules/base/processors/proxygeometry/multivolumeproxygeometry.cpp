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

#include "multivolumeproxygeometry.h"

#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"
#include "voreen/core/datastructures/geometry/geometrysequence.h"

namespace voreen {

MultiVolumeProxyGeometry::MultiVolumeProxyGeometry()
    : Processor()
    , inport_(Port::INPORT, "volumehandle.volumehandle", "Volume Input", true)
    , outport_(Port::OUTPORT, "proxygeometry.geometry", "Proxy Geometry Output")
{
    addPort(inport_);
    addPort(outport_);
}

MultiVolumeProxyGeometry::~MultiVolumeProxyGeometry() {
}

Processor* MultiVolumeProxyGeometry::create() const {
    return new MultiVolumeProxyGeometry();
}

void MultiVolumeProxyGeometry::process() {
    GeometrySequence* seq = new GeometrySequence(true);

    std::vector<const VolumeBase*> data = inport_.getAllData();
    for(size_t d=0; d<data.size(); ++d) {
        if(!data[d])
            continue;

        const VolumeBase* volume = data[d];

        tgt::vec3 coordLlf(0, 0, 0);
        tgt::vec3 coordUrb(1, 1, 1);

        TriangleMeshGeometryColorNormal* mesh = new TriangleMeshGeometryColorNormal();
        mesh->addCube(VertexNormal(coordLlf, coordLlf), VertexNormal(coordUrb, coordUrb));
        mesh->transform(volume->getTextureToWorldMatrix());
        seq->addGeometry(mesh);
    }

    outport_.setData(seq);
}

} // namespace
