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

#include "clipregiongeometrycreator.h"

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

namespace voreen {

ClipRegionGeometryCreator::ClipRegionGeometryCreator()
    : Processor()
    , inport_(Port::INPORT, "volumeinport", "Volume Inport")
    , outport_(Port::OUTPORT, "geometryoutport", "Mesh Outport")
    , boundingBox_("boundingbox", "Bounding Box (Clip Region)", tgt::IntBounds(tgt::ivec3(0), tgt::ivec3(1)), tgt::ivec3(0), tgt::ivec3(1))
{
    addPort(inport_);
    addPort(outport_);

    addProperty(boundingBox_);
}

Processor* ClipRegionGeometryCreator::create() const {
    return new ClipRegionGeometryCreator();
}

void ClipRegionGeometryCreator::process() {
    outport_.clear();

    GlMeshGeometryUInt16Simple* mesh = new GlMeshGeometryUInt16Simple();

    // get corners in voxel space
    tgt::vec3 llf(tgt::vec3(boundingBox_.get().getLLF()) - tgt::vec3(0.5f));
    tgt::vec3 urb(tgt::vec3(boundingBox_.get().getURB()) + tgt::vec3(0.5f));

    tgt::vec3 lrf(urb.x, llf.y, llf.z);
    tgt::vec3 lrb(urb.x, llf.y, urb.z);
    tgt::vec3 llb(llf.x, llf.y, urb.z);

    tgt::vec3 ulb(llf.x, urb.y, urb.z);
    tgt::vec3 ulf(llf.x, urb.y, llf.z);
    tgt::vec3 urf(urb.x, urb.y, llf.z);

    // create the actual geometry
    mesh->addVertex(urb);
    mesh->addVertex(urf);
    mesh->addVertex(ulf);
    mesh->addVertex(urb);
    mesh->addVertex(ulf);
    mesh->addVertex(ulb);

    mesh->addVertex(llf);
    mesh->addVertex(ulf);
    mesh->addVertex(urf);
    mesh->addVertex(llf);
    mesh->addVertex(urf);
    mesh->addVertex(lrf);

    mesh->addVertex(llf);
    mesh->addVertex(llb);
    mesh->addVertex(ulb);
    mesh->addVertex(llf);
    mesh->addVertex(ulb);
    mesh->addVertex(ulf);

    mesh->addVertex(urb);
    mesh->addVertex(ulb);
    mesh->addVertex(llb);
    mesh->addVertex(urb);
    mesh->addVertex(llb);
    mesh->addVertex(lrb);

    mesh->addVertex(urb);
    mesh->addVertex(lrb);
    mesh->addVertex(lrf);
    mesh->addVertex(urb);
    mesh->addVertex(lrf);
    mesh->addVertex(urf);

    mesh->addVertex(llf);
    mesh->addVertex(lrf);
    mesh->addVertex(lrb);
    mesh->addVertex(llf);
    mesh->addVertex(lrb);
    mesh->addVertex(llb);

    // set transformation from voxel to world space
    mesh->setTransformationMatrix(inport_.getData()->getVoxelToWorldMatrix());

    outport_.setData(mesh, true);
}

void ClipRegionGeometryCreator::adjustPropertiesToInput() {
    if (inport_.getData()){
        tgt::svec3 oldVolumeDimensions = tgt::svec3(boundingBox_.getMaxValue()) + tgt::svec3::one;

        tgt::svec3 volumeDimensions = inport_.getData()->getDimensions();

        // reset clip region if the volume dimensions have changed, else keep the value (e.g. for keeping the clipping for multiple time steps)
        if (volumeDimensions != oldVolumeDimensions) {
            // set max values and reset clip region
            boundingBox_.setMaxValue(tgt::ivec3(volumeDimensions) - tgt::ivec3::one);
            boundingBox_.set(tgt::IntBounds(tgt::ivec3(0), tgt::ivec3(volumeDimensions) - tgt::ivec3::one));
        }
    }
}

}   // namespace
