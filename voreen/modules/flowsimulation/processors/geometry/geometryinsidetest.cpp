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

#include "geometryinsidetest.h"

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include <olb3D.h>
#ifndef OLB_PRECOMPILED
#include "olb3D.hh"
#endif
using namespace olb;
typedef double T;

namespace voreen {

const std::string GeometryInsideTest::loggerCat_("voreen.flowsimulation.GeometryInsideTest");

GeometryInsideTest::GeometryInsideTest()
    : Processor()
    , inport_(Port::INPORT, "geometryinsidetest.inport", "")
    , outport_(Port::OUTPORT, "geometryinsidetest.outport", "ID Volume Output", false)
    , method_("method", "Method")
    , dimensions_("dimensions", "Dimensions", 256, 32, 1024)
    , path_("path", "STL geometry", "Path", "", "STL (*.stl)", FileDialogProperty::OPEN_FILE)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(method_);
        method_.addOption("fast", "Fast", 0);
        method_.addOption("accurate", "Accurate", 1);
    addProperty(dimensions_);
    dimensions_.setTracking(false);
    addProperty(path_);
}

Processor* GeometryInsideTest::create() const {
    return new GeometryInsideTest();
}

bool GeometryInsideTest::isReady() const {
    if(!isInitialized()) {
        setNotReadyErrorMessage("Not initialized");
        return false;
    }

    if(!dynamic_cast<const GlMeshGeometryBase*>(inport_.getData())) {
        setNotReadyErrorMessage("Invalid input");
        return false;
    }

    return true;
}

void GeometryInsideTest::process() {

    const GlMeshGeometryBase* inputGeometry = dynamic_cast<const GlMeshGeometryBase*>(inport_.getData());
    tgtAssert(inputGeometry, "Invalid input");

    try {
        std::ofstream file(path_.get());
        inputGeometry->exportAsStl(file);
        file.close();
    }
    catch (std::exception& e) {
        LERROR("Failed to export mesh: " << e.what());
        outport_.setData(nullptr);
        return;
    }

    const int OUTSIDE = 0;
    const int INSIDE = 255;

    T spacing = tgt::max(inputGeometry->getBoundingBox(true).diagonal()) / dimensions_.get();
    STLreader<T> stlReader(path_.get(), spacing, 1.0, method_.getValue());
    Cuboid3D<T> cuboid(stlReader, spacing);
    BlockGeometry3D<T> geometry(cuboid, 0);
    geometry.rename(OUTSIDE, INSIDE, stlReader);

    tgt::svec3 dim(geometry.getNx(), geometry.getNy(), geometry.getNz());
    VolumeRAM_UInt8* idVolume = new VolumeRAM_UInt8(dim);

#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for
#endif
    for(int z=0; z<geometry.getNz(); z++) {
        for(int y=0; y<geometry.getNy(); y++) {
            for(int x=0; x<geometry.getNx(); x++) {
                idVolume->voxel(x, y, z) = geometry.getMaterial(x, y, z);
            }
        }
    }

    tgt::Vector3<T> offset = tgt::Vector3<T>::fromPointer(cuboid.getOrigin().data());
    Volume* outputVolume = new Volume(idVolume, tgt::vec3(spacing), offset);
    outport_.setData(outputVolume);
}

}