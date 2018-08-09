/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "geometrytransformationvolume.h"

#include "voreen/core/datastructures/geometry/geometry.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

using tgt::mat4;

namespace voreen {

const std::string GeometryTransformationVolume::loggerCat_("voreen.base.GeometryTransformationVolume");

GeometryTransformationVolume::GeometryTransformationVolume()
  : Processor()
  , geometryInport_(Port::INPORT, "geometry.inport", "Geometry Input")
  , volumeInport_(Port::INPORT, "volume.inport", "Volume Input")
  , geometryOutport_(Port::OUTPORT, "geometry.outport", "Transformed Geometry Output")
  , enableProcessing_("enableProcessing", "Enable", true)
  , sourceCoordinateSystem_("sourceCoordinateSystem", "Source Coordinate System")
  , targetCoordinateSystem_("targetCoordinateSystem", "Dest Coordinate System")
  , forceUpdate_(true)
{
    addPort(geometryInport_);
    addPort(volumeInport_);
    addPort(geometryOutport_);

    addProperty(enableProcessing_);
    enableProcessing_.onChange(MemberFunctionCallback<GeometryTransformationVolume>(this, &GeometryTransformationVolume::forceUpdate));

    sourceCoordinateSystem_.addOption("voxel-coordinates", "Voxel Coordinates");
    sourceCoordinateSystem_.addOption("volume-coordinates", "Volume Coordinates");
    sourceCoordinateSystem_.addOption("world-coordinates", "World Coordinates");
    sourceCoordinateSystem_.addOption("texture-coordinates", "Texture Coordinates");
    sourceCoordinateSystem_.onChange(MemberFunctionCallback<GeometryTransformationVolume>(this, &GeometryTransformationVolume::forceUpdate));
    addProperty(sourceCoordinateSystem_);

    targetCoordinateSystem_.addOption("voxel-coordinates", "Voxel Coordinates");
    targetCoordinateSystem_.addOption("volume-coordinates", "Volume Coordinates");
    targetCoordinateSystem_.addOption("world-coordinates", "World Coordinates");
    targetCoordinateSystem_.addOption("texture-coordinates", "Texture Coordinates");
    targetCoordinateSystem_.select("world-coordinates");
    targetCoordinateSystem_.onChange(MemberFunctionCallback<GeometryTransformationVolume>(this, &GeometryTransformationVolume::forceUpdate));
    addProperty(targetCoordinateSystem_);

}

Processor* GeometryTransformationVolume::create() const {
    return new GeometryTransformationVolume();
}

void GeometryTransformationVolume::process() {
    // pass-through input geometry, if processing is disabled
    if (!enableProcessing_.get()) {
        geometryOutport_.setData(geometryInport_.getData(), false);
        forceUpdate_ = false;
        return;
    }

    // do nothing, if neither input data nor coordinate system selection has changed
    if (!geometryInport_.hasChanged() && !volumeInport_.hasChanged() && !forceUpdate_)
        return;

    // return if no input data present
    if (!geometryInport_.hasData() || !volumeInport_.hasData())
        return;

    // retrieve geometry from inport
    const Geometry* inputGeometry = geometryInport_.getData();
    tgtAssert(inputGeometry, "no input geometry");

    // retrieve volume from inport
    const VolumeBase* volumeHandle = volumeInport_.getData();
    tgtAssert(volumeHandle, "No volume");
    //tgtAssert(volumeHandle->getRepresentation<VolumeRAM>(), "No volume");

    //
    // determine transformation matrix
    //
    tgt::mat4 transformation = tgt::mat4::identity;
    if (targetCoordinateSystem_.isSelected("voxel-coordinates")) {
        if (sourceCoordinateSystem_.isSelected("voxel-coordinates")) {
            // voxel to voxel coordinates (no transformation necessary)
        } else if (sourceCoordinateSystem_.isSelected("volume-coordinates")) {
            transformation = volumeHandle->getPhysicalToVoxelMatrix();
        } else if (sourceCoordinateSystem_.isSelected("world-coordinates")) {
            transformation = volumeHandle->getWorldToVoxelMatrix();
        } else if (sourceCoordinateSystem_.isSelected("texture-coordinates")) {
            transformation = volumeHandle->getTextureToVoxelMatrix();
        }
    }
    // volume coordinates
    else if (targetCoordinateSystem_.isSelected("volume-coordinates")) {
        if (sourceCoordinateSystem_.isSelected("voxel-coordinates")) {
            transformation = volumeHandle->getVoxelToPhysicalMatrix();
        } else if (sourceCoordinateSystem_.isSelected("volume-coordinates")) {
            // No transformation necessary
        } else if (sourceCoordinateSystem_.isSelected("world-coordinates")) {
            transformation = volumeHandle->getWorldToPhysicalMatrix();
        } else if (sourceCoordinateSystem_.isSelected("texture-coordinates")) {
            transformation = volumeHandle->getTextureToPhysicalMatrix();
        }
    }
    // world coordinates
    else if (targetCoordinateSystem_.isSelected("world-coordinates")) {
        if (sourceCoordinateSystem_.isSelected("voxel-coordinates")) {
            transformation = volumeHandle->getVoxelToWorldMatrix();
        } else if (sourceCoordinateSystem_.isSelected("volume-coordinates")) {
            transformation = volumeHandle->getPhysicalToWorldMatrix();
        } else if (sourceCoordinateSystem_.isSelected("world-coordinates")) {
            // world to world coordinates (no transformation necessary)
        } else if (sourceCoordinateSystem_.isSelected("texture-coordinates")) {
            transformation = volumeHandle->getTextureToWorldMatrix();
        }
    }
    // texture coordinates: [0:1.0]^3
    else if (targetCoordinateSystem_.isSelected("texture-coordinates")) {
        if (sourceCoordinateSystem_.isSelected("voxel-coordinates")) {
            transformation = volumeHandle->getVoxelToTextureMatrix();
        } else if (sourceCoordinateSystem_.isSelected("volume-coordinates")) {
            transformation = volumeHandle->getPhysicalToTextureMatrix();
        } else if (sourceCoordinateSystem_.isSelected("world-coordinates")) {
            transformation = volumeHandle->getWorldToTextureMatrix();
        } else if (sourceCoordinateSystem_.isSelected("texture-coordinates")) {
            // no transformation necessary
        }
    }

    // create transformed output geometry
    std::unique_ptr<Geometry> outputGeometry = inputGeometry->clone();
    outputGeometry->transform(transformation);

    geometryOutport_.setData(outputGeometry.release());
    forceUpdate_ = false;
}

void GeometryTransformationVolume::forceUpdate() {
    forceUpdate_ = true;
}


} // namespace voreen
