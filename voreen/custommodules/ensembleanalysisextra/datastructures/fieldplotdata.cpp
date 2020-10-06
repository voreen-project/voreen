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

#include "fieldplotdata.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumedisk.h"

#include <cmath>
#include <algorithm>

namespace voreen {

FieldPlotData::FieldPlotData(size_t width, size_t height, size_t numSlices)
{
    try {
        representation_ = new VolumeRAM_Float(tgt::svec3(width, height, numSlices));
    } catch(std::bad_alloc& e) {
        LERRORC("voreen.fieldplotdata", "Failed to allocate plot data memory");
        throw;
    }

    representation_->clear();
    plotData_ = new Volume(representation_, tgt::vec3::one, tgt::vec3::zero); // takes ownership for representation
}

FieldPlotData::FieldPlotData(Volume* volume) {
    plotData_ = volume;
    tgtAssert(plotData_, "Volume was null");
    representation_ = dynamic_cast<VolumeRAM_Float*>(plotData_->getWritableRepresentation<VolumeRAM>());
    tgtAssert(representation_, "No RAM representation available");
}

FieldPlotData::~FieldPlotData() {
    delete plotData_;
}

void FieldPlotData::drawConnection(size_t x1, size_t x2, float v1, float v2, size_t sliceNumber) {

    // TODO: figure out why this is necessary.
    x1 = std::min(x1, getWidth() - 1);
    x2 = std::min(x2, getWidth() - 1);

    long y1 = tgt::clamp<long>(v1 * (getHeight()-1), 0, getHeight() - 1);
    long y2 = tgt::clamp<long>(v2 * (getHeight()-1), 0, getHeight() - 1);

    // Bresenham line algorithm.
    long dx = x2-x1;
    long dy = -std::abs(y2-y1);
    long err = dx+dy;
    long sy = y1<y2 ? 1 : -1;

    while (x1!=x2 || y1!=y2) {
        representation_->voxel(x1, y1, sliceNumber)++;
        long e2 = 2*err;
        if (e2 > dy) { err += dy; x1 += 1;  } /* e_xy+e_x > 0 */
        if (e2 < dx) { err += dx; y1 += sy; } /* e_xy+e_y < 0 */
    }

    // Note: the last x coord wont't be drawn.
}


Volume* FieldPlotData::getVolume() const {
    return plotData_;
}

size_t FieldPlotData::getWidth()  const {
    return plotData_->getDimensions().x;
}
size_t FieldPlotData::getHeight() const {
    return plotData_->getDimensions().y;
}

}
