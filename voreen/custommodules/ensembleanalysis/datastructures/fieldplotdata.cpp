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

    long y1 = static_cast<long>(v1 * (getHeight()-1));
    long y2 = static_cast<long>(v2 * (getHeight()-1));

    long dx = x2-x1;
    long dy = -std::abs(y2-y1), sy = y1<y2 ? 1 : -1;
    long err = dx+dy;

    while (true) {
        if (x1==x2 && y1==y2) break;
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
