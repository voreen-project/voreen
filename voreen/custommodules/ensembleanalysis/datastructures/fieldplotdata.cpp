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
    } catch(std::bad_alloc e) {
        LERRORC("voreen.fieldplotdata", "Failed to allocate plot data memory");
        throw;
    }

    representation_->clear();
    plotData_ = new Volume(representation_, tgt::vec3::one, tgt::vec3::zero); // takes ownership for representation
}

FieldPlotData::FieldPlotData(Volume* volume) {
    const VolumeRAM* representation = volume->getRepresentation<VolumeRAM>();
    tgtAssert(representation, "No RAM representation available");
    representation_ = dynamic_cast<VolumeRAM_Float*>(representation->clone());
    plotData_ = new Volume(representation_, volume);

}

FieldPlotData::~FieldPlotData() {
    delete plotData_;
}

void FieldPlotData::putSingleMass(size_t x, float v, size_t sliceNumber) {

    float y = v * (getHeight() - 2.0f) + 0.5f;

    size_t xPos = x;
    float yPos = y;

    size_t yPos1 = static_cast<size_t>(yPos - 0.5f);
    size_t yPos2 = static_cast<size_t>(yPos + 0.5f);

    if (xPos >= 0 && xPos <= getWidth()-1 && yPos1 >= 0 && yPos2 <= getHeight()-1) {
        representation_->voxel(xPos, yPos1, sliceNumber) += std::abs(yPos - yPos2);
        representation_->voxel(xPos, yPos2, sliceNumber) += std::abs(yPos - yPos1);
    }

}

void FieldPlotData::drawConnection(size_t x1, size_t x2, float v1, float v2, size_t sliceNumber) {

    float y1 = v1 * (getHeight() - 2.0f) + 0.5f;
    float y2 = v2 * (getHeight() - 2.0f) + 0.5f;

    float gradient = (y2-y1) / (x2-x1);
    size_t lineLength = x2 - x1;

    for (size_t i=0; i<lineLength; i++)
    {
        size_t xPos = x1 + i;
        float yPos = y1 + gradient*i;

        size_t yPos1 = static_cast<int>(yPos - 0.5f);
        size_t yPos2 = static_cast<int>(yPos + 0.5f);

        if (xPos >= 0 && xPos <= getWidth()-1 && yPos1 >= 0 && yPos2 <= getHeight()-1)
        {
            representation_->voxel(xPos, yPos1, sliceNumber) += std::abs(yPos - yPos2);
            representation_->voxel(xPos, yPos2, sliceNumber) += std::abs(yPos - yPos1);
        }
    }
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

void FieldPlotData::serialize(Serializer& s) const {

}

void FieldPlotData::deserialize(Deserializer& s) {

}

}
