#include "fieldplotdata.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumedisk.h"

#include <cmath>
#include <algorithm>

namespace voreen {

FieldPlotData::FieldPlotData(int width, int height, int numSlices)
{
    try {
        representation_ = new VolumeRAM_Float(tgt::svec3(width, height, numSlices));
    } catch(std::bad_alloc e) {
        LERRORC("voreen.fieldplotdata", "Failed to allocate plot data memory");
        throw;
    }

    memset(representation_->voxel(), 0, representation_->getNumBytes());
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

void FieldPlotData::putSingleMass (int x, float v, float minValue, float maxValue, int sliceNumber) {

    float y = (v-minValue)/(maxValue-minValue)*(getHeight() - 2.0f) + 0.5f;

    int   xPos = x;
    float yPos = y;

    int yPos1 = static_cast<int>(yPos - 0.5f);
    int yPos2 = static_cast<int>(yPos + 0.5f);

    if (xPos >= 0 && xPos <= getWidth()-1 && yPos1 >= 0 && yPos2 <= getHeight()-1) {
        representation_->voxel(xPos, yPos1, sliceNumber) += std::abs(yPos - yPos2);
        representation_->voxel(xPos, yPos2, sliceNumber) += std::abs(yPos - yPos1);
    }

}

void FieldPlotData::drawConnection(int x1, int x2, float v1, float v2, float minValue, float maxValue, int sliceNumber) {

    float y1 = (v1 - minValue) / (maxValue - minValue) * (getHeight() - 2.0f) + 0.5f;
    float y2 = (v2 - minValue) / (maxValue - minValue) * (getHeight() - 2.0f) + 0.5f;

    float gradient = (y2-y1) / (x2-x1);
    int lineLength = x2 - x1;

    for (int i=0; i<lineLength; i++)
    {
        int xPos = x1 + i;
        float yPos = y1 + gradient*i;

        int yPos1 = static_cast<int>(yPos - 0.5f);
        int yPos2 = static_cast<int>(yPos + 0.5f);

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

int FieldPlotData::getWidth()  const {
    return static_cast<int>(plotData_->getDimensions().x);
}
int FieldPlotData::getHeight() const {
    return static_cast<int>(plotData_->getDimensions().y);
}

void FieldPlotData::serialize(Serializer& s) const {

}

void FieldPlotData::deserialize(Deserializer& s) {

}

}
