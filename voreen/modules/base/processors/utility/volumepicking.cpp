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

#include "volumepicking.h"

namespace voreen {

const std::string VolumePicking::loggerCat_("voreen.base.VolumePicking");

VolumePicking::VolumePicking()
    : RenderProcessor(),
      inportFirstHit_(Port::INPORT, "inport.firsthit", "First-hit Image"),
      inportVolume_(Port::INPORT, "inport.volume", "Volume Input"),
      outportPicked_(Port::OUTPORT, "outport.picked", "Picked Points Geometry"),
      outCoordinateSystem_("outCoordinateSystem", "Output Coordinate System", INVALID_PARAMETERS),
      clearLast_("clearLast", "Clear Last"),
      clearAll_("clearAll", "Clear All"),
      numPicked_("numPicked", "Number of Points", 0, 0, 1000, Processor::VALID),
      lastPicked_("lastPicked", "Last Picked", tgt::vec3(0.f), tgt::vec3(-10000.f), tgt::vec3(10000.f),
                  Processor::VALID),
      pickEvent_("pickEvent", "Picking", this, &VolumePicking::mouseEvent, tgt::MouseEvent::MOUSE_BUTTON_LEFT,
                 tgt::MouseEvent::PRESSED, tgt::Event::CTRL),
      lastMousePosition_(tgt::vec2(-1))
{
    addPort(inportFirstHit_);
    addPort(inportVolume_);
    addPort(outportPicked_);

    outCoordinateSystem_.addOption("voxel",     "Voxel Coordinates");
    outCoordinateSystem_.addOption("volume",    "Volume Coordinates");
    outCoordinateSystem_.addOption("world",     "World Coordinates");
    outCoordinateSystem_.addOption("texture",   "Texture Coordinates");
    outCoordinateSystem_.select("voxel");

    clearLast_.onClick(MemberFunctionCallback<VolumePicking>(this, &VolumePicking::clearLast));
    clearAll_.onClick(MemberFunctionCallback<VolumePicking>(this, &VolumePicking::clearAll));
    numPicked_.setReadOnlyFlag(true);
    lastPicked_.setReadOnlyFlag(true);

    addProperty(outCoordinateSystem_);
    addProperty(numPicked_);
    addProperty(lastPicked_);
    addProperty(clearLast_);
    addProperty(clearAll_);
    addEventProperty(pickEvent_);
}

Processor* VolumePicking::create() const {
    return new VolumePicking();
}

void VolumePicking::process() {

    if (lastMousePosition_ != tgt::ivec2(-1)) {
        tgt::vec3 pickedPos = getPickedPosition(lastMousePosition_);
        pickedPositions_.push_back(pickedPos);
        LINFO("Picked pos: " << pickedPos);
        lastMousePosition_ = tgt::ivec2(-1);
    }

    PointListGeometryVec3* transformedPositions = getTransformedPositions();

    if (transformedPositions->getNumPoints() > 0)
        lastPicked_.set(transformedPositions->getData().back());
    else
        lastPicked_.set(tgt::vec3(0.f));

    outportPicked_.setData(transformedPositions);

    updateProperties();
}

void VolumePicking::invalidate(int inv) {
    if (inv == INVALID_RESULT)
        return;

    RenderProcessor::invalidate(inv);
}

void VolumePicking::serialize(Serializer& s) const  {
    RenderProcessor::serialize(s);

    s.serialize("pickedPositions", pickedPositions_);
}

void VolumePicking::deserialize(Deserializer& d) {
    RenderProcessor::deserialize(d);

    try {
        d.deserialize("pickedPositions", pickedPositions_);
    }
    catch (SerializationNoSuchDataException&) {
        d.removeLastError();
    }
}

tgt::vec3 VolumePicking::getPickedPosition(tgt::ivec2 mousePosition) {
    tgtAssert(inportVolume_.hasData() && inportFirstHit_.hasRenderTarget(), "inports not ready");

    tgt::ivec2 dims = inportFirstHit_.getRenderTarget()->getSize();
    if (tgt::hor(tgt::lessThan(mousePosition, tgt::ivec2(0))) || tgt::hor(tgt::greaterThanEqual(mousePosition, dims)) ) {
        LWARNING("Invalid mouse position" << mousePosition);
        return tgt::vec3(-1.f);
    }
    mousePosition.y = dims.y - mousePosition.y;

    // read first hit position from input texture
    tgt::vec3 pickedPos = inportFirstHit_.getRenderTarget()->getColorAtPos(mousePosition).xyz();

    return pickedPos;
}

PointListGeometryVec3* VolumePicking::getTransformedPositions() const {
    tgtAssert(inportVolume_.hasData() && inportFirstHit_.hasRenderTarget(), "inports not ready");

    std::vector<tgt::vec3> transformedPoints = pickedPositions_;

    // voxel coordinates
    if (outCoordinateSystem_.isSelected("voxel")) {
        tgt::vec3 volDims = tgt::vec3(inportVolume_.getData()->getDimensions() - tgt::svec3(1));
        for (size_t i=0; i<transformedPoints.size(); i++)
            transformedPoints[i] *= volDims;
    }
    // volume coordinates
    else if (outCoordinateSystem_.isSelected("volume")) {
        tgt::mat4 textureToPhysicalTrafo = inportVolume_.getData()->getTextureToPhysicalMatrix();
        for (size_t i=0; i<transformedPoints.size(); i++)
            transformedPoints[i] = textureToPhysicalTrafo * transformedPoints[i];
    }
    // world coordinates
    else if (outCoordinateSystem_.isSelected("world")) {
        tgt::vec3 volDims = tgt::vec3(inportVolume_.getData()->getDimensions()-tgt::svec3(1));
        tgt::mat4 voxelToWorldTrafo = inportVolume_.getData()->getVoxelToWorldMatrix();
        for (size_t i=0; i<transformedPoints.size(); i++) {
            transformedPoints[i] = voxelToWorldTrafo * (transformedPoints[i] * volDims);
        }
    }
    // texture coordinates: [0:1.0]^3 (no transformation necessary)
    else if (outCoordinateSystem_.isSelected("texture")) {
    }

    PointListGeometryVec3* result = new PointListGeometryVec3();
    result->setData(transformedPoints);

    return result;
}

void VolumePicking::mouseEvent(tgt::MouseEvent* e) {
    lastMousePosition_ = e->coord();
    invalidate(INVALID_RESULT+1);
}

void VolumePicking::clearLast() {
    if (pickedPositions_.size() > 0) {
        pickedPositions_.pop_back();
        updateProperties();
        invalidate(INVALID_RESULT+1);
    }
}

void VolumePicking::clearAll() {
    if (pickedPositions_.size() > 0) {
        pickedPositions_.clear();
        updateProperties();
        invalidate(INVALID_RESULT+1);
    }
}

void VolumePicking::updateProperties() {
    numPicked_.set(static_cast<int>(pickedPositions_.size()));
}

} // voreen namespace
