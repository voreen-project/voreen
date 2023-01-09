/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "boundingboxsampler.h"

#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

namespace voreen {

const std::string BoundingBoxSampler::loggerCat_("voreen.sciviscontest2022.BoundingBoxSampler");

BoundingBoxSampler::BoundingBoxSampler()
	: Processor()
	, inport_(Port::INPORT, "inport", "Geometry Input", false)
	, outport_(Port::OUTPORT, "outport", "PointSegmentList Output", false)
	, numLines_("numLines", "Number of Lines:", 1, 1, 99)
	, numSegments_("numSegments", "Number of Points on each line:", 2, 2, 999)
	, width_("width", "Line length:", 1.0f, 0.1f, 99.9f)
	, distance_("distance", "Distance from face::", 0.0f, -9999.9, 9999.9f)
	, height_("height", "Sampling region:", tgt::vec2(0.0f, 1.0f))
	, offset_("offset", "Offset:", tgt::vec3::zero, tgt::vec3(-9999.9f), tgt::vec3(9999.9f))
	, side_("side", "Bounding box side")
{
	Processor::addPort(this->inport_);
	Processor::addPort(this->outport_);

	Processor::addProperty(this->numLines_);
	Processor::addProperty(this->numSegments_);
	Processor::addProperty(this->width_);
	Processor::addProperty(this->distance_);
	Processor::addProperty(this->height_);
	Processor::addProperty(this->offset_);

	ON_CHANGE_LAMBDA(this->height_, [&]() {
		auto& value = this->height_.get();
		auto corrected = value;
		corrected.x = std::min(value.x, value.y);
		corrected.y = std::max(value.x, value.y);

		if (value != corrected) {
			this->height_.set(corrected);
		}
	});

	this->side_.addOption("Left", "Sample from the left", BoundingBoxSide::Left);
	this->side_.addOption("Right", "Sample from the right", BoundingBoxSide::Right);
	this->side_.addOption("Front", "Sample from the front", BoundingBoxSide::Front);
	this->side_.addOption("Back", "Sample from the back", BoundingBoxSide::Back);
	this->side_.addOption("Top", "Sample from the top", BoundingBoxSide::Top);
	this->side_.addOption("Bottom", "Sample from the bottom", BoundingBoxSide::Bottom);
	Processor::addProperty(this->side_);
}

void BoundingBoxSampler::process() {
	auto inGeometry = this->inport_.getData();
	if (!inGeometry) {
		LERROR("Invalid input geometry");
	}

	float width;
	float height;
	tgt::vec3 up = tgt::vec3::zero;
	tgt::vec3 back = tgt::vec3::zero;
	tgt::vec3 right = tgt::vec3::zero;

	auto boundingBox = inGeometry->getBoundingBox();
	tgt::vec3 faceCenter = boundingBox.center();

	auto side = this->side_.getValue();
	switch (side) {
	case BoundingBoxSide::Left:
		width = boundingBox.getURB().y - boundingBox.getLLF().y;
		height = boundingBox.getURB().z - boundingBox.getLLF().z;

		faceCenter.x = boundingBox.getLLF().x;

		up = tgt::vec3(0.0f, 0.0f, 1.0f);
		back = tgt::vec3(-1.0f, 0.0f, 0.0f);
		right = tgt::vec3(0.0f, 1.0f, 0.0f);
		break;
	case BoundingBoxSide::Right:
		width = boundingBox.getURB().y - boundingBox.getLLF().y;
		height = boundingBox.getURB().z - boundingBox.getLLF().z;

		faceCenter.x = boundingBox.getURB().x;

		up = tgt::vec3(0.0f, 0.0f, 1.0f);
		back = tgt::vec3(1.0f, 0.0f, 0.0f);
		right = tgt::vec3(0.0f, -1.0f, 0.0f);
		break;
	case BoundingBoxSide::Front:
		width = boundingBox.getURB().x - boundingBox.getLLF().x;
		height = boundingBox.getURB().z - boundingBox.getLLF().z;

		faceCenter.y = boundingBox.getLLF().y;

		up = tgt::vec3(0.0f, 0.0f, 1.0f);
		back = tgt::vec3(0.0f, -1.0f, 0.0f);
		right = tgt::vec3(1.0f, 0.0f, 0.0f);
		break;
	case BoundingBoxSide::Back:
		width = boundingBox.getURB().x - boundingBox.getLLF().x;
		height = boundingBox.getURB().z - boundingBox.getLLF().z;

		faceCenter.y = boundingBox.getURB().y;

		up = tgt::vec3(0.0f, 0.0f, 1.0f);
		back = tgt::vec3(0.0f, 1.0f, 0.0f);
		right = tgt::vec3(-1.0f, 0.0f, 0.0f);
		break;
	case BoundingBoxSide::Top:
		width = boundingBox.getURB().x - boundingBox.getLLF().x;
		height = boundingBox.getURB().y - boundingBox.getLLF().y;

		faceCenter.z = boundingBox.getURB().z;

		up = tgt::vec3(0.0f, 1.0f, 0.0f);
		back = tgt::vec3(0.0f, 0.0f, 1.0f);
		right = tgt::vec3(1.0f, 0.0f, 0.0f);
		break;
	case BoundingBoxSide::Bottom:
		width = boundingBox.getURB().x - boundingBox.getLLF().x;
		height = boundingBox.getURB().y - boundingBox.getLLF().y;

		faceCenter.z = boundingBox.getLLF().z;

		up = tgt::vec3(0.0f, 1.0f, 0.0f);
		back = tgt::vec3(0.0f, 0.0f, -1.0f);
		right = tgt::vec3(1.0f, 0.0f, 0.0f);
		break;
	}

	float regionHeight = (this->height_.get().y - this->height_.get().x) * height;
	tgt::vec3 faceBottom = faceCenter - ((height / 2.0f) * up) + (this->height_.get().x * height * up) + this->offset_.get();

	tgt::vec3 startLeft = faceBottom - (this->width_.get() * width / 2.0f * right) + (back * this->distance_.get());
	tgt::vec3 startRight = faceBottom + (this->width_.get() * width / 2.0f * right) + (back * this->distance_.get());

	tgt::vec3 endLeft = startLeft + (regionHeight * up);
	tgt::vec3 endRight = startRight + (regionHeight * up);

	auto numLines = this->numLines_.get();
	auto segmentVector = numLines == 1 ? tgt::vec3::zero : ((endLeft - startLeft) / ((float)numLines - 1.0f));

	auto currLeft = startLeft;
	auto currRight = startRight;
	std::unique_ptr<PointSegmentListGeometryVec3> output{ new PointSegmentListGeometryVec3() };
	for (int i = 0; i < numLines; ++i) {
		std::vector<tgt::vec3> segment {};
		auto lineSegmentVector = (currRight - currLeft) / (float)(this->numSegments_.get() - 1);
		for (int j = 0; j < this->numSegments_.get(); ++j) {
			auto point = currLeft + ((float)j * lineSegmentVector);
			segment.push_back(point);
		}

		output->addSegment(std::move(segment));
		currLeft += segmentVector;
		currRight += segmentVector;
	}
	this->outport_.setData(output.release());
}

}