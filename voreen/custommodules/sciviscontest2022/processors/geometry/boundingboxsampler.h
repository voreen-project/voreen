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

#ifndef VRN_BOUNDINGBOXSAMPLER_H
#define VRN_BOUNDINGBOXSAMPLER_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/intproperty.h"

namespace voreen {

class VRN_CORE_API BoundingBoxSampler : public Processor {
public:
	enum class BoundingBoxSide {
		Left,
		Right,
		Front,
		Back,
		Top,
		Bottom,
	};

	BoundingBoxSampler();
	Processor* create() const override { return new BoundingBoxSampler(); }

	std::string getClassName() const override { return "BoundingBoxSampler"; }
	std::string getCategory() const override { return "Geometry"; }
	CodeState getCodeState() const override   { return CODE_STATE_EXPERIMENTAL;  }

protected:
	void process() override;

private:
	GeometryPort inport_;
	GeometryPort outport_;

	IntProperty numLines_;
	IntProperty numSegments_;
	FloatProperty width_;
	FloatProperty distance_;
	FloatVec2Property height_;
	FloatVec3Property offset_;
	OptionProperty<BoundingBoxSide> side_;

	static const std::string loggerCat_;
};

static std::ostream& operator<<(std::ostream& stream, BoundingBoxSampler::BoundingBoxSide mode) {
	switch (mode)
	{
	case BoundingBoxSampler::BoundingBoxSide::Left:
		stream << "Left";
		break;
	case BoundingBoxSampler::BoundingBoxSide::Right:
		stream << "Right";
		break;
	case BoundingBoxSampler::BoundingBoxSide::Front:
		stream << "Front";
		break;
	case BoundingBoxSampler::BoundingBoxSide::Back:
		stream << "Back";
		break;
	case BoundingBoxSampler::BoundingBoxSide::Top:
		stream << "Top";
		break;
	case BoundingBoxSampler::BoundingBoxSide::Bottom:
		stream << "Bottom";
		break;
	default:
		stream << "Invalid Interpolation";
		break;
	}
	return stream;
}

}

#endif // VRN_BOUNDINGBOXSAMPLER_H
