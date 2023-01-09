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

#ifndef VRN_CURVECREATER_H
#define VRN_CURVECREATOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/processors/asynccomputeprocessor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

namespace voreen {

using CurveCreatorInput = const PointSegmentListGeometryVec3*;
using CurveCreatorOutput = std::unique_ptr<PointSegmentListGeometryVec3>;

/**
 * @brief The curvecreator is a processor that takes point segments (lists of lists of points) as
 * input and computes evenly spaced point segments on a curve through the input point segments. 
 * Every input segment corresponds to one output segment. Input sampling rate (how many of the input points
 * will be used) and output sampling rate (how dense should the output be) can be adjusted using the properties.
 * The input must be in voxel space. The output will also be in voxel space. 
 * 
 */
class VRN_CORE_API CurveCreator : public AsyncComputeProcessor<CurveCreatorInput, CurveCreatorOutput> {
public:
    enum class InterpolationMode : uint32_t {
        INTERPOLATE_CUBIC,
        INTERPOLATE_BEZIER,
        INTERPOLATE_RESAMPLE
    };

    CurveCreator();

    virtual Processor* create() const override { return new CurveCreator(); }

    virtual std::string getClassName() const override { return "CurveCreator"; }
    virtual std::string getCategory() const override { return "Geometry"; }
protected:
    ComputeInput prepareComputeInput() override;
	ComputeOutput compute(ComputeInput input, ProgressReporter &progressReporter) const override;
	void processComputeOutput(ComputeOutput output) override;
private:
    GeometryPort outport_;
    GeometryPort inport_;

    FloatProperty sampleLengthPorperty_;
    FloatProperty inputSamplingDensityProperty_;
    OptionProperty<InterpolationMode> interpolationModeProperty_;
};

static std::ostream& operator<<(std::ostream& stream, CurveCreator::InterpolationMode mode) {
    switch (mode)
    {
    case CurveCreator::InterpolationMode::INTERPOLATE_CUBIC:
        stream << "Cubic Interpolation";
        break;
    case CurveCreator::InterpolationMode::INTERPOLATE_BEZIER:
        stream << "Bezier Interpolation";
        break;
    case CurveCreator::InterpolationMode::INTERPOLATE_RESAMPLE:
        stream << "Resampling";
        break;
    default:
        stream << "Invalid Interpolation";
        break;
    }
    return stream;
}

}

#endif