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
#include "curvecreator.h"
#include "../../datastructures/curve.h"

#include <random>

namespace voreen {

CurveCreator::CurveCreator() :
    AsyncComputeProcessor()
    , outport_(Port::OUTPORT, "curvecreator.out", "PointList Output", false)
    , inport_(Port::INPORT, "curvecreator.in", "PointList Input", false)
    , sampleLengthPorperty_("curvecreator.sampleLength", "Samplingrate: ", 5.f, 0.05f, 100.f)
    , inputSamplingDensityProperty_("curvecreator.inputSamplindDensity", "Input Sampling Density (%): ", 12.5f, 1.f, 100.f)
    , interpolationModeProperty_("curvecreator.interpolationMode", "Interpolation Mode: ") {

    addPort(outport_);
    addPort(inport_);

    interpolationModeProperty_.addOption("Cubic", "Interpolate using cubic splines", InterpolationMode::INTERPOLATE_CUBIC);
    interpolationModeProperty_.addOption("Resample", "Simply resample the input points evenly", InterpolationMode::INTERPOLATE_RESAMPLE);
    interpolationModeProperty_.setDefaultValue("Cubic");

    ON_CHANGE_LAMBDA(inport_, [&]() { if (isReady()) process(); });
    
    addProperty(sampleLengthPorperty_);
    addProperty(inputSamplingDensityProperty_);
    addProperty(interpolationModeProperty_);
}

CurveCreatorInput CurveCreator::prepareComputeInput() 
{
    return static_cast<const PointSegmentListGeometryVec3*>(inport_.getData());
}

CurveCreatorOutput CurveCreator::compute(ComputeInput input, ProgressReporter &progressReporter) const 
{
    CurveCreatorOutput output{ new CurveCreatorOutput::element_type() };
    const auto samplingLength = 1.f / sampleLengthPorperty_.get();
    // turn [0,100] percentages to [0,1] percentages:
    const auto percentageOfInputPoints = inputSamplingDensityProperty_.get() / 100.f;
    // use only every 1/percentage point from inpute:
    auto samplingDensity = (1 / percentageOfInputPoints);

    auto mode = interpolationModeProperty_.getValue();

    auto progressPerSegment = 1.f / input->getNumSegments();

    for (auto i = 0u; i < static_cast<unsigned int>(input->getNumSegments()); ++i) {
        const auto& segment = input->getSegment(i);
        std::vector<tgt::vec3> sampledPoints;

        // sample input points
        for (auto j = 0u; j < segment.size(); ++j) {
            auto index = j * samplingDensity;
            if (index >= segment.size())
                break;

            sampledPoints.emplace_back(segment.at(index));
        }
        
        std::unique_ptr<Curve> curve;
        if (mode == InterpolationMode::INTERPOLATE_CUBIC) 
        {
            curve.reset(new CubicCurve(std::move(sampledPoints)));
        }
        if (mode == InterpolationMode::INTERPOLATE_RESAMPLE)
        {
            sampledPoints = Curve::resamplePoints(sampledPoints);
            output->addSegment(std::move(sampledPoints));
            continue;
        }
        
        auto totalLength = curve->getTotalLength();
        std::vector<tgt::vec3> resultSegment;
        for (auto t = 0.f; t < totalLength; t += samplingLength) {
            resultSegment.emplace_back(curve->sample(t));
            progressReporter.setProgress(i * progressPerSegment + (t / totalLength) * progressPerSegment);
        }
        resultSegment = Curve::resamplePoints(resultSegment);
        output->addSegment(std::move(resultSegment));
    }
    progressReporter.setProgress(1.f);
    return output;
}

void CurveCreator::processComputeOutput(ComputeOutput output) {
    outport_.setData(output.release());
}

}
