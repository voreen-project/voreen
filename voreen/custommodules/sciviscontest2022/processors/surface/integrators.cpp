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

#include "integrators.h"

namespace voreen {

//////////////////////// Streamline Integration ////////////////////////

StreamIntegratorData StreamIntegrator::initLine(const SpatialSampler &sampler,
												const StreamIntegratorInput &input,
												Streamline &line,
												const tgt::vec3 &seedPoint,
												bool &finished) const {
	tgtAssert(finished==false, "finished must be set to false");

	// Sample the velocity of the flow at the seed point.
	auto velocity = sampler.sample(seedPoint);

	// Return an empty line in case the initial velocity is
	// zero and the line point is not affected by a flow.
	if (velocity==tgt::vec3::zero) {
		finished = true;
		return StreamIntegratorData{
			seedPoint,
			velocity
		};
	}

	// Add initial element.
	line.addElementAtEnd(Streamline::StreamlineElement(seedPoint, velocity));
	return StreamIntegratorData{
		seedPoint,
		velocity
	};
}

void StreamIntegrator::integrateStep(const SpatialSampler &sampler,
									 const StreamIntegratorInput &input,
									 Streamline &line,
									 StreamIntegratorData &stepInput,
									 bool &finished) const {
	tgtAssert(finished==false, "finished must be set to false");
	constexpr float epsilon = 1e-5f; // std::numeric_limits<float>::epsilon() is not enough.

	// Execute 4th order Runge-Kutta step.
	tgt::vec3 k1 = stepInput.velocity * input.stepSize * input.velocityUnitConversion;
	tgt::vec3 k2 = sampler.sample(stepInput.position + (k1/2.0f)) * input.stepSize * input.velocityUnitConversion;
	tgt::vec3 k3 = sampler.sample(stepInput.position + (k2/2.0f)) * input.stepSize * input.velocityUnitConversion;
	tgt::vec3 k4 = sampler.sample(stepInput.position + k3) * input.stepSize * input.velocityUnitConversion;
	tgt::vec3 dPos = ((k1/6.0f) + (k2/3.0f) + (k3/3.0f) + (k4/6.0f));

	// Update position.
	stepInput.position += dPos;

	// Check constrains.
	auto clamped = tgt::clamp(stepInput.position, input.boundsMin, input.boundsMax);
	bool hasStep = (stepInput.position==clamped); // Is in bounds.
	hasStep &= (stepInput.position!=line.getLastElement().position_); // Position has changed.

	stepInput.velocity = sampler.sample(stepInput.position);
	float magnitude = tgt::length(stepInput.velocity);
	hasStep &= (stepInput.velocity!=tgt::vec3::zero)
		&& (magnitude > input.absoluteMagnitudeThreshold.x - epsilon)
		&& (magnitude < input.absoluteMagnitudeThreshold.y + epsilon)
		&& std::acos(std::abs(tgt::dot(line.getLastElement().velocity_, stepInput.velocity))/
			(tgt::length(line.getLastElement().velocity_)*magnitude)) <= input.stopIntegrationAngleThreshold;

	if (hasStep) {
		line.addElementAtEnd(Streamline::StreamlineElement(stepInput.position, stepInput.velocity));
		if (line.getNumElements() >= input.upperLengthThreshold) {
			finished = true;
		}
	} else {
		finished = true;
	}
}

}