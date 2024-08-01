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

#ifndef VRN_INTEGRATORS_H
#define VRN_INTEGRATORS_H

#include "modules/flowanalysis/datastructures/streamline.h"
#include "modules/flowanalysis/utils/flowutils.h"

namespace voreen {

template<typename S, typename I, typename P>
class IntegratorBase {
public:
	IntegratorBase() = default;
	virtual ~IntegratorBase() = default;

	/**
	 * Integrates a line to completion starting from a seed point.
	 * @param seedPoint Seed point in world space.
	 */
	virtual Streamline integrate(const S &sampler,
								 const I &input,
								 const tgt::vec3 &seedPoint) const {
		Streamline line{};
		bool finished = false;

		P stepInput = this->initLine(sampler, input, line, seedPoint, finished);
		while (!finished) {
			this->integrateStep(sampler, input, line, stepInput, finished);
		}

		return line;
	}

	/**
	 * Initializes the integration of a line starting from a seed point.
	 * @param seedPoint Seed point in world space.
	 */
	virtual P initLine(const S &sampler,
									  const I &input,
									  Streamline &line,
									  const tgt::vec3 &seedPoint,
									  bool &finished) const = 0;
	virtual void integrateStep(const S &sampler,
							   const I &input,
							   Streamline &line,
							   P &stepInput,
							   bool &finished) const = 0;
};

//////////////////////// Streamline Integration ////////////////////////

struct StreamIntegratorInput {
	tgt::vec3 boundsMin;
	tgt::vec3 boundsMax;
	float stepSize;
	size_t upperLengthThreshold;
	tgt::vec2 absoluteMagnitudeThreshold;
	float stopIntegrationAngleThreshold;
	float velocityUnitConversion;
};

struct StreamIntegratorData {
	tgt::vec3 position;
	tgt::vec3 velocity;
};

class StreamIntegrator : public IntegratorBase<SpatialSampler,
											   StreamIntegratorInput,
											   StreamIntegratorData> {
public:
	StreamIntegrator() = default;
	~StreamIntegrator() noexcept override = default;

	StreamIntegratorData initLine(const SpatialSampler &sampler,
								  const StreamIntegratorInput &input,
								  Streamline &line,
								  const tgt::vec3 &seedPoint,
								  bool &finished) const override;
	void integrateStep(const SpatialSampler &sampler,
					   const StreamIntegratorInput &input,
					   Streamline &line,
					   StreamIntegratorData &stepInput,
					   bool &finished) const override;
};

}

#endif // VRN_INTEGRATORS_H
