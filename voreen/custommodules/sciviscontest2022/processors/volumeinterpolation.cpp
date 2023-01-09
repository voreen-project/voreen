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

#include "volumeinterpolation.h"

#include "voreen/core/ports/conditions/portconditionvolumetype.h"

#include <memory>

namespace voreen {

const std::string VolumeInterpolation::loggerCat_("voreen.sciviscontest2022.VolumeInterpolation");

VolumeInterpolation::VolumeInterpolation()
	: Processor()
	, inputVolumes_(Port::INPORT, "inputVolumes", "Volume List Input")
	, outputVolume_(Port::OUTPORT, "outputVolume", "Output")
	, timeStep_("timeStep", "Interpolation step", 0.0f, 0.0f, 1.0f)
{
	inputVolumes_.onChange(MemberFunctionCallback<VolumeInterpolation>(this, &VolumeInterpolation::updateTimeRange));
	timeStep_.onChange(MemberFunctionCallback<VolumeInterpolation>(this, &VolumeInterpolation::forceReload));

	Processor::addPort(inputVolumes_);
	Processor::addPort(outputVolume_);

	Processor::addProperty(timeStep_);
}

Processor* VolumeInterpolation::create() const {
	return new VolumeInterpolation();
}

std::unique_ptr<VolumeBase> interpolate(float timeStep, const VolumeList* input) {
	std::size_t firstIdx = std::min(static_cast<std::size_t>(std::floor(timeStep)), input->size());
	std::size_t secondIdx = std::min(static_cast<std::size_t>(std::ceil(timeStep)), input->size());

	float alpha = timeStep - firstIdx;

	const VolumeBase* firstVolume = input->at(firstIdx);
	const VolumeBase* secondVolume = input->at(secondIdx);

	if (firstIdx == secondIdx) {
		return std::unique_ptr<VolumeBase>(firstVolume->clone());
	}

	VolumeRAMRepresentationLock firstVolumeData(firstVolume);
	VolumeRAMRepresentationLock secondVolumeData(secondVolume);

	tgt::Vector3<long> dimensions = tgt::min(firstVolumeData->getDimensions(), secondVolumeData->getDimensions());
	RealWorldMapping firstRWM = firstVolume->getRealWorldMapping();
	RealWorldMapping secondRWM = secondVolume->getRealWorldMapping();

	auto outputVolume = std::unique_ptr<VolumeAtomic<tgt::Vector3<float>>>{ new VolumeAtomic<tgt::Vector3<float>>(dimensions) };

#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for default(none) shared(dimensions, firstVolumeData, secondVolumeData, firstRWM, secondRWM, outputVolume, alpha)
#endif
	for (long z = 0; z < dimensions.z; z++) {
		for (long y = 0; y < dimensions.y; y++) {
			for (long x = 0; x < dimensions.x; x++) {
				tgt::svec3 pos(x, y, z);

				tgt::Vector3<float> value {};
				for (std::size_t channel = 0; channel < 3; ++channel) {
					float first = firstVolumeData->getVoxelNormalized(x, y, z, channel);
					float second = secondVolumeData->getVoxelNormalized(x, y, z, channel);

					first = firstRWM.normalizedToRealWorld(first);
					second = secondRWM.normalizedToRealWorld(second);

					float val = ((1.0f - alpha) * first) + (alpha * second);
					value[channel] = val;
				}

				outputVolume->voxel(pos) = value;
			}
		}
	}

	std::unique_ptr<Volume> volume { new Volume(outputVolume.release(), firstVolume) };
	volume->setRealWorldMapping(RealWorldMapping()); // Override to default rwm.
	volume->setModality(firstVolume->getModality());
	return volume;
}

void VolumeInterpolation::process() {
	auto* volumeList = this->inputVolumes_.getData();
	if (volumeList) {
		auto volume = interpolate(this->timeStep_.get(), volumeList);
		this->outputVolume_.setData(volume.release());
	}
}

void VolumeInterpolation::forceReload() {
	invalidate();
}

void VolumeInterpolation::updateTimeRange() {
	auto volumeList = this->inputVolumes_.getData();
	if (volumeList && !volumeList->empty()) {
		this->timeStep_.setMaxValue(static_cast<float>(volumeList->size() - 1));
	} else {
		this->timeStep_.setMaxValue(1.0f);
	}
}

}
