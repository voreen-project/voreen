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

#include "streamsurfaceviewer.h"

#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"

namespace voreen {

template <typename T>
void populateColors(VolumeRAM::Filter filter,
					BlendFunction blendFunction,
					const tgt::vec4& baseColor,
					const std::vector<PortDataPointer<VolumeBase>>& volumes,
					const TriangleMeshGeometry<T>* geometry,
					TriangleMeshGeometry<T>* output,
					const std::vector<std::unique_ptr<TransFunc1DKeys>>& transferFunctions);

StreamSurfaceViewer::StreamSurfaceViewer()
	: AsyncComputeProcessor()
	, volumeInport_(Port::INPORT, "volumeInport", "Volume Input", true)
	, geometryInport_(Port::INPORT, "geometryInport", "Geometry input")
	, geometryOutport_(Port::OUTPORT, "geometryOutport", "Geometry output")
	, blendFunction_("blendFunction", "Blending:", Processor::INVALID_RESULT, false)
	, filterMode_("filterMode", "Filter:", Processor::INVALID_RESULT, false)
	, baseColor_("baseColor", "Base color:")
	, volumeSelector_("volumeSelector", "Volume:", Processor::INVALID_RESULT, false)
	, transferFunction_("transferFunction_", "Color Map:")
	, volumeSettings_()
{
	AsyncComputeProcessor::addPort(this->volumeInport_);
	AsyncComputeProcessor::addPort(this->geometryInport_);
	AsyncComputeProcessor::addPort(this->geometryOutport_);
	ON_CHANGE(this->volumeInport_, StreamSurfaceViewer, onVolumeChange);

	AsyncComputeProcessor::addProperty(this->blendFunction_);
	this->blendFunction_.addOption("min", "Minimum", BlendFunction::MIN);
	this->blendFunction_.addOption("max", "Maximum", BlendFunction::MAX);
	this->blendFunction_.addOption("add", "Add", BlendFunction::ADD);
	this->blendFunction_.addOption("sub", "Subtract", BlendFunction::SUBTRACT);
	this->blendFunction_.addOption("mul", "Multiply", BlendFunction::MULTIPLY);
	this->blendFunction_.addOption("avg", "Average", BlendFunction::AVERAGE);
	this->blendFunction_.setGroupID("surface");

	AsyncComputeProcessor::addProperty(this->filterMode_);
	this->filterMode_.addOption("nearest", "Nearest", VolumeRAM::Filter::NEAREST);
	this->filterMode_.addOption("linear", "Linear", VolumeRAM::Filter::LINEAR);
	this->filterMode_.addOption("cubic", "Cubic", VolumeRAM::Filter::CUBIC);
	this->filterMode_.setGroupID("surface");

	AsyncComputeProcessor::addProperty(this->baseColor_);
	this->baseColor_.setGroupID("surface");

	AsyncComputeProcessor::addProperty(this->volumeSelector_);
	ON_PROPERTY_CHANGE(this->volumeSelector_, StreamSurfaceViewer, onVolumeSelected);
	this->volumeSelector_.setGroupID("surface");

	AsyncComputeProcessor::addProperty(this->transferFunction_);
	ON_PROPERTY_CHANGE(this->transferFunction_, StreamSurfaceViewer, onVolumeSettingsChange);
	this->transferFunction_.setGroupID("surface");
}

StreamSurfaceViewerInput StreamSurfaceViewer::prepareComputeInput() {
	auto flowVolumes = this->volumeInport_.getThreadSafeAllData();
	if(flowVolumes.empty()) {
		throw InvalidInputException("No volume", InvalidInputException::S_ERROR);
	}

	auto geometry = this->geometryInport_.getThreadSafeData();
	auto geometryPtr = dynamic_cast<const TriangleMeshGeometryBase*>(geometry.get());

	if (geometryPtr == nullptr) {
		throw InvalidInputException("Geometry does not inherit from TriangleMeshGeometryBase", InvalidInputException::S_ERROR);
	}

	std::unique_ptr<Geometry> output;
	switch (geometryPtr->getVertexLayout()) {
	case TriangleMeshGeometryBase::COLOR:
		output = std::unique_ptr<Geometry>{new TriangleMeshGeometryColor()};
		break;
	case TriangleMeshGeometryBase::COLOR_NORMAL:
		output = std::unique_ptr<Geometry>{new TriangleMeshGeometryColorNormal()};
		break;
	case TriangleMeshGeometryBase::COLOR_TEXCOORD:
		output = std::unique_ptr<Geometry>{new TriangleMeshGeometryColorTexCoord()};
		break;
	case TriangleMeshGeometryBase::COLOR_NORMAL_TEXCOORD:
		output = std::unique_ptr<Geometry>{new TriangleMeshGeometryColorNormalTexCoord()};
		break;
	case TriangleMeshGeometryBase::SIMPLE:
	case TriangleMeshGeometryBase::NORMAL:
	case TriangleMeshGeometryBase::TEXCOORD:
	case TriangleMeshGeometryBase::NORMAL_TEXCOORD:
		throw InvalidInputException("Geometry has no color component", InvalidInputException::S_ERROR);
	}

	std::vector<std::unique_ptr<TransFunc1DKeys>> transferFunctions {};
	for (auto& volume : flowVolumes) {
		auto modality = volume->getModality().getName();
		auto& settings = this->volumeSettings_.at(modality);
		transferFunctions.emplace_back(settings.transFuncKeys.clone());
	}

	return StreamSurfaceViewerInput {
		this->filterMode_.getValue(),
		this->blendFunction_.getValue(),
		this->baseColor_.get(),
		std::move(output),
		std::move(geometry),
		std::move(flowVolumes),
		std::move(transferFunctions)
	};
}

StreamSurfaceViewerOutput StreamSurfaceViewer::compute(StreamSurfaceViewerInput input, ProgressReporter &progressReporter) const {
	auto output = std::move(input.output);
	auto geometry = static_cast<const TriangleMeshGeometryBase*>(input.geometry.get());
	auto filter = input.filter;
	auto blendFunc = input.blendFunction;
	auto& baseColor = input.baseColor;

	switch (geometry->getVertexLayout()) {
	case TriangleMeshGeometryBase::COLOR: {
		auto geometryPtr = static_cast<const TriangleMeshGeometry<VertexColor>*>(geometry);
		auto outputPtr = static_cast<TriangleMeshGeometry<VertexColor>*>(output.get());
		populateColors(filter, blendFunc, baseColor,  input.volumes, geometryPtr, outputPtr, input.transferFunctions);
		break;
	}
	case TriangleMeshGeometryBase::COLOR_NORMAL: {
		auto geometryPtr = static_cast<const TriangleMeshGeometry<VertexColorNormal>*>(geometry);
		auto outputPtr = static_cast<TriangleMeshGeometry<VertexColorNormal>*>(output.get());
		populateColors(filter, blendFunc, baseColor,  input.volumes, geometryPtr, outputPtr, input.transferFunctions);
		break;
	}
	case TriangleMeshGeometryBase::COLOR_TEXCOORD: {
		auto geometryPtr = static_cast<const TriangleMeshGeometry<VertexColorTexCoord>*>(geometry);
		auto outputPtr = static_cast<TriangleMeshGeometry<VertexColorTexCoord>*>(output.get());
		populateColors(filter, blendFunc, baseColor,  input.volumes, geometryPtr, outputPtr, input.transferFunctions);
		break;
	}
	case TriangleMeshGeometryBase::COLOR_NORMAL_TEXCOORD: {
		auto geometryPtr = static_cast<const TriangleMeshGeometry<VertexColorNormalTexCoord>*>(geometry);
		auto outputPtr = static_cast<TriangleMeshGeometry<VertexColorNormalTexCoord>*>(output.get());
		populateColors(filter, blendFunc, baseColor,  input.volumes, geometryPtr, outputPtr, input.transferFunctions);
		break;
	}
	case TriangleMeshGeometryBase::SIMPLE:
	case TriangleMeshGeometryBase::NORMAL:
	case TriangleMeshGeometryBase::TEXCOORD:
	case TriangleMeshGeometryBase::NORMAL_TEXCOORD:
		throw InvalidInputException("Geometry has no color component", InvalidInputException::S_ERROR);
	}

	return StreamSurfaceViewerOutput {
		std::move(output)
	};
}

void StreamSurfaceViewer::processComputeOutput(StreamSurfaceViewerOutput output) {
	this->geometryOutport_.setData(output.geometry.release());
}

void StreamSurfaceViewer::onVolumeChange() {
	// Remove settings of old volumes.
	this->volumeSelector_.setOptions({});
	this->volumeSelector_.reset();

	auto volumes = this->volumeInport_.getAllData();
	if (volumes.empty()) {
		this->volumeSettings_.clear();
		return;
	}

	for (auto it = this->volumeSettings_.begin(); it != this->volumeSettings_.end();) {
		bool found = false;
		for (auto& volume : volumes) {
			auto modality = volume->getModality().getName();
			if (it->first == modality) {
				found = true;
			}
		}

		if (!found && this->volumeSelector_.hasKey(it->first)) {
			this->volumeSelector_.removeOption(it->first);
			it = this->volumeSettings_.erase(it);
		} else {
			++it;
		}
	}

	// Add settings for new volumes.
	for (auto& volume : volumes) {
		auto modality = volume->getModality().getName();
		if (this->volumeSettings_.find(modality) == this->volumeSettings_.end()) {
			this->volumeSettings_.insert({modality, StreamSurfaceViewerVolumeSettings {
				volume,
				{}
			}});
		}
		this->volumeSelector_.addOption(modality, modality);
	}
}

void StreamSurfaceViewer::onVolumeSelected() {
	// Load the settings.
	auto& selected = this->volumeSelector_.get();
	auto setting_it = this->volumeSettings_.find(selected);
	if (setting_it == this->volumeSettings_.end()) {
		return;
	}

	auto& settings = setting_it->second;

	// Apply the loaded settings.
	this->transferFunction_.set1DKeys(settings.transFuncKeys.clone());
	this->transferFunction_.setVolume(settings.volume);
}

void StreamSurfaceViewer::onVolumeSettingsChange() {
	// Load the settings.
	auto& selected = this->volumeSelector_.get();
	auto setting_it = this->volumeSettings_.find(selected);
	if (setting_it == this->volumeSettings_.end()) {
		return;
	}
	auto& settings = setting_it->second;

	// Apply the new settings.
	settings.transFuncKeys.clearKeys();
	settings.transFuncKeys.setMemberValuesFrom(this->transferFunction_.get());
}

template <typename T>
void populateColors(VolumeRAM::Filter filter,
					BlendFunction blendFunction,
					const tgt::vec4& baseColor,
					const std::vector<PortDataPointer<VolumeBase>>& volumes,
					const TriangleMeshGeometry<T>* geometry,
					TriangleMeshGeometry<T>* output,
					const std::vector<std::unique_ptr<TransFunc1DKeys>>& transferFunctions) {
	// Copy mesh to output.
	output->addMesh(geometry);

	// Cache colors so that we don't need to continuously swap the loaded volume.
	std::vector<tgt::vec4> baseColors;
	baseColors.resize(volumes.size(), baseColor);

	std::vector<std::vector<tgt::vec4>> vertexColors {};
	vertexColors.resize(output->getNumTriangles() * 3, baseColors);

	// Compute the colors.
	for (size_t i = 0; i < volumes.size(); i++) {
		const auto& volume = volumes[i];
		VolumeRAMRepresentationLock representation {volume};
		RealWorldMapping rwm = volume->getRealWorldMapping();
		auto worldToVoxelMatrix = volume->getWorldToVoxelMatrix();
		auto& transferFunction = transferFunctions[i];
		auto keys = transferFunction->getKeys();

		for (size_t j = 0; j < output->getNumTriangles(); j++) {
			auto& triangle = output->getTriangle(j);
			for (size_t k = 0; k < 3; k++) {
				auto& vertex = triangle.v_[k];
				auto voxelPos = worldToVoxelMatrix * vertex.pos_;

				float sample;
				switch (filter) {
				case VolumeRAM::NEAREST:
					sample = representation->getVoxelNormalized(voxelPos);
					break;
				case VolumeRAM::LINEAR:
					sample = representation->getVoxelNormalizedLinear(voxelPos);
					break;
				case VolumeRAM::CUBIC:
					sample = representation->getVoxelNormalizedCubic(voxelPos);
					break;
				}

				sample = rwm.normalizedToRealWorld(sample);
				sample = transferFunction->realWorldToNormalized(sample);

				auto idx = (j * 3) + k;
				auto& vertexColor = vertexColors[idx][i];

				// Find the pair of keys which span the sampled intensity.
				auto lastKeyIt = std::find_if(keys.begin(), keys.end(),
											   [&](const TransFuncMappingKey* key) { return key->getIntensity() > sample; });

				auto alphaMode = transferFunction->getAlphaMode();
				auto transformAlpha = [&](tgt::vec4 color) {
					switch (alphaMode) {
					case TransFunc1DKeys::AlphaMode::TF_ZERO_ALPHA:
						return tgt::vec4{ color.xyz(), 0.0f };
					case TransFunc1DKeys::AlphaMode::TF_USE_ALPHA:
						return color;
					case TransFunc1DKeys::AlphaMode::TF_ONE_ALPHA:
						return tgt::vec4{ color.xyz(), 1.0f };
					}
				};

				// Interpolate the color using the transfer function and the sampled intensity.
				if (keys.empty()) {
					continue;
				} else if (lastKeyIt == keys.begin()) {
					vertexColor = transformAlpha(tgt::vec4 { keys[0]->getColorL() } / 255.0f);
					continue;
				} else if (lastKeyIt == keys.end()) {
					vertexColor = transformAlpha(tgt::vec4 { (*(lastKeyIt - 1))->getColorR() } / 255.0f);
					continue;
				}

				auto& firstKey = **(lastKeyIt - 1);
				auto& lastKey = **lastKeyIt;

				tgt::vec4 leftColor = transformAlpha(tgt::vec4 { firstKey.getColorR() } / 255.0f);
				tgt::vec4 rightColor = transformAlpha(tgt::vec4 { lastKey.getColorL() } / 255.0f);

				auto leftIntensity = firstKey.getIntensity();
				auto rightIntensity = lastKey.getIntensity();

				auto factor = (sample - leftIntensity) / (rightIntensity - leftIntensity);
				vertexColor = ((1.0f - factor) * leftColor) + (factor * rightColor);
			}
		}
	}

	// Apply the computed color to the vertices.
	for (size_t i = 0; i < vertexColors.size(); i++) {
		size_t triangleIdx = i / 3;
		size_t vertexIdx = i % 3;

		// Blend the colors using the requested blend function.
		tgt::vec4 color = baseColor;
		for (auto& c : vertexColors[i]) {
			switch (blendFunction) {
			case MIN:
				color = tgt::min(color, c);
				break;
			case MAX:
				color = tgt::max(color, c);
				break;
			case ADD:
				color += c;
				break;
			case SUBTRACT:
				color -= c;
				break;
			case MULTIPLY:
				color *= c;
				break;
			case AVERAGE:
				color += c;
				break;
			}
		}

		switch (blendFunction) {
		case MIN:
		case MAX:
		case MULTIPLY:
			break;
		case ADD:
		case SUBTRACT:
			color = tgt::clamp(color, tgt::vec4::zero, tgt::vec4::one);
			break;
		case AVERAGE:
			color /= static_cast<float>(vertexColors[i].size());
			break;
		}

		auto triangle = output->getTriangle(triangleIdx);
		triangle.v_[vertexIdx].setColor(color);
		output->setTriangle(triangle, triangleIdx);
	}
}

}
