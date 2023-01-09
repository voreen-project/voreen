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

#include "streamsurfacecreator.h"

#include "integrators.h"
#include "modules/flowanalysis/datastructures/streamlinelist.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

#include <memory>
#include <random>

namespace voreen {

class RibbonStore {
public:
	struct RibbonBlueprint;

	RibbonStore(const VolumeRAMRepresentationLock& lock,
				const PortDataPointer<VolumeBase>& volume,
				VolumeRAM::Filter filterMode,
				InputSpace inputSpace,
				const std::vector<tgt::vec3>& seedPoints,
				StreamIntegratorInput input,
				int maxNumStreamlines,
				float ripThreshold,
				bool enableRipping,
				bool enableMerging,
				bool enableSplitting);

	const std::vector<Streamline>& getStreamLines() const;
	const std::vector<RibbonBlueprint>& getRibbonBlueprints() const;

	enum class EdgeCase {
		None,	///< No edge case occurred.
		Rip,	///< Rip up a ribbon.
		Split,	///< Split a ribbon into two by adding a new streamline in the middle.
		Merge,	///< Merge two ribbons into one.
	};

	struct RibbonInfo {
		size_t leftLine;	///< Index of the left streamline.
		size_t rightLine;	///< Index of the right streamline.
		size_t leftStart;	///< Index of the element on the left streamline with which we start the operation.
		size_t rightStart;	///< Index of the element on the right streamline with which we start the operation.
		size_t leftEnd;		///< Index one past the end on the left streamline.
		size_t rightEnd;	///< Index one past the end on the right streamline.
		float leftUCoord;	///< U-Coordinate of the left streamline.
		float rightUCoord;	///< U-Coordinate of the right streamline.
		float leftVCoord;	///< V-Coordinate of the left streamline.
		float rightVCoord;	///< V-Coordinate of the right streamline.
	};

	struct NoneCase {};

	struct RipCase {};

	struct SplitCase {
		size_t middleLine;		///< Index of the middle streamline.
		size_t l0;				///< Index of the l0 section on the left streamline.
		size_t r0;				///< Index of the r0 section on the right streamline.
		size_t l1;				///< Index of the l1 section on the left streamline.
		size_t r1;				///< Index of the r1 section on the right streamline.
		float l0VCoord;			///< V-Coordinate of the l0 section.
		float r0VCoord;			///< V-Coordinate of the r0 section.
		tgt::vec2 middleCoord;	///< TexCoord of the middle section.
	};

	struct MergeCase {
		size_t middleLine;		///< Index of the middle streamline.
		size_t l0;				///< Index of the l0 section on the left streamline.
		size_t r0;				///< Index of the r0 section on the right streamline.
		size_t l1;				///< Index of the l1 section on the left streamline.
		size_t r1;				///< Index of the r1 section on the right streamline.
		size_t m0;				///< Index of the m0 section on the middle streamline.
		tgt::vec2 middleCoord;	///< TexCoord of the middle section.
	};

	struct RibbonBuildEdgeCase {
		EdgeCase edgeCase;
		union {
			NoneCase none;
			RipCase rip;
			SplitCase split;
			MergeCase merge;
		} op;
	};

	struct RibbonBlueprint {
		RibbonInfo info;
		RibbonBuildEdgeCase edgeCase;
	};
private:
	std::vector<Streamline> lines_;
	std::vector<RibbonBlueprint> blueprints_;

	static const std::string loggerCat_;
};

const std::string StreamSurfaceCreator::loggerCat_("sciviscontest2022.StreamSurfaceCreator");
const std::string RibbonStore::loggerCat_("sciviscontest2022.StreamSurfaceCreator.RibbonStore");

StreamSurfaceCreator::StreamSurfaceCreator()
    : AsyncComputeProcessor()
    , volumeInport_(Port::INPORT, "volumeInport", "Volume Input")
	, seedingInport_(Port::INPORT, "seedingInport", "Seeding points Input")
	, geometryOutport_(Port::OUTPORT, "geometryOutport", "Geometry Output")
	, streamlineOutport_(Port::OUTPORT, "streamlineOutport", "Streamlines Output")
	, maxNumStreamlines_("maxNumStreamlines", "Maximum number of streamlines", 1, 1, std::numeric_limits<int>::max())
	, streamlineLengthThreshold_("streamlineLengthThreshold", "Restrict streamline length", 1000, 2, 10000)
	, absoluteMagnitudeThreshold_("absoluteMagnitudeThreshold", "Threshold of Magnitude (absolute)", tgt::vec2(0.0f, 1000.0f), 0.0f, 9999.99f)
    , stepSize_("stepSize", "Integration stepSize", 0.25f, 0.001f, 100.0f)
	, fitAbsoluteMagnitudeThreshold_("fitAbsoluteMagnitude", "Fit absolute Threshold to Input", false)
	, rippingThresholdAngle_("rippingThresholdAngle", "Angle at which a ribbon rips apart", 120, 0, 360)
	, filterMode_("filterModeProp", "Filtering:", Processor::INVALID_RESULT, false, Property::LOD_DEVELOPMENT)
	, inputSpace_("inputSpaceProp", "Input Space:", Processor::INVALID_RESULT, false, Property::LOD_DEVELOPMENT)
	, seedPointsOffset_("seedPointsOffset", "Seed point offset:", tgt::vec3::zero, tgt::vec3(-9999.9f), tgt::vec3(9999.9f), Processor::INVALID_RESULT, NumericProperty<tgt::vec3>::STATIC, Property::LOD_DEVELOPMENT)
	, linesOutportEnabled_("linesOutportEnabled", "Output Streamlines", false, Processor::INVALID_RESULT, Property::LOD_DEBUG)
	, rippingEnabled_("rippingEnabled", "Handle ripping edge case", true, Processor::INVALID_RESULT, Property::LOD_DEBUG)
	, mergingEnabled_("mergingEnabled", "Handle merging edge case", true, Processor::INVALID_RESULT, Property::LOD_DEBUG)
	, splittingEnabled_("splittingEnabled", "Handle splitting edge case", true, Processor::INVALID_RESULT, Property::LOD_DEBUG)
	, velocityUnitConversion_("velocityUnitConversion", "Input Velocity Unit")
{
	AsyncComputeProcessor::addPort(this->volumeInport_);
	AsyncComputeProcessor::addPort(this->seedingInport_);
	AsyncComputeProcessor::addPort(this->geometryOutport_);
	this->streamlineOutport_.deinitialize();

	AsyncComputeProcessor::addProperty(maxNumStreamlines_);
	maxNumStreamlines_.setTracking(false);
	maxNumStreamlines_.setGroupID("stream surface");

	AsyncComputeProcessor::addProperty(streamlineLengthThreshold_);
	streamlineLengthThreshold_.setTracking(false);
	streamlineLengthThreshold_.setGroupID("stream surface");

	AsyncComputeProcessor::addProperty(absoluteMagnitudeThreshold_);
	absoluteMagnitudeThreshold_.setTracking(false);
	absoluteMagnitudeThreshold_.setNumDecimals(2);
	absoluteMagnitudeThreshold_.setGroupID("stream surface");

	AsyncComputeProcessor::addProperty(stepSize_);
	stepSize_.setGroupID("stream surface");

	AsyncComputeProcessor::addProperty(fitAbsoluteMagnitudeThreshold_);
	ON_CHANGE(fitAbsoluteMagnitudeThreshold_, StreamSurfaceCreator, adjustPropertiesToInput);
	fitAbsoluteMagnitudeThreshold_.setGroupID("stream surface");

	AsyncComputeProcessor::addProperty(rippingThresholdAngle_);
	rippingThresholdAngle_.setTracking(false);
	rippingThresholdAngle_.setGroupID("stream surface");

	AsyncComputeProcessor::addProperty(filterMode_);
	filterMode_.addOption("linear", "Linear", VolumeRAM::LINEAR);
	filterMode_.addOption("nearest", "Nearest", VolumeRAM::NEAREST);
	filterMode_.setGroupID("stream surface");

	AsyncComputeProcessor::addProperty(inputSpace_);
	inputSpace_.addOption("voxel", "Voxel space", InputSpace::Voxel);
	inputSpace_.addOption("world", "World space", InputSpace::World);
	inputSpace_.setGroupID("stream surface");

	AsyncComputeProcessor::addProperty(seedPointsOffset_);
	seedPointsOffset_.setGroupID("stream surface");

	AsyncComputeProcessor::addProperty(this->linesOutportEnabled_);
	ON_CHANGE(this->linesOutportEnabled_, StreamSurfaceCreator, linesOutportEnabledChange);
	linesOutportEnabled_.setGroupID("stream surface");

	AsyncComputeProcessor::addProperty(this->rippingEnabled_);
	rippingEnabled_.setGroupID("stream surface");

	AsyncComputeProcessor::addProperty(this->mergingEnabled_);
	mergingEnabled_.setGroupID("stream surface");

	AsyncComputeProcessor::addProperty(this->splittingEnabled_);
	splittingEnabled_.setGroupID("stream surface");

	AsyncComputeProcessor::addProperty(this->velocityUnitConversion_);
	velocityUnitConversion_.addOption("km/s", "km/s", 1000000.0f);
	velocityUnitConversion_.addOption("m/s", "m/s", 1000.0f);
	velocityUnitConversion_.addOption("cm/s", "cm/s", 10.0f);
	velocityUnitConversion_.addOption("mm/s", "mm/s", 1.0f);
	velocityUnitConversion_.set("m/s");

	setPropertyGroupGuiName("stream surface", "Stream Surface Settings");
}

StreamSurfaceCreatorInput StreamSurfaceCreator::prepareComputeInput() {
	auto flowVolume = volumeInport_.getThreadSafeData();
	if(!flowVolume) {
		throw InvalidInputException("No volume", InvalidInputException::S_ERROR);
	}

	tgt::mat4 worldToVoxelMatrix = flowVolume->getWorldToVoxelMatrix();
	tgt::Bounds roi = flowVolume->getBoundingBox().getBoundingBox();

	auto inGeometry = dynamic_cast<const PointSegmentListGeometryVec3*>(seedingInport_.getData());
	if (!inGeometry) {
		throw InvalidInputException("Invalid Geometry type", InvalidInputException::S_ERROR);
	}

	std::vector<std::vector<tgt::vec3>> segments;
	segments.reserve(inGeometry->getNumSegments());
	for (auto i = 0u; i < inGeometry->getNumSegments(); ++i) {
		segments.emplace_back(inGeometry->getSegment(i));
	}

	auto offset = this->seedPointsOffset_.get();

	for (auto& segment: segments) {
		for (auto& point: segment) {
			point += offset;
		}
	}

	std::unique_ptr<GeometrySequence> geometry(new GeometrySequence());
	std::unique_ptr<StreamlineListBase> streamlines {};
	if (this->linesOutportEnabled_.get()) {
		streamlines = std::unique_ptr<StreamlineListBase>(new StreamlineList());
	}

	return StreamSurfaceCreatorInput {
		this->rippingEnabled_.get(),
		this->mergingEnabled_.get(),
		this->splittingEnabled_.get(),
		this->maxNumStreamlines_.get(),
		this->rippingThresholdAngle_.get(),
		this->streamlineLengthThreshold_.get(),
        this->stepSize_.get(),
		this->velocityUnitConversion_.getValue(),
		this->absoluteMagnitudeThreshold_.get(),
		this->filterMode_.getValue(),
		this->inputSpace_.getValue(),
		volumeInport_.getThreadSafeData(),
		std::move(segments),
		std::move(geometry),
		std::move(streamlines)
	};
}

StreamSurfaceCreatorOutput StreamSurfaceCreator::compute(StreamSurfaceCreatorInput input, ProgressReporter &progressReporter) const {
	// The algorithm is conceptually pretty simple and composed of three steps:
	// 1. Generate the required streamlines.
	// 2. Analyze the computed streamlines for the three edge cases of
	// 		adding/removing particles and ripping of ribbons, like described in
	//		the original paper.
	// 3. Construct the ribbons (surface).

	auto ripThreshold = static_cast<float>(input.ripThreshold);
	PortDataPointer<VolumeBase> flowVolume = std::move(input.flowVolume);
	VolumeRAMRepresentationLock representation(flowVolume);
	float stepSize = input.stepSize;
	size_t upperLengthThreshold = input.streamlineLengthThreshold;
	tgt::vec3 dimensions = tgt::vec3(representation->getDimensions() - tgt::svec3::one);
	auto voxelToWorldMatrix = flowVolume->getVoxelToWorldMatrix();

	StreamIntegratorInput integratorInput {
		voxelToWorldMatrix * tgt::vec3::zero,
		voxelToWorldMatrix * dimensions,
		stepSize,
		upperLengthThreshold,
		input.absoluteMagnitudeThreshold,
		2 * tgt::PIf,	// Ignore angle
		input.velocityUnitConversion
	};

	std::unique_ptr<StreamlineListBase> streamlines = std::move(input.streamlines);
	std::unique_ptr<GeometrySequence> geometry = std::move(input.geometry);

	// Initialize the vector of sub-geometries for parallel processing.
	std::vector<std::unique_ptr<TriangleMeshGeometryTexCoord>> subGeometries{};
	std::vector<std::vector<Streamline>> subStreamlines{};
	subGeometries.reserve(input.seedPoints.size());
	subStreamlines.reserve(input.seedPoints.size());
	for (auto i = 0u; i < input.seedPoints.size(); ++i) {
		subGeometries.emplace_back(new TriangleMeshGeometryTexCoord());
		subStreamlines.emplace_back();
	}

	ThreadedTaskProgressReporter progress(progressReporter, input.seedPoints.size());
	bool aborted = false;

	// Compute each geometry individually.
#ifdef VRN_MODULE_OPENMP
	#pragma omp parallel for default(none)\
	shared(aborted, input, representation, flowVolume, integratorInput,\
	progress, ripThreshold, subGeometries, subStreamlines)
	for (long i=0; i<static_cast<long>(input.seedPoints.size()); i++) {
		if (aborted) {
			continue;
		}
#else
	for (auto i = 0u; i < input.seedPoints.size(); ++i) {
#endif // VRN_MODULE_OPENMP

		auto& subGeometry = subGeometries[i];
		auto& subStreamline = subStreamlines[i];
		const auto& seedPoints = input.seedPoints[i];

		// The first and second step are implemented in the `RibbonStore`.
		RibbonStore store(representation,
						flowVolume,
						input.filterMode,
						input.inputSpace,
						seedPoints,
						integratorInput,
						input.maxNumStreamlines,
						ripThreshold,
						input.enableRipping,
						input.enableMerging,
						input.enableSplitting);

		// 3. Construct the ribbons by adding triangles along the length of neighboring
		// streamlines.
		auto& lines = store.getStreamLines();
		auto& ribbons = store.getRibbonBlueprints();
		for(auto& ribbon: ribbons) {
			auto info = ribbon.info;
			auto edgeCase = ribbon.edgeCase;

			auto& leftLine = lines[info.leftLine];
			auto& rightLine = lines[info.rightLine];

			size_t l0 = info.leftStart;
			size_t r0 = info.rightStart;
			bool caughtUp = (info.leftEnd == info.leftStart || info.rightEnd == info.rightStart)
				|| (l0 >= info.leftEnd && r0 >= info.rightEnd);

			float leftVCoord = info.leftVCoord;
			float rightVCoord = info.rightVCoord;

			// Construct the ribbon.
			while (!caughtUp) {
				size_t l1 = l0 + 1;
				size_t r1 = r0 + 1;

				auto& l0Elem = leftLine.getElementAt(l0);
				auto& r0Elem = rightLine.getElementAt(r0);

				// Add triangle the fist two vertices of the
				// triangle (L0, R0, _) counter clockwise.
				Triangle<VertexTexCoord> face;
				face.v_[0] = VertexTexCoord(l0Elem.position_, { info.leftUCoord, leftVCoord }, tgt::ivec2::zero);
				face.v_[1] = VertexTexCoord(r0Elem.position_, { info.rightUCoord, rightVCoord }, tgt::ivec2::zero);

				if (l1 >= info.leftEnd && r1 >= info.rightEnd) {
					// We are caught up on both sides, interrupt.
					break;
				} else if (l1 >= info.leftEnd) {
					// We are caught up on the left side but need to catch up on the right.
					auto& r1Elem = rightLine.getElementAt(r1);

					// Add triangle (L0, R0, R1) counter clockwise.
					face.v_[2] = VertexTexCoord(r1Elem.position_, { info.rightUCoord, rightVCoord + 1 }, tgt::ivec2::zero);

					// Advance the right.
					r0++;
					rightVCoord++;
				} else if (r1 >= info.rightEnd) {
					// We are caught up on the right side but need to catch up on the left.
					auto& l1Elem = leftLine.getElementAt(l1);

					// Add triangle (L0, R0, L1) counter clockwise.
					face.v_[2] = VertexTexCoord(l1Elem.position_, { info.leftUCoord, leftVCoord + 1 }, tgt::ivec2::zero);

					// Advance the left.
					l0++;
					leftVCoord++;
				} else {
					// We are not caught up, so we need to determine on which side to place
					// the last vertex. This is done using a greedy algorithm by comparing
					// the length of the lines (L0, R1) and (L1, R0) and selecting the
					// shorter one.
					auto& l1Elem = leftLine.getElementAt(l1);
					auto& r1Elem = rightLine.getElementAt(r1);

					auto l0r1Dist = tgt::distanceSq(l0Elem.position_, r1Elem.position_);
					auto l1r0Dist = tgt::distanceSq(l1Elem.position_, r0Elem.position_);
					auto minDist = std::min(l0r1Dist, l1r0Dist);

					auto advanceLeft = l1r0Dist == minDist;
					if (advanceLeft) {
						l0++;
						leftVCoord++;
						face.v_[2] = VertexTexCoord(l1Elem.position_, { info.leftUCoord, leftVCoord }, tgt::ivec2::zero);
					} else {
						r0++;
						rightVCoord++;
						face.v_[2] = VertexTexCoord(r1Elem.position_, { info.rightUCoord, rightVCoord }, tgt::ivec2::zero);
					}
				}

				subGeometry->addTriangle(face);
				caughtUp = (l0 == info.leftEnd - 1) && (r0 == info.rightEnd - 1);
			}

			// Handle the edge case at the end of the ribbon.
			switch(edgeCase.edgeCase) {
			case RibbonStore::EdgeCase::None:
			case RibbonStore::EdgeCase::Rip:
				break;
			case RibbonStore::EdgeCase::Split:
				{
					auto split = edgeCase.op.split;

					auto& middleLine = lines[split.middleLine];
					auto& l0Elem = leftLine.getElementAt(split.l0);
					auto& r0Elem = rightLine.getElementAt(split.r0);
					auto& l1Elem = leftLine.getElementAt(split.l1);
					auto& r1Elem = rightLine.getElementAt(split.r1);
					auto& middleElem = middleLine.getElementAt(0);

					// Add the triangle (L0, R0, Middle).
					Triangle<VertexTexCoord> face;
					face.v_[0] = VertexTexCoord(l0Elem.position_, { info.leftUCoord, leftVCoord }, tgt::ivec2::zero);
					face.v_[1] = VertexTexCoord(r0Elem.position_, { info.rightUCoord, rightVCoord }, tgt::ivec2::zero);
					face.v_[2] = VertexTexCoord(middleElem.position_, split.middleCoord, tgt::ivec2::zero);
					subGeometry->addTriangle(face);

					// Add the triangle (L0, Middle, L1).
					face.v_[0] = VertexTexCoord(l0Elem.position_, { info.leftUCoord, leftVCoord }, tgt::ivec2::zero);
					face.v_[1] = VertexTexCoord(middleElem.position_, split.middleCoord, tgt::ivec2::zero);
					face.v_[2] = VertexTexCoord(l1Elem.position_, { info.leftUCoord, leftVCoord + 1.0f }, tgt::ivec2::zero);
					subGeometry->addTriangle(face);

					// Add the triangle (Middle, R0, R1).
					face.v_[0] = VertexTexCoord(middleElem.position_, split.middleCoord, tgt::ivec2::zero);
					face.v_[1] = VertexTexCoord(r0Elem.position_, { info.rightUCoord, rightVCoord }, tgt::ivec2::zero);
					face.v_[2] = VertexTexCoord(r1Elem.position_, { info.rightUCoord, rightVCoord + 1.0f }, tgt::ivec2::zero);
					subGeometry->addTriangle(face);
				}
				break;
			case RibbonStore::EdgeCase::Merge:
			{
				auto merge = edgeCase.op.merge;

				auto& middleLine = lines[merge.middleLine];
				auto& l0Elem = leftLine.getElementAt(merge.l0);
				auto& l1Elem = leftLine.getElementAt(merge.l1);
				auto& r0Elem = rightLine.getElementAt(merge.r0);
				auto& r1Elem = rightLine.getElementAt(merge.r1);
				auto& m0Elem = middleLine.getElementAt(merge.m0);

				// Add the triangle (L0, M0, L1).
				Triangle<VertexTexCoord> face;
				face.v_[0] = VertexTexCoord(l0Elem.position_, { info.leftUCoord, leftVCoord }, tgt::ivec2::zero);
				face.v_[1] = VertexTexCoord(m0Elem.position_, merge.middleCoord, tgt::ivec2::zero);
				face.v_[2] = VertexTexCoord(l1Elem.position_, { info.leftUCoord, leftVCoord + 1.0f }, tgt::ivec2::zero);
				subGeometry->addTriangle(face);

				// Add the triangle (M0, R0, R1).
				face.v_[0] = VertexTexCoord(m0Elem.position_, merge.middleCoord, tgt::ivec2::zero);
				face.v_[1] = VertexTexCoord(r0Elem.position_, { info.rightUCoord, rightVCoord }, tgt::ivec2::zero);
				face.v_[2] = VertexTexCoord(r1Elem.position_, { info.rightUCoord, rightVCoord + 1.0f }, tgt::ivec2::zero);
				subGeometry->addTriangle(face);

				// Add the triangle (M0, R1, L1).
				face.v_[0] = VertexTexCoord(m0Elem.position_, merge.middleCoord, tgt::ivec2::zero);
				face.v_[1] = VertexTexCoord(r1Elem.position_, { info.rightUCoord, rightVCoord + 1.0f }, tgt::ivec2::zero);
				face.v_[2] = VertexTexCoord(l1Elem.position_, { info.leftUCoord, leftVCoord + 1.0f }, tgt::ivec2::zero);
				subGeometry->addTriangle(face);
			}
				break;
			}
		}

		subStreamline.clear();
		std::copy(lines.begin(), lines.end(), std::back_inserter(subStreamline));

		// Update the progress.
#ifdef VRN_MODULE_OPENMP
#pragma omp critical
		if (progress.reportStepDone()) {
			aborted = true;
		}
#else
		if (progress.reportStepDone()) {
			aborted = true;
            break;
		}
#endif // VRN_MODULE_OPENMP
	}

	if (aborted) {
		throw boost::thread_interrupted();
	}

	// Insert sub-geometries into the main geometry.
	for(auto& subGeometry: subGeometries) {
		geometry->addGeometry(subGeometry.release());
	}

	if (streamlines) {
		for (auto& lines: subStreamlines) {
			for (auto& line: lines) {
				streamlines->addStreamline(line);
			}
		}
	}

	return StreamSurfaceCreatorOutput {
		std::move(geometry),
		std::move(streamlines)
	};
}

void StreamSurfaceCreator::processComputeOutput(StreamSurfaceCreatorOutput output) {
	this->geometryOutport_.setData(output.geometry.release());

	if (output.streamlines) {
		this->streamlineOutport_.setData(output.streamlines.release());
	}
}

void StreamSurfaceCreator::adjustPropertiesToInput() {

	const VolumeBase* volume = volumeInport_.getData();
	auto inGeometry = dynamic_cast<const PointSegmentListGeometryVec3*>(seedingInport_.getData());

	if(fitAbsoluteMagnitudeThreshold_.get() && volume !=nullptr) {
		auto* data = volume->getDerivedData<VolumeMinMaxMagnitude>();

		absoluteMagnitudeThreshold_.setMinValue(data->getMinMagnitude());
		absoluteMagnitudeThreshold_.setMaxValue(data->getMaxMagnitude());
		absoluteMagnitudeThreshold_.set(tgt::vec2(data->getMinMagnitude(), data->getMaxMagnitude()));
	}
	else {
		absoluteMagnitudeThreshold_.setMinValue(0.0f);
		absoluteMagnitudeThreshold_.setMaxValue(5000.0f);
	}

	if (inGeometry == nullptr || inGeometry->getNumSegments() == 0) {
		maxNumStreamlines_.set(1);
	} else {
		auto& data = inGeometry->getData();
		auto max = std::max_element(data.begin(), data.end(),
									[](const std::vector<tgt::vec3>& l, const std::vector<tgt::vec3>& r) { return l.size() < r.size(); });
		maxNumStreamlines_.set((int)(*max).size() * 20);
	}
}

void StreamSurfaceCreator::linesOutportEnabledChange() {
	auto enabled = this->linesOutportEnabled_.get();
	if (enabled) {
		this->streamlineOutport_.initialize();
		AsyncComputeProcessor::addPort(this->streamlineOutport_);
	} else {
		this->streamlineOutport_.deinitialize();
		AsyncComputeProcessor::removePort(this->streamlineOutport_);
	}
}

RibbonStore::RibbonStore(const VolumeRAMRepresentationLock& lock,
						 const PortDataPointer<VolumeBase>& volume,
						 VolumeRAM::Filter filterMode,
						 InputSpace inputSpace,
						 const std::vector<tgt::vec3>& seedPoints,
						 StreamIntegratorInput input,
						 int maxNumStreamlines,
						 float ripThreshold,
						 bool enableRipping,
						 bool enableMerging,
						 bool enableSplitting)
	: blueprints_()
	, lines_()
{
	tgtAssert(ripThreshold >= 0 && ripThreshold <= 360, "The angle must be in the range [0, 360]")
	constexpr auto tau = 6.283185307179586f;
	auto ripThresholdRadians = ripThreshold * tau / 360.0f;

	// 1. Like described above, we first generate the required streamlines
	// using a streamline integrator, which integrates only in the direction
	// of the flow, like `StreamIntegrator`.
	StreamIntegrator integrator {};
	auto spaceToWorldMatrix = tgt::mat4::identity;
	SpatialSampler sampler(*lock, volume->getRealWorldMapping(), filterMode, volume->getWorldToVoxelMatrix());
	std::vector<Streamline> streamlines(seedPoints.size());
	streamlines.resize(seedPoints.size(), {});

	switch(inputSpace) {
	case InputSpace::Voxel:
		spaceToWorldMatrix = volume->getVoxelToWorldMatrix();
		break;
	case InputSpace::World:
		break;
	}

	// Use the distance of a seed point on the curve as the U-Coordinate.
	float seedCurveLength = 0.0f;
	std::vector<float> seedCurveUVCoords{};
	seedCurveUVCoords.reserve(seedPoints.size());
	seedCurveUVCoords.push_back(0.0f);
	for (auto it = seedPoints.begin() + 1; it != seedPoints.end(); ++it) {
		auto seed = spaceToWorldMatrix * *it;
		auto lastSeed = spaceToWorldMatrix * *(it - 1);
		seedCurveLength += tgt::distance(lastSeed, seed);
		seedCurveUVCoords.push_back(seedCurveLength);
	}

	// Normalize the U-Coordinate from [0, seedCurveLength] to [0, 1].
	for(auto& coord: seedCurveUVCoords) {
		coord /= seedCurveLength;
	}

	// Compute initial streamlines.
#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for default(none) shared(seedPoints, streamlines, integrator, sampler, input, spaceToWorldMatrix)
#endif // VRN_MODULE_OPENMP
	for (auto i = 0u; i < streamlines.size(); ++i) {
		auto& seed = seedPoints[i];
		auto seedWorld = spaceToWorldMatrix * seed;
		streamlines[i] = integrator.integrate(sampler, input, seedWorld);
	}

	std::vector<RibbonBlueprint> blueprints;
	blueprints.resize(streamlines.size(), {});

	// Early exit if we have less than two lines.
	if (streamlines.size() <= 1) {
		this->blueprints_ = std::move(blueprints);
		this->lines_ = std::move(streamlines);
		return;
	}

	// 2. Analyze the resulting streamlines to identify the three edge cases described
	// in the paper. This is done by grouping neighboring pairs of streamlines into
	// ribbons, which are then inserted into a processing queue. Each ribbon can then
	// be analyzed independently of the other ribbons.
	struct RibbonInfo {
		size_t leftIdx;			///< Index of the left line.
		size_t rightIdx;		///< Index of the right line.
		size_t leftStart;		///< Start index on the left line.
		size_t rightStart;		///< Start index on the right line.
		size_t leftLength;		///< Maximum length of the left line.
		size_t rightLength;		///< Maximum length of the right line.
		float leftUCoord;		///< U-Coordinate of the left line.
		float rightUCoord;		///< U-Coordinate of the right line.
		float leftVCoord;		///< V-Coordinate of the left line.
		float rightVCoord;		///< V-Coordinate of the right line.
	};

	std::vector<RibbonInfo> processingQueue;
	std::vector<Streamline> tmpStreamlines {};
	std::vector<RibbonInfo> tmpProcessingQueue {};
	processingQueue.reserve(streamlines.size() - 1);

	// Add the initial ribbons on the processing queue.
	for (size_t left = 0; left < streamlines.size() - 1; ++left) {
		size_t rightIdx = left + 1;
		processingQueue.push_back(RibbonInfo {
			left,
			rightIdx,
			0,
			0,
			std::min(streamlines[left].getNumElements(), input.upperLengthThreshold),
			std::min(streamlines[rightIdx].getNumElements(), input.upperLengthThreshold),
			seedCurveUVCoords[left],
			seedCurveUVCoords[rightIdx],
			0.0f,
			0.0f,
		});
	}

	// Work through the processing queue.
	bool logWarning = true;
	bool warningTriggered = false;
	size_t blueprintOffset = 0u;
	while (!processingQueue.empty()) {
		blueprints.resize(blueprintOffset + processingQueue.size(), {});

#ifdef VRN_MODULE_OPENMP
		#pragma omp parallel for default(none)\
		shared(blueprints, processingQueue, tmpProcessingQueue, streamlines, tmpStreamlines, input, integrator, sampler, warningTriggered)\
		firstprivate(enableRipping, ripThresholdRadians, enableSplitting, maxNumStreamlines, blueprintOffset)
#endif // VRN_MODULE_OPENMP
		for (auto i = 0u; i < processingQueue.size(); ++i) {
			RibbonInfo ribbon = processingQueue[i];

			auto& leftLine = streamlines[ribbon.leftIdx];
			auto& rightLine = streamlines[ribbon.rightIdx];

			size_t l0 = ribbon.leftStart;
			size_t r0 = ribbon.rightStart;
			size_t lEnd = l0 + ribbon.leftLength;
			size_t rEnd = r0 + ribbon.rightLength;
			float lVCoord = ribbon.leftVCoord;
			float rVCoord = ribbon.rightVCoord;
			bool caughtUp = (l0 >= lEnd) && (r0 >= rEnd);

			RibbonBlueprint blueprint {};
			blueprint.edgeCase = RibbonBuildEdgeCase { EdgeCase::None, { NoneCase {} } };

			// Traverse the entire ribbon using the greedy algorithm and check for
			// the presence of the edge cases at each step.
			while (!caughtUp) {
				size_t l1 = l0 + 1;
				size_t r1 = r0 + 1;

				// If we are caught up on either side we assert that no edge case can show up.
				if (l1 >= lEnd) {
					r0++;
					rVCoord++;
				} else if (r1 >= rEnd) {
					l0++;
					lVCoord++;
				} else {
					auto& l0Elem = leftLine.getElementAt(l0);
					auto& l1Elem = leftLine.getElementAt(l1);
					auto& r0Elem = rightLine.getElementAt(r0);
					auto& r1Elem = rightLine.getElementAt(r1);

					auto velocityMagnitude = tgt::length(l0Elem.velocity_) * tgt::length(r0Elem.velocity_);
					auto velocityCos = tgt::dot(l0Elem.velocity_, r0Elem.velocity_) / velocityMagnitude;
					auto velocityAngle = std::acos(velocityCos);

					// Case 1: Rip up a ribbon.
					if (enableRipping && velocityAngle >= ripThresholdRadians) {
						// If the two streamlines are moving in opposing directions we rip up the
						// ribbon. The original paper does not clearly specify how the "ripping"
						// is intended to be implemented, therefore we take the easy route and
						// simply stop the processing of the current ribbon.
						blueprint.edgeCase.edgeCase = EdgeCase::Rip;
						blueprint.edgeCase.op.rip = RipCase {};

						lEnd = l0;
						rEnd = r0;
						break;
					}

					auto l0r0Dist = tgt::distanceSq(l0Elem.position_, r0Elem.position_);
					auto l1r1Dist = tgt::distanceSq(l1Elem.position_, r1Elem.position_);
					auto l0l1Dist = tgt::distanceSq(l0Elem.position_, l1Elem.position_);
					auto r0r1Dist = tgt::distanceSq(r0Elem.position_, r1Elem.position_);
					auto maxWidth = std::max(l0r0Dist, l1r1Dist);
					auto maxHeight = std::max(l0l1Dist, r0r1Dist);

					// Case 2: Splitting a ribbon in two.
					bool stop = false;
					if (enableSplitting && maxWidth / maxHeight >= 4.0f) {

#ifdef VRN_MODULE_OPENMP
						#pragma omp critical
#endif // VRN_MODULE_OPENMP
						if (streamlines.size() + tmpStreamlines.size() >= maxNumStreamlines) {
							warningTriggered = true;
						} else {
							// We detected that we need to split up the ribbon. This is implemented
							// by computing a new streamline in the middle of the existing two streamlines
							// and then adding the two new ribbons to the processing queue.
							auto leftLength = lEnd - l1;
							auto rightLength = rEnd - r1;
							auto middleLength = std::max(leftLength, rightLength);
							auto lineInput = input;
							lineInput.upperLengthThreshold = middleLength;

							auto middlePoint = (l1Elem.position_ + r1Elem.position_)/2.0f;
							auto line = integrator.integrate(sampler, lineInput, middlePoint);

							if (line.getNumElements()!=0 && !(l1 >= lEnd - 1 && r1 >= rEnd - 1)) {
								// We were able to compute the new streamline; all that is left to do
								// is record the edge case, add the new ribbons to the processing queue
								// and stop the processing of the current ribbon.
								size_t middleLineIdx = streamlines.size() + tmpStreamlines.size();
								auto middleUCoord = (ribbon.leftUCoord + ribbon.rightUCoord) / 2.0f;
								auto middleVCoord = (lVCoord + rVCoord + 2.0f) / 2.0f;

								blueprint.edgeCase.edgeCase = EdgeCase::Split;
								blueprint.edgeCase.op.split = SplitCase{
									middleLineIdx,
									l0,
									r0,
									l1,
									r1,
									lVCoord,
									rVCoord,
									{ middleUCoord, middleVCoord },
								};

								middleLength = line.getNumElements();
								tmpStreamlines.push_back(std::move(line));
								tmpProcessingQueue.push_back(RibbonInfo{
									ribbon.leftIdx,
									middleLineIdx,
									l1,
									0,
									leftLength,
									middleLength,
									ribbon.leftUCoord,
									middleUCoord,
									lVCoord,
									middleVCoord
								});
								tmpProcessingQueue.push_back(RibbonInfo{
									middleLineIdx,
									ribbon.rightIdx,
									0,
									r1,
									middleLength,
									rightLength,
									middleUCoord,
									ribbon.rightUCoord,
									middleVCoord,
									rVCoord,
								});

								lEnd = l1;
								rEnd = r1;
								stop = true;
							}
						}
					}

					if (stop) {
						break;
					}

					/// Todo: Implement merging

					// Greedy algorithm for determining on which side to advance by selecting
					// the side with the minimum diagonal length.
					auto l0r1Dist = tgt::distanceSq(l0Elem.position_, r1Elem.position_);
					auto l1r0Dist = tgt::distanceSq(l1Elem.position_, r0Elem.position_);
					auto minDist = std::min(l0r1Dist, l1r0Dist);
					auto advanceLeft = l1r0Dist == minDist;
					if (advanceLeft) {
						l0++;
						lVCoord++;
					} else {
						r0++;
						rVCoord++;
					}
				}

				caughtUp = (l0 >= lEnd || l0 >= lEnd - 1) && (r0 >= rEnd || r0 >= rEnd - 1);
			}

			// At this point we processed the ribbon, and we must therefore record
			// how far we have to follow each streamline.
			blueprint.info = RibbonStore::RibbonInfo {
				ribbon.leftIdx,
				ribbon.rightIdx,
				ribbon.leftStart,
				ribbon.rightStart,
				lEnd,
				rEnd,
				ribbon.leftUCoord,
				ribbon.rightUCoord,
				ribbon.leftVCoord,
				ribbon.rightVCoord,
			};
			blueprints[blueprintOffset + i] = blueprint;
		}

		if (warningTriggered && logWarning) {
			logWarning = false;
			LWARNING(
				"Exceeded the maximum number of streamlines. Restricting to " << maxNumStreamlines);
		}

		blueprintOffset += processingQueue.size();
		processingQueue = tmpProcessingQueue;
		tmpProcessingQueue.clear();

		streamlines.insert(streamlines.end(), tmpStreamlines.begin(), tmpStreamlines.end());
		tmpStreamlines.clear();
	}

	this->lines_ = std::move(streamlines);
	this->blueprints_ = std::move(blueprints);
}

const std::vector<Streamline>& RibbonStore::getStreamLines() const {
	return this->lines_;
}

const std::vector<RibbonStore::RibbonBlueprint>& RibbonStore::getRibbonBlueprints() const {
	return this->blueprints_;
}

}

