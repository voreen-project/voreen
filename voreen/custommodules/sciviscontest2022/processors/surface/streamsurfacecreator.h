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

#ifndef VRN_STREAMSURFACECREATOR_H
#define VRN_STREAMSURFACECREATOR_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/ports/geometryport.h"
#include "voreen/core/ports/volumeport.h"

#include "modules/flowanalysis/ports/streamlinelistport.h"
#include "voreen/core/datastructures/geometry/geometrysequence.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/vectorproperty.h"

namespace voreen {

enum InputSpace {
	Voxel,
	World
};

struct StreamSurfaceCreatorInput {
	bool enableRipping;
	bool enableMerging;
	bool enableSplitting;
	int maxNumStreamlines;
	int ripThreshold;
	int streamlineLengthThreshold;
    float stepSize;
	float velocityUnitConversion;
	tgt::vec2 absoluteMagnitudeThreshold;
	VolumeRAM::Filter filterMode;
	InputSpace inputSpace;
	PortDataPointer<VolumeBase> flowVolume;
	std::vector<std::vector<tgt::vec3>> seedPoints;
	std::unique_ptr<GeometrySequence> geometry;
	std::unique_ptr<StreamlineListBase> streamlines;
};

struct StreamSurfaceCreatorOutput {
	std::unique_ptr<GeometrySequence> geometry;
	std::unique_ptr<StreamlineListBase> streamlines;
};

class VRN_CORE_API StreamSurfaceCreator: public AsyncComputeProcessor<StreamSurfaceCreatorInput, StreamSurfaceCreatorOutput> {
public:
    StreamSurfaceCreator();

    Processor* create() const override { return new StreamSurfaceCreator(); }
	bool isReady() const override {
		return (geometryOutport_.isConnected() || streamlineOutport_.isConnected())
				&& volumeInport_.isConnected()
				&& seedingInport_.isConnected();
	}

    std::string getCategory() const override { return "Stream surface Processing"; }
    std::string getClassName() const override { return "StreamSurfaceCreator"; }
    Processor::CodeState getCodeState() const override { return CODE_STATE_EXPERIMENTAL; }

protected:
    void setDescriptions() override {
		setDescription("This processor is used to create stream surfaces from a vec3 volume. The resulting streamlines can be visualized or modified "
					   "by other processors.");
		this->maxNumStreamlines_.setDescription("Upper threshold for the number of generated streamlines.");
		this->streamlineLengthThreshold_.setDescription("Maximum length of the generated streamlines.");
		this->absoluteMagnitudeThreshold_.setDescription("Flow data points outside the threshold interval will not be used for streamline construction.");
		this->fitAbsoluteMagnitudeThreshold_.setDescription("Automatically adapts the absolute magnitude threshold on input change");

		this->linesOutportEnabled_.setDescription("Enables the streamlines out port");
		this->rippingEnabled_.setDescription("Enable handling of the ripping edge case");
		this->mergingEnabled_.setDescription("Enable handling of the merging edge case");
		this->splittingEnabled_.setDescription("Enable handling of the splitting edge case");
    }

    ComputeInput prepareComputeInput() override;
	ComputeOutput compute(ComputeInput input, ProgressReporter &progressReporter) const override;
	void processComputeOutput(ComputeOutput output) override;
	void adjustPropertiesToInput() override;

private:
	void linesOutportEnabledChange();

	// Ports
	VolumePort volumeInport_;
	GeometryPort seedingInport_;
    GeometryPort geometryOutport_;
	StreamlineListPort streamlineOutport_;

	// Streamsurface settings
	IntProperty maxNumStreamlines_;						///< Upper threshold for the number of generated streamlines.
	IntProperty streamlineLengthThreshold_;     		///< Restrict number of elements
	FloatIntervalProperty absoluteMagnitudeThreshold_;  ///< Only magnitudes in this interval are used
	BoolProperty fitAbsoluteMagnitudeThreshold_;        ///< Fit magnitude on input change?
	IntProperty rippingThresholdAngle_;         		///< Angle at which the ripping edge case is detected
	OptionProperty<VolumeRAM::Filter> filterMode_;      ///< Filtering inside the dataset
	OptionProperty<InputSpace> inputSpace_;				///< Space of the seed points
	FloatVec3Property seedPointsOffset_;				///< Seed points offset
	FloatProperty stepSize_;							///< Size of one integration step.
	FloatOptionProperty velocityUnitConversion_;		///< Scaling factor for the step size.

	// Debug
	BoolProperty linesOutportEnabled_;					///< Enable optional streamline output port
	BoolProperty rippingEnabled_;						///< Enable handling of the ripping edge case
	BoolProperty mergingEnabled_;						///< Enable handling of the merging edge case
	BoolProperty splittingEnabled_;						///< Enable handling of the splitting edge case

    static const std::string loggerCat_;
};

}

#endif // VRN_STREAMSURFACECREATOR_H
