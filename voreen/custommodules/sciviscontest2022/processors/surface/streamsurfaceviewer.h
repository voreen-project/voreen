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

#ifndef VRN_STREAMSURFACEVIEWER_H
#define VRN_STREAMSURFACEVIEWER_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/ports/geometryport.h"
#include "voreen/core/ports/volumeport.h"

#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"

#include <map>
#include <vector>
#include <string>

namespace voreen {

enum BlendFunction {
	MIN,
	MAX,
	ADD,
	SUBTRACT,
	MULTIPLY,
	AVERAGE
};

struct StreamSurfaceViewerInput {
	VolumeRAM::Filter filter;
	BlendFunction blendFunction;
	tgt::vec4 baseColor;
	std::unique_ptr<Geometry> output;
	PortDataPointer<Geometry> geometry;
	std::vector<PortDataPointer<VolumeBase>> volumes;
	std::vector<std::unique_ptr<TransFunc1DKeys>> transferFunctions;
};

struct StreamSurfaceViewerOutput {
	std::unique_ptr<Geometry> geometry;
};

struct StreamSurfaceViewerVolumeSettings {
	const VolumeBase* volume;
	TransFunc1DKeys transFuncKeys;
};

class VRN_CORE_API StreamSurfaceViewer: public AsyncComputeProcessor<StreamSurfaceViewerInput, StreamSurfaceViewerOutput> {
public:
	StreamSurfaceViewer();

	Processor* create() const override { return new StreamSurfaceViewer(); }

	std::string getCategory() const override { return "Stream surface Processing"; }
	std::string getClassName() const override { return "StreamSurfaceViewer"; }
	Processor::CodeState getCodeState() const override { return CODE_STATE_EXPERIMENTAL; }

protected:
	void setDescriptions() override {
		setDescription("This processor is used to map volume data onto the color "
					   "of a stream surface");
	}

	ComputeInput prepareComputeInput() override;
	ComputeOutput compute(ComputeInput input, ProgressReporter &progressReporter) const override;
	void processComputeOutput(ComputeOutput output) override;

private:
	void onVolumeChange();
	void onVolumeSelected();
	void onVolumeSettingsChange();

	// Ports
	VolumePort volumeInport_;
	GeometryPort geometryInport_;
	GeometryPort geometryOutport_;

	// Settings
	OptionProperty<BlendFunction> blendFunction_;		///< Blend function of the mapped colors.
	OptionProperty<VolumeRAM::Filter> filterMode_;		///< Filter modes of the volume.
	ColorProperty baseColor_;							///< Base color blending color.
	StringOptionProperty volumeSelector_;				///< Option for modifying properties pertaining to a specific volume.
	TransFunc1DKeysProperty transferFunction_;          ///< tf for color coding the volume samples.

	// Per Volume settings
	std::map<std::string, StreamSurfaceViewerVolumeSettings> volumeSettings_;
};

}

#endif // VRN_STREAMSURFACEVIEWER_H
