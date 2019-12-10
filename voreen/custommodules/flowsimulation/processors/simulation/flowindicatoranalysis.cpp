/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "flowindicatoranalysis.h"

#include "voreen/core/ports/conditions/portconditionvolumelist.h"

#include "modules/plotting/datastructures/plotcell.h"
#include "modules/plotting/datastructures/plotdata.h"
#include "modules/plotting/datastructures/plotrow.h"

#include "../../utils/utils.h"

namespace voreen {

FlowIndicatorAnalysis::FlowIndicatorAnalysis()
    : AsyncComputeProcessor()
    , parameterPort_(Port::INPORT, "input.parameter", "Parameter Port")
    , volumeListPort_(Port::INPORT, "input.volumelist", "Volume List Port")
    , outport_(Port::OUTPORT,"output.plot", "Plot Port")
{
    addPort(parameterPort_);
    addPort(volumeListPort_);
    volumeListPort_.addCondition(new PortConditionVolumeListEnsemble());
    addPort(outport_);
}

FlowIndicatorAnalysisInput FlowIndicatorAnalysis::prepareComputeInput() {

    auto volumes = volumeListPort_.getThreadSafeData();
    auto parameterSetEnsemble = parameterPort_.getThreadSafeData();

    std::vector<FlowIndicator> indicators;
    for(const FlowIndicator& indicator : parameterSetEnsemble->getFlowIndicators()) {
        // TODO: use only FLUX MEASURE?
        indicators.push_back(indicator);
    }

    if(indicators.empty()) {
        throw InvalidInputException("No indicators", InvalidInputException::S_IGNORE);
    }

    std::unique_ptr<PlotData> output(new PlotData(1, indicators.size()));
    output->setColumnLabel(0, "Time [s]");
    for(size_t i=0; i<indicators.size(); i++) {
        output->setColumnLabel(i+1, "Indicator " + std::to_string(indicators[i].id_));
    }

    return FlowIndicatorAnalysisInput {
        std::move(volumes),
        std::move(indicators),
        std::move(output)
    };
}

FlowIndicatorAnalysisOutput FlowIndicatorAnalysis::compute(FlowIndicatorAnalysisInput input, ProgressReporter& progressReporter) const {

    auto volumes = std::move(input.volumes);
    auto indicators = std::move(input.indicators);
    std::unique_ptr<PlotData> data = std::move(input.output);

    progressReporter.setProgress(0.0f);

    for(size_t t = 0; t < volumes->size(); t++) {
        const VolumeBase* volume = volumes->at(t);
        std::vector<PlotCellValue> values;
        values.push_back(PlotCellValue(volume->getTimestep()));
        for(const FlowIndicator& indicator : indicators) {
            tgt::vec3 v = utils::sampleDisk(volume, indicator.center_, indicator.normal_, indicator.radius_);
            float magnitude = tgt::length(v);
            values.push_back(PlotCellValue(magnitude));
        }
        data->insert(values);

        progressReporter.setProgress((t+1.0f) / volumes->size());
    }

    // If we only have a single row, copy and shift to render the line visible.
    if(data->getRowsCount() == 1) {
        std::vector<PlotCellValue> copy = data->getRow(0).getCells();
        copy[0] = PlotCellValue(copy[0].getValue() + 1.0); // Add one second.
        data->insert(copy);
    }

    progressReporter.setProgress(1.0f);

    return FlowIndicatorAnalysisOutput {
        std::move(data)
    };
}

void FlowIndicatorAnalysis::processComputeOutput(FlowIndicatorAnalysisOutput output) {
    outport_.setData(output.plotData.release(), true);
}

}
