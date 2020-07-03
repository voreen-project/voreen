/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_FLOWINDICATORANALYSIS_H
#define VRN_FLOWINDICATORANALYSIS_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/ports/genericport.h"
#include "modules/plotting/ports/plotport.h"
#include "modules/plotting/datastructures/plotdata.h"
#include "../../ports/flowparametrizationport.h"

namespace voreen {

struct FlowIndicatorAnalysisInput {
    PortDataPointer<VolumeList> volumes;
    std::vector<FlowIndicator> indicators;
    std::unique_ptr<PlotData> output;
    std::function<std::vector<float>(const std::vector<tgt::vec3>&)> outputFunc;
    bool transformSamples;
};

struct FlowIndicatorAnalysisOutput {
    std::unique_ptr<PlotData> plotData;
};

class VRN_CORE_API FlowIndicatorAnalysis : public AsyncComputeProcessor<FlowIndicatorAnalysisInput, FlowIndicatorAnalysisOutput> {
public:
    FlowIndicatorAnalysis();
    virtual Processor* create() const { return new FlowIndicatorAnalysis(); }

    virtual std::string getCategory() const  { return "Plotting"; }
    virtual std::string getClassName() const { return "FlowIndicatorAnalysis"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL; }

protected:

    virtual void setDescriptions() override {
        setDescription("This processors allows to sample and disks defined by FlowIndicatorDetection processors. "
                       "The output can be plotted e.g. using a LinePlot processor. "
                       "It allows to export the curves into simple CSV files which can again be used by the "
                       "FlowIndicatorDetection processor to setup a velocity curve for a CFD simulation.");
    }

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

private:

    static bool isTimeSeries(const VolumeList* list);

    void onParametersChange();

    void exportVelocityCurve();

    FlowParametrizationPort parameterPort_;
    VolumeListPort volumeListPort_;
    PlotPort outport_;

    StringOptionProperty outputQuantity_;
    StringOptionProperty indicator_;
    BoolProperty transformSamples_;

    FileDialogProperty exportCurvePath_;
    ButtonProperty saveButton_;

    static const std::string loggerCat_;
};

} //namespace
#endif

