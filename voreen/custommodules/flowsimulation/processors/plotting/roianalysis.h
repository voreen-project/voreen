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

#ifndef VRN_ROIANALYSIS_H
#define VRN_ROIANALYSIS_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/ports/genericport.h"
#include "modules/plotting/ports/plotport.h"
#include "modules/plotting/datastructures/plotdata.h"
#include "../../ports/flowparametrizationport.h"

namespace voreen {

struct RoiAnalysisInput {
    PortDataPointer<VolumeList> volumes;
    std::vector<tgt::vec3> seedPoints;
    std::unique_ptr<PlotData> output;
    std::function<std::vector<float>(const std::vector<tgt::vec3>&)> outputFunc;
};

struct RoiAnalysisOutput {
    std::unique_ptr<PlotData> plotData;
};

class VRN_CORE_API RoiAnalysis : public AsyncComputeProcessor<RoiAnalysisInput, RoiAnalysisOutput> {
public:
    RoiAnalysis();
    virtual Processor* create() const { return new RoiAnalysis(); }

    virtual std::string getCategory() const  { return "Plotting"; }
    virtual std::string getClassName() const { return "RoiAnalysis"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL; }

protected:

    virtual void setDescriptions() override {
        setDescription("This processors allows to sample a roi defined by FlowIndicatorDetection processors. "
                       "The output can be plotted e.g. using a LinePlot processor.");
    }

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

private:

    static bool isTimeSeries(const VolumeList* list);

    VolumeListPort volumeListPort_;
    VolumePort maskPort_;
    PlotPort outport_;

    StringOptionProperty outputQuantity_;

    static const std::string loggerCat_;
};

} //namespace

#endif

