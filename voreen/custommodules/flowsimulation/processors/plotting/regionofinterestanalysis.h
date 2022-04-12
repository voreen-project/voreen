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

#ifndef VRN_REGIONOFINTERESTANALYSIS_H
#define VRN_REGIONOFINTERESTANALYSIS_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/ports/genericport.h"
#include "modules/plotting/ports/plotport.h"
#include "modules/plotting/datastructures/plotdata.h"
#include "../../ports/flowparametrizationport.h"

namespace voreen {

struct RegionOfInterestAnalysisInput {
    PortDataPointer<VolumeList> volumes;
    std::vector<std::vector<tgt::vec3>> segments;
};

struct RegionOfInterestAnalysisOutput {
    std::vector<std::unique_ptr<PlotData>> output;
};

class VRN_CORE_API RegionOfInterestAnalysis : public AsyncComputeProcessor<RegionOfInterestAnalysisInput, RegionOfInterestAnalysisOutput> {
public:
    RegionOfInterestAnalysis();
    virtual Processor* create() const { return new RegionOfInterestAnalysis(); }

    virtual std::string getCategory() const  { return "Plotting"; }
    virtual std::string getClassName() const { return "RegionOfInterestAnalysis"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL; }

protected:

    virtual void setDescriptions() override {
        setDescription("This processors allows to sample a roi defined by RegionOfInterest2D processor. "
                       "The output can be plotted e.g. using a LinePlot processor.");
    }

    virtual void beforeProcess();

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

private:

    void onInputDataChange();
    void onSelectedSegmentChange();

    VolumeListPort volumeListPort_;
    GeometryPort regionOfInterestPort_;
    PlotPort outport_;

    IntProperty selectedSegment_;

    std::vector<std::unique_ptr<PlotData>> plotData_;

    static const std::string loggerCat_;
};

} //namespace

#endif

