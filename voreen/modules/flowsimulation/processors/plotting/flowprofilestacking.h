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

#ifndef VRN_FLOWPROFILESTACKING_H
#define VRN_FLOWPROFILESTACKING_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/ports/genericport.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "modules/plotting/ports/plotport.h"
#include "modules/plotting/datastructures/plotdata.h"

namespace voreen {

struct FlowProfileStackingInput {
    PortDataPointer<VolumeList> volumes;
    std::vector<std::vector<tgt::vec3>> segments;
    int numSamplesRadius_;
    int numSamplesAngle_;
    int numSamplesZ_;
    tgt::vec2 zRange_;
    bool values_;
};

struct FlowProfileStackingOutput {
    std::vector<std::unique_ptr<PlotData>> output;
};

class VRN_CORE_API FlowProfileStacking : public AsyncComputeProcessor<FlowProfileStackingInput, FlowProfileStackingOutput> {
public:
    FlowProfileStacking();
    virtual Processor* create() const { return new FlowProfileStacking(); }

    virtual std::string getCategory() const  { return "Plotting"; }
    virtual std::string getClassName() const { return "FlowProfileStacking"; }
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

    StringOptionProperty stackingMode_;
    IntProperty selectedSegment_;
    IntProperty selectedSlice_;
    IntProperty numSamplesRadius_;
    IntProperty numSamplesAngle_;
    IntProperty numSamplesZ_;
    FloatIntervalProperty zRange_;

    std::vector<std::unique_ptr<PlotData>> plotData_;

    static const std::string loggerCat_;
};

} //namespace

#endif

