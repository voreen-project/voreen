/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_FLOWINDICATORDETECTION_H
#define VRN_FLOWINDICATORDETECTION_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/matrixproperty.h"
#include "voreen/core/properties/string/stringtableproperty.h"

#include "../../datastructures/flowparameters.h"
#include "../../ports/flowparametrizationport.h"

#include "custommodules/vesseltopology/ports/vesselgraphport.h"

namespace voreen {

/**
 * This processor is being used to select in and out flow.
 */
class VRN_CORE_API FlowIndicatorDetection : public Processor {
public:
    FlowIndicatorDetection();
    virtual Processor* create() const         { return new FlowIndicatorDetection();    }

    virtual std::string getClassName() const  { return "FlowIndicatorDetection";        }
    virtual std::string getCategory() const   { return "Simulation";                    }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL;         }

    virtual bool isReady() const;
    virtual void process();
    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

protected:

    virtual void adjustPropertiesToInput();

    virtual void setDescriptions() {
        setDescription("This processor is being used to select in and out flow.");
    }

private:

    void onSelectionChange();
    void onConfigChange();
    void onInputChange();
    void buildTable();

    VesselGraphPort vesselGraphPort_;
    VolumePort volumePort_;
    FlowParametrizationPort flowParametrizationPort_;

    StringProperty ensembleName_;
    FloatProperty simulationTime_;
    FloatProperty temporalResolution_;
    IntProperty spatialResolution_;
    IntProperty numTimeSteps_;

    OptionProperty<FlowFunction> flowFunction_;
    OptionProperty<FlowDirection> flowDirection_;
    FloatProperty radius_;

    StringTableProperty flowIndicatorTable_;

    IntProperty firstRefNode_;
    IntProperty numRefNodes_;
    IntProperty angleThreshold_;

    std::vector<FlowIndicator> flowIndicators_;

    static const std::string loggerCat_;
};

}   //namespace

#endif
