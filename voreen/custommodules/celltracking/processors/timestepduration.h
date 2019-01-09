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

#ifndef VRN_TIMESTEPDURATION_H
#define VRN_TIMESTEPDURATION_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/textport.h"

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"

namespace voreen {

/**
 * Gets a volume list and a selected volume and outputs the current time for a time series.
 */
class VRN_CORE_API TimestepDuration : public Processor {

public:
    TimestepDuration();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "TimestepDuration";     }
    virtual std::string getCategory() const  { return "Utility";         }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL;  }

protected:
    virtual void setDescriptions() {
        setDescription("Outputs the current time of the selected time step in a time series.");
    }

    virtual void process();
    virtual void adjustPropertiesToInput();

    /// Inport for the volume list.
    VolumeListPort inport_;

    /// Outport for the elapsed time as text
    TextPort outport_;

    /// the first time step (which is defined as elapsed time = 0)
    IntProperty startTimestep_;

    /// how long is the time between two consecutive time steps (in seconds)
    FloatProperty timestepDuration_;

    /// the current time step for which the elapsed time should be computed
    IntProperty currentTimestep_;

    /// elapsed time for the current time step
    FloatProperty elapsedTime_;

    static const std::string loggerCat_;

};

} // namespace

#endif
