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

#ifndef VRN_FLOWCHARACTERISTICS_H
#define VRN_FLOWCHARACTERISTICS_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/volumeport.h"

#ifdef VRN_MODULE_VESSELNETWORKANALYSIS
#include "modules/vesselnetworkanalysis/ports/vesselgraphport.h"
#endif

namespace voreen {

/**
 * Used to extract and set characteristics of a time series.
 */
class VRN_CORE_API FlowCharacteristics : public Processor {
public:
    FlowCharacteristics();

    virtual Processor* create() const { return new FlowCharacteristics(); }
    virtual std::string getCategory() const  { return "Utility"; }
    virtual std::string getClassName() const { return "FlowCharacteristics"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL; }

    virtual bool isReady() const;
    virtual void process();

protected:
    virtual void setDescriptions() {
        setDescription("Used to extract and set characteristics of a flow time series.");
    }

private:

    VolumeListPort inport_;
#ifdef VRN_MODULE_VESSELNETWORKANALYSIS
    VesselGraphPort vesselGraphPort_;
#endif

    FloatProperty temporalResolution_;
    FloatProperty characteristicLength_;
    FloatProperty minVelocity_;
    FloatProperty maxVelocity_;
    ButtonProperty resetButton_;
};

}

#endif

