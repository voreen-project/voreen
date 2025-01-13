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

#ifndef VRN_RELATIVEPRESSUREFROMVESSELGRAPH_H
#define VRN_RELATIVEPRESSUREFROMVESSELGRAPH_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/geometryport.h"

#include "modules/vesselnetworkanalysis/ports/vesselgraphport.h"

namespace voreen {

/**
 * This processor is being used to generate simple synthetic simulation geometries.
 */
class VRN_CORE_API RelativePressureFromVesselGraph : public Processor {
public:
    RelativePressureFromVesselGraph();
    virtual Processor* create() const         { return new RelativePressureFromVesselGraph();    }

    virtual std::string getClassName() const  { return "RelativePressureFromVesselGraph";        }
    virtual std::string getCategory() const   { return "utility";                                }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL;                  }

protected:

    virtual void setDescriptions() {
        setDescription("This processor is used to calculate the relative pressure between two end nodes of "
                       "a vessel graph with a single edge.");
    }

    virtual bool isReady() const;
    virtual void process();

private:

    VesselGraphPort vesselGraphInport_;
    VolumePort pressureVolumeInport_;

    GeometryPort sampleRegionsOutport_;
    VolumePort relativePressureOutport_;

    FloatProperty relativePressure_;

    static const std::string loggerCat_;
};

}   //namespace

#endif
