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

#ifndef VRN_VESSELGRAPHPERTURBATION_H
#define VRN_VESSELGRAPHPERTURBATION_H

#include "voreen/core/processors/processor.h"

#include "../ports/vesselgraphport.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/optionproperty.h"

#include <random>

namespace voreen {

class VesselGraphPerturbation : public Processor {
public:
    VesselGraphPerturbation();
    virtual ~VesselGraphPerturbation();
    virtual std::string getCategory() const { return "Geometry"; }
    virtual std::string getClassName() const { return "VesselGraphPerturbation"; }
    virtual CodeState getCodeState() const { return Processor::CODE_STATE_EXPERIMENTAL; }
    virtual Processor* create() const { return new VesselGraphPerturbation(); }

protected:
    virtual void setDescriptions() {
        setDescription("This processor can be used to perturb VesselGraphs using a number of methods.");
    }

    enum PerturbationMethod {
        ADD_EDGES, // See NetMets: software for quantifying [...] network segmentation, figure 4 (a)
        SPLIT_NODES, // ... (b)
        SUBDIVIDE_EDGES, // ... (c)
        SPLIT_EDGES, // not shown, cut an edge in half, adding two nodes: o-------o -> o--o o--o

        // Not from Netmets
        MOVE_NODES, // move node positions
        CHANGE_PROPERTIES, // move node positions

        COMBINED, // Combines all other perturbations
    };

    virtual void process();

    // ports
    VesselGraphPort inport_;
    VesselGraphPort outport_;

    // properties
    BoolProperty enabled_;
    OptionProperty<PerturbationMethod> perturbationMethod_;
    FloatProperty perturbationAmount_; // \in [0,1]

    BoolProperty usePredeterminedSeed_;
    IntProperty predeterminedSeed_;

    std::random_device randomDevice;

    static const std::string loggerCat_;
};

} // namespace voreen
#endif // VRN_VESSELGRAPHPERTURBATION_H
