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

#ifndef VRN_VESSELGRAPHREFINER_H
#define VRN_VESSELGRAPHREFINER_H

#include "voreen/core/processors/processor.h"

#include "../ports/vesselgraphport.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/optionproperty.h"

namespace voreen {

class VesselGraphRefiner : public Processor {
public:
    VesselGraphRefiner();
    virtual ~VesselGraphRefiner();
    virtual std::string getCategory() const { return "Geometry"; }
    virtual std::string getClassName() const { return "VesselGraphRefiner"; }
    virtual CodeState getCodeState() const { return Processor::CODE_STATE_EXPERIMENTAL; }
    virtual Processor* create() const { return new VesselGraphRefiner(); }

protected:
    virtual void setDescriptions() {
        setDescription("This processor implements a refinement step of the paper \"Scalable Robust Graph and Feature Extraction for Arbitrary Vessel Networks in Volumetric Datasets\" by Drees et al. "
                "Please note that this processor does not (and cannot) recompute the properties of the remaining edges ."
                "It should therefore typically only be used for debugging purposes or as a quick way to find an appropriate bulge size parameter for a full extraction using <b>VesselGraphCrator</b>."
                );
        enabled_.setDescription("If not enabled, the original graph is passed on unmodified.");
        refinementMethod_.setDescription("'End Recursive' is the method used in the above mentioned paper. 'All' also tries to remove central edges in the graph.");
        minBulgeSize_.setDescription("Edges with a bulge size below this threshold will be considered for deletion during the refinement. A bulge size of 1.0 roughly corresponds to hemisphere-shaped bulge.");
    }

    enum RefinementMethod {
        ALL,
        END_RECURSIVE
    };

    virtual void process();

    std::function<bool(const VesselGraphEdge& edge)> createRemovableEdgePredicate() const;

    //void adjustPropertiesToInput();

    VesselGraphPort inport_;
    VesselGraphPort outport_;

    // properties
    BoolProperty enabled_;
    OptionProperty<RefinementMethod> refinementMethod_;
    IntProperty maxIterations_;
    IntProperty minVoxelLength_;
    FloatProperty minElongation_;
    FloatProperty minBulgeSize_;

    static const std::string loggerCat_;
};

} // namespace voreen
#endif // VRN_VESSELGRAPHREFINER_H
