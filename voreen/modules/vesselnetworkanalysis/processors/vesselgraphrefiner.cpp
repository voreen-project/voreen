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

#include "vesselgraphrefiner.h"

#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "../algorithm/vesselgraphrefinement.h"

namespace voreen {



const std::string VesselGraphRefiner::loggerCat_("voreen.vesselgraphrefiner");


VesselGraphRefiner::VesselGraphRefiner()
    : Processor()
    , inport_(Port::INPORT, "graph.input", "Graph Input", false, Processor::INVALID_RESULT)
    , outport_(Port::OUTPORT, "graph.output", "Normalized Graph Output", false, Processor::VALID)
    , enabled_("enabled", "Enabled", true)
    , refinementMethod_("refinementMethod", "Refinement Method")
    , maxIterations_("maxIterations", "Max Iterations", 1000, 1, 1000)
    , minVoxelLength_("minVoxelLength", "Min Voxel Length", 0, 0, 50)
    , minElongation_("minElongation", "Minimum Elongation", 0, 0, 5)
    , minBulgeSize_("minBulgeSize", "Minimum Bulge Size", 0, 0, 10)
{

    addPort(inport_);
    addPort(outport_);

    addProperty(enabled_);
    addProperty(minVoxelLength_);
    addProperty(minElongation_);
    addProperty(minBulgeSize_);

    addProperty(refinementMethod_);
        refinementMethod_.addOption("all", "All", ALL);
        refinementMethod_.addOption("end_recursive", "End Recursive", END_RECURSIVE);
        refinementMethod_.selectByValue(END_RECURSIVE);
    addProperty(maxIterations_);
}

VesselGraphRefiner::~VesselGraphRefiner() {
}

void VesselGraphRefiner::process() {
    const VesselGraph* input = inport_.getData();
    if(!input) {
        outport_.setData(nullptr);
        return;
    }
    if(!enabled_.get()) {
        outport_.setData(input, false);
        return;
    }

    std::unique_ptr<VesselGraph> output(nullptr);
    switch(refinementMethod_.getValue()) {
        case ALL:
            output = VesselGraphRefinement::removeAllEdges(*input, createRemovableEdgePredicate());
            break;
        case END_RECURSIVE:
            output = VesselGraphRefinement::removeEndEdgesRecursively(*input, createRemovableEdgePredicate(), maxIterations_.get());
            break;
    }

    outport_.setData(output.release());
}

std::function<bool(const VesselGraphEdge& edge)> VesselGraphRefiner::createRemovableEdgePredicate() const {
    size_t minVoxelLength = minVoxelLength_.get();
    float minElongation = minElongation_.get();
    float minBulgeSize = minBulgeSize_.get();
    return [minVoxelLength, minElongation, minBulgeSize] (const VesselGraphEdge& edge) {
        return !edge.hasValidData() || edge.getVoxels().size() < minVoxelLength || edge.getElongation() < minElongation || edge.getRelativeBulgeSize() < minBulgeSize;
    };
}


} // namespace voreen
