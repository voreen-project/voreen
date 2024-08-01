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

#ifndef VRN_VESSELGRAPHREFINEMENT_H
#define VRN_VESSELGRAPHREFINEMENT_H


#include "../datastructures/vesselgraph.h"
#include <functional>
#include <memory>

namespace voreen {

struct VesselGraphRefinement {
    typedef std::function<bool(const VesselGraphEdge&)> RemovableEdgeCheck;

    static std::unique_ptr<VesselGraph> removeEndEdgesRecursively(const VesselGraph& input, RemovableEdgeCheck isRemovableEdge, size_t maxIterations = std::numeric_limits<size_t>::max());

    // Remove edges if they are individually deletable (according to isRemoveableEdge)
    static std::unique_ptr<VesselGraph> removeAllEdges(const VesselGraph&, RemovableEdgeCheck);

    static std::unique_ptr<VesselGraph> removeEndEdges(const VesselGraph& input, RemovableEdgeCheck isRemovableEdge, size_t& num_removed_edges);

    // Will also be called after all edge removal functions
    static std::unique_ptr<VesselGraph> removeDregree2Nodes(const VesselGraph&);
};

} // namespace voreen
#endif // VRN_VESSELGRAPHREFINEMENT_H
