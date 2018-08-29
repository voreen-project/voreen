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

#ifndef VRN_FORCEDIRECTEDNODEGRAPH_H
#define VRN_FORCEDIRECTEDNODEGRAPH_H

#include "nodegraph.h"

namespace voreen {

/**
 * Represents a graph of nodes and connections
 **/
class ForceDirectedNodeGraph : public NodeGraph {
public:
    /**
     * Constructs a new empty node graph.
     **/
    ForceDirectedNodeGraph();

    /**
     * Default destructor, deletes all held NodeGraphNodes.
     **/
    virtual ~ForceDirectedNodeGraph();

    void initLayout();

    void applyForceDirectedLayoutStep(double spring, double charge, double damping);

    void calculateForces();

    void moveNodes();

    // - setter methods ---------------------------------------------------------------------------
    void SetSpringForce(double value);
    void SetSpringLength(double value);
    void SetRepulsionForce(double value);
    void SetHorizontalMagneticForce(double value);
    void SetVerticalMagneticForce(double value);
    void SetMagneticAlpha(double value);
    void SetMagneticBeta(double value);
    void SetEdgeAngleEqualization(double value);
    void SetDampingFactor(double value);

private:
    double springForce_;                ///< force that pulls two connected nodes together
    double springLength_;               ///< zero energy length of springs
    double repulsionForce_;             ///< repelling force between two nodes (that they don't get to close to each other)

    double horizontalMagneticForce_;    ///< force of horizontal aligned magnetic field
    double verticalMagneticForce_;      ///< force of vertical aligned magnetic field
    double magneticAlpha_;              ///< power of edge length in magnetic field
    double magneticBeta_;               ///< power of angle length in magnetic field

    double edgeAngleEqualization_;      ///< power of inter-edge angle equalization force

    double dampingFactor_;              ///< damping factor of iterative solver
};

}

#endif // VRN_FORCEDIRECTEDNODEGRAPH_H

