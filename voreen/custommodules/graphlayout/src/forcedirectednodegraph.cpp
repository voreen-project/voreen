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

#include "../include/forcedirectednodegraph.h"
#include "tgt/assert.h"
#include "tgt/tgt_math.h"

#include <algorithm>

namespace voreen {

ForceDirectedNodeGraph::ForceDirectedNodeGraph()
    : NodeGraph()
    , springForce_(0)
    , springLength_(1)
    , repulsionForce_(0)
    , horizontalMagneticForce_(0)
    , verticalMagneticForce_(0)
    , magneticAlpha_(1)
    , magneticBeta_(1)
    , edgeAngleEqualization_(0)
    , dampingFactor_(0.5)
{
}

ForceDirectedNodeGraph::~ForceDirectedNodeGraph() {
}

void ForceDirectedNodeGraph::initLayout() {
    // initial layout inspired by Sugiyama/Misue
    double radius = (springLength_ * nodes_.size()) / (2 * tgt::PI);
    double angle  = (2 * tgt::PI) / nodes_.size();
    std::vector<NodeGraphNode*> shuffledNodes;
    shuffledNodes.reserve(nodes_.size());

    for (std::map<int, NodeGraphNode*>::iterator it = nodes_.begin(); it != nodes_.end(); ++it) {
        shuffledNodes.push_back(it->second);
    }
    std::random_shuffle(shuffledNodes.begin(), shuffledNodes.end());

    double currAngle = angle;
    for (std::vector<NodeGraphNode*>::iterator it = shuffledNodes.begin(); it != shuffledNodes.end(); ++it, currAngle += angle) {
        (*it)->position_ = tgt::dvec2(cos(currAngle) * radius, sin(currAngle) * radius);
    }
}

void ForceDirectedNodeGraph::applyForceDirectedLayoutStep(double spring, double charge, double damping) {
    //applyForces(spring, charge);
    springForce_ = spring;
    repulsionForce_ = charge;
    dampingFactor_ = damping;
    calculateForces();
    moveNodes();
}

void ForceDirectedNodeGraph::calculateForces() {
    std::map<int, NodeGraphNode*>::iterator firstIt, secondIt;

    for (firstIt = nodes_.begin(); firstIt != nodes_.end(); ++firstIt) {
        std::map<double, NodeGraphNode*> edgeAngles;

        for (secondIt = nodes_.begin(); secondIt != nodes_.end(); ++secondIt) {
            // don't apply forces, if nodes are the same
            if (firstIt->first != secondIt->first) {
                tgt::dvec2 direction = secondIt->second->position_ - firstIt->second->position_;
                tgt::dvec2 normDirection = tgt::normalize(direction);

                double distance = tgt::length(direction);
                double force = 0;

                // spring force
                if (firstIt->second->isConnectedTo(secondIt->second)) {
                    force = springForce_ * (distance - springLength_);
                }
                // repulsion force
                else {
                    force = -((firstIt->second->mass_ * secondIt->second->mass_) / (distance * distance)) * repulsionForce_;
                }

                firstIt->second->velocity_ += direction * (force / distance);

                // magnetic field force inspired by Sugiyama/Misue
                if (firstIt->second->isConnectedTo(secondIt->second)) {
                    double angle = atan2(direction.y, direction.x);
                    
                    if (edgeAngleEqualization_ > 0) 
                        edgeAngles.insert(std::make_pair(angle, secondIt->second));

                    bool quadrant24 = (angle >= 0 && angle < tgt::PI/2) || angle < -tgt::PI/2;

                    // calculate minimum positive angle to horizontal field
                    double horAngle = tgt::PI/2 - fabs((fabs(angle) - tgt::PI/2));
                    // calculate minimum positive angle to vertical field
                    double vertAngle = fabs(tgt::PI/2 - (fabs(angle)));

                    // calculate force of magnetic field
                    bool hor = ((horAngle / horizontalMagneticForce_) < (vertAngle / verticalMagneticForce_));
                    if (hor)
                        force = horizontalMagneticForce_ * pow(distance, magneticAlpha_) * pow(horAngle, magneticBeta_);
                    else
                        force = verticalMagneticForce_ * pow(distance, magneticAlpha_) * pow(vertAngle, magneticBeta_);

                    // calculate and apply horizontal field rotational forces if counter clockwise rotation:
                    if ((hor && !quadrant24) || (!hor && quadrant24)) {
                        firstIt->second->velocity_  += tgt::dvec2( normDirection.y * force, -normDirection.x * force);
                        secondIt->second->velocity_ += tgt::dvec2(-normDirection.y * force,  normDirection.x * force);
                    }
                    // calculate and apply vertical field rotational forces if clockwise rotation:
                    else {
                        firstIt->second->velocity_  += tgt::dvec2(-normDirection.y * force,  normDirection.x * force);
                        secondIt->second->velocity_ += tgt::dvec2( normDirection.y * force, -normDirection.x * force);
                    }
                }
            }
        }
        // equal inter-edge angles
        if (edgeAngleEqualization_ > 0 && edgeAngles.size() > 1) {
            // to avoid any special cases add first angle as last angle (+ 2PI)
            std::map<double, NodeGraphNode*>::const_iterator minIt, maxIt, curIt, nextIt;
            minIt = maxIt = curIt = nextIt = edgeAngles.begin();
            edgeAngles.insert(std::make_pair(curIt->first + 2*tgt::PI, curIt->second));

            // find pair of edges with minimum/maximum angle
            ++nextIt;
            double minAngle = nextIt->first - curIt->first;
            double maxAngle = minAngle;
            for (; nextIt != edgeAngles.end(); ++curIt, ++nextIt) {
                double angle = nextIt->first - curIt->first;
                if (angle < minAngle) {
                    minAngle = angle;
                    minIt = curIt;
                }
                if (angle > maxAngle) {
                    maxAngle = angle;
                    maxIt = curIt;
                }
            }

            // move node ...
            // ... towards smallest angle
            double movingAngle = (minIt->first + (++minIt)->first)/2;
            tgt::dvec2 move(cos(movingAngle), sin(movingAngle));
            move *= edgeAngleEqualization_ * (((2.0 * tgt::PI) / firstIt->second->getConnectionCount()) - minAngle);
            firstIt->second->velocity_ += move;
            // ... off largest angle
            movingAngle = (maxIt->first + (++maxIt)->first)/2 + tgt::PI;
            move = tgt::dvec2(cos(movingAngle), sin(movingAngle));
            move *= edgeAngleEqualization_ * (maxAngle - ((2.0 * tgt::PI) / firstIt->second->getConnectionCount()));
            firstIt->second->velocity_ += move;
        }
    }
}

void ForceDirectedNodeGraph::moveNodes() {
    std::map<int, NodeGraphNode*>::iterator firstIt;

    for (firstIt = nodes_.begin(); firstIt != nodes_.end(); ++firstIt) {
        firstIt->second->position_ += firstIt->second->velocity_;
        firstIt->second->velocity_ *= dampingFactor_;
    }
}

// - setter methods ---------------------------------------------------------------------------
void ForceDirectedNodeGraph::SetSpringForce(double value) {
    springForce_ = value;
}

void ForceDirectedNodeGraph::SetSpringLength(double value) {
    springLength_ = value;
}

void ForceDirectedNodeGraph::SetRepulsionForce(double value) {
    repulsionForce_ = value;
}

void ForceDirectedNodeGraph::SetHorizontalMagneticForce(double value) {
    horizontalMagneticForce_ = value;
}

void ForceDirectedNodeGraph::SetVerticalMagneticForce(double value) {
    verticalMagneticForce_ = value;
}

void ForceDirectedNodeGraph::SetMagneticAlpha(double value) {
    magneticAlpha_ = value;
}

void ForceDirectedNodeGraph::SetMagneticBeta(double value) {
    magneticBeta_ = value;
}

void ForceDirectedNodeGraph::SetEdgeAngleEqualization(double value) {
    edgeAngleEqualization_ = value;
}

void ForceDirectedNodeGraph::SetDampingFactor(double value) {
    dampingFactor_ = value;
}


} // namespace voreen
