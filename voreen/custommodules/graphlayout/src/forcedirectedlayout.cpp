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

#include "../include/forcedirectedlayout.h"
#include "voreen/core/processors/processorwidgetfactory.h"

#include "modules/plotting/datastructures/plotdata.h"
#include "modules/plotting/datastructures/plotrow.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

namespace voreen {

const std::string ForceDirectedLayout::loggerCat_("voreen.Plotting.ForceDirectedLayout");


ForceDirectedLayout::ForceDirectedLayout()
    : Processor()
    , nodesIn_(Port::INPORT, "PlotData.IncomingNodes")
    , nodesOut_(Port::OUTPORT, "PlotData.OutgoingNodes")
    , edgesIn_(Port::INPORT, "PlotData.IncomingEdges")
    , edgesOut_(Port::OUTPORT, "PlotData.OutgoingEdges")
    , springForce_("springForce", "Spring Force", 1, 0, 2)
    , springLength_("springLength", "Spring Length", 1, 0, 10)
    , charge_("charge", "Charge", 1, 0, 4)
    , horizontalMagneticForce_("horizontalMagneticForce", "Horizontal Force", 0, 0, 2)
    , verticalMagneticForce_("verticalMagneticForce", "Vertical Force", 0, 0, 2)
    , magneticAlpha_("magneticAlpha", "Magnetic Power of Distance", 1, 0, 2)
    , magneticBeta_("magneticBeta", "Magnetic Power of Angle", 1, 0, 4)
    , edgeAngleEqualization_("EdgeAngleEqualization", "Edge-Angle Equalization Power", 0, 0, 2)
    , damping_("damping", "Damping", .5, 0, 1)
    , threshold_("threshold", "Equilibrium Threshold", .05, 0, 5)
    , initLayout_("initLayout", "Initial Layout", Processor::INVALID_RESULT)
    , applyForcesStep_("applyForcesStep", "Apply Forces Once", Processor::INVALID_RESULT)
    , applyForces_("applyForces", "Apply Forces", Processor::INVALID_RESULT)
{
    addPort(nodesIn_);
    addPort(nodesOut_);
    addPort(edgesIn_);
    addPort(edgesOut_);

    addProperty(springForce_);
    addProperty(springLength_);
    addProperty(charge_);
    addProperty(horizontalMagneticForce_);
    addProperty(verticalMagneticForce_);
    addProperty(magneticAlpha_);
    addProperty(magneticBeta_);
    addProperty(edgeAngleEqualization_);
    addProperty(damping_);
    addProperty(threshold_);
    addProperty(initLayout_);
    addProperty(applyForces_);
    addProperty(applyForcesStep_);

    springForce_.setGroupID("SimpleForces");
    springLength_.setGroupID("SimpleForces");
    charge_.setGroupID("SimpleForces");

    horizontalMagneticForce_.setGroupID("MagneticForces");
    verticalMagneticForce_.setGroupID("MagneticForces");
    magneticAlpha_.setGroupID("MagneticForces");
    magneticBeta_.setGroupID("MagneticForces");
    edgeAngleEqualization_.setGroupID("MagneticForces");

    damping_.setGroupID("Layout");
    threshold_.setGroupID("Layout");
    initLayout_.setGroupID("Layout");
    applyForces_.setGroupID("Layout");
    applyForcesStep_.setGroupID("Layout");
        
    setPropertyGroupGuiName("SimpleForces", "Simple Forces");
    setPropertyGroupGuiName("MagneticForces", "Magnetic Forces");
    setPropertyGroupGuiName("Layout", "Layout");

    initLayout_.onChange(MemberFunctionCallback<ForceDirectedLayout>(this, &ForceDirectedLayout::initLayout));
    applyForces_.onChange(MemberFunctionCallback<ForceDirectedLayout>(this, &ForceDirectedLayout::applyForces));
    applyForcesStep_.onChange(MemberFunctionCallback<ForceDirectedLayout>(this, &ForceDirectedLayout::applyForcesStep));
}

ForceDirectedLayout::~ForceDirectedLayout() {
}

Processor* ForceDirectedLayout::create() const {
    return new ForceDirectedLayout();
}

bool ForceDirectedLayout::isEndProcessor() const {
    return !(nodesOut_.isConnected() && edgesOut_.isConnected());
}

bool ForceDirectedLayout::isReady() const {
    return true;
}

void ForceDirectedLayout::process() {
    if (nodesIn_.isReady() && edgesIn_.isReady() && (nodesIn_.hasChanged() || edgesIn_.hasChanged()))
        readInports();
    writeOutports();
}

void ForceDirectedLayout::initialize() {
    Processor::initialize();
    readInports();
}

void ForceDirectedLayout::deinitialize() {
    nodesOut_.setData(0, false);
    edgesOut_.setData(0, false);
    Processor::deinitialize();
}

void ForceDirectedLayout::readInports() {
    nodeGraph_.clearNodes();

    const int idCol    = 0;
    const int xCol     = 1;
    const int yCol     = 2;
    const int massCol  = 3;
    const int conn1Col = 0;
    const int conn2Col = 1;

    const PlotData* nodes = dynamic_cast<const PlotData*>(nodesIn_.getData());
    const PlotData* edges = dynamic_cast<const PlotData*>(edgesIn_.getData());

    if (nodes == 0 || edges == 0)
        return;

    // first convert nodes
    for (std::vector<PlotRowValue>::const_iterator it = nodes->getRowsBegin(); it != nodes->getRowsEnd(); ++it) {
        nodeGraph_.addNode(static_cast<int>(it->getValueAt(idCol)), it->getValueAt(xCol), it->getValueAt(yCol), it->getValueAt(massCol));
    }

    // then convert edges
    for (std::vector<PlotRowValue>::const_iterator it = edges->getRowsBegin(); it != edges->getRowsEnd(); ++it) {
        nodeGraph_.connectNodes(static_cast<int>(it->getValueAt(conn1Col)), static_cast<int>(it->getValueAt(conn2Col)));
    }

    calculateForces();
}

void ForceDirectedLayout::calculateForces() {
    nodeGraph_.SetSpringForce(springForce_.get());
    nodeGraph_.SetSpringLength(springLength_.get());
    nodeGraph_.SetRepulsionForce(charge_.get());
    nodeGraph_.SetHorizontalMagneticForce(horizontalMagneticForce_.get());
    nodeGraph_.SetVerticalMagneticForce(verticalMagneticForce_.get());
    nodeGraph_.SetMagneticAlpha(magneticAlpha_.get());
    nodeGraph_.SetMagneticBeta(magneticBeta_.get());
    nodeGraph_.SetDampingFactor(damping_.get());
    nodeGraph_.SetEdgeAngleEqualization(edgeAngleEqualization_.get());

    nodeGraph_.calculateForces();
}

void ForceDirectedLayout::moveNodes() {
    nodeGraph_.moveNodes();
}

void ForceDirectedLayout::initLayout() {
    nodeGraph_.SetSpringForce(springForce_.get());
    nodeGraph_.SetSpringLength(springLength_.get());

    nodeGraph_.initLayout();    
    calculateForces();
}

void ForceDirectedLayout::applyForces() {
    double threshold = threshold_.get();
    int step = 0;
    double maxVelocity = 0;
    do {
        moveNodes();
        calculateForces();

        maxVelocity = 0;
        for (std::map<int, NodeGraphNode*>::const_iterator firstIt = nodeGraph_.getNodesBegin(); firstIt != nodeGraph_.getNodesEnd(); ++firstIt) {
            maxVelocity = std::max(maxVelocity, tgt::length(firstIt->second->velocity_));
        }
    } while (++step <= 256 && maxVelocity > threshold);

}

void ForceDirectedLayout::applyForcesStep() {
    moveNodes();
    calculateForces();
    writeOutports();
}

void ForceDirectedLayout::writeOutports() {
    const PlotData* edges = dynamic_cast<const PlotData*>(edgesIn_.getData());
    if (edges == 0)
        return;

    PlotData* newNodes = new PlotData(1, 4);
    newNodes->setColumnLabel(0, "id");
    newNodes->setColumnLabel(1, "x");
    newNodes->setColumnLabel(2, "y");
    newNodes->setColumnLabel(3, "dx");
    newNodes->setColumnLabel(4, "dy");

    PlotData* newEdges = new PlotData(2, 0);
    newEdges->setColumnLabel(0, "first");
    newEdges->setColumnLabel(1, "second");

    for (std::map<int, NodeGraphNode*>::const_iterator it = nodeGraph_.getNodesBegin(); it != nodeGraph_.getNodesEnd(); ++it) {
        std::vector<PlotCellValue> toInsert;
        toInsert.push_back(PlotCellValue(static_cast<plot_t>(it->first)));
        toInsert.push_back(PlotCellValue(it->second->position_.x));
        toInsert.push_back(PlotCellValue(it->second->position_.y));
        toInsert.push_back(PlotCellValue(it->second->velocity_.x));
        toInsert.push_back(PlotCellValue(it->second->velocity_.y));
        newNodes->insert(toInsert);
    }

    for (std::vector<PlotRowValue>::const_iterator it = edges->getRowsBegin(); it != edges->getRowsEnd(); ++it) {
        std::vector<PlotCellValue> toInsert;
        toInsert.push_back(PlotCellValue(static_cast<plot_t>(it->getValueAt(0))));
        toInsert.push_back(PlotCellValue(static_cast<plot_t>(it->getValueAt(1))));
        newEdges->insert(toInsert);
    }

    nodesOut_.setData(newNodes, true);
    edgesOut_.setData(newEdges, true);
}



}
