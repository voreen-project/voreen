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

#include "../include/nodegraphsource.h"

#include "modules/plotting/datastructures/plotdata.h"
#include "modules/plotting/datastructures/plotcell.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

#include "voreen/core/voreenapplication.h"
#include "tgt/assert.h"

#include <ctime>
#include <vector>
#include <set>
#include <list>

namespace voreen {

const std::string NodeGraphSource::loggerCat_("voreen.NodeGraphSource");

NodeGraphSource::NodeGraphSource()
    : Processor()
    , outPortNodes_(Port::OUTPORT, "Nodes")
    , outPortConnections_(Port::OUTPORT, "Connections")
    , nodeData_(0)
    , connectionData_(0)
    , nodeCount_("NodeCount", "Number of Nodes", 10, 1, 200)
    , connectionCount_("ConnectionCount", "Number of Connections", 6, 1, 100)
    , omitNotConnectedNodes_("OmitNotConnectedNodes_", "Remove not connected Nodes", true)
    , generateTree_("GenerateTree", "Generate Tree", false)
    , minimumChildCount_("MinChildCount", "Minimum Child Count", 1, 1, 10)
    , maximumChildCount_("MaxChildCount", "MaximumChild Count", 4, 1, 10)
{
    addProperty(nodeCount_);
    addProperty(generateTree_);

    addProperty(connectionCount_);
    addProperty(omitNotConnectedNodes_);
    addProperty(maximumChildCount_);
    maximumChildCount_.setVisibleFlag(false);

    connectionCount_.setGroupID("GenerateGraph");
    omitNotConnectedNodes_.setGroupID("GenerateGraph");

    maximumChildCount_.setGroupID("GenerateTree");

    addPort(outPortNodes_);
    addPort(outPortConnections_);

    setPropertyGroupGuiName("GenerateGraph", "Generate Graph");
    setPropertyGroupGuiName("GenerateTree", "Generate Tree");

    generateTree_.onChange(MemberFunctionCallback<NodeGraphSource>(this, &NodeGraphSource::updateProps));
}

NodeGraphSource::~NodeGraphSource() {
    // nothing to do here
}

Processor* NodeGraphSource::create() const {
    return new NodeGraphSource();
}

void NodeGraphSource::process() {
    tgtAssert(nodeData_ != 0 && connectionData_ != 0, "no plotdata object");
    
    if (! generateTree_.get())
        createRandomGraph();
    else 
        createRandomTree();

    Processor::setProgress(1.f);
}

void NodeGraphSource::initialize() {
    Processor::initialize();

    nodeData_ = new PlotData(0, 0);
    connectionData_ = new PlotData(0, 0);
    // initialize random seed:
    srand(static_cast<unsigned int>(time(0)));
}

void NodeGraphSource::deinitialize() {
    if (nodeData_ != outPortNodes_.getData())
        delete nodeData_;
    if (connectionData_ != outPortConnections_.getData())
        delete connectionData_;

    outPortNodes_.setData(0, true);
    outPortConnections_.setData(0, true);
    nodeData_ = 0;
    connectionData_ = 0;

    Processor::deinitialize();
}

void NodeGraphSource::createRandomGraph() {
    // step 0: some initialization
    PlotData* oldNodes = nodeData_;
    PlotData* oldConnections = connectionData_;

    PlotData* newNodes = new PlotData(1, 3);
    PlotData* newConnections = new PlotData(2, 0);

    plot_t rmax = static_cast<plot_t>(RAND_MAX);
    std::set<int> connectedNodes;

    // step 1: calculate maximum possible connection count for n nodes
    size_t sum = 0;
    for (int i = 1; i < nodeCount_.get(); ++i)
        sum += i;

    // step 2: gather all possible connections
    std::vector< std::vector<PlotCellValue> > allEdges;
    allEdges.resize(sum);
    int row = 0;
    for (int i = 0; i < nodeCount_.get(); ++i) {
        for (int j = i+1; j < nodeCount_.get(); ++j) {
            if (i != j) {
                PlotCellValue value(i);
                allEdges[row].push_back(value);
                allEdges[row].push_back(value);
                ++row;
            }
        }
    }

    // step 3: randomly pick m connections
    for (int i = 0; i < connectionCount_.get() && !allEdges.empty(); ++i) {
        std::vector< std::vector<PlotCellValue> >::iterator position = allEdges.begin();
        position += rand() % allEdges.size();
        newConnections->insert(*position);
        if (omitNotConnectedNodes_.get()) {
            connectedNodes.insert(static_cast<int>((*position)[0].getValue()));
            connectedNodes.insert(static_cast<int>((*position)[1].getValue()));
        }
        allEdges.erase(position);
    }

    // step 4: remove not connected nodes if neccesairy
    if (! omitNotConnectedNodes_.get()) {
        for (int i = 0; i < nodeCount_.get(); ++i)
            connectedNodes.insert(i);
    }

    // step 5: create n random nodes
    for (std::set<int>::iterator it = connectedNodes.begin(); it != connectedNodes.end(); ++it) {
        std::vector<PlotCellValue> node;
        node.push_back(PlotCellValue(*it));
        node.push_back(PlotCellValue(static_cast<plot_t>(rand())/rmax * 10.0));
        node.push_back(PlotCellValue(static_cast<plot_t>(rand())/rmax * 10.0));
        node.push_back(PlotCellValue(1.0));
        newNodes->insert(node);
    }

    // step 6: send data to PlotPorts
    nodeData_ = newNodes;
    connectionData_ = newConnections;

    outPortNodes_.setData(nodeData_);
    outPortConnections_.setData(connectionData_);

    delete oldNodes;
    delete oldConnections;
}


void NodeGraphSource::createRandomTree() {
    // step 0: some initialization
    PlotData* oldNodes = nodeData_;
    PlotData* oldConnections = connectionData_;

    PlotData* newNodes = new PlotData(1, 3);
    PlotData* newConnections = new PlotData(2, 0);

    plot_t rmax = static_cast<plot_t>(RAND_MAX);
    int currentNode = 0;
    std::list<int> nodesToProcess;

    // step 1: create initial node (root)
    std::vector<PlotCellValue> node;
    node.push_back(currentNode);
    node.push_back(PlotCellValue(static_cast<plot_t>(rand())/rmax * 10.0));
    node.push_back(PlotCellValue(static_cast<plot_t>(rand())/rmax * 10.0));
    node.push_back(PlotCellValue(1.0));
    newNodes->insert(node);
    nodesToProcess.push_back(currentNode);
    ++currentNode;

    // step 2: create tree structure
    while (currentNode < nodeCount_.get()) {
        int currentParent = nodesToProcess.front();
        nodesToProcess.pop_front();
        int numChildren = (rand() % maximumChildCount_.get()) + 1;

        for (int i = 0; i < numChildren; ++i) {
            std::vector<PlotCellValue> child;
            child.push_back(currentNode);
            child.push_back(PlotCellValue(static_cast<plot_t>(rand())/rmax * 10.0));
            child.push_back(PlotCellValue(static_cast<plot_t>(rand())/rmax * 10.0));
            child.push_back(PlotCellValue(1.0));
            newNodes->insert(child);
            nodesToProcess.push_back(currentNode);

            std::vector<PlotCellValue> edge;
            edge.push_back(currentParent);
            edge.push_back(currentNode);
            newConnections->insert(edge);

            ++currentNode;
        }
    }

    // step 3: send data to PlotPorts
    nodeData_ = newNodes;
    connectionData_ = newConnections;

    outPortNodes_.setData(nodeData_);
    outPortConnections_.setData(connectionData_);

    delete oldNodes;
    delete oldConnections;
}

void NodeGraphSource::updateProps() {
    if (generateTree_.get()) {
        connectionCount_.setVisibleFlag(false);
        omitNotConnectedNodes_.setVisibleFlag(false);
        maximumChildCount_.setVisibleFlag(true);
    }
    else {
        connectionCount_.setVisibleFlag(true);
        omitNotConnectedNodes_.setVisibleFlag(true);
        maximumChildCount_.setVisibleFlag(false);
    }
}

} // namespace voreen
