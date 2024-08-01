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

#include "voreen/core/network/networkgraph.h"

#include <algorithm>
#include <deque>
#include <iostream>
#include <stack>
#include <vector>

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/coprocessorport.h"

namespace voreen {

NetworkGraph::GraphNode::GraphNode(Processor* const processor, const int id)
    : processor_(processor),
    id_(id),
    successors_(),
    inDegree_(0),
    outDegree_(0),
    marked_(false),
    minDistanceToRoot_(0),
    tDiscovered_(0),
    tFinished_(0)
{
}

NetworkGraph::GraphNode::GraphNode(Processor* const processor)
    : processor_(processor),
    id_(-1),
    successors_(),
    inDegree_(0),
    outDegree_(0),
    marked_(false),
    minDistanceToRoot_(0),
    tDiscovered_(0),
    tFinished_(0)
{
}

NetworkGraph::GraphNode::GraphNode(const int id)
    : processor_(0),
    id_(id),
    successors_(),
    inDegree_(0),
    outDegree_(0),
    marked_(false),
    minDistanceToRoot_(0),
    tDiscovered_(0),
    tFinished_(0)
{
}

NetworkGraph::GraphNode::GraphNode(const NetworkGraph::GraphNode& other) {
    copy(other);
}

NetworkGraph::GraphNode::~GraphNode() {
}

bool NetworkGraph::GraphNode::operator<(const NetworkGraph::GraphNode& other) const {
    return (id_ < other.id_);
}

bool NetworkGraph::GraphNode::operator==(const NetworkGraph::GraphNode& other) const {
    return ((id_ == other.id_) && (processor_ == other.processor_));
}

bool NetworkGraph::GraphNode::operator!=(const NetworkGraph::GraphNode& other) const {
    return !(*this == other);
}

NetworkGraph::GraphNode& NetworkGraph::GraphNode::operator=(const NetworkGraph::GraphNode& other) {
    copy(other);
    return *this;
}

void NetworkGraph::GraphNode::addSuccessor(NetworkGraph::GraphNode* const successor)
{
    if (successor == 0)
        return;

    successors_.push_back(successor);
    ++(successor->inDegree_);
    ++outDegree_;
}

void NetworkGraph::GraphNode::clearSuccessors() {
    for (size_t i = 0; i < successors_.size(); ++i) {
        if (successors_[i]->inDegree_ > 0)
            --(successors_[i]->inDegree_);
    }
    successors_.clear();
    outDegree_ = 0;
}

// private methods
//

void NetworkGraph::GraphNode::copy(const NetworkGraph::GraphNode& other) {
    processor_ = other.processor_;  // shallow copy only!
    id_ = other.id_;
    successors_ = std::vector<GraphNode*>();    // do not copy pointers to successors!
    inDegree_ = 0;
    outDegree_ = 0;
    marked_ = other.marked_;
    minDistanceToRoot_ = other.minDistanceToRoot_;
    tDiscovered_ = other.tDiscovered_;
    tFinished_= other.tFinished_;
}

// ============================================================================

NetworkGraph::CollectSuccessorsVisitor::CollectSuccessorsVisitor() {}

NetworkGraph::CollectSuccessorsVisitor::~CollectSuccessorsVisitor() {}

const std::set<NetworkGraph::GraphNode*>& NetworkGraph::CollectSuccessorsVisitor::getSuccessors() const {
    return successors_;
}

const std::vector<NetworkGraph::GraphNode*>& NetworkGraph::CollectSuccessorsVisitor::getSuccessorsOrdered() const {
    return successorOrdered_;
}

bool NetworkGraph::CollectSuccessorsVisitor::visit(NetworkGraph::GraphNode* const graphNode) {
    if (graphNode) {
        successors_.insert(graphNode);
        successorOrdered_.push_back(graphNode);
        return true;
    }
    return false;
}

// ============================================================================

const std::string NetworkGraph::loggerCat_("voreen.NetworkGraph");

NetworkGraph::NetworkGraph(const std::vector<Processor*>& processors)
    : nodes_(),
      nodeCounter_(0),
      networkGraphTransposed_(0)
{
    buildGraph(processors, LoopPortTypeCheck(true), LoopPortTypeCheck(), true);
}

NetworkGraph::NetworkGraph(const std::vector<Processor*>& processors,
                           const PortTypeCheck& type, const PortTypeCheck& loopType, const bool keepNodes)
    : nodes_(),
      nodeCounter_(0),
      networkGraphTransposed_(0)
{
    buildGraph(processors, type, loopType, keepNodes);
}

NetworkGraph::NetworkGraph(const NetworkGraph& other) {
    networkGraphTransposed_ = 0;
    copy(other);
}

NetworkGraph::~NetworkGraph() {
    clearGraph();
}

NetworkGraph& NetworkGraph::operator=(const NetworkGraph& other) {
    networkGraphTransposed_ = 0;
    copy(other);
    return *this;
}

bool NetworkGraph::operator==(const NetworkGraph& other) const {

    // number of nodes equal?
    if (nodeCounter_ != other.nodeCounter_)
        return false;

    // compare nodes
    NodeSet::const_iterator iter;
    for (iter = nodes_.begin(); iter != nodes_.end(); ++iter) {
        GraphNode* node = *iter;
        // node contained by other graph?
        if (other.nodes_.find(node) == other.nodes_.end())
            return false;
        // other node equal?
        GraphNode* otherNode = *(other.nodes_.find(node));
        if ( *node != *otherNode)
            return false;
        // in/out degrees equal?
        if (node->inDegree_ != otherNode->inDegree_)
            return false;
        if (node->outDegree_ != otherNode->outDegree_)
            return false;
        // successor counts equal?
        if (node->getSuccessors().size() != otherNode->getSuccessors().size())
            return false;
        // successors equal (consider that order of successors might vary)?
        for (size_t i = 0; i<node->getSuccessors().size(); ++i) {
            bool found = false;
            for (size_t j = 0; j<otherNode->getSuccessors().size() && !found; ++j) {
                found = (*(node->getSuccessors().at(i)) == *(otherNode->getSuccessors().at(j)));
            }
            if (!found)
                return false;
        }
        for (size_t i = 0; i<otherNode->getSuccessors().size(); ++i) {
            bool found = false;
            for (size_t j = 0; j<node->getSuccessors().size() && !found; ++j) {
                found = (*(otherNode->getSuccessors().at(i)) == *(node->getSuccessors().at(j)));
            }
            if (!found)
                return false;
        }
    }

    return true;
}

bool NetworkGraph::operator!=(const NetworkGraph& other) const {
    return !(*this == other);
}

bool NetworkGraph::containsProcessor(Processor* const processor) const {
    return (!findNodes(processor).empty());
}

bool NetworkGraph::isSuccessor(Processor* predecessor, Processor* successor) const {

    std::set<Processor*> predecessors;
    predecessors.insert(predecessor);
    const std::set<Processor*> successors = getSuccessors(predecessors);

    return (successors.find(successor) != successors.end());
}

bool NetworkGraph::isSuccessor(Port* predecessorPort, Processor* successor) const {
    tgtAssert(predecessorPort, "null pointer passed");
    tgtAssert(successor, "null pointer passed");
    if (predecessorPort->isOutport())
        return (getSuccessors(predecessorPort).count(successor));
    else if (predecessorPort->isInport())
        return (isSuccessor(predecessorPort->getProcessor(), successor));
    else {
        LERROR("isSuccessor(): passed port is neither inport nor outport");
        return false;
    }
}

bool NetworkGraph::isSuccessor(Port* predecessor, Port* successor) const {
    tgtAssert(predecessor, "null pointer passed");
    tgtAssert(successor, "null pointer passed");

    if (predecessor == successor)
        return true;
    else if (predecessor->isInport() && successor->isInport())
        return (isSuccessor(predecessor->getProcessor(), successor));
    else if (predecessor->isOutport() && successor->isInport())
        return (isSuccessor(predecessor, successor->getProcessor()) && isSuccessor(predecessor->getProcessor(), successor));
    else if (predecessor->isOutport() && successor->isOutport())
        return (isSuccessor(predecessor, successor->getProcessor()));
    else if (predecessor->isInport() && successor->isOutport())
        return (isSuccessor(predecessor->getProcessor(), successor->getProcessor()));
    else {
        LERROR("isSuccessor(): passed ports are neither inports not outports");
        return false;
    }
}

bool NetworkGraph::isSuccessor(Processor* predecessor, Port* successorPort) const {
    tgtAssert(predecessor, "null pointer passed");
    tgtAssert(successorPort, "null pointer passed");

    if (successorPort->isInport())
        return (getPredecessors(successorPort).count(predecessor) > 0);
    else if (successorPort->isOutport())
        return (isSuccessor(predecessor, successorPort->getProcessor()));
    else {
        LERROR("isSuccessor(): passed port is neither inport nor outport");
        return false;
    }
}

bool NetworkGraph::isPathElement(Processor* processor, Processor* pathRoot, Processor* pathEnd) const {
    return (isSuccessor(pathRoot, processor) && isSuccessor(processor, pathEnd));
}

bool NetworkGraph::isPathElement(Processor* processor, Port* pathRoot, Port* pathEnd) const {
    tgtAssert(processor, "null pointer passed");
    tgtAssert(pathRoot, "null pointer passed");
    tgtAssert(pathEnd, "null pointer passed");

    if (pathRoot == pathEnd)
        return (pathRoot->getProcessor() == processor);
    else if (pathRoot->isInport() && pathEnd->isInport())
        return (isSuccessor(pathRoot->getProcessor(), processor) && isSuccessor(processor, pathEnd));
    else if (pathRoot->isOutport() && pathEnd->isInport())
        return (isSuccessor(pathRoot, processor) && isSuccessor(processor, pathEnd));
    else if (pathRoot->isOutport() && pathEnd->isOutport())
        return (isSuccessor(pathRoot, processor) && isSuccessor(processor, pathEnd->getProcessor()));
    else if (pathRoot->isInport() && pathEnd->isOutport())
        return (isSuccessor(pathRoot->getProcessor(), processor) && isSuccessor(processor, pathEnd->getProcessor()));
    else {
        LERROR("isSuccessor(): passed ports are neither inports not outports");
        return false;
    }
}

bool NetworkGraph::isPathElement(Port* port, Port* pathRoot, Port* pathEnd) const {
        return (isSuccessor(pathRoot, port) && isSuccessor(port, pathEnd));
}

void NetworkGraph::fullTraverseBreadthFirst(GraphVisitor* const visitor) const {
    setAllNodesMarked(false);
    std::vector<GraphNode*> roots = identifyRoots();
    for (size_t i = 0; i < roots.size(); ++i) {
        if (roots[i]->isMarked() == false)
            traverseBreadthFirstInternal(roots[i], visitor);
    }
}

void NetworkGraph::fullTraverseDepthFirst(GraphVisitor* const onDiscoveryVisitor,
                                          GraphVisitor* const onFinishVisitor) const
{
    setAllNodesDiscovered(0);
    setAllNodesFinished(0);

    int startTime = 0;      // start timestamp is 0 for first run
    std::vector<GraphNode*> roots = identifyRoots();
    for (size_t i = 0; i < roots.size(); ++i) {
        if (roots[i]->isDiscovered() == false)
            startTime = traverseDepthFirstInternal(roots[i], startTime, onDiscoveryVisitor, onFinishVisitor);
    }
}


std::vector<Processor*> NetworkGraph::getProcessors(const bool sortByID) const {
    std::vector<Processor*> processors;

    // If no sorting is required, simply copy the processor pointers to a vector.
    //
    if (sortByID == false) {
        for (NodeSet::const_iterator it = nodes_.begin(); it != nodes_.end(); ++it)
            processors.push_back((*it)->getProcessor());
    } else {
        // Insert all nodes into a vector
        //
        std::vector<GraphNode*> nodes;
        for (NodeSet::const_iterator it = nodes_.begin(); it != nodes_.end(); ++it)
            nodes.push_back((*it));

        // Sort that vector according to the IDs of the nodes.
        //
        std::sort(nodes.begin(), nodes.end(), NodeIDComparator());

        // Take the processors from the nodes in the sorted order.
        //
        for (size_t i = 0; i < nodes.size(); ++i)
            if (std::find(processors.begin(), processors.end(), nodes[i]->getProcessor()) == processors.end())
                processors.push_back(nodes[i]->getProcessor());
    }
    return processors;
}

std::vector<NetworkGraph::GraphNode*> NetworkGraph::getGraphNodes() const {
    std::vector<NetworkGraph::GraphNode*> nodes;

    for (NodeSet::const_iterator it = nodes_.begin(); it != nodes_.end(); ++it)
            nodes.push_back((*it));

    return nodes;
}

NetworkGraph* NetworkGraph::getSubGraph(const PortTypeCheck& type, const PortTypeCheck& loopType,
                                       const bool keepNodes) const {
    NetworkGraph* subGraph = new NetworkGraph(getProcessors(), type, loopType, keepNodes);
    return subGraph;
}

NetworkGraph* NetworkGraph::getTransposed() const {
    // Copy all nodes to a new graph. The successors will not be set for
    // new nodes yet.
    //
    NetworkGraph* transposed = new NetworkGraph(nodes_);
    for (NodeSet::const_iterator it = nodes_.begin(); it != nodes_.end(); ++it) {
        const GraphNode* const n = *it;
        const std::vector<GraphNode*>& successors = n->getSuccessors();

        // Find that current node in the transposed graph and add itself as
        // successor to all of its succesors...
        //
        GraphNode* const nt = transposed->findNode(n->getID());
        if (!nt) {
            LWARNING("Graph node with id " << n->getID() << " not present in the transposed graph");
            continue;
        }

        for (size_t s = 0; s < successors.size(); ++s) {
            const GraphNode* succ = successors[s];
            GraphNode* const st = transposed->findNode(succ->getID());
            if (st)
                st->addSuccessor(nt);
            else {
                LWARNING("Graph node with id " << succ->getID() << " not present in the transposed graph");
            }
        }
    }
    return transposed;
}

std::set<Processor*> NetworkGraph::getPredecessors(const std::set<Processor*>& processors) const {

    // generate transposed graph, if not present
    if (!networkGraphTransposed_)
        networkGraphTransposed_ = getTransposed();
    tgtAssert(networkGraphTransposed_, "Transposed network graph not generated");

    // perform successor search on transposed graph
    return networkGraphTransposed_->getSuccessors(processors);
}

std::set<Processor*> NetworkGraph::getPredecessors(const std::vector<Processor*>& processors) const {
    std::set<Processor*> proSet;
    std::copy(processors.begin(), processors.end(), std::inserter(proSet,proSet.begin()));
    return getPredecessors(proSet);
}

std::set<Processor*> NetworkGraph::getPredecessors(Processor* processor) const {
    tgtAssert(processor, "null pointer passed");
    std::set<Processor*> procSet;
    procSet.insert(processor);
    return getPredecessors(procSet);
}

std::set<Processor*> NetworkGraph::getPredecessors(Port* inport) const {
    tgtAssert(inport, "null pointer passed");
    tgtAssert(inport->isInport(), "passed port is not an inport");
    Processor* ownerProc = inport->getProcessor();
    tgtAssert(containsProcessor(ownerProc), "passed port's processor is not part of the network graph");

    // get all direct predecessors of the port
    std::set<Processor*> directPortPredecessors;
    const std::vector<const Port*> connectedPorts = inport->getConnected();
    for (size_t i=0; i<connectedPorts.size(); i++)
        directPortPredecessors.insert(connectedPorts.at(i)->getProcessor());

    // remove all direct predecessors that are not connected in the network graph
    std::set<Processor*> procPredecessors = getPredecessors(ownerProc);
    std::set<Processor*>::iterator it = directPortPredecessors.begin();
    while (it != directPortPredecessors.end()) {
        if (!procPredecessors.count(*it)) {
            directPortPredecessors.erase(it++);
        } else {
            ++it;
        }
    }

    // get predecessors of remaining direct port predecessors
    return getPredecessors(directPortPredecessors);
}

std::set<Processor*> NetworkGraph::getSuccessors(const std::set<Processor*>& processors) const {

    std::set<Processor*> result;
    CollectSuccessorsVisitor csv;

    setAllNodesMarked(false);
    std::set<Processor*>::const_iterator procIter;
    for (procIter = processors.begin(); procIter != processors.end(); ++procIter) {
        traverseBreadthFirst(*procIter, &csv);
    }

    std::set<NetworkGraph::GraphNode*>::const_iterator nodeIter;
    for (nodeIter = csv.getSuccessors().begin(); nodeIter != csv.getSuccessors().end(); ++nodeIter) {
        if (result.find((*nodeIter)->getProcessor()) == result.end())
            result.insert((*nodeIter)->getProcessor());
    }

    return result;
}

std::set<Processor*> NetworkGraph::getSuccessors(Processor* processor) const {
    tgtAssert(processor, "null pointer passed");
    std::set<Processor*> procSet;
    procSet.insert(processor);
    return getSuccessors(procSet);
}

std::set<Processor*> NetworkGraph::getSuccessors(Port* outport) const {
    tgtAssert(outport, "null pointer passed");
    tgtAssert(outport->isOutport(), "passed port is not an outport");
    Processor* ownerProc = outport->getProcessor();
    tgtAssert(containsProcessor(ownerProc), "passed port's processor is not part of the network graph");

    // get all direct successors of the port
    std::set<Processor*> directPortSuccessors;
    const std::vector<const Port*> connectedPorts = outport->getConnected();
    for (size_t i=0; i<connectedPorts.size(); i++)
        directPortSuccessors.insert(connectedPorts.at(i)->getProcessor());

    // remove all direct successors that are not connected in the network graph
    std::set<Processor*> procSuccessors = getSuccessors(ownerProc);
    std::set<Processor*>::iterator it = directPortSuccessors.begin();
    while (it != directPortSuccessors.end()) {
        if (!procSuccessors.count(*it)) {
            directPortSuccessors.erase(it++);
        } else {
            ++it;
        }
    }

    // get successors of remaining direct port predecessors
    return getSuccessors(directPortSuccessors);
}

std::vector<Processor*> NetworkGraph::sortTopologically(const std::set<Processor*>& processorSubset) const {
    // Perfrom depth-first search at first
    //
    fullTraverseDepthFirst();

    // Copy the nodes to a vector so that it can be sorted.
    //
    std::vector<GraphNode*> nodesVector;
    for (NodeSet::const_iterator it = nodes_.begin(); it != nodes_.end(); ++it)
        nodesVector.push_back(*it);

    // Sort the nodes according to their 'finish' timestep in descending order.
    //
    std::sort(nodesVector.begin(), nodesVector.end(), NodeSortTopologicalComparator());

    std::vector<Processor*> processors;

    // If the processor subset is not empty, perform some addtional tests.
    //
    if (processorSubset.empty() == false) {
        for (size_t i = 0; i < nodesVector.size(); ++i) {
            // Find the processor of the current node within the subset. If it is present,
            // add that processor to the vector which is about to be returned.
            //
            std::set<Processor*>::const_iterator it = processorSubset.find(nodesVector[i]->getProcessor());
            if (it != processorSubset.end())
                processors.push_back(nodesVector[i]->getProcessor());
        }
    } else {
        for (size_t i = 0; i < nodesVector.size(); ++i)
            processors.push_back(nodesVector[i]->getProcessor());
    }
    return processors;
}

bool NetworkGraph::traverseBreadthFirst(Processor* processor, GraphVisitor* const visitor) const {

    std::set<GraphNode*> roots = findNodes(processor);
    if (roots.empty())
        return false;

    setAllNodesMarked(false);
    std::set<GraphNode*>::iterator iter;
    for (iter = roots.begin(); iter != roots.end(); ++iter)
        traverseBreadthFirstInternal(*iter, visitor);
    return true;
}

bool NetworkGraph::traverseDepthFirst(Processor* processor,
                                      GraphVisitor* const onDiscoveryVisitor,
                                      GraphVisitor* const onFinishVisitor) const
{
    std::set<GraphNode*> roots = findNodes(processor);
    if (roots.empty())
        return false;

    setAllNodesMarked(false);
    std::set<GraphNode*>::iterator iter;
    setAllNodesDiscovered(0);
    setAllNodesFinished(0);
    for (iter = roots.begin(); iter != roots.end(); ++iter)
        traverseDepthFirstInternal(*iter, 0, onDiscoveryVisitor, onFinishVisitor);
    return true;
}

// private methods
//

NetworkGraph::NetworkGraph()
    : nodes_(),
      nodeCounter_(0),
      networkGraphTransposed_(0)
{
}

NetworkGraph::NetworkGraph(const NodeSet& nodes)
    : nodes_(),
      nodeCounter_(0),
      networkGraphTransposed_(0)
{
    for (NodeSet::const_iterator it = nodes.begin(); it != nodes.end(); ++it) {
        // Copy that node in the node set by making a 'deep' copy. Note that
        // the new node will have no successors set yet!
        //
        const GraphNode* const n = *it;
        std::pair<NodeSet::iterator, bool> res =
            nodes_.insert(new GraphNode(n->getProcessor(), n->getID()));

        // Increase counter if a new node was created.
        //
        if (res.second == true)
            ++nodeCounter_;
    }
}

NetworkGraph::GraphNode* NetworkGraph::addNode(Processor* const processor, bool* const added)
{

    if (networkGraphTransposed_) {
        LWARNING("Transposed network graph already generated");
        delete networkGraphTransposed_;
        networkGraphTransposed_ = 0;
    }

    if (processor == 0) {
        if (added != 0)
            *added = false;
        return 0;
    }

    GraphNode* node = new GraphNode(processor, static_cast<int>(nodeCounter_++));
    std::pair<NodeSet::iterator, bool> res = nodes_.insert(node);

    if (added != 0)
        *added = res.second;

    // If that node was already present in this graph, delete this new
    // and superfluous node.
    //
    if (res.second == false) {
        --nodeCounter_;
        delete node;
        return (*(res.first));
    }

    return node;
}

void NetworkGraph::connectNodes(const PortTypeCheck& ptc) {

    if (networkGraphTransposed_) {
        LWARNING("Transposed network graph already generated");
        delete networkGraphTransposed_;
        networkGraphTransposed_ = 0;
    }

    // clear successors and in/out degrees of all nodes
    for (NodeSet::const_iterator it = nodes_.begin(); it != nodes_.end(); ++it) {
        GraphNode* const n = *it;
        n->clearSuccessors();
        n->inDegree_ = 0;
    }

    // connect nodes
    for (NodeSet::const_iterator it = nodes_.begin(); it != nodes_.end(); ++it) {
        GraphNode* const n1 = *it;
        Processor* const p1 = n1->getProcessor();
        tgtAssert(p1, "GraphNode without processor");

        // Get all outports, including co-processor outports
        //
        std::vector<Port*> outports = p1->getOutports();
        const std::vector<CoProcessorPort*>& coOutports = p1->getCoProcessorOutports();
        outports.insert(outports.end(), coOutports.begin(), coOutports.end());

        int outportCounter = 0;
        for (size_t p = 0; p < outports.size(); ++p) {
            if (!ptc.isA(outports[p]))
                continue;

            ++outportCounter;
            const std::vector<const Port*> connectedPorts = outports[p]->getConnected();
            for (size_t c = 0; c < connectedPorts.size(); ++c) {
                const Port* const outport = connectedPorts[c];
                if (!ptc.isA(outport))
                    continue;

                // Find the corresponding node in the graph for that succeeding processor
                //
                Processor* p2 = outport->getProcessor();
                std::set<GraphNode*> nodes2 = findNodes(p2);
                if (nodes2.empty())
                    LERROR("buildGraph(): processor '" << p2->getID() << "' should already have a node, too!");
                else {
                    if (nodes2.size() > 1)
                        LWARNING("connect nodes: more than one graph node for processor \"" << p2->getID() << "\"");
                    std::set<GraphNode*>::iterator iter;
                    for (iter = nodes2.begin(); iter != nodes2.end(); ++iter)
                        n1->addSuccessor(*iter);
                }
            }   // for (c
        }   // for (p
    }   // for (it  (processors)
}

void NetworkGraph::unrollLoops(const PortTypeCheck& loopType) {

    if (networkGraphTransposed_) {
        LWARNING("Transposed network graph already generated");
        delete networkGraphTransposed_;
        networkGraphTransposed_ = 0;
    }

    // vector containing the loop initiator and loop finalizer
    // of all valid loops in the graph
    std::vector<std::pair<const Port*, Port*> > loops;

    // iterate over all graph nodes
    for (NodeSet::const_iterator it = nodes_.begin(); it != nodes_.end(); ++it) {
        GraphNode* const n1 = *it;
        Processor* const p1 = n1->getProcessor();
        tgtAssert(p1, "GraphNode without processor");

        // Get all outports, including co-processor outports
        //
        std::vector<Port*> outports = p1->getOutports();
        const std::vector<CoProcessorPort*>& coOutports = p1->getCoProcessorOutports();
        outports.insert(outports.end(), coOutports.begin(), coOutports.end());

        for (size_t p = 0; p < outports.size(); ++p) {

            if (!loopType.isA(outports[p]))
                continue;

            if (!outports[p]->isConnected()) {
                continue;
            }

            Port* destPort = outports[p];
            const Port* srcPort = destPort->getConnected().front();
            Processor* srcProc = srcPort->getProcessor();
            Processor* destProc = destPort->getProcessor();

            // check loop ports for multiple connections
            if (destPort->getConnected().size() > 1) {
                LWARNING("Loop port \"" << destPort->getID() << "\" of processor \"" <<
                         destProc->getID() << "\" has multiple connections. Skipping.");
                continue;
            }
            if (srcPort->getConnected().size() > 1) {
                LWARNING("Loop port \"" << srcPort->getID() << "\" of processor \"" <<
                         srcProc->getID() << "\" has multiple connections. Skipping.");
                continue;
            }

            // check loop processors for multiple loop ports
            std::vector<Port*> srcPorts = srcProc->getPorts();
            for (size_t i=0; i<srcPorts.size(); ++i) {
                if (loopType.isA(srcPorts[i]) && (srcPorts[i] != srcPort)) {
                    LWARNING("Processor \"" << srcProc->getID() << "\" has multiple loop ports. Skipping.");
                    continue;
                }
            }
            std::vector<Port*> destPorts = destProc->getPorts();
            for (size_t i=0; i<destPorts.size(); ++i) {
                if (loopType.isA(destPorts[i]) && (destPorts[i] != destPort)) {
                    LWARNING("Processor \"" << destProc->getID() << "\" has multiple loop ports. Skipping.");
                    continue;
                }
            }

            // check whether loop end proc is successor of loop start proc
            if (!isSuccessor(srcProc, destProc)) {
                LWARNING("Skipping loop: \"" << destProc->getID() << "\" is not a successor of \"" << srcProc->getID() << "\"");
                continue;
            }

            loops.push_back(std::pair<const Port*, Port*>(srcPort, destPort));

        }   // for (p
    }   // for (it  (processors)

    // discard overlapping loops (nested loops are allowed)
    std::set<int> overlappings;
    for (size_t i=0; i<loops.size(); ++i) {
        if (overlappings.find(static_cast<const int>(i)) != overlappings.end())
            continue;
        for (size_t j=0; j<loops.size(); ++j) {
            if (i==j)
                continue;
            if (isPathElement(loops[j].first->getProcessor(), loops[i].first->getProcessor(), loops[i].second->getProcessor()) !=
                isPathElement(loops[j].second->getProcessor(), loops[i].first->getProcessor(), loops[i].second->getProcessor())  ) {
                LWARNING("Loops \"" << loops[i].first->getProcessor()->getID() << "\" -> \"" << loops[i].second->getProcessor()->getID() << "\" and \"" <<
                          loops[j].first->getProcessor()->getID() << "\" -> \"" << loops[j].second->getProcessor()->getID() << "\" overlap. Skipping.");
                overlappings.insert(static_cast<int>(i));
                overlappings.insert(static_cast<int>(j));
            }
        }
    }
    std::vector<std::pair<const Port*, Port*> > validLoops;
    for (size_t i=0; i<loops.size(); ++i) {
        if (overlappings.find(static_cast<int>(i)) == overlappings.end())
            validLoops.push_back(loops[i]);
    }

    //
    // now we are left with valid loops => unroll
    //
    bool finished = false;
    while (!finished) {
        finished = true;
        for (size_t i=0; i<validLoops.size(); ++i) {

            Processor* loopSrc = validLoops[i].first->getProcessor();
            Processor* loopDest = validLoops[i].second->getProcessor();

            // todo check nested
            bool hasNested = false;
            for (size_t j=0; j<validLoops.size() && !hasNested; ++j) {
                if (i==j)
                    continue;
                if (isPathElement(validLoops[j].first->getProcessor(), loopSrc, loopDest))
                    hasNested = true;
            }
            if (hasNested)
                continue;

            int numIterations = loopType.getNumIterations(validLoops[i].second);

            LDEBUG("Unrolling loop \"" << loopSrc->getID() << "\" -> \"" << loopDest->getID() << "\" (" << numIterations << " iterations)");
            validLoops.erase(validLoops.begin() + i);

            // detect path nodes
            std::set<GraphNode*> loopSrcNodes = findNodes(loopSrc);
            if (loopSrcNodes.size() != 1) {
                LERROR("Exactly one loop source node expected. Found: " << loopSrcNodes.size() << ". Skipping");
                continue;
            }
            GraphNode* loopSrcNode = *loopSrcNodes.begin();
            std::set<GraphNode*> loopDestNodes = findNodes(loopDest);
            if (loopDestNodes.size() != 1) {
                LERROR("Exactly one loop dest node expected. Found: " << loopDestNodes.size() << ". Skipping");
                continue;
            }
            GraphNode* loopDestNode = *loopDestNodes.begin();

            std::set<GraphNode*> pathNodes = getPathNodes(*loopSrcNodes.begin(), *loopDestNodes.begin());
            tgtAssert(pathNodes.find(loopSrcNode) != pathNodes.end(), "Loop source node not contained by path nodes");
            tgtAssert(pathNodes.find(loopDestNode) != pathNodes.end(), "Loop dest node not contained by path nodes");

            GraphNode* prevLoopDestNode = loopDestNode;
            for (int iteration=1; iteration<numIterations; ++iteration) {
                std::set<GraphNode*> duplicateNodes;
                std::map<GraphNode*,GraphNode*> duplicateMap;
                for (std::set<GraphNode*>::iterator iter = pathNodes.begin(); iter!=pathNodes.end(); ++iter) {
                    GraphNode* duplicate = new GraphNode(*(*iter));
                    duplicateMap[*iter] = duplicate;
                    duplicate->id_ = static_cast<int>(nodeCounter_++);
                    nodes_.insert(duplicate);
                    duplicateNodes.insert(duplicate);
                }
                for (std::set<GraphNode*>::iterator iter = pathNodes.begin(); iter!=pathNodes.end(); ++iter) {
                    GraphNode* original = *iter;
                    GraphNode* duplicate = duplicateMap[original];
                    for (size_t succ=0; succ<original->successors_.size(); ++succ) {
                        if (duplicateMap.find(original->successors_[succ]) != duplicateMap.end())
                            duplicate->addSuccessor(duplicateMap[original->successors_[succ]]);
                        else
                            duplicate->addSuccessor(original->successors_[succ]);
                    }
                }
                tgtAssert(duplicateMap.find(loopSrcNode) != duplicateMap.end(), "Loop source node has no key in duplicate map");
                prevLoopDestNode->addSuccessor(duplicateMap[loopSrcNode]);
                prevLoopDestNode = duplicateMap[loopDestNode];
            }

            finished = false;
            break;
        }
    }
    if (!validLoops.empty())
        LWARNING("Valid loops vector should be empty at this point");

}

void NetworkGraph::buildGraph(const std::vector<Processor*>& processors, const PortTypeCheck& ptc,
                              const PortTypeCheck& loopPtc, const bool keepNodes)
{
    clearGraph();

    // Add new nodes for all given processors to the graph
    //
    for (size_t i = 0; i < processors.size(); ++i) {
        if ((keepNodes == false) && !ptc.hasA(processors[i]))
            continue;

        addNode(processors[i]);
    }

    connectNodes(ptc);

    unrollLoops(loopPtc);
}

void NetworkGraph::clearGraph() {

    delete networkGraphTransposed_;
    networkGraphTransposed_ = 0;

    for (NodeSet::iterator it = nodes_.begin(); it != nodes_.end(); ++it)
        delete (*it);
    nodes_.clear();
    nodeCounter_ = 0;
}

void NetworkGraph::copy(const NetworkGraph& other) {
    clearGraph();

    for (NodeSet::const_iterator it = other.nodes_.begin(); it != other.nodes_.end(); ++it)
        nodes_.insert( new GraphNode(*(*it)) );

    // The successors are NOT copied by GraphNode, because they would have the address
    // in the other graph. So the successors must be re-inserted with the addresses they
    // have obtained by being copied to the new graph.
    //
    for (NodeSet::const_iterator it = other.nodes_.begin(); it != other.nodes_.end(); ++it) {
        const std::vector<GraphNode*>& successors = (*it)->getSuccessors();
        GraphNode* n = findNode((*it)->getID());
        if (n == 0) {
            LWARNING("Graph node with id " << (*it)->getID() << " not present in the copied graph");
            continue;
        }
        for (size_t i = 0; i < successors.size(); ++i) {
            const GraphNode* succ = successors[i];
            GraphNode* s = findNode(succ->getID());
            if (s != 0) {
                n->addSuccessor(s);
            }
            else {
                LWARNING("Graph node with id " << succ->getID() << " not present in the copied graph");
            }
        }
    }

    // finally copy over node counter
    nodeCounter_ = other.nodeCounter_;
}

void NetworkGraph::pruneGraph(const std::set<GraphNode*>& nodes) {

    // collect and remove all graph nodes that are not an element of 'nodes'
    std::vector<GraphNode*> toDelete;
    NodeSet::iterator iter = nodes_.begin();
    while (iter != nodes_.end()) {
        if (nodes.find(*iter) == nodes.end()) {
            NodeSet::iterator tIter = iter;
            iter++;   // increment before erasing!
            toDelete.push_back(*tIter);
            nodes_.erase(tIter);
        }
        else
            iter++;
    }

    // remove all successor pointers that reference an removed node
    for (iter = nodes_.begin(); iter != nodes_.end(); ++iter) {
        GraphNode* node = *iter;
        std::vector<GraphNode*>::iterator sucIter = node->successors_.begin();
        while (sucIter != node->successors_.end()) {
            if (!findNode((*sucIter)->getID()))
                sucIter = node->successors_.erase(sucIter);
            else
                sucIter++;
        }
    }

    // update nodeCounter_
    nodeCounter_ = nodes_.size();

    // update outDegree and reset inDegree of all remaining nodes
    for (iter = nodes_.begin(); iter != nodes_.end(); ++iter) {
        GraphNode* node = *iter;
        node->inDegree_ = 0;
        node->outDegree_ = static_cast<int>(node->successors_.size());
    }
    // update inDegree
    for (iter = nodes_.begin(); iter != nodes_.end(); ++iter) {
        GraphNode* node = *iter;
        for (size_t i=0; i<node->successors_.size(); ++i)
            node->successors_[i]->inDegree_++;
    }

    // finally delete all pruned nodes
    for (size_t i=0; i<toDelete.size(); ++i)
        delete toDelete[i];

}

std::set<NetworkGraph::GraphNode*> NetworkGraph::findNodes(Processor* processor) const {

    std::set<NetworkGraph::GraphNode*> result;

    NodeSet::const_iterator it;
    for (it = nodes_.begin(); it != nodes_.end(); ++it){
        if ((*it)->getProcessor() == processor)
            result.insert(*it);
    }

    return result;
}

NetworkGraph::GraphNode* NetworkGraph::findNode(const int id) const {
    if (id < 0)
        return 0;

    GraphNode tmp(id);      // node with no processor, i.e. NULL pointer
    NodeSet::const_iterator it = nodes_.find(&tmp);
    if (it == nodes_.end())
        return 0;
    else
        return (*it);
}

std::set<NetworkGraph::GraphNode*> NetworkGraph::getSuccessorNodes(const std::set<NetworkGraph::GraphNode*>& nodes) const {

    CollectSuccessorsVisitor csv;

    setAllNodesMarked(false);
    std::set<NetworkGraph::GraphNode*>::const_iterator iter;
    for (iter = nodes.begin(); iter != nodes.end(); ++iter) {
        // detect node in graph with equal id
        GraphNode* node = findNode((*iter)->getID());
        if (node)
            traverseBreadthFirstInternal(node, &csv);
    }

    return csv.getSuccessors();
}

std::set<NetworkGraph::GraphNode*> NetworkGraph::getPredecessorNodes(const std::set<NetworkGraph::GraphNode*>& nodes) const {

    // generate transposed graph, if not present
    if (!networkGraphTransposed_)
        networkGraphTransposed_ = getTransposed();
    tgtAssert(networkGraphTransposed_, "Transposed network graph not generated");

    // perform successor search on transposed graph
    std::set<GraphNode*> nodesTransposed;
    for (std::set<GraphNode*>::const_iterator iter=nodes.begin(); iter != nodes.end(); ++iter) {
        GraphNode* node = networkGraphTransposed_->findNode((*iter)->getID());
        if (node)
            nodesTransposed.insert(node);
    }
    std::set<GraphNode*> successorsTransposed = networkGraphTransposed_->getSuccessorNodes(nodesTransposed);

    // transform result back into original graph
    std::set<GraphNode*> result;
    for (std::set<GraphNode*>::iterator iter=successorsTransposed.begin(); iter != successorsTransposed.end(); ++iter) {
        GraphNode* node = findNode((*iter)->getID());
        if (node)
            result.insert(node);
    }

    return result;
}

std::set<NetworkGraph::GraphNode*> NetworkGraph::getPathNodes(const GraphNode* predecessor, const GraphNode* successor) const {

    NetworkGraph pathGraph(*this);

    // detect successor nodes of predecessor
    std::set<GraphNode*> successorNodes;
    successorNodes.insert(const_cast<GraphNode*>(predecessor));
    successorNodes = pathGraph.getSuccessorNodes(successorNodes);

    // detect successor nodes
    std::set<GraphNode*> predecessorNodes;
    predecessorNodes.insert(const_cast<GraphNode*>(successor));
    predecessorNodes = pathGraph.getPredecessorNodes(predecessorNodes);

    // reduce graph to intersection of predecessor and successor nodes
    pathGraph.pruneGraph(successorNodes);
    pathGraph.pruneGraph(predecessorNodes);

    // returns nodes with sames ids as the nodes contained by the pruned graph
    std::set<GraphNode*> result;
    for (NodeSet::iterator iter = pathGraph.nodes_.begin(); iter != pathGraph.nodes_.end(); ++iter) {
        GraphNode* node = findNode((*iter)->getID());
        if (node)
            result.insert(node);
    }
    return result;
}

std::vector<NetworkGraph::GraphNode*> NetworkGraph::identifyRoots() const {
    std::vector<GraphNode*> roots;
    for (NodeSet::const_iterator it = nodes_.begin(); it != nodes_.end(); ++it) {
        if ((*it)->getInDegree() == 0)
            roots.push_back(*it);
    }
    return roots;
}

void NetworkGraph::setAllNodesMarked(const bool marked) const {
    for (NodeSet::const_iterator it = nodes_.begin(); it != nodes_.end(); ++it) {
        GraphNode* const node = const_cast<GraphNode*>(*it);
        node->setMarked(marked);
    }
}

void NetworkGraph::setAllNodesDiscovered(const int value) const {
    for (NodeSet::const_iterator it = nodes_.begin(); it != nodes_.end(); ++it) {
        GraphNode* const node = const_cast<GraphNode*>(*it);
        node->setDiscovered(value);
    }
}

void NetworkGraph::setAllNodesFinished(const int value) const {
    for (NodeSet::const_iterator it = nodes_.begin(); it != nodes_.end(); ++it) {
        GraphNode* const node = const_cast<GraphNode*>(*it);
        node->setFinished(value);
    }
}

void NetworkGraph::traverseBreadthFirstInternal(GraphNode* const root,
                                                GraphVisitor* const visitor) const
{
    if (root == 0)
        return;

    std::deque<GraphNode*> q;
    root->mark();
    q.push_back(root);

    while (q.empty() == false) {
        GraphNode* const n = q.front();
        q.pop_front();

        // Apply visitor if present
        //
        if (visitor != 0)
            visitor->visit(n);

        const std::vector<GraphNode*>& successors = n->getSuccessors();
        const int level = n->minDistanceToRoot_;
        for (size_t i = 0; i < successors.size(); ++i) {
            GraphNode* const n = successors[i];
            if (n->isMarked() == false) {
                n->mark();
                q.push_back(n);
                n->minDistanceToRoot_ = (level + 1);
            }
        }
    }   // while
}

int NetworkGraph::traverseDepthFirstInternal(GraphNode* const root, const int startTime,
                                             GraphVisitor* const onDiscoveryVisitor,
                                             GraphVisitor* const onFinishVisitor) const
{
    int time = startTime;
    std::stack<GraphNode*> stack;
    stack.push(root);
    root->setDiscovered(++time);
    if (onDiscoveryVisitor != 0)
        onDiscoveryVisitor->visit(root);

    while (stack.empty() == false) {
        GraphNode* const n = stack.top();
        const std::vector<GraphNode*>& successors = n->getSuccessors();
        bool newDiscovery = false;
        for (size_t j = 0; j < successors.size(); ++j) {
            GraphNode* succ = successors[j];
            if (succ->isDiscovered() == false) {
                succ->setDiscovered(++time);

                // Apply discovery visitor, if present.
                //
                if (onDiscoveryVisitor != 0)
                    onDiscoveryVisitor->visit(succ);

                stack.push(succ);
                newDiscovery = true;
                break;
            }
        }

        if (newDiscovery == false) {
            stack.pop();
            n->setFinished(++time);

            // Apply finish visitor, if present.
            //
            if (onFinishVisitor != 0)
                onFinishVisitor->visit(n);
        }
    }   // while
    return time;
}

}   // namespace
