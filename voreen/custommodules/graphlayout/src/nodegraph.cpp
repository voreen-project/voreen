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

#include "../include/nodegraph.h"

#include "tgt/assert.h"

namespace voreen {

// - NodeGraphNode implementation -----------------------------------------------------------------

NodeGraphNode::NodeGraphNode(plot_t x, plot_t y, plot_t mass)
    : position_(x, y)
    , velocity_(0, 0)
    , mass_(mass)
{
}

void NodeGraphNode::connect(const NodeGraphNode* otherNode) {
    if (otherNode != 0) {
        connections_.insert(otherNode);
    }
}

void NodeGraphNode::disconnect(const NodeGraphNode* otherNode) {
    if (otherNode != 0) {
        connections_.erase(otherNode);
    }
}

bool NodeGraphNode::isConnectedTo(const NodeGraphNode* otherNode) const {
    return ((otherNode != 0) && (connections_.find(otherNode) != connections_.end()));
}

std::set<const NodeGraphNode*>::const_iterator NodeGraphNode::getConnectionsBegin() const {
    return connections_.begin();
}

std::set<const NodeGraphNode*>::const_iterator NodeGraphNode::getConnectionsEnd() const {
    return connections_.end();
}

size_t NodeGraphNode::getConnectionCount() const {
    return connections_.size();
}

// - NodeGraph implementation ---------------------------------------------------------------------

NodeGraph::NodeGraph() {
}

NodeGraph::~NodeGraph() {
    clearNodes();
}


NodeGraphNode* NodeGraph::addNode(int id, plot_t x, plot_t y, plot_t mass) {
    std::map<int, NodeGraphNode*>::iterator it = nodes_.find(id);
    if (it != nodes_.end()) {
        it->second->position_.x = x;
        it->second->position_.y = y;
        it->second->mass_ = mass;
        return it->second;
    }
    else {
        NodeGraphNode* newNode = new NodeGraphNode(x, y, mass);
        nodes_.insert(std::make_pair(id, newNode));
        return newNode;
    }
}

void NodeGraph::connectNodes(int first, int second) {
    std::map<int, NodeGraphNode*>::iterator it1 = nodes_.find(first);
    std::map<int, NodeGraphNode*>::iterator it2 = nodes_.find(second);
    if (it1 != nodes_.end() && it2 != nodes_.end()) {
        it1->second->connect(it2->second);
        it2->second->connect(it1->second);
    }
}

void NodeGraph::connectNodes(NodeGraphNode* first, NodeGraphNode* second) {
    if (first != 0 && second != 0) {
        first->connect(second);
        second->connect(first);
    }
}

NodeGraphNode* NodeGraph::getNode(int id) {
    std::map<int, NodeGraphNode*>::iterator it = nodes_.find(id);
    if (it != nodes_.end()) {
        return it->second;
    }
    else {
        return 0;
    }
}

void NodeGraph::clearNodes() {
    for (std::map<int, NodeGraphNode*>::iterator it = nodes_.begin(); it != nodes_.end(); ++it) {
        delete it->second;
    }
    nodes_.clear();
}

size_t NodeGraph::getNodesCount() const {
    return nodes_.size();
}

std::map<int, NodeGraphNode*>::const_iterator NodeGraph::getNodesBegin() const {
    return nodes_.begin();
}

std::map<int, NodeGraphNode*>::const_iterator NodeGraph::getNodesEnd() const {
    return nodes_.end();
}

} // namespace voreen
