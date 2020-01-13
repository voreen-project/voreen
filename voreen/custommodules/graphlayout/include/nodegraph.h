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

#ifndef VRN_NODEGRAPH_H
#define VRN_NODEGRAPH_H

#include "modules/plotting/datastructures/interval.h"
#include "modules/plotting/datastructures/plotbase.h"
#include "tgt/vector.h"

#include <map>
#include <set>

namespace voreen {

/**
 * Represents a single node of a NodeGraph.
 **/
struct NodeGraphNode {
public:
    /**
     * Constructs a new node at given coordiantes
     *
     * \param   x       x coordinate
     * \param   y       y coordinate
     * \param   mass    mass
     **/
    NodeGraphNode(plot_t x, plot_t y, plot_t mass);

    /**
     * Connects this node with \a otherNode.
     *
     * \param   \otherNode  other node to connect with
     **/
    void connect(const NodeGraphNode* otherNode);

    /**
     * Disconnects this node from \a otherNode.
     *
     * \param   \otherNode  other node to disconnect from
     **/
    void disconnect(const NodeGraphNode* otherNode);

    /**
     * Returns true, if this node is connected to \a otherNode
     *
     * \param   \otherNode  other node check for connection
     **/
    bool isConnectedTo(const NodeGraphNode* otherNode) const;

    /**
     * Returns an iterator to the begin of all connections.
     */
    std::set<const NodeGraphNode*>::const_iterator getConnectionsBegin() const;
    /**
     * Returns an iterator to the end of all connections.
     */
    std::set<const NodeGraphNode*>::const_iterator getConnectionsEnd() const;
    /**
     * Returns the number of nodes connected to this node
     **/
    size_t getConnectionCount() const;

    tgt::dvec2 position_;   ///< position of the node
    tgt::dvec2 velocity_;   ///< velocity/force of the node

    plot_t mass_;       ///< mass of the node

private:
    std::set<const NodeGraphNode*> connections_;      ///< set of all connections to other graph nodes
};


/**
 * Represents a graph of nodes and connections
 **/
class NodeGraph {
public:
    /**
     * Constructs a new empty node graph.
     **/
    NodeGraph();

    /**
     * Default destructor, deletes all held NodeGraphNodes.
     **/
    virtual ~NodeGraph();

    /**
     * Adds a new node at given coordinates to the graph
     *
     * \param   id      id of the node
     * \param   x       x coordinate
     * \param   y       y coordinate
     * \param   mass    mass
     **/
    NodeGraphNode* addNode(int id, plot_t x, plot_t y, plot_t mass);

    /**
     * Connects the to nodes given by ids.
     *
     * \param   first   id of the first node to connect
     * \param   second  id of the second node to connect
     **/
    void connectNodes(int first, int second);

    /**
     * Connects the to nodes given by pointers.
     *
     * \param   first   pointer to the first node to connect
     * \param   second  pointer to the second node to connect
     **/
    void connectNodes(NodeGraphNode* first, NodeGraphNode* second);

    /**
     * Returns a pointer to the node with the given id.
     *
     * \param   id  id of the node to return
     **/
    NodeGraphNode* getNode(int id);

    /**
     * Removes all nodes and connections from this graph. All pointers and iterators will be invalidated!
     */
    void clearNodes();

    /**
     * Returns the number of nodes in this node graph.
     **/
    size_t getNodesCount() const;

    /**
     * Returns an iterator to the begin of all nodes.
     */
    std::map<int, NodeGraphNode*>::const_iterator getNodesBegin() const;

    /**
     * Returns an iterator to the end of all nodes.
     */
    std::map<int, NodeGraphNode*>::const_iterator getNodesEnd() const;

protected:
    std::map<int, NodeGraphNode*> nodes_;       ///< map of all contained nodes associated to their ids

};


}

#endif // VRN_NODEGRAPH_H

