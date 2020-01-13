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

#ifndef VRN_NETWORKGRAPH_H
#define VRN_NETWORKGRAPH_H

#include <set>
#include <typeinfo>
#include <vector>

#include "voreen/core/voreencoreapi.h"
#include "voreen/core/ports/port.h"

namespace voreen {

class VRN_CORE_API PortTypeCheck {
    public:
        virtual ~PortTypeCheck() {}
        virtual bool isA(const Port* p) const = 0;
        virtual bool hasA(const Processor* p) const = 0;
        virtual int getNumIterations(const Port* p) const = 0;
};

// ----------------------------------------------------------------------------

class VRN_CORE_API PortTypeCheckReject : public PortTypeCheck {
    public:
        virtual bool isA(const Port* /*p*/) const {
            return false;
        }
        virtual bool hasA(const Processor* /*p*/) const {
            return false;
        }
        int getNumIterations(const Port* /*p*/) const {
            return -1;
        }
};

// ----------------------------------------------------------------------------

class VRN_CORE_API LoopPortTypeCheck : public PortTypeCheck {
    public:
        LoopPortTypeCheck() :
            inverse_(false)
        {}

        LoopPortTypeCheck(bool inverse) :
            inverse_(inverse)
        {}
        virtual bool isA(const Port* p) const {
            if (inverse_)
                return !(p->isLoopPort());
            else
                return p->isLoopPort();
        }
        virtual bool hasA(const Processor* p) const {
            bool hasPort = false;
            for (size_t i=0; i<p->getPorts().size(); ++i) {
                if (inverse_) {
                    if (!p->getPorts()[i]->isLoopPort())
                        hasPort = true;
                }
                else {
                    if (p->getPorts()[i]->isLoopPort())
                        hasPort = true;
                }
            }
            return hasPort;
        }
        int getNumIterations(const Port* p) const {
            return std::max(p->getNumLoopIterations(), 1);
        }

    private:
        bool inverse_;
};

// ----------------------------------------------------------------------------

template <typename T>
class GenericPortTypeCheck : public PortTypeCheck {
    public:
        GenericPortTypeCheck() :
            inverse_(false)
        {}

        GenericPortTypeCheck(bool inverse) :
            inverse_(inverse)
        {}

        virtual bool isA(const Port* p) const {
            if (dynamic_cast<const T*>(p))
                return !inverse_;
            else
                return inverse_;
            //return dynamic_cast<T*>(p);
        }

        bool hasA(const Processor* p) const {
            std::vector<Port*> ports = p->getPorts();
            for (size_t i = 0; i < ports.size(); ++i) {
                if (isA(ports[i]))
                    return true;
            }
            return false;
        }
        int getNumIterations(const Port* p) const {
            return std::max(p->getNumLoopIterations(), 1);
        }

    private:
        bool inverse_;
};

// ----------------------------------------------------------------------------

/**
 * Class for representing processor networks as directed graphs. This class performs
 * all kinds of topological analysis, sorting and traversing within the networks.
 * It is primarily intended to be used by <code>NetworkEvaluator</code> for
 * determining the processing order and traversing the networks.
 * The visitor design-pattern is used to perform certain tasks during the graph
 * traversal. See also <code>GraphVisitor</code> class.
 *
 * @author  Dirk Feldmann, July/August 2009
 */
class VRN_CORE_API NetworkGraph {
public:

    /**
     * Class for representing a node within the Graph. The nodes store the related
     * pointer to the <code>Processor</code> object and their direct successors so
     * that they form adjacency lists. Each node yields an ID to identify them with
     * a graph.
     */
    class GraphNode {
    friend class NetworkGraph;
    public:
        /**
         * Comparator class allowing to compare the content of the passed GraphNode
         * pointers instead of only comporing the pointers with each other.
         */
        struct GraphNodeComparator {
            /**
             * Operator returning false if one or more of the passed pointers are NULL
             * or the result of (*n1 < *n2) otherwise.
             */
            bool operator()(const GraphNode* const n1, const GraphNode* const n2) const {
                if ((n1 == 0) || (n2 == 0))
                    return false;
                return (*n1 < *n2);
            }
        };

    public:
        /**
         * Ctor which uses the passed ID for the node.
         * This ctor is primary intended to be used by NetworkGraph class only
         * for transposing graphs.
         *
         * @param   processor   The pointer to the Processor object which shall
         *                      be wrapped by this GraphNode. This pointer should
         *                      not be NULL.
         * @param   id  ID of the node. Make sure that this ID is really unique
         *              within the graph containing this node.
         */
        explicit GraphNode(Processor* const processor, const int id);

        /**
         * Ctor creating a new node for the passed Processor objects. The
         * pointer shall not be NULL and the ID of the new node is set to the
         * invalid value -1.
         * This ctor is primary used by the NetworkGraph to search for nodes
         * by the processor pointer.
         *
         * @param   processor   The pointer to the Processor object which shall
         *                      be wrapped by this GraphNode. This pointer should
         *                      not be NULL.
         */
        explicit GraphNode(Processor* const processor);

        /**
         * Ctor creating a new node which holds only and ID, but no Processor pointer.
         * The ID should be unique in the graph using this node.
         * This ctor is primary used by the NetworkGraph to search for nodes
         * by ID.
         *
         * @param   id  ID of the node.
         */
        explicit GraphNode(int id);

        /**
         * Copy-ctor performing shallow copies (pointer to Processor objects)
         * and does not copy the successors of this node (see comment on copy-
         * ctor of <code>NetworkGraph</code>).
         */
        GraphNode(const GraphNode& other);

        /**
         * Dtor. Has nothing to do.
         */
        ~GraphNode();

        /**
         * Assignment operator performing shallow copies only and does not
         * copy the successors of this node (see comment on copy-
         * ctor of <code>NetworkGraph</code>).
         */
        GraphNode& operator =(const GraphNode& other);

        /**
         * A graph node is defined to be "less" than another,
         * if its id is less than the other's.
         */
        bool operator<(const GraphNode& other) const;

        /**
         * Two graph nodes are defined to be equal, if both
         * have the same id and wrap the same processor.
         * Thus, the operator can be used to identify nodes
         * across graphs.
         */
        bool operator==(const GraphNode& other) const;

        /**
         * Convenience operator. Returns !(*this == other).
         */
        bool operator!=(const GraphNode& other) const;

        /**
         * Adds the passed GraphNode pointer to its successor vector, if
         * it is not NULL. Otherwise nothing will happen.
         * By adding a successor, this method increase the out-degree of this
         * node the in-degree of the successor.
         *
         * @param   successor   Successor to be added. This pointer shall not be NULL.
         */
        void addSuccessor(GraphNode* const successor);

        /**
         * Removes all successors from this node. Its out-degree will be set to 0
         * and the in-degrees of all former successor will be decreased by 1.
         */
        void clearSuccessors();

        /**
         * Return the ID of this node.
         */
        int getID() const { return id_; }

        /**
         * Returns this in-degree of this node. This is the number of direct predecessors
         * and does necessarily correspond to the number of inports to the related
         * Processor object!
         */
        int getInDegree() const { return inDegree_; }

        /**
         * Returns the minimal distance to the root which has bee calculated during
         * the last breadth-first search. The term 'root' herein also referes to that
         * particular root used by the last breadth-first search. The default value
         * is 0, so if no breadth-first search has been performed in the graph containing
         * this node, all nodes will yield zero.
         */
        int getMinDistanceToRoot() const { return minDistanceToRoot_; }

        /**
         * Returns the out-degree of this node. This is the number of direct successors
         * nodes and does not necessarily correspond to the number of outports on the
         * related Processor object!
         */
        int getOutDegree() const { return outDegree_; }

        /**
         * Return the processor object wrapped by this node.
         */
        Processor* getProcessor() const { return processor_; }

        /**
         * Returns all successor nodes of this one.
         */
        const std::vector<GraphNode*>& getSuccessors() const { return successors_; }

        /**
         * Returns the value of an internal bool marker which is used by breadth-first
         * search.
         */
        bool isMarked() const { return marked_; }

        /**
         * Sets the values of the internal bool marker, which is e.g. used by
         * breadth-first search.
         */
        void setMarked(const bool marked) { marked_ = marked; }

        /**
         * Convenience method to set the internal bool marker to true.
         */
        void mark() { marked_ = true; }

        /**
         * Convenience method to set the internal bool marker to false.
         */
        void unmark() { marked_ = false; }

        /**
         * Return the 'timestamp' of discovery of this node during depth-first search.
         */
        int getDiscovered() const { return tDiscovered_; }

        /**
         * Return the 'timestamp' of finalization of this node during depth-first search.
         */
        int getFinished() const { return tFinished_; }

        /**
         * Convenience method to determine whether this node has already been discovered
         * during depth-first search.
         */
        bool isDiscovered() const { return (tDiscovered_ > 0); }

        /**
         * Convenience method to determine whether this node has already been finished
         * during depth-first search.
         */
        bool isFinished() const { return (tFinished_ > 0); }

        /**
         * Sets the timestamp for discovery of this node during depth-first search to
         * the given value. The timestamps need to be strictly monotonic increasing
         * and must be > 0.
         */
        void setDiscovered(const int value) { tDiscovered_ = value; }

        /**
         * Sets the timestamp for discovery of this node during depth-first search to
         * the given value. The timestamps need to be strictly monotonic increasing
         * and must be > 0.
         */
        void setFinished(const int value) { tFinished_ = value; }

    private:
        /**
         * Performs internal shallow copying for copy-ctor and operator=.
         * The successors of this node will not be copied. See comment on
         * <code>NetworkGraph::copy()</code>.
         */
        void copy(const GraphNode& other);

    private:
        /** Pointer to the wrapped processor object. */
        Processor* processor_;

        /** ID of this node. Must be unique within graph using it. */
        int id_;

        /** All direct successor nodes of this one. */
        std::vector<GraphNode*> successors_;

        /** Number of direct predecessor nodes. */
        int inDegree_;

        /** Number of direct successor nodes. */
        int outDegree_;

        /** Marker used by breadth-first search. */
        bool marked_;

        /** Minimal distance to the root determined during last breadth-first search. */
        int minDistanceToRoot_;

        /** Timestamp of discovery of this node during depth-first search. */
        int tDiscovered_;

        /** Timestamp of finalization of this node during depth-first search. */
        int tFinished_;
    };  // class GraphNode

    /**
     * Pure virtual base class for visitor design-pattern used by <code>NetworkGraph</code> class.
     * Visitors used by NetworkGraph have to implement this interface.
     */
    class GraphVisitor {
    public:
        virtual ~GraphVisitor() {}

        /**
         * Intended to perfom the visit on the given GraphNode object. The return
         * value shall indicate possible failures during the visit, e.g. a passed
         * NULL-pointer to a GrapNode object, etc..
         *
         * @return  true shall be returned if the visit has been successful, false
         *          otherwise.
         */
        virtual bool visit(NetworkGraph::GraphNode* const graphNode) = 0;
    };

    /**
     * GraphVisitor subclass designed to collect all visited processors during
     * graph traversing.
     */
    class CollectSuccessorsVisitor : public GraphVisitor {
    public:
        CollectSuccessorsVisitor();
        virtual ~CollectSuccessorsVisitor();
        const std::set<GraphNode*>& getSuccessors() const;
        const std::vector<GraphNode*>& getSuccessorsOrdered() const;
        virtual bool visit(GraphNode* const graphNode);
    private:
        std::set<GraphNode*> successors_;
        std::vector<GraphNode*> successorOrdered_;
    };  // class CollectSuccessorsVisitor

public:
    /**
     * Ctor contruction NetworkGraph with according to the given processors and
     * the connections between them.
     */
    explicit NetworkGraph(const std::vector<Processor*>& processors);

    /**
     * Ctor contruction NetworkGraph with according to the given processors and
     * the connections between them, but confines them to processors and connections
     * which use the given PortTypeCheck. Additionally, connections matching the
     * passed loop type criterion are considered as static loops and unrolled.
     *
     * @param   type    Type of ports and connections to which the graph shall be
     *                  confined. Acts like a kind of mask.
     * @param   loopType    Ports and connections of this type are considered as static loops
     *                      and are unrolled, i.e., all nodes in the graph representing
     *                      processors between the source and the destination of the loop edge
     *                      are duplicated and inserted after the loop destination. The number
     *                      of duplications of the loop body is also determined by loopType.
     * @param   keepNodes   Determines whether processors/nodes which are not affected
     *                      by the given PortTypeCheck shall be kept in the graph or
     *                      disposed.
     */
    explicit NetworkGraph(const std::vector<Processor*>& processors, const PortTypeCheck& type,
        const PortTypeCheck& loopType, bool keepNodes = false);

    /**
     * Copy-ctor. This methods calls <code>copy()</code>
     */
    NetworkGraph(const NetworkGraph& other);

    /**
     * Dtor. Calls <code>clearGraph()</code>.
     */
    ~NetworkGraph();

    /**
     * Assignment-operator internally calling <code>copy()</code>.
     */
    NetworkGraph& operator=(const NetworkGraph& other);

    /**
     * Two network graphs are defined to be equal, if they
     * have the same graph nodes and all graph nodes have
     * the same successors.
     *
     * @note The graph node comparison is done by using
     *       the GraphNode class' equality operator.
     *       Therefore, equal graphs do not necessarily
     *       share GraphNode objects.
     */
    bool operator==(const NetworkGraph& other) const;

    /**
     * Convenience operator. Returns !(*this == other).
     */
    bool operator!=(const NetworkGraph& other) const;

    /**
     * Indicates whether the given processor is containing within this graph.
     *
     * @return  true if the processor is contained, or false otherwise of if
     *          the given pointer is NULL.
     */
    bool containsProcessor(Processor* const processor) const;

    /**
     * Returns all processors from all nodes in this graph. The order may be
     * chosen to be unsorted or the processors can be ordered according to the
     * ID of their nodes. The ID of the nodes corresponds to the order of insertion
     * into the graph.
     *
     * @param   sortByID    Determines whether to sort the processor according to
     *                      the IDs of their corresponding nodes or not.
     * @return  Vector containing pointer to all Processor objects which are
     *          represented by nodes within this graph. Does not contain duplicates.
     */
    std::vector<Processor*> getProcessors(const bool sortByID = false) const;

     /**
     * Returns all graphNodes from this graph
     *
     * @return  Vector containing pointer to all GraphNode objects within this graph.
     */

    std::vector<NetworkGraph::GraphNode*> getGraphNodes() const;

    /**
     * Returns a subset of this graph containing only processors with ports and connections
     * between them, which match the given PortTypeCheck criteria. All other processors
     * having different ports and connections can be chosen to be disposed.
     * The caller has to take care of deleting the returned graph.
     *
     * @param   type        PortTypeCheck specifying, which type of port and connections shall
     *                      be maintained in the new sub-graph.
     * @param   loopType    Determines the type of ports and connections to be treated as static loops
     *                      and unrolled. \sa buildGraph
     * @param   keepNodes   Determines whether to keep nodes, which do not match the
     *                      criteria specified via the given PortTypeCheck.
     * @return  A new NetworkGraph which contains only nodes with processors having ports
     *          and connections of the type specified by the parameter type.
     */
    NetworkGraph* getSubGraph(const PortTypeCheck& type = GenericPortTypeCheck<Port>(),
                              const PortTypeCheck& loopType = PortTypeCheckReject(),
                              const bool keepNodes = false) const;

    /**
     * Returns the transposed graph. A transposed directed graph is a graph
     * where all edges (u, v) are replaced by the edges (v, u). This is done by
     * inverting the role of predecessor and successor for all nodes within
     * this graph.
     * The caller has to take care of deleting the returned graph.
     *
     * @return  Pointer to the transposed graph. Caller has to free it via
     *          delete.
     */
    NetworkGraph* getTransposed() const;

    /**
     * Returns true, if \p successor is a direct or indirect successor of \p predecessor
     * in the network graph.
     */
    bool isSuccessor(Processor* predecessor, Processor* successor) const;

    /**
     * Returns true, if \p successor is a direct or indirect successor of \p predecessorPort
     * in the network graph. The passed port may be an inport or an outport.
     */
    bool isSuccessor(Port* predecessorPort, Processor* successor) const;

    /**
     * Returns true, if \p successorPort is a direct or indirect successor of \p predecessor
     * in the network graph. The passed port may be an inport or an outport.
     */
    bool isSuccessor(Processor* predecessor, Port* successorPort) const;

    /**
     * Returns true, if \p successor is a direct or indirect successor of \p predecessor
     * in the network graph. Both parameters may be inports or outports.
     */
    bool isSuccessor(Port* predecessor, Port* successor) const;

    /**
     * Returns true, if \p processor lies on a path between \p pathRoot and \p pathEnd.
.     */
    bool isPathElement(Processor* processor, Processor* pathRoot, Processor* pathEnd) const;

    /**
     * Returns true, if \p processor lies on a path between \p pathRootPort and \p pathEndPort.
     * Both ports may either be an inport or an outport.
     */
    bool isPathElement(Processor* processor, Port* pathRootPort, Port* pathEndPort) const;

    /**
     * Returns true, if \p port lies on a path between \p pathRootPort and \p pathEndPort.
     * All passed ports may either be an inport or an outport.
     * If \p port is identical to \p pathRootPort or \p pathEndPort and a path between root and end exists,
     * \p port is defined to be an path element.
.    */
    bool isPathElement(Port* port, Port* pathRootPort, Port* pathEndPort) const;

    /**
     * Returns all predecessing processors of the passed processors.
     *
     * @param  processors The processors whose predecessors are to be determined.
     * @return Set of predecessing processors without duplicates, including the passed ones.
     */
    std::set<Processor*> getPredecessors(const std::set<Processor*>& processors) const;

    /// @overload
    std::set<Processor*> getPredecessors(const std::vector<Processor*>& processors) const;

    /// @overload
    std::set<Processor*> getPredecessors(Processor* processor) const;

    /**
     * Returns all predecessing processors of the passed port,
     * which must be an inport.
     *
     * @param port The port whose predecessors are to be determined.
     * @return Set of predecessing processors without duplicates.
     */
    std::set<Processor*> getPredecessors(Port* inport) const;

    /**
     * Returns all successing processors of the passed processors.
     *
     * @param  processors The processors whose successors are to be determined.
     * @return Set of successing processors without duplicates, including the passed ones.
     */
    std::set<Processor*> getSuccessors(const std::set<Processor*>& processors) const;

    /// @overload
    std::set<Processor*> getSuccessors(Processor* processor) const;

    /**
     * Returns all successing processors of the passed port,
     * which must be an outport.
     *
     * @param port The port whose successors are to be determined.
     * @return Set of successing processors without duplicates.
     */
    std::set<Processor*> getSuccessors(Port* outport) const;

    /**
     * Traverses the entire graph in breadth-first order. This means,
     * that every node within the graph is a potential root. The roots
     * for the search are determined by <code>identifyRoots()</code>.
     * Thereby, all connected components will be found. This method
     * internally calls <code>traverseBreadthFirstInternal()</code>.
     * After traversing, all bool marks of the nodes are set to true.
     * Before traversing, this method sets all bool markers of the nodes
     * to false. Furthermore, the breadth-first search determines the
     * minimal distance to the root for this node, for each connected
     * component of the graph.
     *
     * @param   visitor A visitor according to the visitor design-pattern
     *                  which is applied to each node when it is discovered
     *                  by breadth-first search, and if it is not NULL.
     */
    void fullTraverseBreadthFirst(GraphVisitor* const visitor = 0) const;

    /**
     * Traverses the entire graph in depth-first order. This means, that
     * every node within the graph is a potential root. The roots for the
     * search are determined by <code>identifyRoots()</code>. Thereby,
     * all connected components will be found. This method internally calls
     * <code>traverseBreadthFirstInternal()</code>.
     * After traversing, all discovery timestamps and all finalization
     * timestamps are > 0. Before traversing, all these timestamps are set
     * to 0 again.
     *
     * @param   onDiscoveryVisitor  A vistor according to the visitor desgin-
     * pattern which will be applied to all nodes when they are discovered by
     * the depth-first search, and if it is not NULL.
     * @param   onFinishVisitor  A vistor according to the visitor desgin-
     * pattern which will be applied to all nodes when they are finished by
     * the depth-first search, and if it is not NULL.
     */
    void fullTraverseDepthFirst(GraphVisitor* const onDiscoveryVisitor = 0,
        GraphVisitor* const onFinishVisitor = 0) const;

    /**
     * Sorts the graph topolocial. The topological sorting is used to determine the
     * order of evaluation for the processors by the NetworkEvaluator.
     *
     * NOTE: Topological sorting works fine in acyclic directed graphs, but I don't know
     * whether this will work in all cases, if cycles are present. No check is for cycles
     * in the graph is performed, because they might be wanted.
     *
     * @param   processorSubset If the given set is not empty, only the processors which are
     * are also in this set will be included in the returned vector. Otherwise all processors
     * will be included.
     * @return  Topological sorting of processors which can be used as rendering order.
     */
    std::vector<Processor*> sortTopologically(const std::set<Processor*>& processorSubset = std::set<Processor*>()) const;

    /**
     * Performs breadth-first search in the graph, starting at the node which
     * contains the given processor, if present. Otherwise, nothing will happen.
     * This method internally calls <code>traverseBreadthFirstInternal()</code>.
     * After traversing, all bool marks of the nodes are set to true.
     * Before traversing, this method sets all bool markers of the nodes
     * to false. Furthermore, the breadth-first search determines the
     * minimal distance to the given root for this node.
     *
     * @param   root    The node containing this processor will act as the root
     *                  for the breadth-first search.
     * @param   visitor A visitor according to the visitor design-pattern
     *                  which is applied to each node when it is discovered
     *                  by breadth-first search, and if it is not NULL.
     * @return  true if the search was perfomed, false if it was not, e.g. if
     *          no node with the given processor could be found.
     */
    bool traverseBreadthFirst(Processor* root, GraphVisitor* const visitor = 0) const;

    /**
     * Performs breadth-first search in the graph, starting at the node which
     * contains the given processor, if present. Otherwise, nothing will happen.
     * This method internally calls <code>traverseBreadthFirstInternal()</code>.
     * After traversing, all discovery timestamps and all finalization
     * timestamps are > 0. Before traversing, all these timestamps are set
     * to 0 again.
     *
     * @param   root    The node containing this processor will act as the root
     *                  for the breadth-first search.
     * @param   onDiscoveryVisitor  A vistor according to the visitor desgin-
     *      pattern which will be applied to all nodes when they are discovered by
     *      the depth-first search, and if it is not NULL.
     * @param   onFinishVisitor  A vistor according to the visitor desgin-
     *      pattern which will be applied to all nodes when they are finished by
     *      the depth-first search, and if it is not NULL.
     * @return  true if the search was perfomed, false if it was not, e.g. if
     *          no node with the given processor could be found.
     */
    bool traverseDepthFirst(Processor* root, GraphVisitor* const onDiscoveryVisitor = 0,
        GraphVisitor* const onFinishVisitor = 0) const;

private:
    typedef std::set<GraphNode*, GraphNode::GraphNodeComparator> NodeSet;

    /**
     * Comparator used to sort the nodes topological. This requires a previous
     * depth-first search for the affected nodes.
     */
    struct NodeSortTopologicalComparator {
        /**
         * This operator returns the result from the comparion of the 'finish' timestamps
         * of the nodes. The nodes have to sorted in descending order of their 'finish'
         * timestamps.
         */
        bool operator ()(const GraphNode* const n1, const GraphNode* const n2) const {
            if ((n1 == 0) || (n2 == 0))
                return false;
            return (n1->getFinished() >= n2->getFinished());
        }
    };

    /**
     * Comparator class to compare pointers to GraphNode objects by their ID only, rather
     * than taking into account the pointers to the internal Processor objects as well.
     */
    struct NodeIDComparator {
        /**
         * Operator which defines n1 to be smaller than n2 if the ID of n1 is less
         * than the one of n2.
         */
        bool operator ()(const GraphNode* const n1, const GraphNode* const n2) const {
            if ((n1 == 0) || (n2 == 0))
                return false;
            return (n1->getID() < n2->getID());
        }
    };

private:
    /**
     * Ctor creating an empty graph.
     */
    NetworkGraph();

    /**
     * Ctor construction a graph with deep copies of the nodes in the given set, but
     * without making any connection between the nodes and without making any
     * determination of successors for the nodes.
     * This is used for instance by creating transposed graphs.
     */
    NetworkGraph(const NodeSet& nodes);

    /**
     * Adds a new GraphNode fot the given processor to this graph and returns it.
     * The new GraphNode objects gets an ID which corresponds to its number on
     * insertion. If the node already exists in the graph, the existing node will
     * be returned and no new node will be created.
     *
     * @param   processor   Pointer to the Processor object for which a node shall be
     *                      added. This pointer must not be NULL.
     * @param   added   Pointer to a bool which may hold the information, if a new node
     *                  has been created or not. If this pointer is NULL, the information
     *                  about possible new creations will not be stored to it.
     * @return  The new GraphNode object, if it was created newly, or the exisiting one,
     *          containing the given processor.
     */
    GraphNode* addNode(Processor* const processor, bool* const added = 0);

    /**
     * Connects the nodes within the graph. The given PortTypeCheck can confine
     * the connections to be made.
     *
     * @param   ptc PortTypeCheck which determines, whether a connection between
     *              the nodes will be established, according to the type of the
     *              port on the internal Processor pointer of the node.
     */
    void connectNodes(const PortTypeCheck& ptc = GenericPortTypeCheck<Port>());

    /**
     * Unrolls all loops within the graph by duplicating all graph nodes between
     * the source and the destination of the loop connection and inserting them after
     * the loop connection destination.
     *
     * @param loopType Determines the type of ports and connections that are to be
     *                 considered loop edges and additionally retrieves the number
     *                 of iterations of a loop connection.
     */
    void unrollLoops(const PortTypeCheck& loopType = GenericPortTypeCheck<Port>());

    /**
     * Acutally creates the graph from the given processors. The PortTypeCheck may be
     * used to determine what kind of processor and connections may be contained in
     * the graph. The LoopPortTypeCheck determines what connections are to be considered
     * as static loop and unrolled.
     * Before creating the graph, the previous graph will be cleared by a call to
     * <code>clearGraph()</code>.
     *
     * @param   processors  Vector of all processors which shall be contained.
     * @param   ptc         PortTypeCheck to confine the nodes/processors and connections
     *                      to those, which fulfill its criterion.
     * @param   loopPtc     Ports and connections of this type are considered as static loops
     *                      and are unrolled, i.e., all nodes in the graph representing
     *                      processors between the source and the destination of the loop edge
     *                      are duplicated and inserted after the loop destination. The number
     *                      of duplications of the loop body is also determined by loopPtc.
     * @param   keepNodes   Determines whether nodes for processors not matching the
     *                      the criterion by parameter ptc shall be kept in the graph
     *                      or not.
     */
    void buildGraph(const std::vector<Processor*>& processors, const PortTypeCheck& ptc,
                    const PortTypeCheck& loopPtc, const bool keepNodes = true);

    /**
     * Clears the entire graph be deleting all nodes and resetting the node counter.
     */
    void clearGraph();

    /**
     * Copies the graph. This is done by making deep copies of the nodes from the other
     * graph, but without successors. When all nodes are copied, all connections for
     * the successor in the other graph will be made for the copies. This is necessary
     * because if the successors would have been copied from the other graph, the nodes
     * in this new graph would know their memory addresses in the other graph! But a
     * corresponding successor S for a node N in this graph has a diffrent address S' in
     * the other graph.
     */
    void copy(const NetworkGraph& other);

    /**
     * Removes all nodes from the graph that are not contained
     * by \p nodes.
     */
    void pruneGraph(const std::set<GraphNode*>& nodes);

    /**
     * Returns a set of all GraphNodes wrapping the passed processor.
     */
    std::set<GraphNode*> findNodes(Processor* processor) const;

    /**
     * Returns the pointer to the GraphNode object with the given ID, if
     * existing, or NULL otherwise.
     */
    GraphNode* findNode(const int id) const;

    /**
     * Returns all successing graph nodes of the passed nodes.
     *
     * @param  nodes The graph nodes whose successors are to be determined.
     * @return Set of successing graph nodes, including the passed ones.
     *         The result set does not contain duplicate graph nodes, but a
     *         single processor object might be wrapped by multiple graph nodes.
     */
    std::set<GraphNode*> getSuccessorNodes(const std::set<GraphNode*>& nodes) const;

    /**
     * Returns all predecessing graph nodes of the passed nodes.
     *
     * @param  nodes The graph nodes whose predecessors are to be determined.
     * @return Set of predecessing graph nodes, including the passed ones.
     *         The result set does not contain duplicate graph nodes, but a
     *         single processor object might be wrapped by multiple graph nodes.
     */
    std::set<GraphNode*> getPredecessorNodes(const std::set<GraphNode*>& nodes) const;

    /**
     * Returns all nodes that lie on a path between \p predecessor and \p successor.
     */
    std::set<GraphNode*> getPathNodes(const GraphNode* predecessor, const GraphNode* successor) const;

    /**
     * Identifies roots for full travering the graph via breadth-first or depth-first
     * search. A node is assumed to be a root, if it has a in-degree of 0.
     */
    std::vector<GraphNode*> identifyRoots() const;

    /**
     * Sets the 'discovery' timestamps of all nodes to the given value. This is used to
     * intialize depth-first search.
     */
    void setAllNodesDiscovered(const int value) const;

    /**
     * Sets the 'finish' timestamps of all nodes to the given value. This is used to
     * intialize depth-first search.
     */
    void setAllNodesFinished(const int value) const;

    /**
     * Sets the bool mark of all nodes to the given value. This is used to
     * intialize breadth-first search.
     */
    void setAllNodesMarked(const bool marked) const;

    /**
     * Performs the actual breadth-first rraversin in this graph by
     * starting at the given GraphNode.
     * After traversing, all bool marks of the nodes which were found
     * during the search are set to true.
     * Before traversing, this method sets all bool markers of the nodes
     * to false. Furthermore, the breadth-first search determines the
     * minimal distance to the given root for this node.
     *
     * @param   root    The node which will act as the root for the
     *                  breadth-first search.
     * @param   visitor A visitor according to the visitor design-pattern
     *                  which is applied to each node when it is discovered
     *                  by breadth-first search, and if it is not NULL.
     */
    void traverseBreadthFirstInternal(GraphNode* const root, GraphVisitor* const visitor = 0) const;

    /**
     * After the transversing, every node which coudl be reaced from the given root
     * has a timestamp for its discovery and a timestamp for its finalization,
     * > 0 each and strictly monotic increasing.
     * The timestamps indicate the order of traversal. A node is finished, if the
     * depth-first search has visited all nodes which are reachable form that node
     * and it is about to 'backtrack' to its predecessor.
     *
     * @param   root    The node which will act as the root for the depth-first search.
     * @param   startTime   Timestamp at which the depth-first search starts. This
     *      should be 0 for the first connected component, and T1 for the second,
     *      with T1 denoting the timestamp returned by this method from the for run,
     *      etc..
     * @param   onDiscoveryVisitor  A vistor according to the visitor desgin-
     *      pattern which will be applied to all nodes when they are discovered by
     *      the depth-first search, and if it is not NULL.
     * @param   onFinishVisitor  A vistor according to the visitor desgin-
     *      pattern which will be applied to all nodes when they are finished by
     *      the depth-first search, and if it is not NULL.
     * @return  Return the timestamp of the finalization of the root. This value
     *      is shall be used for a successive call to this method as startTime when
     *      traversing another connected component of this graph with a different root.
     */
    int traverseDepthFirstInternal(GraphNode* const root, const int startTime = 0,
        GraphVisitor* const onDiscoveryVisitor = 0, GraphVisitor* const onFinishVisitor = 0) const;

private:

    /** The set containing all nodes of this graph. */
    NodeSet nodes_;

    /**
     * Counter for the nodes. This counter is used to assign IDs to new nodes
     * according to their order of insertion.
     */
    size_t nodeCounter_;

    /**
     * Transposed graph. Is used for certain operations like
     * the detection of predecessors. Is generated on-demand.
     */
    mutable NetworkGraph* networkGraphTransposed_;

    static const std::string loggerCat_;
};

}   // namespace

#endif
