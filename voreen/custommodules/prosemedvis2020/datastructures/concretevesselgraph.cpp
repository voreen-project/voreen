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

#include "concretevesselgraph.h"
#include <vector>
#include "tgt/vector.h"

namespace voreen {

    // Implementations of ConcreteVesselGraph

    std::vector<ConcreteVesselGraphNode> ConcreteVesselGraph::getNodes() const {
        return nodes;
    }

    std::vector<ConcreteVesselGraphEdge> ConcreteVesselGraph::getEdges() const {
        return edges;
    }

	ConcreteVesselGraphEdge& ConcreteVesselGraph::getEdge(size_t id) {
		for (ConcreteVesselGraphEdge& edge : edges) {
			if (edge.getID() == id) {
				return edge;
			}
		}
	}

    void ConcreteVesselGraph::addNode(ConcreteVesselGraphNode node) {
        for (ConcreteVesselGraphNode& n : nodes) {
            if (n.getId() == node.getId()) {
                return;
            }
        }
        nodes.emplace_back(node);
    }

    void ConcreteVesselGraph::addEdge(ConcreteVesselGraphEdge edge) {
        for (ConcreteVesselGraphEdge& e : edges) {
            if (e.getID() == edge.getID()) {
                return;
            }
        }
        edges.push_back(edge);
    }

    ConcreteVesselGraphNode& ConcreteVesselGraph::getNode(size_t id) {
		for (int i = 0; i < nodes.size(); i++) {
            if (nodes[i].getId() == id) {
                return nodes[i];
            }
        }
    }
    void ConcreteVesselGraph::removeEdge(size_t edgeID) {
        for (int i = 0;i < edges.size();i++) {
            if (edges[i].getID() == edgeID) {
                edges.erase(edges.begin()+i);
                break;
            }
        }
    }
    void ConcreteVesselGraph::changeNodeLabel(size_t id, std::string newLabel) {
		for (int i = 0; i < nodes.size(); i++) {
            if (nodes[i].getId() == id) {
                nodes[i].label = newLabel;
				return;
            }
        }
    }
    void ConcreteVesselGraph::changeEdgeLabel(size_t id, std::string newLabel) {
		for (int i = 0; i < edges.size(); i++) {
            if (edges[i].getID() == id) {
                edges[i].setLabel(newLabel);
				return;
            }
        }
    }

	void ConcreteVesselGraph::serialize(Serializer& s) const {

		s.serialize("nodes", nodes);

		s.serialize("edges", edges);

	}
	void ConcreteVesselGraph::deserialize(Deserializer& s) {

		s.deserialize("nodes", nodes);

		s.deserialize("edges", edges);
	}

    // Implementations of ConcreteVesselGraphNode

    ConcreteVesselGraphNode::ConcreteVesselGraphNode(size_t id, std::string label, NodeType type,
        tgt::vec3 position) : id(id), label(label), type(type), pos(position) {
	}

	size_t ConcreteVesselGraphNode::getId() {
        return id;
    }
    void ConcreteVesselGraphNode::setId(size_t id) {
        this->id=id;
    }
    tgt::vec3 ConcreteVesselGraphNode::getPos() const {
        return pos;
    }

    NodeType ConcreteVesselGraphNode::getType() const {
        return type;
    }

    int ConcreteVesselGraphNode::getDegree() const {
        return edgeIDs.size();
    }

    void ConcreteVesselGraphNode::addEdge(size_t newEdge) {
        edgeIDs.emplace_back(newEdge);
    }

    bool ConcreteVesselGraphNode::hasEdge(size_t edge) {
        for (size_t e : edgeIDs) {
            if (e == edge) {
                return true;
            }
        }
        return false;
    }
    std::vector<size_t> ConcreteVesselGraphNode::getEdgeIDs() {
        return edgeIDs;
    }

	std::string ConcreteVesselGraphNode::getLabel() {
		return label;
	}

	void ConcreteVesselGraphNode::setLabel(std::string newLabel) {
		label = newLabel;
	}

	void ConcreteVesselGraphNode::eraseEdge(size_t id) {
		auto iter = std::find(edgeIDs.begin(), edgeIDs.end(), id);
		if (iter != edgeIDs.end()) {
			edgeIDs.erase(iter);
		}
	}

	void ConcreteVesselGraphNode::serialize(Serializer& s) const {
		s.serialize("id", id);
		s.serialize("label", label);
		int iType = (int)type;
		s.serialize("type", iType);
		s.serialize("pos", pos);
		s.serialize("edges", edgeIDs);
	}
	
	void ConcreteVesselGraphNode::deserialize(Deserializer& s) {
		s.deserialize("id", id);
		s.deserialize("label", label);
		int iType;
		s.deserialize("type", iType);
		type = (NodeType) iType;
		s.deserialize("pos", pos);
		s.deserialize("edges", edgeIDs);

	}

    // Implementations of ConcreteVesselGraphEdge
    ConcreteVesselGraphEdge::ConcreteVesselGraphEdge(size_t id, std::string label, float radius,
		size_t startNode, size_t endNode)
        : id(id), label(label), radius(radius), startNode(startNode), endNode(endNode) {}

	size_t ConcreteVesselGraphEdge::getID() {
		return id;
	}

	void ConcreteVesselGraphEdge::setID(size_t newID) {
		id = newID;
	}

    std::string ConcreteVesselGraphEdge::getLabel() const {
        return label;
    }
    void ConcreteVesselGraphEdge::setLabel(std::string label) {
        this->label = label;
    }
    float ConcreteVesselGraphEdge::getRadius() const {
        return radius;
    }
    void ConcreteVesselGraphEdge::setRadius(float radius) {
        this->radius = radius;
    }

	size_t ConcreteVesselGraphEdge::getStart()  {
        return startNode;
    }

	size_t ConcreteVesselGraphEdge::getEnd()   {
        return endNode;
    }

    void ConcreteVesselGraphEdge::setEdgePoints(std::vector<ConcreteVesselGraphEdgePoint> newEdgePoints) {
        edgePoints = newEdgePoints;
    }


    std::tuple<tgt::vec3, tgt::vec3,float> ConcreteVesselGraphEdge::interpolate(float t) {
        float position = (edgePoints.size() - 1) * t;
        int edgepoint = (int)position;
        float  k= position-edgepoint;

        float radius;
        tgt::vec3 normal;
        tgt::vec3 pos;
        if (edgepoint != edgePoints.size() - 1) {
            radius = edgePoints[edgepoint].getRadius()*(1-k) + edgePoints[edgepoint+1].getRadius()*k;
            tgt::vec3 firstPoint = edgePoints[edgepoint].getPos();
            tgt::vec3 secondPoint = edgePoints[edgepoint+1].getPos();
            normal = secondPoint - firstPoint;
            pos = firstPoint+ k*normal;
            tgt::normalize(normal);
        }
        else {
            radius =  edgePoints[edgepoint].getRadius();
            tgt::vec3 firstPoint = edgePoints[edgepoint-1].getPos();
            tgt::vec3 secondPoint = edgePoints[edgepoint].getPos();
            normal = secondPoint - firstPoint;
            pos = edgePoints[edgepoint].getPos();
            tgt::normalize(normal);
        }
        return std::make_tuple(pos,normal,radius);
    }

	std::vector<ConcreteVesselGraphEdgePoint> ConcreteVesselGraphEdge::getEdgePoints() {
		return edgePoints;
	}

	void ConcreteVesselGraphEdge::serialize(Serializer& s) const {
		s.serialize("id", id);
		s.serialize("startnode", startNode);
		s.serialize("endnode", endNode);
		s.serialize("edgePoints", edgePoints);
		s.serialize("label", label);
		s.serialize("radius", radius);
	}

	void ConcreteVesselGraphEdge::deserialize(Deserializer& s) {

		s.deserialize("id", id);
		s.deserialize("startnode", startNode);
		s.deserialize("endnode", endNode);
		s.deserialize("edgePoints", edgePoints);
		s.deserialize("label", label);
		s.deserialize("radius", radius);
	}



    // Implementations of ConcreteVesselGraphEdgePoint

    tgt::vec3 ConcreteVesselGraphEdgePoint::getPos() const {
        return pos;
    }

    float ConcreteVesselGraphEdgePoint::getRadius() const {
        return radius;
    }
    ConcreteVesselGraphEdgePoint::ConcreteVesselGraphEdgePoint(tgt::vec3 pos, float radius) :pos(pos), radius(radius) {};

	void ConcreteVesselGraphEdgePoint::serialize(Serializer& s) const {
		s.serialize("pos", pos);
		s.serialize("radius", radius);
	}
	
	void ConcreteVesselGraphEdgePoint::deserialize(Deserializer& s) {
		s.deserialize("pos", pos);
		s.deserialize("radius", radius);
	}
}