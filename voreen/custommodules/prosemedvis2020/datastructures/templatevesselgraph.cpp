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

#include "templatevesselgraph.h"
#include <vector>

namespace voreen {

    std::vector<TemplateVesselGraphNode> TemplateVesselGraph::getNodes() const {
        return nodes;
    }
    
    std::vector<TemplateVesselGraphEdge> TemplateVesselGraph::getEdges() const {
        return edges;
    }

    void TemplateVesselGraph::addNode(TemplateVesselGraphNode& node) {
        for (TemplateVesselGraphNode n : nodes) {
            if (n.getID() == node.getID()) {
                return;
            }
        }
        nodes.push_back(node);
    }
    
    void TemplateVesselGraph::addEdge(TemplateVesselGraphEdge& edge) {
        for (TemplateVesselGraphEdge e : edges) {
            if (e.getLabel() == edge.getLabel()) {
                return;
            }
        }
        edges.push_back(edge);
    }

	const TemplateVesselGraphNode& TemplateVesselGraph::getNode(size_t id) const{
		for (int i = 0; i < nodes.size(); i++) {
			if (nodes[i].getID() == id) {
				return nodes[i];
			}
		}
	}

	const TemplateVesselGraphEdge& TemplateVesselGraph::getEdge(size_t id) const{
		for (int i = 0; i < edges.size(); i++) {
			if (edges[i].getID() == id) {
				return edges[i];
			}
		}
	}


    TemplateVesselGraphNode::TemplateVesselGraphNode(size_t id, std::string label, NodeType type) : id(id), label(label), type(type) {}

    void TemplateVesselGraphNode::addEdge(size_t newEdge) {
        edgeIDs.push_back(newEdge);
    }

    bool TemplateVesselGraphNode::hasEdge(size_t edge) {
        for (size_t e : edgeIDs) {
            if (e == edge) {
                return true;
            }
        }
        return false;
    }

    int TemplateVesselGraphNode::getDegree() const {
        return edgeIDs.size();
    }

    size_t TemplateVesselGraphNode::getID() const {
        return id;
    }

    NodeType TemplateVesselGraphNode::getType() const {
        return type;
    }
    std::vector<size_t> TemplateVesselGraphNode::getEdgeIDs() const {
        return edgeIDs;
    };

	void TemplateVesselGraphNode::setID(size_t newID) {
		id = newID;
	}

	std::string TemplateVesselGraphNode::getLabel() const{
		return label;
	}

    TemplateVesselGraphEdge::TemplateVesselGraphEdge(size_t id, std::string label, float relativeRadius,
        size_t startNode, size_t endNode) 
        : id(id), label(label), relativeRadius(relativeRadius), startNode(startNode), endNode(endNode) {}

    size_t TemplateVesselGraphEdge::getStart(){
        return startNode;
    }

    size_t TemplateVesselGraphEdge::getEnd(){
        return endNode;
    }

    std::string TemplateVesselGraphEdge::getLabel() const {
        return label;
    }
    float TemplateVesselGraphEdge::getRadius() const{
        return relativeRadius;
    }

	void TemplateVesselGraphEdge::setID(size_t newID) {
		id = newID;
	}

	size_t TemplateVesselGraphEdge::getID() const{
		return id;
	}
}
