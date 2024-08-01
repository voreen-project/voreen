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

#ifndef VRN_CONCRETEVESSELGRAPH_H
#define VRN_CONCRETEVESSELGRAPH_H

#include <iostream>
#include <vector>
#include "tgt/vector.h"
#include "templatevesselgraph.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/io/serialization/serializer.h"
#include "voreen/core/io/serialization/deserializer.h"


namespace voreen {


	struct ConcreteVesselGraph;
	struct ConcreteVesselGraphNode;
	struct ConcreteVesselGraphEdge;
	struct ConcreteVesselGraphEdgePoint;

	struct ConcreteVesselGraph :public Serializable {

		std::vector<ConcreteVesselGraphNode> nodes;
		std::vector<ConcreteVesselGraphEdge> edges;


	public:
		ConcreteVesselGraphEdge& getEdge(size_t id);
		std::vector<ConcreteVesselGraphNode> getNodes() const;
		std::vector<ConcreteVesselGraphEdge> getEdges() const;
		void addNode(ConcreteVesselGraphNode node);
		void addEdge(ConcreteVesselGraphEdge edge);
		ConcreteVesselGraphNode& getNode(size_t id);
		void removeEdge(size_t edgeID);
		void changeNodeLabel(size_t id, std::string newLabel);
		void changeEdgeLabel(size_t id, std::string newLabel);

		#ifndef VRN_VESSELNETWORKANALYSIS_MINIMAL_VESSELGRAPH
		virtual void serialize(Serializer& s) const;
		virtual void deserialize(Deserializer& s);
		#endif
	};

	struct ConcreteVesselGraphNode : public Serializable {

		size_t id;
		std::string label;
		tgt::vec3 pos;
		NodeType type;
		std::vector<size_t> edgeIDs;

	public:
		ConcreteVesselGraphNode(size_t id, std::string label, NodeType type, tgt::vec3 position);
		ConcreteVesselGraphNode() = default;

		size_t getId();
		void setId(size_t id);
		std::string getLabel();
		void setLabel(std::string label);
		tgt::vec3 getPos() const;
		NodeType getType() const;

		int getDegree() const;
		void addEdge(size_t newEdge);
		bool hasEdge(size_t edge);

		std::vector<size_t> getEdgeIDs();
		void eraseEdge(size_t id);

		#ifndef VRN_VESSELNETWORKANALYSIS_MINIMAL_VESSELGRAPH
		virtual void serialize(Serializer& s) const;
		virtual void deserialize(Deserializer& s);
		#endif

	};


	struct ConcreteVesselGraphEdge : public Serializable {

		size_t id;
		std::string label;
		float radius;

		size_t startNode;
		size_t endNode;
		std::vector<ConcreteVesselGraphEdgePoint> edgePoints;

	public:
		ConcreteVesselGraphEdge(size_t id, std::string label, float radius,
			size_t startNode, size_t endNode);
		ConcreteVesselGraphEdge() = default;

		size_t getID();
		void setID(size_t newID);
		std::string getLabel() const;
		void setLabel(std::string label);
		float getRadius() const;
		void setRadius(float radius);
		size_t getStart();
		size_t getEnd();
		std::vector<ConcreteVesselGraphEdgePoint> getEdgePoints();
		void setEdgePoints(std::vector<ConcreteVesselGraphEdgePoint> newEdgePoints);
		std::tuple<tgt::vec3, tgt::vec3,float> interpolate(float t);

#ifndef VRN_VESSELNETWORKANALYSIS_MINIMAL_VESSELGRAPH
		virtual void serialize(Serializer& s) const;
		virtual void deserialize(Deserializer& s);
#endif
	};



	struct ConcreteVesselGraphEdgePoint :public Serializable {
		tgt::vec3 pos;
		float radius;

	public:
		tgt::vec3 getPos() const;
		float getRadius() const;
		ConcreteVesselGraphEdgePoint(tgt::vec3 pos, float radius);
		ConcreteVesselGraphEdgePoint() = default;

		#ifndef VRN_VESSELNETWORKANALYSIS_MINIMAL_VESSELGRAPH
		virtual void serialize(Serializer& s) const;
		virtual void deserialize(Deserializer& s);
		#endif
	};
}
#endif
