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

#ifndef VRN_TEMPLATEVESSELGRAPH_H
#define VRN_TEMPLATEVESSELGRAPH_H

#include <iostream>
#include <vector>
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/io/serialization/serializer.h"
#include "voreen/core/io/serialization/deserializer.h"


namespace voreen {

    struct TemplateVesselGraph;
    struct TemplateVesselGraphNode;
    struct TemplateVesselGraphEdge;

    enum NodeType { END_NODE=1, BOUNDARY_NODE, BIFURCATION_NODE };


    struct TemplateVesselGraph {

        private:
            std::vector<TemplateVesselGraphNode> nodes;
            std::vector<TemplateVesselGraphEdge> edges;

        public:
            std::vector<TemplateVesselGraphNode> getNodes() const;
            std::vector<TemplateVesselGraphEdge> getEdges() const;
            void addNode(TemplateVesselGraphNode& node);		
            void addEdge(TemplateVesselGraphEdge& edge);
			const TemplateVesselGraphNode& getNode(size_t id) const;
			const TemplateVesselGraphEdge& getEdge(size_t id) const;
    };

    struct TemplateVesselGraphNode {

        size_t id;
		std::string label;
        NodeType type;

        std::vector<size_t> edgeIDs;

        public:
            TemplateVesselGraphNode(size_t id, std::string label, NodeType type);
    
            void addEdge(size_t newEdge);
            bool hasEdge(size_t edge);
            std::vector<size_t> getEdgeIDs() const;

            int getDegree() const;
			void setID(size_t newID);
            size_t getID() const;
			std::string getLabel() const;
            NodeType getType() const;
    };

    struct TemplateVesselGraphEdge {

		size_t id;
        std::string label;
        float relativeRadius;

        size_t startNode;
        size_t endNode;

        public:
            TemplateVesselGraphEdge(size_t id, std::string label, float relativeRadius, size_t startNode, size_t endNode);

            size_t getStart();
            size_t getEnd();
            std::string getLabel() const;
            float getRadius() const;
			size_t getID() const;
			void setID(size_t newID);
    };



}

#endif
