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

#include "concretevesselgraphcreator.h"
#include "modules/vesselnetworkanalysis/datastructures/vesselgraph.h"
#include "voreen/core/datastructures/diskarraystorage.h"
#include "custommodules/flowsimulation/utils/utils.h"
namespace voreen {
    
    ConcreteVesselGraphCreator::ConcreteVesselGraphCreator() : Processor(),
        flowInport(Port::INPORT, "flow", "Flow Input"),
        templateVesselGraphInport(Port::INPORT, "TemplateVesselGraph", "Template Vessel Graph Input"),
        vesselGraphInport(Port::INPORT, "VesselGraph", "Vessel Graph Input"),
        concreteVesselGraphOutport(Port::OUTPORT, "ConcreteVesselGraph", "Concrete Vessel Graph Output"),
		filePath("filePath", "Path for saved file", "Choose Path", "", ".json", FileDialogProperty::SAVE_FILE, Processor::VALID), 
		saveButton("saveButton", "Save", Processor::VALID)
    {
        addPort(flowInport);
        addPort(templateVesselGraphInport);
        addPort(vesselGraphInport);
        addPort(concreteVesselGraphOutport);

		addProperty(filePath);
			ON_CHANGE(filePath, ConcreteVesselGraphCreator, saveCurrentConcreteVesselGraph);
		addProperty(saveButton);
			ON_CHANGE(saveButton, ConcreteVesselGraphCreator, saveCurrentConcreteVesselGraph);
			
    }

    Processor* ConcreteVesselGraphCreator::create() const {
        return new ConcreteVesselGraphCreator();
    }

	void ConcreteVesselGraphCreator::saveCurrentConcreteVesselGraph() {
		const std::string path = filePath.get();
		if (path.empty()) {
			return;
		}

		const ConcreteVesselGraph* input = concreteVesselGraphOutport.getData();
		if (!input) {
			return;
		}

		JsonSerializer serializer;
		try {
			serializer.serialize("graph", *input);
		}
		catch (SerializationException& s) {
			LERROR("Could not serialize graph: " << s.what());
			return;
		}
		try {
			std::fstream f(path, std::ios::out);

			serializer.write(f, true, true);
		}
		catch (...) {
			LERROR("Could not save graph " << path);
			return;
		}
		LINFO("Saved graph to " << path << ".");
	}


    void ConcreteVesselGraphCreator::process() {
		size_t nodeCounter = 0;
		size_t edgeCounter = 1;
        ConcreteVesselGraph* cvg = new ConcreteVesselGraph();

        const VesselGraph* vg = vesselGraphInport.getData();

        const VolumeBase* flow = flowInport.getData();
        VolumeRAMRepresentationLock volumeData(flow);

        const TemplateVesselGraph* tvg = templateVesselGraphInport.getData();

        //add nodes from vesselgraph to concretevesselgraph
        DiskArray<VesselGraphNode> nodeArray = vg->getNodes();
		std::map<uint32_t, size_t> nodeMaps;
        for (int i = 0;i < nodeArray.size();i++) {
            int degree = nodeArray.at(i).getDegree();
            NodeType type;

            if (degree == 1) {
                type = END_NODE;
            }
            else {
                type = BIFURCATION_NODE;
            }
            //add node
            ConcreteVesselGraphNode cn = ConcreteVesselGraphNode(nodeCounter, std::string("not matched"), type, nodeArray.at(i).pos_);
			nodeMaps[nodeArray.at(i).getID().raw()] = nodeCounter;
			nodeCounter++;
            cvg->addNode(cn);
        }
        

       //add edges from vesselgraph to concretevesselgraph
        DiskArray<VesselGraphEdge> edgeArray = vg->getEdges();
        for (int i = 0;i < edgeArray.size();i++) {

            //nodes from current edge
            VesselGraphNode& node1 = edgeArray.at(i).getNode1();
            VesselGraphNode& node2 = edgeArray.at(i).getNode2();

            //connection vector between nodes
            tgt::vec3 diff = node2.pos_ - node1.pos_;

            //flow vector
			std::vector<tgt::vec3> samples = utils::sampleDisk(flowInport.getData(), node1.pos_, tgt::normalize(diff), edgeArray.at(i).getMinRadiusAvg(), false);

			tgt::vec3 accumDirection = tgt::vec3::zero;
			for (const auto& sample : samples) {
				accumDirection += sample;
			}

            tgt::vec3 voxel = tgt::vec3::zero;
            for (size_t channel = 0; channel < 3; channel++) {
                voxel[channel] = volumeData->getVoxelNormalized(node1.pos_, channel);
            }

            float dotproduct = tgt::dot(diff, accumDirection);
			
            //new concretevesselgraphnodes
            ConcreteVesselGraphNode& c1 = cvg->getNode(nodeMaps[node1.getID().raw()]);
            ConcreteVesselGraphNode& c2 = cvg->getNode(nodeMaps[node2.getID().raw()]);

            const DiskArray<VesselSkeletonVoxel>& voxelArray = edgeArray.at(i).getVoxels();

            //dotproduct>=0: connection vector and flow vector alligned => direction of edge c1->c2
            if (dotproduct >= 0) {

                //new concretevesselgraphedge, actual radius added later
                ConcreteVesselGraphEdge ce = ConcreteVesselGraphEdge(edgeCounter, std::string("not matched"), -1, c1.getId(), c2.getId());
				edgeCounter++;

                //add edge to node and concretevesselgraph
                c1.addEdge(ce.getID());

                //convert vesselskeletonvoxel to concretevesselgraphedgepoints
                std::vector<ConcreteVesselGraphEdgePoint> points = std::vector<ConcreteVesselGraphEdgePoint>();

                //add startnode as concretevesselgraphedgepoint
                if (vg->getNode(VGNodeID(c1.getId())).getRadius() < 0.1) {
                    points.push_back(ConcreteVesselGraphEdgePoint(c1.getPos(), voxelArray.at(0).avgDistToSurface_));
                }
                else {
                    points.push_back(ConcreteVesselGraphEdgePoint(c1.getPos(), vg->getNode(VGNodeID(c1.getId())).getRadius()));
                }

                //determine average radius of edge
                float radius = 0;
                for (int i = 0;i < voxelArray.size();i++) {
                    points.push_back(ConcreteVesselGraphEdgePoint(voxelArray.at(i).pos_, voxelArray.at(i).avgDistToSurface_));

                    //update average radius
                    radius += voxelArray.at(i).avgDistToSurface_;
                }
				radius = radius/voxelArray.size();

                //add endnode as concretevesselgraphedgepoint
                if (vg->getNode(VGNodeID(c2.getId())).getRadius() < 0.1) {
                    points.push_back(ConcreteVesselGraphEdgePoint(c2.getPos(), voxelArray.at(voxelArray.size() - 1).avgDistToSurface_));
                }
                else {
                    points.push_back(ConcreteVesselGraphEdgePoint(c2.getPos(), vg->getNode(VGNodeID(c2.getId())).getRadius()));
                }
                ce.setRadius(radius);
                ce.setEdgePoints(points);

				cvg->addEdge(ce);
            }
            else {
                //new concretevesselgraphedge, actual radius added later
                ConcreteVesselGraphEdge ce = ConcreteVesselGraphEdge(edgeCounter, std::string("not matched"), -1, c2.getId(), c1.getId());
				edgeCounter++;
                //add edge to node and concretevesselgraph
                c2.addEdge(ce.getID());

                //convert vesselskeletonvoxel to concretevesselgraphedgepoints
                std::vector<ConcreteVesselGraphEdgePoint> points = std::vector<ConcreteVesselGraphEdgePoint>();

                //add startnode as concretevesselgraphedgepoint
                if (vg->getNode(VGNodeID(c2.getId())).getRadius() < 0.1) {
                    points.push_back(ConcreteVesselGraphEdgePoint(c2.getPos(), voxelArray.at(voxelArray.size() - 1).avgDistToSurface_));
                }
                else {
                    points.push_back(ConcreteVesselGraphEdgePoint(c2.getPos(), vg->getNode(VGNodeID(c2.getId())).getRadius()));
                }

                //determine average radius of edge
                float radius = 0;//average Radius
                for (int i = voxelArray.size() -1;i>-1;i--) {
                    points.push_back(ConcreteVesselGraphEdgePoint(voxelArray.at(i).pos_, voxelArray.at(i).avgDistToSurface_));

                    //update average radius
                    radius +=  voxelArray.at(i).avgDistToSurface_;
                }
				radius = radius / voxelArray.size();

                //add endnode as concretevesselgraphedgepoint
                if (vg->getNode(VGNodeID(c1.getId())).getRadius() < 0.1) {
                    points.push_back(ConcreteVesselGraphEdgePoint(c1.getPos(), voxelArray.at(0).avgDistToSurface_));
                }
                else {
                    points.push_back(ConcreteVesselGraphEdgePoint(c1.getPos(), vg->getNode(VGNodeID(c1.getId())).getRadius()));
                }
                ce.setRadius(radius);
                ce.setEdgePoints(points);

				cvg->addEdge(ce);
            }



        }
		/*
		tgtAssert(cvg->getNodes().size() == 4, "nodes");
		tgtAssert(cvg->getEdges().size() == 3, "edges");
		for (ConcreteVesselGraphNode& node : cvg->getNodes()) {
			tgtAssert(node.getLabel().compare("not matched") == 0, "node label");
			std::cout << node.getDegree() << std::endl;
		}
		for (ConcreteVesselGraphEdge& edge : cvg->getEdges()) {
			tgtAssert(edge.getLabel().compare("not matched") == 0, "edge label");
		}
		*/
        //Insert Boundarynodes
        std::vector<ConcreteVesselGraphNode>nodes = cvg->getNodes();
        std::vector<ConcreteVesselGraphNode>bifurcationnodes = std::vector<ConcreteVesselGraphNode>();
        for (ConcreteVesselGraphNode& node : nodes) {
            if (node.getType() == BIFURCATION_NODE) {
                bifurcationnodes.push_back(node);
            }
        }

        //edges that need a boundary node
        std::vector<size_t>boundaryedges = std::vector<size_t>();
        for (ConcreteVesselGraphNode& bn : bifurcationnodes) {
            //boundaryedges.insert(boundaryedges.end(), bn.getEdges().begin(), bn.getEdges().end());
			for (size_t be : bn.getEdgeIDs()) {
				boundaryedges.push_back(be);
			}
        }
        
        //insert boundarynodes where radiusratio significantly changes
        for (size_t be : boundaryedges) {

            //get edgepoints from current edge
            std::vector<ConcreteVesselGraphEdgePoint>edgepoints = cvg->getEdge(be).getEdgePoints();

			
            float prevRadius=edgepoints[0].getRadius();
            ConcreteVesselGraphEdgePoint boundarynodeedgepoint;
            int i = 0;
			int splitIndex = 0;
			float currentMinRatio = std::numeric_limits<float>::max();
            for (ConcreteVesselGraphEdgePoint& ep : edgepoints) {
				// if radius = 0 edge has no more interesting points
				if (ep.getRadius() == 0) continue;
				float ratio = ep.getRadius()/prevRadius;
                if (ratio < currentMinRatio) { //check radiusratio
					currentMinRatio = ratio;
                    boundarynodeedgepoint = ep;
					splitIndex = i;
                }
                prevRadius = ep.getRadius();
                i++;
            }

            //split edgepoints at boundarynode
            std::vector<ConcreteVesselGraphEdgePoint> edgepoints1(edgepoints.begin(), edgepoints.begin() +splitIndex+1);
            std::vector<ConcreteVesselGraphEdgePoint> edgepoints2(edgepoints.begin()+splitIndex, edgepoints.end());

            //new node and egdes 
            ConcreteVesselGraphNode boundarynode = ConcreteVesselGraphNode(nodeCounter, "not matched",BOUNDARY_NODE,boundarynodeedgepoint.getPos());
			nodeCounter++;

			//calculate new radius
			float radius1 = 0.f;
			for (int i = 1; i < edgepoints1.size(); i++) {
				radius1 += edgepoints1[i].getRadius();
			}
			radius1 /= edgepoints1.size() - 1;
            ConcreteVesselGraphEdge edge1 = ConcreteVesselGraphEdge(edgeCounter, "not matched", radius1, cvg->getEdge(be).getStart(),boundarynode.getId());
			edgeCounter++;

			//calculate new radius
			float radius2 = 0.f;
			for (int i = 1; i < edgepoints2.size(); i++) {
				radius2 += edgepoints2[i].getRadius();
			}
			radius2 /= edgepoints2.size() - 1;
            ConcreteVesselGraphEdge edge2 = ConcreteVesselGraphEdge(edgeCounter, "not matched", radius2, boundarynode.getId(), cvg->getEdge(be).getEnd());
			edgeCounter++;

            edge1.setEdgePoints(edgepoints1);
            edge2.setEdgePoints(edgepoints2);
			boundarynode.addEdge(edge2.getID());
			cvg->getNode(cvg->getEdge(be).getStart()).eraseEdge(be);
			cvg->getNode(cvg->getEdge(be).getStart()).addEdge(edge1.getID());
            cvg->addNode(boundarynode);
            cvg->addEdge(edge1);
            cvg->addEdge(edge2);
            cvg->removeEdge(be);
        }

		//tgtAssert(cvg->getNodes().size() == 6, "nodes");
		//tgtAssert(cvg->getEdges().size() == 5, "edges");

        //Matching
        //Finding startnode in TemplateVesselGraph
        size_t templateStartNode = tvg->getNodes()[0].getID();
        bool flag=true;
        while (flag) {
            flag = false;
            for (TemplateVesselGraphEdge& edge : tvg->getEdges()) {
                if (edge.getEnd() == templateStartNode) {
                    templateStartNode = edge.getStart();
                    flag = true;
                }
            }
        }
        //Finding startnode in ConcreteVesselGraph
        size_t concreteStartNode = cvg->getNodes()[0].getId();
        flag = true;
        while (flag) {
            flag = false;
            for (ConcreteVesselGraphEdge& edge : cvg->getEdges()) {
                if (edge.getEnd() == concreteStartNode) {
                    concreteStartNode = edge.getStart();
                    flag = true;
                }
            }
        }

        std::vector<size_t> currentTemplateNodes = std::vector<size_t>();//Queue
        currentTemplateNodes.push_back(templateStartNode);
        std::vector<size_t> currentConcreteNodes = std::vector<size_t>();//Queue
        currentConcreteNodes.push_back(concreteStartNode);
        std::vector<TemplateVesselGraphEdge> currentTemplateEdgesToMatch = std::vector<TemplateVesselGraphEdge>();
        std::vector<ConcreteVesselGraphEdge> currentConcreteEdgesToMatch = std::vector<ConcreteVesselGraphEdge>();
		while (currentTemplateNodes.size()) {

			//Match Node IDs
			cvg->changeNodeLabel(currentConcreteNodes[0], tvg->getNode(currentTemplateNodes[0]).getLabel());

			//get outgoing edges from current Node
			TemplateVesselGraphNode templateNode = tvg->getNode(currentTemplateNodes[0]);
			std::vector<size_t> templateEdgeIDs = templateNode.getEdgeIDs();
			for (size_t edge: templateEdgeIDs)
			{
				currentTemplateEdgesToMatch.push_back(tvg->getEdges()[edge]);
			}

			ConcreteVesselGraphNode concreteNode = cvg->getNode(currentConcreteNodes[0]);
			std::vector<size_t> concreteEdgeIDs = concreteNode.getEdgeIDs();
			for (size_t i = 0; i < concreteEdgeIDs.size(); i++)
			{
				currentConcreteEdgesToMatch.push_back(cvg->getEdge(concreteEdgeIDs[i]));
			}

			
            //sort edges using radius
			std::sort(currentConcreteEdgesToMatch.begin(), currentConcreteEdgesToMatch.end(), [&](const ConcreteVesselGraphEdge& e1,
				const ConcreteVesselGraphEdge& e2) {return e1.getRadius() < e2.getRadius();});
            std::sort(currentTemplateEdgesToMatch.begin(), currentTemplateEdgesToMatch.end(), [&](const TemplateVesselGraphEdge& e1,
                const TemplateVesselGraphEdge& e2) {return e1.getRadius() < e2.getRadius();});

            //match edge labels 
            for (int i = 0;i < currentTemplateEdgesToMatch.size();i++) {
				if (i >= currentConcreteEdgesToMatch.size()) {
					LWARNINGC("ConcreteVesselGraphCreator", "Matching Error. Check inverts and mirroring of the velocity volumes (VolumeListMultiChannelAdapter).");
					break;
				}
                cvg->changeEdgeLabel(currentConcreteEdgesToMatch[i].getID(), currentTemplateEdgesToMatch[i].getLabel());
                currentTemplateNodes.push_back(currentTemplateEdgesToMatch[i].getEnd());
                currentConcreteNodes.push_back(currentConcreteEdgesToMatch[i].getEnd());
            }

            //first node done, remove from queue
            currentTemplateNodes.erase(currentTemplateNodes.begin());
            currentConcreteNodes.erase(currentConcreteNodes.begin());
			currentTemplateEdgesToMatch = std::vector<TemplateVesselGraphEdge>();
			currentConcreteEdgesToMatch = std::vector<ConcreteVesselGraphEdge>();
        }
		/*
		for (ConcreteVesselGraphNode& node : cvg->getNodes()) {
			tgtAssert(node.getLabel().compare("not matched") != 0, "node label");
			std::cout << node.getDegree() << std::endl;
		}
		for (ConcreteVesselGraphEdge& edge : cvg->getEdges()) {
			tgtAssert(edge.getLabel().compare("not matched") != 0, "edge label");
		}
		*/
        concreteVesselGraphOutport.setData(cvg);

		
    }
}