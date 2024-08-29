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

#include "abstractvisualizationcreator.h"
#include "voreen/core/processors/asynccomputeprocessor.h"
#include "modules/plotting/datastructures/plotcell.h"
#include "modules/plotting/datastructures/plotdata.h"
#include "modules/plotting/datastructures/plotrow.h"
#include "modules/flowsimulation/utils/utils.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
namespace voreen {

    AbstractVisualizationCreator::AbstractVisualizationCreator() : Processor(),
        templateVesselGraphInport(Port::INPORT, "TemplateVesselGraph", "Template Vessel Graph Input"),
        controlVolumeInport(Port::INPORT, "Control Volume", "control Volume List Port"),
        controlConcreteVesselGraphInport(Port::INPORT, "Control Concrete Vesselgraph", "control Concrete Vesselgraph Port"),
        diseasedVolumeInport(Port::INPORT, "Diseased Volume", "diseased Volume List Port"),
        diseasedConcreteVesselGraphInport(Port::INPORT, "Diseased Concrete Vesselgraph", "diseased Concrete Vesselgraph Port"),
        treatedVolumeInport(Port::INPORT, "Treated Volume", "treated Volume List Port"),
        treatedConcreteVesselGraphInport(Port::INPORT, "Treated Concrete Vesselgraph", "treated Concrete Vesselgraph Port"),
        plotOutport(Port::OUTPORT, "Plot", "Plot Port"),
        controlPointsOutport(Port::OUTPORT, "Control Points", "control Points Port"),
        diseasedPointsOutport(Port::OUTPORT, "diseased Points", "diseased Points Port"),
        treatedPointsOutport(Port::OUTPORT, "treated Points", "treated Points Port"),
        controlParameterOutport(Port::OUTPORT, "control parametrization outport", "Control Flow Parametrization Output"),
        diseasedParameterOutport(Port::OUTPORT, "diseased parametrization outport", "Diseased Flow Parametrization Output"),
        treatedParameterOutport(Port::OUTPORT, "treated parametrization outport", "Treated Flow Parametrization Output"),
        edgeLabel("edge label","Edge Label"),
        edgeSlider("edge slider", "Edge Slider"),
        calculateButton("calculate button","Calculate max StdDeviation"),
		controlEdgeIds("controlEdgeIds", "Hidden Control EdgeIds"),
		diseasedEdgeIds("diseasedEdgeIds", "Hidden Diseased EdgeIds"),
		treatedEdgeIds("treatedEdgeIds", "Hidden Treated EdgeIds")

    {
        addPort(templateVesselGraphInport);
        addPort(controlVolumeInport);
        addPort(controlConcreteVesselGraphInport);
        addPort(diseasedVolumeInport);
        addPort(diseasedConcreteVesselGraphInport);
        addPort(treatedVolumeInport);
        addPort(treatedConcreteVesselGraphInport);
        addPort(plotOutport);
        addPort(controlPointsOutport);
        addPort(diseasedPointsOutport);
        addPort(treatedPointsOutport);
        addPort(controlParameterOutport);
        addPort(diseasedParameterOutport);
        addPort(treatedParameterOutport);
        ON_CHANGE(templateVesselGraphInport, AbstractVisualizationCreator, setLabels);
        addProperty(edgeLabel);
        ON_CHANGE(edgeLabel, AbstractVisualizationCreator, remember);
        addProperty(edgeSlider);
        edgeSlider.set(0.5f);
        addProperty(calculateButton);
        ON_CHANGE(calculateButton,AbstractVisualizationCreator,calculate);
		addProperty(controlEdgeIds);
		controlEdgeIds.setVisibleFlag(false);
		addProperty(diseasedEdgeIds);
		diseasedEdgeIds.setVisibleFlag(false);
		addProperty(treatedEdgeIds);
		treatedEdgeIds.setVisibleFlag(false);
    }

    Processor* AbstractVisualizationCreator::create() const {
        return new AbstractVisualizationCreator();
    }

    void AbstractVisualizationCreator::process() {
        //load volumelists and concretevesselgraphs
        auto controlVolumes = controlVolumeInport.getThreadSafeData();
        if (!controlVolumes || controlVolumes->empty()) {
            throw InvalidInputException("No volumes", InvalidInputException::S_IGNORE);
        }
        const ConcreteVesselGraph* controlCVG = controlConcreteVesselGraphInport.getData();

        auto diseasedVolumes = diseasedVolumeInport.getThreadSafeData();
        if (!diseasedVolumes || diseasedVolumes->empty()) {
            throw InvalidInputException("No volumes", InvalidInputException::S_IGNORE);
        }
        const ConcreteVesselGraph* diseasedCVG = diseasedConcreteVesselGraphInport.getData();

        auto treatedVolumes = treatedVolumeInport.getThreadSafeData();
        if (!treatedVolumes || treatedVolumes->empty()) {
            throw InvalidInputException("No volumes", InvalidInputException::S_IGNORE);
        }
        const ConcreteVesselGraph* treatedCVG = treatedConcreteVesselGraphInport.getData();

        //init PlotData
        PlotData* output = new PlotData(1, 3);
        output->setColumnLabel(0, "Timesteps");
        output->setColumnLabel(1, "control velocity");
        output->setColumnLabel(2, "diseased velocity");
        output->setColumnLabel(3, "treated velocity");

        std::vector<PlotCellValue> values;

        //selected edgelabel
        std::string label = edgeLabel.getValue();

        //selected position on the edge
        float interpolatePosition = edgeSlider.get();
        rememberValues[label] = interpolatePosition;
        
        //positions 
        std::tuple<tgt::vec3, tgt::vec3,float> controlInterpolation = std::make_tuple(tgt::vec3::zero, tgt::vec3::zero,0.0);
        std::tuple<tgt::vec3, tgt::vec3, float>   diseasedInterpolation = std::make_tuple(tgt::vec3::zero, tgt::vec3::zero, 0.0);
        std::tuple<tgt::vec3, tgt::vec3, float>   treatedInterpolation = std::make_tuple(tgt::vec3::zero, tgt::vec3::zero, 0.0);
        //radius 
        float controlAvgRadius, diseasedAvgRadius, treatedAvgRadius;

        //Control
        std::vector<ConcreteVesselGraphEdge> edges = controlCVG->getEdges();

		//add edgeIds to invisible property
		controlEdgeIds.reset();
		int index = 0;
		int maxId = 0;
		//find out max, to adapt to connectedcomponentselector
		for (int i = 0; i < edges.size(); i++) {
			if (edges[i].getID() > maxId) { maxId = edges[i].getID(); }
		}
		//add rows to property and find clicked edge
		for (int i = 1; i <= maxId; i++) {
			controlEdgeIds.addRow(std::to_string(i));
			for (int j = 0; j < edges.size(); j++) {
				if (i == edges[j].getID()) {
					if (edges[j].getLabel().compare(label) == 0) { index = i-1; }
				}
				
			}
		}

		//select the id from the edge chosen by user 
		std::vector<int> select = { index };
		controlEdgeIds.setSelectedRowIndices(select);

		//Diseased
		edges = diseasedCVG->getEdges();

		//add rows to property and find clicked edge 
		for (int i = 1; i <= maxId; i++) {
			diseasedEdgeIds.addRow(std::to_string(i));
			for (int j = 0; j < edges.size(); j++) {
				if (i == edges[j].getID()) {
					if (edges[j].getLabel().compare(label) == 0) { index = i - 1; }
				}

			}
		}
		//select the id from the edge chosen by user
		select = { index };
		diseasedEdgeIds.setSelectedRowIndices(select);

		//Treated
		edges = treatedCVG->getEdges();

		//add rows to property and find clicked edge 
		for (int i = 1; i <= maxId; i++) {
			treatedEdgeIds.addRow(std::to_string(i));
			for (int j = 0; j < edges.size(); j++) {
				if (i == edges[j].getID()) {
					if (edges[j].getLabel().compare(label) == 0) { index = i - 1; }
				}

			}
		}
		//select the id from the edge chosen by user
		select = { index };
		treatedEdgeIds.setSelectedRowIndices(select);


		
		edges = controlCVG->getEdges();
        //find edge in concretevesselgraph and interpolate
        for (int i = 0;i < edges.size();i++) {
            if (edges[i].getLabel().compare(label) == 0) {
                controlInterpolation = edges[i].interpolate(interpolatePosition);
                break;
            }
        }

        edges = diseasedCVG->getEdges();
        for (int i = 0;i < edges.size();i++) {
            if (edges[i].getLabel().compare(label) == 0) {
                diseasedInterpolation =edges[i].interpolate(interpolatePosition);
                break;
            }
        }

        edges = treatedCVG->getEdges();
        for (int i = 0;i < edges.size();i++) {
            if (edges[i].getLabel().compare(label) == 0) {
                treatedInterpolation = edges[i].interpolate(interpolatePosition);
                break;
            }
        }

        //Set FlowParameter
        FlowSimulationConfig* controlFlowParameterSetEnsemble = new FlowSimulationConfig("Control");
        FlowIndicator indicator;
        indicator.type_ = FlowIndicatorType::FIT_MEASURE;
        indicator.center_ = std::get<0>(controlInterpolation);
        indicator.normal_ = std::get<1>(controlInterpolation);
        indicator.radius_ = std::get<2>(controlInterpolation);
        controlFlowParameterSetEnsemble->addFlowIndicator(indicator);
        controlParameterOutport.setData(controlFlowParameterSetEnsemble);

        FlowSimulationConfig* diseasedFlowParameterSetEnsemble = new FlowSimulationConfig("Diseased");
        indicator.type_ = FlowIndicatorType::FIT_MEASURE;
        indicator.center_ = std::get<0>(diseasedInterpolation);
        indicator.normal_ = std::get<1>(diseasedInterpolation);
        indicator.radius_ = std::get<2>(diseasedInterpolation);
        diseasedFlowParameterSetEnsemble->addFlowIndicator(indicator);
        diseasedParameterOutport.setData(diseasedFlowParameterSetEnsemble);

        FlowSimulationConfig* treatedFlowParameterSetEnsemble = new FlowSimulationConfig("Treated");
        indicator.type_ = FlowIndicatorType::FIT_MEASURE;
        indicator.center_ = std::get<0>(treatedInterpolation);
        indicator.normal_ = std::get<1>(treatedInterpolation);
        indicator.radius_ = std::get<2>(treatedInterpolation);
        treatedFlowParameterSetEnsemble->addFlowIndicator(indicator);
        treatedParameterOutport.setData(treatedFlowParameterSetEnsemble);

        size_t minSize = std::min(controlVolumes->size(),diseasedVolumes->size());
        minSize = std::min(minSize, treatedVolumes->size());
        //sample data
        for (int i = 0;i < minSize; i++) {
            const VolumeBase* controlVolume = controlVolumes->at(i);
            const VolumeBase* diseasedVolume = diseasedVolumes->at(i);
            const VolumeBase* treatedVolume = treatedVolumes->at(i);
           
            std::vector<tgt::vec3> controlSamples = utils::sampleDisk(controlVolume,std::get<0>(controlInterpolation), std::get<1>(controlInterpolation),std::get<2>(controlInterpolation), false);
            std::vector<tgt::vec3> diseasedSamples = utils::sampleDisk(diseasedVolume, std::get<0>(diseasedInterpolation), std::get<1>(diseasedInterpolation), std::get<2>(diseasedInterpolation), false);
            std::vector<tgt::vec3> treatedSamples = utils::sampleDisk(treatedVolume, std::get<0>(treatedInterpolation), std::get<1>(treatedInterpolation), std::get<2>(treatedInterpolation), false);
            tgt::vec3 controlAccumVelocity = tgt::vec3::zero;
            tgt::vec3 diseasedAccumVelocity = tgt::vec3::zero;
            tgt::vec3 treatedAccumVelocity = tgt::vec3::zero;
                
            for (const auto& sample : controlSamples) {
                controlAccumVelocity += sample/(float)controlSamples.size();
            }
            for (const auto& sample : diseasedSamples) {
                diseasedAccumVelocity += sample / (float)diseasedSamples.size();
            }
            for (const auto& sample : treatedSamples) {
                treatedAccumVelocity += sample / (float)treatedSamples.size();
            }

                //add data
            values.push_back(PlotCellValue(i));
            values.push_back(PlotCellValue(tgt::length(controlAccumVelocity)));
            values.push_back(PlotCellValue(tgt::length(diseasedAccumVelocity)));
            values.push_back(PlotCellValue(tgt::length(treatedAccumVelocity)));
            output->insert(values);
            values.clear();

            }
        //Set pointsegmentlists

        PointSegmentListGeometryVec3* controlEdgeGeom = new PointSegmentListGeometryVec3();
        for (auto& edge : controlCVG->getEdges()) {
            std::vector<tgt::vec3> segment;
            for (auto& edgePoint : edge.getEdgePoints()) {
                segment.push_back(edgePoint.getPos());
            }
            controlEdgeGeom->addSegment(segment);
        }
        controlPointsOutport.setData(controlEdgeGeom);
        
        PointSegmentListGeometryVec3* diseasedEdgeGeom = new PointSegmentListGeometryVec3();
        for (auto& edge : diseasedCVG->getEdges()) {
            std::vector<tgt::vec3> segment;
            for (auto& edgePoint : edge.getEdgePoints()) {
                segment.push_back(edgePoint.getPos());
            }
            diseasedEdgeGeom->addSegment(segment);
        }
        diseasedPointsOutport.setData(diseasedEdgeGeom);

        PointSegmentListGeometryVec3* treatedEdgeGeom = new PointSegmentListGeometryVec3();
        for (auto& edge : treatedCVG->getEdges()) {
            std::vector<tgt::vec3> segment;
            for (auto& edgePoint : edge.getEdgePoints()) {
                segment.push_back(edgePoint.getPos());
            }
            treatedEdgeGeom->addSegment(segment);
        }
        treatedPointsOutport.setData(treatedEdgeGeom);
    
    plotOutport.setData(output);
    }

    void AbstractVisualizationCreator::setLabels() {
        const TemplateVesselGraph* tvg=templateVesselGraphInport.getData();
        if (tvg != NULL) {
            std::vector<TemplateVesselGraphEdge> edges = tvg->getEdges();
            for (int i = 0;i < edges.size();i++) {
                rememberValues[edges[i].getLabel()] = 0.5f;
                edgeLabel.addOption(edges[i].getLabel(), edges[i].getLabel());
            }
        }
    }
    void AbstractVisualizationCreator::remember() {
        edgeSlider.set(rememberValues[edgeLabel.getValue()]);
    }

	void AbstractVisualizationCreator::selectIdInProperty(std::string edgelabel) {
		
	}
    void AbstractVisualizationCreator::calculate() {
        if (templateVesselGraphInport.hasData() &&
            controlVolumeInport.hasData() &&
            controlConcreteVesselGraphInport.hasData() &&
            diseasedVolumeInport.hasData() &&
            diseasedConcreteVesselGraphInport.hasData() &&
            treatedVolumeInport.hasData() &&
            treatedConcreteVesselGraphInport.hasData()) {

            const TemplateVesselGraph* tvg = templateVesselGraphInport.getData();
            auto controlVolumes = controlVolumeInport.getThreadSafeData();
            if (!controlVolumes || controlVolumes->empty()) {
                throw InvalidInputException("No volumes", InvalidInputException::S_IGNORE);
            }
            const ConcreteVesselGraph* controlCVG = controlConcreteVesselGraphInport.getData();

            auto diseasedVolumes = diseasedVolumeInport.getThreadSafeData();
            if (!diseasedVolumes || diseasedVolumes->empty()) {
                throw InvalidInputException("No volumes", InvalidInputException::S_IGNORE);
            }
            const ConcreteVesselGraph* diseasedCVG = diseasedConcreteVesselGraphInport.getData();

            auto treatedVolumes = treatedVolumeInport.getThreadSafeData();
            if (!treatedVolumes || treatedVolumes->empty()) {
                throw InvalidInputException("No volumes", InvalidInputException::S_IGNORE);
            }
            const ConcreteVesselGraph* treatedCVG = treatedConcreteVesselGraphInport.getData();
            std::vector<TemplateVesselGraphEdge> templateEdges = tvg->getEdges();
            for (int i = 0;i < templateEdges.size();i++) {
                //edge to inspect
                std::string label = templateEdges[i].getLabel();
                float value,maxDeviation=-1;
                for (int j = 0;j < 11;j++) {

                    //selected position on the edge
                    float interpolatePosition = 0.1*j;

                    //positions 
                    std::tuple<tgt::vec3, tgt::vec3, float> controlInterpolation = std::make_tuple(tgt::vec3::zero, tgt::vec3::zero, 0.0);
                    std::tuple<tgt::vec3, tgt::vec3, float>   diseasedInterpolation = std::make_tuple(tgt::vec3::zero, tgt::vec3::zero, 0.0);
                    std::tuple<tgt::vec3, tgt::vec3, float>   treatedInterpolation = std::make_tuple(tgt::vec3::zero, tgt::vec3::zero, 0.0);
                    //radius 
                    float controlAvgRadius, diseasedAvgRadius, treatedAvgRadius;

                    //find edge in concretevesselgraph and interpolate
                    std::vector<ConcreteVesselGraphEdge> edges = controlCVG->getEdges();
                    for (int i = 0;i < edges.size();i++) {
                        if (edges[i].getLabel().compare(label) == 0) {
                            controlInterpolation = edges[i].interpolate(interpolatePosition);
                            break;
                        }
                    }

                    edges = diseasedCVG->getEdges();
                    for (int i = 0;i < edges.size();i++) {
                        if (edges[i].getLabel().compare(label) == 0) {
                            diseasedInterpolation = edges[i].interpolate(interpolatePosition);
                            break;
                        }
                    }

                    edges = treatedCVG->getEdges();
                    for (int i = 0;i < edges.size();i++) {
                        if (edges[i].getLabel().compare(label) == 0) {
                            treatedInterpolation = edges[i].interpolate(interpolatePosition);
                            break;
                        }
                    }

                    size_t minSize = std::min(controlVolumes->size(), diseasedVolumes->size());
                    minSize = std::min(minSize, treatedVolumes->size());
                    //sample data
                    float stdDeviation = 0;
                    for (int i = 0;i < minSize; i++) {
                        const VolumeBase* controlVolume = controlVolumes->at(i);
                        const VolumeBase* diseasedVolume = diseasedVolumes->at(i);
                        const VolumeBase* treatedVolume = treatedVolumes->at(i);

                        std::vector<tgt::vec3> controlSamples = utils::sampleDisk(controlVolume, std::get<0>(controlInterpolation), std::get<1>(controlInterpolation), std::get<2>(controlInterpolation), false);
                        std::vector<tgt::vec3> diseasedSamples = utils::sampleDisk(diseasedVolume, std::get<0>(diseasedInterpolation), std::get<1>(diseasedInterpolation), std::get<2>(diseasedInterpolation), false);
                        std::vector<tgt::vec3> treatedSamples = utils::sampleDisk(treatedVolume, std::get<0>(treatedInterpolation), std::get<1>(treatedInterpolation), std::get<2>(treatedInterpolation), false);
                        tgt::vec3 controlAccumVelocity = tgt::vec3::zero;
                        tgt::vec3 diseasedAccumVelocity = tgt::vec3::zero;
                        tgt::vec3 treatedAccumVelocity = tgt::vec3::zero;

                        for (const auto& sample : controlSamples) {
                            controlAccumVelocity += sample / (float)controlSamples.size();
                        }
                        for (const auto& sample : diseasedSamples) {
                            diseasedAccumVelocity += sample / (float)diseasedSamples.size();
                        }
                        for (const auto& sample : treatedSamples) {
                            treatedAccumVelocity += sample / (float)treatedSamples.size();
                        }
                        //calculate standard Deviation
                        float controlValue = tgt::length(controlAccumVelocity);
                        float diseasedValue = tgt::length(diseasedAccumVelocity);
                        float treatedValue = tgt::length(treatedAccumVelocity);
                        float mean = (controlValue + diseasedValue + treatedValue) / 3;
                        stdDeviation += std::sqrt(((controlValue - mean) * (controlValue - mean) +
                                (diseasedValue - mean) * (diseasedValue - mean) +
                                (treatedValue - mean) * (treatedValue - mean)) / 3);
                    }
                    //update maxDeviation
                    if (stdDeviation > maxDeviation) {
                        value = interpolatePosition;
                        maxDeviation = stdDeviation;
                    }
                   
                } 
                rememberValues[label] = value;
            }
            remember();
        }

    }
}