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

#include "cvglabelswapper.h"

namespace voreen {

	CvgLabelSwapper::CvgLabelSwapper() : Processor(),
		edgeLabel1("firstSwapTarget", "First swap target"),
		edgeLabel2("secondSwapTarget", "Second swap target"),
		concreteVesselGraphInport(Port::INPORT, "ConcreteVesselGraphInport", "Concrete Vessel Graph Input"),
		concreteVesselGraphOutport(Port::OUTPORT, "ConcreteVesselGraphOutport", "Concrete Vessel Graph Output")
	{
		addProperty(edgeLabel1);
		addProperty(edgeLabel2);
		addPort(concreteVesselGraphInport);
		addPort(concreteVesselGraphOutport);
	}

	Processor* CvgLabelSwapper::create() const {
		return new CvgLabelSwapper();
	}

	void CvgLabelSwapper::process() {
		const ConcreteVesselGraph* cvg = concreteVesselGraphInport.getData();

		ConcreteVesselGraph* newCvg = new ConcreteVesselGraph();
		for (ConcreteVesselGraphEdge edge : cvg->getEdges()) {
			newCvg->addEdge(edge);
		}
		for (ConcreteVesselGraphNode node : cvg->getNodes()) {
			newCvg->addNode(node);
		}

		//init label options
		for (int i = 0; i < newCvg->getEdges().size(); i++) {
			edgeLabel1.addOption(newCvg->getEdges()[i].getLabel(), newCvg->getEdges()[i].getLabel());
			edgeLabel2.addOption(newCvg->getEdges()[i].getLabel(), newCvg->getEdges()[i].getLabel());
		}

		std::string label1 = edgeLabel1.getValue();
		std::string label2 = edgeLabel2.getValue();

		for (int i = 0; i < cvg->getEdges().size(); i++) {
			if (cvg->getEdges()[i].getLabel().compare(label1) == 0) {
				newCvg->changeEdgeLabel(cvg->getEdges()[i].getID(),label2);
			}else if (cvg->getEdges()[i].getLabel().compare(label2) == 0) {
				newCvg->changeEdgeLabel(cvg->getEdges()[i].getID(), label1);
			}
		}

		concreteVesselGraphOutport.setData(newCvg);
	}
}