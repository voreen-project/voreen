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

#include "templatevesselgraphloader.h"

namespace voreen {

	TemplateVesselGraphLoader::TemplateVesselGraphLoader() : Processor(),
		vesselGraphOutput_(Port::OUTPORT, "vesselGraph", "TemplateVesselGraph Output")
	{
		addPort(vesselGraphOutput_);
	}

	Processor* TemplateVesselGraphLoader::create() const {
		return new TemplateVesselGraphLoader();
	}

	void TemplateVesselGraphLoader::process() {
		TemplateVesselGraph *templateVesselGraphOutput = new TemplateVesselGraph();

		TemplateVesselGraphNode mpa_begin = TemplateVesselGraphNode(0, "mpa_begin", END_NODE);
		TemplateVesselGraphNode mpa_bifurcation = TemplateVesselGraphNode(1, "mpa_bifurcation", BIFURCATION_NODE);
		TemplateVesselGraphNode rpa_begin = TemplateVesselGraphNode(2, "rpa_begin", BOUNDARY_NODE);
		TemplateVesselGraphNode rpa_end = TemplateVesselGraphNode(3, "rpa_end", END_NODE);
		TemplateVesselGraphNode lpa_begin = TemplateVesselGraphNode(4, "lpa_begin", BOUNDARY_NODE);
		TemplateVesselGraphNode lpa_end = TemplateVesselGraphNode(5, "lpa_end", END_NODE);

		TemplateVesselGraphEdge mpa = TemplateVesselGraphEdge(0, std::string("mpa"), 1.0, mpa_begin.getID(), mpa_bifurcation.getID());
		TemplateVesselGraphEdge mpa_rpa = TemplateVesselGraphEdge(1, std::string("mpa_rpa"), 0.75, mpa_bifurcation.getID(), rpa_begin.getID());
		TemplateVesselGraphEdge rpa = TemplateVesselGraphEdge(2, std::string("rpa"), 0.75, rpa_begin.getID(), rpa_end.getID());
		TemplateVesselGraphEdge mpa_lpa = TemplateVesselGraphEdge(3, std::string("mpa_lpa"), 0.55, mpa_bifurcation.getID(), lpa_begin.getID());
		TemplateVesselGraphEdge lpa = TemplateVesselGraphEdge(4, std::string("lpa"), 0.55, lpa_begin.getID(), lpa_end.getID());

		mpa_begin.addEdge(mpa.getID());
		mpa_bifurcation.addEdge(mpa_rpa.getID());
		mpa_bifurcation.addEdge(mpa_lpa.getID());
		rpa_begin.addEdge(rpa.getID());
		lpa_begin.addEdge(lpa.getID());

		// Build TemplateVesselGraph
		templateVesselGraphOutput->addNode(mpa_begin);
		templateVesselGraphOutput->addNode(mpa_bifurcation);
		templateVesselGraphOutput->addNode(rpa_begin);
		templateVesselGraphOutput->addNode(rpa_end);
		templateVesselGraphOutput->addNode(lpa_begin);
		templateVesselGraphOutput->addNode(lpa_end);

		templateVesselGraphOutput->addEdge(mpa);
		templateVesselGraphOutput->addEdge(mpa_rpa);
		templateVesselGraphOutput->addEdge(rpa);
		templateVesselGraphOutput->addEdge(mpa_lpa);
		templateVesselGraphOutput->addEdge(lpa);

		// Return TemplateVesselGraph
		vesselGraphOutput_.setData(templateVesselGraphOutput);
		return;
	}
}
