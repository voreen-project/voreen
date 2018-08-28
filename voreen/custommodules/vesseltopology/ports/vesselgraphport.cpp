/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "vesselgraphport.h"

namespace voreen {

VesselGraphPort::VesselGraphPort(PortDirection direction, const std::string& name, const std::string& guiName, bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel)
    : GenericPort<VesselGraph>(direction, name, guiName, allowMultipleConnections, invalidationLevel)
{
}

VesselGraphPort::~VesselGraphPort() {
}
Port* VesselGraphPort::create(PortDirection direction, const std::string& id, const std::string& guiName) const
{
    return new VesselGraphPort(direction,id,guiName);
}
tgt::col3 VesselGraphPort::getColorHint() const {
    return tgt::col3(128, 255, 0);
}

std::string VesselGraphPort::getContentDescription() const {
    std::stringstream strstr;
    //port values
    strstr  << Port::getContentDescription();

    if (hasData()) {
        const VesselGraph* graph = getData();

        strstr << std::endl << "Nodes: " << graph->getNodes().size();
        strstr << std::endl << "Edges: " << graph->getEdges().size();
        strstr << std::endl << "Bounding Box: " << graph->getBounds().getLLF() << " x " << graph->getBounds().getURB();
    }
    return strstr.str();
}

std::string VesselGraphPort::getContentDescriptionHTML() const {
    return getContentDescription();
}
} // namespace voreen
