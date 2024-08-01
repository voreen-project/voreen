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

#include "flowsimulationconfigport.h"

namespace voreen {

FlowSimulationConfigPort::FlowSimulationConfigPort(PortDirection direction, const std::string& id, const std::string& guiName, bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel)
    : GenericPort<FlowSimulationConfig>(direction, id, guiName, allowMultipleConnections, invalidationLevel) {
}

std::string FlowSimulationConfigPort::getClassName() const {
    return "FlowParametrizationPort"; // For backwards compatibility.
}

Port* FlowSimulationConfigPort::create(PortDirection direction, const std::string& id, const std::string& guiName) const {
    return new FlowSimulationConfigPort(direction,id,guiName);
}

tgt::col3 FlowSimulationConfigPort::getColorHint() const {
    return tgt::col3(64, 64, 255);
}

std::string FlowSimulationConfigPort::getContentDescription() const {
    std::stringstream strstr;
    strstr << Port::getContentDescription();
    if(hasData()) {
        strstr << std::endl << "Name: " << getData()->getName();
        strstr << std::endl << "Members: " << getData()->getFlowParameterSets().size();
        strstr << std::endl << "Indicators: " << getData()->getFlowIndicators().size();
    }
    return strstr.str();
}

std::string FlowSimulationConfigPort::getContentDescriptionHTML() const {
    std::stringstream strstr;
    strstr << Port::getContentDescriptionHTML();
    if(hasData()) {
        strstr << "<br>" << "Name: " << getData()->getName();
        strstr << "<br>" << "Members: " << getData()->getFlowParameterSets().size();
        strstr << "<br>" << "Indicators: " << getData()->getFlowIndicators().size();
    }
    return strstr.str();
}

} // namespace
