/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "voreen/core/ports/coprocessorport.h"

namespace voreen {

CoProcessorPort::CoProcessorPort(PortDirection direction, const std::string& id, const std::string& guiName, const bool allowMultipleConnections,
                Processor::InvalidationLevel invalidationLevel)
    : Port(direction, id, guiName, allowMultipleConnections, invalidationLevel)
{}

std::vector<Processor*> CoProcessorPort::getConnectedProcessors() const {
    if (isOutport()) {
        tgtAssert(false, "Called CoProcessorPort::getConnectedProcessors on outport!");
        LERRORC("voreen.coprocessorport", "Called CoProcessorPort::getConnectedProcessors on outport!");
        return std::vector<Processor*>();
    }
    else {
        std::vector<Processor*> processors;
        for (size_t i = 0; i < connectedPorts_.size(); ++i) {
            processors.push_back(connectedPorts_[i]->getProcessor());
        }
        return processors;
    }
}

Processor* CoProcessorPort::getConnectedProcessor() const {
    if (isOutport()) {
        tgtAssert(false, "Called CoProcessorPort::getConnectedProcessor on outport!");
        LERRORC("voreen.coprocessorport", "Called CoProcessorPort::getConnectedProcessor on outport!");
        return 0;
    }
    else {
        if (isConnected())
            return connectedPorts_[0]->getProcessor();
        else
            return 0;
    }
}

bool CoProcessorPort::isReady() const {
    return ((getConnectedProcessor() != 0) && getConnectedProcessor()->isReady());
}

tgt::col3 CoProcessorPort::getColorHint() const {
    return tgt::col3(0, 255, 0);
}

std::string CoProcessorPort::getContentDescription() const {
    std::stringstream strstr;
    strstr << getGuiName() << std::endl
           << "Type: " << getClassName();
    return strstr.str();
}

std::string CoProcessorPort::getContentDescriptionHTML() const {
    std::stringstream strstr;
    strstr << "<center><font><b>" << getGuiName() << "</b></font></center>"
           << "Type: " << getClassName();
    return strstr.str();
}

} // namespace
