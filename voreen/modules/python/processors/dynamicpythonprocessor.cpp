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

#include "dynamicpythonprocessor.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/ports/renderport.h"

namespace voreen {

const std::string DynamicPythonProcessor::loggerCat_("voreen.DynamicPythonProcessor");

DynamicPythonProcessor::DynamicPythonProcessor()
    : RenderProcessor()
    , portList_("portList", "Port List", true)
{
    // Add available ports.
    addPortItem(new VolumePort(Port::INPORT, ""));
    addPortItem(new VolumePort(Port::OUTPORT, ""));
    //addPortItem(new GeometryPort(Port::INPORT, ""));
    //addPortItem(new GeometryPort(Port::OUTPORT, ""));
    addPortItem(new RenderPort(Port::INPORT, ""));
    addPortItem(new RenderPort(Port::OUTPORT, ""));

    // Add properties.
    addProperty(portList_);
    ON_CHANGE(portList_, DynamicPythonProcessor, onPortListChange);

    /*
    // TODO: define naming schema!
    // Name Generator spitting out "<PortType>(<instance>)"
    InteractiveListProperty::NameGenerator nameGenerator =
            [this] (const InteractiveListProperty::Instance& instance) {
                std::string name = portList_.getItems()[instance.itemId_];
                name += "(" + std::to_string(portInstances_[name].size()) + ")";
                return name;
            };

    portList_.setNameGenerator(nameGenerator);
    */
}

DynamicPythonProcessor::~DynamicPythonProcessor() {
    portInstances_.clear();
    while(!getPorts().empty()) {
        Port* port = getPorts().front();
        removePort(port);
        delete port;
    }
}

Processor* DynamicPythonProcessor::create() const {
    return new DynamicPythonProcessor();
}

void DynamicPythonProcessor::initialize() {
    RenderProcessor::initialize();
}

void DynamicPythonProcessor::deinitialize() {
    RenderProcessor::deinitialize();
}

void DynamicPythonProcessor::process() {
}

void DynamicPythonProcessor::serialize(Serializer& s) const {
    RenderProcessor::serialize(s);
    //portList_.serialize(s);
}
void DynamicPythonProcessor::deserialize(Deserializer& s) {
    //portList_.deserialize(s);
    RenderProcessor::deserialize(s);
    onPortListChange();
}

// private methods
//

void DynamicPythonProcessor::onPortListChange() {

    // Gather orphaned ports.
    std::vector<Port*> deletedPorts;
    std::vector<InteractiveListProperty::Instance> newPorts = portList_.getInstances();
    for(Port* port : getPorts()) {
        bool found = false;
        for(auto iter = newPorts.begin(); iter != newPorts.end(); iter++) {
            if(portList_.getInstanceName(*iter) == port->getID()) {
                found = true;
                newPorts.erase(iter);
                break;
            }
        }

        if(!found) {
            deletedPorts.push_back(port);
        }
    }

    // Delete orphaned ports.
    for(Port* port : deletedPorts) {
        //std::vector<Port*>& ports = portInstances_[port->getClassName() + (port->isInport() ? "-in" : "-out")]; // TODO: remove!!
        //ports.erase(std::find(ports.begin(), ports.end(), port));
        port->deinitialize();
        removePort(port);
        delete port;
    }

    // Create new Ports.
    for (auto& newPort : newPorts) {
        Port* item = portItems_[newPort.itemId_].get();
        Port* port = item->create(item->isInport() ? Port::INPORT : Port::OUTPORT, portList_.getInstanceName(newPort));
        addPort(port);
        port->initialize();
        //portInstances_[portList_.getItems()[newPort.itemId_]].push_back(port);
    }
}

void DynamicPythonProcessor::addPortItem(Port* port) {
    std::string name = port->getClassName() + (port->isInport() ? "-in" : "-out");
    portList_.addItem(name);
    portItems_.push_back(std::unique_ptr<Port>(port));
}



} // namespace voreen
