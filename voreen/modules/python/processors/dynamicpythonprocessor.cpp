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
#include "voreen/core/ports/renderport.h"
#include "voreen/core/processors/processorwidget.h"

namespace voreen {

const std::string DynamicPythonProcessor::loggerCat_("voreen.DynamicPythonProcessor");

DynamicPythonProcessor::DynamicPythonProcessor()
    : RenderProcessor()
    , portList_("portList", "Port List", true)
    , enabled_("enabled", "Enabled", true)
    , pythonScript_("pythonScript", "Python Script")
{
    // Add available ports.
    addPortItem(new VolumePort(Port::INPORT, ""));
    addPortItem(new VolumePort(Port::OUTPORT, ""));
    addPortItem(new RenderPort(Port::INPORT, ""));
    addPortItem(new RenderPort(Port::OUTPORT, "")); // Add properties to be able to set format?

    // Add properties.
    addProperty(portList_);
    ON_CHANGE(portList_, DynamicPythonProcessor, onPortListChange);

    addProperty(enabled_);
    //addProperty(pythonScript_); // Don't add property here, since the editor is included as processor widget!
    pythonScript_.setOwner(this); // ..but override owner!
    ON_CHANGE(pythonScript_, DynamicPythonProcessor, onScriptChange);

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
    portList_.clear(); // Will trigger onPortListChange an delete remaining ports.
}

Processor* DynamicPythonProcessor::create() const {
    return new DynamicPythonProcessor();
}

void DynamicPythonProcessor::initialize() {
    RenderProcessor::initialize();
    getProcessorWidget()->updateFromProcessor();
}

void DynamicPythonProcessor::deinitialize() {
    RenderProcessor::deinitialize();
}

bool DynamicPythonProcessor::isReady() const {

    if(enabled_.get() && !valid_) {
        setNotReadyErrorMessage("Script contains invalid code!");
        return false;
    }

    return RenderProcessor::isReady();
}

void DynamicPythonProcessor::process() {
    if(enabled_.get()) {
        PythonScript script = pythonScript_.get();
        script.run(true);
    }
}

void DynamicPythonProcessor::serialize(Serializer& s) const {
    bool usePointerContentSerialization = s.getUsePointerContentSerialization();
    try {
        std::vector<const Property*> pyVec;
        pyVec.push_back(&pythonScript_);
        s.setUsePointerContentSerialization(true);
        s.serialize("PriorityProperty", pyVec, "Property");
    } catch (SerializationException &e) {
        LWARNING(std::string("Python Script serialization failed: ") + e.what());
    }
    s.setUsePointerContentSerialization(usePointerContentSerialization);

    RenderProcessor::serialize(s);
}
void DynamicPythonProcessor::deserialize(Deserializer& s) {
    RenderProcessor::deserialize(s);

    bool usePointerContentSerialization = s.getUsePointerContentSerialization();
    try {
        std::vector<Property*> pyVec;
        pyVec.push_back(&pythonScript_);
        s.setUsePointerContentSerialization(true);
        s.deserialize("PriorityProperty", pyVec, "Property");
    } catch (SerializationNoSuchDataException&) {
        s.removeLastError();
    }
    s.setUsePointerContentSerialization(usePointerContentSerialization);

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
    }
}

void DynamicPythonProcessor::onScriptChange() {

    // Reset valid flag.
    valid_ = false;

    const std::string functionName = /*"s/g"*/"etPortData";

    // Parse script and search for global modifications.
    const std::string& source = pythonScript_.get().getSource();
    size_t portDataPos = source.find(functionName);
    while(portDataPos != std::string::npos) {

        size_t openQuote = source.find_first_of('"', portDataPos) + 1;
        size_t endQuote  = source.find_first_of('"', openQuote);

        std::string processorArgumentString = source.substr(openQuote, endQuote-openQuote);
        if(processorArgumentString != getGuiName()) {
            LERROR("Within a processor-local python script only the processor's own ports can be modified!");
            return;
        }

        portDataPos = source.find(functionName, portDataPos+1);
    }

    valid_ = true;
    /*
    PythonScript script = pythonScript_.get();
    script.compile(true);

    if(valid_) {
        LINFO("Script successfully saved!");
    }
    else {
        LERROR(script.getLog());
    }
    */

    invalidate(Processor::INVALID_PROGRAM);
}

void DynamicPythonProcessor::addPortItem(Port* port) {
    std::string name = port->getClassName() + (port->isInport() ? "-in" : "-out");
    portList_.addItem(name);
    portItems_.push_back(std::unique_ptr<Port>(port));
}



} // namespace voreen
