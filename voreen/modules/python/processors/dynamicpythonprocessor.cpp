/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include <regex>

namespace voreen {

const std::string DynamicPythonProcessor::loggerCat_("voreen.DynamicPythonProcessor");

DynamicPythonProcessor::DynamicPythonProcessor()
    : RenderProcessor()
    , portList_("portList", "Port List", true)
    , enabled_("enabled", "Enabled", true)
    , pythonProperty_("pythonScript", "Python Script")
    , valid_(false)
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
    //addProperty(pythonProperty_); // Don't add property here, since the editor is included as processor widget!
    pythonProperty_.setOwner(this); // ..but override owner!
    ON_CHANGE(pythonProperty_, DynamicPythonProcessor, onScriptChange);
}

DynamicPythonProcessor::~DynamicPythonProcessor() {
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

        // The following call will clear the console output.
        // Additionally, it set's the source code. Hence, we don't call it here.
        //getProcessorWidget()->updateFromProcessor();

        // Run script.
        pythonScript_.run(true);
    }
}

void DynamicPythonProcessor::serialize(Serializer& s) const {
    bool usePointerContentSerialization = s.getUsePointerContentSerialization();
    try {
        std::vector<const Property*> pyVec;
        pyVec.push_back(&pythonProperty_);
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
        pyVec.push_back(&pythonProperty_);
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
            if(iter->getName() == port->getID()) {
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
        Port* item = portItems_[newPort.getItemId()].get();
        Port* port = item->create(item->isInport() ? Port::INPORT : Port::OUTPORT, newPort.getName());
        addPort(port);
        port->initialize();
    }
}

void DynamicPythonProcessor::onScriptChange() {

    // Reset valid flag.
    valid_ = false;

    // Copy script to keep its ID.
    pythonScript_ = pythonProperty_.get();

    // Parse script and search for global modifications.
    std::string source = pythonScript_.getSource();

    const std::regex getPortDataRegex(R"(getPortData\(\s*(\"[^\"]*\")\s*\))");
    source = std::regex_replace(source, getPortDataRegex, "getPortData(\"" + getGuiName() + "\", $1)");

    const std::regex setPortDataRegex(R"(setPortData\(\s*(\"[^\"]*\"\s*,\s*[a-zA-Z0-9]+\s*)\))");
    source = std::regex_replace(source, setPortDataRegex, "setPortData(\"" + getGuiName() + "\", $1)");

    // Apply modification and compile.
    pythonScript_.setSource(source);
    valid_ = pythonScript_.compile();
}

void DynamicPythonProcessor::addPortItem(Port* port) {
    std::string name = port->getClassName() + (port->isInport() ? "-in" : "-out");
    portList_.addItem(name);
    portItems_.push_back(std::unique_ptr<Port>(port));
}



} // namespace voreen
