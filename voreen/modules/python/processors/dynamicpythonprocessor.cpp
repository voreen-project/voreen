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

namespace voreen {

const std::string DynamicPythonProcessor::loggerCat_("voreen.DynamicPythonProcessor");

DynamicPythonProcessor::DynamicPythonProcessor()
    : Processor()
    , portList_("portList", "Port List")
{
    // Add available ports.
    addPortItem(new VolumePort(Port::INPORT, ""));
    addPortItem(new VolumePort(Port::OUTPORT, ""));
    addPortItem(new GeometryPort(Port::INPORT, ""));
    addPortItem(new GeometryPort(Port::OUTPORT, ""));

    // Add properties.
    addProperty(portList_);
}

DynamicPythonProcessor::~DynamicPythonProcessor() {
}

Processor* DynamicPythonProcessor::create() const {
    return new DynamicPythonProcessor();
}

void DynamicPythonProcessor::initialize() {
    Processor::initialize();
}

void DynamicPythonProcessor::deinitialize() {
    Processor::deinitialize();
}

void DynamicPythonProcessor::process() {
}

// private methods
//

void DynamicPythonProcessor::onPortListChange() {

    // Check if instance was deleted.
    bool numInstancesChanged = portList_.getInstances().size() != numInstances_;
    if(numInstancesChanged) {
        // Handle removal.
        if(numInstances_ > portList_.getInstances().size() && selectedInstance_) {
            // Assumes that only the selected item can be removed!
            tgtAssert(numInstances_ == portList_.getInstances().size() + 1, "Only single instance removal allowed!");
            Port* port = getPort(portList_.getInstanceName(*selectedInstance_));
            if(port) {
                removePort(port);
                port->deinitialize();
                delete port;
            }
            selectedInstance_.reset();
        } else if(numInstances_ < portList_.getInstances().size() && selectedInstance_) {
            // Assumes that only the selected item can be added!
            tgtAssert(numInstances_ == portList_.getInstances().size() - 1, "Only single instance addition allowed!");
            Port* item = portItems_[selectedInstance_->itemId_].get();
            Port* port = item->create(item->isInport() ? Port::INPORT : Port::OUTPORT, portList_.getInstanceName(*selectedInstance_));
            port->initialize();
            addPort(port);
        }
        numInstances_ = portList_.getInstances().size();
    }

    // Hide old group.
    if(selectedInstance_) {
        // We need to reset here, because otherwise onFilterPropertyChange
        // will be triggered while the current instance is restored.
        selectedInstance_.reset();
    }

    // Show new group.
    boost::optional<InteractiveListProperty::Instance> currentInstance;
    if(portList_.getSelectedInstance() != -1) {
        currentInstance = portList_.getInstances()[portList_.getSelectedInstance()];
    }

    selectedInstance_ = currentInstance;
}

void DynamicPythonProcessor::addPortItem(Port* port) {
    std::string name = port->getClassName() + (port->isInport() ? " Inport" : " Outport");
    portList_.addItem(name);
    portItems_.push_back(std::unique_ptr<Port>(port));
}



} // namespace voreen
