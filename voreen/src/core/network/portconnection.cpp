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

#include "voreen/core/network/portconnection.h"

#include "voreen/core/ports/port.h"

namespace voreen {

PortConnection::PortEntry::PortEntry(Port* port)
    : port_(port)
{}

void PortConnection::PortEntry::serialize(Serializer& s) const {
    s.serialize("name", port_->getID());
    s.serialize("Processor", port_->getProcessor());
}

void PortConnection::PortEntry::deserialize(Deserializer& s) {
    std::string name;
    Processor* processor;

    s.deserialize("name", name);
    s.deserialize("Processor", processor);

    // Was processor not deserialized?
    if (!processor)
        // Cancel port search, since we are not able to proceed...
        return;

    std::vector<Port*> ports = processor->getPorts();
    for (size_t i=0; i<ports.size(); ++i) {
        if (ports[i]->getID() == name) {
            port_ = ports[i];
            break;
        }
    }
}

Port* PortConnection::PortEntry::getPort() const {
    return port_;
}

PortConnection::PortConnection(Port* outport, Port* inport)
    : outport_(outport)
    , inport_(inport)
{}

PortConnection::PortConnection()
    : outport_(0)
    , inport_(0)
{}

void PortConnection::serialize(Serializer& s) const {
    s.serialize("Outport", outport_);
    s.serialize("Inport", inport_);
}

void PortConnection::deserialize(Deserializer& s) {
    s.deserialize("Outport", outport_);
    s.deserialize("Inport", inport_);
}

void PortConnection::setOutport(Port* value) {
    outport_ = PortEntry(value);
}

Port* PortConnection::getOutport() const {
    return outport_.getPort();
}

void PortConnection::setInport(Port* value) {
    inport_ = PortEntry(value);
}

Port* PortConnection::getInport() const {
    return inport_.getPort();
}

} // namespace
