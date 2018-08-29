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

#include "voreen/core/ports/textport.h"
#include "voreen/core/datastructures/rendertarget/rendertarget.h"
#include "voreen/core/processors/renderprocessor.h"

namespace voreen {

const std::string TextPort::loggerCat_("voreen.TextPort");

TextPort::TextPort(PortDirection direction, const std::string& name, const std::string& guiName,
                   bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel)
    : Port(direction, name, guiName, allowMultipleConnections, invalidationLevel)
{
}

TextPort::~TextPort() {
}

std::string TextPort::getContentDescription() const {
    std::stringstream strstr;
    strstr << Port::getContentDescription() << std::endl
            << "Content: \"" << getData() << "\"";
    return strstr.str();
}

std::string TextPort::getContentDescriptionHTML() const {
    std::stringstream strstr;
    strstr  << Port::getContentDescriptionHTML() << "<br>"
            << "Content: &quot;" << getData() << "&quot;";
    return strstr.str();
}

void TextPort::forwardData() const{
    for(std::vector<Port*>::const_iterator it = forwardPorts_.begin(); it != forwardPorts_.end(); ++it){
        dynamic_cast<TextPort*>(*it)->setData(getData());
    }
}

void TextPort::setData(const std::string& data) {
    tgtAssert(isOutport(), "called setData on inport!");
    notifyDataWillChange();
    portData_ = data;
    notifyDataHasChanged();
    invalidatePort();
}

std::string TextPort::getData() const {
    if (isOutport())
        return portData_;
    else {
        for (size_t i = 0; i < connectedPorts_.size(); ++i) {
            if (!connectedPorts_[i]->isOutport())
                continue;

            TextPort* p = static_cast< TextPort* >(connectedPorts_[i]);
            if (p)
                return p->getData();
        }
    }
    return "";
}

std::vector<std::string> TextPort::getAllData() const {
    if (isOutport()) {
        std::vector<std::string> allData;
        allData.push_back(portData_);
        return allData;
    }
    else {
        std::vector<std::string> allData;
        for (size_t i = 0; i < connectedPorts_.size(); ++i) {
            if (connectedPorts_[i]->isOutport() == false)
                continue;
            TextPort* p = static_cast<TextPort*>(connectedPorts_[i]);
            allData.push_back(p->getData());
        }
        return allData;
    }
}

std::vector<const TextPort*> TextPort::getConnected() const {
    std::vector<const TextPort*> ports;
    for (size_t i = 0; i < connectedPorts_.size(); ++i) {
        TextPort* p = static_cast<TextPort*>(connectedPorts_[i]);

        ports.push_back(p);
    }
    return ports;
}

bool TextPort::isReady() const {
    return Port::isReady();
}

void TextPort::clear() {
    if (!isOutport())
        LERROR("clear() called on inport");
    else
        portData_ = "";
}

bool TextPort::hasData() const {
    return true;
}

tgt::col3 TextPort::getColorHint() const {
    return tgt::col3(0, 255, 255);
}

bool TextPort::isDataThreadSafe() const {
    return true;
}

} // namespace
