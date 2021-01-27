/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "fieldplotdataport.h"

namespace voreen {

FieldPlotDataPort::FieldPlotDataPort(PortDirection direction, const std::string& id, const std::string& guiName, bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel)
    : GenericPort<FieldPlotData>(direction, id, guiName, allowMultipleConnections, invalidationLevel) {
}

std::string FieldPlotDataPort::getClassName() const {
    return "FieldPlotDataPort";
}

Port* FieldPlotDataPort::create(PortDirection direction, const std::string& id, const std::string& guiName) const {
    return new FieldPlotDataPort(direction,id,guiName);
}

tgt::col3 FieldPlotDataPort::getColorHint() const {
    return tgt::col3(150, 72, 4);
}

std::string FieldPlotDataPort::getContentDescription() const {
    std::stringstream strstr;
    strstr << Port::getContentDescription();
    if(hasData()) {
        strstr << std::endl << "Plot width/height: " << getData()->getWidth() << "/" << getData()->getHeight();
        strstr << std::endl << "Plot representation dim." << getData()->getVolume()->getDimensions().x << "/" << getData()->getVolume()->getDimensions().y << "/" << getData()->getVolume()->getDimensions().z;
    }
    return strstr.str();
}

std::string FieldPlotDataPort::getContentDescriptionHTML() const {
    std::stringstream strstr;
    strstr << Port::getContentDescriptionHTML();
    if(hasData()) {
        strstr << "<br>" << "Plot width/height: " << getData()->getWidth() << "/" << getData()->getHeight();
        strstr << "<br>" << "Plot representation dim." << getData()->getVolume()->getDimensions().x << "/" << getData()->getVolume()->getDimensions().y << "/" << getData()->getVolume()->getDimensions().z;
    }
    return strstr.str();
}



} // namespace
