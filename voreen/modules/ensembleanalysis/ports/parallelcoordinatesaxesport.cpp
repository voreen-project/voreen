/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "parallelcoordinatesaxesport.h"
#include <sstream>

namespace voreen {

ParallelCoordinatesAxesPort::ParallelCoordinatesAxesPort( PortDirection direction, const std::string& id, const std::string& guiName )
    : GenericPort<ParallelCoordinatesAxes>( direction, id, guiName, false, direction == INPORT ? Processor::INVALID_RESULT : Processor::VALID )
{}

Port* ParallelCoordinatesAxesPort::create( PortDirection direction, const std::string& id, const std::string& guiName ) const {
    return new ParallelCoordinatesAxesPort( direction, id, guiName );
}

std::string ParallelCoordinatesAxesPort::getClassName() const {
    return "ParallelCoordinatesAxesPort";
}

std::string ParallelCoordinatesAxesPort::getContentDescription() const {
    auto stream = std::stringstream();
    stream << Port::getContentDescriptionHTML();

    if( this->hasData() )
    {
        const auto axes = this->getData();
        stream << std::endl << "Runs: " << axes->runs();
        stream << std::endl << "Timesteps: " << axes->timesteps();
        stream << std::endl << "Fields: " << axes->fields();
        stream << std::endl << "Samples: " << axes->samples();
        stream << std::endl << "Memory: " << axes->memorySize() / 1000000.0f << " MB";
    }

    return stream.str();
}

std::string ParallelCoordinatesAxesPort::getContentDescriptionHTML() const {
    auto stream = std::stringstream();
    stream << Port::getContentDescriptionHTML();

    if( this->hasData() )
    {
        const auto axes = this->getData();
        stream << "<br/>" << "Runs: " << axes->runs();
        stream << "<br/>" << "Timesteps: " << axes->timesteps();
        stream << "<br/>" << "Fields: " << axes->fields();
        stream << "<br/>" << "Samples: " << axes->samples();
        stream << "<br/>" << "Memory: " << axes->memorySize() / 1000000.0f << " MB";
    }

    return stream.str();
}

}