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

#include "vortexport.h"
#include <sstream>

namespace voreen {
// ------------------- //
// --- Vortex Port --- //
// ------------------- //

VortexPort::VortexPort( PortDirection direction, const std::string& id, const std::string& guiName, bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel )
    : GenericPort<Vortex>( direction, id, guiName, allowMultipleConnections, invalidationLevel )
{}

Port* VortexPort::create( PortDirection direction, const std::string& id, const std::string& guiName ) const
{
    return new VortexPort( direction, id, guiName );
}
std::string VortexPort::getClassName() const
{
    return "VortexPort";
}
std::string VortexPort::getContentDescriptionHTML() const
{
    auto stream = std::stringstream();
    stream << Port::getContentDescriptionHTML();

    if( this->hasData() )
    {
        const auto vortex = this->getData();
        stream << "<br/>Orientation: " << to_string( vortex->getOrientation() );
    }

    return stream.str();
}

// ------------------------------ //
// --- Vortex Collection Port --- //
// ------------------------------ //

VortexCollectionPort::VortexCollectionPort( PortDirection direction, const std::string& id, const std::string& guiName, bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel )
    : GenericPort<VortexCollection>( direction, id, guiName, allowMultipleConnections, invalidationLevel )
{}

Port* VortexCollectionPort::create( PortDirection direction, const std::string& id, const std::string& guiName ) const
{
    return new VortexPort( direction, id, guiName );
}
std::string VortexCollectionPort::getClassName() const
{
    return "VortexCollectionPort";
}
std::string VortexCollectionPort::getContentDescriptionHTML() const
{
    auto stream = std::stringstream();
    stream << Port::getContentDescriptionHTML();

    if( this->hasData() )
    {
        const auto collection = this->getData();
        stream << "<br/>Runs: " << std::to_string( collection->runs() );
        stream << "<br/>Timesteps: " << std::to_string( collection->timesteps() );
        stream << "<br/>Vortices: " << std::to_string( collection->totalNumVortices() );
    }

    return stream.str();
}

// ------------------------ //
// --- Vortex List Port --- //
// ------------------------ //

VortexListPort::VortexListPort( PortDirection direction, const std::string& id, const std::string& guiName, bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel )
    : GenericPort<std::vector<Vortex>>( direction, id, guiName, allowMultipleConnections, invalidationLevel )
{}

Port* VortexListPort::create( PortDirection direction, const std::string& id, const std::string& guiName ) const
{
    return new VortexListPort( direction, id, guiName );
}
std::string VortexListPort::getClassName() const
{
    return "VortexListPort";
}
std::string VortexListPort::getContentDescriptionHTML() const
{
    auto stream = std::stringstream();
    stream << Port::getContentDescriptionHTML();

    if( this->hasData() )
    {
        const auto vortices = this->getData();
        stream << "<br/>Number of vortices: " << vortices->size();
    }

    return stream.str();
}

}
