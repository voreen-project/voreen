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

#include "vortex.h"

namespace voreen {

// -------------- //
// --- Vortex --- //
// -------------- //

Vortex::Vortex( Orientation orientation, std::vector<tgt::vec3> coreline ) : _orientation( orientation ), _coreline( std::move( coreline ) )
{}
Vortex::Vortex( std::istream& stream )
{
    stream.read( reinterpret_cast<char*>( &_orientation ), sizeof( int32_t ) );

    size_t corelineSize;
    stream.read( reinterpret_cast<char*>( &corelineSize ), sizeof( size_t ) );

    _coreline.resize( corelineSize );
    stream.read( reinterpret_cast<char*>( _coreline.data() ), corelineSize * sizeof( tgt::vec3 ) );
}

void Vortex::serialize( std::ostream& stream ) const
{
    const auto corelineSize = _coreline.size();
    stream.write( reinterpret_cast<const char*>( &_orientation ), sizeof( int32_t ) );
    stream.write( reinterpret_cast<const char*>( &corelineSize ), sizeof( size_t ) );
    stream.write( reinterpret_cast<const char*>( _coreline.data() ), corelineSize * sizeof( tgt::vec3 ) );
}


void Vortex::setOrientation( Orientation orientation ) noexcept
{
    _orientation = orientation;
}
Vortex::Orientation Vortex::getOrientation() const noexcept
{
    return _orientation;
}

void Vortex::setCoreline( std::vector<tgt::vec3> coreline ) noexcept
{
    _coreline = std::move( coreline );
}
const std::vector<tgt::vec3>& Vortex::coreline() const noexcept
{
    return _coreline;
}

// ------------------------- //
// --- Vortex Collection --- //
// ------------------------- //

VortexCollection::VortexID::VortexID( size_t run, size_t timestep, size_t index ) : run( run ), timestep( timestep ), index( index )
{}

bool VortexCollection::VortexID::operator==( const VortexCollection::VortexID& other ) const noexcept
{
    return run == other.run && timestep == other.timestep && index == other.index;
}
bool VortexCollection::VortexID::operator!=( const VortexCollection::VortexID& other ) const noexcept
{
    return run != other.run || timestep != other.timestep || index != other.index;
}

VortexCollection::VortexID VortexCollection::VortexID::Invalid = VortexCollection::VortexID( ~0u, ~0u, ~0u );



VortexCollection::VortexCollection( size_t runs, size_t timesteps ) : _runs( runs ), _timesteps( timesteps ), _vortices( runs* timesteps ), _matches( runs* timesteps )
{}
VortexCollection::VortexCollection( std::istream& stream )
{
    stream.read( reinterpret_cast<char*>( &_runs ), sizeof( size_t ) );
    stream.read( reinterpret_cast<char*>( &_timesteps ), sizeof( size_t ) );

    _vortices.resize( _runs * _timesteps );
    _matches.resize( _runs * _timesteps );
    for( size_t i = 0; i < _vortices.size(); ++i )
    {
        size_t vortexCount;
        stream.read( reinterpret_cast<char*>( &vortexCount ), sizeof( size_t ) );

        _vortices[i].resize( vortexCount );
        _matches[i].resize( vortexCount );
        for( size_t j = 0; j < vortexCount; ++j )
        {
            _vortices[i][j] = Vortex( stream );

            size_t matchCount;
            stream.read( reinterpret_cast<char*>( &matchCount ), sizeof( size_t ) );

            _matches[i][j].resize( matchCount );
            for( size_t k = 0; k < matchCount; ++k )
            {
                stream.read( reinterpret_cast<char*>( &_matches[i][j][k].run ), sizeof( size_t ) );
                stream.read( reinterpret_cast<char*>( &_matches[i][j][k].timestep ), sizeof( size_t ) );
                stream.read( reinterpret_cast<char*>( &_matches[i][j][k].index ), sizeof( size_t ) );
            }
        }
    }
}

void VortexCollection::serialize( std::ostream& stream ) const
{
    stream.write( reinterpret_cast<const char*>( &_runs ), sizeof( size_t ) );
    stream.write( reinterpret_cast<const char*>( &_timesteps ), sizeof( size_t ) );

    for( size_t i = 0; i < _vortices.size(); ++i )
    {
        const auto vortexCount = _vortices[i].size();
        stream.write( reinterpret_cast<const char*>( &vortexCount ), sizeof( size_t ) );

        for( size_t j = 0; j < _vortices[i].size(); ++j )
        {
            _vortices[i][j].serialize( stream );

            const auto matchCount = _matches[i][j].size();
            stream.write( reinterpret_cast<const char*>( &matchCount ), sizeof( size_t ) );
            for( size_t k = 0; k < matchCount; ++k )
            {
                stream.write( reinterpret_cast<const char*>( &_matches[i][j][k].run ), sizeof( size_t ) );
                stream.write( reinterpret_cast<const char*>( &_matches[i][j][k].timestep ), sizeof( size_t ) );
                stream.write( reinterpret_cast<const char*>( &_matches[i][j][k].index ), sizeof( size_t ) );
            }
        }
    }
}

size_t VortexCollection::runs() const noexcept
{
    return _runs;
}
size_t VortexCollection::timesteps() const noexcept
{
    return _timesteps;
}
size_t VortexCollection::totalNumVortices() const
{
    size_t count = 0;
    for( const auto& vortices : _vortices )
        count += vortices.size();
    return count;
}

const std::vector<Vortex>& VortexCollection::vortices( size_t run, size_t timestep ) const
{
    return _vortices[run * _timesteps + timestep];
}
void VortexCollection::setVortices( size_t run, size_t timestep, std::vector<Vortex> vortices )
{
    const auto index = run * _timesteps + timestep;
    _matches[index].resize( vortices.size() );
    _vortices[index] = std::move( vortices );
}

const std::vector<VortexCollection::VortexID>& VortexCollection::matches( size_t run, size_t timestep, size_t index ) const
{
    return _matches[run * _timesteps + timestep][index];
}
const std::vector<VortexCollection::VortexID>& VortexCollection::matches( VortexID vortexID ) const
{
    return this->matches( vortexID.run, vortexID.timestep, vortexID.index );
}

void VortexCollection::addMatch( VortexCollection::VortexID first, VortexCollection::VortexID second )
{
    _matches[first.run * _timesteps + first.timestep][first.index].push_back( second );
    _matches[second.run * _timesteps + second.timestep][second.index].push_back( first );
}

// -------------- //
// --- Global --- //
// -------------- //

std::string to_string( Vortex::Orientation orientation ) {
    switch( orientation )
    {
    case Vortex::Orientation::eUnknown:
        return "unknown";
    case Vortex::Orientation::eClockwise:
        return "clockwise";
    case Vortex::Orientation::eCounterClockwise:
        return "counter-clockwise";
    default:
        return "Vortex::Orientation";
    }
}

}
