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

#ifndef VRN_VORTEX_H
#define VRN_VORTEX_H

#include "tgt/vector.h"
#include <vector>

namespace voreen {

class Vortex {
public:
    enum class Orientation : int32_t
    {
        eUnknown = 0x0,
        eClockwise = 0x1,
        eCounterClockwise = 0x2
    };

    Vortex() = default;
    Vortex( Orientation orientation, std::vector<tgt::vec3> coreline );
    Vortex( std::istream& stream );

    void serialize( std::ostream& stream ) const;

    void setOrientation( Orientation orientation ) noexcept;
    Orientation getOrientation() const noexcept;

    void setCoreline( std::vector<tgt::vec3> coreline ) noexcept;
    const std::vector<tgt::vec3>& coreline() const noexcept;

private:
    Orientation _orientation;
    std::vector<tgt::vec3> _coreline;
};

class VortexCollection {
public:
    struct VortexID {
        size_t member, timestep, index;

        VortexID() noexcept = default;
        VortexID( size_t member, size_t timestep, size_t index );

        bool operator==( const VortexID& other ) const noexcept;
        bool operator!=( const VortexID& other ) const noexcept;

        static VortexID Invalid;
    };

    VortexCollection() = default;
    VortexCollection( size_t members, size_t timesteps );
    VortexCollection( std::istream& stream );

    void serialize( std::ostream& stream ) const;

    size_t members() const noexcept;
    size_t timesteps() const noexcept;
    size_t totalNumVortices() const;

    const std::vector<Vortex>& vortices( size_t members, size_t timestep ) const;
    void setVortices( size_t members, size_t timestep, std::vector<Vortex> vortices );

    const std::vector<VortexID>& matches( size_t members, size_t timestep, size_t index ) const;
    const std::vector<VortexID>& matches( VortexID vortexID ) const;

    void addMatch( VortexID first, VortexID second );

private:
    size_t _members, _timesteps;
    std::vector<std::vector<Vortex>> _vortices;
    std::vector<std::vector<std::vector<VortexID>>> _matches;
};

std::string to_string( Vortex::Orientation orientation );
}

#endif // VRN_VORTEX_H