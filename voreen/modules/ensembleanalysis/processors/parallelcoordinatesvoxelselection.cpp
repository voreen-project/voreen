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

#include "parallelcoordinatesvoxelselection.h"

namespace voreen {
ParallelCoordinatesVoxelSelection::ParallelCoordinatesVoxelSelection()
    : Processor(),
      ensembleport_(Port::INPORT, "port_ensemble", "Ensemble Input" ),
      volumeport_(Port::OUTPORT, "port_volume", "Binary Volume Output" ),
      propertySections_("property_sections", "Sections", Processor::INVALID_RESULT )
{
    // --- Initialize Ports --- //
    this->addPort(ensembleport_ );
    this->addPort(volumeport_ );

    // --- Initialize Properties --- //
    this->addProperty(propertySections_ );
    propertySections_.setVisibleFlag(false );
}

void ParallelCoordinatesVoxelSelection::process()
{
    const auto ensemble = ensembleport_.getData();

    // --- Find Member --- //
    const EnsembleMember* member = nullptr;
    for( const auto& m : ensemble->getMembers() )
    {
        if(m.getName() == propertySections_.get().member )
        {
            member = &m;
            break;
        }
    }
    if( !member || member->getTimeSteps().size() <= propertySections_.get().timestep )
    {
        volumeport_.clear();
        return;
    }

    // --- Collect Volumes --- //
    auto volumes = std::vector<const VolumeRAM*>();
    for( const auto& field : propertySections_.get().fields )
    {
        const VolumeBase* volume = member->getTimeSteps()[propertySections_.get().timestep].getVolume(field);
        if(volume)
            volumes.push_back( volume->getRepresentation<VolumeRAM>() );
    }
    if(volumes.size() != propertySections_.get().fields.size() )
    {
        volumeport_.clear();
        return;
    }

    // --- Fill Output --- //
    auto volume = std::unique_ptr<VolumeRAM_UInt8>( new VolumeRAM_UInt8( volumes.front()->getDimensions() ) );

#pragma omp parallel for
    for( long i = 0; i < static_cast<long>( volume->getNumVoxels() ); ++i )
    {
        volume->voxel( i ) = std::numeric_limits<uint8_t>::max();

        for( size_t j = 0; j < volumes.size(); ++j )
        {
            const auto& sections = propertySections_.get().sections[j];
            const auto value = volumes[j]->getVoxelNormalized( i );

            auto selected = sections.empty();
            for( const auto& section : sections )
            {
                if( value >= section.first && value <= section.second )
                {
                    selected = true;
                    break;
                }
            }

            if( !selected )
            {
                volume->voxel( i ) = 0;
                break;
            }
        }
    }

    volumeport_.setData(new Volume(volume.release(), tgt::vec3::one, tgt::vec3::zero ) );
}

}