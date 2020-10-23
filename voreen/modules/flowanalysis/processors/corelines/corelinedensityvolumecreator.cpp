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

#include "corelinedensityvolumecreator.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

namespace voreen {

CorelineDensityVolumeCreator::CorelineDensityVolumeCreator() : Processor(),
    _inCorelines( Port::INPORT, "inportCorelines", "Corelines", true ),
    _inVolume( Port::INPORT, "inportVolume", "One of the volumes used for parallel vectors and coreline creator. Only used for the output dimensions." ),
    _out( Port::OUTPORT, "outport", "Binary volume of coreline presence in voxel" )
{
    this->addPort( _inCorelines );
    this->addPort( _inVolume );
    this->addPort( _out );
}

void CorelineDensityVolumeCreator::Process( const std::vector<std::vector<tgt::vec3>>& corelines, VolumeRAM_Float& outVolume )
{
    outVolume.clear();

    for( long i = 0; i < static_cast<long>( corelines.size() ); ++i ) {
        for( const auto point : corelines[i] ) {
            outVolume.voxel( std::lround( point.x ), std::lround( point.y ), std::lround( point.z ) ) += 1.0f;
        }
    }
}

void CorelineDensityVolumeCreator::process() {

    auto corelines = std::vector<std::vector<tgt::vec3>>();
    for( const auto* geometry : _inCorelines.getAllData() )
    {
        const auto* current = dynamic_cast<const PointSegmentListGeometryVec3*>( geometry );
        if(current) {
            const auto& data = current->getData();
            corelines.insert(corelines.end(), data.begin(), data.end());
        }
    }

    const auto dim = _inVolume.getData()->getDimensions();
    auto binVol = std::unique_ptr<VolumeRAM_Float>( new VolumeRAM_Float( dim ) );

    CorelineDensityVolumeCreator::Process( corelines, *binVol );

    _out.setData( new Volume( binVol.release(), _inVolume.getData()->getSpacing(), _inVolume.getData()->getOffset() ) );
}

} // namespace voreen
