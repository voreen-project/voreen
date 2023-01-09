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

#include "corelinelengthvolumecreator.h"
#include "voreen/core/datastructures/geometry/geometrysequence.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

#include <tuple>
#include <set>

namespace voreen {

CorelineLengthVolumeCreator::CorelineLengthVolumeCreator() : Processor(),
    _inCorelines( Port::INPORT, "inportCorelines", "Corelines", true ),
    _inVolume( Port::INPORT, "inportVolume", "One of the volumes used for parallel vectors and coreline creator. Only used for the output dimensions." ),
    _out( Port::OUTPORT, "outport", "Binary volume of coreline presence in voxel" )
{
    this->addPort( _inCorelines );
    this->addPort( _inVolume );
    this->addPort( _out );
}

void CorelineLengthVolumeCreator::Process( const std::vector<std::vector<tgt::vec3>>& corelines, VolumeRAM_Float& outVolume )
{
    outVolume.clear();

	auto countVol = std::unique_ptr<VolumeRAM_Float>( new VolumeRAM_Float( outVolume.getDimensions() ) );
	countVol->clear();

	std::set<std::tuple<std::size_t, std::size_t, std::size_t>> traversedVoxels {};

	for (const auto& coreline : corelines) {
		if (coreline.size() < 2) {
			continue;
		}

		float corelineLength = 0.0f;
		for (std::size_t i = 1; i < coreline.size(); ++i) {
			corelineLength += tgt::distance(coreline[i - 1], coreline[i]);
		}

		traversedVoxels.clear();

		for (const auto& point : coreline) {
			std::size_t x = std::lround(point.x);
			std::size_t y = std::lround(point.y);
			std::size_t z = std::lround(point.z);

			auto result = traversedVoxels.emplace(x, y, z);
			if (result.second) {    // inserted ?
				outVolume.voxel(x, y, z) += corelineLength;
				countVol->voxel(x, y, z) += 1.0f;
			}
		}
	}

	// Average voxels by the number of corelines passing through it.
	auto numVoxels = outVolume.getNumVoxels();
	for (std::size_t i = 0; i < numVoxels; ++i) {
		auto count = std::max(1.0f, countVol->voxel(i));
		outVolume.voxel(i) /= count;
	}
}

void CorelineLengthVolumeCreator::process() {

    auto corelines = std::vector<std::vector<tgt::vec3>>();
    for( const auto* geometry : _inCorelines.getAllData() )
    {
        const auto* current = dynamic_cast<const PointSegmentListGeometryVec3*>( geometry );
        if(current) {
            const auto& data = current->getData();
            corelines.insert(corelines.end(), data.begin(), data.end());
        } else {
			const auto* sequence = dynamic_cast<const GeometrySequence*>(geometry);
			if (!sequence) {
				continue;
			}

			auto numGeometries = sequence->getNumGeometries();
			for (auto i = 0u; i < numGeometries; ++i) {
				const auto* subGeometry = dynamic_cast<const PointSegmentListGeometryVec3*>(sequence->getGeometry(i));
				if (subGeometry) {
					const auto& data = subGeometry->getData();
					corelines.insert(corelines.end(), data.begin(), data.end());
				}
			}
		}
    }

    const auto dim = _inVolume.getData()->getDimensions();
    auto binVol = std::unique_ptr<VolumeRAM_Float>( new VolumeRAM_Float( dim ) );

	CorelineLengthVolumeCreator::Process( corelines, *binVol );

    _out.setData( new Volume( binVol.release(), _inVolume.getData()->getSpacing(), _inVolume.getData()->getOffset() ) );
}

} // namespace voreen
