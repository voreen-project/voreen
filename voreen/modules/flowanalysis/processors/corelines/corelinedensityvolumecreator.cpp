#include <chrono>

#include "corelinedensityvolumecreator.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

namespace voreen
{
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
		const auto timeBegin = std::chrono::high_resolution_clock::now();

#pragma omp parallel for
		for( long i = 0; i < static_cast<long>( outVolume.getNumVoxels() ); ++i )
		{
			outVolume.voxel( i ) = false;
		}


		for( long i = 0; i < static_cast<long>( corelines.size() ); ++i )
		{
			for( const auto point : corelines[i] )
			{
				outVolume.voxel( std::lround( point.x ), std::lround( point.y ), std::lround( point.z ) ) += 1.0f;
			}
		}

		const auto timeEnd = std::chrono::high_resolution_clock::now();
		const auto time = std::chrono::duration_cast<std::chrono::milliseconds>( timeEnd - timeBegin ).count();
		std::cout << "[CorelineDensityVolumeCreator]: Finished in " << time << " ms." << std::endl;
	}

	void CorelineDensityVolumeCreator::process()
	{
		if( !_inCorelines.hasData() || !_inVolume.hasData() )
			return;

		auto corelines = std::vector<std::vector<tgt::vec3>>();
		for( const auto geometry : _inCorelines.getAllData() )
		{
			const auto& current = dynamic_cast<const PointSegmentListGeometryVec3*>( geometry )->getData();
			corelines.insert( corelines.end(), current.begin(), current.end() );
		}

		const auto dim = _inVolume.getData()->getDimensions();
		auto binVol = std::unique_ptr<VolumeRAM_Float>( new VolumeRAM_Float( dim ) );

		CorelineDensityVolumeCreator::Process( corelines, *binVol );

		_out.setData( new Volume( binVol.release(), _inVolume.getData()->getSpacing(), _inVolume.getData()->getOffset() ) );
	}
} // namespace voreen
