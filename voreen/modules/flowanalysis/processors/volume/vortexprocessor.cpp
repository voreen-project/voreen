#include "vortexprocessor.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "Eigen/Eigenvalues"

#include <chrono>

namespace voreen
{
	VortexProcessor::VortexProcessor() : Processor(),
		_inportVelocity( Port::INPORT, "inport_velociy", "Velocity" ),
		_inportMask( Port::INPORT, "inport_mask", "Mask" ),
		_outportJacobi( Port::OUTPORT, "outport_jacobi", "Jacobi" ),
		_outportDelta( Port::OUTPORT, "outport_delta", "Delta-Criterion" ),
		_outportLambda2( Port::OUTPORT, "outport_lambda2", "Lambda2-Criterion" ),
		_outportQ( Port::OUTPORT, "outport_q", "Q-Criterion" )
	{
		this->addPort( _inportVelocity );
		this->addPort( _inportMask );
		this->addPort( _outportJacobi );
		this->addPort( _outportDelta );
		this->addPort( _outportLambda2 );
		this->addPort( _outportQ );
	}

	void VortexProcessor::process()
	{
		const auto velocityVolume = _inportVelocity.getData();
		const auto maskVolume = _inportMask.getData();

		if( !velocityVolume || !maskVolume || velocityVolume->getDimensions() != maskVolume->getDimensions() )
		{
			_outportJacobi.clear();
			_outportDelta.clear();
			_outportLambda2.clear();
			_outportQ.clear();
			return;
		}

		const auto velocity = velocityVolume->getRepresentation<VolumeRAM>();
		const auto mask = maskVolume->getRepresentation<VolumeRAM>();
		const auto dim = velocity->getDimensions();

		auto volumeJacobi = std::unique_ptr<VolumeRAM_Mat3Float>( new VolumeRAM_Mat3Float( dim ) );
		auto volumeDelta = std::unique_ptr<VolumeRAM_Float>( new VolumeRAM_Float( dim ) );
		auto volumeLambda2 = std::unique_ptr<VolumeRAM_Float>( new VolumeRAM_Float( dim ) );
		auto volumeQ = std::unique_ptr<VolumeRAM_Float>( new VolumeRAM_Float( dim ) );

		VortexProcessor::Process( *velocity, *mask, *volumeJacobi, volumeDelta.get(), volumeLambda2.get(), volumeQ.get() );

		// --- Update Outports --- //
		auto volume = new Volume( volumeJacobi.release(), tgt::vec3( 1.0f, 1.0f, 1.0f ), tgt::vec3::zero );
		volume->getMetaDataContainer().addMetaData( "name", new StringMetaData( "JACOBI" ) );
		_outportJacobi.setData( volume );

		volume = new Volume( volumeDelta.release(), tgt::vec3( 1.0f, 1.0f, 1.0f ), tgt::vec3::zero );
		volume->getMetaDataContainer().addMetaData( "name", new StringMetaData( "DELTA" ) );
		_outportDelta.setData( volume );

		volume = new Volume( volumeLambda2.release(), tgt::vec3( 1.0f, 1.0f, 1.0f ), tgt::vec3::zero );
		volume->getMetaDataContainer().addMetaData( "name", new StringMetaData( "LAMBDA2" ) );
		_outportLambda2.setData( volume );

		volume = new Volume( volumeQ.release(), tgt::vec3( 1.0f, 1.0f, 1.0f ), tgt::vec3::zero );
		volume->getMetaDataContainer().addMetaData( "name", new StringMetaData( "Q" ) );
		_outportQ.setData( volume );
	}

	void VortexProcessor::Process( const VolumeRAM& velocity, const VolumeRAM& mask, VolumeRAM_Mat3Float& outJacobi, VolumeRAM_Float* outDelta, VolumeRAM_Float* outLambda2, VolumeRAM_Float* outQ )
	{
		const auto timeBegin = std::chrono::high_resolution_clock::now();

		const std::array<tgt::svec3, 7> offsets { tgt::svec3( 0, 0, 0 ),
			tgt::svec3( -1, 0, 0 ), tgt::svec3( 1, 0, 0 ),
			tgt::svec3( 0, -1, 0 ), tgt::svec3( 0, 1, 0 ),
			tgt::svec3( 0, 0, -1 ), tgt::svec3( 0, 0, 1 )
		};

		const auto dim = velocity.getDimensions();

		auto voxels = std::vector<tgt::ivec3>();
		voxels.reserve( mask.getNumVoxels() );
		for( auto x = 3; x < dim.x - 3; ++x ) for( auto y = 3; y < dim.y - 3; ++y ) for( auto z = 3; z < dim.z - 3; ++z )
		{
			auto addVoxel = true;
			for( const auto offset : offsets )
			{
				if( mask.getVoxelNormalized( tgt::svec3( x, y, z ) + offset ) <= 0.0f )
				{
					addVoxel = false;
					break;
				}
			}
			if( addVoxel ) voxels.push_back( tgt::ivec3( x, y, z ) );
		}
		voxels.shrink_to_fit();

#pragma omp parallel for
		for( long i = 0; i < static_cast<long>( voxels.size() ); ++i )
		{
			const auto position = static_cast<tgt::svec3>( voxels[i] );

			auto jacobi = Eigen::Matrix3f();
			for( int i = 0; i < 3; ++i )
			{
				for( int j = 0; j < 3; ++j )
				{
					const auto b3 = velocity.getVoxelNormalized( position + offsets[( j + 1 ) * 2] * size_t( 3 ), i );
					const auto b2 = velocity.getVoxelNormalized( position + offsets[( j + 1 ) * 2] * size_t( 2 ), i );
					const auto b1 = velocity.getVoxelNormalized( position + offsets[( j + 1 ) * 2], i );
					const auto c = velocity.getVoxelNormalized( position, i );
					const auto f1 = velocity.getVoxelNormalized( position + offsets[( j + 1 ) * 2 - 1], i );
					const auto f2 = velocity.getVoxelNormalized( position + offsets[( j + 1 ) * 2 - 1] * size_t( 2 ), i );
					const auto f3 = velocity.getVoxelNormalized( position + offsets[( j + 1 ) * 2 - 1] * size_t( 3 ), i );

					const auto db2 = ( b1 - b3 ) / 2.0f;
					const auto db1 = ( c - b2 ) / 2.0f;
					const auto dc = ( f1 - b1 ) / 2.0f;
					const auto df1 = ( f2 - c ) / 2.0f;
					const auto df2 = ( f3 - f1 ) / 2.0f;

					const auto ddb1 = ( dc - db2 ) / 2.0f;
					const auto ddf1 = ( df2 - dc ) / 2.0f;

					const auto t = 0.5f;
					//jacobi( i, j ) = ( 2.0f * t * t * t - 3.0f * t * t + 1.0f ) * db1
					//	+ ( -2.0f * t * t * t + 3.0f * t * t ) * df1
					//	+ ( t * t * t - 2.0f * t * t + t ) * ddb1
					//	+ ( t * t * t - t * t ) * ddf1;
					jacobi( i, j ) = 3.0f * ( 2.0f * b1 - 2.0f * f1 + db1 + df1 ) * t * t +
						2.0f * ( -3.0f * b1 + 3.0f * f1 - 2.0f * db1 - df1 ) * t
						+ db1;
					// jacobi( i, j ) = dc;
				}
			}
			for( int i = 0; i < 3; ++i ) for( int j = 0; j < 3; ++j )
				outJacobi.voxel( position )[i][j] = jacobi( i, j );

			// --- Compute Delta Criterion --- //
			auto eigenvalues = Eigen::Vector3cf();
			if( outDelta )
			{
				eigenvalues = jacobi.eigenvalues();
				outDelta->voxel( position ) = std::max( { eigenvalues.x().imag(), eigenvalues.y().imag(), eigenvalues.z().imag() } );
			}

			// --- Compute Strain Rate and Vorticity --- //
			const auto jacobiT = jacobi.transpose();
			const auto S = 0.5f * ( jacobi + jacobiT );
			const auto O = 0.5f * ( jacobi - jacobiT );

			// --- Compute Lambda2 Criterion --- //
			if( outLambda2 )
			{
				eigenvalues = ( S * S + O * O ).eigenvalues();
				if( eigenvalues.x().real() < eigenvalues.y().real() )
					if( eigenvalues.y().real() < eigenvalues.z().real() ) outLambda2->voxel( position ) = eigenvalues.y().real();
					else if( eigenvalues.x().real() < eigenvalues.z().real() ) outLambda2->voxel( position ) = eigenvalues.z().real();
					else outLambda2->voxel( position ) = eigenvalues.x().real();
				else
					if( eigenvalues.x().real() < eigenvalues.z().real() ) outLambda2->voxel( position ) = eigenvalues.x().real();
					else if( eigenvalues.y().real() < eigenvalues.z().real() ) outLambda2->voxel( position ) = eigenvalues.z().real();
					else outLambda2->voxel( position ) = eigenvalues.y().real();
			}

			// --- Compute Q Criterion --- //
			if( outQ )
			{
				eigenvalues = ( S.transpose() * S ).eigenvalues();
				const auto spectralNormS = std::max( { eigenvalues.x().real(), eigenvalues.y().real(), eigenvalues.z().real() } );
				eigenvalues = ( O.transpose() * O ).eigenvalues();
				const auto spectralNormO = std::max( { eigenvalues.x().real(), eigenvalues.y().real(), eigenvalues.z().real() } );
				outQ->voxel( position ) = 0.5f * ( spectralNormO - spectralNormS );
			}
		}

		const auto timeEnd = std::chrono::high_resolution_clock::now();
		const auto time = std::chrono::duration_cast<std::chrono::milliseconds>( timeEnd - timeBegin ).count();
		// std::cout << "[VortexProcessor]: Finished in " << time / 1000.0f << " seconds." << std::endl;
	}
}