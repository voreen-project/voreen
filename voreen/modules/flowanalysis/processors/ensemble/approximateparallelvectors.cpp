#include "approximateparallelvectors.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"

#include "modules/flowanalysis/processors/corelines/parallelvectors.h"
#include "modules/flowanalysis/processors/volume/accelerationprocessor.h"
#include "modules/flowanalysis/processors/volume/vortexprocessor.h"

#include "Eigen/Eigenvalues"
#include <chrono>

namespace voreen
{
	ApproximateParallelVectors::ApproximateParallelVectors() : Processor(),
		_inportEnsemble( Port::INPORT, "inport_ensemble", "Ensemble" ),
		_inportMask( Port::INPORT, "inport_mask", "Mask" ),
		_outportA( Port::OUTPORT, "outport_a", "Volume A" ),
		_outportB( Port::OUTPORT, "outport_b", "Volume B" ),
		_outportEps( Port::OUTPORT, "outport_eps", "Volume Epsilon" ),

		_propertyTimestep( "property_timestep", "Timestep", 0, 0, std::numeric_limits<int>::max(), Processor::VALID ),
		_propertyEpsilon( "property_epsilon", "Epsilon Threshold", -1.0f, -1.0f, 1.0f, Processor::VALID ),
		_propertyUseAcceleration( "property_use_aceleration", "Use Acceleration", false, Processor::VALID ),
		_propertyUpdateButton( "property_update_button", "Update", Processor::VALID )
	{
		this->addPort( _inportEnsemble );
		this->addPort( _inportMask );
		this->addPort( _outportA );
		this->addPort( _outportB );
		this->addPort( _outportEps );

		this->addProperty( _propertyTimestep );
		this->addProperty( _propertyEpsilon );
		this->addProperty( _propertyUseAcceleration );
		this->addProperty( _propertyUpdateButton );

		_inportEnsemble.onNewData( LambdaFunctionCallback( [this]
		{
			_propertyTimestep.setMaxValue( std::max( 0, static_cast<int>( _inportEnsemble.getData()->getMinNumTimeSteps() ) - 1 ) );
		} ) );
		_propertyUpdateButton.onClick( MemberFunctionCallback<ApproximateParallelVectors>( this, &ApproximateParallelVectors::updateButton ) );
	}

	Processor* ApproximateParallelVectors::create() const
	{
		return new ApproximateParallelVectors();
	}
	std::string ApproximateParallelVectors::getClassName() const
	{
		return "ApproximateParallelVectors";
	}
	std::string ApproximateParallelVectors::getCategory() const
	{
		return "Vortex Processing";
	}

	void ApproximateParallelVectors::process()
	{}
	void ApproximateParallelVectors::updateButton()
	{
		const auto ensemble = _inportEnsemble.getData();
		const auto maskVolumeBase = _inportMask.getData();
		const auto timestep = _propertyTimestep.get();

		if( !ensemble || !maskVolumeBase )
		{
			_outportA.clear();
			_outportB.clear();
			_outportEps.clear();
			return;
		}

		// --- Gather Voxels --- //
		const auto volumeMask = VolumeRAMRepresentationLock( maskVolumeBase );
		const auto dim = volumeMask->getDimensions();
		const auto& runs = ensemble->getMembers();

		auto voxels = std::vector<size_t>();
		voxels.reserve( volumeMask->getNumVoxels() );
		for( size_t i = 0; i < volumeMask->getNumVoxels(); ++i )
			if( volumeMask->getVoxelNormalized( i ) > 0.0 )
				voxels.push_back( i );
		voxels.shrink_to_fit();

		const auto getVelocity = [timestep, dim, &runs, &voxels] ( size_t index )
		{
			const auto volumeLockU = VolumeRAMRepresentationLock( runs[index].getTimeSteps()[timestep].getVolume( "U" ) );
			const auto volumeLockV = VolumeRAMRepresentationLock( runs[index].getTimeSteps()[timestep].getVolume( "V" ) );
			const auto volumeLockW = VolumeRAMRepresentationLock( runs[index].getTimeSteps()[timestep].getVolume( "W" ) );

			const auto volumeU = dynamic_cast<const VolumeRAM_Double*>( volumeLockU.operator->() );
			const auto volumeV = dynamic_cast<const VolumeRAM_Double*>( volumeLockV.operator->() );
			const auto volumeW = dynamic_cast<const VolumeRAM_Double*>( volumeLockW.operator->() );

			auto volumeVelocity = std::unique_ptr<VolumeRAM_3xDouble>( new VolumeRAM_3xDouble( dim ) );
			volumeVelocity->fill( tgt::dvec3::zero );
#pragma omp parallel for
			for( long i = 0; i < static_cast<long>( voxels.size() ); ++i )
			{
				const auto voxelIndex = voxels[i];
				volumeVelocity->voxel( voxelIndex ) = tgt::dvec3( volumeU->voxel( voxelIndex ), volumeV->voxel( voxelIndex ), volumeW->voxel( voxelIndex ) );
			}

			return std::unique_ptr<Volume>( new Volume( volumeVelocity.release(), tgt::vec3::one, tgt::vec3::zero ) );
		};

		// --- Create Volumes --- //
		auto volumeA = std::unique_ptr<VolumeRAM_3xDouble>( new VolumeRAM_3xDouble( dim ) );
		auto volumeC = std::unique_ptr<VolumeRAM_Mat3Float>( new VolumeRAM_Mat3Float( dim ) );
		auto volumeB = std::unique_ptr<VolumeRAM_3xDouble>( new VolumeRAM_3xDouble( dim ) );
		auto volumeEps = std::unique_ptr<VolumeRAM_Double>( new VolumeRAM_Double( dim ) );

		volumeA->fill( tgt::dvec3::zero );
		volumeC->fill( tgt::dmat3::zero );
		volumeB->fill( tgt::dvec3::zero );
		volumeEps->fill( 0.0 );

		// --- Calculate A --- //
		for( size_t i = 0; i < runs.size(); ++i )
		{
			const auto volumeBaseVelocity = getVelocity( i );
			const auto volumeBaseVelocityConv = std::unique_ptr<Volume>( VolumeOperatorConvert().apply<tgt::dvec3>( volumeBaseVelocity.get() ) );

			const auto volumeVelocity = dynamic_cast<const VolumeRAM_3xDouble*>( volumeBaseVelocity->getRepresentation<VolumeRAM>() );
			const auto volumeVelocityConv = dynamic_cast<const VolumeRAM_3xDouble*>( volumeBaseVelocityConv->getRepresentation<VolumeRAM>() );

#pragma omp parallel for
			for( long i = 0; i < static_cast<long>( voxels.size() ); ++i )
			{
				const auto voxelIndex = voxels[i];
				volumeA->voxel( voxelIndex ) += volumeVelocityConv->voxel( voxelIndex ) / static_cast<double>( runs.size() );
			}

			if( _propertyUseAcceleration.get() )
			{
				auto volumeJacobi = std::unique_ptr<VolumeRAM_Mat3Float>( new VolumeRAM_Mat3Float( dim ) );
				VortexProcessor::Process( *volumeVelocity, *volumeMask.operator->(), *volumeJacobi, nullptr, nullptr, nullptr );

				auto volumeAcceleration = std::unique_ptr<VolumeRAM_3xDouble>( new VolumeRAM_3xDouble( dim ) );
				AccelerationProcessor::Process( *volumeJacobi, *volumeVelocity, *volumeAcceleration );

#pragma omp parallel for
				for( long i = 0; i < static_cast<long>( voxels.size() ); ++i )
				{
					const auto voxelIndex = voxels[i];
					volumeA->voxel( voxelIndex ) += volumeAcceleration->voxel( voxelIndex ) / static_cast<double>( runs.size() );
				}
			}
		}

		// --- Calculate C --- //
		for( size_t i = 0; i < runs.size(); ++i )
		{
			const auto volumeBaseVelocity = getVelocity( i );
			const auto volumeBaseVelocityConv = std::unique_ptr<Volume>( VolumeOperatorConvert().apply<tgt::dvec3>( volumeBaseVelocity.get() ) );

			const auto volumeVelocity = dynamic_cast<const VolumeRAM_3xDouble*>( volumeBaseVelocity->getRepresentation<VolumeRAM>() );
			const auto volumeVelocityConv = dynamic_cast<const VolumeRAM_3xDouble*>( volumeBaseVelocityConv->getRepresentation<VolumeRAM>() );

#pragma omp parallel for
			for( long i = 0; i < static_cast<long>( voxels.size() ); ++i )
			{
				const auto voxelIndex = voxels[i];
				const auto v = volumeVelocityConv->voxel( voxelIndex );
				const auto a = volumeA->voxel( voxelIndex );

				auto& covariance = volumeC->voxel( voxelIndex );
				for( int i = 0; i < 3; ++i ) for( int j = 0; j < 3; ++j )
					covariance[i][j] += ( v[i] - a[i] ) * ( v[j] - a[j] );
			}

			if( _propertyUseAcceleration.get() )
			{
				auto volumeJacobi = std::unique_ptr<VolumeRAM_Mat3Float>( new VolumeRAM_Mat3Float( dim ) );
				VortexProcessor::Process( *volumeVelocity, *volumeMask.operator->(), *volumeJacobi, nullptr, nullptr, nullptr );

				auto volumeAcceleration = std::unique_ptr<VolumeRAM_3xDouble>( new VolumeRAM_3xDouble( dim ) );
				AccelerationProcessor::Process( *volumeJacobi, *volumeVelocity, *volumeAcceleration );

#pragma omp parallel for
				for( long i = 0; i < static_cast<long>( voxels.size() ); ++i )
				{
					const auto voxelIndex = voxels[i];
					const auto acc = volumeAcceleration->voxel( voxelIndex );
					const auto a = volumeA->voxel( voxelIndex );

					auto& covariance = volumeC->voxel( voxelIndex );
					for( int i = 0; i < 3; ++i ) for( int j = 0; j < 3; ++j )
						covariance[i][j] += ( acc[i] - a[i] ) * ( acc[j] - a[j] );
				}
			}
		}

		// --- Calculate B --- //
		auto positiveEpsCount = 0;
		for( long i = 0; i < static_cast<long>( voxels.size() ); ++i )
		{
			const auto voxelIndex = voxels[i];
			const auto a = volumeA->voxel( voxelIndex );
			const auto c = volumeC->voxel( voxelIndex );
			const auto b = volumeB->voxel( voxelIndex ) = tgt::dvec3( c * a );

			auto eigenMatrix = Eigen::Matrix3d();
			for( int i = 0; i < 3; ++i ) for( int j = 0; j < 3; ++j )
				eigenMatrix( i, j ) = c[i][j];
			auto eigenvalues = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>( eigenMatrix ).eigenvalues();
			std::sort( eigenvalues.data(), eigenvalues.data() + 3 );

			const auto chi = std::abs( tgt::length( b ) / tgt::length( a ) - std::abs( eigenvalues.z() ) ) < 0.00001;
			volumeEps->voxel( voxelIndex ) = ( eigenvalues.z() - eigenvalues.y() ) / eigenvalues.sum() * ( chi ? 1.0 : -1.0 );
		}

		_outportA.setData( new Volume( volumeA.release(), tgt::vec3::one, tgt::vec3::zero ) );
		_outportB.setData( new Volume( volumeB.release(), tgt::vec3::one, tgt::vec3::zero ) );
		_outportEps.setData( new Volume( volumeEps.release(), tgt::vec3::one, tgt::vec3::zero ) );
	}
}