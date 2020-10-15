#pragma once

#ifndef VRN_FLOWMAPPROCESSOR_H
#define VRN_FLOWMAPPROCESSOR_H

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/ports/volumeport.h"
#include "Eigen/Eigenvalues"

namespace voreen
{
	class FlowMapProcessor : public Processor
	{
	public:
		FlowMapProcessor() : Processor(),
			_inportX( Port::INPORT, "inport_x", "Velocity Field List (x)" ),
			_inportY( Port::INPORT, "inport_y", "Velocity Field List (y)" ),
			_inportZ( Port::INPORT, "inport_z", "Velocity Field List (z)" ),
			_outportFlowmap( Port::OUTPORT, "outport_flowmap", "Flow Map" ),
			_outportFTLE( Port::OUTPORT, "outport_ftle", "FTLE Field" ),
			_propertyTimesteps( "property_timesteps", "Timesteps", 0, 0, std::numeric_limits<int>::max(), 0, std::numeric_limits<int>::max(), Processor::VALID ),
			_propertyStepsize( "property_stepsize", "Stepsize", 1.0f, 0.0f, 100.0f, Processor::VALID ),
			_propertyComputeFlowmap( "property_comute_flowmap", "Compute Flowmap", Processor::VALID ),
			_propertyComputeFTLEField( "property_comute_ftle_field", "Compute FTLE Field", Processor::VALID )
		{
			this->addPort( _inportX );
			this->addPort( _inportY );
			this->addPort( _inportZ );
			this->addPort( _outportFlowmap );
			this->addPort( _outportFTLE );

			this->addProperty( _propertyTimesteps );
			this->addProperty( _propertyStepsize );
			this->addProperty( _propertyComputeFlowmap );
			this->addProperty( _propertyComputeFTLEField );

			_inportX.onNewData( MemberFunctionCallback<FlowMapProcessor>( this, &FlowMapProcessor::onNewInportData ) );
			_inportY.onNewData( MemberFunctionCallback<FlowMapProcessor>( this, &FlowMapProcessor::onNewInportData ) );
			_inportZ.onNewData( MemberFunctionCallback<FlowMapProcessor>( this, &FlowMapProcessor::onNewInportData ) );

			_propertyComputeFlowmap.onClick( MemberFunctionCallback<FlowMapProcessor>( this, &FlowMapProcessor::computeFlowmap ) );
			_propertyComputeFTLEField.onClick( MemberFunctionCallback<FlowMapProcessor>( this, &FlowMapProcessor::computeFTLEField ) );
		}
		Processor* create() const override
		{
			return new FlowMapProcessor();
		}
		std::string getClassName() const override
		{
			return "FlowMapProcessor";
		}
		std::string getCategory() const override
		{
			return "Flow";
		}

		bool isReady() const override
		{
			return _inportX.isReady() && _inportY.isReady() && _inportY.isReady();
		}

	private:
		void process() override {}
		void onNewInportData()
		{
			const auto volumeListX = _inportX.getData();
			const auto volumeListY = _inportY.getData();
			const auto volumeListZ = _inportZ.getData();
			if( !volumeListX || !volumeListY || !volumeListZ ) return;

			const auto timesteps = std::min( volumeListX->size(), std::min( volumeListY->size(), volumeListZ->size() ) );
			_propertyTimesteps.setMaxValue( static_cast<int>( timesteps - 1 ) );
			_propertyStepsize.setMaxValue( timesteps - 2 );
		}
		void computeFlowmap()
		{
			const auto volumeListX = _inportX.getData();
			const auto volumeListY = _inportY.getData();
			const auto volumeListZ = _inportZ.getData();

			if( !volumeListX || !volumeListY || !volumeListZ ||
				!( volumeListX->size() == volumeListY->size() && volumeListY->size() == volumeListZ->size() ) )
			{
				_outportFlowmap.clear();
				_outportFTLE.clear();
				return;
			}

			const auto dim = volumeListX->first()->getDimensions();
			const auto stepsize = static_cast<float>( _propertyStepsize.get() );
			const auto interval = _propertyTimesteps.get();

			std::cout << "[FlowMapProcessor]: Flow Map Computation --> Interval = " << interval << ", Stepsize = " << stepsize << std::endl;

			auto volumes = std::vector<std::array<const VolumeRAM_Double*, 3>>( interval.y - interval.x + 1, std::array<const VolumeRAM_Double*, 3> {} );
			const auto velocity = [dim, &volumes] ( tgt::dvec3 p, double t )
			{
				const auto sample = [dim, &volumes] ( tgt::vec3 p, size_t t )
				{
					if( p.x >= 0.0 && p.x < dim.x && p.y >= 0.0 && p.y < dim.y && p.z >= 0.0 && p.z < dim.z )
						return tgt::vec3( volumes[t][0]->voxel( p ), volumes[t][1]->voxel( p ), volumes[t][2]->voxel( p ) );
					return tgt::vec3( 0.0, 0.0, 0.0 );
				};
				const auto trilinearInterpolation = [&sample] ( tgt::vec3 p, size_t t )
				{
					const auto lx = std::floor( p.x ), hx = std::ceil( p.x );
					const auto ly = std::floor( p.y ), hy = std::ceil( p.y );
					const auto lz = std::floor( p.z ), hz = std::ceil( p.z );
					const auto wx = p.x - lx, wy = p.y - ly, wz = p.z - lz;

					const tgt::vec3 corners[8] {
						sample( tgt::vec3( lx, ly, lz ), t ), sample( tgt::vec3( lx, ly, hz ), t ),
						sample( tgt::vec3( lx, hy, lz ), t ), sample( tgt::vec3( lx, hy, hz ), t ),
						sample( tgt::vec3( hx, ly, lz ), t ), sample( tgt::vec3( hx, ly, hz ), t ),
						sample( tgt::vec3( hx, hy, lz ), t ), sample( tgt::vec3( hx, hy, hz ), t )
					};

					const tgt::vec3 xinterpolation[4] {
						( 1.0f - wx ) * corners[0] + wx * corners[4],
						( 1.0f - wx ) * corners[1] + wx * corners[5],
						( 1.0f - wx ) * corners[2] + wx * corners[6],
						( 1.0f - wx ) * corners[3] + wx * corners[7]
					};

					const tgt::vec3 yinterpolation[2] {
						( 1.0f - wy ) * xinterpolation[0] + wy * xinterpolation[2],
						( 1.0f - wy ) * xinterpolation[1] + wy * xinterpolation[3]
					};

					return ( 1.0f - wz ) * yinterpolation[0] + wz * yinterpolation[1];
				};

				const auto first = trilinearInterpolation( p, static_cast<size_t>( std::floor( t ) ) );
				const auto second = trilinearInterpolation( p, static_cast<size_t>( std::ceil( t ) ) );
				const auto w = static_cast<float>( t - std::floor( t ) );

				return ( 1.0f - w ) * first + w * second;
			};

			auto flowmap = std::unique_ptr<VolumeRAM_3xFloat>( new VolumeRAM_3xFloat( dim ) );
			flowmap->fill( tgt::dvec3( 0.0, 0.0, 0.0 ) );

			for( auto t = 0; t <= interval.y - interval.x - std::ceil( stepsize ); ++t )
			{
				std::cout << "[FlowMapProcessor]: Flow Map Computation --> Timestep = " << t << std::endl;
				for( size_t i = t; i <= t + std::ceil( stepsize ); ++i ) if( volumes[i][0] == nullptr )
				{
					volumes[i][0] = dynamic_cast<const VolumeRAM_Double*>( volumeListX->at( interval.x + i )->getRepresentation<VolumeRAM>() );
					volumes[i][1] = dynamic_cast<const VolumeRAM_Double*>( volumeListY->at( interval.x + i )->getRepresentation<VolumeRAM>() );
					volumes[i][2] = dynamic_cast<const VolumeRAM_Double*>( volumeListZ->at( interval.x + i )->getRepresentation<VolumeRAM>() );
				}

				// const auto volumeX = dynamic_cast<const VolumeRAM_Double*>( volumeListX->at( interval.x + t )->getRepresentation<VolumeRAM>() );
				// const auto volumeY = dynamic_cast<const VolumeRAM_Double*>( volumeListY->at( interval.x + t )->getRepresentation<VolumeRAM>() );
				// const auto volumeZ = dynamic_cast<const VolumeRAM_Double*>( volumeListZ->at( interval.x + t )->getRepresentation<VolumeRAM>() );

#pragma omp parallel for
				for( long x = 0; x < static_cast<long>( dim.x ); ++x ) for( long y = 0; y < dim.y; ++y ) for( long z = 0; z < dim.z; ++z )
				{
					auto& flow = flowmap->voxel( x, y, z );
					const auto p = tgt::vec3( x + 0.5f, y + 0.5f, z + 0.5f ) + flow;
					if( p.x < 0.0 || p.x >= dim.x || p.y < 0.0 || p.y >= dim.y || p.z < 0.0 || p.z >= dim.z ) continue;

					const auto k1 = velocity( p, t );
					const auto k2 = velocity( p + k1 * stepsize / 2.0f, t + stepsize / 2.0f );
					const auto k3 = velocity( p + k2 * stepsize / 2.0f, t + stepsize / 2.0f );
					const auto k4 = velocity( p + k3 * stepsize, t + stepsize );
					flow += stepsize * ( k1 + 2.0f * k2 + 2.0f * k3 + k4 ) / 6.0f;

					// flow += tgt::dvec3( volumeX->voxel( p ), volumeY->voxel( p ), volumeZ->voxel( p ) );
				}
			}

			auto output = new Volume( flowmap.release(), tgt::vec3( 1, 1, 1 ), tgt::vec3::zero );
			output->getMetaDataContainer().addMetaData( "name", new StringMetaData( "FlowMap" ) );
			_outportFlowmap.setData( output );
		}
		void computeFTLEField()
		{
			if( !_outportFlowmap.hasData() ) this->computeFlowmap();
			const auto flowmap = dynamic_cast<const VolumeRAM_3xFloat*>( _outportFlowmap.getData()->getRepresentation<VolumeRAM>() );

			const auto dim = flowmap->getDimensions();
			const auto interval = _propertyTimesteps.get();

			auto ftleField = std::unique_ptr<VolumeRAM_Float>( new VolumeRAM_Float( dim ) );
			auto totalMin = std::numeric_limits<float>::max(), totalMax = std::numeric_limits<float>::min();

#pragma omp parallel for
			for( long x = 0; x < static_cast<long>( dim.x ); ++x )
			{
				std::cout << "[FlowMapProcessor]: FTLE Field Computation -> x = " << x << std::endl;

				for( long y = 0; y < dim.y; ++y )
				{
					for( long z = 0; z < dim.z; ++z )
					{
						auto spacing = tgt::vec3( 2.0, 2.0, 2.0 );
						auto xl = x - 1, xh = x + 1;
						auto yl = y - 1, yh = y + 1;
						auto zl = z - 1, zh = z + 1;

						if( xl < 0 ) xl = 0, spacing.x = 1.0;
						if( xh >= dim.x ) xh = dim.x - 1, spacing.x = 1.0;
						if( yl < 0 ) yl = 0, spacing.y = 1.0;
						if( yh >= dim.y ) yh = dim.y - 1, spacing.y = 1.0;
						if( zl < 0 ) zl = 0, spacing.z = 1.0;
						if( zh >= dim.z ) zh = dim.z - 1, spacing.z = 1.0;

						const tgt::vec3 gradients[3] {
							( flowmap->voxel( xl, y, z ) - flowmap->voxel( xh, y, z ) ) / spacing.x,
							( flowmap->voxel( x, yl, z ) - flowmap->voxel( x, yh, z ) ) / spacing.y,
							( flowmap->voxel( x, y, zl ) - flowmap->voxel( x, y, zh ) ) / spacing.z
						};

						auto jacobi = Eigen::Matrix3f();
						for( int i = 0; i < 3; ++i ) jacobi.col( i ) = Eigen::Vector3f( gradients[i].x, gradients[i].y, gradients[i].z );

						const auto eigenvalues = Eigen::EigenSolver<Eigen::Matrix3f>( jacobi * jacobi.transpose(), false ).eigenvalues();
						const auto max = std::max( eigenvalues.x().real(), std::max( eigenvalues.y().real(), eigenvalues.z().real() ) );
						const auto ftle = std::log( std::sqrt( max ) ) / ( interval.y - interval.x );
						ftleField->voxel( x, y, z ) = ftle;

						if( max != 0.0f ) totalMin = std::min( totalMin, ftle );
						totalMax = std::max( totalMax, ftle );
					}
				}
			}

			std::cout << "[FlowMapProcessor] Finished FTLE Map Computation --> Range: Min = " << totalMin << ", Max = " << totalMax << std::endl;

			auto output = new Volume( ftleField.release(), tgt::vec3( 1, 1, 1 ), tgt::vec3::zero );
			output->getMetaDataContainer().addMetaData( "name", new StringMetaData( "FTLE" ) );
			_outportFTLE.setData( output );
		}

		/*void processGyre()
		{
		const auto volumeXHandle = _inport1.getData()->at( 0 );
		const auto volumeYHandle = _inport2.getData()->at( 0 );
		auto dimensions = volumeXHandle->getDimensions();

		if( !volumeXHandle || !volumeYHandle )
		{
		std::cout << "[FTLEProcessor]: Invalid Input" << std::endl;
		_outport.clear();
		return;
		}


		// --- Fill Flow Map --- //
		const auto volumeX = dynamic_cast<const VolumeRAM_Float*>( volumeXHandle->getRepresentation<VolumeRAM>() );
		const auto volumeY = dynamic_cast<const VolumeRAM_Float*>( volumeYHandle->getRepresentation<VolumeRAM>() );
		auto flowmap = std::unique_ptr<VolumeRAM_2xDouble>( new VolumeRAM_2xDouble( tgt::svec3( dimensions.x, dimensions.y, 1 ) ) );
		flowmap->fill( tgt::dvec2( 0.0, 0.0 ) );

		for( size_t t = 0; t < dimensions.z - std::ceil( _propertyStepsize.get() ); ++t )
		{
		std::cout << "[FTLEProcessor]: Filling Flow Field -> t = " << t << std::endl;

		#pragma omp parallel for
		for( long x = 0; x < static_cast<long>( dimensions.x ); ++x )
		{
		for( long y = 0; y < dimensions.y; ++y )
		{
		auto& flow = flowmap->voxel( x, y, 0 );
		const auto position = tgt::dvec2( x, y ) + flow;
		if( position.x >= 0.0 && position.x < dimensions.x && position.y >= 0.0 && position.y < dimensions.y )
		{
		const auto velocity = [dimensions, volumeX, volumeY] ( tgt::dvec2 p, double t )
		{
		const auto sample = [dimensions, volumeX, volumeY] ( tgt::dvec2 p, size_t t )
		{
		if( p.x >= 0.0 && p.x < dimensions.x && p.y >= 0.0 && p.y < dimensions.y )
		return tgt::dvec2( volumeX->voxel( p.x, p.y, t ), volumeY->voxel( p.x, p.y, t ) );
		return tgt::dvec2( 0.0, 0.0 );
		};

		const auto trilinear = [&sample] ( tgt::dvec2 p, size_t t )
		{
		const auto lx = std::floor( p.x ), hx = std::ceil( p.x );
		const auto ly = std::floor( p.y ), hy = std::ceil( p.y );
		const auto wx = p.x - lx, wy = p.y - ly;

		const auto c1 = sample( tgt::dvec2( lx, hy ), t );
		const auto c2 = sample( tgt::dvec2( hx, hy ), t );
		const auto c3 = sample( tgt::dvec2( lx, ly ), t );
		const auto c4 = sample( tgt::dvec2( hx, ly ), t );

		const auto i1 = c1 * wx + c2 * ( 1.0 - wx );
		const auto i2 = c3 * wx + c4 * ( 1.0 - wx );

		return i1 * wy + i2 * ( 1.0 - wy );
		};

		const auto i1 = trilinear( p, static_cast<size_t>( std::floor( t ) ) );
		const auto i2 = trilinear( p, static_cast<size_t>( std::ceil( t ) ) );
		const auto w = std::ceil( t ) - t;

		return i1 * w + i2 * ( 1.0 - w );
		};

		const auto h = static_cast<double>( _propertyStepsize.get() );
		const auto k1 = velocity( position, t );
		const auto k2 = velocity( position + k1 * h / 2.0, t + h / 2.0 );
		const auto k3 = velocity( position + k2 * h / 2.0, t + h / 2.0 );
		const auto k4 = velocity( position + k3 * h, t + h );
		flow += h * ( k1 + 2.0 * k2 + 2.0 * k3 + k4 ) / 6.0;
		}
		}
		}
		}

		// --- Fill FTLE Field --- //
		auto ftleField = std::unique_ptr<VolumeRAM_Float>( new VolumeRAM_Float( tgt::svec3( dimensions.x, dimensions.y, 1 ) ) );
		ftleField->fill( 0.0f );
		auto totalMin = std::numeric_limits<double>::max(), totalMax = std::numeric_limits<double>::min();

		#pragma omp parallel for
		for( long x = 0; x < static_cast<long>( dimensions.x ); ++x )
		{
		std::cout << "[FTLEProcessor]: Filling FTLE Field -> x=" << x << std::endl;

		for( long y = 0; y < dimensions.y; ++y )
		{
		auto spacing = tgt::dvec2( 2.0, 2.0 );
		auto xl = x - 1, xh = x + 1;
		auto yl = y - 1, yh = y + 1;

		if( xl < 0 ) xl = 0, spacing.x = 1.0;
		if( xh >= dimensions.x ) xh = dimensions.x - 1, spacing.x = 1.0;
		if( yl < 0 ) yl = 0, spacing.y = 1.0;
		if( yh >= dimensions.y ) yh = dimensions.y - 1, spacing.y = 1.0;

		auto jacobi = Eigen::Matrix2d();
		jacobi.col( 0 ) = reinterpret_cast<Eigen::Vector2d&&>( ( flowmap->voxel( xl, y, 0 ) - flowmap->voxel( xh, y, 0 ) ) / spacing.x );
		jacobi.col( 1 ) = reinterpret_cast<Eigen::Vector2d&&>( ( flowmap->voxel( x, yl, 0 ) - flowmap->voxel( x, yh, 0 ) ) / spacing.y );

		const auto eigenvalues = Eigen::EigenSolver<Eigen::Matrix2d>( jacobi.transpose() * jacobi, false ).eigenvalues();
		const auto max = std::max( eigenvalues.x().real(), eigenvalues.y().real() );
		const auto ftle = std::log( std::sqrt( max ) ) / dimensions.z;

		totalMin = std::min( totalMin, max );
		totalMax = std::max( totalMax, max );
		ftleField->voxel( x, y, 0 ) = ftle;
		}
		}

		std::cout << "[FTLEProcessor]: Finished with eigenvalue between " << totalMin << " and " << totalMax << std::endl;

		auto output = new Volume( ftleField.release(), tgt::vec3( 1, 1, 1 ), tgt::vec3::zero );
		output->getMetaDataContainer().addMetaData( "name", new StringMetaData( "FTLE" ) );
		_outport.setData( output );
		}*/

		VolumeListPort _inportX, _inportY, _inportZ;
		VolumePort _outportFlowmap, _outportFTLE;

		IntIntervalProperty _propertyTimesteps;
		FloatProperty _propertyStepsize;
		ButtonProperty _propertyComputeFlowmap, _propertyComputeFTLEField;
	};
}

#endif // VRN_FLOWMAPPROCESSOR_H