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

#include "approximateparallelvectors.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"

#include "modules/flowanalysis/processors/volume/acceleration.h"
#include "custommodules/sciviscontest2020/processors/vortexprocessor.h"

#include "Eigen/Eigenvalues"

namespace voreen {
ApproximateParallelVectors::ApproximateParallelVectors()
    : Processor()
    , _inportEnsemble( Port::INPORT, "inport_ensemble", "Ensemble" )
    , _inportMask( Port::INPORT, "inport_mask", "Mask" )
    , _outportA( Port::OUTPORT, "outport_a", "Volume A" )
    , _outportB( Port::OUTPORT, "outport_b", "Volume B" )
    , _outportEps( Port::OUTPORT, "outport_eps", "Volume Epsilon" )
    , _propertyTime( "property_time", "Time", 0, 0, std::numeric_limits<float>::max() )
    , _propertyEpsilon( "property_epsilon", "Epsilon Threshold", -1.0f, -1.0f, 1.0f )
    , _propertyUseAcceleration( "property_use_aceleration", "Use Acceleration", false )
{
    this->addPort( _inportEnsemble );
    this->addPort( _inportMask );
    this->addPort( _outportA );
    this->addPort( _outportB );
    this->addPort( _outportEps );

    this->addProperty( _propertyTime );
    this->addProperty( _propertyEpsilon );
    this->addProperty( _propertyUseAcceleration );

    _inportEnsemble.onNewData( LambdaFunctionCallback( [this] {
        const auto* ensemble = _inportEnsemble.getData();
        _propertyTime.setMinValue(ensemble->getStartTime());
        _propertyTime.setMaxValue(ensemble->getEndTime());
    } ) );
}

Processor* ApproximateParallelVectors::create() const {
    return new ApproximateParallelVectors();
}
std::string ApproximateParallelVectors::getClassName() const {
    return "ApproximateParallelVectors";
}
std::string ApproximateParallelVectors::getCategory() const {
    return "Vortex Processing";
}

void ApproximateParallelVectors::process() {

    const auto ensemble = _inportEnsemble.getData();
    const auto maskVolumeBase = _inportMask.getData();
    const auto time = _propertyTime.get();

    // --- Gather Voxels --- //
    const auto volumeMask = VolumeRAMRepresentationLock( maskVolumeBase );
    const auto dim = volumeMask->getDimensions();
    const auto& members = ensemble->getMembers();

    auto voxels = std::vector<size_t>();
    voxels.reserve( volumeMask->getNumVoxels() );
    for( size_t i = 0; i < volumeMask->getNumVoxels(); ++i )
        if( volumeMask->getVoxelNormalized( i ) > 0.0 )
            voxels.push_back( i );
    voxels.shrink_to_fit();

    const auto getVelocity = [time, dim, &members, &voxels] (size_t index)
    {
        size_t timeStep = members[index].getTimeStep(time);
        const auto volumeLockU = VolumeRAMRepresentationLock(members[index].getTimeSteps()[timeStep].getVolume("U" ) );
        const auto volumeLockV = VolumeRAMRepresentationLock(members[index].getTimeSteps()[timeStep].getVolume("V" ) );
        const auto volumeLockW = VolumeRAMRepresentationLock(members[index].getTimeSteps()[timeStep].getVolume("W" ) );

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
    auto volumeA = std::unique_ptr<VolumeRAM_3xFloat>( new VolumeRAM_3xFloat( dim ) );
    auto volumeC = std::unique_ptr<VolumeRAM_Mat3Float>( new VolumeRAM_Mat3Float( dim ) );
    auto volumeB = std::unique_ptr<VolumeRAM_3xFloat>( new VolumeRAM_3xFloat( dim ) );
    auto volumeEps = std::unique_ptr<VolumeRAM_Float>( new VolumeRAM_Float( dim ) );

    volumeA->fill( tgt::dvec3::zero );
    volumeC->fill( tgt::dmat3::zero );
    volumeB->fill( tgt::dvec3::zero );
    volumeEps->fill( 0.0 );

    // --- Calculate A --- //
    for(size_t i = 0; i < members.size(); ++i )
    {
        const auto volumeBaseVelocity = getVelocity( i );
        const auto volumeVelocity = dynamic_cast<const VolumeRAM_3xFloat*>( volumeBaseVelocity->getRepresentation<VolumeRAM>() );

#pragma omp parallel for
        for( long i = 0; i < static_cast<long>( voxels.size() ); ++i )
        {
            const auto voxelIndex = voxels[i];
            volumeA->voxel( voxelIndex ) += volumeVelocity->voxel( voxelIndex ) / static_cast<float>( members.size() );
        }

        if( _propertyUseAcceleration.get() )
        {
            auto volumeJacobi = std::unique_ptr<VolumeRAM_Mat3Float>( new VolumeRAM_Mat3Float( dim ) );
            VortexProcessor::Process( *volumeVelocity, *volumeMask.operator->(), *volumeJacobi, nullptr, nullptr, nullptr );

            auto volumeAcceleration = std::unique_ptr<VolumeRAM_3xFloat>( new VolumeRAM_3xFloat( dim ) );
            Acceleration::Process(*volumeJacobi, *volumeVelocity, *volumeAcceleration );

#pragma omp parallel for
            for( long i = 0; i < static_cast<long>( voxels.size() ); ++i )
            {
                const auto voxelIndex = voxels[i];
                volumeA->voxel( voxelIndex ) += volumeAcceleration->voxel( voxelIndex ) / static_cast<float>( members.size() );
            }
        }
    }

    // --- Calculate C --- //
    for(size_t i = 0; i < members.size(); ++i )
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

            auto volumeAcceleration = std::unique_ptr<VolumeRAM_3xFloat>( new VolumeRAM_3xFloat( dim ) );
            Acceleration::Process(*volumeJacobi, *volumeVelocity, *volumeAcceleration );

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