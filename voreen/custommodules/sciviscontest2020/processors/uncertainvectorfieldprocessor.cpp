/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "uncertainvectorfieldprocessor.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "Eigen/Eigenvalues"

#include <chrono>
#include <random>

namespace voreen {

UncertainVectorFieldProcessor::UncertainVectorFieldProcessor() : Processor(),
    _inportEnsemble( Port::INPORT, "inport_ensemble", "Ensemble Input" ),
    _inportMask( Port::INPORT, "inport_mask", "Volume Mask Input" ),
    _outportQ( Port::OUTPORT, "outport_q", "Q" ),
    _outportL( Port::OUTPORT, "ourport_lambda_2", "Lambda2" ),

    _propertySelectedMembers( "property_selected_members", "Members", Processor::VALID ),
    _propertyTimestep( "property_timestep", "Timestep", 0, 0, std::numeric_limits<int>::max(), Processor::VALID ),
    _propertySampleCount( "property_samplecount", "Sample Count", 64, 1, 512, Processor::VALID ),
    _propertyThresholdQ( "property_threshold_q", "Threshold Q", 0.0f, -1.0f, 1.0f, Processor::VALID ),
    _propertyThresholdL( "property_threshold_l", "Threshold Lambda2", 0.0f, -1.0f, 1.0f, Processor::VALID ),
    _propertyUpdate( "property_update", "Update", Processor::INVALID_RESULT )
{
    this->addPort( _inportEnsemble );
    this->addPort( _inportMask );
    this->addPort( _outportQ );
    this->addPort( _outportL );

    this->addProperty( _propertySelectedMembers );
    this->addProperty( _propertyTimestep );
    this->addProperty( _propertySampleCount );
    this->addProperty( _propertyThresholdQ );
    this->addProperty( _propertyThresholdL );
    this->addProperty( _propertyUpdate );

    _inportEnsemble.onNewData( LambdaFunctionCallback( [this]
    {
        _propertySelectedMembers.reset();
        for( const auto& member : _inportEnsemble.getData()->getMembers() )
            _propertySelectedMembers.addRow( member.getName() );
        _propertySelectedMembers.invalidate();

        _propertyTimestep.setMaxValue( std::max( 0, static_cast<int>( _inportEnsemble.getData()->getMinNumTimeSteps() ) - 1 ) );
    } ) );
    _propertyUpdate.onClick( MemberFunctionCallback<UncertainVectorFieldProcessor>( this, &UncertainVectorFieldProcessor::update ) );
}

void UncertainVectorFieldProcessor::update()
{
    const auto ensemble = _inportEnsemble.getData();
    const auto maskVolume = _inportMask.getData();
    if( !ensemble || !maskVolume )
    {
        _outportQ.clear();
        _outportL.clear();
        return;
    }

    const auto timestep = _propertyTimestep.get();
    const auto sampleCount = _propertySampleCount.get();
    const auto thresholdQ = _propertyThresholdQ.get();
    const auto thresholdL = _propertyThresholdL.get();

    const std::array<tgt::svec3, 7> offsets { tgt::svec3( 0, 0, 0 ),
        tgt::svec3( -1, 0, 0 ), tgt::svec3( 1, 0, 0 ),
        tgt::svec3( 0, -1, 0 ), tgt::svec3( 0, 1, 0 ),
        tgt::svec3( 0, 0, -1 ), tgt::svec3( 0, 0, 1 )
    };

    const auto& members = ensemble->getMembers();
    const auto& memberIndices = _propertySelectedMembers.get();
    const auto mask = maskVolume->getRepresentation<VolumeRAM>();
    const auto dim = mask->getDimensions();

    auto voxels = std::vector<tgt::ivec3>();
    voxels.reserve( mask->getNumVoxels() );
    for( auto x = 1; x < dim.x - 1; ++x ) for( auto y = 1; y < dim.y - 1; ++y ) for( auto z = 1; z < dim.z - 1; ++z )
    {
        auto addVoxel = true;
        for( const auto offset : offsets )
        {
            if( mask->getVoxelNormalized( tgt::svec3( x, y, z ) + offset ) <= 0.0f )
            {
                addVoxel = false;
                break;
            }
        }
        if( addVoxel ) voxels.push_back( tgt::ivec3( x, y, z ) );
    }
    voxels.shrink_to_fit();

    std::cout << "[UncertainVectorFieldProcessor]: Starting Computation -> timestep = " << timestep << ", sampleCount = " << sampleCount << ", thresholdQ = " << thresholdQ << ", thresholdL = " << thresholdL << ", voxelCount = " << voxels.size() << std::endl;

    auto volumeQ = std::unique_ptr<VolumeRAM_Float>( new VolumeRAM_Float( dim ) );
    auto volumeL = std::unique_ptr<VolumeRAM_Float>( new VolumeRAM_Float( dim ) );

    auto velocities = std::vector<tgt::vec3>( voxels.size() * offsets.size() * memberIndices.size() );
    for( size_t i = 0; i < memberIndices.size(); ++i )
    {
        const auto memberIndex = memberIndices[i];
        auto volumeX = VolumeRAMRepresentationLock( members[memberIndex].getTimeSteps()[timestep].getVolume( "U" ) );
        auto volumeY = VolumeRAMRepresentationLock( members[memberIndex].getTimeSteps()[timestep].getVolume( "V" ) );
        auto volumeZ = VolumeRAMRepresentationLock( members[memberIndex].getTimeSteps()[timestep].getVolume( "W" ) );

#pragma omp parallel for
        for( long j = 0; j < static_cast<long>( voxels.size() ); ++j )
        {
            for( size_t k = 0; k < offsets.size(); ++k )
            {
                const auto index = offsets.size() * ( ( memberIndices.size() * j ) + i ) + k;
                const auto position = tgt::svec3( voxels[j] ) + offsets[k];
                velocities[index] = tgt::vec3(
                    volumeX->getVoxelNormalized( position ),
                    volumeY->getVoxelNormalized( position ),
                    volumeZ->getVoxelNormalized( position )
                );
            }
        }
    }

    std::cout << "[UncertainVectorFieldProcessor]: Finished gathering velocities -> count = " << velocities.size() << std::endl;

#pragma omp parallel for
    for( long i = 0; i < static_cast<long>( voxels.size() ); ++i )
    {
        const auto timeBegin = std::chrono::high_resolution_clock::now();
        const auto x = voxels[i].x, y = voxels[i].y, z = voxels[i].z;

        // --- Gather velocity vectors --- //
        auto vectors = std::array<std::vector<tgt::vec3>, 7>();
        for( size_t i = 0; i < vectors.size(); ++i ) vectors[i].resize( memberIndices.size() );
        for( size_t j = 0; j < memberIndices.size(); ++j )
        {
            for( size_t k = 0; k < offsets.size(); ++k )
            {
                const auto index = offsets.size() * ( ( memberIndices.size() * i ) + j ) + k;
                vectors[k][j] = velocities[index];
            }
        }

        // --- Calculate 21D-Mean --- //
        auto mean = Eigen::Matrix<float, 21, 1>();
        for( size_t i = 0; i < 21; ++i )
        {
            mean[i] = 0.0f;
            for( size_t j = 0; j < memberIndices.size(); ++j )
                mean[i] += vectors[i / 3][j][i % 3];
            mean[i] /= memberIndices.size();
        }

        // --- Calculate 21x21 Covariance Matrix --- //
        auto covariance = Eigen::Matrix<float, 21, 21>();
        for( size_t i = 0; i < 21; ++i )
        {
            for( size_t j = 0; j < 21; ++j )
            {
                covariance( i, j ) = 0.0f;
                for( size_t k = 0; k < memberIndices.size(); ++k )
                    covariance( i, j ) += ( vectors[i / 3][k][i % 3] - mean[i] ) * ( vectors[j / 3][k][j % 3] - mean[j] );
                covariance( i, j ) /= memberIndices.size();
            }
        }

        const auto cholesky = Eigen::LLT<Eigen::Matrix<float, 21, 21>>( covariance ).matrixL().toDenseMatrix();
        auto engine = std::mt19937( std::random_device()( ) );
        auto normalDistribution = std::normal_distribution<float>( 0.0f, 1.0f );

        auto countQ = 0, countL = 0;
        for( size_t i = 0; i < sampleCount; ++i )
        {
            // --- Generate Sample 21D-Vector --- //
            auto standardNormal = Eigen::Matrix<float, 21, 1>();
            for( size_t j = 0; j < 21; ++j )
                standardNormal[j] = normalDistribution( engine );
            const auto sample = mean + cholesky * standardNormal;

            // --- Calculate Jacobian --- //
            auto jacobian = Eigen::Matrix3f();
            jacobian <<
                sample[6] - sample[3], sample[12] - sample[9], sample[18] - sample[15],
                sample[7] - sample[4], sample[13] - sample[10], sample[19] - sample[16],
                sample[8] - sample[5], sample[14] - sample[11], sample[20] - sample[17];
            jacobian /= 2.0f;

            // --- Calculate S and Omega --- //
            const auto transpose = jacobian.transpose();
            const auto S = ( jacobian + transpose ) / 2.0f;
            const auto O = ( jacobian - transpose ) / 2.0f;

            // --- Q Criterium --- //
            auto eigenvalues = ( S.transpose() * S ).eigenvalues();
            const auto spectralNormS = std::max( { eigenvalues.x().real(), eigenvalues.y().real(), eigenvalues.z().real() } );
            eigenvalues = ( O.transpose() * O ).eigenvalues();
            const auto spectralNormO = std::max( { eigenvalues.x().real(), eigenvalues.y().real(), eigenvalues.z().real() } );
            if( 0.5f * ( spectralNormO - spectralNormS ) > thresholdQ ) ++countQ;

            // --- Lambda 2 Criterium --- //
            auto lambda2 = 0.0f;
            eigenvalues = ( S * S + O * O ).eigenvalues();
            if( eigenvalues.x().real() < eigenvalues.y().real() )
            {
                if( eigenvalues.y().real() < eigenvalues.z().real() ) lambda2 = eigenvalues.y().real();
                else if( eigenvalues.x().real() < eigenvalues.z().real() ) lambda2 = eigenvalues.z().real();
                else lambda2 = eigenvalues.x().real();
            }
            else
            {
                if( eigenvalues.x().real() < eigenvalues.z().real() ) lambda2 = eigenvalues.x().real();
                else if( eigenvalues.y().real() < eigenvalues.z().real() ) lambda2 = eigenvalues.z().real();
                else lambda2 = eigenvalues.y().real();
            }
            if( lambda2 < thresholdL ) ++countL;
        }

        volumeQ->voxel( x, y, z ) = static_cast<float>( countQ ) / sampleCount;
        volumeL->voxel( x, y, z ) = static_cast<float>( countL ) / sampleCount;

        const auto timeEnd = std::chrono::high_resolution_clock::now();
        //if( omp_get_thread_num() == 0 && i % ( voxels.size() / 100 / omp_get_num_threads() ) == 0 )
        //{
        //    const auto time = std::chrono::duration_cast<std::chrono::microseconds>( timeEnd - timeBegin ).count();
        //    const auto remaining = voxels.size() / omp_get_num_threads() - i;
        //    const auto seconds = time / 1000000.0 * remaining;
        //    const auto minutes = static_cast<int>( seconds / 60.0 );
        //
        //    std::cout << "[UncertainVectorFieldProcessor]: Remaining Time: " << minutes << " min " << static_cast<int>( seconds ) % 60 << " s" << std::endl;
        //}
    }

    // --- Update Outports --- //
    auto volume = new Volume( volumeQ.release(), tgt::vec3( 1.0f, 1.0f, 1.0f ), tgt::vec3::zero );
    volume->getMetaDataContainer().addMetaData( "name", new StringMetaData( "Q" ) );
    _outportQ.setData( volume );

    volume = new Volume( volumeL.release(), tgt::vec3( 1.0f, 1.0f, 1.0f ), tgt::vec3::zero );
    volume->getMetaDataContainer().addMetaData( "name", new StringMetaData( "LAMBDA2" ) );
    _outportL.setData( volume );
}

}