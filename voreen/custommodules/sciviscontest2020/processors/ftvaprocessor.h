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

#ifndef VRN_FTVAPROCESSOR_H
#define VRN_FTVAPROCESSOR_H

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/ports/volumeport.h"
#include "Eigen/Eigenvalues"

namespace voreen {
class FTVAProcessor : public Processor {
public:
    FTVAProcessor() : Processor(),
        _inport( Port::INPORT, "inport", "List of Flow Maps" ),
        _outport( Port::OUTPORT, "outport", "FTVA Field" ),
        _propertyTimesteps( "property_timesteps", "Timesteps", 0, 0, std::numeric_limits<int>::max(), Processor::VALID ),
        _propertyComputeFTVAField( "property_compute_ftva_field", "Compute FTVA Field", Processor::VALID )
    {
        this->addPort( _inport );
        this->addPort( _outport );

        this->addProperty( _propertyTimesteps );
        this->addProperty( _propertyComputeFTVAField );

        _propertyComputeFTVAField.onClick( MemberFunctionCallback<FTVAProcessor>( this, &FTVAProcessor::computeFTVAField ) );
    }
    Processor* create() const override
    {
        return new FTVAProcessor();
    }
    std::string getClassName() const override
    {
        return "FTVAProcessor";
    }
    std::string getCategory() const override
    {
        return "Flow";
    }

private:
    void process() override {}
    void computeFTVAField()
    {
        const auto volumes = _inport.getData();
        if( !volumes )
        {
            _outport.clear();
            return;
        }

        const auto dim = volumes->first()->getDimensions();
        const std::array<tgt::ivec3, 7 > offsets { tgt::ivec3( 0, 0, 0 ),
            tgt::ivec3( -1, 0, 0 ), tgt::ivec3( 1, 0, 0 ),
            tgt::ivec3( 0, -1, 0 ), tgt::ivec3( 0, 1, 0 ),
            tgt::ivec3( 0, 0, -1 ), tgt::ivec3( 0, 0, 1 )
        };

        auto countField = std::unique_ptr<VolumeRAM_Int32>( new VolumeRAM_Int32( dim ) );
        auto meanField = std::unique_ptr<VolumeRAM_3xFloat>( new VolumeRAM_3xFloat( dim ) );
        auto covarianceField = std::unique_ptr<VolumeRAM_Mat3Float>( new VolumeRAM_Mat3Float( dim ) );
        auto ftvaField = std::unique_ptr<VolumeRAM_Float>( new VolumeRAM_Float( dim ) );

        countField->fill( 0 );
        meanField->fill( tgt::vec3::zero );
        covarianceField->fill( tgt::mat3::zero );

        // --- Calculate Means --- //
        for( size_t i = 0; i < volumes->size(); ++i )
        {
            const auto volume = dynamic_cast<const VolumeRAM_3xFloat*>( volumes->at( i )->getRepresentation<VolumeRAM>() );
            std::cout << "[FTVAProcessor]: Computing means -> volume = " << i << std::endl;

#pragma omp parallel for
            for( long x = 0; x < dim.x; ++x ) for( long y = 0; y < dim.y; ++y ) for( long z = 0; z < dim.z; ++z )
            {
                auto& count = countField->voxel( x, y, z );
                auto& mean = meanField->voxel( x, y, z );

                for( const auto offset : offsets )
                {
                    const auto position = tgt::ivec3( x, y, z ) + offset;
                    if( position.x >= 0 && position.x < dim.x && position.y >= 0 && position.y < dim.y && position.z >= 0 && position.z < dim.z )
                    {
                        count += 1;
                        mean += volume->voxel( position.x, position.y, position.z );
                    }
                }
            }
        }

        // --- Calculate Covariance Matrices --- //
        for( size_t i = 0; i < volumes->size(); ++i )
        {
            const auto volume = dynamic_cast<const VolumeRAM_3xFloat*>( volumes->at( i )->getRepresentation<VolumeRAM>() );
            std::cout << "[FTVAProcessor]: Computing covariance matrices -> volume = " << i << std::endl;

#pragma omp parallel for
            for( long x = 0; x < dim.x; ++x ) for( long y = 0; y < dim.y; ++y ) for( long z = 0; z < dim.z; ++z )
            {
                const auto count = countField->voxel( x, y, z );
                const auto mean = meanField->voxel( x, y, z ) / static_cast<float>( count );
                auto& covariance = covarianceField->voxel( x, y, z );

                for( int j = 0; j < 3; ++j )
                {
                    for( int k = 0; k < 3; ++k )
                    {
                        for( const auto offset : offsets )
                        {
                            const auto position = tgt::ivec3( x, y, z ) + offset;
                            if( position.x >= 0 && position.x < dim.x && position.y >= 0 && position.y < dim.y && position.z >= 0 && position.z < dim.z )
                            {
                                const auto flow = volume->voxel( position.x, position.y, position.z );
                                covariance[j][k] += ( flow[j] - mean[j] ) * ( flow[k] - mean[k] );
                            }
                        }
                    }
                }
            }
        }

        // --- Calculate FTVA Field --- //
#pragma omp parallel for
        for( long x = 0; x < dim.x; ++x )
        {
            std::cout << "[FTVAProcessor]: Computing FTVA field -> x = " << x << std::endl;
            for( long y = 0; y < dim.y; ++y ) for( long z = 0; z < dim.z; ++z )
            {
                const auto count = countField->voxel( x, y, z );
                const auto covariance = covarianceField->voxel( x, y, z ) / static_cast<float>( count - 1 );

                using MatrixType = Eigen::Matrix<float, 3, 3, Eigen::StorageOptions::RowMajor | Eigen::StorageOptions::AutoAlign>;
                const auto eigenvalues = Eigen::EigenSolver<MatrixType>( reinterpret_cast<const MatrixType&>( covariance ), false ).eigenvalues();
                const auto max = std::max( eigenvalues.x().real(), std::max( eigenvalues.y().real(), eigenvalues.z().real() ) );

                ftvaField->voxel( x, y, z ) = std::log( std::sqrt( max ) ) / _propertyTimesteps.get();
            }
        }

        auto output = new Volume( ftvaField.release(), tgt::vec3( 1, 1, 1 ), tgt::vec3::zero );
        output->getMetaDataContainer().addMetaData( "name", new StringMetaData( "FTVA" ) );
        _outport.setData( output );

        /*
        const auto volumes = _inport.getData();
        const auto dim = volumes->first()->getDimensions();

        auto ftvaField = std::unique_ptr<VolumeRAM_Float>( new VolumeRAM_Float( dim ) );

        auto pointsVolume = std::vector<std::vector<std::vector<std::vector<tgt::vec3>>>>( dim.x, std::vector<std::vector<std::vector<tgt::vec3>>>( dim.y, std::vector<std::vector<tgt::vec3>>( dim.z, std::vector<tgt::vec3>( volumes->size() ) ) ) ); // A lot of memory
        for( size_t i = 0; i < volumes->size(); ++i )
        {
            std::cout << "[FTVAProcessor]: Gathering points -> volume: " << i << std::endl;
            const auto volume = dynamic_cast<const VolumeRAM_3xFloat*>( volumes->at( i )->getRepresentation<VolumeRAM>() );

#pragma omp parallel for
            for( long x = 0; x < static_cast<long>( dim.x ); ++x ) for( long y = 0; y < dim.y; ++y ) for( long z = 0; z < dim.z; ++z )
            {
                pointsVolume[x][y][z][i] = volume->voxel( x, y, z );
            }
        }

#pragma omp parallel for
        for( long x = 0; x < static_cast<long>( dim.x ); ++x )
        {
            std::cout << "[FTVAProcessor]: Computing FTVA -> x: " << x << std::endl;
            for( long y = 0; y < dim.y; ++y ) for( long z = 0; z < dim.z; ++z )
            {
                auto points = std::vector<tgt::vec3>();
                points.reserve( 7 * volumes->size() );
                for( size_t i = 0; i < volumes->size(); ++i )
                {
                    points.push_back( pointsVolume[x][y][z][i] );
                    if( x > 0 ) points.push_back( pointsVolume[x - 1][y][z][i] );
                    if( x < dim.x - 1 ) points.push_back( pointsVolume[x + 1][y][z][i] );
                    if( y > 0 ) points.push_back( pointsVolume[x][y - 1][z][i] );
                    if( y < dim.y - 1 ) points.push_back( pointsVolume[x][y + 1][z][i] );
                    if( z > 0 ) points.push_back( pointsVolume[x][y][z - 1][i] );
                    if( z < dim.z - 1 ) points.push_back( pointsVolume[x][y][z + 1][i] );
                }

                auto mean = tgt::vec3();
                for( const auto point : points ) mean += point;
                mean /= static_cast<float>( points.size() );

                auto covariance = Eigen::Matrix3f();
                for( int i = 0; i < 3; ++i )
                {
                    for( int j = 0; j < 3; ++j )
                    {
                        covariance( i, j ) = 0.0f;
                        for( size_t k = 0; k < points.size(); ++k )
                        {
                            covariance( i, j ) += ( points[k][i] - mean[i] ) * ( points[k][j] - mean[j] );
                        }
                        covariance( i, j ) /= points.size() - 1;
                    }
                }

                const auto eigenvalues = Eigen::EigenSolver<Eigen::Matrix3f>( covariance, false ).eigenvalues();
                const auto max = std::max( eigenvalues.x().real(), std::max( eigenvalues.y().real(), eigenvalues.z().real() ) );

                ftvaField->voxel( x, y, z ) = std::log( std::sqrt( max ) ) / _propertyTimesteps.get();
            }
        }

        auto output = new Volume( ftvaField.release(), tgt::vec3( 1, 1, 1 ), tgt::vec3::zero );
        output->getMetaDataContainer().addMetaData( "name", new StringMetaData( "FTVA" ) );
        _outport.setData( output );*/
    }

    VolumeListPort _inport;
    VolumePort _outport;

    IntProperty _propertyTimesteps;
    ButtonProperty _propertyComputeFTVAField;
};

}

#endif // VRN_FTVAPROCESSOR_H