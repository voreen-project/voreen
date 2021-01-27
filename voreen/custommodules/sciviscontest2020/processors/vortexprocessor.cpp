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

#include "vortexprocessor.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorgradient.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"


#include "Eigen/Eigenvalues"

namespace voreen {

const std::string VortexProcessor::loggerCat_("voreen.flowanalysis.VortexProcessor");

VortexProcessor::VortexProcessor()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , inputVolume_(Port::INPORT, "vortexprocessor.inputVolume", "Volume Input")
    , outputVolumeJacobi_(Port::OUTPORT, "vortexprocessor.outputVolumeJacobi", "Jacobi matrix")
    , outputVolumeDelta_(Port::OUTPORT, "vortexprocessor.outputVolumeDelta", "Delta-criterion")
    , outputVolumeQ_(Port::OUTPORT, "vortexprocessor.outputVolumeQ", "Q-criterion")
    , outputVolumeLamda2_(Port::OUTPORT, "vortexprocessor.outputVolumeLambda2", "Lambda2-criterion")
{
    addPort(inputVolume_);
    inputVolume_.addCondition(new PortConditionVolumeChannelCount(3));

    addPort(outputVolumeJacobi_);
    addPort(outputVolumeDelta_);
    addPort(outputVolumeQ_);
    addPort(outputVolumeLamda2_);
}

Processor* VortexProcessor::create() const {
    return new VortexProcessor();
}

bool VortexProcessor::isReady() const {
    if(!inputVolume_.isReady()) {
        setNotReadyErrorMessage("No input");
        return false;
    }

    if(!outputVolumeJacobi_.isReady() &&
        !outputVolumeDelta_.isReady() &&
        !outputVolumeQ_.isReady() &&
        !outputVolumeLamda2_.isReady())
    {
        setNotReadyErrorMessage("No outport connected");
        return false;
    }

    return true;
}

VortexProcessorIO VortexProcessor::prepareComputeInput() {

    auto inputVolume = inputVolume_.getThreadSafeData();
    if(!inputVolume) {
        throw InvalidInputException("No input", InvalidInputException::S_ERROR);
    }

    const tgt::svec3 dim = inputVolume->getDimensions();

    std::unique_ptr<VolumeRAM_Mat3Float> volumeJacobi;
    if(outputVolumeJacobi_.isConnected()) {
        volumeJacobi.reset(new VolumeRAM_Mat3Float( dim ));
        volumeJacobi->fill(tgt::mat3::identity);
    }
    std::unique_ptr<VolumeRAM_Float> volumeDelta;
    if(outputVolumeDelta_.isConnected()) {
        volumeDelta.reset(new VolumeRAM_Float( dim ));
        volumeDelta->fill(0.0f);
    }
    std::unique_ptr<VolumeRAM_Float> volumeQ;
    if(outputVolumeQ_.isConnected()) {
        volumeQ.reset(new VolumeRAM_Float( dim ));
        volumeQ->fill(0.0f);
    }
    std::unique_ptr<VolumeRAM_Float> volumeLambda2;
    if(outputVolumeLamda2_.isConnected()) {
        volumeLambda2.reset(new VolumeRAM_Float( dim ));
        volumeLambda2->fill(0.0f);
    }

    return VortexProcessorIO{
            std::move(inputVolume),
            std::move(volumeJacobi),
            std::move(volumeDelta),
            std::move(volumeQ),
            std::move(volumeLambda2)
    };
}

VortexProcessorIO VortexProcessor::compute(VortexProcessorIO input, ProgressReporter& progressReporter) const {

    auto& inputVolume = input.inputVolume;
    RealWorldMapping rwm = inputVolume->getRealWorldMapping();
    VolumeRAMRepresentationLock inputVolumeData(inputVolume);

    auto& volumeJacobi = input.outputJacobi;
    auto& volumeDelta = input.outputDelta;
    auto& volumeQ = input.outputQ;
    auto& volumeLambda2 = input.outputLambda2;

    tgt::Vector3<long> dimensions = inputVolumeData->getDimensions();
    tgt::vec3 spacing = inputVolume->getSpacing();

    const std::array<tgt::svec3, 7> offsets { tgt::svec3( 0, 0, 0 ),
                                              tgt::svec3( -1, 0, 0 ), tgt::svec3( 1, 0, 0 ),
                                              tgt::svec3( 0, -1, 0 ), tgt::svec3( 0, 1, 0 ),
                                              tgt::svec3( 0, 0, -1 ), tgt::svec3( 0, 0, 1 )
    };

    ThreadedTaskProgressReporter progress(progressReporter, dimensions.z);
    bool aborted = false;

#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for
#endif
    for (long z = 3; z < dimensions.z-3; z++) {
        if (aborted) {
            continue;
        }
        for (long y = 3; y < dimensions.y-3; y++) {
            for (long x = 3; x < dimensions.x-3; x++) {
                tgt::svec3 pos(x, y, z);

                // --- Compute Jacobian --- //
                auto jacobi = Eigen::Matrix3f();
                for (size_t i = 0; i < 3; i++) {
                    for (size_t j = 0; j < 3; j++) {
                        const auto b3 = rwm.normalizedToRealWorld(inputVolumeData->getVoxelNormalized(pos + offsets[(j + 1) * 2] * size_t(3), i));
                        const auto b2 = rwm.normalizedToRealWorld(inputVolumeData->getVoxelNormalized(pos + offsets[(j + 1) * 2] * size_t(2), i));
                        const auto b1 = rwm.normalizedToRealWorld(inputVolumeData->getVoxelNormalized(pos + offsets[(j + 1) * 2], i));
                        const auto c  = rwm.normalizedToRealWorld(inputVolumeData->getVoxelNormalized(pos, i));
                        const auto f1 = rwm.normalizedToRealWorld(inputVolumeData->getVoxelNormalized(pos + offsets[(j + 1) * 2 - 1], i));
                        const auto f2 = rwm.normalizedToRealWorld(inputVolumeData->getVoxelNormalized(pos + offsets[(j + 1) * 2 - 1] * size_t(2), i));
                        const auto f3 = rwm.normalizedToRealWorld(inputVolumeData->getVoxelNormalized(pos + offsets[(j + 1) * 2 - 1] * size_t(3), i));

                        const auto db2 = (b1 - b3) / 2.0f;
                        const auto db1 = (c - b2) / 2.0f;
                        const auto dc  = (f1 - b1) / 2.0f;
                        const auto df1 = (f2 - c) / 2.0f;
                        const auto df2 = (f3 - f1) / 2.0f;

                        const auto t = 0.5f;
                        jacobi(i, j) = 3.0f * (2.0f * b1 - 2.0f * f1 + db1 + df1) * t * t +
                                       2.0f * (-3.0f * b1 + 3.0f * f1 - 2.0f * db1 - df1) * t
                                       + db1;
                    }
                }

                if(volumeJacobi) {
                    for (size_t i = 0; i < 3; i++) {
                        for (size_t j = 0; j < 3; j++) {
                            volumeJacobi->voxel(pos)[i][j] = jacobi(i, j);
                        }
                    }
                }

                // --- Compute Delta Criterion --- //
                auto eigenvalues = Eigen::Vector3cf();
                if (volumeDelta) {
                    eigenvalues = jacobi.eigenvalues();
                    volumeDelta->voxel(pos) = std::max({eigenvalues.x().imag(), eigenvalues.y().imag(), eigenvalues.z().imag()});
                }

                // --- Compute Strain Rate and Vorticity --- //
                const auto jacobiT = jacobi.transpose();
                const auto S = 0.5f * (jacobi + jacobiT);
                const auto O = 0.5f * (jacobi - jacobiT);

                // --- Compute Q Criterion --- //
                if (volumeQ) {
                    eigenvalues = (S.transpose() * S).eigenvalues();
                    const auto spectralNormS = std::max({eigenvalues.x().real(), eigenvalues.y().real(), eigenvalues.z().real()});
                    eigenvalues = (O.transpose() * O).eigenvalues();
                    const auto spectralNormO = std::max({eigenvalues.x().real(), eigenvalues.y().real(), eigenvalues.z().real()});
                    volumeQ->voxel(pos) = 0.5f * (spectralNormO - spectralNormS);
                }

                // --- Compute Lambda2 Criterion --- //
                if (volumeLambda2) {
                    eigenvalues = (S * S + O * O).eigenvalues();
                    if (eigenvalues.x().real() < eigenvalues.y().real())
                        if (eigenvalues.y().real() < eigenvalues.z().real()) {
                            volumeLambda2->voxel(pos) = eigenvalues.y().real();
                        }
                        else if (eigenvalues.x().real() < eigenvalues.z().real()) {
                            volumeLambda2->voxel(pos) = eigenvalues.z().real();
                        }
                        else {
                            volumeLambda2->voxel(pos) = eigenvalues.x().real();
                        }
                    else if (eigenvalues.x().real() < eigenvalues.z().real()) {
                        volumeLambda2->voxel(pos) = eigenvalues.x().real();
                    }
                    else if (eigenvalues.y().real() < eigenvalues.z().real()) {
                        volumeLambda2->voxel(pos) = eigenvalues.z().real();
                    }
                    else {
                        volumeLambda2->voxel(pos) = eigenvalues.y().real();
                    }

                    if(volumeLambda2->voxel(pos) > 0.0f) {
                        volumeLambda2->voxel(pos) = 0.0f;
                    }
                    else {
                        volumeLambda2->voxel(pos) = -volumeLambda2->voxel(pos);
                    }
                }
            }
        }

        if (progress.reportStepDone()) {
#ifdef VRN_MODULE_OPENMP
            #pragma omp critical
            aborted = true;
#else
            aborted = true;
            break;
#endif
        }
    }

    if (aborted) {
        throw boost::thread_interrupted();
    }

    progressReporter.setProgress(1.0f);

    return input;
}

void VortexProcessor::processComputeOutput(VortexProcessorIO output) {

    tgt::vec3 spacing = output.inputVolume->getSpacing();
    tgt::vec3 offset = output.inputVolume->getOffset();

    if(output.outputJacobi) {
        auto* volume = new Volume(output.outputJacobi.release(), spacing, offset);
        volume->setMetaDataValue<StringMetaData>("name", "jacobi");
        outputVolumeJacobi_.setData(volume);
        // No need to set real world mapping, the matrix contains real world values already.
    }

    if(output.outputDelta) {
        auto* volume = new Volume(output.outputDelta.release(), spacing, offset);
        volume->setMetaDataValue<StringMetaData>("name", "delta");
        outputVolumeDelta_.setData(volume);
    }

    if(output.outputQ) {
        auto* volume = new Volume(output.outputQ.release(), spacing, offset);
        volume->setMetaDataValue<StringMetaData>("name", "Q");
        outputVolumeQ_.setData(volume);
    }

    if(output.outputLambda2) {
        auto* volume = new Volume(output.outputLambda2.release(), spacing, offset);
        volume->setMetaDataValue<StringMetaData>("name", "lambda2");
        outputVolumeLamda2_.setData(volume);
    }
}


void VortexProcessor::Process( const VolumeRAM& velocity, const VolumeRAM& mask, VolumeRAM_Mat3Float& outJacobi, VolumeRAM_Float* outDelta, VolumeRAM_Float* outLambda2, VolumeRAM_Float* outQ )
{
    const std::array<tgt::svec3, 7> offsets { tgt::svec3( 0, 0, 0 ),
                                              tgt::svec3( -1, 0, 0 ), tgt::svec3( 1, 0, 0 ),
                                              tgt::svec3( 0, -1, 0 ), tgt::svec3( 0, 1, 0 ),
                                              tgt::svec3( 0, 0, -1 ), tgt::svec3( 0, 0, 1 )
    };

    const auto dim = velocity.getDimensions();

    auto voxels = std::vector<tgt::ivec3>();
    voxels.reserve( mask.getNumVoxels() );
    for( auto x = 3; x < dim.x - 3; ++x )
        for( auto y = 3; y < dim.y - 3; ++y )
            for( auto z = 3; z < dim.z - 3; ++z )
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
}


} // namespace voreen
