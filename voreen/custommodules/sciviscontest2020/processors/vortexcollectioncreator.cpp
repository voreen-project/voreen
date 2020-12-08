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

#include "vortexcollectioncreator.h"

#include "modules/flowanalysis/processors/geometry/parallelvectors.h"
#include "modules/flowanalysis/processors/geometry/corelinecreator.h"
#include "modules/flowanalysis/processors/volume/acceleration.h"

#include "custommodules/sciviscontest2020/processors/curlprocessor.h"
#include "custommodules/sciviscontest2020/processors/vortexprocessor.h"
#include "custommodules/sciviscontest2020/processors/rotationaldirectionprocessor.h"
#include "custommodules/sciviscontest2020/processors/vortextracking.h"

#include "voreen/core/memorymanager/volumememorymanager.h"

#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"

namespace voreen{

VortexCollectionCreator::VortexCollectionCreator() : Processor(),
    _inportEnsemble( Port::INPORT, "inport_ensemble", "Ensemble Input" ),
    _outportVortexCollection( Port::OUTPORT, "outport_vortex_collection", "Vortex Collection Output", false ),
    _propertySelectedMembers( "property_selected_members", "Members", Processor::VALID ),
    _propertyTimestepInterval( "property_timestep_interval", "Timesteps", 0, 0, std::numeric_limits<int>::max(), 0, std::numeric_limits<int>::max(), Processor::VALID ),
    _propertyCorelineLength( "property_coreline_length", "Minimal Coreline Length", 64, 2, 256, Processor::VALID ),
    _propertyUpdateButton( "property_update_button", "Update", Processor::VALID ),
    _propertyFileDialog( "property_file_dialog", "Save File", "Select File...", "", "Vortex Collection (*.vc)", FileDialogProperty::FileMode::SAVE_FILE, Processor::VALID ),
    _propertySaveButton( "property_save_button", "Save File", Processor::VALID )
{
    this->addPort( _inportEnsemble );
    this->addPort( _outportVortexCollection );

    this->addProperty( _propertySelectedMembers );
    this->addProperty( _propertyTimestepInterval );
    this->addProperty( _propertyCorelineLength );
    this->addProperty( _propertyUpdateButton );
    this->addProperty( _propertyFileDialog );
    this->addProperty( _propertySaveButton );

    _inportEnsemble.onNewData( LambdaFunctionCallback( [this]
    {
        _propertySelectedMembers.blockCallbacks( true );
        _propertySelectedMembers.reset();
        for( const auto& member : _inportEnsemble.getData()->getMembers() )
            _propertySelectedMembers.addRow( member.getName(), member.getColor() );
        _propertySelectedMembers.blockCallbacks( false );
        _propertySelectedMembers.invalidate();
    } ) );
    _propertySelectedMembers.onChange( LambdaFunctionCallback( [this]
    {
        if( !_inportEnsemble.hasData() ) return;

        auto ensemble = EnsembleDataset();
        for( const auto member : _propertySelectedMembers.get() )
            ensemble.addMember( _inportEnsemble.getData()->getMembers()[member] );

        _propertyTimestepInterval.blockCallbacks( true );
        _propertyTimestepInterval.setMaxValue( static_cast<int>( ensemble.getMinNumTimeSteps() - 1 ) );
        _propertyTimestepInterval.blockCallbacks( false );
        _propertyTimestepInterval.invalidate();
    } ) );
    _propertyUpdateButton.onClick( MemberFunctionCallback<VortexCollectionCreator>( this, &VortexCollectionCreator::updateButton ) );
    _propertySaveButton.onClick( LambdaFunctionCallback( [this]
    {
        if( _outportVortexCollection.hasData() && _propertyFileDialog.get() != "" )
        {
            auto stream = std::ofstream( _propertyFileDialog.get(), std::ios::out | std::ios::binary );
            _outportVortexCollection.getData()->serialize( stream );
        }
    } ) );
}

void VortexCollectionCreator::updateButton()
{
    const auto* ensemble = _inportEnsemble.getData();
    if( !ensemble ) return;

    const auto& members = _propertySelectedMembers.get();
    const auto interval = _propertyTimestepInterval.get();

    auto collection = std::unique_ptr<VortexCollection>( new VortexCollection( members.size(), interval.y - interval.x + 1 ) );

    for( size_t i = 0; i < members.size(); ++i )
    {
        for( size_t j = 0; j <= interval.y - interval.x; ++j )
        {
            const auto timeStart = std::chrono::high_resolution_clock::now();
            //const auto& fields = ensemble->getMembers()[members[i]].getTimeSteps()[interval.x + j].getFieldNames();
            const auto& timeStep = ensemble->getMembers()[members[i]].getTimeSteps()[interval.x + j];

            const auto maskLock = VolumeRAMRepresentationLock( timeStep.getVolume( "SALT" ) );
            const auto mask = dynamic_cast<const VolumeRAM_Float*>( maskLock.operator->() );
            const auto dim = mask->getDimensions();

            // --- Channel Merger --- //
            const auto volumeULock = VolumeRAMRepresentationLock( timeStep.getVolume( "U" ) );
            const auto volumeVLock = VolumeRAMRepresentationLock( timeStep.getVolume( "V" ) );
            const auto volumeWLock = VolumeRAMRepresentationLock( timeStep.getVolume( "W" ) );

            const auto volumeU = dynamic_cast<const VolumeRAM_Float*>( *volumeULock );
            const auto volumeV = dynamic_cast<const VolumeRAM_Float*>( *volumeVLock );
            const auto volumeW = dynamic_cast<const VolumeRAM_Float*>( *volumeWLock );

            auto velocityRAM = std::unique_ptr<VolumeRAM_3xFloat>( new VolumeRAM_3xFloat( dim ) );
#pragma omp parallel for
            for( long x = 0; x < dim.x; ++x ) {
                for (long y = 0; y < dim.y; ++y) {
                    for (long z = 0; z < dim.z; ++z) {
                        velocityRAM->voxel(x, y, z) = tgt::dvec3(volumeU->voxel(x, y, z),
                                                                 volumeV->voxel(x, y, z),
                                                                 volumeW->voxel(x, y, z));
                    }
                }
            }

            const auto velocityVolume = std::unique_ptr<Volume>( new Volume( velocityRAM.release(), timeStep.getVolume( "U" ) ) );
            const auto velocityLock = VolumeRAMRepresentationLock( velocityVolume.get() );
            const auto velocity = dynamic_cast<const VolumeRAM_3xFloat*>( velocityLock.operator->() );

            const auto velocityFloatVolume = std::unique_ptr<Volume>( VolumeOperatorConvert().apply<tgt::dvec3>( velocityVolume.get() ) );
            const auto velocityFloatLock = VolumeRAMRepresentationLock( velocityFloatVolume.get() );
            const auto velocityFloat = dynamic_cast<const VolumeRAM_3xFloat*>( velocityFloatLock.operator->() );

            // --- Vortex Processor --- //
            auto jacobi = std::unique_ptr<VolumeRAM_Mat3Float>( new VolumeRAM_Mat3Float( dim ) );
            VortexProcessor::Process( *velocity, *mask, *jacobi, nullptr, nullptr, nullptr );

            // --- Acceleration Processor --- //
            auto acceleration = std::unique_ptr<VolumeRAM_3xFloat>(new VolumeRAM_3xFloat( dim ) );
            Acceleration::Process(*jacobi, *velocity, *acceleration);

            // --- Parallel Vectors --- //
            auto parallelVectors = std::unique_ptr<ParallelVectorSolutions>( new ParallelVectorSolutions() );
            ParallelVectors::Process( *velocityFloat, *acceleration, jacobi.get(), mask, *parallelVectors );
            std::cout << "parallel vectors solutions: " << parallelVectors->solutions.size() << std::endl;

            // --- Coreline Creator --- //
            auto corelines = std::vector<std::vector<tgt::vec3>>();
            CorelineCreator::Process( *parallelVectors, _propertyCorelineLength.get(), corelines );

            // --- Curl Processor --- //
            auto curl = std::unique_ptr<VolumeRAM_3xFloat>( new VolumeRAM_3xFloat( dim ) );
            CurlProcessor::Process( *jacobi, *curl );

            auto vortices = std::vector<Vortex>( corelines.size() );
            for( size_t i = 0; i < vortices.size(); ++i )
            {
                vortices[i].setCoreline( std::move( corelines[i] ) );
                RotationalDirectionProcessor::Process( *curl, vortices[i] );
            }
            collection->setVortices( i, j, std::move( vortices ) );
        }
    }

    for( size_t member = 0; member < collection->members(); ++member )
    {
        for( size_t timestep = 0; timestep < collection->timesteps() - 1; ++timestep )
        {
            const auto& vortices = collection->vortices( member, timestep );
            auto vorticesNextTimestep = collection->vortices( member, timestep + 1 );

            for( size_t i = 0; i < vortices.size(); ++i )
            {
                size_t index;
                VortexTracking::Process( vortices[i], vorticesNextTimestep, 6.5f, index );
                if( index != std::numeric_limits<size_t>::max() )
                {
                    collection->addMatch( VortexCollection::VortexID( member, timestep, i ), VortexCollection::VortexID( member, timestep + 1, index ) );
                }
            }
        }
    }

    _outportVortexCollection.setData( collection.release() );
}

}
