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

#include "vortexlistselector.h"

#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

namespace voreen {

VortexListSelector::VortexListSelector() : Processor(),
    _inportVortexCollection( Port::INPORT, "inport_vortex_collection", "Vortex Collection" ),
    _outportVortexList( Port::OUTPORT, "outport_vortex_list", "Vortex List" ),
    _outportGeometry( Port::OUTPORT, "outport_geometry", "Coreline List" ),
    _propertyMembers( "property_members", "Members", Processor::VALID ),
    _propertyTimesteps( "property_timesteps", "Timesteps", 0, 0, std::numeric_limits<int>::max(), 0, std::numeric_limits<int>::max(), Processor::VALID ),
    _propertyCorelineLength( "property_coreline_length", "Minimal Coreline Length", 40, 2, std::numeric_limits<int>::max(), Processor::VALID ),
    _Rotation("Rotation", "Direction of Rotation:")
{
    this->addPort( _inportVortexCollection );
    this->addPort( _outportVortexList );
    this->addPort( _outportGeometry );

    this->addProperty( _propertyMembers );
    this->addProperty( _propertyTimesteps );
    this->addProperty( _propertyCorelineLength );

    this->addProperty(_Rotation);
    _Rotation.reset();
    _Rotation.addOption("B", "Both", OPTION_B);
    _Rotation.addOption("CW", "Clockwise", OPTION_CW);
    _Rotation.addOption("CCW", "Counter-Clockwise", OPTION_CCW);
    _Rotation.selectByValue(OPTION_B);

    _inportVortexCollection.onNewData( LambdaFunctionCallback( [this]
    {
        _propertyMembers.blockCallbacks( true );
        _propertyMembers.reset();
        for( size_t i = 0; i < _inportVortexCollection.getData()->members(); ++i )
            _propertyMembers.addRow( std::to_string( i ) );
        _propertyMembers.blockCallbacks( false );

        _propertyTimesteps.blockCallbacks( true );
        _propertyTimesteps.setMinValue( 0 );
        _propertyTimesteps.setMaxValue( static_cast<int>( _inportVortexCollection.getData()->timesteps() ) - 1 );
        _propertyTimesteps.blockCallbacks( false );

        this->updatePropertyCorelineLength();
    } ) );

    _propertyMembers.onChange( LambdaFunctionCallback( [this]
    {
        if( !_inportVortexCollection.hasData() ) return;

        this->updatePropertyCorelineLength();
        this->invalidate();
    } ) );
    _propertyTimesteps.onChange( LambdaFunctionCallback( [this]
    {
        if( !_inportVortexCollection.hasData() ) return;

        this->updatePropertyCorelineLength();
        this->invalidate();
    } ) );
    _propertyCorelineLength.onChange( LambdaFunctionCallback( [this]
    {
        this->invalidate();
    } ) );
}

void VortexListSelector::Process( const VortexCollection& vortexCollection, const std::vector<int>& members, int firstTimestep, int lastTimestep, int minLength, RotationOptions rot, std::vector<Vortex>& outVortexList )
{
    outVortexList.clear();
    for( const auto member : members )
        for( int timestep = firstTimestep; timestep <= lastTimestep; ++timestep )
            for( const auto& vortex : vortexCollection.vortices( member, timestep ) )
                if( vortex.coreline().size() >= minLength )
                    switch(rot){
                    case OPTION_B:
                        outVortexList.push_back(vortex);
                        break;
                    case OPTION_CW:
                        if (vortex.getOrientation() == Vortex::Orientation::eClockwise) { outVortexList.push_back(vortex); }//falls clockwise dreht
                        break;
                    case OPTION_CCW:
                        if (vortex.getOrientation() == Vortex::Orientation::eCounterClockwise) { outVortexList.push_back(vortex); }//falls counterclockwise dreht
                        break;
                    default:
                        break;
                    }
    outVortexList.shrink_to_fit();
}

void VortexListSelector::process()
{
    const auto collection = _inportVortexCollection.getData();
    if( !collection )
    {
        _outportGeometry.clear();
        return;
    }

    auto outVortexList = std::unique_ptr<std::vector<Vortex>>( new std::vector<Vortex>());
    VortexListSelector::Process( *collection, _propertyMembers.get(), _propertyTimesteps.get().x, _propertyTimesteps.get().y, _propertyCorelineLength.get(), _Rotation.getValue(), *outVortexList );

    auto corelines = std::vector<std::vector<tgt::vec3>>();
    corelines.reserve(outVortexList->size());

    for( const auto& vortex : *outVortexList )
        corelines.push_back( vortex.coreline() );

    auto geometry = new PointSegmentListGeometryVec3();
    geometry->setData( std::move( corelines ) );
    _outportVortexList.setData( outVortexList.release() );
    _outportGeometry.setData( geometry );
}

void VortexListSelector::updatePropertyCorelineLength()
{
    auto minLength = std::numeric_limits<int>::max();
    auto maxLength = 0;
    for( const auto member : _propertyMembers.get() )
    {
        for( auto timestep = _propertyTimesteps.get().x; timestep <= _propertyTimesteps.get().y; ++timestep )
        {
            for( const auto& vortex : _inportVortexCollection.getData()->vortices( member, timestep ) )
            {
                minLength = std::min( minLength, static_cast<int>( vortex.coreline().size() ) );
                maxLength = std::max( maxLength, static_cast<int>( vortex.coreline().size() ) );
            }
        }
    }

    _propertyCorelineLength.blockCallbacks( true );
    _propertyCorelineLength.setMinValue( minLength );
    _propertyCorelineLength.setMaxValue( maxLength );
    _propertyCorelineLength.blockCallbacks( false );
}

}
