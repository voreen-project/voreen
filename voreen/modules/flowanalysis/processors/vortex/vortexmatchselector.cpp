/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "vortexmatchselector.h"

#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "modules/flowanalysis/processors/vortex/vortextracking.h"

namespace voreen {

VortexMatchSelector::VortexMatchSelector() : RenderProcessor(),
    _inportVortexCollection( Port::INPORT, "inport_vortex_collection", "Vortex Collection" ),
    _inportVortexCollectionMatch( Port::INPORT, "inport_vortex_collection_match", "Vortex Collection for Matching" ),
    _outportRender( Port::OUTPORT, "outport_render", "Merge Tree", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER ),
    _outportGeometryClockwise( Port::OUTPORT, "outport_geometry_clockwise", "Corelines Clockwise" ),
    _outportGeometryCounterClockwise( Port::OUTPORT, "outport_geometry_counter_clockwise", "Corelines Counter-Clockwise" ),

    _eventPropertyHover( new EventProperty<VortexMatchSelector>( "event_property_hover", "Hover Event Property", this, &VortexMatchSelector::onHoverEvent, tgt::MouseEvent::MOUSE_BUTTON_NONE, tgt::MouseEvent::MouseAction::MOTION, tgt::Event::MODIFIER_NONE, false, true ) ),
    _eventPropertyMouse( new EventProperty<VortexMatchSelector>( "event_property_mouse", "Mouse Event Property", this, &VortexMatchSelector::onMouseEvent, static_cast<tgt::MouseEvent::MouseButtons>( tgt::MouseEvent::MOUSE_BUTTON_LEFT | tgt::MouseEvent::MOUSE_BUTTON_RIGHT ) ) ),
    _eventPropertyWheel( new EventProperty<VortexMatchSelector>( "event_property_wheel", "Wheel Event Property", this, &VortexMatchSelector::onWheelEvent, tgt::MouseEvent::MOUSE_BUTTON_ALL, tgt::MouseEvent::MouseAction::WHEEL, tgt::Event::MODIFIER_NONE, false, true ) ),

    _propertyMaximumMatchDistance( "property_maximum_match_distance", "Maximum Match Distance", 5.0f, 0.0f, 32.0f, Processor::VALID ),
    _propertySelectedRuns( "property_selected_runs", "Runs", Processor::VALID ),
    _propertyMinimumGroupSize( "property_minimum_group_size", "Minimum Group Size", 1, 1, std::numeric_limits<int>::max(), Processor::VALID ),
    _propertyOrientation( "property_orientaton", "Orientation", Processor::VALID, true ),
    _propertyTimestepInterval( "property_timestep_interval", "Timesteps", 0, 0, std::numeric_limits<int>::max(), 0, std::numeric_limits<int>::max(), Processor::VALID ),
    _propertyCorelineLength( "property_coreline_length", "Minimum Coreline Points", 0, 0, std::numeric_limits<int>::max(), Processor::VALID ),
    _propertySelectAll( "property_select_all", "Select All", Processor::VALID ),
    _propertyClearSelection( "property_clear_selection", "Clear Selection", Processor::VALID )
{
    this->addPort( _inportVortexCollection );
    this->addPort( _inportVortexCollectionMatch );
    this->addPort( _outportRender );
    this->addPort( _outportGeometryClockwise );
    this->addPort( _outportGeometryCounterClockwise );

    this->addEventProperty( _eventPropertyHover.get() );
    this->addEventProperty( _eventPropertyMouse.get() );
    this->addEventProperty( _eventPropertyWheel.get() );

    this->addProperty( _propertyMaximumMatchDistance );
    this->addProperty( _propertySelectedRuns );
    this->addProperty( _propertyMinimumGroupSize );
    this->addProperty( _propertyOrientation );
    this->addProperty( _propertyTimestepInterval );
    this->addProperty( _propertyCorelineLength );
    this->addProperty( _propertySelectAll );
    this->addProperty( _propertyClearSelection );

    _propertyOrientation.setOptions( std::deque<Option<int>> {
        Option<int>( "Both", "Both", 0 ),
            Option<int>( "Clockwise", "Clockwise", static_cast<int>( Vortex::Orientation::eClockwise ) ),
            Option<int>( "Counter-Clockwise", "Counter-Clockwise", static_cast<int>( Vortex::Orientation::eCounterClockwise ) )
    } );

    _inportVortexCollection.onNewData( MemberFunctionCallback<VortexMatchSelector>( this, &VortexMatchSelector::onNewData ) );
    _inportVortexCollectionMatch.onChange( MemberFunctionCallback<VortexMatchSelector>( this, &VortexMatchSelector::updateProbabilities ) );

    _propertyMaximumMatchDistance.onChange( LambdaFunctionCallback( [this]
    {
        this->updateProbabilities();
        this->invalidate();
    } ) );
    _propertySelectedRuns.onChange( LambdaFunctionCallback( [this]
    {
        this->updateVortexInfos();
        this->invalidate();
    } ) );
    _propertyMinimumGroupSize.onChange( LambdaFunctionCallback( [this]
    {
        this->updateVortexInfos();
        this->invalidate();
    } ) );
    _propertyOrientation.onChange( LambdaFunctionCallback( [this]
    {
        this->updateVortexInfos();
        this->invalidate();
    } ) );
    _propertyTimestepInterval.onChange( LambdaFunctionCallback( [this]
    {
        this->updateGeometry();
        this->invalidate();
    } ) );
    _propertyCorelineLength.onChange( LambdaFunctionCallback( [this]
    {
        this->updateGeometry();
        this->invalidate();
    } ) );
    _propertySelectAll.onChange( LambdaFunctionCallback( [this]
    {
        if( !_inportVortexCollection.hasData() ) return;

        for( auto& vortexInfos : _vortexInfos ) for( auto& vortexInfo : vortexInfos )
            vortexInfo.selected = true;
        this->updateGeometry();
        this->invalidate();
    } ) );
    _propertyClearSelection.onChange( LambdaFunctionCallback( [this]
    {
        if( !_inportVortexCollection.hasData() ) return;

        for( auto& vortexInfos : _vortexInfos ) for( auto& vortexInfo : vortexInfos )
            vortexInfo.selected = false;
        this->updateGeometry();
        this->invalidate();
    } ) );
}

void VortexMatchSelector::process()
{
    const auto collection = _inportVortexCollection.getData();
    if( !collection )
    {
        _outportRender.clear();
        _outportGeometryClockwise.clear();
        _outportGeometryCounterClockwise.clear();
        return;
    }

    _outportRender.activateTarget();

    // --- Clear Background --- //
    glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

    // --- Draw Boundary --- //
    const auto renderSize = _outportRender.getSize();
    glViewport( 0, 0, renderSize.x, renderSize.y );
    glBegin( GL_LINES );
    glColor4f( 0.0f, 0.0f, 0.0, 1.0f );
    glVertex2f( -1.0f, 0.0f );
    glVertex2f( 1.0f, 0.0f );
    glEnd();

    // --- Draw Timestep Interval --- //
    const auto dim = renderSize.y / 2;
    glViewport( 0, dim, renderSize.x, dim );
    glBegin( GL_LINES );
    glColor4f( 0.0f, 0.0f, 0.0, 1.0f );

    const auto interval = _propertyTimestepInterval.get();
    auto x = this->treeToClip( tgt::vec2( 7.5f + interval.x * 15.0f, 0.0f ) ).x;
    glVertex2f( x, -1.0f );
    glVertex2f( x, 0.99f );

    x = this->treeToClip( tgt::vec2( 7.5f + ( interval.y + 1 ) * 15.0f, 0.0f ) ).x;
    glVertex2f( x, -1.0f );
    glVertex2f( x, 0.99f );
    glEnd();
    this->render( Visualization::eTree );

    glViewport( ( renderSize.x - dim ) / 2, 0, dim, dim );
    this->render( Visualization::eWorld );

    glClearColor( 0.0f, 0.0f, 0.0f, 0.0f );
    _outportRender.deactivateTarget();
}
void VortexMatchSelector::render( Visualization visualization )
{
    const auto collection = _inportVortexCollection.getData();
    const auto size =  _outportRender.getSize();

    // --- Draw Vortices --- //
    for( size_t timestep = 0; timestep < collection->timesteps(); ++timestep )
    {
        for( size_t run = 0; run < collection->runs(); ++run )
        {
            const auto& vortices = collection->vortices( run, timestep );
            const auto& vortexInfos = _vortexInfos[run * collection->timesteps() + timestep];

            for( size_t vortexIndex = 0; vortexIndex < vortices.size(); ++vortexIndex )
            {
                const auto& vortex = vortices[vortexIndex];
                const auto& vortexInfo = vortexInfos[vortexIndex];
                if( !vortexInfo.visible ) continue;

                float depth = 0.1f;
                float size = vortexInfo.predecessor ? 7.0f : 10.0f;
                tgt::vec4 color;

                if( visualization == Visualization::eWorld ) size -= 3.0f;
                if( vortex.getOrientation() == Vortex::Orientation::eCounterClockwise ) color = tgt::vec4( 0.0f, 0.0f, 1.0f, 1.0f );
                else if( vortex.getOrientation() == Vortex::Orientation::eClockwise ) color = tgt::vec4( 1.0f, 0.0f, 0.0f, 1.0f );
                else color = tgt::vec4( 0.2f, 0.2f, 0.2f, 1.0f );
                if( vortexInfo.highlighted || vortexInfo.selected ) color *= 0.6f;
                color.w = 0.25f + vortexInfo.probability * 0.75f;

                if( vortexInfo.highlighted ) depth = 0.0f, size += 3.0f, color = tgt::vec4( 0.0f, 1.0f, 0.0f, 1.0f );

                const auto pos = visualization == Visualization::eTree ?
                    this->treeToClip( vortexInfo.positionTree ) : this->worldToClip( vortexInfo.positionWorld );

                glPointSize( size );
                glColor4f( color.x, color.y, color.z, color.w );
                glBegin( GL_POINTS );
                glVertex3f( pos.x, pos.y, depth );
                glEnd();

            }
        }
    }

    // --- Draw Lines --- //
    glBegin( GL_LINES );
    glColor4f( 0.5f, 0.5f, 0.5f, 1.0f );
    for( size_t run = 0; run < collection->runs(); ++run )
    {
        for( size_t timestep = 0; timestep < collection->timesteps(); ++timestep )
        {
            const auto& vortices = collection->vortices( run, timestep );
            const auto& vortexInfos = _vortexInfos[run * collection->timesteps() + timestep];
            for( size_t vortexIndex = 0; vortexIndex < vortices.size(); ++vortexIndex )
            {
                const auto& vortexInfo = vortexInfos[vortexIndex];
                if( !vortexInfo.visible ) continue;

                const auto begin = visualization == Visualization::eTree ?
                    this->treeToClip( vortexInfo.positionTree ) : this->worldToClip( vortexInfo.positionWorld );
                const auto& matches = collection->matches( run, timestep, vortexIndex );
                for( const auto& match : matches )
                {
                    const auto& matchInfo = _vortexInfos[match.run * collection->timesteps() + match.timestep][match.index];
                    if( !matchInfo.visible ) continue;

                    const auto end = visualization == Visualization::eTree ?
                        this->treeToClip( matchInfo.positionTree ) : this->worldToClip( matchInfo.positionWorld );
                    glVertex3f( begin.x, begin.y, 0.2f );
                    glVertex3f( end.x, end.y, 0.2f );
                }
            }
        }
    }
    glEnd();
}

tgt::vec2 VortexMatchSelector::treeToClip( tgt::vec2 pos ) const
{
    const auto size = tgt::vec2( _outportRender.getSize().x, _outportRender.getSize().y / 2 );
    return ( _offsetTree + pos ) / size * 2.0f - 1.0f;
}
tgt::vec2 VortexMatchSelector::worldToClip( tgt::vec2 pos ) const
{
    return _offsetWorld + _zoomWorld * ( pos * 2.0f - 1.0f );
}

void VortexMatchSelector::onNewData()
{
    const auto collection = _inportVortexCollection.getData();

    // --- Initialize Vortex Infos --- //
    auto minCorelineLength = std::numeric_limits<size_t>::max();
    auto maxCorelineLength = std::numeric_limits<size_t>::min();
    auto positionWorldMin = tgt::vec2( std::numeric_limits<float>::max(), std::numeric_limits<float>::max() );
    auto positionWorldMax = tgt::vec2( std::numeric_limits<float>::min(), std::numeric_limits<float>::min() );

    _vortexInfos = std::vector<std::vector<VortexInfo>>( collection->runs() * collection->timesteps() );
    for( size_t run = 0; run < collection->runs(); ++run ) for( size_t timestep = 0; timestep < collection->timesteps(); ++timestep )
    {
        const auto& vortices = collection->vortices( run, timestep );
        auto& vortexInfos = _vortexInfos[run * collection->timesteps() + timestep];
        vortexInfos.resize( vortices.size() );

        for( size_t vortexIndex = 0; vortexIndex < vortices.size(); ++vortexIndex )
        {
            const auto& vortex = vortices[vortexIndex];
            auto& vortexInfo = vortexInfos[vortexIndex];

            minCorelineLength = std::min( minCorelineLength, vortex.coreline().size() );
            maxCorelineLength = std::max( maxCorelineLength, vortex.coreline().size() );

            // --- Calculate World Position --- //
            auto& positionWorld = vortexInfo.positionWorld = tgt::vec2::zero;
            for( const auto point : vortex.coreline() )
                positionWorld += point.xy();
            positionWorld /= static_cast<float>( vortex.coreline().size() );
            positionWorld.y = 1.0f - positionWorld.y;

            positionWorldMin.x = std::min( positionWorldMin.x, positionWorld.x );
            positionWorldMin.y = std::min( positionWorldMin.y, positionWorld.y );
            positionWorldMax.x = std::max( positionWorldMax.x, positionWorld.x );
            positionWorldMax.y = std::max( positionWorldMax.y, positionWorld.y );
        }
    }
    for( auto& vortexInfos : _vortexInfos ) for( auto& vortexInfo : vortexInfos )
        vortexInfo.positionWorld = ( vortexInfo.positionWorld - positionWorldMin ) / ( positionWorldMax - positionWorldMin );

    // --- Find Vortex Groups --- //
    int32_t currentGroup = 0;
    for( size_t run = 0; run < collection->runs(); ++run ) for( size_t timestep = 0; timestep < collection->timesteps(); ++timestep )
        for( size_t vortexIndex = 0; vortexIndex < collection->vortices( run, timestep ).size(); ++vortexIndex )
            this->fillGroup( currentGroup++, VortexCollection::VortexID( run, timestep, vortexIndex ) );
    _groupInfos = std::vector<GroupInfo>( currentGroup );


    // --- Update Properties --- //
    _propertySelectedRuns.blockCallbacks( true );
    _propertySelectedRuns.reset();
    for( size_t run = 0; run < collection->runs(); ++run )
        _propertySelectedRuns.addRow( std::to_string( run ) );
    _propertySelectedRuns.blockCallbacks( false );

    _propertyTimestepInterval.blockCallbacks( true );
    const auto timestepMax = std::max( 0, static_cast<int>( collection->timesteps() ) - 1 );
    _propertyTimestepInterval.setMaxValue( timestepMax );
    _propertyTimestepInterval.set( tgt::ivec2( 0, 0 ) );
    _propertyTimestepInterval.blockCallbacks( false );

    _propertyCorelineLength.blockCallbacks( true );
    _propertyCorelineLength.setMinValue( static_cast<int>( minCorelineLength ) );
    _propertyCorelineLength.setMaxValue( static_cast<int>( maxCorelineLength ) );
    _propertyCorelineLength.blockCallbacks( false );

    this->updateProbabilities();
    this->updateVortexInfos();
    this->invalidate();
}
void VortexMatchSelector::onHoverEvent( tgt::MouseEvent* event )
{
    const auto collection = _inportVortexCollection.getData();
    if( !collection ) return;

    const auto renderSize = tgt::vec2( _outportRender.getSize() );
    const auto cursor = tgt::vec2( event->x(), renderSize.y - event->y() );
    const auto dim = renderSize.y / 2;

    struct { VortexCollection::VortexID vortex; float distance; } current { VortexCollection::VortexID::Invalid, 7.5f };
    for( size_t run = 0; run < collection->runs(); ++run ) 	for( size_t timestep = 0; timestep < collection->timesteps(); ++timestep )
    {
        auto& vortexInfos = _vortexInfos[run * collection->timesteps() + timestep];

        for( size_t vortexIndex = 0; vortexIndex < vortexInfos.size(); ++vortexIndex )
        {
            auto& vortexInfo = vortexInfos[vortexIndex];
            if( !vortexInfo.visible ) continue;
            vortexInfo.highlighted = false;

            const auto pixelTree = tgt::vec2( vortexInfo.positionTree.x, vortexInfo.positionTree.y + _offsetTree.y + dim );
            auto distance = tgt::length( cursor - pixelTree );
            if( distance < current.distance ) current = { VortexCollection::VortexID( run, timestep, vortexIndex ), distance };

            // const auto pixelWorld = tgt::vec2( ( renderSize.x - dim ) / 2.0f + vortexInfo.positionWorld.x * dim, vortexInfo.positionWorld.y * dim );
            const auto pixelWorld = tgt::vec2( ( renderSize.x - dim ) / 2.0f, 0.0f ) + ( this->worldToClip( vortexInfo.positionWorld ) + 1.0f ) / 2.0f * tgt::vec2( dim, dim );
            distance = tgt::length( cursor - pixelWorld );
            if( distance < current.distance ) current = { VortexCollection::VortexID( run, timestep, vortexIndex ), distance };
        }
    }

    if( _hoveredVortex != current.vortex )
    {
        _hoveredVortex = current.vortex;
        if( _hoveredVortex != VortexCollection::VortexID::Invalid )
        {
            auto& vortex = _vortexInfos[_hoveredVortex.run * collection->timesteps() + _hoveredVortex.timestep][_hoveredVortex.index];
            vortex.highlighted = true;
            std::cout << "[VortexMatchProcessor] Probability of Hovered Vortex: " << vortex.probability << std::endl;
        }

        this->process();
    }
}
void VortexMatchSelector::onMouseEvent( tgt::MouseEvent* event )
{
    const auto collection = _inportVortexCollection.getData();
    if( !collection ) return;

    if( _hoveredVortex != VortexCollection::VortexID::Invalid )
    {
        if( event->action() == tgt::MouseEvent::MouseAction::RELEASED )
        {
            auto& vortexInfo = _vortexInfos[_hoveredVortex.run * collection->timesteps() + _hoveredVortex.timestep][_hoveredVortex.index];

            if( event->button() == tgt::MouseEvent::MouseButtons::MOUSE_BUTTON_LEFT )
                this->fillSelection( !vortexInfo.selected, _hoveredVortex );
            else if( event->button() == tgt::MouseEvent::MouseButtons::MOUSE_BUTTON_RIGHT )
                vortexInfo.selected = !vortexInfo.selected;

            this->updateGeometry();
            this->process();
        }
    }
}
void VortexMatchSelector::onWheelEvent( tgt::MouseEvent* event )
{
    if( !_inportVortexCollection.hasData() ) return;

    const auto renderSize = _outportRender.getSize();
    const auto cursor = tgt::vec2( event->x(), renderSize.y - event->y() );
    const auto dim = renderSize.y / 2;

    if( cursor.y >= dim )
    {
        if( event->button() == tgt::MouseEvent::MOUSE_WHEEL_DOWN ) _offsetTree.y += 20.0f;
        else if( event->button() == tgt::MouseEvent::MOUSE_WHEEL_UP ) _offsetTree.y -= 20.0f;
        _offsetTree.y = std::min( 0.0f, _offsetTree.y );
    }
    else
    {
        auto point = tgt::vec2( cursor.x - ( renderSize.x - dim ) / 2.0f, cursor.y ) / tgt::vec2( dim, dim );
        point = point * 2.0f - 1.0f;
        point = ( point - _offsetWorld ) / _zoomWorld;
        const auto before = _offsetWorld + _zoomWorld * point;

        if( event->button() == tgt::MouseEvent::MOUSE_WHEEL_DOWN ) _zoomWorld /= 1.05f;
        else if( event->button() == tgt::MouseEvent::MOUSE_WHEEL_UP ) _zoomWorld *= 1.05f;
        _zoomWorld = std::max( _zoomWorld, 0.95f );

        const auto after = _zoomWorld * point;
        _offsetWorld = ( _zoomWorld == 0.95f ) ? tgt::vec2::zero : before - after;
    }

    this->process();
}

void VortexMatchSelector::updateProbabilities()
{
    const auto collection = _inportVortexCollection.getData();
    const auto collectionMatch = _inportVortexCollectionMatch.getData();
    if( !collection || !collectionMatch || collection->timesteps() != collectionMatch->timesteps() || _vortexInfos.size() == 0 )
    {
        for( auto& vortexInfos : _vortexInfos ) for( auto& vortexInfo : vortexInfos )
            vortexInfo.probability = 1.0f;
        return;
    }

    auto count = 0;
#pragma omp parallel for shared( count )
    for( long timestep = 0; timestep < static_cast<long>( collection->timesteps() ); ++timestep )
    {
#pragma omp critical
        std::cout << "[VortexMatchSelector] Computing Probabilities: timesteps = " << ++count << std::endl;

        for( size_t run = 0; run < collection->runs(); ++run )
        {
            const auto& vortices = collection->vortices( run, timestep );
            auto& vortexInfos = _vortexInfos[run * collection->timesteps() + timestep];

            for( size_t vortexIndex = 0; vortexIndex < vortices.size(); ++vortexIndex )
            {
                const auto& vortex = vortices[vortexIndex];
                auto& vortexInfo = vortexInfos[vortexIndex];

                size_t count = 0;
                for( size_t runMatching = 0; runMatching < collectionMatch->runs(); ++runMatching )
                {
                    const auto& vorticesMatching = collectionMatch->vortices( runMatching, timestep );
                    size_t index;
                    VortexTracking::Process( vortex, vorticesMatching, _propertyMaximumMatchDistance.get(), index );
                    if( index != std::numeric_limits<size_t>::max() ) ++count;
                }

                vortexInfo.probability = static_cast<float>( count ) / collectionMatch->runs();
            }
        }
    }
}
void VortexMatchSelector::updateVortexInfos()
{
    const auto collection = _inportVortexCollection.getData();
    const auto orientation = _propertyOrientation.getValue();
    if( !collection ) return;

    // --- Update Vortex Visibility --- //
    const auto selectedRuns = [this, collection]
    {
        const auto& selectedRuns = _propertySelectedRuns.get();
        auto result = std::vector<bool>( collection->runs() );
        for( const auto run : selectedRuns ) result[run] = true;
        return result;
    }( );

    for( size_t run = 0; run < collection->runs(); ++run )
    {
#pragma omp parallel for
        for( long timestep = 0; timestep < static_cast<long>( collection->timesteps() ); ++timestep )
        {
            const auto& vortices = collection->vortices( run, timestep );
            auto& vortexInfos = _vortexInfos[run * collection->timesteps() + timestep];

            for( size_t vortexIndex = 0; vortexIndex < vortices.size(); ++vortexIndex )
            {
                const auto& vortex = vortices[vortexIndex];
                auto& vortexInfo = vortexInfos[vortexIndex];
                vortexInfo.visible = !orientation || orientation == static_cast<int>( vortex.getOrientation() );
                vortexInfo.visible &= selectedRuns[run];

                // --- Check for Predecessor --- //
                vortexInfo.predecessor = false;
                for( const auto& match : collection->matches( run, timestep, vortexIndex ) ) if( match.run == run && match.timestep == timestep - 1 )
                {
                    const auto& matchInfo = _vortexInfos[match.run * collection->timesteps() + match.timestep][match.index];
                    vortexInfo.predecessor = true;
                }
            }
        }
    }

    // --- Update Group Infos --- //
#pragma omp parallel for
    for( long i = 0; i < static_cast<long>( _groupInfos.size() ); ++i ) _groupInfos[i] = GroupInfo {};
    for( size_t timestep = 0; timestep < collection->timesteps(); ++timestep )
    {
        auto groupWidths = std::vector<uint32_t>( _groupInfos.size() );
        for( size_t run = 0; run < collection->runs(); ++run )
        {
            const auto& vortexInfos = _vortexInfos[run * collection->timesteps() + timestep];

            for( size_t vortexIndex = 0; vortexIndex < vortexInfos.size(); ++vortexIndex )
            {
                const auto& vortexInfo = vortexInfos[vortexIndex];
                if( vortexInfo.visible )
                {
                    ++_groupInfos[vortexInfo.group].members;
                    ++groupWidths[vortexInfo.group];
                }
            }
        }

        for( size_t i = 0; i < _groupInfos.size(); ++i )
            _groupInfos[i].width = std::max( _groupInfos[i].width, groupWidths[i] );
    }

    // --- Update Vortex Visibility --- //
    for( size_t run = 0; run < collection->runs(); ++run )
    {
#pragma omp parallel for
        for( long timestep = 0; timestep < static_cast<long>( collection->timesteps() ); ++timestep )
        {
            auto& vortexInfos = _vortexInfos[run * collection->timesteps() + timestep];
            for( size_t vortexIndex = 0; vortexIndex < vortexInfos.size(); ++vortexIndex )
            {
                auto& vortexInfo = vortexInfos[vortexIndex];
                vortexInfo.visible &= _groupInfos[vortexInfo.group].members >= _propertyMinimumGroupSize.get();
            }
        }
    }

    // --- Update Properties --- //
    uint32_t maximumMembers = 1;
    for( const auto& groupInfo : _groupInfos ) maximumMembers = std::max( maximumMembers, groupInfo.members );
    _propertyMinimumGroupSize.blockCallbacks( true );
    _propertyMinimumGroupSize.setMaxValue( maximumMembers );
    _propertyMinimumGroupSize.blockCallbacks( false );

    // --- Update Vortex Tree Positions --- //
    const auto padding = tgt::vec2( 15.0f, 20.0f );
    const auto step = tgt::vec2( 15.0f, 15.0f );

    const auto groupHeights = [this, padding, step]
    {
        auto heights = std::vector<float>( _groupInfos.size() );
        for( size_t i = 0; i < heights.size(); ++i )
        {
            if( i == 0 ) heights[i] = padding.y;
            else
            {
                const auto previous = _groupInfos[i - 1].width * ( _groupInfos[i - 1].members >= _propertyMinimumGroupSize.get() );
                heights[i] = heights[i - 1] + previous * step.y;
            }
        }
        return heights;
    }( );

#pragma omp parallel for
    for( long timestep = 0; timestep < static_cast<long>( collection->timesteps() ); ++timestep )
    {
        auto heights = groupHeights;
        for( size_t run = 0; run < collection->runs(); ++run )
        {
            auto& vortexInfos = _vortexInfos[run * collection->timesteps() + timestep];
            for( size_t vortexIndex = 0; vortexIndex < vortexInfos.size(); ++vortexIndex )
            {
                auto& vortexInfo = vortexInfos[vortexIndex];
                if( !vortexInfo.visible ) continue;

                vortexInfo.positionTree.x = padding.x + timestep * step.x;

                auto& height = heights[vortexInfo.group];
                vortexInfo.positionTree.y = heights[vortexInfo.group];
                height += step.y;
            }
        }
    }

    this->updateGeometry();
}
void VortexMatchSelector::updateGeometry()
{
    const auto collection = _inportVortexCollection.getData();
    const auto interval = _propertyTimestepInterval.get();
    if( !collection ) return;

    // --- Gather Corelines of Selected Vortices --- //
    auto segmentsClockwise = std::vector<std::vector<tgt::vec3>>();
    auto segmentsCounterClockwise = std::vector<std::vector<tgt::vec3>>();
    for( size_t run = 0; run < collection->runs(); ++run )
    {
        for( size_t timestep = interval.x; timestep <= interval.y; ++timestep )
        {
            const auto& vortices = collection->vortices( run, timestep );
            const auto& vortexInfos = _vortexInfos[run * collection->timesteps() + timestep];
            for( size_t vortexIndex = 0; vortexIndex < vortices.size(); ++vortexIndex )
            {
                const auto& vortex = vortices[vortexIndex];
                const auto& vortexInfo = vortexInfos[vortexIndex];

                if( vortexInfo.visible && vortexInfo.selected && vortex.coreline().size() >= _propertyCorelineLength.get() )
                {
                    if( vortex.getOrientation() == Vortex::Orientation::eClockwise )
                        segmentsClockwise.push_back( vortex.coreline() );
                    else if( vortex.getOrientation() == Vortex::Orientation::eCounterClockwise )
                        segmentsCounterClockwise.push_back( vortex.coreline() );
                }
            }
        }
    }

    segmentsClockwise.shrink_to_fit();
    auto geometryClockwise = std::unique_ptr<PointSegmentListGeometryVec3>( new PointSegmentListGeometryVec3() );
    geometryClockwise->setData( std::move( segmentsClockwise ) );
    _outportGeometryClockwise.setData( geometryClockwise.release() );

    segmentsCounterClockwise.shrink_to_fit();
    auto geometryCounterClockwise = std::unique_ptr<PointSegmentListGeometryVec3>( new PointSegmentListGeometryVec3() );
    geometryCounterClockwise->setData( std::move( segmentsCounterClockwise ) );
    _outportGeometryCounterClockwise.setData( geometryCounterClockwise.release() );
}

void VortexMatchSelector::fillGroup( const int32_t group, const VortexCollection::VortexID& id )
{
    const auto collection = _inportVortexCollection.getData();
    auto& vortexInfo = _vortexInfos[id.run * collection->timesteps() + id.timestep][id.index];

    if( vortexInfo.group == -1 )
    {
        vortexInfo.group = group;
        for( const auto& match : collection->matches( id ) )
            this->fillGroup( group, match );
    }
}
void VortexMatchSelector::fillSelection( const bool selected, const VortexCollection::VortexID& id )
{
    const auto collection = _inportVortexCollection.getData();
    auto& vortexInfo = _vortexInfos[id.run * collection->timesteps() + id.timestep][id.index];

    vortexInfo.selected = selected;
    for( const auto& match : collection->matches( id ) ) if( match.timestep == id.timestep + 1 )
        this->fillSelection( selected, match );
}

} // namespace voreen
