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

#include "particlerenderer.h"

namespace voreen {

ParticleRenderer::ParticleRenderer() : RenderProcessor(),
    _inportEnsemble( Port::INPORT, "inport_ensemble", "Ensemble" ),
    _outportImage( Port::OUTPORT, "outport_image", "Image", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER ),

    _propertySelectedMember( "property_selected_member", "Member", Processor::VALID, true ),
    _propertyFlowComponentX( "property_flow_component_x", "Flow Component X", Processor::VALID, true ),
    _propertyFlowComponentY( "property_flow_component_y", "Flow Component Y", Processor::VALID, true ),
    _propertyFlowComponentZ( "property_flow_component_z", "Flow Component Z", Processor::VALID, true ),
    _propertySelectedField( "property_selected_field", "Field", Processor::VALID, true ),
    _propertySeedPoints( "property_seed_points", "Seed Points", 50000, 1, 1000000, Processor::VALID ),
    _propertyUpdateMesh( "property_update_mesh", "Update", Processor::VALID ),

    _propertyTransferFunction( "property_transfer_function", "Transfer Function", Processor::INVALID_RESULT ),
    _propertyResetDomain( "property_reset_domain", "Reset Domain" ),
    _propertyLineWidth( "property_line_width", "Line Width", 1.0f, 1.0f, 4.0f ),
    _propertyTimesteps( "property_timesteps", "Timesteps", 0, 0, std::numeric_limits<int>::max() ),
    _propertyShader( "property_shader", "Shader", "particlerenderer.frag", "particlerenderer.vert", "" ),
    _propertyCamera( "property_camera", "Camera", tgt::Camera(), true, true, 500.0f, Processor::INVALID_RESULT )
{
    this->addPort( _inportEnsemble );
    this->addPort( _outportImage );

    this->addProperty( _propertySelectedMember );
    this->addProperty( _propertyFlowComponentX );
    this->addProperty( _propertyFlowComponentY );
    this->addProperty( _propertyFlowComponentZ );
    this->addProperty( _propertySelectedField );
    this->addProperty( _propertySeedPoints );
    this->addProperty( _propertyUpdateMesh );

    this->addProperty( _propertyTransferFunction );
    this->addProperty( _propertyResetDomain );
    this->addProperty( _propertyLineWidth );
    this->addProperty( _propertyTimesteps );
    this->addProperty( _propertyShader );
    this->addProperty( _propertyCamera );

    _interactionHandlerCamera.reset( new CameraInteractionHandler( "interaction_handler_camera", "Camera Interaction", &_propertyCamera ) );
    this->addInteractionHandler( _interactionHandlerCamera.get() );

    _inportEnsemble.onNewData( LambdaFunctionCallback( [this]
    {
        const auto ensemble = _inportEnsemble.getData();
        const auto& members = ensemble->getMembers();
        const auto& fields = ensemble->getCommonFieldNames();

        auto options = std::deque<Option<int>>();
        for( auto i = 0; i < members.size(); ++i )
        {
            const auto& member = members[i];
            options.push_back( Option<int>( member.getName(), member.getName(), i, tgt::vec4( member.getColor(), 1.0f ) ) );
        }
        _propertySelectedMember.setOptions( std::move( options ) );

        options.clear();
        for( auto i = 0; i < fields.size(); ++i )
        {
            const auto& field = fields[i];
            options.push_back( Option<int>( field, field, i ) );
        }
        _propertyFlowComponentX.setOptions( options );
        _propertyFlowComponentY.setOptions( options );
        _propertyFlowComponentZ.setOptions( options );
        _propertySelectedField.setOptions( std::move( options ) );

        _propertySelectedMember.invalidate();
        _propertyFlowComponentX.invalidate();
        _propertyFlowComponentY.invalidate();
        _propertyFlowComponentZ.invalidate();
        _propertySelectedField.invalidate();

        _propertyTimesteps.setMaxValue( std::max( 0, static_cast<int>( ensemble->getMinNumTimeSteps() ) - 1 ) );
    } ) );
    _propertyUpdateMesh.onChange( MemberFunctionCallback<ParticleRenderer>( this, &ParticleRenderer::updateMesh ) );
    _propertyResetDomain.onClick( LambdaFunctionCallback( [this]
    {
        _propertyTransferFunction.get()->setDomain( _domain );
        _propertyTransferFunction.invalidate();
    } ) );
}

Processor* ParticleRenderer::create() const
{
    return new ParticleRenderer();
}
std::string ParticleRenderer::getCategory() const
{
    return "Flow Visualization";
}
std::string ParticleRenderer::getClassName() const
{
    return "ParticleRenderer";
}

void ParticleRenderer::initialize()
{
    RenderProcessor::initialize();

    // --- Build Shader --- //
    auto header = RenderProcessor::generateHeader( &tgt::GpuCapabilities::GlVersion::SHADER_VERSION_330 );
    header += _propertyTransferFunction.get()->getShaderDefines();
    _propertyShader.setHeader( header );
    _propertyShader.rebuild();
}
void ParticleRenderer::deinitialize()
{
    RenderProcessor::deinitialize();
}
void ParticleRenderer::process()
{
    _outportImage.activateTarget();
    _outportImage.clearTarget();

    if( auto shader = _propertyShader.getShader() )
    {
        shader->activate();
        this->setGlobalShaderParameters( shader, &_propertyCamera.get(), _outportImage.getSize() );

        auto textureUnit = tgt::TextureUnit();
        textureUnit.activate();
        _propertyTransferFunction.get()->getTexture()->bind();
        _propertyTransferFunction.get()->setUniform( shader, "transFunParam_", "transFuncTex_", textureUnit.getUnitNumber() );

        shader->setUniform( "domain", _propertyTransferFunction.get()->getDomain() );
        shader->setUniform( "timesteps", _propertyTimesteps.get() );

        glEnable( GL_BLEND );
        glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
        if( _propertyTransferFunction.get()->getAlphaMode() == TransFuncBase::TF_USE_ALPHA )
            glDisable( GL_DEPTH_TEST );

        glEnable( GL_LINE_SMOOTH );
        glLineWidth( _propertyLineWidth.get() );
        if( _mesh ) _mesh->render();
        glDisable( GL_LINE_SMOOTH );

        glBlendEquation( GL_FUNC_ADD );
        glBlendFunc( GL_ONE, GL_ZERO );
        glDisable( GL_BLEND );
        glEnable( GL_DEPTH_TEST );
        shader->deactivate();
    }

    _outportImage.deactivateTarget();
}

void ParticleRenderer::updateMesh()
{
    const auto ensemble = _inportEnsemble.getData();
    if( !ensemble ) return;

    const auto& member = ensemble->getMembers()[_propertySelectedMember.getValue()];
    const auto& timesteps = member.getTimeSteps();
    if( timesteps.size() < 2 ) return;

    const auto& flowComponentNameX = _propertyFlowComponentX.get();
    const auto& flowComponentNameY = _propertyFlowComponentY.get();
    const auto& flowComponentNameZ = _propertyFlowComponentZ.get();
    const auto& fieldName = _propertySelectedField.get();
    const auto seedPoints = _propertySeedPoints.get();

    const auto volumeMask = VolumeRAMRepresentationLock( timesteps.front().getVolume( fieldName ) );
    const auto dim = volumeMask->getDimensions();

    // --- Initialize Seed Points --- //
    auto lines = std::vector<std::vector<VertexColor>>( seedPoints );
    auto finished = std::vector<bool>( lines.size(), false );
    {
        auto voxels = std::vector<tgt::vec3>();
        voxels.reserve( volumeMask->getNumVoxels() );
        for( size_t x = 0; x < dim.x; ++x ) for( size_t y = 0; y < dim.y; ++y ) for( size_t z = 0; z < dim.z; ++z )
        {
            if( volumeMask->getVoxelNormalized( x, y, z ) > 0.0 )
                voxels.push_back( tgt::vec3( x, y, z ) + tgt::vec3( 0.5f, 0.5f, 0.5f ) );
        }
        std::random_shuffle( voxels.begin(), voxels.end() );
        if( voxels.size() < lines.size() ) lines.resize( voxels.size() );

#pragma omp parallel for
        for( long i = 0; i < static_cast<long>( lines.size() ); ++i )
            lines[i] = { VertexColor( voxels[i], tgt::vec4::zero ) };
    }

    _domain = tgt::vec2( std::numeric_limits<float>::max(), std::numeric_limits<float>::min() );
    for( size_t timestep = 0; true; ++timestep )
    {
        std::cout << "[ParticleRenderer] Updating Mesh: timestep = " << timestep << std::endl;
        const auto volumeValue = VolumeRAMRepresentationLock( timesteps[timestep].getVolume( fieldName ) );

#pragma omp parallel for
        for( long lineIndex = 0; lineIndex < static_cast<long>( lines.size() ); ++lineIndex )
        {
            if( finished[lineIndex] ) continue;

            auto& current = lines[lineIndex].back();
            const auto value = volumeValue->getVoxelNormalized( current.pos_ );

            lines[lineIndex].back().color_ = tgt::vec4( value, static_cast<float>( timestep ), 0.0f, 0.0f );
            _domain.x = std::min( _domain.x, value );
            _domain.y = std::max( _domain.y, value );
        }
        if( timestep ) const_cast<VolumeBase*>( timesteps[timestep].getVolume( fieldName ) )->removeRepresentation<VolumeRAM>();

        if( timestep < timesteps.size() - 2 )
        {
            const auto volumeX = VolumeRAMRepresentationLock( timesteps[timestep].getVolume( flowComponentNameX ) );
            const auto volumeY = VolumeRAMRepresentationLock( timesteps[timestep].getVolume( flowComponentNameY ) );
            const auto volumeZ = VolumeRAMRepresentationLock( timesteps[timestep].getVolume( flowComponentNameZ ) );

            const auto volumeXNext = VolumeRAMRepresentationLock( timesteps[timestep + 1].getVolume( flowComponentNameX ) );
            const auto volumeYNext = VolumeRAMRepresentationLock( timesteps[timestep + 1].getVolume( flowComponentNameY ) );
            const auto volumeZNext = VolumeRAMRepresentationLock( timesteps[timestep + 1].getVolume( flowComponentNameZ ) );

            const auto getVelocityRK4 = [&] ( tgt::vec3 p, float t )
            {
                const auto sample = [&] ( tgt::vec3 p, bool t )
                {
                    if( p.x >= 0.0 && p.x < dim.x && p.y >= 0.0 && p.y < dim.y && p.z >= 0.0 && p.z < dim.z && volumeMask->getVoxelNormalized( p ) > 0.0 )
                    {
                        return t ?
                            tgt::vec3( volumeXNext->getVoxelNormalized( p ), volumeYNext->getVoxelNormalized( p ), volumeZNext->getVoxelNormalized( p ) )
                            : tgt::vec3( volumeX->getVoxelNormalized( p ), volumeY->getVoxelNormalized( p ), volumeZ->getVoxelNormalized( p ) );
                    }
                    return tgt::vec3( 0.0, 0.0, 0.0 );
                };
                const auto trilinearInterpolation = [&] ( tgt::vec3 p, bool t )
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

                const auto first = trilinearInterpolation( p, false );
                const auto second = trilinearInterpolation( p, true );

                return ( 1.0f - t ) * first + t * second;
            };

#pragma omp parallel for
            for( long lineIndex = 0; lineIndex < static_cast<long>( lines.size() ); ++lineIndex )
            {
                if( finished[lineIndex] ) continue;

                auto& line = lines[lineIndex];
                auto& current = line.back();
                auto next = tgt::vec3();

                const auto k1 = getVelocityRK4( current.pos_, 0.0f );
                const auto k2 = getVelocityRK4( current.pos_ + 0.5f * k1, 0.5f );
                const auto k3 = getVelocityRK4( current.pos_ + 0.5f * k2, 0.5f );
                const auto k4 = getVelocityRK4( current.pos_ + k3, 1.0f );
                next = current.pos_ + ( k1 + 2.0f * k2 + 2.0f * k3 + k4 ) / 6.0f;

                if( next.x < 0.0f || next.x >= dim.x || next.y < 0.0f || next.y >= dim.y || next.z < 0.0f || next.z >= dim.z || volumeMask->getVoxelNormalized( next ) <= 0.0 )
                {
                    finished[lineIndex] = true;
                    continue;
                }

                line.push_back( VertexColor( next, tgt::vec4::zero ) );
            }
        }
        else break;
    }

    // --- Build Mesh --- //
    auto mesh = std::unique_ptr<GlMeshGeometryUInt32Color>( new GlMeshGeometryUInt32Color() );
    mesh->setPrimitiveType( GL_LINE_STRIP );
    mesh->enablePrimitiveRestart();

    for( const auto& line : lines )
    {
        for( const auto point : line )
        {
            mesh->addIndex( static_cast<uint32_t>( mesh->getVertices().size() ) );
            mesh->addVertex( point );
        }
        mesh->addIndex( mesh->getPrimitiveRestartIndex() );
    }

    _propertyTransferFunction.get()->setDomain( _domain );
    _propertyTransferFunction.invalidate();

    _mesh.reset( mesh.release() );
}

}