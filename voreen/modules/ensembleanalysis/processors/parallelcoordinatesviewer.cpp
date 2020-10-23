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

#include "parallelcoordinatesviewer.h"
#include "voreen/core/properties/fontproperty.h"

namespace util {

tgt::vec2 screenToPoint( tgt::ivec2 screen, tgt::ivec2 viewport ) {
    const auto x = static_cast<float>( screen.x ) / viewport.x * 2.0f - 1.0f;
    const auto y = -1.0f * ( static_cast<float>( screen.y ) / viewport.y * 2.0f - 1.0f );
    return tgt::vec2( x, y );
}
tgt::ivec2 pointToScreen( tgt::vec2 point, tgt::ivec2 viewport ) {
    const auto x = static_cast<int>( ( point.x + 1.0f ) / 2.0f * viewport.x );
    const auto y = static_cast<int>( ( point.y + 1.0f ) / 2.0f * viewport.y );
    return tgt::ivec2( x, y );
}

}

namespace voreen {

const float ParallelCoordinatesViewer::X_LIMIT = 0.95f;
const float ParallelCoordinatesViewer::Y_LIMIT = 0.9f;

ParallelCoordinatesViewer::ParallelCoordinatesViewer()
    : RenderProcessor()
    , _axesport( Port::INPORT, "port_axes", "Parallel Coordinates Axes" )
    , _renderport( Port::OUTPORT, "port_render", "RenderPort", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , _eventPropertyHover( new EventProperty<ParallelCoordinatesViewer>( "event_property_hover", "Hover Event Property", this,
        &ParallelCoordinatesViewer::onHoverEvent, tgt::MouseEvent::MOUSE_BUTTON_NONE, tgt::MouseEvent::MouseAction::MOTION,
        tgt::Event::MODIFIER_NONE, false, true ) )
    , _eventPropertyMouse( new EventProperty<ParallelCoordinatesViewer>( "event_property_mouse",
        "Mouse Event Property", this, &ParallelCoordinatesViewer::onMouseEvent,
        static_cast<tgt::MouseEvent::MouseButtons>( tgt::MouseEvent::MOUSE_BUTTON_LEFT | tgt::MouseEvent::MOUSE_BUTTON_RIGHT ),
        tgt::MouseEvent::ACTION_ALL, tgt::Event::MODIFIER_NONE, false, true ) )
    , _propertySelectedMember( "property_selected_member", "Member", Processor::VALID, true )
    , _propertyVisualizationMode( "property_visualization_mode", "Visualization", Processor::VALID, true )
    , _propertySelectedTimestep( "property_selected_timestep", "Timestep", 0, 0, std::numeric_limits<int>::max(), Processor::VALID )
    , _propertySelectedField( "property_selected_field", "Field", Processor::VALID, true )
    , _propertyFieldInterval( "property_field_interval", "Interval", 0.0f, 0.0f, 0.0f, std::numeric_limits<float>::min(), std::numeric_limits<float>::max(), Processor::VALID )
    , _propertyDensityBlending( "property_density_blending", "Blend Densities", false, Processor::VALID )
    , _propertyDensityVisibleSamples( "property_density_visible_samples", "Unselected Samples for Maximum Density", 1, 1, std::numeric_limits<int>::max(), Processor::VALID )
    , _propertyDensitySelectedSamples( "property_density_selected_samples", "Selected Samples for Maximum Density", 1, 1, std::numeric_limits<int>::max(), Processor::VALID )
    , _propertyTransFuncField { // These have to be initialized here or the program crashes on shutdown
    IntOptionProperty( "property_trans_func_field0", "Transfer Function Field", Processor::VALID, true ),
    IntOptionProperty( "property_trans_func_field1", "Transfer Function Field", Processor::VALID, true ),
    IntOptionProperty( "property_trans_func_field2", "Transfer Function Field", Processor::VALID, true ),
    IntOptionProperty( "property_trans_func_field3", "Transfer Function Field", Processor::VALID, true ) }
    , _propertySelectedSamples( "property_selected_samples", "Selected Samples", std::vector<int>(), Processor::VALID )
    , _propertySections( "property_section", "Sections", Processor::VALID )
    , _shaderProgram(0)
    , _vertexArray(0)
    , _indexBuffer(0)
{
    // --- Initialize Ports --- //
    this->addPort( _axesport );
    this->addPort( _renderport );

    // --- Initialize Event Properties --- //
    this->addEventProperty( _eventPropertyHover.get() );
    this->addEventProperty( _eventPropertyMouse.get() );

    // --- Add Properties --- //
    this->addProperty( _propertySelectedMember );
    this->addProperty( _propertyVisualizationMode );

    this->addProperty( _propertySelectedTimestep );
    this->addProperty( _propertySelectedField );
    this->addProperty( _propertyFieldInterval );

    this->addProperty( _propertyDensityBlending );
    this->addProperty( _propertyDensityVisibleSamples );
    this->addProperty( _propertyDensitySelectedSamples );

    for( size_t i = 0; i < _propertyTransFuncField.size(); ++i )
    {
        // _propertyTransFuncField[i] = IntOptionProperty( "property_trans_func_field" + std::to_string( i ), "Transfer Function Field", Processor::VALID );
        this->addProperty( _propertyTransFuncField[i] );

        _propertyTransFunc[i] = TransFunc1DKeysProperty( "property_trans_func" + std::to_string( i ), "Transfer Function", Processor::VALID );
        this->addProperty( _propertyTransFunc[i] );

        _propertyTransFuncIntern[i] = TransFunc1DKeysProperty( "property_trans_func_intern" + std::to_string( i ), "Transfer Function (Linking)", Processor::VALID );
        this->addProperty( _propertyTransFuncIntern[i] );
    }

    this->addProperty( _propertySelectedSamples );
    this->addProperty( _propertySections );

    // --- Initialize Properties --- //
    _propertyVisualizationMode.setOptions( std::deque<Option<int>> {
        Option<int>( "fields", "Single Time Step", 0 ),
        Option<int>( "timesteps", "All Time Steps", 1 )
    });
    for( size_t i = 0; i < _propertyTransFuncIntern.size(); ++i )
        _propertyTransFuncIntern[i].setVisibleFlag( false );
    _propertySelectedSamples.setVisibleFlag( false );
    _propertySections.setVisibleFlag( false );

    // --- Initialize Callbacks --- //
    _axesport.onNewData( MemberFunctionCallback<ParallelCoordinatesViewer>( this, &ParallelCoordinatesViewer::onNewInportData ) );
    _propertySelectedMember.onChange( LambdaFunctionCallback( [this]
    {
        if( !_axesport.hasData() ) return;
        this->updateSampleStates();
        this->invalidate();
    } ) );
    _propertyVisualizationMode.onChange( LambdaFunctionCallback( [this]
    {
        if( !_axesport.hasData() ) return;
        this->onNewInportData();
        this->invalidate();
    } ) );
    _propertySelectedTimestep.onChange( LambdaFunctionCallback( [this]
    {
        if( !_axesport.hasData() ) return;
        this->updateSampleStates();
        this->invalidate();
    } ) );
    _propertySelectedField.onChange( LambdaFunctionCallback( [this]
    {
        if( !_axesport.hasData() ) return;
        const auto axes = _axesport.getData();
        const auto axisIndex = _propertySelectedField.getValue();
        const auto range = axes->getRange( axisIndex );

        _propertyFieldInterval.blockCallbacks( true );
        _propertyFieldInterval.setMinValue( range.first );
        _propertyFieldInterval.setMaxValue( range.second );
        _propertyFieldInterval.blockCallbacks( false );

        if( _propertyVisualizationMode.getValue() )
        {
            _propertyFieldInterval.blockCallbacks( true );
            _propertyFieldInterval.set( tgt::vec2( range.first, range.second ) );
            _propertyFieldInterval.blockCallbacks( false );

            // Update uniform buffer
            for( auto& uniform : _uniformBufferVec )
            {
                uniform.y = range.first;
                uniform.z = range.second;
            }
            glUseProgram( _shaderProgram );
            glUniform3fv( glGetUniformLocation( _shaderProgram, "axes" ), static_cast<GLsizei>( _uniformBufferVec.size() ), reinterpret_cast<const GLfloat*>( _uniformBufferVec.data() ) );
            glUseProgram( 0 );

            for( auto& sections : _sections ) sections.clear();
            this->updateSampleStates();

            this->invalidate();
        }
        else
        {
            const auto uniform = _uniformBufferVec[axisIndex];
            _propertyFieldInterval.blockCallbacks( true );
            _propertyFieldInterval.set( tgt::vec2( uniform.y, uniform.z ) );
            _propertyFieldInterval.blockCallbacks( false );
        }


    } ) );
    _propertyFieldInterval.onChange( LambdaFunctionCallback( [this]
    {
        if( !_axesport.hasData() ) return;

        const auto interval = _propertyFieldInterval.get();
        if( _propertyVisualizationMode.getValue() )
        {
            for( auto& uniform : _uniformBufferVec )
            {
                uniform.y = interval.x;
                uniform.z = interval.y;
            }
        }
        else
        {
            const auto axisIndex = _propertySelectedField.getValue();
            auto& uniform = _uniformBufferVec[axisIndex];
            uniform.y = interval.x;
            uniform.z = interval.y;

            for( size_t i = 0; i < _propertyTransFuncField.size(); ++i )
                if( _propertyTransFuncField[i].getValue() == axisIndex )
                    this->updateTransferFunction( i );
        }

        // Update uniform buffer
        glUseProgram( _shaderProgram );
        glUniform3fv( glGetUniformLocation( _shaderProgram, "axes" ), static_cast<GLsizei>( _uniformBufferVec.size() ), reinterpret_cast<const GLfloat*>( _uniformBufferVec.data() ) );
        glUseProgram( 0 );

        this->updateSampleStates();
        this->invalidate();
    } ) );
    _propertyDensityBlending.onChange( LambdaFunctionCallback( [this] { if( _axesport.hasData() ) this->invalidate(); } ) );
    _propertyDensityVisibleSamples.onChange( LambdaFunctionCallback( [this] { if( _axesport.hasData() ) this->invalidate(); } ) );
    _propertyDensitySelectedSamples.onChange( LambdaFunctionCallback( [this] { if( _axesport.hasData() ) this->invalidate(); } ) );
    for( size_t i = 0; i < _propertyTransFuncField.size(); ++i )
    {
        _propertyTransFuncField[i].onChange( LambdaFunctionCallback( [this, i]
        {
            if( !_axesport.hasData() ) return;
            if( _propertyTransFuncField[i].getValue() >= 0 )
            {
                this->updateTransferFunction( i );
            }
            else
            {
                auto begin = new TransFuncMappingKey( 0.0f, tgt::col4( 0, 0, 0, 0 ) );
                auto end = new TransFuncMappingKey( 1.0f, tgt::col4( 0, 0, 0, 0 ) );
                _propertyTransFuncIntern[i].get()->setDomain( 0.0f, 1.0f );
                _propertyTransFuncIntern[i].get()->setKeys( { begin, end } );
                _propertyTransFuncIntern[i].invalidate();
            }
        } ) );

        _propertyTransFunc[i].onChange( LambdaFunctionCallback( [this, i]
        {
            if( !_axesport.hasData() ) return;
            this->updateTransferFunction( i );
        } ) );
    }
    _propertySelectedSamples.onChange( LambdaFunctionCallback( [this]
    {
        if( !_axesport.hasData() ) return;

        const auto& values = _propertySelectedSamples.get();
        if( values.size() != _samplesSelection.size() ) return;

        _samplesSelection = values;
        for( auto& sections : _sections ) sections.clear();
        this->invalidate();
    } ) );
}

void ParallelCoordinatesViewer::setDescriptions() {
    setDescription("Viewer for parallel coordinates created by ParallelCoordinatesCreator");
    _propertySelectedMember.setDescription("Currently displayed member");
    _propertyVisualizationMode.setDescription("Visualization mode for currently displayed member:<br>"
                                              "<strong>Single Time Step</strong>: Each axis represents a field<br>"
                                              "<strong>All Time Steps</strong>: Each axis represents a time step<br>");
    _propertySelectedTimestep.setDescription("If single Time Step visualization mode is selected, "
                                             "this property determines the time step to be displayed");
    _propertySelectedTimestep.setDescription("If all Time Steps are visualized, this property determines<br>"
                                             "the field to be displayed");
    _propertyDensityBlending.setDescription("Enable to use a density visualization for rendered lines");
    //_propertyDensityVisibleSamples.setDescription("Determines the mini");
    //_propertyDensitySelectedSamples.setDescription("");

    _propertyTransFuncField.front().setDescription("Assigns the respective transfer function to the selected field");
    _propertyTransFunc.front().setDescription("Transfer function to be used for the selected field");
}

void ParallelCoordinatesViewer::initialize()
{
    RenderProcessor::initialize();

    // --- Helper function to compile a shader --- //
    const auto compileShader = [] ( GLenum type, const char* shaderCode )
    {
        GLuint shader = glCreateShader( type );
        glShaderSource( shader, 1, &shaderCode, nullptr );
        glCompileShader( shader );

        int success;
        glGetShaderiv( shader, GL_COMPILE_STATUS, &success );
        if( !success )
        {
            char infoLog[512];
            glGetShaderInfoLog( shader, 512, nullptr, infoLog );
            LERROR("Shader compilation failed\n" << infoLog);
            shader = 0;
        }

        return shader;
    };

    // --- Create shader program --- //
    const auto vertexShader = compileShader( GL_VERTEX_SHADER,
        "#version 330 core\n"
        "\n"
        "uniform int operator;\n"
        "uniform int operand;\n"
        "uniform vec3 axes[128];\n"
        "layout( location = 0 ) in float inYCoord;\n"
        "\n"
        "void main()\n"
        "{\n"
        "	int index = int( operator == 0? ( gl_VertexID % operand ) : ( gl_VertexID / operand ) );"
        "	vec3 axis = axes[index];\n"
        "	float y = ( inYCoord - axis.y ) / ( axis.z - axis.y );\n"
        "   gl_Position = vec4( axis.x, 0.9 * ( y * 2.0 - 1.0 ), 0.0, 1.0 );\n"
        "}"
    );
    const auto fragmentShader = compileShader( GL_FRAGMENT_SHADER,
        "#version 330 core\n"
        "\n"
        "uniform vec4 color;\n"
        "\n"
        "out vec4 outColor;\n"
        "\n"
        "void main()\n"
        "{\n"
        "   outColor = color; // vec4( vec3( 0.8 ), 0.95 / 110683 );\n"
        "}"
    );
    if( !vertexShader || !fragmentShader ) return;

    _shaderProgram = glCreateProgram();
    glAttachShader( _shaderProgram, vertexShader );
    glAttachShader( _shaderProgram, fragmentShader );
    glLinkProgram( _shaderProgram );

    int success;
    glGetProgramiv( _shaderProgram, GL_LINK_STATUS, &success );
    if( !success )
    {
        char infoLog[512];
        glGetProgramInfoLog( _shaderProgram, 512, nullptr, infoLog );
        LERROR("Shader program linking failed\n" << infoLog);
    }

    glDeleteShader( vertexShader );
    glDeleteShader( fragmentShader );

    // --- Generate VAO & VBO --- //
    glGenVertexArrays( 1, &_vertexArray );
    glBindVertexArray( _vertexArray );

    glGenBuffers( 1, &_indexBuffer );
    glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, _indexBuffer );
}
void ParallelCoordinatesViewer::deinitialize()
{
    RenderProcessor::deinitialize();

    glDeleteVertexArrays( 1, &_vertexArray );
    glDeleteBuffers( 1, &_indexBuffer );
    glDeleteProgram( _shaderProgram );
}
void ParallelCoordinatesViewer::process()
{
    if( !_axesport.hasData() ) return;
    const auto axes = _axesport.getData();
    const auto axisCount = _propertyVisualizationMode.getValue() ? axes->timesteps() : axes->fields();

    // --- Clear Background --- //
    _renderport.activateTarget();
    glClearColor( 0.95f, 0.95f, 0.95f, _propertyDensityVisibleSamples.get() == 1 ? 1.0f : 0.2f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    // --- Draw Lines --- //
    glUseProgram( _shaderProgram );
    glDisable( GL_DEPTH_TEST );
    glEnable( GL_BLEND );

    const auto offset = _propertyVisualizationMode.getValue() ? reinterpret_cast<const void*>( _propertySelectedMember.getValue() * axes->getStrideMember() ) : reinterpret_cast<const void*>( _propertySelectedMember.getValue() * axes->getStrideMember() + _propertySelectedTimestep.get() * axes->getStrideTimestep() );
    glBindVertexArray( _vertexArray );
    glEnableVertexAttribArray( 0 );
    glBindBuffer( GL_ARRAY_BUFFER, axes->getVertexBuffer() );
    glVertexAttribPointer( 0, 1, GL_FLOAT, GL_FALSE, 0, offset );

    // Draw unselected lines
    glBlendFuncSeparate( GL_DST_ALPHA, GL_ZERO, GL_ONE, GL_ONE );
    glUniform4f( glGetUniformLocation( _shaderProgram, "color" ), 0.259f, 0.529f, 0.961f, 0.8f / _propertyDensityVisibleSamples.get() );
    for( size_t i = 0; i < axes->samples(); ++i )  if( _samplesVisiblity[i] && !_samplesSelection[i] )
    {
        const auto baseVertex = _propertyVisualizationMode.getValue() ? static_cast<GLint>( _propertySelectedField.getValue() + i * axes->fields() ) : static_cast<GLint>( i * axes->fields() );
        glDrawElementsBaseVertex( GL_LINE_STRIP, static_cast<GLsizei>( axisCount ), GL_UNSIGNED_INT, nullptr, baseVertex );
    }

    // Reset background alpha
    if( !_propertyDensityBlending.get() )
    {
        glBindVertexArray( 0 );
        glUseProgram( 0 );

        // Reset alpha channel
        glBegin( GL_TRIANGLE_STRIP );
        glBlendFuncSeparate( GL_ZERO, GL_ONE, GL_ONE, GL_ZERO );
        glColor4f( 1.0f, 0.0f, 0.0f, _propertyDensitySelectedSamples.get() == 1 ? 1.0f : 0.2f );
        glVertex2f( -1.0f, 1.0f );
        glVertex2f( -1.0f, -1.0f );
        glVertex2f( 1.0f, 1.0f );
        glVertex2f( 1.0f, -1.0f );
        glEnd();

        glUseProgram( _shaderProgram );
        glBindVertexArray( _vertexArray );
    }

    // Draw selected lines
    glBlendFuncSeparate( GL_DST_ALPHA, GL_ZERO, GL_ONE, GL_ONE );
    glUniform4f( glGetUniformLocation( _shaderProgram, "color" ), 0.922f, 0.251f, 0.204f, 0.8f / _propertyDensitySelectedSamples.get() );
    for( size_t i = 0; i < axes->samples(); ++i )  if( _samplesVisiblity[i] && _samplesSelection[i] )
    {
        const auto baseVertex = _propertyVisualizationMode.getValue() ? static_cast<GLint>( _propertySelectedField.getValue() + i * axes->fields() ) : static_cast<GLint>( i * axes->fields() );
        glDrawElementsBaseVertex( GL_LINE_STRIP, static_cast<GLsizei>( axisCount ), GL_UNSIGNED_INT, nullptr, baseVertex );
    }

    glBindVertexArray( 0 );
    glUseProgram( 0 );
    glDisable( GL_BLEND );

    // --- Draw Axes & Top Boundary --- //
    glBegin( GL_LINES );
    glEnable( GL_LINE_SMOOTH );
    glLineWidth( 3.0f );

    glColor4f( 0.0f, 0.0f, 0.0f, 1.0f );
    for( size_t i = 0; i < axisCount; ++i )
    {
        const auto x = _uniformBufferVec[i].x;
        glVertex2f( x, -Y_LIMIT );
        glVertex2f(x, Y_LIMIT );
    }
    glVertex2f(-1.0f, Y_LIMIT );
    glVertex2f(1.0f, Y_LIMIT );

    // --- Draw sections --- //
    glColor4f( 0.988f, 0.729f, 0.012f, 1.0f );
    for( size_t i = 0; i < axisCount; ++i )
    {
        const auto& sections = _sections[i];
        const auto uniform = _uniformBufferVec[i];

        for( auto& section : sections )
        {
            const auto low = Y_LIMIT * ((section.first - uniform.y ) / (uniform.z - uniform.y ) * 2.0f - 1.0f );
            const auto high = Y_LIMIT * ((section.second - uniform.y ) / (uniform.z - uniform.y ) * 2.0f - 1.0f );
            if( high < -Y_LIMIT || low > Y_LIMIT ) continue;

            glVertex2f( uniform.x, std::max(-Y_LIMIT, low ) );
            glVertex2f( uniform.x, std::min(Y_LIMIT, high ) );
        }
    }

    // --- Setup font --- //
    const auto axisNamesY = util::pointToScreen(tgt::vec2(0.0f, Y_LIMIT ), _renderport.getSize() ).y + 10;
    auto font = tgt::Font( VoreenApplication::app()->getFontPath( "VeraMono.ttf" ), _renderport.getSize().y - axisNamesY - 10, 0.0f );
    font.setTextAlignment( tgt::Font::TextAlignment::TopRight );
    font.setFontColor( tgt::vec4( 0.05f, 0.05f, 0.05f, 1.0f ) );

    // --- Draw active section --- //
    if( _interaction == Interaction::eSelection )
    {
        const auto index = _propertyVisualizationMode.getValue() ? _hoveredScreenAxis : _indexBufferVec[_hoveredScreenAxis];
        const auto x = this->screenAxisToCoordinate( _hoveredScreenAxis );
        const auto uniform = _uniformBufferVec[index];

        glColor4f( 0.741f, 0.547f, 0.009f, 1.0f );
        auto y = ( _activeSection.first - uniform.y ) / ( uniform.z - uniform.y );
        glVertex2f(x, Y_LIMIT * (y * 2.0f - 1.0f ) );
        y = ( _activeSection.second - uniform.y ) / ( uniform.z - uniform.y );
        glVertex2f(x, Y_LIMIT * (y * 2.0f - 1.0f ) );

        const auto minmax = std::minmax( _activeSection.first, _activeSection.second );

        auto str = std::to_string( minmax.first );
        if( minmax.second != minmax.first ) str += " - " + std::to_string( minmax.second );

        const auto textWidth = font.getSize( tgt::vec3( 512.0f, 512.0f, 0.0f ), str, _renderport.getSize() ).x;
        auto offset = textWidth / 2;
        if( x - offset < 10 ) offset += 10 - ( x - offset );
        else if( x + offset > _renderport.getSize().x - 10 ) offset -= ( x + offset ) - ( _renderport.getSize().x - 10 );

        font.render( tgt::vec3( x + offset, 10.0f, 0.0f ), str, _renderport.getSize() );
    }

    glEnd();
    glDisable( GL_LINE_SMOOTH );
    glLineWidth( 1.0f );

    // --- Draw axes names --- //
    if( _propertyVisualizationMode.getValue() ) font.setFontSize( 16 );

    for( size_t i = 0; i < axisCount; ++i )
    {
        const auto text = _propertyVisualizationMode.getValue() ? std::to_string( i ) : axes->getFieldName( i );
        const auto x = util::pointToScreen( tgt::vec2( _uniformBufferVec[i].x, 0.0f ), _renderport.getSize() ).x;
        const auto textWidth = font.getSize( tgt::vec3( 512.0f, 512.0f, 0.0f ), text, _renderport.getSize() ).x;

        auto offset = textWidth / 2;
        if( x - offset < 10 ) offset += 10 - ( x - offset );
        else if( x + offset > _renderport.getSize().x - 10 ) offset -= ( x + offset ) - ( _renderport.getSize().x - 10 );

        font.setFontColor( _hoveredScreenAxis != ~0u && i == ( _propertyVisualizationMode.getValue() ? _hoveredScreenAxis : _indexBufferVec[_hoveredScreenAxis] ) ? tgt::vec4( 0.3f, 0.3f, 0.3f, 1.0f ) : tgt::vec4( 0.05f, 0.05f, 0.05f, 1.0f ) );
        font.render( tgt::vec3( x + offset, static_cast<float>( axisNamesY ), 0.0f ), text, _renderport.getSize() );
    }

    // --- Draw axes intervals --- //
    {
        font.setFontColor( tgt::vec4( 0.05f, 0.05f, 0.05f, 1.0f ) );
        font.setTextAlignment( tgt::Font::TopLeft );
        font.setFontSize( 16 );

        const std::array<float, 3> ycoords {
            static_cast<float>( util::pointToScreen(tgt::vec2( 0.0f, -Y_LIMIT ), _renderport.getSize() ).y ),
            static_cast<float>( util::pointToScreen( tgt::vec2( 0.0f, 0.0f ), _renderport.getSize() ).y ),
            static_cast<float>( util::pointToScreen(tgt::vec2(0.0f, Y_LIMIT ), _renderport.getSize() ).y )
        };

        const std::array<tgt::Font::TextAlignment, 4> alignments {
            tgt::Font::TopLeft, tgt::Font::TopRight, tgt::Font::BottomLeft, tgt::Font::BottomRight
        };

        for( size_t i = 0; i < ( _propertyVisualizationMode.getValue() ? 1 : axisCount ); ++i )
        {
            const auto uniform = _uniformBufferVec[i];
            const auto x = util::pointToScreen( tgt::vec2( uniform.x, 0.0f ), _renderport.getSize() ).x;

            const std::array<std::string, 3> texts {
                std::to_string( uniform.y ),
                std::to_string( ( uniform.y + uniform.z ) / 2.0f ),
                std::to_string( uniform.z )
            };

            for( size_t j = 0; j < 3; ++j )
            {
                const auto textWidth = font.getSize( tgt::vec3( 512.0f, 512.0f, 0.0f ), texts[j], _renderport.getSize() ).x;

                auto xoffset = 5.0f;
                auto alignment = ( j == 2 ) ? tgt::Font::BottomLeft : tgt::Font::TopLeft;
                if( x + textWidth > _renderport.getSize().x - 5 )
                {
                    alignment = static_cast<tgt::Font::TextAlignment>( alignment + 2 );
                    xoffset = -5.0f;
                }

                font.setTextAlignment( alignment );
                font.render( tgt::vec3( x + xoffset, ycoords[j], 0.0f ), texts[j], _renderport.getSize() );
            }
        }
    }

    glClearColor( 0.0f, 0.0f, 0.0f, 0.0f );
    glBlendFunc(GL_ONE, GL_ZERO);
    glDisable( GL_BLEND );

    _renderport.deactivateTarget();
}

void ParallelCoordinatesViewer::onNewInportData()
{
    const auto axes = _axesport.getData();
    const auto axisCount = _propertyVisualizationMode.getValue() ? axes->timesteps() : axes->fields();

    auto memberOptions = std::deque<Option<int>>();
    for( size_t i = 0; i < axes->members(); ++i )
    {
        const auto name = axes->getMemberName( i );
        memberOptions.push_back( Option<int>( name, name, static_cast<int>( i ) ) );
    }

    // Selected member
    _propertySelectedMember.blockCallbacks( true );
    _propertySelectedMember.setOptions( memberOptions );
    _propertySelectedMember.invalidate();
    _propertySelectedMember.blockCallbacks( false );

    // Selected timestep
    _propertySelectedTimestep.blockCallbacks( true );
    _propertySelectedTimestep.setMaxValue( axes->timesteps() - 1 );
    _propertySelectedTimestep.setVisibleFlag( !_propertyVisualizationMode.getValue() );
    _propertySelectedTimestep.blockCallbacks( false );

    // Selected field
    auto fieldOptions = std::deque<Option<int>>();
    for( size_t i = 0; i < axes->fields(); ++i )
    {
        const auto name = axes->getFieldName( i );
        fieldOptions.push_back( Option<int>( name, name, static_cast<int>( i ) ) );
    }
    _propertySelectedField.blockCallbacks( true );
    _propertySelectedField.setOptions( fieldOptions );
    _propertySelectedField.invalidate();
    _propertySelectedField.blockCallbacks( false );

    // Field interval
    const auto fieldInterval = axes->getRange( _propertySelectedField.getValue() );
    _propertyFieldInterval.blockCallbacks( true );
    _propertyFieldInterval.setMinValue( fieldInterval.first );
    _propertyFieldInterval.setMaxValue( fieldInterval.second );
    _propertyFieldInterval.set( tgt::vec2( fieldInterval.first, fieldInterval.second ) );
    _propertyFieldInterval.blockCallbacks( false );

    // Density visible samples
    _propertyDensityVisibleSamples.blockCallbacks( true );
    _propertyDensityVisibleSamples.setMaxValue( axes->samples() );
    _propertyDensityVisibleSamples.blockCallbacks( false );

    // Density selected samples
    _propertyDensitySelectedSamples.blockCallbacks( true );
    _propertyDensitySelectedSamples.setMaxValue( axes->samples() );
    _propertyDensitySelectedSamples.blockCallbacks( false );

    // Transfer function fields
    fieldOptions.push_front( Option<int>( "None", "None", -1 ) );
    for( size_t i = 0; i < _propertyTransFuncField.size(); ++i )
    {
        _propertyTransFuncField[i].blockCallbacks( true );
        _propertyTransFuncField[i].setOptions( fieldOptions );
        _propertyTransFuncField[i].setVisibleFlag( !_propertyVisualizationMode.getValue() );
        _propertyTransFuncField[i].invalidate();
        _propertyTransFuncField[i].blockCallbacks( false );
    }

    // Transfer functions
    for( size_t i = 0; i < _propertyTransFunc.size(); ++i )
    {
        _propertyTransFunc[i].blockCallbacks( true );
        _propertyTransFunc[i].setVisibleFlag( !_propertyVisualizationMode.getValue() );
        _propertyTransFunc[i].blockCallbacks( false );
    }

    // Selected samples
    _propertySelectedSamples.blockCallbacks( true );
    _propertySelectedSamples.set( std::vector<int>( axes->samples(), false ) );
    _propertySelectedSamples.blockCallbacks( false );

    // Index buffer
    _indexBufferVec.resize( axisCount );
    for( size_t i = 0; i < axisCount; ++i )
    {
        _indexBufferVec[i] = _propertyVisualizationMode.getValue() ?
            static_cast<GLuint>( i * axes->fields() * axes->samples() ) : static_cast<GLuint>( i );
    }
    glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, _indexBuffer );
    glBufferData( GL_ELEMENT_ARRAY_BUFFER, sizeof( GLuint ) * axisCount, _indexBufferVec.data(), GL_STATIC_DRAW );

    // Uniform buffer
    _uniformBufferVec.resize( axisCount );
    for( size_t i = 0; i < axisCount; ++i )
    {
        const auto range = axes->getRange( _propertyVisualizationMode.getValue() ? _propertySelectedField.getValue() : i );
        _uniformBufferVec[i] = tgt::vec3( this->screenAxisToCoordinate( i ), range.first, range.second );
    }
    glUseProgram( _shaderProgram );
    glUniform1i( glGetUniformLocation( _shaderProgram, "operator" ), _propertyVisualizationMode.getValue() ? 1 : 0 );
    glUniform1i( glGetUniformLocation( _shaderProgram, "operand" ), static_cast<GLint>( _propertyVisualizationMode.getValue() ? axes->samples() * axes->fields() : axes->fields() ) );
    glUniform3fv( glGetUniformLocation( _shaderProgram, "axes" ), static_cast<GLsizei>( _uniformBufferVec.size() ), reinterpret_cast<const GLfloat*>( _uniformBufferVec.data() ) );
    glUseProgram( 0 );

    // Interaction data
    _interaction = Interaction::eNone;
    _hoveredScreenAxis = ~0u;
    _activeSection = std::pair<float, float>( 0.0f, 0.0f );
    _sections = std::vector<std::list<std::pair<float, float>>>( axisCount );
    _samplesVisiblity = std::vector<bool>( axes->samples(), true );
    _samplesSelection = std::vector<int>( axes->samples(), false );
}
void ParallelCoordinatesViewer::onHoverEvent( tgt::MouseEvent* event )
{
    if( !_axesport.hasData() ) return;
    const auto axes = _axesport.getData();
    const auto axisCount = _propertyVisualizationMode.getValue() ? axes->timesteps() : axes->fields();

    const auto point = util::screenToPoint( event->coord(), event->viewport() );

    struct { size_t axis; int distance; } current { ~0u, point.y > Y_LIMIT ? 50 : 10 };
    for( size_t i = 0; i < axisCount; ++i )
    {
        const auto screen = util::pointToScreen( tgt::vec2( this->screenAxisToCoordinate( i ), 0.0f ), event->viewport() );
        const auto distance = std::abs( screen.x - event->x() );
        if( distance < current.distance ) current = { i, distance };
    }

    if( current.axis != _hoveredScreenAxis )
    {
        _hoveredScreenAxis = current.axis;
        this->process();
    }
}
void ParallelCoordinatesViewer::onMouseEvent( tgt::MouseEvent* event )
{
    if( !_axesport.hasData() || _hoveredScreenAxis == ~0u ) return;

    const auto axes = _axesport.getData();
    const auto axisCount = _propertyVisualizationMode.getValue() ? axes->timesteps() : axes->fields();
    const auto axisIndex = _propertyVisualizationMode.getValue() ? _hoveredScreenAxis : _indexBufferVec[_hoveredScreenAxis];

    const auto point = util::screenToPoint( event->coord(), event->viewport() );
    auto updateUniformBuffer = false, render = false;

    if( event->action() == tgt::MouseEvent::MouseAction::PRESSED )
    {
        if(point.y > Y_LIMIT && !_propertyVisualizationMode.getValue() ) // Axes in visualization mode 'timesteps' cant be moved
        {
            _interaction = Interaction::eMoving;
        }
        else
        {
            const auto uniform = _uniformBufferVec[axisIndex];
            const auto x = std::max( 0.0f, std::min( 1.0f, ((point.y / Y_LIMIT ) + 1.0f ) / 2.0f ) );
            _activeSection.first = _activeSection.second = uniform.y + x * ( uniform.z - uniform.y );
            _interaction = Interaction::eSelection;
            render = true;
        }
    }
    else if( event->action() == tgt::MouseEvent::MouseAction::RELEASED )
    {
        if( _interaction == Interaction::eMoving ) // Only occurs in visualization mode 'fields'
        {
            _uniformBufferVec[axisIndex].x = this->screenAxisToCoordinate( _hoveredScreenAxis );
            updateUniformBuffer = render = true;
        }
        else if( _interaction == Interaction::eSelection )
        {
            // Swap value so that first < second
            if( _activeSection.first > _activeSection.second )
                std::swap( _activeSection.first, _activeSection.second );

            // Add or remove section
            auto& sections = _sections[axisIndex];
            if( event->button() == tgt::MouseEvent::MOUSE_BUTTON_LEFT ) // Add section
            {
                sections.insert( std::upper_bound( sections.begin(), sections.end(), _activeSection ), _activeSection );
                for( auto it = sections.begin(); std::next( it ) != sections.end(); )
                {
                    auto next = std::next( it );
                    if( it->second >= next->first )
                    {
                        it->second = std::max( it->second, next->second );
                        sections.erase( next );
                    }
                    else ++it;
                }
            }
            else // Remove section
            {
                for( auto it = sections.begin(); it != sections.end(); )
                {
                    if( it->first < _activeSection.first )
                    {
                        if( it->second <= _activeSection.first ) ++it;
                        else
                        {
                            if( it->second <= _activeSection.second )
                            {
                                it->second = _activeSection.first;
                                ++it;
                            }
                            else
                            {
                                sections.insert( std::next( it ), std::pair<float, float>( _activeSection.second, it->second ) );
                                it->second = _activeSection.first;
                                break;
                            }
                        }
                    }
                    else if( it->first < _activeSection.second )
                    {
                        if( it->second <= _activeSection.second ) it = sections.erase( it );
                        else
                        {
                            it->first = _activeSection.second;
                            break;
                        }
                    }
                    else break;
                }
            }

            this->updateSampleStates();
            if( !_propertyVisualizationMode.getValue() )
                for( size_t i = 0; i < _propertyTransFuncField.size(); ++i )
                    if( static_cast<size_t>(_propertyTransFuncField[i].getValue()) == axisIndex )
                        this->updateTransferFunction( i );


            _activeSection = std::pair<float, float>( 0.0f, 0.0f );
            render = true;
        }

        _interaction = Interaction::eNone;
    }
    else if( event->action() == tgt::MouseEvent::MouseAction::MOTION )
    {
        if( _interaction == Interaction::eMoving ) // Only occurs in visualization mode 'fields'
        {
            const auto x = std::max(-X_LIMIT - 0.01f, std::min(X_LIMIT + 0.01f, point.x ) );
            const auto prevx = this->screenAxisToCoordinate( _hoveredScreenAxis - 1 );
            const auto nextx = this->screenAxisToCoordinate( _hoveredScreenAxis + 1 );

            // Check if axes have to be swapped
            size_t swap = ~0u;
            if( _hoveredScreenAxis > 0 && x < prevx ) swap = _hoveredScreenAxis - 1;
            else if( _hoveredScreenAxis < axes->fields() - 1 && x > nextx ) swap = _hoveredScreenAxis + 1;

            if( swap != ~0u )
            {
                _uniformBufferVec[_indexBufferVec[_hoveredScreenAxis]].x = this->screenAxisToCoordinate( _hoveredScreenAxis );
                std::swap( _indexBufferVec[_hoveredScreenAxis], _indexBufferVec[swap] );
                std::swap( _uniformBufferVec[_indexBufferVec[_hoveredScreenAxis]].x, _uniformBufferVec[_indexBufferVec[swap]].x );

                glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, _indexBuffer );
                glBufferData( GL_ELEMENT_ARRAY_BUFFER, sizeof( GLuint ) * axisCount, _indexBufferVec.data(), GL_STATIC_DRAW );
                glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, 0 );

                _hoveredScreenAxis = swap;
            }
            else _uniformBufferVec[_indexBufferVec[_hoveredScreenAxis]].x = x;
            updateUniformBuffer = render = true;
        }
        else if( _interaction == Interaction::eSelection )
        {
            const auto uniform = _uniformBufferVec[axisIndex];
            const auto x = std::max( 0.0f, std::min( 1.0f, ((point.y / Y_LIMIT ) + 1.0f ) / 2.0f ) );
            _activeSection.second = uniform.y + x * ( uniform.z - uniform.y );
            render = true;
        }
    }

    // --- Perform Updates If Requested --- //
    if( updateUniformBuffer )
    {
        glUseProgram( _shaderProgram );
        glUniform3fv( glGetUniformLocation( _shaderProgram, "axes" ), static_cast<GLsizei>( _uniformBufferVec.size() ), reinterpret_cast<const GLfloat*>( _uniformBufferVec.data() ) );
        glUseProgram( 0 );
    }
    if( render ) this->process();
}

void ParallelCoordinatesViewer::updateSampleStates()
{
    const auto axes = _axesport.getData();
    const auto axisCount = _propertyVisualizationMode.getValue() ? axes->timesteps() : axes->fields();
    const auto sectionsCount = std::accumulate( _sections.begin(), _sections.end(), size_t( 0 ),
        [] ( size_t val, const std::list<std::pair<float, float>>& sections ) { return val + sections.size(); } );

    std::fill( _samplesSelection.begin(), _samplesSelection.end(), false );
    for( size_t i = 0; i < axes->samples(); ++i )
    {
        _samplesVisiblity[i] = true;

        for( size_t j = 0; j < axisCount; ++j )
        {
            const auto uniform = _uniformBufferVec[j];
            const auto field = _propertyVisualizationMode.getValue() ? _propertySelectedField.getValue() : j;
            const auto timestep = _propertyVisualizationMode.getValue() ? j : _propertySelectedTimestep.get();
            const auto value = axes->getValue( field, i, timestep, _propertySelectedMember.getValue() );

            if( value < uniform.y || value > uniform.z )
            {
                _samplesVisiblity[i] = false;
                break;
            }
        }


        if( sectionsCount && _samplesVisiblity[i] )
        {
            _samplesSelection[i] = true;
            for( size_t j = 0; j < axisCount; ++j ) if( _sections[j].size() )
            {
                const auto field = _propertyVisualizationMode.getValue() ? _propertySelectedField.getValue() : j;
                const auto timestep = _propertyVisualizationMode.getValue() ? j : _propertySelectedTimestep.get();
                const auto value = axes->getValue( field, i, timestep, _propertySelectedMember.getValue() );

                bool select = false;
                for( const auto section : _sections[j] ) if( value >= section.first && value <= section.second )
                {
                    select = true;
                    break;
                }

                if( !select )
                {
                    _samplesSelection[i] = false;
                    break;
                }
            }
        }
    }

    _propertySelectedSamples.blockCallbacks( true );
    _propertySelectedSamples.set( _samplesSelection );
    _propertySelectedSamples.blockCallbacks( false );

    _propertySections.set( ParallelCoordinatesSectionsPropertyData( _propertySelectedMember.getDescription(), static_cast<size_t>( _propertySelectedTimestep.get() ), axes->getFieldNames(), _sections ) );
}
void ParallelCoordinatesViewer::updateTransferFunction( size_t index )
{
    if( _propertyTransFuncField[index].getValue() < 0 ) return;

    const auto inputKeys = _propertyTransFunc[index].get()->getKeys();
    const auto interpolationFromIntensity = [&inputKeys] ( float intensity )
    {
        for( size_t i = 0; i < inputKeys.size(); ++i )
        {
            if( inputKeys[i]->getIntensity() >= intensity )
            {
                auto lowerIntensity = 0.0f, upperIntensity = inputKeys[i]->getIntensity();
                auto lowerColor = tgt::col4(), upperColor = inputKeys[i]->getColorL();

                if( i == 0 ) lowerIntensity = 0.0f, lowerColor = tgt::col4( 0, 0, 0, 0 );
                else lowerIntensity = inputKeys[i - 1]->getIntensity(), lowerColor = inputKeys[i - 1]->getColorR();


                const auto x = ( intensity - lowerIntensity ) / ( upperIntensity - lowerIntensity );

                const auto r = static_cast<uint8_t>( x * upperColor.r + ( 1.0f - x ) * lowerColor.r );
                const auto g = static_cast<uint8_t>( x * upperColor.g + ( 1.0f - x ) * lowerColor.g );
                const auto b = static_cast<uint8_t>( x * upperColor.b + ( 1.0f - x ) * lowerColor.b );
                const auto a = static_cast<uint8_t>( x * upperColor.a + ( 1.0f - x ) * lowerColor.a );

                return tgt::col4( r, g, b, a );
            }
        }

        return tgt::col4( 0, 0, 0, 0 );
    };

    const auto axisIndex = _propertyTransFuncField[index].getValue();
    const auto& sections = _sections[axisIndex];
    auto transferFunctionKeys = std::unique_ptr<TransFunc1DKeys>( _propertyTransFunc[index].get()->clone() );
    auto keys = std::vector<TransFuncMappingKey*>();

    const auto range =  _uniformBufferVec[axisIndex].yz(); // tgt::vec2( sections.front().first, sections.back().second );
    transferFunctionKeys->setDomain( range.x, range.y );

    if( sections.size() )
    {
        for( const auto section : sections )
        {
            auto intensity = ( section.first - range.x ) / ( range.y - range.x );
            auto interpolation = interpolationFromIntensity( intensity );
            auto begin = new TransFuncMappingKey( intensity, interpolation );
            begin->setSplit( true );
            begin->setColorL( tgt::col4( 0, 0, 0, 0 ) );
            keys.push_back( begin );

            intensity = ( section.second - range.x ) / ( range.y - range.x );
            interpolation = interpolationFromIntensity( intensity );
            auto end = new TransFuncMappingKey( intensity, interpolation );
            end->setSplit( true );
            end->setColorR( tgt::col4( 0, 0, 0, 0 ) );
            keys.push_back( end );
        }
    }
    transferFunctionKeys->setKeys( std::move( keys ) );

    for( const auto key : inputKeys )
    {
        for( const auto section : sections )
        {
            const auto lower = ( section.first - range.x ) / ( range.y - range.x );
            const auto upper = ( section.second - range.x ) / ( range.y - range.x );
            if( key->getIntensity() >= lower && key->getIntensity() <= upper )
            {
                transferFunctionKeys->addKey( key->clone() );
                break;
            }
        }
    }

    _propertyTransFuncIntern[index].set1DKeys( transferFunctionKeys.release() );
    _propertyTransFuncIntern[index].invalidate();
}

float ParallelCoordinatesViewer::screenAxisToCoordinate( size_t axis ) const
{
    const auto axes = _axesport.getData();
    const auto axisCount = _propertyVisualizationMode.getValue() ? axes->timesteps() : axes->fields();
    return -X_LIMIT + axis * ((2.0f * X_LIMIT ) / (axisCount - 1 ) );
}

}
