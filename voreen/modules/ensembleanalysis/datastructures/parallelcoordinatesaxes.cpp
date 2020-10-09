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

#include "parallelcoordinatesaxes.h"
#include <fstream>

namespace voreen
{

ParallelCoordinatesAxes::ParallelCoordinatesAxes( std::vector<std::string> members, std::vector<std::string> fields, std::vector<std::pair<float, float>> ranges, std::vector<float> values, size_t numTimesteps, size_t numSamples )
    : timesteps_(numTimesteps ), samples_(numSamples ), members_(std::move(members ) ), fields_(std::move(fields ) ), ranges_(std::move(ranges ) ), values_(std::move(values ) ), vertexBuffer_(0 )
{}

ParallelCoordinatesAxes::ParallelCoordinatesAxes( const std::string& filepath ) : vertexBuffer_(0 )
{
    auto stream = std::ifstream( filepath, std::ios::in | std::ios::binary );
    if( !stream ) {
        std::cout << "[ERROR]: ParallelCoordinatesAxes::ParallelCoordinatesAxes --> Failed to open file " << filepath << std::endl;
        return;
    }

    size_t members, fields;
    stream.read( reinterpret_cast<char*>( &members ), sizeof( size_t ) );
    stream.read(reinterpret_cast<char*>( &timesteps_ ), sizeof( size_t ) );
    stream.read( reinterpret_cast<char*>( &fields ), sizeof( size_t ) );
    stream.read(reinterpret_cast<char*>( &samples_ ), sizeof( size_t ) );

    members_.resize(members );
    for( size_t i = 0; i < members; ++i ) {
        size_t length;
        stream.read( reinterpret_cast<char*>( &length ), sizeof( size_t ) );

        members_[i].resize(length );
        stream.read(reinterpret_cast<char*>( &members_[i][0] ), length );
    }

    fields_.resize(fields );
    for( size_t i = 0; i < fields; ++i ) {
        size_t length;
        stream.read( reinterpret_cast<char*>( &length ), sizeof( size_t ) );

        fields_[i].resize(length );
        stream.read(reinterpret_cast<char*>( &fields_[i][0] ), length );
    }

    ranges_.resize(fields );
    stream.read(reinterpret_cast<char*>( ranges_.data() ), ranges_.size() * sizeof( std::pair<float, float> ) );

    values_.resize(members * timesteps_ * fields * samples_ );
    stream.read(reinterpret_cast<char*>( values_.data() ), values_.size() * sizeof( float ) );

    if( !stream ) std::cout << "[ERROR]: ParallelCoordinatesAxes::ParallelCoordinatesAxes --> Failed read file " << filepath << std::endl;
}
ParallelCoordinatesAxes::~ParallelCoordinatesAxes() {
    if( vertexBuffer_ ) glDeleteBuffers(1, &vertexBuffer_ );
}

void ParallelCoordinatesAxes::serialize( const std::string& filepath ) const {

    auto stream = std::ofstream( filepath, std::ios::out | std::ios::binary );
    if( !stream ) {
        std::cout << "[ERROR]: ParallelCoordinatesAxes::serialize --> Failed to open file " << filepath << std::endl;
        return;
    }

    const auto members = this->members(), timesteps = this->timesteps(), fields = this->fields(), samples = this->samples();
    stream.write( reinterpret_cast<const char*>( &members ), sizeof( size_t ) );
    stream.write( reinterpret_cast<const char*>( &timesteps ), sizeof( size_t ) );
    stream.write( reinterpret_cast<const char*>( &fields ), sizeof( size_t ) );
    stream.write( reinterpret_cast<const char*>( &samples ), sizeof( size_t ) );

    for( size_t i = 0; i < members; ++i ) {
        const auto length = members_[i].size();
        stream.write( reinterpret_cast<const char*>( &length ), sizeof( size_t ) );
        stream.write(reinterpret_cast<const char*>( members_[i].data() ), length );
    }

    for( size_t i = 0; i < fields; ++i ) {
        const auto length = fields_[i].size();
        stream.write( reinterpret_cast<const char*>( &length ), sizeof( size_t ) );
        stream.write(reinterpret_cast<const char*>( fields_[i].data() ), length );
    }

    stream.write(reinterpret_cast<const char*>( ranges_.data() ), ranges_.size() * sizeof( std::pair<float, float> ) );
    stream.write(reinterpret_cast<const char*>( values_.data() ), values_.size() * sizeof( float ) );

    if( !stream ) std::cout << "[ERROR]: ParallelCoordinatesAxes::serialize --> Failed write file " << filepath << std::endl;
}

size_t ParallelCoordinatesAxes::members() const noexcept {
    return members_.size();
}
size_t ParallelCoordinatesAxes::timesteps() const noexcept {
    return timesteps_;
}
size_t ParallelCoordinatesAxes::fields() const noexcept {
    return fields_.size();
}
size_t ParallelCoordinatesAxes::samples() const noexcept {
    return samples_;
}

const std::string& ParallelCoordinatesAxes::getMemberName( size_t i ) const {
    return members_[i];
}
const std::string& ParallelCoordinatesAxes::getFieldName( size_t i ) const {
    return fields_[i];
}

const std::vector<std::string>& ParallelCoordinatesAxes::getMemberNames() const noexcept {
    return members_;
}
const std::vector<std::string>& ParallelCoordinatesAxes::getFieldNames() const noexcept {
    return fields_;
}

std::pair<float, float> ParallelCoordinatesAxes::getRange( size_t field ) const {
    return ranges_[field];
}
const std::vector<std::pair<float, float>>& ParallelCoordinatesAxes::getRanges() const noexcept {
    return ranges_;
}

float ParallelCoordinatesAxes::getValue( size_t field, size_t sample, size_t timestep, size_t member ) const {
    return values_[(this->timesteps() * this->samples() * this->fields() * member ) + (this->samples() * this->fields() * timestep ) + (this->fields() * sample ) + field];
}
const std::vector<float>& ParallelCoordinatesAxes::getValues() const noexcept {
    return values_;
}

size_t ParallelCoordinatesAxes::getStrideMember() const noexcept {
    return timesteps_ * fields_.size() * samples_ * sizeof( float );
}

size_t ParallelCoordinatesAxes::getStrideTimestep() const noexcept {
    return fields_.size() * samples_ * sizeof( float );
}

size_t ParallelCoordinatesAxes::memorySize() const noexcept {
    return values_.size() * sizeof( float ) + ranges_.size() * sizeof( std::pair<float, float> );
}

GLuint ParallelCoordinatesAxes::getVertexBuffer() const {
    if( !vertexBuffer_ ) {
        glGenBuffers( 1, &vertexBuffer_ );
        glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer_ );
        glBufferData(GL_ARRAY_BUFFER, values_.size() * sizeof( float ), values_.data(), GL_STATIC_DRAW );
    }

    return vertexBuffer_;
}

}