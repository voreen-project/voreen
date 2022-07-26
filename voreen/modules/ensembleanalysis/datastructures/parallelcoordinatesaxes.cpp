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

#include "parallelcoordinatesaxes.h"
#include <fstream>

#include "voreen/core/io/serialization/deserializer.h"
#include "voreen/core/io/serialization/serializer.h"

#include "voreen/core/io/serialization/jsondeserializer.h"
#include "voreen/core/io/serialization/jsonserializer.h"

#include "voreen/core/utils/voreenfilepathhelper.h"

#include "tgt/filesystem.h"

namespace voreen
{

ParallelCoordinatesAxes::ParallelCoordinatesAxes( std::string ensembleHash,
                                                  std::vector<std::string> members,
                                                  std::vector<std::pair<std::string, int>> fields,
                                                  std::vector<std::string> axesLabels,
                                                  std::vector<tgt::vec2> ranges,
                                                  std::vector<float> values,
                                                  tgt::Bounds bounds,
                                                  size_t numTimesteps,
                                                  size_t numSamples )
    : timesteps_(numTimesteps)
    , samples_(numSamples)
    , ensembleHash_(std::move(ensembleHash))
    , members_(std::move(members))
    , fields_(std::move(fields))
    , axesLabels_(std::move(axesLabels))
    , ranges_(std::move(ranges))
    , values_(std::move(values))
    , bounds_(bounds)
    , vertexBuffer_(0)
{}

ParallelCoordinatesAxes::ParallelCoordinatesAxes( const std::string& filepath )
    : vertexBuffer_(0)
{
    auto stream = std::ifstream( filepath );
    if( !stream ) {
        LERRORC("voreen.ParallelCoordinateAxes", "Failed to open file " << filepath);
        return;
    }

    try {
        JsonDeserializer json;
        Deserializer s(json);
        json.read(stream, false);
        s.deserialize("ensembleHash", ensembleHash_);
        s.deserialize("members", members_);
        s.deserialize("fields", fields_);
        s.deserialize("axesLabels", axesLabels_);
        s.deserialize("timesteps", timesteps_);
        s.optionalDeserialize("bounds", bounds_, tgt::Bounds()); // Have been added later.
        s.deserialize("samples", samples_);
        s.deserialize("ranges", ranges_);
        try {
            s.deserialize("values", values_);
        } catch (SerializationException& e) {
            s.removeLastError();
            try {
                s.deserializeBinaryBlob("values_binary", values_);
            } catch (SerializationException& e) {
                s.removeLastError();

                VoreenFilePathHelper tmp;
                s.deserialize("rawDataPath", tmp);
                size_t numValues = members_.size() * fields_.size() * timesteps_ * samples_;
                values_.resize(numValues);
                size_t numbytes = values_.size() * sizeof(float);

                std::fstream rawin(tmp.getPath().c_str(), std::ios::in | std::ios::binary);
                rawin.read(reinterpret_cast<char*>(values_.data()), numbytes);
            }
        }
        LINFOC("voreen.ParallelCoordinateAxes", "Loaded successfully");
    } catch(tgt::Exception& e) {
        LERRORC("voreen.ParallelCoordinateAxes", e.what());
    }
}
ParallelCoordinatesAxes::~ParallelCoordinatesAxes() {
    if( vertexBuffer_ ) glDeleteBuffers(1, &vertexBuffer_ );
}

void ParallelCoordinatesAxes::serialize( const std::string& filepath, bool binary) const {

    auto stream = std::ofstream( filepath );
    if( !stream ) {
        LERRORC("voreen.ParallelCoordinateAxes", "Failed to open file " << filepath);
        return;
    }

    try {
        JsonSerializer json;
        Serializer s(json);
        s.serialize("ensembleHash", ensembleHash_);
        s.serialize("members", members_);
        s.serialize("fields", fields_);
        s.serialize("axesLabels", axesLabels_);
        s.serialize("timesteps", timesteps_);
        s.serialize("bounds", bounds_);
        s.serialize("samples", samples_);
        s.serialize("ranges", ranges_);
        if (!binary) {
            s.serialize("values", values_);
        }
        else {
            s.serialize("rawDataPath", VoreenFilePathHelper(filepath));

            std::string rawname = tgt::FileSystem::fullBaseName(filepath) + ".raw";
            const char* data = reinterpret_cast<const char*>(values_.data());
            size_t numbytes = values_.size() * sizeof(float);

            std::fstream rawout(rawname.c_str(), std::ios::out | std::ios::binary);
            rawout.write(data, numbytes);
        }
        json.write(stream, true, false);
        LINFOC("voreen.ParallelCoordinateAxes", "Saved successfully");
    } catch(tgt::Exception& e) {
        LERRORC("voreen.ParallelCoordinateAxes", e.what());
    }
}

size_t ParallelCoordinatesAxes::members() const noexcept {
    return members_.size();
}
size_t ParallelCoordinatesAxes::fields() const noexcept {
    return fields_.size();
}
size_t ParallelCoordinatesAxes::timesteps() const noexcept {
    return timesteps_;
}
size_t ParallelCoordinatesAxes::samples() const noexcept {
    return samples_;
}

const std::string& ParallelCoordinatesAxes::getEnsembleHash() const {
    return ensembleHash_;
}

const std::string& ParallelCoordinatesAxes::getMemberName( size_t i ) const {
    return members_[i];
}
const std::string& ParallelCoordinatesAxes::getFieldName( size_t i ) const {
    return fields_[i].first;
}
int ParallelCoordinatesAxes::getChannel( size_t i) const {
    return fields_[i].second;
}

const std::vector<std::string>& ParallelCoordinatesAxes::getMemberNames() const noexcept {
    return members_;
}
const std::vector<std::pair<std::string, int>>& ParallelCoordinatesAxes::getFields() const noexcept {
    return fields_;
}
const std::vector<std::string>& ParallelCoordinatesAxes::getAxesLabels() const noexcept {
    return axesLabels_;
}

tgt::vec2 ParallelCoordinatesAxes::getRange( size_t field ) const {
    return ranges_[field];
}
const std::vector<tgt::vec2>& ParallelCoordinatesAxes::getRanges() const noexcept {
    return ranges_;
}

float ParallelCoordinatesAxes::getValue( size_t field, size_t sample, size_t timestep, size_t member ) const {
    return values_[(this->timesteps() * this->samples() * this->fields() * member ) + (this->samples() * this->fields() * timestep ) + (this->fields() * sample ) + field];
}
const std::vector<float>& ParallelCoordinatesAxes::getValues() const noexcept {
    return values_;
}

const tgt::Bounds& ParallelCoordinatesAxes::getBounds() const {
    return bounds_;
}

size_t ParallelCoordinatesAxes::getStrideMember() const noexcept {
    return timesteps_ * fields_.size() * samples_ * sizeof( float );
}

size_t ParallelCoordinatesAxes::getStrideTimestep() const noexcept {
    return fields_.size() * samples_ * sizeof( float );
}

size_t ParallelCoordinatesAxes::memorySize() const noexcept {
    return values_.size() * sizeof( float ) + ranges_.size() * sizeof( tgt::vec2 );
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
