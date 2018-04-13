/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "voreen/core/io/serialization/jsondeserializer.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/voreenmodule.h"

#include "voreen/core/io/serialization/deserializer.h"

#include "rapidjson/istreamwrapper.h"
#include "rapidjson/error/en.h"

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#if WIN32 // HACK: Needed due to bug #13164 (described boost bugtracker)
#define BOOST_IOSTREAMS_DYN_LINK
#endif
#include <boost/iostreams/filter/gzip.hpp>

namespace voreen {

const std::string JsonDeserializer::loggerCat_ = "voreen.JsonDeserializer";

JsonDeserializer::JsonDeserializer(std::string documentPath)
    : documentPath_(documentPath)
    , document_()
    , documentBuffer_()
    , node_(&document_) //TODO: we might want to do that later in read... (or just read in the constructor?)
{
    /*
     * TODO: do we need this?
    // register application (as proxy for modules)
    if (VoreenApplication::app()) {
        registerFactory(VoreenApplication::app());
        registerFactories(VoreenApplication::app()->getSerializerFactories());
    }
    else {
        LWARNING("VoreenApplication not instantiated");
    }
    */
}

JsonDeserializer::~JsonDeserializer() {
}

std::string JsonDeserializer::getDocumentPath() const {
    return documentPath_;
}

void JsonDeserializer::deserialize(const std::string& key, bool& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, char& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, signed char& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, unsigned char& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, uint16_t& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, int16_t& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, uint32_t& data) {
    deserializeFromMember(key, data);
}

#ifdef __APPLE__
void JsonDeserializer::deserialize(const std::string& key, long unsigned int& data) {
    deserializeFromMember(key, data);
}
#endif

void JsonDeserializer::deserialize(const std::string& key, int32_t& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, uint64_t& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, int64_t& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, float& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, double& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, long double& data) {
    throw SerializationException("long double is not supported");
}

void JsonDeserializer::deserialize(const std::string& key, std::string& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, tgt::vec2& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, tgt::vec3& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, tgt::vec4& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, tgt::dvec2& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, tgt::dvec3& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, tgt::dvec4& data) {
    deserializeFromMember(key, data);
}


void JsonDeserializer::deserialize(const std::string& key, tgt::ivec2& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, tgt::ivec3& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, tgt::ivec4& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, tgt::col3& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, tgt::col4& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, tgt::mat2& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, tgt::mat3& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, tgt::mat4& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, tgt::Matrix2d& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, tgt::Matrix3d& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, tgt::Matrix4d& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::deserialize(const std::string& key, Serializable& data) {
    deserializeFromMember(key, data);
}

void JsonDeserializer::read(std::istream& stream)
{
    using namespace boost::iostreams;

    // Prepare gzip compressing stream
    filtering_istream decompressingStream;
    decompressingStream.push(gzip_decompressor());
    decompressingStream.push(stream);

    documentBuffer_.clear();
    std::copy(std::istream_iterator<char>(decompressingStream), std::istream_iterator<char>(), std::back_inserter(documentBuffer_));
    documentBuffer_.push_back(0); //append string-terminating zero byte.
    rapidjson::ParseResult parse_result = document_.ParseInsitu<rapidjson::kParseNanAndInfFlag>(documentBuffer_.data());

    //rapidjson::IStreamWrapper s(decompressingStream);
    //rapidjson::ParseResult parse_result = document_.ParseStream<rapidjson::kParseNanAndInfFlag>(s);
    if (parse_result.IsError()) {
        throw SerializationException(std::string("Could not parse json: ") + rapidjson::GetParseError_En(parse_result.Code()) + " at " + std::to_string(parse_result.Offset()));
    }
    // Not sure if this is necessary, but it shouldn't hurt:
    node_ = &document_;
}

void JsonDeserializer::readFromValue(const rapidjson::Value& val, Serializable& data) {
    if(!val.IsObject()) {
        throw SerializationException("Not an object");
    } else {
        TemporaryNodeChanger nodeChanger(*this, val);
        Deserializer deserializer(*this);
        data.deserialize(deserializer);
    }
}

} // namespace
