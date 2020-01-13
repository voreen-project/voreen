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

#include "voreen/core/io/serialization/jsonserializer.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/voreenmodule.h"
#include "voreen/core/animation/animation.h"

#include "rapidjson/prettywriter.h"
#include "rapidjson/writer.h"
#include "rapidjson/ostreamwrapper.h"

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#if WIN32 // HACK: Needed due to bug #13164 (described boost bugtracker)
#define BOOST_IOSTREAMS_DYN_LINK
#endif
#include <boost/iostreams/filter/gzip.hpp>

#include "voreen/core/io/serialization/serializer.h"

namespace voreen {

const std::string JsonSerializer::loggerCat_ = "voreen.JsonSerializer";

JsonSerializer::JsonSerializer(std::string documentPath)
    : documentPath_(documentPath)
    , document_()
    , allocator_(document_.GetAllocator())
    , node_(&document_)
{
    node_->SetObject();
}

JsonSerializer::~JsonSerializer() {
}

std::string JsonSerializer::getDocumentPath() const {
    return documentPath_;
}

void JsonSerializer::resolveUnresolvedReferences() {
    // TODO ?
}

void JsonSerializer::serialize(const std::string& key, const bool& data) {
    serializeSimpleTypes(key, data);
}

void JsonSerializer::serialize(const std::string& key, const char& data) {
    serializeSimpleTypes(key, data);
}

void JsonSerializer::serialize(const std::string& key, const signed char& data) {
    serializeSimpleTypes(key, data);
}

void JsonSerializer::serialize(const std::string& key, const unsigned char& data) {
    serializeSimpleTypes(key, data);
}

void JsonSerializer::serialize(const std::string& key, const uint16_t& data) {
    serializeSimpleTypes(key, data);
}

void JsonSerializer::serialize(const std::string& key, const int16_t& data) {
    serializeSimpleTypes(key, data);
}

void JsonSerializer::serialize(const std::string& key, const uint32_t& data) {
    serializeSimpleTypes(key, data);
}

void JsonSerializer::serialize(const std::string& key, const int32_t& data) {
    serializeSimpleTypes(key, data);
}

#ifdef __APPLE__
void JsonSerializer::serialize(const std::string& key, const long unsigned int& data) {
    serializeSimpleTypes(key, data);
}
#endif

void JsonSerializer::serialize(const std::string& key, const uint64_t& data) {
    serializeSimpleTypes(key, data);
}

void JsonSerializer::serialize(const std::string& key, const int64_t& data) {
    serializeSimpleTypes(key, data);
}

void JsonSerializer::serialize(const std::string& key, const float& data) {
    serializeSimpleTypes(key, data);
}

void JsonSerializer::serialize(const std::string& key, const double& data) {
    serializeSimpleTypes(key, data);
}

void JsonSerializer::serialize(const std::string& key, const long double& data) {
    //serializeSimpleTypes(key, data);
    throw SerializationException("long double serialization is not supported");
}

void JsonSerializer::serialize(const std::string& key, const std::string& data) {
    serializeSimpleTypes(key, data);
}

void JsonSerializer::serialize(const std::string& key, const tgt::vec2& data) {
    serializeTgtVector(key, data);
}

void JsonSerializer::serialize(const std::string& key, const tgt::vec3& data) {
    serializeTgtVector(key, data);
}

void JsonSerializer::serialize(const std::string& key, const tgt::vec4& data) {
    serializeTgtVector(key, data);
}

void JsonSerializer::serialize(const std::string& key, const tgt::dvec2& data) {
    serializeTgtVector(key, data);
}

void JsonSerializer::serialize(const std::string& key, const tgt::dvec3& data) {
    serializeTgtVector(key, data);
}

void JsonSerializer::serialize(const std::string& key, const tgt::dvec4& data) {
    serializeTgtVector(key, data);
}

void JsonSerializer::serialize(const std::string& key, const tgt::ivec2& data) {
    serializeTgtVector(key, data);
}

void JsonSerializer::serialize(const std::string& key, const tgt::ivec3& data) {
    serializeTgtVector(key, data);
}

void JsonSerializer::serialize(const std::string& key, const tgt::ivec4& data) {
    serializeTgtVector(key, data);
}

void JsonSerializer::serialize(const std::string& key, const tgt::svec2& data) {
    serializeTgtVector(key, data);
}

void JsonSerializer::serialize(const std::string& key, const tgt::svec3& data) {
    serializeTgtVector(key, data);
}

void JsonSerializer::serialize(const std::string& key, const tgt::svec4& data) {
    serializeTgtVector(key, data);
}

void JsonSerializer::serialize(const std::string& key, const tgt::col3& data) {
    serializeTgtVector(key, data);
}

void JsonSerializer::serialize(const std::string& key, const tgt::col4& data) {
    serializeTgtVector(key, data);
}

void JsonSerializer::serialize(const std::string& key, const tgt::mat2& data) {
    serializeTgtVector(key+".row0", data[0]);
    serializeTgtVector(key+".row1", data[1]);
}

void JsonSerializer::serialize(const std::string& key, const tgt::mat3& data) {
    serializeTgtVector(key+".row0", data[0]);
    serializeTgtVector(key+".row1", data[1]);
    serializeTgtVector(key+".row2", data[2]);
}

void JsonSerializer::serialize(const std::string& key, const tgt::mat4& data) {
    serializeTgtVector(key+".row0", data[0]);
    serializeTgtVector(key+".row1", data[1]);
    serializeTgtVector(key+".row2", data[2]);
    serializeTgtVector(key+".row3", data[3]);
}

void JsonSerializer::serialize(const std::string& key, const tgt::Matrix2d& data) {
    serializeTgtVector(key+".row0", data[0]);
    serializeTgtVector(key+".row1", data[1]);
}

void JsonSerializer::serialize(const std::string& key, const tgt::Matrix3d& data) {
    serializeTgtVector(key+".row0", data[0]);
    serializeTgtVector(key+".row1", data[1]);
    serializeTgtVector(key+".row2", data[2]);
}

void JsonSerializer::serialize(const std::string& key, const tgt::Matrix4d& data) {
    serializeTgtVector(key+".row0", data[0]);
    serializeTgtVector(key+".row1", data[1]);
    serializeTgtVector(key+".row2", data[2]);
    serializeTgtVector(key+".row3", data[3]);
}

void JsonSerializer::serialize(const std::string& key, const Serializable& data) {
    serializeSimpleTypes(key, data);
}

rapidjson::Value JsonSerializer::serializeToValue(const Serializable& data) {
    rapidjson::Value val;
    val.SetObject();

    {
        TemporaryNodeChanger nodeChanger(*this, val);
        Serializer serializer(*this);
        data.serialize(serializer);
    }

    return val;
}

void JsonSerializer::write(std::ostream& stream, bool pretty, bool compressed) {
    using namespace boost::iostreams;

    resolveUnresolvedReferences();

    // Prepare gzip compressing stream
    // We have to do it in the outside scope (that lives longer than the OStreamWrapper sb
    // because sb references this stream if compression is used.
    filtering_ostream compressingStream;
    compressingStream.push(gzip_compressor(
                gzip_params(
                    gzip::best_compression
                    )
                ));

    // Wrap compressing stream for rapidjson
    std::unique_ptr<rapidjson::OStreamWrapper> sb;
    if(compressed) {
        compressingStream.push(stream);
        sb.reset(new rapidjson::OStreamWrapper(compressingStream));
    } else {
        sb.reset(new rapidjson::OStreamWrapper(stream));
    }

    if(pretty) {
        rapidjson::PrettyWriter<rapidjson::OStreamWrapper, rapidjson::UTF8<>, rapidjson::UTF8<>, rapidjson::CrtAllocator, rapidjson::kWriteNanAndInfFlag> writer(*sb);
        document_.Accept(writer);
    } else {
        rapidjson::Writer<rapidjson::OStreamWrapper, rapidjson::UTF8<>, rapidjson::UTF8<>, rapidjson::CrtAllocator, rapidjson::kWriteNanAndInfFlag> writer(*sb);
        document_.Accept(writer);
    }
}

} // namespace
