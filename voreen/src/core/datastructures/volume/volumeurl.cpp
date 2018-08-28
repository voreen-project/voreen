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

#include "voreen/core/datastructures/volume/volumeurl.h"

#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/utils/voreenfilepathhelper.h"

#include "tgt/filesystem.h"

#include <string>

namespace {
int lower_case(int c) {
    return tolower(c);
}
}

using std::string;

namespace voreen {

const std::string VolumeURL::loggerCat_("voreen.VolumeURL");

VolumeURL::VolumeURL()
    : protocol_("")
    , path_("")
{}

VolumeURL::VolumeURL(const VolumeURL& rhs)
    : Serializable()
{
    protocol_ = rhs.getProtocol();
    path_ = rhs.getPath();
    searchParameterMap_ = rhs.searchParameterMap_;

    metaDataContainer_ = rhs.getMetaDataContainer();
}

VolumeURL::VolumeURL(const std::string& url) {
    parseURL(url, protocol_, path_, searchParameterMap_);
    path_ = tgt::FileSystem::cleanupPath(path_);

    // fixes problems with symbolic links under linux
    if (tgt::FileSystem::isAbsolutePath(path_))
        path_ = tgt::FileSystem::absolutePath(path_);
}

VolumeURL::~VolumeURL() {
}

VolumeURL& VolumeURL::operator=(const VolumeURL& rhs) {
    protocol_ = rhs.getProtocol();
    path_ = rhs.getPath();
    searchParameterMap_ = rhs.searchParameterMap_;
    metaDataContainer_ = rhs.getMetaDataContainer();
    return *this;
}

bool VolumeURL::operator==(const VolumeURL& rhs) const {
    return (getURL() == rhs.getURL());
}


VolumeURL::VolumeURL(const std::string& protocol, const std::string& filepath, const std::string& searchString) {
    // TODO: validate parameters
    protocol_ = protocol;
    path_ = tgt::FileSystem::cleanupPath(filepath);
    // fixes problems with symbolic links under linux
    if (tgt::FileSystem::isAbsolutePath(path_))
        path_ = tgt::FileSystem::absolutePath(path_);

    searchParameterMap_ = parseSearchString(searchString);
}

void VolumeURL::serialize(Serializer& s) const {
    s.serialize("path", VoreenFilePathHelper(path_));
    s.serialize("protocol", protocol_);
    s.serialize("searchParameters", searchParameterMap_);

    // Determine if path is in relative format.
    // Since this is unlikely to happen, it is serialized optionally.
    if (!tgt::FileSystem::isAbsolutePath(tgt::FileSystem::cleanupPath(path_)))
        s.serialize("absolute", false);
}

void VolumeURL::deserialize(Deserializer& s) {

    try {
        VoreenFilePathHelper tmp;
        s.deserialize("path", tmp);
        path_ = tgt::FileSystem::cleanupPath(tmp.getPath());
        s.deserialize("protocol", protocol_);
        s.deserialize("searchParameters", searchParameterMap_);

        // Convert into relative path, if necessary.
        bool absolute;
        s.optionalDeserialize("absolute", absolute, true);
        if (!absolute && !s.getDocumentPath().empty())
            path_ = tgt::FileSystem::relativePath(path_, tgt::FileSystem::dirName(s.getDocumentPath()));

        return; // We are done here.
    } catch(SerializationNoSuchDataException&) {
        s.removeLastError();
    }

    // look for alternative attribute name 'url' (legacy)
    bool found = false;
    try {
        std::string url;
        s.deserialize("url", url);
        parseURL(url, protocol_, path_, searchParameterMap_);
        found = true;
    }
    catch (SerializationNoSuchDataException&) {
        s.removeLastError();
    }

    // look for alternative attribute name 'filename' (legacy)
    if (!found) {
        try {
            std::string url;
            s.deserialize("filename", url);
            parseURL(url, protocol_, path_, searchParameterMap_);
            found = true;
        }
        catch (SerializationNoSuchDataException&) {
            s.removeLastError();
        }
    }
    if (!found)
        throw SerializationNoSuchDataException("VolumeURL: neither attribute 'filename' nor attribute 'url' found");

    std::string basePath = tgt::FileSystem::dirName(s.getDocumentPath());
    if (!basePath.empty()) {
        VolumeURL originConv;
        try {
            originConv = VolumeSerializerPopulator().getVolumeSerializer()->convertOriginToAbsolutePath(*this, basePath);
        }
        catch (tgt::UnsupportedFormatException& e) {
            throw SerializationException(std::string(e.what()));
        }
        path_ = originConv.getPath();
    }

    path_ = tgt::FileSystem::cleanupPath(path_);
}

std::string VolumeURL::getURL() const {
    return constructURL(protocol_, path_, searchParameterMap_);
}

std::string VolumeURL::getPath() const {
    return path_;
}

std::string VolumeURL::getFilename() const {
    return tgt::FileSystem::fileName(getPath());
}

std::string VolumeURL::getSearchString() const {
    return constructSearchString(searchParameterMap_);
}

std::string VolumeURL::getProtocol() const {
    return protocol_;
}

void VolumeURL::addSearchParameter(const std::string& key, const std::string& value) {
    if (key.empty() || value.empty()) {
        LWARNING("Search key or value empty.");
        return;
    }

    searchParameterMap_[key] = value;
}

void VolumeURL::removeSearchParameter(const std::string& key) {
    if (key == "")
        return;

    std::map<std::string, std::string>::iterator it = searchParameterMap_.find(key);
    if (it != searchParameterMap_.end())
        searchParameterMap_.erase(it);
}

std::string VolumeURL::getSearchParameter(const std::string& k, bool caseSensitive) const {

    std::string key = k;
    if (caseSensitive) {
        std::map<std::string, std::string>::const_iterator it = searchParameterMap_.find(key);
        if (it != searchParameterMap_.end())
            return it->second;
        else
            return "";
    }
    else {
        transform(key.begin(), key.end(), key.begin(), lower_case);
        for (std::map<std::string, std::string>::const_iterator it = searchParameterMap_.begin();
                it != searchParameterMap_.end(); ++it) {
            std::string curKey = it->first;
            transform(curKey.begin(), curKey.end(), curKey.begin(), lower_case);
            if (key == curKey)
                return it->second;
        }
        return "";
    }
}

MetaDataContainer& VolumeURL::getMetaDataContainer() {
    return metaDataContainer_;
}

const MetaDataContainer& VolumeURL::getMetaDataContainer() const {
    return metaDataContainer_;
}

std::string VolumeURL::convertURLToRelativePath(const std::string& url, const std::string& basePath) {
    if (basePath.empty())
        return url;

    VolumeURL origin(url);
    VolumeURL originConv;
    try {
        originConv = VolumeSerializerPopulator().getVolumeSerializer()->convertOriginToRelativePath(origin, basePath);
    }
    catch (tgt::UnsupportedFormatException& e) {
        LWARNING(std::string(e.what()));
        originConv = origin;
    }

    // serialize with unix path separators for platform consistency
    std::string result = constructURL(originConv.getProtocol(), tgt::FileSystem::cleanupPath(originConv.getPath(), false), originConv.searchParameterMap_);
    return result;
}

std::string VolumeURL::convertURLToAbsolutePath(const std::string& url, const std::string& basePath) {
    if (basePath.empty())
        return url;

    VolumeURL origin(url);
    VolumeURL originConv;
    std::string result;
    try {
        originConv = VolumeSerializerPopulator().getVolumeSerializer()->convertOriginToAbsolutePath(origin, basePath);
        result = originConv.getURL();
    }
    catch (tgt::UnsupportedFormatException& e) {
        LWARNING(std::string(e.what()));
        result = url;
    }

    return result;
    //return tgt::FileSystem::cleanupPath(result);
}

std::string VolumeURL::constructURL(const std::string& protocol, const std::string& path, const std::map<std::string, std::string>& searchParameters) {
    std::string url = path;

    if (!protocol.empty())
        url = protocol + "://" + url;

    if (!searchParameters.empty())
        url += "?" + constructSearchString(searchParameters);

    return url;
}

std::string VolumeURL::constructSearchString(const std::map<std::string, std::string>& searchParameters) {
    std::string searchString;
    for (std::map<std::string, std::string>::const_iterator it = searchParameters.begin(); it != searchParameters.end(); ++it) {
        if (!searchString.empty())
            searchString += "&";
        searchString += escapeString(it->first) + "=" + escapeString(it->second);
    }
    return searchString;
}

void VolumeURL::parseURL(const std::string& url, std::string& protocol, std::string& path, std::map<std::string, std::string>& searchParameters) {
    // protocol
    string::size_type sep_pos = url.find("://");
    if (sep_pos == std::string::npos) {
        // URL does not contain protocol specifier
        protocol = "";
    }
    else {
        // return substring before protocol separator
        protocol = url.substr(0, sep_pos);
    }


    // path and searchString
    std::string fullPath;
    if (sep_pos == std::string::npos) {
        fullPath = url;
    }
    else {
        fullPath = url.substr(sep_pos + 3);
    }

    // get search string
    sep_pos = fullPath.find("?");
    std::string searchString;
    if (sep_pos == std::string::npos) {
        // URL does not contain search string
        path = fullPath;
        searchString = "";
    }
    else {
        // separate path from searchstring
        path = fullPath.substr(0, sep_pos);
        searchString = fullPath.substr(sep_pos + 1);
    }

    searchParameters = parseSearchString(searchString);
}

std::map<std::string, std::string> VolumeURL::parseSearchString(std::string searchString)  {
    std::map<std::string, std::string> searchParameters;

    // temporarily replace escaped '&' by '#1#' and '=' by '#2#' in order to allow splitting
    // of the search string by these chars
    searchString = strReplaceAll(searchString, "\\&", "#1#");
    searchString = strReplaceAll(searchString, "\\=", "#2#");

    std::vector<std::string> searchElems = strSplit(searchString, '&');
    for (size_t i=0; i<searchElems.size(); i++) {
        std::vector<std::string> keyValuePair = strSplit(searchElems.at(i), '=');
        if (keyValuePair.size() == 2) {
            std::string key = keyValuePair.at(0);
            std::string value = keyValuePair.at(1);

            key = strReplaceAll(key, "#1#", "\\&");
            key = strReplaceAll(key, "#2#", "\\=");
            key = unescapeString(key);

            value = strReplaceAll(value, "#1#", "\\&");
            value = strReplaceAll(value, "#2#", "\\=");
            value = unescapeString(value);

            searchParameters.insert(std::make_pair(key, value));
        }
    }

    return searchParameters;
}

std::string VolumeURL::escapeString(const std::string& str) {
    std::string result = str;

    result = strReplaceAll(result, "\\", "\\\\");
    result = strReplaceAll(result, "?", "\\?");
    result = strReplaceAll(result, "&", "\\&");
    result = strReplaceAll(result, "=", "\\=");
    result = strReplaceAll(result, " ", "\\ ");

    return result;
}

std::string VolumeURL::unescapeString(const std::string& str) {
    std::string result = str;

    result = strReplaceAll(result, "\\?", "?");
    result = strReplaceAll(result, "\\&", "&");
    result = strReplaceAll(result, "\\=", "=");
    result = strReplaceAll(result, "\\ ", " ");
    result = strReplaceAll(result, "\\\\", "\\");

    return result;
}

} // namespace
