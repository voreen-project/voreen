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

#include "voreen/core/io/volumeserializer.h"

#include <cstring>

#include "voreen/core/io/volumereader.h"
#include "voreen/core/io/volumewriter.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/utils/stringutils.h"
#include "tgt/filesystem.h"

namespace voreen {

//------------------------------------------------------------------------------

const std::string VolumeSerializer::loggerCat_("voreen.VolumeSerializer");

VolumeSerializer::VolumeSerializer() {
}

VolumeSerializer::~VolumeSerializer() {
}

std::vector<VolumeReader*> VolumeSerializer::getReaders(const std::string& url) const {
    if (url.empty())
        return readers_;

    VolumeURL origin(url);

    // check if a reader for the URL's resource type is available
    std::string protocol = origin.getProtocol();
    if (!protocol.empty()) {
        if (readersProtocolMap_.find(protocol) != readersProtocolMap_.end())
            return readersProtocolMap_.find(protocol)->second;
        else
            throw tgt::UnsupportedFormatException(protocol, url);
    }

    // check if a reader for the filename is available
    std::string filename = tgt::FileSystem::fileName(origin.getPath());
    if (readersFilenameMap_.find(filename) != readersFilenameMap_.end())
        return readersFilenameMap_.find(filename)->second;

    // check if a reader for the filename extension is available
    std::string extension = tgt::FileSystem::fileExtension(origin.getPath(), true);

    // If no file extension is present, check for readers which could propably handle them.
    if (extension.empty()) {
        std::vector<VolumeReader*> result;
        for (VolumeReader* reader : readers_) {
            if (reader->canHandleMissingExtension())
                result.push_back(reader);
        }
        return result;
    }

    // construct extension also from last two filename components, e.g. ome.tiff
    std::string doubleExtension;
    std::vector<std::string> filenameSplit = strSplit(origin.getFilename(), ".");
    if (filenameSplit.size() > 2)
        doubleExtension = filenameSplit.at(filenameSplit.size()-2) + "." + filenameSplit.at(filenameSplit.size()-1);

    if (readersExtensionMap_.find(doubleExtension) != readersExtensionMap_.end()) //< prefer double extension over single extension
        return readersExtensionMap_.find(doubleExtension)->second;
    if (readersExtensionMap_.find(extension) != readersExtensionMap_.end())
        return readersExtensionMap_.find(extension)->second;
    else
        throw tgt::UnsupportedFormatException(extension, url);
}

VolumeReader* VolumeSerializer::getReaderByName(const std::string& className) const {
    for (size_t i=0; i<readers_.size(); i++) {
        if (readers_.at(i)->getClassName() == className)
            return readers_.at(i);
    }
    return 0;
}

std::vector<VolumeWriter*> VolumeSerializer::getWriters(const std::string& url) const {
    if (url.empty())
        return writers_;

    VolumeURL origin(url);

    // check if a writer for the resource type is available
    std::string protocol = origin.getProtocol();
    if (!protocol.empty()) {
        if (writersProtocolMap_.find(protocol) != writersProtocolMap_.end())
            return writersProtocolMap_.find(protocol)->second;
        else
            throw tgt::UnsupportedFormatException(protocol, url);
    }

    // check if a reader for the filename is available
    std::string filename = tgt::FileSystem::fileName(origin.getPath());
    if (writersFilenameMap_.find(filename) != writersFilenameMap_.end())
        return writersFilenameMap_.find(filename)->second;

    // check if a reader for the filename extension is available
    std::string extension = tgt::FileSystem::fileExtension(origin.getPath(), true);
    if (writersExtensionMap_.find(extension) != writersExtensionMap_.end())
        return writersExtensionMap_.find(extension)->second;
    else
        throw tgt::UnsupportedFormatException(extension, url);
}

VolumeWriter* VolumeSerializer::getWriterByName(const std::string& className) const {
    for (size_t i=0; i<writers_.size(); i++) {
        if (writers_.at(i)->getClassName() == className)
            return writers_.at(i);
    }
    return 0;
}

VolumeList* VolumeSerializer::read(const std::string& url) const {
    std::vector<VolumeReader*> matchingReaders = getReaders(url);
    tgtAssert(!matchingReaders.empty(), "readers vector empty (exception expected)");
    if (matchingReaders.size() == 1) {
        VolumeList* collection = matchingReaders.front()->read(url);
        if (collection)
            appendPreferredReaderToOriginURLs(collection, matchingReaders.front());
        return collection;
    }
    else {
         std::vector<std::string> errors;

        // if preferred reader is specified and available, use it for loading the data set
        VolumeReader* prefReader = 0;
        VolumeURL urlOrigin(url);
        std::string prefReaderStr = urlOrigin.getSearchParameter("preferredReader");
        if (prefReaderStr != "")
            prefReader = getReaderByName(prefReaderStr);
        if (prefReader) {
            try {
                VolumeList* collection = prefReader->read(url);
                if (collection)
                    appendPreferredReaderToOriginURLs(collection,  prefReader);
                return collection;
            }
            catch (const tgt::FileException& e) {
                errors.push_back(e.what());
            }
        }

        // if no preferred reader or preferred reader could not load the data set,
        // iterate over all possibly matching readers and try to load data set, collect error messages
        for (size_t i=0; i<matchingReaders.size(); i++) {
            if (matchingReaders.at(i) == prefReader)
                continue;
            try {
                VolumeList* collection = matchingReaders.at(i)->read(url);
                if (collection)
                    appendPreferredReaderToOriginURLs(collection,  matchingReaders.at(i));
                return collection;
            }
            catch (const tgt::FileException& e) {
                errors.push_back(e.what());
            }
        }
        std::string errorMsg = "No VolumeReader was able to read the data set: " + strJoin(errors, "; ");
        throw tgt::FileException(errorMsg);

    }
}

VolumeList* VolumeSerializer::readSlices(const std::string& url, size_t firstSlice, size_t lastSlice) const {
    std::vector<VolumeReader*> matchingReaders = getReaders(url);
    tgtAssert(!matchingReaders.empty(), "readers vector empty (exception expected)");
    if (matchingReaders.size() == 1) {
        VolumeList* collection = matchingReaders.front()->readSlices(url, firstSlice, lastSlice);
        if (collection)
            appendPreferredReaderToOriginURLs(collection, matchingReaders.front());
        return collection;
    }
    else {
        // iterate over all possibly matching readers try to load data set, collect error messages
        std::vector<std::string> errors;
        for (size_t i=0; i<matchingReaders.size(); i++) {
            try {
                VolumeList* collection = matchingReaders.at(i)->readSlices(url, firstSlice, lastSlice);
                if (collection)
                    appendPreferredReaderToOriginURLs(collection, matchingReaders.at(i));
                return collection;
            }
            catch (const tgt::FileException& e) {
                errors.push_back(e.what());
            }
        }
        std::string errorMsg = "No VolumeReader was able to read the data set: " + strJoin(errors, "; ");
        throw tgt::FileException(errorMsg);
    }
}

VolumeList* VolumeSerializer::readBrick(const std::string& url, tgt::ivec3 brickStartPos,
    int brickSize) const {
    std::vector<VolumeReader*> matchingReaders = getReaders(url);
    tgtAssert(!matchingReaders.empty(), "readers vector empty (exception expected)");
    if (matchingReaders.size() == 1) {
        VolumeList* collection = matchingReaders.front()->readBrick(url, brickStartPos, brickSize);
        if (collection)
            appendPreferredReaderToOriginURLs(collection, matchingReaders.front());
        return collection;
    }
    else {
        // iterate over all possibly matching readers try to load data set, collect error messages
        std::vector<std::string> errors;
        for (size_t i=0; i<matchingReaders.size(); i++) {
            try {
                VolumeList* collection = matchingReaders.at(i)->readBrick(url, brickStartPos, brickSize);
                if (collection)
                    appendPreferredReaderToOriginURLs(collection, matchingReaders.at(i));
                return collection;
            }
            catch (const tgt::FileException& e) {
                errors.push_back(e.what());
            }
        }
        std::string errorMsg = "No VolumeReader was able to read the data set: " + strJoin(errors, "; ");
        throw tgt::FileException(errorMsg);
    }
}

VolumeBase* VolumeSerializer::read(const VolumeURL& origin) const {
    std::vector<VolumeReader*> matchingReaders = getReaders(origin.getURL());
    tgtAssert(!matchingReaders.empty(), "readers vector empty (exception expected)");
    if (matchingReaders.size() == 1) {
        VolumeBase* handle = matchingReaders.front()->read(origin);
        appendPreferredReaderToOriginURLs(handle, matchingReaders.front());
        return handle;
    }
    else {
        std::vector<std::string> errors;

        // if preferred reader is specified and available, use it for loading the data set
        VolumeReader* prefReader = 0;
        std::string prefReaderStr = origin.getSearchParameter("preferredReader");
        if (prefReaderStr != "")
            prefReader = getReaderByName(prefReaderStr);
        if (prefReader) {
            try {
                VolumeBase* handle = prefReader->read(origin);
                if (handle)
                    appendPreferredReaderToOriginURLs(handle,  prefReader);
                return handle;
            }
            catch (const tgt::FileException& e) {
                errors.push_back(e.what());
            }
        }

        // if no preferred reader or preferred reader could not load the data set,
        // iterate over all possibly matching readers and try to load data set, collect error messages
        for (size_t i=0; i<matchingReaders.size(); i++) {
            if (matchingReaders.at(i) == prefReader)
                continue;
            try {
                VolumeBase* handle = matchingReaders.at(i)->read(origin);
                if (handle)
                    appendPreferredReaderToOriginURLs(handle,  matchingReaders.at(i));
                return handle;
            }
            catch (const tgt::FileException& e) {
                errors.push_back(e.what());
            }
        }
        std::string errorMsg = "No VolumeReader was able to read the data set: " + strJoin(errors, "; ");
        throw tgt::FileException(errorMsg);
    }
}

void VolumeSerializer::write(const std::string& url, const VolumeBase* volumeHandle) const {
    std::vector<VolumeWriter*> matchingWriters = getWriters(url);
    tgtAssert(!matchingWriters.empty(), "writers vector empty (exception expected)");
    if (matchingWriters.size() == 1) {
        return matchingWriters.front()->write(url, volumeHandle);
    }
    else {
        // iterate over all possibly matching writers try to load data set, collect error messages
        std::vector<std::string> errors;
        for (size_t i=0; i<matchingWriters.size(); i++) {
            try {
                return matchingWriters.at(i)->write(url, volumeHandle);
            }
            catch (const tgt::FileException& e) {
                errors.push_back(e.what());
            }
        }
        std::string errorMsg = "No VolumeWriter was able to write the data set: " + strJoin(errors, "; ");
        throw tgt::FileException(errorMsg);
    }
}

VolumeURL VolumeSerializer::convertOriginToRelativePath(const VolumeURL& origin, const std::string& basePath) const {
    std::vector<VolumeReader*> readers = getReaders(origin.getURL());
    tgtAssert(!readers.empty(), "readers vector must not be empty");
    return readers.front()->convertOriginToRelativePath(origin, basePath); // TODO: how to to select appropriate reader?
}

VolumeURL VolumeSerializer::convertOriginToAbsolutePath(const VolumeURL& origin, const std::string& basePath) const {
    std::vector<VolumeReader*> readers = getReaders(origin.getURL());
    tgtAssert(!readers.empty(), "readers vector must not be empty");
    return readers.front()->convertOriginToAbsolutePath(origin, basePath); // TODO: how to to select appropriate reader?
}

void VolumeSerializer::registerReader(VolumeReader* vr) {
    tgtAssert(vr, "null pointer passed");
    readers_.push_back(vr);

    // update extensions map
    const std::vector<std::string>& extensions = vr->getSupportedExtensions();
    for (size_t i = 0; i < extensions.size(); ++i) {
        // insert new entry (empty vector), if extension is not already mapped
        if (readersExtensionMap_.find(extensions[i]) == readersExtensionMap_.end()) {
            std::vector<VolumeReader*> readersVec;
            readersExtensionMap_.insert(std::make_pair(extensions[i], readersVec));
        }
        // add reader to mapped vector
        tgtAssert(readersExtensionMap_.find(extensions[i]) != readersExtensionMap_.end(), "map entry should exist here");
        readersExtensionMap_.find(extensions[i])->second.push_back(vr);
    }

    // update filename map
    const std::vector<std::string>& filenames = vr->getSupportedFilenames();
    for (size_t i = 0; i < filenames.size(); ++i) {
        // insert new entry (empty vector), if filename is not already mapped
        if (readersFilenameMap_.find(filenames[i]) == readersFilenameMap_.end()) {
            std::vector<VolumeReader*> readersVec;
            readersFilenameMap_.insert(std::make_pair(filenames[i], readersVec));
        }
        // add reader to mapped vector
        tgtAssert(readersFilenameMap_.find(filenames[i]) != readersFilenameMap_.end(), "map entry should exist here");
        readersFilenameMap_.find(filenames[i])->second.push_back(vr);
    }

    // update protocols map
    const std::vector<std::string>& formats = vr->getSupportedProtocols();
    for (size_t i = 0; i < formats.size(); ++i) {
        // insert new entry (empty vector), if extension is not already mapped
        if (readersProtocolMap_.find(formats[i]) == readersProtocolMap_.end()) {
            std::vector<VolumeReader*> readersVec;
            readersProtocolMap_.insert(std::make_pair(formats[i], readersVec));
        }
        // add reader to mapped vector
        tgtAssert(readersProtocolMap_.find(formats[i]) != readersProtocolMap_.end(), "map entry should exist here");
        readersProtocolMap_.find(formats[i])->second.push_back(vr);
    }
}

void VolumeSerializer::registerWriter(VolumeWriter* vw) {
    tgtAssert(vw, "null pointer passed");
    writers_.push_back(vw);

    // update extensions map
    const std::vector<std::string>& extensions = vw->getSupportedExtensions();
    for (size_t i = 0; i < extensions.size(); ++i) {
        // insert new entry (empty vector), if extension is not already mapped
        if (writersExtensionMap_.find(extensions[i]) == writersExtensionMap_.end()) {
            std::vector<VolumeWriter*> writersVec;
            writersExtensionMap_.insert(std::make_pair(extensions[i], writersVec));
        }
        // add writer to mapped vector
        tgtAssert(writersExtensionMap_.find(extensions[i]) != writersExtensionMap_.end(), "map entry should exist here");
        writersExtensionMap_.find(extensions[i])->second.push_back(vw);
    }

    // update filename map
    const std::vector<std::string>& filenames = vw->getSupportedFilenames();
    for (size_t i = 0; i < filenames.size(); ++i) {
        // insert new entry (empty vector), if filename is not already mapped
        if (writersFilenameMap_.find(filenames[i]) == writersFilenameMap_.end()) {
            std::vector<VolumeWriter*> writersVec;
            writersFilenameMap_.insert(std::make_pair(filenames[i], writersVec));
        }
        // add writer to mapped vector
        tgtAssert(writersFilenameMap_.find(filenames[i]) != writersFilenameMap_.end(), "map entry should exist here");
        writersFilenameMap_.find(filenames[i])->second.push_back(vw);
    }

    // update protocols map
    const std::vector<std::string>& formats = vw->getSupportedProtocols();
    for (size_t i = 0; i < formats.size(); ++i) {
        // insert new entry (empty vector), if protocol is not already mapped
        if (writersProtocolMap_.find(formats[i]) == writersProtocolMap_.end()) {
            std::vector<VolumeWriter*> writersVec;
            writersProtocolMap_.insert(std::make_pair(formats[i], writersVec));
        }
        // add writer to mapped vector
        tgtAssert(writersProtocolMap_.find(formats[i]) != writersProtocolMap_.end(), "map entry should exist here");
        writersProtocolMap_.find(formats[i])->second.push_back(vw);
    }
}

void VolumeSerializer::setProgressBar(ProgressBar* progressBar) {
    for (size_t i=0; i<readers_.size(); i++)
        readers_.at(i)->setProgressBar(progressBar);

    for (size_t i=0; i<writers_.size(); i++)
        writers_.at(i)->setProgressBar(progressBar);
}

void VolumeSerializer::appendPreferredReaderToOriginURLs(const VolumeList* collection, const VolumeReader* volumeReader) const {
    tgtAssert(collection, "null pointer passed");
    tgtAssert(volumeReader, "null pointer passed");

    for (size_t i=0; i<collection->size(); i++) {
        if (collection->at(i) == 0) {
            LERROR("collection contains null pointer");
            return;
        }
        appendPreferredReaderToOriginURLs(collection->at(i), volumeReader);
    }
}

void VolumeSerializer::appendPreferredReaderToOriginURLs(VolumeBase* handle, const VolumeReader* volumeReader) const {
    tgtAssert(handle, "null pointer passed");
    tgtAssert(volumeReader, "null pointer passed");

    handle->getOrigin().addSearchParameter("preferredReader", volumeReader->getClassName());
}


} // namespace voreen
