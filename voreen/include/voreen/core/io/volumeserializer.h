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

#ifndef VRN_VOLUMESERIALIZER_H_
#define VRN_VOLUMESERIALIZER_H_

#include <map>
#include <vector>
#include <string>

#include "tgt/exception.h"

#include "voreen/core/voreencoreapi.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumelist.h"

namespace voreen {

// forward declarations
class Volume;
class VolumeReader;
class VolumeWriter;
class VolumeSerializer;
class ProgressBar;

/**
 * You can register several VolumeReader and VolumeWriter instances in this class.
 * A call to \a read or \a write will then automatically select the proper VolumeReader /
 * VolumeWriter and delegate the loading / saving to it.
 */
class VRN_CORE_API VolumeSerializer {
public:

    VolumeSerializer();
    ~VolumeSerializer();

    /**
     * Finds suitable VolumeReaders for the specified URL.
     * If an empty string is passed, all registered readers are returned.
     *
     * @param url URL of the file to load. May be a plain filename.
     * @return VolumeReaders that are able to read from the URL.
     *
     * @throw UnsupportedFormatException if no suitable reader has been found.
     */
    std::vector<VolumeReader*> getReaders(const std::string& url = "") const;

    /**
     * Returns the reader with the passed className, or 0 if
     * no reader with that className is available.
     */
    VolumeReader* getReaderByName(const std::string& className) const;

    /**
     * Finds suitable VolumeWriters for the specified URL.
     * If an empty string is passed, all registered writers are returned.
     *
     * @param url URL of the file to write. May be a plain filename.
     * @return VolumeWriters that are able to read from the filename.
     *
     * @throw UnsupportedFormatException if no suitable writer has been found.
     */
    std::vector<VolumeWriter*> getWriters(const std::string& url = "") const;

    /**
     * Returns the writer with the passed className, or 0 if
     * no writer with that className is available.
     */
    VolumeWriter* getWriterByName(const std::string& className) const;

    /**
     * Loads one or multiple volumes from the specified URL.
     *
     * @param url The URL the data set is to be read from.
     *      This may be a file path, e.g. /file/to/volume.dat,
     *      or a complete URL with resource type and inner path, e.g.
     *      zip://path/to/archive.zip/volume.dat
     * @return VolumeList containing all volumes read from the url.
     *      the caller is responsible for freeing the memory.
     */
    VolumeList* read(const std::string& url) const;

    /**
     * Instead of reading the whole dataset, only some slices are read from file and written
     * to the newly built volume. This isn't supported by all file formats.
     *
     * @return new VolumeList, the caller is responsible for freeing the memory.
     */
    VolumeList* readSlices(const std::string& url, size_t firstSlice, size_t lastSlice) const;

    /**
     * Instead of reading the whole dataset, only a brick of volume data, specified by
     * its starting location and its dimensions, is read. This isn't supported by all file formats.
     *
     * @return new VolumeList, the caller is responsible for freeing the memory.
     */
    VolumeList* readBrick(const std::string& url, tgt::ivec3 brickStartPos, int brickSize) const;

    /**
     * Reads a single volume from the specified origin.
     *
     * This function is mainly used for the deserialization
     * of VolumeHandles.
     *
     * @return Volume encapsulating the loaded volume.
     *      The caller is responsible for freeing the memory.
     */
    VolumeBase* read(const VolumeURL& origin) const;

    /**
     * Converts the file path of the passed origin to a path relative to the passed base path.
     * The conversion is delegated to the appropriate reader.
     *
     * @return origin with relative file path. If the passed origin already has a relative path or
     *      if it does not contain a file path at all, an identical copy is returned.
     */
    VolumeURL convertOriginToRelativePath(const VolumeURL& origin, const std::string& basePath) const;
    /**
     * Converts the passed origin's file path from a path relative to the passed base path
     * to an absolute path. The conversion is delegated to the appropriate reader.
     *
     * @return origin with absolute file path. If the passed origin already has an absolute path or
     *      if it does not contain a file path at all, an identical copy is returned.
     */
    VolumeURL convertOriginToAbsolutePath(const VolumeURL& origin, const std::string& basePath) const;

    /**
     * Saves a Volume to the given file.
     *
     * @param url The URL of the file where the data should be written to. May be a plain filename.
     * @param volume The Volume that should be saved.
     */
    void write(const std::string& url, const VolumeBase* volume) const;

    /**
     * Use this method to register a VolumeReader.
     */
    void registerReader(VolumeReader* vr);

    /**
     * Use this method to register a VolumeWriter.
     */
    void registerWriter(VolumeWriter* vw);

    /**
     * Assigns the passed progress bar to all registered readers and writers.
     * May be null.
     */
    void setProgressBar(ProgressBar* progressBar);

private:
    void appendPreferredReaderToOriginURLs(const VolumeList* collection, const VolumeReader* volumeReader) const;
    void appendPreferredReaderToOriginURLs(VolumeBase* handle, const VolumeReader* volumeReader) const;

    /// all registered readers
    std::vector<VolumeReader*> readers_;

    /// all registered readers
    std::vector<VolumeWriter*> writers_;

    /// maps from filename extensions to the appropriate reader
    std::map<std::string, std::vector<VolumeReader*> > readersExtensionMap_;

    /// maps from filenames to the appropriate reader
    std::map<std::string, std::vector<VolumeReader*> > readersFilenameMap_;

    /// maps from protocol strings (e.g. zip, dicom) to the appropriate reader
    std::map<std::string, std::vector<VolumeReader*> > readersProtocolMap_;

    /// maps from filename extensions to the appropriate writer
    std::map<std::string, std::vector<VolumeWriter*> > writersExtensionMap_;

    /// maps from filenames to the appropriate writer
    std::map<std::string, std::vector<VolumeWriter*> > writersFilenameMap_;

    /// maps from protocol strings (e.g. zip, dicom) to the appropriate writer
    std::map<std::string, std::vector<VolumeWriter*> > writersProtocolMap_;

    static const std::string loggerCat_;

};

} // namespace voreen

#endif // VRN_VOLUMESERIALIZER_H_
