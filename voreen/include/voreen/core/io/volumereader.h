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

#ifndef VRN_VOLUMEREADER_H
#define VRN_VOLUMEREADER_H

#include <string>
#include <vector>

#include "voreen/core/voreencoreapi.h"
#include "voreen/core/datastructures/volume/volumelist.h"
#include "voreen/core/io/progressbar.h"

#include "tgt/exception.h"

namespace voreen {

// forward declarations
class Volume;

/**
 * Reads one or multiple volume data sets from a passed location.
 * Implement this class in order to support a new volume data format.
 */
class VRN_CORE_API VolumeReader {
public:
    VolumeReader(ProgressBar* progress = 0);
    virtual ~VolumeReader() {}

    /**
     * Returns the name of this class as a string.
     * Necessary due to the lack of code reflection in C++.
     */
    virtual std::string getClassName() const = 0;

    /**
     * Returns a description of the file format the reader supports.
     *
     * @note The description is to be shown in a file dialog
     *      and should therefore be short.
     */
    virtual std::string getFormatDescription() const = 0;

    /**
     * Virtual constructor.
     */
    virtual VolumeReader* create(ProgressBar* progress = 0) const = 0;

    /**
     * Loads one or multiple volumes from the specified URL.
     *
     * @param url The URL the data set is to be read from.
     *      This may be a file path, e.g. /file/to/volume.dat,
     *      or a complete URL with resource type and inner path, e.g.
     *      zip://path/to/archive.zip/volume.dat
     *
     * @return VolumeList containing all volumes read from the url.
     *      the caller is responsible for freeing the memory.
     *
     * @throw tgt::FileException if the data set could not be loaded
     */
    virtual VolumeList* read(const std::string& url) = 0;

    /**
     * Instead of reading a whole dataset, only some slices are read from file and written
     * to the newly built volume.
     *
     * Override this function in order to provide a brick-wise loading routine.
     * The default implementation throws an exception.
     *
     * @throw tgt::FileException if the data set could not be loaded
     */
    virtual VolumeList* readSlices(const std::string& url, size_t firstSlice = 0, size_t lastSlice = 0);

    /**
     * Instead of reading the whole dataset, only a brick of volume data, specified by
     * its starting location and its dimensions, is read.
     *
     * Override this function in order to provide a brick-wise loading routine.
     * The default implementation throws an exception.
     *
     * @throw tgt::FileException if the data set could not be loaded
     */
    virtual VolumeList* readBrick(const std::string& url, tgt::ivec3 start, int dimensions);

    /**
     * Loads a single volume from the passed origin.
     *
     * The default implementation calls the read() function with the
     * origin's path. Override it in order to access the VolumeURL directly.
     *
     * @param origin The origin the data set is to be read from.
     *
     * @return Volume encapsulating the loaded volume. The caller is responsible for freeing the memory.
     *
     * @throw tgt::FileException if the data set could not be loaded
     */
    virtual VolumeBase* read(const VolumeURL& origin);

    /**
     * Returns a list of all volumes that are found at the passed location. It is not guaranteed
     * that the listed volumes can actually be loaded by the reader.
     *
     * The default implementation returns a single VolumeURL encapsulating the passed URL,
     * if the referenced file exists, or an empty vector otherwise. Subclasses that support
     * container files, i.e. files that may store multiple volumes, should override this function.
     * Information that is not necessary for distinctly identifying a single volume but may still
     * be helpful for the user should be stored in the MetaDataContainer of the respective VolumeURL.
     */
    virtual std::vector<VolumeURL> listVolumes(const std::string& url) const;

    /**
     * Converts the file path of the passed origin to a path relative to the passed base path.
     *
     * @note The default implementation performs the conversion under the assumption that the origin's path
     *      refers to a physical file. Override this function in order to perform a custom conversion.
     *
     * @return origin with relative file path. If the passed origin already has a relative path or
     *      if it does not contain a file path at all, an identical copy is returned.
     */
    virtual VolumeURL convertOriginToRelativePath(const VolumeURL& origin, const std::string& basePath) const;

    /**
     * Converts the passed origin's file path from a path relative to the passed base path
     * to an absolute path.
     *
     * @note The default implementation performs the conversion under the assumption that the origin's path
     *      refers to a physical file. Override this function in order to perform a custom conversion.
     *
     * @return origin with absolute file path. If the passed origin already has an absolute path or
     *      if it does not contain a file path at all, an identical copy is returned.
     */
    virtual VolumeURL convertOriginToAbsolutePath(const VolumeURL& origin, const std::string& basePath) const;

    /**
     * Returns true if the reader is able to read files without extension, false otherwise.
     */
    virtual bool canHandleMissingExtension() const;

    /**
     * Returns true if the format supported by this reader supports file watching.
     */
    virtual bool canSupportFileWatching() const;

    /**
     * Returns a list of filename extensions that are supported by this reader.
     */
    const std::vector<std::string>& getSupportedExtensions() const;

    /**
     * Returns a list of filenames that are supported by this reader.
     */
    const std::vector<std::string>& getSupportedFilenames() const;

    /**
     * Returns a list of protocols that are supported by this reader.
     */
    const std::vector<std::string>& getSupportedProtocols() const;

    /**
     * Assigns a progress bar to the reader. May be null.
     */
    void setProgressBar(ProgressBar* progressBar);

    /**
     * Returns the assigned progress bar. May be null.
     */
    ProgressBar* getProgressBar() const;

protected:
    void read(VolumeRAM* volume, FILE* fin);

    /**
     * Reverses the order of the slice in x-direction. This method
     * is called, when the .dat file contains "SliceOrder: -x" in one line.
     */
    void reverseXSliceOrder(VolumeRAM* const volume) const;

    /**
     * Reverses the order of the slice in y-direction. This method
     * is called, when the .dat file contains "SliceOrder: -y" in one line.
     */
    void reverseYSliceOrder(VolumeRAM* const volume) const;

    /**
     * Reverses the order of the slice in z-direction. This method
     * is called, when the .dat file contains "SliceOrder: -z" in one line.
     */
    void reverseZSliceOrder(VolumeRAM* const volume) const;

    /// List of filename extensions supported by the reader.
    std::vector<std::string> extensions_;

    /// List of filenames supported by the reader.
    std::vector<std::string> filenames_;

    /// List of protocols supported by the reader.
    std::vector<std::string> protocols_;

    static const std::string loggerCat_;

private:
    ProgressBar* progress_;
};

} // namespace voreen

#endif // VRN_VOLUMEREADER_H
