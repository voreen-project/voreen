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

#ifndef VRN_VOLUMEWRITER_H
#define VRN_VOLUMEWRITER_H

#include "voreen/core/voreencoreapi.h"
#include <string>
#include <vector>

#include "tgt/exception.h"

namespace voreen {

// forward declarations
class ProgressBar;
class VolumeBase;

/**
 * Writes a volume dataset.
 * Implement this class in order to support a new format.
 */
class VRN_CORE_API VolumeWriter {
public:
    VolumeWriter(ProgressBar* progress = 0);
    virtual ~VolumeWriter() {}

    /**
     * Virtual constructor.
     */
    virtual VolumeWriter* create(ProgressBar* progress = 0) const = 0;

    /**
     * Returns the name of this class as a string.
     * Necessary due to the lack of code reflection in C++.
     */
    virtual std::string getClassName() const = 0;

    /**
     * Returns a description of the file format the writer supports.
     *
     * @note The description is to be shown in a file dialog
     *      and should therefore be short.
     */
    virtual std::string getFormatDescription() const = 0;

    /**
     * Saves a Volume to the given file.
     *
     * @param fileName The file name where the data should go.
     * @param Volume The volume which should be saved.
     */
    virtual void write(const std::string& fileName, const VolumeBase* volumeHandle) = 0;

    /**
     * Returns the filename extensions that are supported by the writer.
     */
    const std::vector<std::string>& getSupportedExtensions() const;

    /**
     * Returns a list of filenames that are supported by this writer.
     */
    const std::vector<std::string>& getSupportedFilenames() const;

    /**
     * Returns a list of protocols that are supported by this writer.
     */
    const std::vector<std::string>& getSupportedProtocols() const;

    /**
     * Use this method as a helper to get a file name without its extension
     *
     * @param fileName File name to be examined.
     * @return File name without extension.
     */
    static std::string getFileNameWithoutExtension(const std::string& filename);

    /**
     * Use this method as a helper to get the extension of a given file name.
     *
     * @param fileName File name to be examined.
     * @return The extension of the file name.
     */
    static std::string getExtension(const std::string& filename);

    /**
     * Assigns a progress bar to the writer. May be null.
     */
    void setProgressBar(ProgressBar* progressBar);

    /**
     * Returns the assigned progress bar. May be null.
     */
    ProgressBar* getProgressBar() const;

protected:
    /// List of filename extensions supported by the writer.
    std::vector<std::string> extensions_;

    /// List of filenames supported by the writer.
    std::vector<std::string> filenames_;

    /// List of protocols supported by the writer.
    std::vector<std::string> protocols_;

    static const std::string loggerCat_;

private:
    ProgressBar* progress_;

};

} // namespace voreen

#endif // VRN_VOLUMEWRITER_H
