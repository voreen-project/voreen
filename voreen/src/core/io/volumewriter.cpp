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

#include "voreen/core/io/volumewriter.h"
#include "voreen/core/io/progressbar.h"

namespace voreen {

const std::string VolumeWriter::loggerCat_("voreen.io.VolumeWriter");

VolumeWriter::VolumeWriter(ProgressBar* progress)
    : progress_(progress)
{}

const std::vector<std::string>& VolumeWriter::getSupportedExtensions() const {
    return extensions_;
}

const std::vector<std::string>& VolumeWriter::getSupportedFilenames() const {
    return filenames_;
}

const std::vector<std::string>& VolumeWriter::getSupportedProtocols() const {
    return protocols_;
}

std::string VolumeWriter::getFileNameWithoutExtension(const std::string& filename) {
    return filename.substr(0, filename.rfind("."));
}

std::string VolumeWriter::getExtension(const std::string& filename) {
    return filename.substr(filename.rfind(".") + 1, filename.length());
}

void VolumeWriter::setProgressBar(ProgressBar* progressBar) {
    progress_ = progressBar;
}

ProgressBar* VolumeWriter::getProgressBar() const {
    return progress_;
}



} // namespace voreen
