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

#include "devilmodule.h"

#include <IL/il.h>
#include <IL/ilu.h>

namespace voreen {

std::vector<std::string> DevILModule::readExtensions_;
std::vector<std::string> DevILModule::writeExtensions_;

DevILModule::DevILModule(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("DevIL");
    setGuiName("DevIL");

    readExtensions_.push_back("png");
    readExtensions_.push_back("jpg");
    readExtensions_.push_back("jpeg");
    readExtensions_.push_back("tif");
    readExtensions_.push_back("tiff");
    readExtensions_.push_back("bmp");
    readExtensions_.push_back("gif");
    readExtensions_.push_back("tga");
    readExtensions_.push_back("pnm");

    writeExtensions_.push_back("png");
    writeExtensions_.push_back("jpg");
    writeExtensions_.push_back("tif");
    writeExtensions_.push_back("bmp");
}

std::string DevILModule::getDevILError() {
    ILenum error = ilGetError();
    return std::string(iluErrorString(error));
}

void DevILModule::initialize() {
    VoreenModule::initialize();

    // initialize DevIL
    ilInit();
    iluInit();
    ilEnable(IL_ORIGIN_SET); //< flip images
}

std::vector<std::string> DevILModule::getSupportedReadExtensions() {
    return readExtensions_;
}

std::vector<std::string> DevILModule::getSupportedWriteExtensions() {
    return writeExtensions_;
}

} // namespace
