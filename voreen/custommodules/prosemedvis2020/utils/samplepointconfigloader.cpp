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

#include "samplepointconfigloader.h"

#include "voreen/core/io/serialization/serializer.h"
#include "voreen/core/io/serialization/deserializer.h"
#include "voreen/core/utils/stringutils.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

namespace voreen {


std::vector<SamplePointConfig> SamplePointConfigLoader::loadSamplepointConfigFile(const std::string& filepath){
    if (filepath == "") {
        return std::vector<SamplePointConfig>();
    }

    std::vector<SamplePointConfig> result;

    std::ifstream infile(filepath);
    if (!infile.is_open()) {
        // File could not be opened --> return empty vector.
        //LWARNING("Failed to open inputfile " << filepath);
        return result;
    }

    // Read inputfile
    std::string line;
    while (std::getline(infile, line)) {
        // Split current line 
        std::stringstream ls(line);
        
        // Read name
        std::string name;
        std::getline(ls, name, ',');
        if (ls.bad()) { continue; }
        name = trim(name);

        // Read xPos
        std::string xPos;
        std::getline(ls, xPos, ',');
        if (ls.bad()) { continue; }
        xPos = trim(xPos);

        // Read yPos
        std::string yPos;
        std::getline(ls, yPos, ',');
        if (ls.bad()) { continue; }
        yPos = trim(yPos);

        // Read zPos
        std::string zPos;
        std::getline(ls, zPos, ',');
        if (ls.bad()) { continue; }
        zPos = trim(zPos);

        // Read Texture-Filename:
        std::string texFilename;
        std::getline(ls, texFilename);
        if (ls.bad()) { continue; }
        texFilename = trim(texFilename);

        SamplePointConfig currentConfig;
        currentConfig.name_ = name;
        currentConfig.defaultPos_.x = std::stof(xPos);
        currentConfig.defaultPos_.y = std::stof(yPos);
        currentConfig.defaultPos_.z = std::stof(zPos);
        currentConfig.textureFilename_ = texFilename;

        result.push_back(currentConfig);
    }

    infile.close();

    return result;
}

void SamplePointConfig::serialize(Serializer& s) const {
    s.serialize("name_", name_);
    s.serialize("defaultPos_", defaultPos_);
    s.serialize("textureFilename_", textureFilename_);
}

void SamplePointConfig::deserialize(Deserializer& s) {
    s.deserialize("name_", name_);
    s.deserialize("defaultPos_", defaultPos_);
    s.deserialize("textureFilename_", textureFilename_);
}

}
