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

#ifndef VRN_COMMANDS_MODIFY_H
#define VRN_COMMANDS_MODIFY_H

//#include "voreen/core/utils/cmdparser/command.h"

#include <vector>
#include <string>

namespace voreen {

struct BrickingInformation;
class Volume;

class CommandCutToPieces /*: public Command*/ {
public:
    CommandCutToPieces();
    bool checkParameters(const std::vector<std::string>& parameters);
    bool execute(const std::vector<std::string>& parameters);

    static std::string loggerCat_;
};

class CommandScale /*: public Command */ {
public:
    CommandScale();
    bool checkParameters(const std::vector<std::string>& parameters);
    bool execute(const std::vector<std::string>& parameters);

    static std::string loggerCat_;
};

class CommandMirrorZ /*: public Command*/ {
public:
    CommandMirrorZ();
    bool execute(const std::vector<std::string>& parameters);

    static std::string loggerCat_;
};

class CommandSubSet /*: public Command*/ {
public:
    CommandSubSet();
    bool checkParameters(const std::vector<std::string>& parameters);
    bool execute(const std::vector<std::string>& parameters);

    static std::string loggerCat_;
};

class CommandBrick /*: public Command*/ {
public:
    CommandBrick();
    bool execute(const std::vector<std::string>& parameters);
    bool checkParameters(const std::vector<std::string>& parameters);
    /**
    * Gets information like spacing etc from a volume. This function is used to fill
    * the brickingInformation struct with the necessary information to brick the volume.
    */
    void getVolumeInformation(BrickingInformation& brickingInformation, const VolumeRAM* volume);

    static std::string loggerCat_;
};


class CommandScaleTexCoords /*: public Command*/ {
public:
    CommandScaleTexCoords();
    bool execute(const std::vector<std::string>& parameters);

    static std::string loggerCat_;
};


}   //namespace voreen

#endif //VRN_COMMANDS_CONVERT_H
