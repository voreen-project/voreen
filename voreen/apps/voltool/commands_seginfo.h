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

#ifndef VRN_COMMANDS_SEGINFO_H
#define VRN_COMMANDS_SEGINFO_H

//#include "voreen/core/utils/cmdparser/command.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

class CommandSegInfo /*: public Command*/ {
public:
    CommandSegInfo();
    bool checkParameters(const std::vector<std::string>& parameters);
    bool execute(const std::vector<std::string>& parameters);

    static std::string loggerCat_;
};

class CommandCutToSegs /*: public Command*/ {
public:
    CommandCutToSegs();
    bool checkParameters(const std::vector<std::string>& parameters);
    bool execute(const std::vector<std::string>& parameters);

    static std::string loggerCat_;
};

class CommandCubify /*: public Command*/ {
public:
    CommandCubify();
    bool checkParameters(const std::vector<std::string>& parameters);
    bool execute(const std::vector<std::string>& parameters);

    static std::string loggerCat_;
};

}   //namespace voreen

#endif //VRN_COMMANDS_SEGINFO_H
