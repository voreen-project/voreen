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

#include <iostream>

#include "raw2flowconverter.h"
#include "stringutility.h"

void printUsage() {
    std::cout << "\nraw2flow:\n---------\n";
    std::cout << "Converts .raw files containing 3D flow data in single\n";
    std::cout << "precision float format to .flow files readable by Voreen module\n";
    std::cout << "Flowreen. The converted file simply has an additional header added.\n\n";
    std::cout << "Usage: raw2flow --dimensions=XxYxZ\n";
    std::cout << "\t[--voxel-orientation={xyz | xzy | yxz | yzx | zxy | zyx}]\n";
    std::cout << "\t[--slice-reversal={x|y|z}]\n";
    std::cout << "\tfilename1 filename2 ... filenameN\n";
}

int main(int argc, char** argv) {
    if (argc < 3) {
        printUsage();
        return -1;
    }

    std::string argDimensions, argOrientation, argReversal;
    int offsetFilenames = argc - 1;
    for (int i = 1; i < (argc - 1); ++i) {
        std::string arg(argv[i]);
        if (arg.substr(0, 2) != "--") {
            offsetFilenames = (i < offsetFilenames) ? i : offsetFilenames;
            continue;
        }

        std::vector<std::string> args = StringUtility::explode(arg, '=');
        if (args.size() != 2) {
            std::cout << "ignoring invalid argument '" << argv[i] << "'...\n";
            continue;
        }

        if (args[0] == "--dimensions")
            argDimensions = args[1];
        else if (args[0] == "--voxel-orientation")
            argOrientation = args[1];
        else if (args[0] == "--slice-reversal")
            argReversal = args[1];
        else
            std::cout << "ignoring unknown option '" << args[0] << "'...\n";
    }

    if (argDimensions.empty() == true) {
        std::cout << "Error: required argument --dimensions is missing or invalid!\n";
        printUsage();
        return -2;
    }

    Raw2FlowConverter converter(3);
    if (converter.setup(argDimensions, argOrientation, argReversal) == true) {
        for (int i = offsetFilenames; i < argc; ++i)
            converter.convert(argv[i]);
    } else {
        std::cout << "Error: failed to setup converter - probably due to wrong dimensions settings.\n";
        return -3;
    }

    return 0;
}
