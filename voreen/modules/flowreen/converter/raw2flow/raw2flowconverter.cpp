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

#include "raw2flowconverter.h"
#include "stringutility.h"
#include <iostream>
#include <fstream>
#include <sstream>

char Raw2FlowConverter::magicNumber_[MAGIC_NUMBER_SIZE] = "VOREENFLOW";
int Raw2FlowConverter::fileVersion_ = 2;

Raw2FlowConverter::Raw2FlowConverter(const unsigned int dimensions)
    : voxelOrientation_(XYZ),
    dimensions_(dimensions),
    x_(0), y_(0), z_(0),
    dataSize_(0),
    sliceReversal_('\0'),
    setup_(false)
{
}

bool Raw2FlowConverter::convert(const std::string& filename) const {
    if (filename.empty() == true) {
        std::cout << "convert(): Error: given file name is empty!\n";
        return false;
    }

    if (setup_ == false) {
        std::cout << "convert(): Error: setup must be called first before converting!\n";
        return false;
    }

    std::string outFile(filename + ".flow");
    std::string fileExtension = StringUtility::getFileExtension(filename);
    if (fileExtension.empty() == false) {
        size_t nameLength = filename.size() - (fileExtension.size() + 1);
        outFile = filename.substr(0, nameLength) + ".flow";
    }

    std::ifstream ifs(filename.c_str(), std::ios_base::in | std::ios_base::binary);
    if (ifs.fail() == true) {
        std::cout << "convert(): Error: cannot open file '" << filename << "' for reading!\n";
        return false;
    }

    std::ofstream ofs(outFile.c_str(), std::ios_base::out | std::ios_base::binary);
    if (ofs.fail() == true) {
        std::cout << "convert(): Error: cannot open file '" << outFile << "' for writing!\n";
        return false;
    }

    std::cout << "writing converted file to '" << outFile << "'...";
    writeHeader(ofs);

    const std::streamsize BUFFER_SIZE = 4092; // 341 elements consisting of 3 floats @ 4 bytes take 4092 bytes
    char* buffer = new char[BUFFER_SIZE];
    while ((ifs.eof() == false) && (ifs.fail() == false)) {
        ifs.read(buffer, BUFFER_SIZE);
        ofs.write(buffer, ifs.gcount());
    }
    delete [] buffer;

    ifs.close();
    ofs.close();
    std::cout << "\tdone.\n";
    return true;
}

bool Raw2FlowConverter::setup(const std::string& sizes, const std::string& voxelOrientation,
                              const std::string& reverseSlices)
{
    if (dimensions_ > 3) {
        std::cout << "setup(): Error: dimension must be less or equal 3! Cannot proceed.\n";
        return false;
    }

    if (setup_ == true) {
        std::cout << "setup(): Error: setup() has already been called! Cannot proceed.\n";
        return false;
    }

    std::vector<unsigned int> sizesVec = parseSizes(sizes);
    if (sizesVec.size() != dimensions_)
        return false;

    for (size_t i = 0; i < dimensions_; ++i)
        sizes_[i] = (sizesVec[i] > 0) ? sizesVec[i] : 1;

    dataSize_ = x_ * y_ * z_ * sizeof(float) * dimensions_;
    voxelOrientation_ = parseVoxelOrientation(voxelOrientation);
    sliceReversal_ = parseSliceReversal(reverseSlices);
    setup_ = true;
    return (dataSize_ > 0);
}

void Raw2FlowConverter::writeHeader(std::ofstream& ofs) const {
    if (setup_ == false) {
        std::cout << "writeHeader(): Error: setup() must be called first before writing anything!\n";
        return;
    }

    ofs.write(magicNumber_, sizeof(char) * MAGIC_NUMBER_SIZE);
    ofs.write(reinterpret_cast<const char*>(&fileVersion_), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&dimensions_), sizeof(unsigned int));
    ofs.write(reinterpret_cast<const char*>(&voxelOrientation_), sizeof(char));
    ofs.write(reinterpret_cast<const char*>(&sliceReversal_), sizeof(char));
    ofs.write(reinterpret_cast<const char*>(&x_), sizeof(unsigned int));
    ofs.write(reinterpret_cast<const char*>(&y_), sizeof(unsigned int));
    ofs.write(reinterpret_cast<const char*>(&z_), sizeof(unsigned int));
    ofs.write(reinterpret_cast<const char*>(&dataSize_), sizeof(unsigned int));
}

std::vector<unsigned int> Raw2FlowConverter::parseSizes(const std::string& input) const {
    std::vector<unsigned int> result;
    std::vector<std::string> parts = StringUtility::explode(input, 'x');

    if (parts.size() != dimensions_)
        return result;

    for (size_t i = 0; i < parts.size(); ++i)
    {
        std::istringstream iss(parts[i]);
        unsigned int aux = 0;
        iss >> aux;
        result.push_back(aux);
    }
    return result;
}

char Raw2FlowConverter::parseSliceReversal(const std::string& input) const {
    if (input.size() != 1)
        return '\0';

    switch (input[0]) {
        case 'x':
        case 'X':
            return 'x';
        case 'y':
        case 'Y':
            return 'y';
        case 'z':
        case 'Z':
            return 'z';
        default:
            break;
    }
    return '\0';
}

Raw2FlowConverter::VoxelOrientation Raw2FlowConverter::parseVoxelOrientation(const std::string& input) const {
    VoxelOrientation orientation = XYZ;

    if (input.size() != dimensions_)
        return orientation;

    // Encode the given character string as number which can be interpreted as a
    // number at base 3.
    //
    int check = 0;
    for (size_t i = 0, digit = 1; i < dimensions_; ++i, digit *= 10) {
        switch (input[i]) {
            case 'x':
            case 'X':
                check += (0 * digit);
                break;
            case 'y':
            case 'Y':
                check += (1 * digit);
                break;
            case 'z':
            case 'Z':
                check += (2 * digit);
                break;
        }
    }

    switch (check) {
        case 210:
            return XYZ;
        case 120:
            return XZY;
        case 201:
            return YXZ;
        case 21:
            return YZX;
        case 102:
            return ZXY;
        case 12:
            return ZYX;
        default:
            break;
    }

    return orientation;
}
