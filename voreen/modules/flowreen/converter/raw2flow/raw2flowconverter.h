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

#ifndef RAW2FLOWCONVERTER_H
#define RAW2FLOWCONVERTER_H

#include <string>
#include <vector>

class Raw2FlowConverter {
public:
    enum VoxelOrientation { XYZ = 1, XZY = 0, YXZ = 2, YZX = 3, ZXY = 4, ZYX = 5};

public:
    Raw2FlowConverter(const unsigned int dimensions = 3);
    ~Raw2FlowConverter() {}
    bool convert(const std::string& filename) const;
    bool setup(const std::string& sizes, const std::string& voxelOrientation = "xyz", const std::string& reverseSlices = "");

private:
    enum { MAGIC_NUMBER_SIZE = 11 };

private:
    std::vector<unsigned int> parseSizes(const std::string& input) const;
    char parseSliceReversal(const std::string& input) const;
    VoxelOrientation parseVoxelOrientation(const std::string& input) const;
    void writeHeader(std::ofstream& ofs) const;

private:
    VoxelOrientation voxelOrientation_;
    static char magicNumber_[MAGIC_NUMBER_SIZE];
    static int fileVersion_;
    unsigned int dimensions_;
    union {
        struct { unsigned int x_, y_, z_; };
        unsigned int sizes_[3];
    };
    unsigned int dataSize_;
    char sliceReversal_;
    bool setup_;
};

#endif
