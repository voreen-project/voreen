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

#ifndef VRN_COMMANDS_MOTION_H
#define VRN_COMMANDS_MOTION_H

//#include "voreen/core/utils/cmdparser/command.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

class CommandCluster /*: public Command*/ {
public:
    CommandCluster();
    bool checkParameters(const std::vector<std::string>& parameters);
    bool execute(const std::vector<std::string>& parameters);
};

class CommandCreateMotion /*: public Command*/ {
public:
    CommandCreateMotion();
    bool checkParameters(const std::vector<std::string>& parameters);
    bool execute(const std::vector<std::string>& parameters);

    /**
    * Creates a dataset with starlike motion. The motion points from the center to the borders.
    */
    VolumeRAM_4xUInt8* createMotionDS(int size);
    /**
    * Creates a dataset with starlike motion vectors. It divides the box in 8 subvolumes. In
    * each subvolume the motion points from the center to the borders.
    */
    VolumeRAM_4xUInt8* createMotionDS2(int size);
    /**
    * Creates a dataset with two spheres of motion. The motion in each sphere points in different directions.
    */
    VolumeRAM_4xUInt8* createMotionDS3(int size);
    VolumeRAM_4xUInt8* createMotionDS4(int size);

    /**
    * Creates motion vectors for each voxel in the box. The motion points starlike from the center
    * border.
    */
    void createMotionBox(VolumeRAM_4xUInt8* dataset, tgt::ivec3 center, int size, int minValue=0, int maxValue=127);
    /**
    * Creates a sphere with motion vectors within the box.
    */
    void createMotionBox2(VolumeRAM_4xUInt8* dataset, tgt::ivec3 center, int size, int radius, tgt::ivec3 direction, int minValue=0, int maxValue=127);
    /**
    * Creates a sphere with motion vectors within the box. The motion turns from slice to slice
    */
    void createMotionBox3(VolumeRAM_4xUInt8* dataset, tgt::ivec3 center, int size, int radius, int minValue=0, int maxValue=127);
};

}   //namespace voreen

#endif //VRN_COMMANDS_MOTION_H
