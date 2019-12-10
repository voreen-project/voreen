/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_BRICKEDVOLUMEBASE_H
#define VRN_BRICKEDVOLUMEBASE_H
#include <tgt/vector.h>
#include "memory.h"


namespace voreen{

/**
 * Baseclass for the bricked volume. 
 */
class BrickedVolumeBase{
public:
    /**
     * Get the dimension of the bricked volume
     */
    virtual tgt::vec3     getDim() = 0;

    /**
     * Get the bricked volume
     * \return The bricked volume bytes
     */
    virtual const void*   getData() = 0;

    /**
     * Offset table for x component
     */
    virtual const size_t* getOffsetX() = 0;
    
    /**
     * Offset table for y component
     */
    virtual const size_t* getOffsetY() = 0;
    
    /**
     * Offset table for y component
     */
    virtual const size_t* getOffsetZ() = 0;
    virtual ~BrickedVolumeBase(){};
};

}
#endif
