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

#ifndef VRN_BRICKEDVOLUME_H
#define VRN_BRICKEDVOLUME_H

#include "brickedvolumebase.h"
#include "memory.h"


#include <tgt/vector.h>


namespace voreen{
/**
 * Template implementation of bricked volume for different 
 * base data types
 */
template<typename T>
class BrickedVolume :public BrickedVolumeBase{
public:
    virtual tgt::vec3     getDim() { return dim_;};
    /**
     * Get the bytes of the bricked volume
     * \return The memory of the bricked volumeRangeMinMax
     *         aligned to 16 bytes to enable use of 
     *         instructions for aligned memory
     */
    virtual const void*   getData() {return data_;};
    virtual const size_t* getOffsetX() { return offsetx_; }
    virtual const size_t* getOffsetY() { return offsety_; }
    virtual const size_t* getOffsetZ() { return offsetz_; }

    /**
     * Constructor of the bricked volumeRangeMinMax
     * \param volume_data The linear volume
     * \param dim Dimensions of the volumeRangeMinMax
     * \param bricksize The desired brick size
     */
    BrickedVolume(T* volume_data, tgt::svec3 dim, int bricksize);
    ~BrickedVolume();

private:
    tgt::svec3 dim_;
    T *        data_;
    size_t *   offsetx_;
    size_t *   offsety_;
    size_t *   offsetz_;

    BrickedVolume (const BrickedVolume &);
    BrickedVolume & operator = (const BrickedVolume &);
};
}

template<typename T>
voreen::BrickedVolume<T>::BrickedVolume(T* volume_data, tgt::svec3 dim, int bricksize){

    // round up sizes of dimensions to next bricksize
    size_t sx = (dim.x+bricksize-1)/bricksize*bricksize;
    size_t sy = (dim.y+bricksize-1)/bricksize*bricksize;
    size_t sz = (dim.z+bricksize-1)/bricksize*bricksize;

    size_t size = sx*sy*sz;

    dim_     = dim;
    data_    = aligned_malloc<T>(size);
    offsetx_ = aligned_malloc<size_t>(dim.x);
    offsety_ = aligned_malloc<size_t>(dim.y);
    offsetz_ = aligned_malloc<size_t>(dim.z);

    size_t bricksize2 = bricksize*bricksize;
    size_t bricksize3 = bricksize*bricksize*bricksize;

    // build offset tables

    // brick id in 1d
#define brick_id(x) ((x)/bricksize)

    // offset in brick
#define brick_off(x) ((x)%bricksize)

    for(size_t i = 0; i != dim.x; i++){
        offsetx_[i]  = brick_id(i)*bricksize3       + brick_off(i);
    }
    for(size_t i = 0; i != dim.y; i++){
        offsety_[i] = brick_id(i)*bricksize3*(sx/bricksize)    + brick_off(i)*bricksize;
    }
    for(size_t i = 0; i != dim.z; i++){
        offsetz_[i] = brick_id(i)*bricksize3*(sx/bricksize)*(sy/bricksize) + brick_off(i)*bricksize2;
    }
#undef brick_id
#undef brick_off

    // copy volume data to new schema
    for(size_t z = 0; z != dim.z; z++){
        for(size_t y = 0; y != dim.y; y++){
            for(size_t x = 0; x != dim.x; x++){
                // x dimension most inner, because it has
                // most cache locality for reading
                size_t idx = offsetx_[x]+offsety_[y]+offsetz_[z];
                size_t old_idx = x+dim.x*y+dim.x*dim.y*z;
                assert(idx < size);
                data_[idx] = volume_data[old_idx];
            }
        }
    }

    // check if the new volume has the correct values
    // and the addressing works
#if 0
    for(int z = 0; z != dim.z; z++){
        for(int y = 0; y != dim.y; y++){
            for(int x = 0; x != dim.x; x++){
                // x dimension most inner, because it has
                // most cache locality for reading
                size_t idx = offsetx_[x]+offsety_[y]+offsetz_[z];
                size_t old_idx = x+dim.x*y+dim.x*dim.y*z;
                assert(idx < size);
                assert(data_[idx] == volume_data[old_idx]);
            }
        }
    }
#endif

}

template<typename T>
voreen::BrickedVolume<T>::~BrickedVolume(){
    aligned_free(data_);
    aligned_free(offsetx_);
    aligned_free(offsety_);
    aligned_free(offsetz_);
}
#endif 
