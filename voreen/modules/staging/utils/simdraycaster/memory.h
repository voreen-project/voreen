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

#pragma once
#include <stdlib.h>
namespace voreen{
    /**
     * Allocate aligned memory. The memory is always big enough for any simd read,
     * even if it would span over the end of the allocation
     * \param bytes The size in bytes
     * \return The allocated memory. Id needs to be freed with the aligned_free* functions
     */
    void* aligned_malloc_bytes(size_t bytes);
    /**
     * Frees aligned memory
     */
    void  aligend_free_bytes(void*);

    /**
     * Typed version of aligned_free_bytes
     * \param size The number of elements, for which memory is allocated
     */
    template<typename T>
    T* aligned_malloc(size_t size){
        return reinterpret_cast<T*>(aligned_malloc_bytes(size*sizeof(T)));
    }

    /**
     * Typed version of aligned_free_bytes
     */
    template<typename T>
    void aligned_free(T* data){
        aligend_free_bytes(data);
    }



};
