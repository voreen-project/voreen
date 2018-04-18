/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
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

#ifndef VRN_COLOR_POOL_H
#define VRN_COLOR_POOL_H

#include "voreen/core/voreencoreapi.h"

#include "tgt/vector.h"

#include <vector>

namespace voreen {

class EnsembleDataset;

class VRN_CORE_API ColorPool {
public:

    static const tgt::vec3 FADE_OUT_COLOR;

    ColorPool(const EnsembleDataset& dataset);
    ~ColorPool();

    /**
     * Returns a unique color for the specified run and time step.
     * Run determines the general color and timeStep the saturation.
     */
    tgt::vec3 getColor(const std::string& run, size_t timeStep = -1);
    tgt::vec3 getColor(size_t run, size_t timeStep = -1);

    /**
     * Grants access to a bunch of colors being quite different from each other.
     * The same index always maps to the same color.
     * If there is no color mapped to the index, the color white will be returned.
     */
    static tgt::vec3 getDistinctColor(size_t index);

private:

    const EnsembleDataset& dataset_;
};

}

#endif
