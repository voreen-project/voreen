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

#include "voreen/core/datastructures/volume/modality.h"

namespace voreen {

const Modality Modality::MODALITY_UNKNOWN("unknown");
const Modality Modality::MODALITY_ANY("any");
const Modality Modality::MODALITY_CT("ct");
const Modality Modality::MODALITY_PET("pet");
const Modality Modality::MODALITY_MR("mr");
const Modality Modality::MODALITY_US("us");
const Modality Modality::MODALITY_SEGMENTATION("segmentation");
const Modality Modality::MODALITY_MASKING("masking");
const Modality Modality::MODALITY_LENSEVOLUME("lense volume");
const Modality Modality::MODALITY_AMBIENTOCCLUSION("ambient occlusion");
const Modality Modality::MODALITY_DYNAMICAMBIENTOCCLUSION("dynamic occlusion");
const Modality Modality::MODALITY_DIRECTION_X("direction x");
const Modality Modality::MODALITY_DIRECTION_Y("direction y");
const Modality Modality::MODALITY_DIRECTION_Z("direction z");
const Modality Modality::MODALITY_DIRECTIONS("directions");
const Modality Modality::MODALITY_NORMALS("normals");
const Modality Modality::MODALITY_GRADIENTS("gradients");
const Modality Modality::MODALITY_GRADIENT_MAGNITUDES("gradient_magnitudes");
const Modality Modality::MODALITY_2ND_DERIVATIVES("2nd_derivatives");
const Modality Modality::MODALITY_FLOW("flow");
const Modality Modality::MODALITY_INDEX_VOLUME("index volume");
const Modality Modality::MODALITY_BRICKED_VOLUME("bricked volume");
const Modality Modality::MODALITY_EEP_VOLUME("eep volume");

Modality::Modality(const std::string& name)
  : name_(name)
{
    if(name.empty()) {
        return;
    }
    for (size_t i = 0; i < modalityNames_().size(); ++i) {
        if (name == modalityNames_()[i])
            return;
    }
    modalityNames_().push_back(name_);
    // Create a copy of the current modality in order to avoid saving dangling pointers of non const objects:
    modalities().push_back(Modality(*this));
}

} // namespace voreen
