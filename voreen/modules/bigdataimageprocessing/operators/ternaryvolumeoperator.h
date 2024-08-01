/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_TERNARYVOLUMEOPERATOR_H
#define VRN_TERNARYVOLUMEOPERATOR_H

#include "voreen/core/datastructures/volume/volumeoperator.h"

namespace voreen {

//Maybe merge this into volumeoperator.h?

// Base interface for all binary volume operators:
class TernaryVolumeOperatorBase {
public:
    virtual bool isCompatible(const VolumeBase* volume1, const VolumeBase* volume2, const VolumeBase* volume3) const = 0;
};

//-----------------------------------------------------------------

// factory-like (does not really produce objects)
template<typename BASE_TYPE>
class UniversalTernaryVolumeOperatorGeneric {
public:
    static const BASE_TYPE* get(const VolumeBase* vh1, const VolumeBase* vh2, const VolumeBase* vh3);

    static void addInstance(BASE_TYPE* inst);
private:
    static UniversalTernaryVolumeOperatorGeneric<BASE_TYPE>* instance_;
    static UniversalTernaryVolumeOperatorGeneric<BASE_TYPE>* getInstance() {
        if(!instance_)
            instance_ = new UniversalTernaryVolumeOperatorGeneric<BASE_TYPE>();

        return instance_;
    }

    std::vector<BASE_TYPE*> instances_;
};

template<typename BASE_TYPE>
UniversalTernaryVolumeOperatorGeneric<BASE_TYPE>* UniversalTernaryVolumeOperatorGeneric<BASE_TYPE>::instance_ = 0;

template<typename BASE_TYPE>
const BASE_TYPE* UniversalTernaryVolumeOperatorGeneric<BASE_TYPE>::get(const VolumeBase* vh1, const VolumeBase* vh2, const VolumeBase* vh3) {
    for(size_t i=0; i<getInstance()->instances_.size(); i++) {
        if(getInstance()->instances_[i]->isCompatible(vh1, vh2, vh3))
            return getInstance()->instances_[i];
    }
    throw VolumeOperatorUnsupportedTypeException();
    return 0;
}

template<typename BASE_TYPE>
void UniversalTernaryVolumeOperatorGeneric<BASE_TYPE>::addInstance(BASE_TYPE* inst) {
    getInstance()->instances_.push_back(inst);
}
} // namespace


#endif
