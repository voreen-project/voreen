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

#ifndef VRN_VOLUMEOPERATOR_H
#define VRN_VOLUMEOPERATOR_H

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/io/progressreporter.h"
#include "voreen/core/utils/exception.h"
#include "tgt/vector.h"

#ifdef VRN_MODULE_FLOWREEN
#include "modules/flowreen/datastructures/deprecated/volumeflow3d.h"
#endif

#include <vector>
#include <limits>

namespace voreen {

/**
 * Thrown by a VolumeOperator is the type of the passed Volume is not supported.
 */
class VRN_CORE_API VolumeOperatorUnsupportedTypeException : public VoreenException {
public:
    VolumeOperatorUnsupportedTypeException(const std::string& type = "unknown")
        : VoreenException("Volume type '" + type + "'not supported by this VolumeOperator")
    {}
};

/**
 * Thrown by a VolumeOperator is the type of the passed Volume is not supported.
 */
class VRN_CORE_API BinaryVolumeOperatorUnsupportedTypeException : public VoreenException {
public:
    BinaryVolumeOperatorUnsupportedTypeException(const std::string& type1 = "unknown", const std::string& type2 = "unknown")
        : VoreenException("Volume types '" + type1 + "' and '" + type2 + "' are either not supported by this VolumeOperator or not identical")
    {}
};

// Base interface for all unary volume operators:
class UnaryVolumeOperatorBase {
public:
    virtual bool isCompatible(const VolumeBase* volume) const = 0;
};

// macro to implement isCompatible in templates:
#define IS_COMPATIBLE \
bool isCompatible(const VolumeBase* volume) const { \
    const VolumeRAM* v = volume->getRepresentation<VolumeRAM>(); \
    if(!v) \
        return false; \
    const VolumeAtomic<T>* va = dynamic_cast<const VolumeAtomic<T>*>(v); \
    if(!va) \
        return false; \
    return true; \
}

//Unary: -----------------------------------------------------------------

// factory-like (does not really produce objects)
template<typename BASE_TYPE>
class UniversalUnaryVolumeOperatorGeneric {
public:
    static const BASE_TYPE* get(const VolumeBase* vh);

    static void addInstance(BASE_TYPE* inst);
private:
    static UniversalUnaryVolumeOperatorGeneric<BASE_TYPE>* instance_;
    static UniversalUnaryVolumeOperatorGeneric<BASE_TYPE>* getInstance() {
        if(!instance_)
            instance_ = new UniversalUnaryVolumeOperatorGeneric<BASE_TYPE>();

        return instance_;
    }

    std::vector<BASE_TYPE*> instances_;
};

template<typename BASE_TYPE>
UniversalUnaryVolumeOperatorGeneric<BASE_TYPE>* UniversalUnaryVolumeOperatorGeneric<BASE_TYPE>::instance_ = 0;

template<typename BASE_TYPE>
const BASE_TYPE* UniversalUnaryVolumeOperatorGeneric<BASE_TYPE>::get(const VolumeBase* vh) {
    for(size_t i=0; i<getInstance()->instances_.size(); i++) {
        if(getInstance()->instances_[i]->isCompatible(vh))
            return getInstance()->instances_[i];
    }
    throw VolumeOperatorUnsupportedTypeException(vh->getFormat());
    return 0;
}

template<typename BASE_TYPE>
void UniversalUnaryVolumeOperatorGeneric<BASE_TYPE>::addInstance(BASE_TYPE* inst) {
    getInstance()->instances_.push_back(inst);
}

// Binary:--------------------------------------------------------------------------------------

// Base interface for all binary volume operators:
class BinaryVolumeOperatorBase {
public:
    virtual bool isCompatible(const VolumeBase* volume1, const VolumeBase* volume2) const = 0;
};

//-----------------------------------------------------------------

// factory-like (does not really produce objects)
template<typename BASE_TYPE>
class UniversalBinaryVolumeOperatorGeneric {
public:
    static const BASE_TYPE* get(const VolumeBase* vh1, const VolumeBase* vh2);

    static void addInstance(BASE_TYPE* inst);
private:
    static UniversalBinaryVolumeOperatorGeneric<BASE_TYPE>* instance_;
    static UniversalBinaryVolumeOperatorGeneric<BASE_TYPE>* getInstance() {
        if(!instance_)
            instance_ = new UniversalBinaryVolumeOperatorGeneric<BASE_TYPE>();

        return instance_;
    }

    std::vector<BASE_TYPE*> instances_;
};

template<typename BASE_TYPE>
UniversalBinaryVolumeOperatorGeneric<BASE_TYPE>* UniversalBinaryVolumeOperatorGeneric<BASE_TYPE>::instance_ = 0;

template<typename BASE_TYPE>
const BASE_TYPE* UniversalBinaryVolumeOperatorGeneric<BASE_TYPE>::get(const VolumeBase* vh1, const VolumeBase* vh2) {
    for(size_t i=0; i<getInstance()->instances_.size(); i++) {
        if(getInstance()->instances_[i]->isCompatible(vh1, vh2))
            return getInstance()->instances_[i];
    }
    throw BinaryVolumeOperatorUnsupportedTypeException(vh1->getFormat(), vh2->getFormat());
    return 0;
}

template<typename BASE_TYPE>
void UniversalBinaryVolumeOperatorGeneric<BASE_TYPE>::addInstance(BASE_TYPE* inst) {
    getInstance()->instances_.push_back(inst);
}

// some macros for easy instanciation (see below):

#define INST_SCALAR_TYPES(univ_type, type) \
    univ_type::addInstance(new type<uint8_t>()); \
    univ_type::addInstance(new type<int8_t>()); \
    univ_type::addInstance(new type<uint16_t>()); \
    univ_type::addInstance(new type<int16_t>()); \
    univ_type::addInstance(new type<uint32_t>()); \
    univ_type::addInstance(new type<int32_t>()); \
    univ_type::addInstance(new type<uint64_t>()); \
    univ_type::addInstance(new type<int64_t>()); \
    univ_type::addInstance(new type<float>()); \
    univ_type::addInstance(new type<double>());

#define INST_VECTOR_TYPES(univ_type, type) \
    univ_type::addInstance(new type<tgt::Vector2<uint8_t> >()); \
    univ_type::addInstance(new type<tgt::Vector2<int8_t> >()); \
    univ_type::addInstance(new type<tgt::Vector2<uint16_t> >()); \
    univ_type::addInstance(new type<tgt::Vector2<int16_t> >()); \
    univ_type::addInstance(new type<tgt::Vector2<uint32_t> >()); \
    univ_type::addInstance(new type<tgt::Vector2<int32_t> >()); \
    univ_type::addInstance(new type<tgt::Vector2<uint64_t> >()); \
    univ_type::addInstance(new type<tgt::Vector2<int64_t> >()); \
    univ_type::addInstance(new type<tgt::Vector2<float> >()); \
    univ_type::addInstance(new type<tgt::Vector2<double> >()); \
    \
    univ_type::addInstance(new type<tgt::Vector3<uint8_t> >()); \
    univ_type::addInstance(new type<tgt::Vector3<int8_t> >()); \
    univ_type::addInstance(new type<tgt::Vector3<uint16_t> >()); \
    univ_type::addInstance(new type<tgt::Vector3<int16_t> >()); \
    univ_type::addInstance(new type<tgt::Vector3<uint32_t> >()); \
    univ_type::addInstance(new type<tgt::Vector3<int32_t> >()); \
    univ_type::addInstance(new type<tgt::Vector3<uint64_t> >()); \
    univ_type::addInstance(new type<tgt::Vector3<int64_t> >()); \
    univ_type::addInstance(new type<tgt::Vector3<float> >()); \
    univ_type::addInstance(new type<tgt::Vector3<double> >()); \
    \
    univ_type::addInstance(new type<tgt::Vector4<uint8_t> >()); \
    univ_type::addInstance(new type<tgt::Vector4<int8_t> >()); \
    univ_type::addInstance(new type<tgt::Vector4<uint16_t> >()); \
    univ_type::addInstance(new type<tgt::Vector4<int16_t> >()); \
    univ_type::addInstance(new type<tgt::Vector4<uint32_t> >()); \
    univ_type::addInstance(new type<tgt::Vector4<int32_t> >()); \
    univ_type::addInstance(new type<tgt::Vector4<uint64_t> >()); \
    univ_type::addInstance(new type<tgt::Vector4<int64_t> >()); \
    univ_type::addInstance(new type<tgt::Vector4<float> >()); \
    univ_type::addInstance(new type<tgt::Vector4<double> >());

#define INST_TENSOR_TYPES(univ_type, type) \
    univ_type::addInstance(new type<Tensor2<float> >());

#define APPLY_OP(vh, ...) get(vh)->apply(vh, ## __VA_ARGS__)
#define APPLY_B_OP(vh1, vh2, ...) get(vh1, vh2)->apply(vh1, vh2, ## __VA_ARGS__)

/**
 * Variation of VRN_FOR_EACH_VOXEL that updates a progress reporter
 * for each z-slice that has been processed. The PROGRESS parameter
 * may be null.
 */
#define VRN_FOR_EACH_VOXEL_WITH_PROGRESS(INDEX, POS, SIZE, PROGRESS) \
    for (tgt::svec3 INDEX = (POS); INDEX.z < (SIZE).z; ++INDEX.z, PROGRESS ? PROGRESS->setProgress(static_cast<float>(INDEX.z) / static_cast<float>(SIZE.z)) : void(0)) \
        for (INDEX.y = (POS).y; INDEX.y < (SIZE).y; ++INDEX.y)\
            for (INDEX.x = (POS).x; INDEX.x < (SIZE).x; ++INDEX.x)

/**
 * Variation of VRN_FOR_EACH_VOXEL that updates a progress reporter
 * for each z-slice that has been processed. In contrast to
 * VRN_FOR_EACH_VOXEL_WITH_PROGRESS, this macro is supposed
 * to be used for tasks that require to traverse the volume
 * multiple times (e.g. one for each axial direction). Therefore,
 * it additionally expects a scale factor as well as an offset
 * for the progress value. These two parameters have to be
 * passed as normalized values with regard to the range [0;1].
 * The PROGRESS parameter may be null.
 */
#define VRN_FOR_EACH_VOXEL_WITH_PROGRESS_SUB_TASK(INDEX, POS, SIZE, PROGRESS, PROGRESS_OFFSET, PROGRESS_SCALE) \
    for (tgt::svec3 INDEX = (POS); INDEX.z < (SIZE).z; ++INDEX.z, PROGRESS ? PROGRESS->setProgress(PROGRESS_OFFSET + ((static_cast<float>(INDEX.z) / static_cast<float>(SIZE.z))*PROGRESS_SCALE)) : void(0)) \
        for (INDEX.y = (POS).y; INDEX.y < (SIZE).y; ++INDEX.y)\
            for (INDEX.x = (POS).x; INDEX.x < (SIZE).x; ++INDEX.x)

} // namespace

#endif // VRN_VOLUMEOPERATOR_H
