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

#ifndef VRN_VOLUMEDERIVEDDATA_H
#define VRN_VOLUMEDERIVEDDATA_H

#include "voreen/core/voreenobject.h"
#include "tgt/assert.h"

#include <boost/thread.hpp>

namespace voreen {

class VolumeBase;

/**
 * Marker interface for data that is derived from a volume dataset,
 * such as histograms or preview images.
 *
 * @note each concrete subclass must provide a default constructor
 * @note each concrete subclass must be handled by a SerializableFactory
 *
 * @see DerivedDataFactory
 * @see VoreenModule::registerSerializerFactory
 */
class VRN_CORE_API VolumeDerivedData : public VoreenSerializableObject {
public:
    VolumeDerivedData();
    virtual ~VolumeDerivedData() {}

    /**
     * Virtual constructor.
     * Check for boo
     *
     * @param handle the volume to derive from
     * @return the constructed derived data object, or null if construction is not possible
     */
    virtual VolumeDerivedData* createFrom(const VolumeBase* handle) const = 0;

    virtual void serialize(Serializer& s) const = 0;
    virtual void deserialize(Deserializer& s) = 0;
};

} // namespace voreen

#endif
