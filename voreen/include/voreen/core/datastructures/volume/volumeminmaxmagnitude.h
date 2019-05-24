/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_VOLUMEMINMAXMAGNITUDE_H
#define VRN_VOLUMEMINMAXMAGNITUDE_H

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumederiveddata.h"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

namespace voreen {

class VRN_CORE_API VolumeMinMaxMagnitude : public VolumeDerivedData {
public:
    /// Empty default constructor required by VolumeDerivedData interface.
    VolumeMinMaxMagnitude();
    VolumeMinMaxMagnitude(const VolumeMinMaxMagnitude& other);
    VolumeMinMaxMagnitude(float minMag, float maxMag);
    virtual std::string getClassName() const { return "VolumeMinMaxMagnitude"; }

    virtual VolumeDerivedData* create() const;

    virtual VolumeDerivedData* createFrom(const VolumeBase* handle) const;

    /// @see VolumeDerivedData
    virtual void serialize(Serializer& s) const;

    /// @see VolumeDerivedData
    virtual void deserialize(Deserializer& s);

    /// Minimum (RealWorld)
    float getMinMagnitude() const;

    /// Maximum (RealWorld)
    float getMaxMagnitude() const;
protected:
    float minMagnitude_;
    float maxMagnitude_;
};

} // namespace voreen

#endif
