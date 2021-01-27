/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_VOLUMEPREVIEW_H
#define VRN_VOLUMEPREVIEW_H

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumederiveddata.h"

#include <string>
#include <iostream>
#include <fstream>

namespace voreen {

class VRN_CORE_API VolumePreview : public VolumeDerivedData {
public:
    /// Empty default constructor required by VolumeDerivedData interface.
    VolumePreview();
    VolumePreview(int height, const std::vector<unsigned char>& data);
    virtual std::string getClassName() const { return "VolumePreview"; }

    virtual VolumeDerivedData* create() const;
    virtual VolumeDerivedData* createFrom(const VolumeBase* handle) const;

    /// @see VolumeDerivedData
    virtual void serialize(Serializer& s) const;

    /// @see VolumeDerivedData
    virtual void deserialize(Deserializer& s);

    int getHeight() const {
        return height_;
    }

    const std::vector<unsigned char>& getData() const {
        return prevData_;
    }

protected:

    int height_;
    std::vector<unsigned char> prevData_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif
