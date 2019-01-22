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

#ifndef VRN_FILTERPROPERTIES_H
#define VRN_FILTERPROPERTIES_H

#include "voreen/core/io/serialization/serializable.h"
#include "voreen/core/properties/property.h"
#include "../volumefiltering/volumefilter.h"

#include <vector>
#include <map>

namespace voreen {

class VolumeBase;
class VolumeFilter;

class FilterProperties : public Serializable {
public:

    static const int DEFAULT_SETTINGS;

    virtual ~FilterProperties();

    const std::vector<Property*> getProperties() const;

    void storeVisibility();
    void restoreVisibility();

    virtual std::string getVolumeFilterName() const = 0;
    virtual void adjustPropertiesToInput(const VolumeBase& input) = 0;
    virtual VolumeFilter* getVolumeFilter(const VolumeBase& volume, int instanceId) const = 0;
    virtual void storeInstance(int instanceId) = 0;
    virtual void restoreInstance(int instanceId) = 0;
    virtual void removeInstance(int instanceId) = 0;

protected:

    virtual void addProperties() = 0;

    std::string getId(const std::string& id) const;

    std::vector<Property*> properties_;
    std::map<Property*, bool> visibilityMap_;

    static const std::string loggerCat_;
};

}

#endif