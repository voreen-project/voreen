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

#ifndef VRN_BINARYMEDIANFILTERPROPERTIES_H
#define VRN_BINARYMEDIANFILTERPROPERTIES_H

#include "filterproperties.h"

#include "../volumefiltering/binarymedianfilter.h"

namespace voreen {

class BinaryMedianFilterProperties : public FilterProperties {
public:
    BinaryMedianFilterProperties();

    virtual std::string getVolumeFilterName() const;

    virtual void adjustPropertiesToInput(const VolumeBase& input);

    virtual VolumeFilter* getVolumeFilter(const VolumeBase& volume, int instanceId) const;
    virtual void restoreInstance(int instanceId);
    virtual void storeInstance(int instanceId);
    virtual void removeInstance(int instanceId);
    virtual void addProperties();
    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

private:

    void updateObjectVoxelThreshold();

    struct Settings {
        bool useUniformExtent_;
        int extentX_;
        int extentY_;
        int extentZ_;
        float binarizationThreshold_;
        SamplingStrategyType samplingStrategyType_;
        float outsideVolumeValue_;
        bool forceMedian_;
        int objectVoxelThreshold_;
    };
    std::map<int, Settings> instanceSettings_;

    BoolProperty useUniformExtent_;
    IntProperty extentX_;
    IntProperty extentY_;
    IntProperty extentZ_;
    FloatProperty binarizationThreshold_;
    OptionProperty<SamplingStrategyType> samplingStrategyType_;
    FloatProperty outsideVolumeValue_;
    BoolProperty forceMedian_;
    IntProperty objectVoxelThreshold_;
};

}

#endif