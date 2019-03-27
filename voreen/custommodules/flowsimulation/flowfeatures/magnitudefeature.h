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

#ifndef VRN_MAGNITUREFEATURE_H
#define VRN_MAGNITUREFEATURE_H

#include "flowfeature.h"

namespace voreen {

class MagnitudeFeature : public FlowFeature {
public:

    MagnitudeFeature(const std::string& sliceBaseType);
    virtual ~MagnitudeFeature();

    virtual std::string getName() const { return "MagnitudeFeature"; }
    virtual VolumeFilter* create() const { return new MagnitudeFeature(sliceBaseType_); }

    int zExtent() const;
    const std::string& getSliceBaseType() const;

    // For now we only support single channel volumes
    size_t getNumInputChannels() const { return 3; };
    size_t getNumOutputChannels() const { return 1; };

    std::unique_ptr<VolumeRAM> getFilteredSlice(const CachingSliceReader* src, int z) const;

private:

    const std::string sliceBaseType_;
};

} // namespace voreen

#endif // VRN_THRESHOLDINGFILTER_H
