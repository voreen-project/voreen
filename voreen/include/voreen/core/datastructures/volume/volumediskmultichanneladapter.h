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

#ifndef VRN_VOLUMEDISKMULTICHANNELADAPTER_H
#define VRN_VOLUMEDISKMULTICHANNELADAPTER_H

#include "volumedisk.h"

namespace voreen {

class VolumeBase;

/**
 * Adapter class for a combining multiple volume (disk or RAM) representations (of same dimensions)
 * into a single, multi-channel volume.
 */
class VolumeDiskMultiChannelAdapter : public VolumeDisk {
public:

    /**
     * Constructor.
     * @param channels Each volume represents one channel of the new representation.
     * @param mirror optional mirror flags for each x, y, z dimension respectively.
     * @param swizzle optional swizzle arguments.
     * @param negate optional negation arguments.
     */
    VolumeDiskMultiChannelAdapter(const std::vector<const VolumeBase*>& channels,
                                  const tgt::bvec3& mirror = tgt::bvec3(false),
                                  const std::vector<size_t>& swizzle = {},
                                  const std::vector<bool>& negate = {});

    std::string getHash() const override;

    VolumeRAM* loadVolume() const override;

    VolumeRAM* loadSlices(const size_t firstZSlice, const size_t lastZSlice) const override;

    VolumeRAM* loadBrick(const tgt::svec3& offset, const tgt::svec3& dimensions) const override;

private:

    std::vector<const VolumeBase*> channels_;
    tgt::bvec3 mirror_;
    std::vector<size_t> swizzle_;
    std::vector<bool> negate_;
};

}

#endif