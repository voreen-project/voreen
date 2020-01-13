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

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

namespace voreen {

VolumeProcessor::VolumeProcessor()
    : Processor()
{}

VolumeProcessor::~VolumeProcessor() {}

tgt::mat4 VolumeProcessor::computeConversionMatrix(const VolumeBase* originVolume, const VolumeBase* destinationVolume) const {
    if (originVolume == destinationVolume)
        return tgt::mat4::identity;
    else {
        tgt::mat4 voxelToWorldOrigin = originVolume->getVoxelToWorldMatrix();
        tgt::mat4 worldToVoxelDestination = destinationVolume->getWorldToVoxelMatrix();
        return worldToVoxelDestination * voxelToWorldOrigin;
    }
}

//------------------------------------------------------------------------

CachingVolumeProcessor::CachingVolumeProcessor()
    : VolumeProcessor()
    , useCaching_("useCaching", "Use Cache", true, VALID)
    , clearCache_("clearCache", "Clear Cache", VALID)
    , cache_(this)
{
    addProperty(useCaching_);

    clearCache_.onChange(MemberFunctionCallback<CachingVolumeProcessor>(this, &CachingVolumeProcessor::clearCache));
    addProperty(clearCache_);

    useCaching_.setGroupID("caching");
    clearCache_.setGroupID("caching");
}

CachingVolumeProcessor::~CachingVolumeProcessor() {}

void CachingVolumeProcessor::beforeProcess() {
    VolumeProcessor::beforeProcess();

    if (useCaching_.get() && VoreenApplication::app()->useCaching()) {
        if (cache_.restore()) {
            setValid();
        }
    }
}

void CachingVolumeProcessor::afterProcess() {
    VolumeProcessor::afterProcess();

    if (useCaching_.get() && VoreenApplication::app()->useCaching()) {
        cache_.store();
    }
}

void CachingVolumeProcessor::initialize() {
    VolumeProcessor::initialize();

    cache_.addAllInports();
    cache_.addAllOutports();
    cache_.addAllProperties();

    cache_.initialize();
}

void CachingVolumeProcessor::clearCache() {
    cache_.clearCache();
}

}   // namespace
