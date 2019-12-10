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

#ifndef VRN_OCTREEQUANTIFICATIONTHREAD
#define VRN_OCTREEQUANTIFICATIONTHREAD

#include "voreen/core/voreencoreapi.h"

#include "voreen/core/utils/backgroundthread.h"

#include "octreequantificationnodequeue.h"
#include "../processors/octreesegmentationquantification.h"
#include "octreequantificationresults.h"

namespace voreen {

class VRN_CORE_API OctreeQuantificationThread : public BackgroundThread { 

public:

    OctreeQuantificationThread(OctreeSegmentationQuantification* quantificationProcessor, OctreeQuantificationNodeQueue& queue, 
            const VolumeOctree* octree, std::vector<const VolumeRAM*> rois, std::vector<tgt::vec3> coordinateRatios, 
            OctreeSegmentationQuantification::SegmentationInterpolationMode simMode, tgt::svec3 volumeURB, size_t numChannels,
            boost::mutex& globalValueMutex); 

protected:

    virtual void threadMain();

    virtual void clearLocalData();

    virtual void pushLocalData();

    virtual void handleInterruption();

    OctreeSegmentationQuantification* quantificationProcessor_;
    OctreeQuantificationNodeQueue& queue_;
    const VolumeOctree* volumeOctree_;
    std::vector<const VolumeRAM*> rois_;
    std::vector<tgt::vec3> coordinateRatios_;
    OctreeSegmentationQuantification::SegmentationInterpolationMode simMode_;
    size_t numChannels_;
    tgt::svec3 volumeURB_;
    boost::mutex& globalValueMutex_;

    OctreeQuantificationResults localResults_;

private:

    // do not use the default constructor
    OctreeQuantificationThread();
};

} // namespace

#endif // VRN_OCTREEQUANTIFICATIONTHREAD
