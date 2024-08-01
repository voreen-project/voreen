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

#ifndef VRN_CLUSTERSPLITTINGTHREAD
#define VRN_CLUSTERSPLITTINGTHREAD

#include "voreen/core/voreencoreapi.h"

#include "voreen/core/utils/backgroundthread.h"

#include "connectedcomponentqueue.h"
#include "../processors/nucleiclustersplitting.h"

#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

namespace voreen {

class VRN_CORE_API ClusterSplittingThread : public BackgroundThread {

public:

    ClusterSplittingThread(NucleiClusterSplitting* processor, ConnectedComponentQueue& queue, boost::mutex& globalValueMutex, bool smoothForSeedDetection, size_t gaussianKernelSize, float gaussianSigma, bool maskBeforeSmoothing, bool expandMarkers, size_t minSeedRank);

protected:

    // do not use the default constructor
    ClusterSplittingThread();

    /**
     * struct for specification of a single seed, i.e., local maximum, within the cluster
     */
    struct SeedDescriptor {
        tgt::svec3 CenterOfMaxima;
        float DistanceToClosestMaximum;
        uint16_t identifier;
    };

    virtual void threadMain();

    virtual void clearLocalData();

    virtual void pushLocalData();

    virtual void handleInterruption();

    NucleiClusterSplitting* processor_;         ///< processor which started the tread an receives results on termination
    ConnectedComponentQueue& queue_;            ///< queue where the thread will get the next connected component for processing

    PointSegmentListGeometryVec3 localResult_;  ///< result for the connected components processed by this thread
    boost::mutex& globalValueMutex_;            ///< mutex for accessing variables (i.e., progress or global results) in the processor

    bool smoothForSeedDetection_;
    size_t gaussianKernelSize_;
    float gaussianSigma_;
    bool maskBeforeSmoothing_;
    bool expandMarkers_;
    size_t minSeedRank_;
};

} // namespace

#endif // VRN_CLUSTERSPLITTING_THREAD
