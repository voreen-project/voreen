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

#ifndef VRN_CONNECTEDCOMPONENTQUEUE
#define VRN_CONNECTEDCOMPONENTQUEUE

#include "voreen/core/voreencoreapi.h"

#include "tgt/vector.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include <queue>

#include <boost/thread.hpp>
#include <boost/thread/condition_variable.hpp>

namespace voreen {

/**
 * Provides a threadsafe data structure for queueing connected components that have to be processed individually for segmenting individual cell nuclei. Use with ClusterSplittingThread class.
 */
class VRN_CORE_API ConnectedComponentQueue {

public:

    /**
     * struct containing the information for a single connected component in the input data
     */
    struct ConnectedComponent {
        size_t id_;                 ///< ID of the connected component
        size_t volume_;             ///< volume of the connected component (in voxels)
        tgt::svec3 llf_;            ///< offset of the input volume's LLF voxel
        tgt::svec3 urb_;

        VolumeRAM_UInt16* image_;    ///< brick of the original image containing this component
        VolumeRAM_UInt32* labels_;   ///< brick of the connected component labeling image containing this component

        tgt::vec3 spacing_;
        tgt::vec3 offset_;
    };

    /**
     * Creates a queue with a fized maximum size. 
     *
     * @param maxSize determines the maximum number of connected components that can simultaneously be queued. If set to 0, the queue does not have a maximum size. 
     */
    ConnectedComponentQueue(size_t maxSize = 0);
                    
    /**
     * Pushes a connected component onto the queue. If the queue is full the pushing thread has to wait until the size of the queue is lower than maxSize.
     */
    void push(ConnectedComponent component);

    /**
     * Returns if the queue is empty.
     */
    bool empty() const;

    /**
     * Returns the next connected component and removes it from the queue. If the queue is empty the popping thread has to wait.
     */
    ConnectedComponent  pop();

private:

    std::queue<ConnectedComponent> queue_;

    size_t maxSize_;

    mutable boost::mutex accessMutex_;
    boost::condition_variable conditionVariableEmpty_;
    boost::condition_variable conditionVariableFull_;

};

}   // namespace

#endif // VRN_CONNECTEDCOMPONENTQUEUE 
