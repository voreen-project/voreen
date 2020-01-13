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

#ifndef VRN_OCTREEQUANTIFICATIONNODEQUEUE
#define VRN_OCTREEQUANTIFICATIONNODEQUEUE

#include "voreen/core/voreencoreapi.h"

#include <queue>

#include <boost/thread.hpp>
#include <boost/thread/condition_variable.hpp>

#include "voreen/core/datastructures/octree/volumeoctree.h"

namespace voreen {

/**
 * Provides a threadsafe data structure for queueing octree nodes that have to be quantified. Use with OctreeQuantificationThread class.
 */
class VRN_CORE_API OctreeQuantificationNodeQueue {

public:

    /**
     * Creates a queue with a fized maximum size. 
     *
     * @param maxSize determines the maximum number of nodes that can simultaneously be queued. If set to 0, the queue does not have a maximum size. 
     */
    OctreeQuantificationNodeQueue(size_t maxSize = 0);
                    

    /**
     * Pushes a node onto the queue. If the queue is full the pushing thread has to wait until the size of the queue is lower than maxSize.
     */
    void push(std::pair<tgt::IntBounds, const VolumeOctreeNode*> node);

    /**
     * Returns if the queue is empty.
     */
    bool empty() const;

    /**
     * Returns the next node and removes it from the queue. If the queue is empty the popping thread has to wait.
     */
    std::pair<tgt::IntBounds, const VolumeOctreeNode*>  pop();

private:

    std::queue<std::pair<tgt::IntBounds, const VolumeOctreeNode*> > queue_;

    size_t maxSize_;

    mutable boost::mutex accessMutex_;
    boost::condition_variable conditionVariableEmpty_;
    boost::condition_variable conditionVariableFull_;

};

}   // namespace

#endif // VRN_OCTREEQUANTIFICATIONNODEQUEUE 
