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

#include "octreequantificationnodequeue.h"

namespace voreen {

OctreeQuantificationNodeQueue::OctreeQuantificationNodeQueue(size_t maxSize)
        : maxSize_(maxSize)
{ }

void OctreeQuantificationNodeQueue::push(std::pair<tgt::IntBounds, const VolumeOctreeNode*> node) {
    boost::mutex::scoped_lock lock(accessMutex_);
    while (maxSize_ && (queue_.size() == maxSize_)) {
        conditionVariableFull_.wait(lock);
    }
    queue_.push(node);
    lock.unlock();
    conditionVariableEmpty_.notify_all();
}

bool OctreeQuantificationNodeQueue::empty() const {
    boost::mutex::scoped_lock lock(accessMutex_);
    return queue_.empty();
}

std::pair<tgt::IntBounds, const VolumeOctreeNode*> OctreeQuantificationNodeQueue::pop() {
    boost::mutex::scoped_lock lock(accessMutex_);
    while(queue_.empty()) {
        conditionVariableEmpty_.wait(lock);
    }                                     
    std::pair<tgt::IntBounds, const VolumeOctreeNode*> node = queue_.front();
    queue_.pop();
    lock.unlock();
    conditionVariableFull_.notify_all();
    return node;
}

}   // namespace
