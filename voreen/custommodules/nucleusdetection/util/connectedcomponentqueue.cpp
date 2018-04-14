/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "connectedcomponentqueue.h"

namespace voreen {

ConnectedComponentQueue::ConnectedComponentQueue(size_t maxSize)
        : maxSize_(maxSize)
{ }

void ConnectedComponentQueue::push(ConnectedComponentQueue::ConnectedComponent component) {
    boost::mutex::scoped_lock lock(accessMutex_);
    while (maxSize_ && (queue_.size() == maxSize_)) {
        conditionVariableFull_.wait(lock);
    }
    queue_.push(component);
    lock.unlock();
    conditionVariableEmpty_.notify_all();
}

bool ConnectedComponentQueue::empty() const {
    boost::mutex::scoped_lock lock(accessMutex_);
    return queue_.empty();
}

ConnectedComponentQueue::ConnectedComponent ConnectedComponentQueue::pop() {
    boost::mutex::scoped_lock lock(accessMutex_);
    while(queue_.empty()) {
        conditionVariableEmpty_.wait(lock);
    }                                     
    ConnectedComponent component = queue_.front();
    queue_.pop();
    lock.unlock();
    conditionVariableFull_.notify_all();
    return component;
}

}   // namespace
