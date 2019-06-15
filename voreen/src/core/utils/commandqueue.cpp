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

#include "voreen/core/utils/commandqueue.h"

#include "voreen/core/datastructures/callback/callback.h"

#include "tgt/assert.h"

#include "boost/thread/lock_guard.hpp"

namespace voreen {

CommandQueue::CommandQueueElement::CommandQueueElement(const void* owner, std::unique_ptr<Callback>&& callback)
    : owner_(owner)
    , callback_(std::move(callback))
    , mutex_()
{
}
CommandQueue::CommandQueueElement::CommandQueueElement(CommandQueueElement&& other)
    : owner_(other.owner_)
    , callback_(std::move(other.callback_))
    , mutex_()
{
}
CommandQueue::CommandQueueElement& CommandQueue::CommandQueueElement::operator=(CommandQueueElement&& other) {
    owner_ = other.owner_;
    callback_ = std::move(other.callback_);
    return *this;
}

CommandQueue::CommandQueue()
{}

CommandQueue::~CommandQueue() {
    boost::lock_guard<boost::shared_mutex> lock(queueMutex_);
    boost::lock_guard<boost::shared_mutex> nlock(nextMutex_);
    for(auto& command : commands_) {
        tgtAssert(!command.callback_, "Unexecuted commands are left in command queue");
    }
    for(auto& command : nextCommands_) {
        tgtAssert(!command.callback_, "Unexecuted commands are left in next command queue");
    }
}

void CommandQueue::executeAll() {
    // Read access to queue in order to read (and execute) all callbacks
    boost::upgrade_lock<boost::shared_mutex> rlock(queueMutex_);

    for(auto& elm : commands_) {
        boost::lock_guard<boost::mutex> lock(elm.mutex_);
        if(elm.callback_) {
            elm.callback_->exec();
        }
    }

    // Upgrade to write access to queue (and grab write access for nextQueue) in order to
    // clear and move elements from next to the current queue
    boost::upgrade_to_unique_lock<boost::shared_mutex> wlock(rlock);
    boost::unique_lock<boost::shared_mutex> nlock(nextMutex_);
    std::swap(commands_, nextCommands_);
    nextCommands_.clear();
}

void CommandQueue::enqueue(const void* owner, const Callback& command) {
    // Write access to next queue in order to push a new element
    boost::unique_lock<boost::shared_mutex> lock(nextMutex_);

    nextCommands_.push_back(CommandQueueElement(owner, std::unique_ptr<Callback>(command.clone())));
}

void CommandQueue::removeAll(const void* owner) {
    // Read access to current and next queue in order to clear callbacks associated with owner
    boost::shared_lock<boost::shared_mutex> lock(queueMutex_);
    boost::shared_lock<boost::shared_mutex> nlock(nextMutex_);

    if (owner) {
        for(auto& elm : commands_) {
            if (elm.owner_ == owner) {
                boost::lock_guard<boost::mutex> lock(elm.mutex_);
                elm.callback_.reset();
            }
        }
        for(auto& elm : nextCommands_) {
            if (elm.owner_ == owner) {
                boost::lock_guard<boost::mutex> lock(elm.mutex_);
                elm.callback_.reset();
            }
        }
    } else {
        for(auto& elm : commands_) {
            boost::lock_guard<boost::mutex> lock(elm.mutex_);
            elm.callback_.reset();
        }
        for(auto& elm : nextCommands_) {
            boost::lock_guard<boost::mutex> lock(elm.mutex_);
            elm.callback_.reset();
        }
    }
}

} // namespace
