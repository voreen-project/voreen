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

#include "voreen/core/utils/commandqueue.h"

#include "voreen/core/datastructures/callback/callback.h"

#include "tgt/assert.h"

#include "boost/thread/lock_guard.hpp"

namespace voreen {

CommandQueue::CommandQueue()
{}

CommandQueue::~CommandQueue() {
    boost::lock_guard<boost::recursive_mutex> lock(mutex_);
    tgtAssert(commands_.empty(), "Some unexecuted commands are left in command queue");
}

void CommandQueue::execute() {
    boost::lock_guard<boost::recursive_mutex> lock(mutex_);

    if (!commands_.empty()) {

        // Execute callback.
        Callback* callback = commands_.front().second.get();
        callback->exec();
        commands_.pop_front();
    }
}

void CommandQueue::executeAll() {
    boost::lock_guard<boost::recursive_mutex> lock(mutex_);

    while (!commands_.empty()) {
        Callback* callback = commands_.front().second.get();
        callback->exec();
        commands_.pop_front();
    }
}

void CommandQueue::enqueue(void* owner, const Callback& command) {
    boost::lock_guard<boost::recursive_mutex> lock(mutex_);

    commands_.push_back(std::make_pair(owner, std::unique_ptr<Callback>(command.clone())));
}

void CommandQueue::removeAll(void* owner) {
    boost::lock_guard<boost::recursive_mutex> lock(mutex_);

    if (owner) {
        for (auto iter = commands_.begin(); iter != commands_.end(); ) {
            if (iter->first == owner)
                iter = commands_.erase(iter);
            else
                iter++;
        }
    }
    else
        commands_.clear();
}

} // namespace
