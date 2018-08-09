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

#ifndef VRN_COMMANDQUEUEPARSER_H
#define VRN_COMMANDQUEUEPARSER_H

#include "voreen/core/voreencoreapi.h"

#include <boost/thread.hpp>

#include <vector>
#include <map>

namespace voreen {

class Callback;

struct CommandQueueElement {
    CommandQueueElement(const void* owner, std::unique_ptr<Callback>&& callback);
    CommandQueueElement(CommandQueueElement&&);
    CommandQueueElement& operator=(CommandQueueElement&&);
    const void* owner_;
    std::unique_ptr<Callback> callback_;
    boost::mutex mutex_; // Mutual exclusion for execution of callback and deletion of the element
};

/**
 * Command queue for executing deferred tasks.
 * This can come in handy for tasks being started in another than the main thread.
 * These tasks are going to be executed each timer event, right before any network evaluation.
 *
 * @see VoreenApplication
 */
class VRN_CORE_API CommandQueue {
public:

    /// Constructor.
    CommandQueue();

    /// Destructor.
    ~CommandQueue();

    /**
     * Executes all commands in the queue in the order of insertion and removes them.
     */
    void executeAll();

    /**
     * Enqueues a single command and registers it's owner.
     *
     * @Note: Owners need to keep track of these tasks!
     */
    void enqueue(const void* owner, const Callback& command);

    /**
     * Removes All tasks owned by the given owner.
     *
     * @param owner if null, every command is removed
     */
    void removeAll(const void* owner = nullptr);

protected:

    boost::shared_mutex queueMutex_; // Either: read (shared) or modify (exclusive)
    boost::mutex nextMutex_; // standard mutex for the commands that will be added in the _next_ executeAll
    std::vector<CommandQueueElement> commands_; ///< Actual command queue, defines order of execution
    std::vector<CommandQueueElement> nextCommands_; ///< Commands that will be added to the command queue soon

private:

    // Delete copy constructor and assignment operator, which enables use of unique_ptr for members.
    CommandQueue(const CommandQueue&);
    CommandQueue& operator=(const CommandQueue&);

};


} // namespace

#endif // VRN_COMMANDQUEUE_H
