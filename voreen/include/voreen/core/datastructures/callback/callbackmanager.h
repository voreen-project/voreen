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

#ifndef VRN_CALLBACKMANAGER_H
#define VRN_CALLBACKMANAGER_H

#include <vector>

#include "voreen/core/datastructures/callback/callback.h"

namespace voreen {


/**
 * A CallbackManager manages a list of callbacks, it can execute later.
 *
 * It is not possible to to remove callbacks, so you have to make sure,
 * the CallbackManager is not executed after data referenced by
 * registered callbacks have been deleted.
 */
class CallbackManager {
public:

    CallbackManager();

    ~CallbackManager();

    /**
     * Register a Callback for later executate.
     */
    void registerCallback(const Callback& a);

    /**
     * Execute all registered Callbacks in the order they were registered for
     * their sideffects
     */
    void execute();

    /** If true, execute does nothing. */
    void blockCallbacks(bool block);

private:
    std::vector<Callback*> callbacks_;  ///< vector of registered callbacks
    bool blockCallbacks_;               ///< if true, execute() does nothing
};

// ----------------------------------------------------------------------------


} // namespace voreen

#endif // VRN_CALLBACKMANAGER_H
