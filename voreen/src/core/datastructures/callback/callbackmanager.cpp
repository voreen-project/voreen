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

#include "voreen/core/datastructures/callback/callbackmanager.h"
#include <cstdlib>

namespace voreen {

CallbackManager::CallbackManager()
    : blockCallbacks_(false)
{
}

CallbackManager::~CallbackManager(){
    for(size_t i = 0; i != callbacks_.size(); i++){
        delete callbacks_[i];
}

}
void CallbackManager::registerCallback( const Callback& a ) {
    callbacks_.push_back(a.clone());
}

void CallbackManager::execute() {
    //do nothing, if blocked
    if(blockCallbacks_) return;

    for(size_t i = 0; i != callbacks_.size(); i++){

        callbacks_[i]->exec();
    }
}

void CallbackManager::blockCallbacks(bool block) {
    blockCallbacks_ = block;
}

}   // namespace
