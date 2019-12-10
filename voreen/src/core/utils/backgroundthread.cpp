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

#include "voreen/core/utils/backgroundthread.h"

namespace voreen {

    /*
     * BackgroundThread
     * ----------------*/

    BackgroundThread::~BackgroundThread() {
        internalThread_.join(); //make sure the internal thread is finished before destroying other class data
    }

    void BackgroundThread::run() {
        //make sure that the background thread is not currently running
        if (isRunning()) {
            internalThread_.interrupt();
            internalThread_.join();
            running_ = false;
        }
        finished_ = false;
        running_ = true;
        //start computation
        internalThread_ = boost::thread(boost::bind(&BackgroundThread::execute, this));
    }

    void BackgroundThread::join() {
        internalThread_.join();
        running_ = false;
    }

    void BackgroundThread::interrupt() {
        internalThread_.interrupt();
    }

    void BackgroundThread::interruptAndJoin() {
        internalThread_.interrupt();
        internalThread_.join();
        running_ = false;
    }

    void BackgroundThread::execute() {
        try {
            threadMain();
            running_ = false;
            finished_ = true;
        }
        catch (boost::thread_interrupted& /*interruption*/) {
            //thread has been interrupted: handle the interruption
            handleInterruption();
            running_ = false;
            return;
        }
    }

} //namespace
