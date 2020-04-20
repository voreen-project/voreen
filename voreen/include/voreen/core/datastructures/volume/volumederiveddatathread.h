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

#ifndef VRN_VOLUMEDERIVEDDATATHREAD_H
#define VRN_VOLUMEDERIVEDDATATHREAD_H

#include "voreen/core/datastructures/volume/volumederiveddata.h"

namespace voreen {

class VolumeBase;

/// Thread calculating a VolumeDerivedData, started by a Volume when calling getDerivedDataThreaded<T>().
class VolumeDerivedDataThreadBase {
public:
    VolumeDerivedDataThreadBase() : result_(0), volume_(0) {}

    virtual ~VolumeDerivedDataThreadBase() {
        if(isRunning()) {
            interrupt();
            join();
        }
    }

    void startThread(const VolumeBase* vb) {
        if(isRunning()) {
            tgtAssert(false, "Thread is already running!\n");
            return;
        }

        volume_ = vb;
        workerThread_ = boost::thread(&VolumeDerivedDataThreadBase::run, this);
    }

    void join() {
        workerThread_.join();
    }

    void interrupt() {
        workerThread_.interrupt();
    }

    bool isRunning() {
        if(workerThread_.joinable() && !workerThread_.timed_join(boost::posix_time::seconds(0)))
            return true;
        else
            return false;
    }

    VolumeDerivedData* getResult() {
        resultMutex_.lock();
        VolumeDerivedData* tmp = result_;
        resultMutex_.unlock();
        return tmp;
    }

protected:
    virtual void run() = 0;

    VolumeDerivedData* result_;
    boost::mutex resultMutex_;

    const VolumeBase* volume_;

    boost::thread workerThread_;
};

// Template class ----------------------------------------------------------------------------------

template <class T>
class VolumeDerivedDataThread : public VolumeDerivedDataThreadBase {
public:
    virtual void run() {
        T dummy;
        VolumeDerivedData* tmp = 0;
        try {
            tmp = dummy.createFrom(volume_);
        }
        catch(boost::thread_interrupted&)
        {
            tmp = 0;
            resultMutex_.lock();
            result_ = 0;
            resultMutex_.unlock();
            return;
        }
        resultMutex_.lock();
        result_ = tmp;
        resultMutex_.unlock();

        notifyVolume();
    }

protected:
    void notifyVolume();
};

} // namespace

// template definitions ---------------------------------------------------------------------------------

#include "voreen/core/datastructures/volume/volumebase.h"

namespace voreen {

template <class T>
void VolumeDerivedDataThread<T>::notifyVolume() {
    volume_->derivedDataThreadFinished<T>(this);
}

} // namespace

#endif
