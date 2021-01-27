/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "jobqueue.h"

#include <algorithm>
#include <iostream>

using namespace voreen;
using std::min;
using boost::thread;
using boost::unique_lock;
using boost::mutex;
using std::swap;
using std::int64_t;

// define global variables
JobQueue voreen::globalJobQueue;

voreen::JobQueue::JobQueue(){
    int cpucount = boost::thread::hardware_concurrency();
    init(cpucount);
}

voreen::JobQueue::JobQueue( int threads ){
    init(threads);
}

voreen::JobQueue::~JobQueue(){
    {
        unique_lock<mutex> ul(lock_);
        jobs_ = std::queue<jobrecord>();// empty queue
        should_exit_ = true;
        wakeThreads_.notify_all();
    }
    waitAll();
    for(auto &t: threadPool_){
        t.join();
    }
}

void voreen::JobQueue::addJob( Job* j ){
    unique_lock<mutex> ul(lock_);
    jobrecord r;
    r.job = j;
    r.id = current_id_++;
    jobs_.push(r);
    wakeThreads_.notify_all();
}

voreen::JobQueueSync voreen::JobQueue::getSyncPoint( ){
    unique_lock<mutex> ul(lock_);
    return current_id_++;
}

voreen::JobQueueSync voreen::JobQueue::getBarrierSyncPoint(){
    unique_lock<mutex> ul(lock_);
    barriers_.push(current_id_++);
    return current_id_++;
}

void voreen::JobQueue::waitSyncPoint( JobQueueSync id ){
    unique_lock<mutex> ul(lock_);
    std::int64_t minid = getLastFinishedJob();

    // Job already finished
    if (id < minid) return;
    
    // Search sync point
    for(int i = 0; i !=waitingSyncs_.size(); i++){
        syncrecord* sr = waitingSyncs_[i];
        if (id == sr->where){
            // syncpoint found
            sr->waiting++;
            // wait for finishing of job
            while(id > getLastFinishedJob())
                sr->cond.wait(ul);
            sr->waiting--;
            if (!sr->waiting) delete sr;
            return;
        }
    }
    // syncpoint not found. Add a new one to the queue
    syncrecord* sr = new syncrecord;
    sr->where = id;
    sr->waiting = 1;
    waitingSyncs_.push_back(sr);
    while(id > getLastFinishedJob())
        sr->cond.wait(ul);
    sr->waiting--;
    if (!sr->waiting) delete sr;
}

bool voreen::JobQueue::finishedSyncPoint( JobQueueSync id ){
    unique_lock<mutex> ul(lock_);
    std::int64_t minid = getLastFinishedJob();
    return getLastFinishedJob()>id;
}

void voreen::JobQueue::waitAll(){
    waitSyncPoint(getSyncPoint());
}

bool voreen::JobQueue::finishedAll(){
    return finishedSyncPoint(getSyncPoint());
}

void voreen::JobQueue::init( int threads ){
    should_exit_ = false;
    unique_lock<mutex> ul(lock_);
    current_id_ = 0;
    for(int i = 0; i != threads; i++){
        threadPool_.push_back( thread(boost::bind(&JobQueue::execute, this, i)));
        currentJob_.push_back(-1);
    }

}

void voreen::JobQueue::execute(int threadid){
    while(!should_exit_){
        jobrecord rec;
        {
            unique_lock<mutex> ul(lock_);
            int64_t id = getLastFinishedJob();

            // clear finished barriers
            while(!barriers_.empty() && barriers_.front() < id){
                barriers_.pop();
                wakeThreads_.notify_all();
            }

            int64_t barrier = barriers_.empty()? current_id_+1: barriers_.front();

            // Get a job!
            if (jobs_.empty() || jobs_.front().id > barrier){
                wakeThreads_.wait(ul);
                continue;
            }else{
                rec = jobs_.front();
                currentJob_.at(threadid) = rec.id;
                jobs_.pop();
            }
        }

        // Do the job, without holding locks
        rec.job->exec();

        // reacquire locks
        unique_lock<mutex> ul2(lock_);
        // Wake all waiting syncs
        currentJob_.at(threadid) = -1;

        int64_t id = getLastFinishedJob();
        std::vector<syncrecord*> s;
        size_t len = waitingSyncs_.size();
        for(size_t i = 0; i  < len; i++){
            syncrecord *sr = waitingSyncs_[i];
            if (sr->where < id){
                s.push_back(sr);
                swap(waitingSyncs_[i], waitingSyncs_[--len]);
            }
        }
        waitingSyncs_.resize(len);
        for(int i = 0; i != s.size(); i++){
            s[i]->cond.notify_all();
        }
    }
}

std::int64_t voreen::JobQueue::getLastFinishedJob(){
    std::int64_t minid=current_id_;
    if (!jobs_.empty()){
        minid = min(jobs_.front().id, minid);
    }

    for(int i = 0; i != currentJob_.size(); i++){
        std::int64_t id = currentJob_[i];
        if (id != -1){
            minid = min(id, minid);
        }
    }
    return minid;

}


