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

#ifndef VRN_JOBQUEUE_H
#define VRN_JOBQUEUE_H

#include <vector>
#include <queue>
#include <cstdint>

#include <boost/thread.hpp>

namespace voreen{
    /**
     * Interface for jobs
     */
    class Job{
    public:
        /**
         * Function that is called for the job.
         */
        virtual void exec() = 0;

    };

    template<typename T>
    class CMFunctionJob: public Job{
    public:
        typedef void (*func)(T);
        CMFunctionJob(func f, T t){
            fun = f;
            data = t;
        }

        virtual void exec(){
            fun(data);

        }

    private:
        T data;
        func fun;
    };

    /**
     * A sync point for a job queue.
     * There is a default Queue with the name "globalCMJobQueue"
     *
     */
    typedef std::int64_t JobQueueSync;

    /**
     * A parallel Queue of Jobs/Tasks.
     * Added Jobs are executed in parallel in a threadpool
     *
     * This structure is threadsafe. Jobs can be added out of
     * already running jobs on the same queue and from multiple
     * threads.
     *
     */
    class JobQueue{
    public:
        /**
         * This Constructor creates the Jobqueue with the as many threads
         * as processors are availible.
         */
        JobQueue();
        /**
         * This Constructor creates the Jobqueue with a custom number of threads.
         * @param threads The number of threads to use.
         */
        JobQueue(int threads);
        ~JobQueue();

        /**
         * Adds a jobs for concurrent executation. The jobs are usually not executed in
         * in the sequence they are added.
         * @param The job to do. "delete" is called on this, after the job is finished
         */
        void addJob(Job* job);

        template<typename T>
        void addFunctionJob(typename CMFunctionJob<T>::func, T);

        /**
         * Get a syncpoint for the for waiting later if all jobs added before are finished.
         *
         * This syncpoint does NOT block.
         */
        JobQueueSync getSyncPoint();
        /**
         * Get a blocking syncpoint for the for waiting later if all jobs added before are finished.
         *
         * This point blocks the queue, so that jobs added after it are only executed after all
         * jobs added before are finished.
         */
        JobQueueSync getBarrierSyncPoint();

        /**
         * Wait till all jobs added before a syncpoint are finished
         * @param sync The syncpoint
         */
        void         waitSyncPoint(JobQueueSync sync);
        /**
         * Check if all jobs added before a syncpoint are finished
         * @param sync The syncpoint
         * @return true if all jobs added before "sync" are finished
         */
        bool         finishedSyncPoint(JobQueueSync sync);

        /**
         * Wait till all jobs are finished.
         */
        void         waitAll();
        /**
         * Check if all jobs are finished
         */
        bool         finishedAll();
    private:
        struct jobrecord{
            std::int64_t id;
            Job *job;
        };

        struct syncrecord{
            std::int64_t where;
            boost::condition_variable cond;
            int waiting;
        };


        std::int64_t current_id_; // = 0
        bool should_exit_; // = false;
        std::vector<boost::thread> threadPool_;
        boost::mutex lock_;
        boost::condition_variable wakeThreads_;
        std::queue<jobrecord> jobs_;
        std::vector<std::int64_t> currentJob_;
        std::queue<std::int64_t> barriers_;

        std::vector<syncrecord*> waitingSyncs_;

        void init(int threads);
        void execute(int id);
        std::int64_t getLastFinishedJob(); // needs lock

        JobQueue (const JobQueue &);
        JobQueue & operator = (const JobQueue &);
    };

    template<typename T>
    void voreen::JobQueue::addFunctionJob( typename CMFunctionJob<T>::func f, T t)
    {
        addJob(new CMFunctionJob<T>(f, t));
    }

    extern JobQueue globalJobQueue;
};
#endif //VRN_JOBQUEUE_H
