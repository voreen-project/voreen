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

#ifndef VRN_ASYNCCOMPUTEPROCESSOR_H
#define VRN_ASYNCCOMPUTEPROCESSOR_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/utils/backgroundthread.h"
#include "voreen/core/utils/commandqueue.h"
#include "voreen/core/ports/port.h"
#include "voreen/core/ports/coprocessorport.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/datastructures/volume/volumeobserver.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/voreenapplication.h"

#include "boost/optional.hpp"

#include <exception>
#include <chrono>
#include <mutex>
#include <cstdlib>

namespace voreen {

class PropertyDisabler {
public:
    PropertyDisabler(PropertyOwner& owner)
        : owner_(owner)
        , propertyWritablilityMap_()
    {
    }

    // Save the current state of the properties as the state to be restored when calling `restore`.
    //
    // This function should be called at the very end of a `PropertyOwner`-constructor  in most cases.
    void saveState(std::function<bool(Property*)> shouldBeIgnored) {
        for(Property* property : owner_.getProperties()) {
            if(!shouldBeIgnored(property)) {
                propertyWritablilityMap_[property] = property->isReadOnlyFlagSet();
            }
        }
    }

    void disable() {
        for(auto pair : propertyWritablilityMap_) {
            pair.first->setReadOnlyFlag(true);
        }
    }

    void restore() {
        for(auto pair : propertyWritablilityMap_) {
            pair.first->setReadOnlyFlag(pair.second);
        }
    }

private:
    PropertyOwner& owner_;
    std::map<Property*, bool> propertyWritablilityMap_;
};

/**
 * Can be thrown to indicate invalid input configurations.
 * See AsyncComputeProcessor::prepareComputeInput()
 */
class InvalidInputException : public std::exception {
public:
    /**
     * Indicates the severity of the error.
     * Both S_ERROR and S_WARNING will inform the user about the problem
     * using the message provided in the constructor.
     * S_IGNORE will cause the warning to be ignored silently.
     */
    enum Severity {
        S_ERROR,
        S_WARNING,
        S_IGNORE
    };
    InvalidInputException(const std::string& msg, Severity);
    virtual ~InvalidInputException() {}
    virtual const char* what() const throw();

    const Severity severity_;
    const std::string msg_;
};

/**
 * This is a processor template for asynchronous computation.
 * Processors whose computation task takes more than a few seconds can (and should!)
 * derive from this class.
 * Instead of performing the heavy computation task in process(), AsyncComputeProcessors
 * work in the compute()-method, which is executed in another thread.
 * Consequently, compute() should only operate on local data and use the provided
 * ProgressReporter instead of using the processors setProgress etc. methods.
 *
 * Both prepareComputeInput() and processComputeOutput() are executed in the main
 * thread and should be used to generate and process the in- and output of compute, respectively.
 *
 * AsyncComputeProcessor also conveniently provides controls and progress displays
 * for the asynchronous computation in a separate property group.
 *
 * ComputeInput and ComputeOutput are template parameter classes that are used to pass
 * data from and to the compute thread. The only requirement to classes used as these is
 * that both most provide a move-constructor.
 */
template<typename I, typename O>
class AsyncComputeProcessor : public Processor,
                              public PortObserver,
                              public DataInvalidationObserver {
protected:

    typedef boost::thread_interrupted InterruptionException;
    static void interruptionPoint();

    typedef O ComputeOutput;
    typedef I ComputeInput;


public:
    AsyncComputeProcessor();
    virtual ~AsyncComputeProcessor();

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    virtual void deserialize(Deserializer& s);

    bool usesExpensiveComputation() const {
        return true;
    }

protected:

    /**
     * Check the input of the processor and generate a ComputeInput object that will be passed
     * to the compute() thread.
     * InvalidInputException can be thrown, if an error occured during input processing.
     * If no exception is thrown, the compute thread will be started with the return value of
     * this function as the input parameter.
     */
    virtual ComputeInput prepareComputeInput() = 0;

    /**
     * Performn the main task of this processor, using input to generate the returned output.
     * The provided ProgressReporter can be used to savely pass progress information to the
     * main thread.
     * The setProgress method of the progress reporter includes an interruption point
     * that is used to interrupt the computation via the main thread.
     * The thread_interrupted exception should only be caught to perform cleanup tasks and
     * should then be rethrown immediately.
     */
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const = 0;

    /**
     * Receives the results of the compute thread, if completed successfully and
     * processes the output accordingly.
     * In this method the computation result can be used to e.g. set outports,
     * properties, etc. using the data generated in compute().
     */
    virtual void processComputeOutput(ComputeOutput output) = 0;

    /**
     * Get references to all critical ports, i.e. all ports whose data has
     * to be valid for the duration of a computation (i.e. during compute).
     * Any change to these ports will stop the computation.
     * Note: This function is expected to always return references to the same ports,
     * when called between initialize() and deinitialize().
     */
    virtual std::vector<std::reference_wrapper<Port>> getCriticalPorts();

    /**
     * Stops the interrupts the background thread, waits until it was stopped
     * and restores property states to the non-running state
     */
    void stopComputeThread() throw();

    /**
     * Interrupts the computation and sets property states accordingly.
     * In contrast to stopComputeThread, this method should be noly called, if the
     * computation did not exit successfully.
     */
    void interruptComputation() throw();

    /// Inherited from PortObserver
    virtual void afterConnectionAdded(const Port* source, const Port* connectedPort);
    virtual void beforeConnectionRemoved(const Port* source, const Port*);
    virtual void dataWillChange(const Port* source);
    virtual void dataHasChanged(const Port* source);

    /// Inherited from DataInvalidationObserver
    virtual void dataAboutToInvalidate(const DataInvalidationObservable* source);
protected:
    void forceComputation();
    // Use this, if the subclass messes with property visibilities themselves
    // and notify AsyncComputeProcessor when the visibility changes.
    void savePropertyEnableStates();

private:

    enum InvalidationBehavior {
        INVALIDATE_RESTART,
        INVALIDATE_ABORT,
        INVALIDATE_ENQUEUE,
    };


    typedef std::chrono::steady_clock Clock;
    typedef std::chrono::time_point<Clock> TimePoint;

    const static std::chrono::milliseconds MINIMUM_PROGRESS_UPDATE_INTERVAL;

    /// Used to manage the visibility of properties during compute()
    void enableRunningState();
    void disableRunningState();
    void enqueueRestart();

    void trySetUpInvalidationBehavior();

    /**
     * Safely passes progress reports through to the processor in the main thread.
     */
    class ComputeProgressReporter : public ProgressReporter {
    public:
        ComputeProgressReporter(AsyncComputeProcessor<ComputeInput, ComputeOutput>& processor);
        virtual ~ComputeProgressReporter();

        virtual void setProgressMessage(const std::string& message);
        virtual std::string getProgressMessage() const;

        virtual void update();
    private:
        AsyncComputeProcessor<ComputeInput, ComputeOutput>& processor_;
        TimePoint startTime_;
        TimePoint lastUpdate_;
        std::string message_;
    };

    /**
     * Ignores all progress reports, e.g. when running non-interactively.
     */
    class NoopProgressReporter : public ProgressReporter {
    public:
        NoopProgressReporter();
        virtual ~NoopProgressReporter();

        virtual void setProgressMessage(const std::string& message);
        virtual std::string getProgressMessage() const;
    private:
        std::string message_;
    };

    /**
     * Abstraction of the compute thread, including input and output data.
     */
    class ComputeThread : public ProcessorBackgroundThread<AsyncComputeProcessor<ComputeInput, ComputeOutput>> {
    public:
        ComputeThread(AsyncComputeProcessor<ComputeInput, ComputeOutput>& processor);
        virtual ~ComputeThread() {}

        /**
         * Set the input to be moved into the computation thread
         */
        void setInput(ComputeInput&& input);

        /**
         * Retrieve (and remove) the output from the (finished) computation thread.
         * If there is no output, e.g. because the thread has not finished computation,
         * a null pointer is returned, indicating "no output".
         */
        std::unique_ptr<ComputeOutput> retrieveOutput();

        /**
         * Remove the output from the (finished) computation thread.
         */
        void clearOutput();

        void threadMain();

    private:
        std::unique_ptr<ComputeInput> input_;
        std::unique_ptr<ComputeOutput> output_;
        std::mutex mutex_;
    };

    BoolProperty synchronousComputation_;
    OptionProperty<InvalidationBehavior> invalidationBehavior_;
    bool invalidationBehaviorPropertyIsSetUp_;
    ButtonProperty manualUpdateButton_;
    ButtonProperty stopUpdateButton_;
    ProgressProperty progressDisplay_;
    StringProperty statusDisplay_;

    /// Used to save/restore property writability before/after expensive computation
    PropertyDisabler propertyDisabler_;

    bool updateForced_;
    bool stopForced_;
    bool restartEnqueued_;
    bool lastInvalidationWasEnqueue_;

    ComputeThread computation_;
    TimePoint computationStartTime_; //Used for final display of required time
};


// --------------------------------------------------------------------------------
// Implementation -----------------------------------------------------------------
// --------------------------------------------------------------------------------

// AsyncComputeProcessor::ComputeThread ------------------------------------

template<class I, class O>
AsyncComputeProcessor<I, O>::ComputeThread::ComputeThread(AsyncComputeProcessor<I, O>& processor)
    : ProcessorBackgroundThread<AsyncComputeProcessor<I, O>>(&processor)
    , input_(nullptr)
    , output_(nullptr)
{
}

template<class I, class O>
void AsyncComputeProcessor<I, O>::ComputeThread::setInput(I&& input) {
    std::lock_guard<std::mutex> guard(mutex_);
    input_.reset(new I(std::move(input)));
}

template<class I, class O>
std::unique_ptr<O> AsyncComputeProcessor<I, O>::ComputeThread::retrieveOutput() {
    std::lock_guard<std::mutex> guard(mutex_);
    return std::move(output_);
}

template<class I, class O>
void AsyncComputeProcessor<I, O>::ComputeThread::clearOutput() {
    std::lock_guard<std::mutex> guard(mutex_);
    return output_.reset(nullptr);
}

template<class I, class O>
void AsyncComputeProcessor<I, O>::ComputeThread::threadMain() {
    tgtAssert(input_, "ComputeThread started without input!");

    mutex_.lock();
    output_.reset();

    std::unique_ptr<ProgressReporter> progress;
    if(ProcessorBackgroundThread<AsyncComputeProcessor<I,O>>::processor_->synchronousComputation_.get()) {
        progress.reset(new NoopProgressReporter());
    } else {
        progress.reset(new ComputeProgressReporter(*ProcessorBackgroundThread<AsyncComputeProcessor<I, O>>::processor_));
    }
    I* input = input_.release();
    mutex_.unlock();

    O* output = new O(std::move(ProcessorBackgroundThread<AsyncComputeProcessor<I,O>>::processor_->compute(std::move(*input), *progress)));

    mutex_.lock();
    output_.reset(output);
    mutex_.unlock();
}

// AsyncComputeProcessor::ComputeProgressReporter ------------------------------------

template<class I, class O>
AsyncComputeProcessor<I,O>::ComputeProgressReporter::ComputeProgressReporter(AsyncComputeProcessor<I, O>& processor)
    : processor_(processor)
    , startTime_(Clock::now())
    , lastUpdate_(Clock::now())
{
}

template<class I, class O>
AsyncComputeProcessor<I,O>::ComputeProgressReporter::~ComputeProgressReporter()
{
    VoreenApplication::app()->getCommandQueue()->removeAll(this);
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::ComputeProgressReporter::setProgressMessage(const std::string& message) {
    message_ = message;
}

template<class I, class O>
std::string AsyncComputeProcessor<I,O>::ComputeProgressReporter::getProgressMessage() const {
    return message_;
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::ComputeProgressReporter::update() {
    float progress = getProgress();
    if(progress != 0) {
        TimePoint now = Clock::now();

        auto timeSinceLastUpdate = now - lastUpdate_;
        auto timeSinceLastUpdateMillis = std::chrono::duration_cast<std::chrono::milliseconds>(timeSinceLastUpdate);
        if(timeSinceLastUpdateMillis < MINIMUM_PROGRESS_UPDATE_INTERVAL) {
            return;
        }
        lastUpdate_ = now;

        auto computeDuration = now - startTime_;
        auto expectedTotalDuration = computeDuration/progress;
        auto remainingDuration = expectedTotalDuration - computeDuration;
        auto remainingDurationMillis = std::chrono::duration_cast<std::chrono::milliseconds>(remainingDuration);
        std::string timeFormat(formatTime(remainingDurationMillis.count()));

        // Enqueue new command for an ui update.
        VoreenApplication::app()->getCommandQueue()->enqueue(this, LambdaFunctionCallback([this, timeFormat, progress] {
                        std::string msg = timeFormat + " remaining";
                        processor_.statusDisplay_.set(msg);
                        processor_.setProgress(progress);
                    }));
    }

    interruptionPoint();
}

// AsyncComputeProcessor::NoopProgressReporter ---------------------------------------
template<class I, class O>
AsyncComputeProcessor<I,O>::NoopProgressReporter::NoopProgressReporter()
    : message_()
{
}

template<class I, class O>
AsyncComputeProcessor<I,O>::NoopProgressReporter::~NoopProgressReporter()
{
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::NoopProgressReporter::setProgressMessage(const std::string& message) {
    message_ = message;
}

template<class I, class O>
std::string AsyncComputeProcessor<I,O>::NoopProgressReporter::getProgressMessage() const {
    return message_;
}


// AsyncComputeProcessor ---------------------------------------------------------
template<class I, class O>
void AsyncComputeProcessor<I,O>::interruptionPoint() {
    boost::this_thread::interruption_point();
}

template<class I, class O>
const std::chrono::milliseconds AsyncComputeProcessor<I,O>::MINIMUM_PROGRESS_UPDATE_INTERVAL(100);

template<class I, class O>
AsyncComputeProcessor<I,O>::AsyncComputeProcessor()
    : Processor()
    , synchronousComputation_("synchronousComputation", "Wait for Result", false, Processor::INVALID_RESULT, Property::LOD_ADVANCED)
    , invalidationBehavior_("invalidationMode", "Invalidation Behavior", Processor::VALID, false, Property::LOD_ADVANCED)
    , invalidationBehaviorPropertyIsSetUp_(false)
    , manualUpdateButton_("manualUpdateButton_", "Start", Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , stopUpdateButton_("stopUpdateButton", "Stop", Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , progressDisplay_("progressDisplay", "Progress")
    , statusDisplay_("statusDisplay", "Status", "Stopped", Processor::VALID)
    , updateForced_(false)
    , stopForced_(false)
    , restartEnqueued_(false)
    , computation_(*this)
    , propertyDisabler_(*this)
{
    addProperty(synchronousComputation_);
        synchronousComputation_.setGroupID("ac_processing");
    addProperty(invalidationBehavior_);
        invalidationBehavior_.setGroupID("ac_processing");
        // Options are added at a later stage! See trySetUpInvalidationBehavior for more info.

    addProperty(manualUpdateButton_);
        manualUpdateButton_.setGroupID("ac_processing");
        ON_CHANGE_LAMBDA(manualUpdateButton_, [this] () {
                if(!isReady()) {
                    std::string msg = "Cannot start computing: " + getNotReadyErrorMessage();
                    VoreenApplication::app()->showMessageBox("Cannot start computing" , msg, true);
                    LERROR(msg);
                    return;
                }
                updateForced_ = true;
                });
    addProperty(stopUpdateButton_);
        stopUpdateButton_.setGroupID("ac_processing");
        stopUpdateButton_.setVisibleFlag(false);
        ON_CHANGE_LAMBDA(stopUpdateButton_, [this] () {
                stopForced_ = true;
                });
    addProperty(progressDisplay_);
        progressDisplay_.setGroupID("ac_processing");
    addProperty(statusDisplay_);
        statusDisplay_.setGroupID("ac_processing");
        statusDisplay_.setEditable(false);
    setPropertyGroupGuiName("ac_processing", "Processing");

    addProgressBar(&progressDisplay_);
}

template<class I, class O>
AsyncComputeProcessor<I,O>::~AsyncComputeProcessor() {
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::savePropertyEnableStates() {
    propertyDisabler_.saveState([this] (Property* p) {
    return p == &manualUpdateButton_
        || p == &stopUpdateButton_
        || p == &progressDisplay_
        || p == &statusDisplay_;
    });
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::enableRunningState() {
    savePropertyEnableStates();
    propertyDisabler_.disable();

    manualUpdateButton_.setVisibleFlag(false);
    stopUpdateButton_.setVisibleFlag(true);
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::disableRunningState() {
    propertyDisabler_.restore();

    manualUpdateButton_.setVisibleFlag(true);
    stopUpdateButton_.setVisibleFlag(false);
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::enqueueRestart() {
    restartEnqueued_ = true;
    lastInvalidationWasEnqueue_ = true;
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::trySetUpInvalidationBehavior() {
    // In this method we dynamically add the options to the
    // invalidationBehavior_ property. Specifically, we _only_ add the enqueue
    // option if the processor has any ports that allow for the enqueue
    // behavior to have any effect at all (which is often not the case!).
    //
    // We need to add the enqueue option before (potential) property
    // deserialization, but can only check ports after they have been added in
    // the child class constructor. AsyncComputeProcessor::deserialize (and
    // specifically before calling Processor::deserialize() is the _only_ place
    // were those two conditions hold (even though this doesn't really have
    // anything to do with deserialization).
    // Another complication is that if a processor is newly created, it is not
    // deserialized, so we dynamically have to decide whether or not we still
    // have to add the enqueue options in invalidate (hence the need for
    // invalidationBehaviorPropertyIsSetUp_).

    if(!invalidationBehaviorPropertyIsSetUp_) {
        bool enqueueMeaningful = false;
        for(auto port : getCriticalPorts()) {
            enqueueMeaningful |= dataSafeToUseAfterInvalidation(&port.get());
        }
        invalidationBehavior_.addOption("invalidateAbort",   "Ignore/Abort", INVALIDATE_ABORT);
        invalidationBehavior_.addOption("invalidateRestart", "Start/Restart", INVALIDATE_RESTART);
        if(enqueueMeaningful) {
            invalidationBehavior_.addOption("invalidateEnqueue", "Start/Enqueue (if possible)", INVALIDATE_ENQUEUE);
        }
        invalidationBehaviorPropertyIsSetUp_ = true;
    }
}

template<class I, class O>
std::vector<std::reference_wrapper<Port>> AsyncComputeProcessor<I,O>::getCriticalPorts() {
    std::vector<std::reference_wrapper<Port>> criticalPorts;
    for(Port* port : getInports()) {
        criticalPorts.push_back(std::ref(*port));
    }
    for(Port* port : getCoProcessorInports()) {
        criticalPorts.push_back(std::ref(*port));
    }
    return criticalPorts;
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::afterConnectionAdded(const Port*, const Port*) {
    // Relevant connection remove/add events also trigger dataHasChanged/dataWillChange callbacks,
    // so no action is required here.
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::beforeConnectionRemoved(const Port*, const Port*) {
    // Relevant connection remove/add events also trigger dataHasChanged/dataWillChange callbacks,
    // so no action is required here.
}

static bool dataSafeToUseAfterInvalidation(const Port* source) {
    // Here the assumption is: We CAN use the port data in a threaded context,
    // BUT we are not able (and thus required) to watch for the invalidation of
    // the data. Hence the only remaining possibility is that data is _copied_
    // (see GenericPort, PortDataPointer, Cloneable) or static and can thus
    // still be used even though the port data itself has changed.
    return !source->isDataInvalidationObservable();
}

template<class I, class O>
void AsyncComputeProcessor<I, O>::dataWillChange(const Port* source) {
    if(invalidationBehavior_.getValue() == INVALIDATE_ENQUEUE && dataSafeToUseAfterInvalidation(source)) {
        enqueueRestart();
    } else {
        interruptComputation();
        // Clear the result, if we have one already, as it will not be valid (i.e., corresponding to the input) afterwards.
        computation_.clearOutput();
    }

    // The port will no longer store it's currently contained data.
    // The data itself might continue existing, therefore we need to remove the observer.
    if(source->isDataInvalidationObservable() && source->hasData()) {
        source->removeDataInvalidationObserver(static_cast<DataInvalidationObserver*>(this));
    }
}

template<class I, class O>
void AsyncComputeProcessor<I, O>::dataHasChanged(const Port* source) {
    if(invalidationBehavior_.getValue() == INVALIDATE_ENQUEUE && dataSafeToUseAfterInvalidation(source)) {
        enqueueRestart();
    } else {
        interruptComputation();
    }

    // The port stores new data.
    // The data might change later, so we register ourselves as an observer if this is possible
    if(source->isDataInvalidationObservable() && source->hasData()) {
        source->addDataInvalidationObserver(static_cast<DataInvalidationObserver*>(this));
    }
}

template<class I, class O>
void AsyncComputeProcessor<I, O>::dataAboutToInvalidate(const DataInvalidationObservable* source) {
    interruptComputation();
}
template<class I, class O>
void AsyncComputeProcessor<I, O>::forceComputation() {
    invalidate();
    updateForced_ = true;
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::initialize() {
    // HACK: (See trySetUpInvalidationBehavior for details)
    trySetUpInvalidationBehavior();

    tgtAssert(invalidationBehaviorPropertyIsSetUp_, "invalidationBehavior_ not set up");

    Processor::initialize();
    // We cannot call getCriticalPorts in constructor, as
    //  1. it is a virtual function
    //  2. the child class did not add any ports so far anyway
    for(auto port : getCriticalPorts()) {
        // We observe all ports
        static_cast<Observable<PortObserver>&>(port.get()).addObserver(static_cast<PortObserver*>(this));

        // The final check tests if the data is thread safe as it. That's true in case of integral types such like std::string, int, ..
        // In that case there is no need to observe the data since it's always passed by value and can't change itself.
        if (!port.get().isDataThreadSafe()) {
            // To get rid of this warning, implement a *PortType*Observer and make AsyncComputeProcessor derive from it
            // and mimic the handling of VolumePorts.
            // In case your port data is thread safe (i.e. by being immutable), simply override Port::isDataThreadSafe() to return true.
            const char* warningMsg = "Detected critical non-thread safe port! Changing input data during computation may result in crashes and/or undefined results!";
            LWARNING(warningMsg);
            tgtAssert(false, warningMsg);
        }
    }

    // Populate propertyWritablilityMap in case stopComputeThread gets called before any computation
    savePropertyEnableStates();

    if(std::getenv("VRN_ASYNCCOMPUTEPROCESSOR_WAIT_FOR_RESULT")) {
        synchronousComputation_.set(true);
    }
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::deinitialize() {
    // We obviously do not want the thread to keep running if the processor is deinitialized,
    // we will never retrieve the result, anyway.
    stopComputeThread();

    for(auto port : getCriticalPorts()) {
        // Do not observe this port anymore
        static_cast<const Observable<PortObserver>&>(port.get()).removeObserver(static_cast<PortObserver*>(this));

        if (port.get().hasData() && port.get().isDataInvalidationObservable()) {
            port.get().removeDataInvalidationObserver(static_cast<DataInvalidationObserver*>(this));
        }
    }

    Processor::deinitialize();
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::deserialize(Deserializer& s) {
    // HACK: (See trySetUpInvalidationBehavior for details)
    trySetUpInvalidationBehavior();

    Processor::deserialize(s);
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::stopComputeThread() throw() {
    computation_.interruptAndJoin();

    disableRunningState();
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::interruptComputation() throw() {
    stopComputeThread();

    // We cannot call anything that causes QT to process events, as e.g. multiple mouse
    // move events would cause the portconnection to be removed recursively.
    //setProgressRange(tgt::vec2(0.0f, 1.0f));
    //setProgress(0.0f);
    statusDisplay_.set("Interrupted");
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::process() {

    bool running = computation_.isRunning();
    bool enqueueRelatedInvalidation = invalidationBehavior_.getValue() == INVALIDATE_ENQUEUE && lastInvalidationWasEnqueue_;
    bool abortNecessary = (updateForced_ || stopForced_
            || (getInvalidationLevel() >= INVALID_RESULT && !enqueueRelatedInvalidation))
        && running;

    if(abortNecessary) {
        interruptComputation();
    }

    std::unique_ptr<O> result = computation_.retrieveOutput();

    bool restartNecessary = updateForced_
        || (!result && getInvalidationLevel() >= INVALID_RESULT && invalidationBehavior_.getValue() != INVALIDATE_ABORT && !enqueueRelatedInvalidation) && !stopForced_
        || restartEnqueued_ && !running && !stopForced_;

    stopForced_ = false;
    updateForced_ = false;
    lastInvalidationWasEnqueue_ = false;

    if(result) {
        processComputeOutput(std::move(*result));
        stopComputeThread();

        TimePoint now = Clock::now();
        auto computeDuration = now - computationStartTime_;
        auto computeDurationMillis = std::chrono::duration_cast<std::chrono::milliseconds>(computeDuration);

        std::string finishedString = "Finished in ";
        finishedString += formatTime(computeDurationMillis.count());
        statusDisplay_.set(finishedString);

        setProgressRange(tgt::vec2(0.0f, 1.0f));
        setProgress(1.0f);
    }

    if(restartNecessary) {
        tgtAssert(!computation_.isRunning(), "Start on running computation");
        setProgressRange(tgt::vec2(0.0f, 1.0f));
        setProgress(0.0f);
        try {
            if(synchronousComputation_.get()) {
                computation_.setInput(std::move(prepareComputeInput()));
                computation_.threadMain();

                std::unique_ptr<O> synchronous_result = computation_.retrieveOutput();

                setProgressRange(tgt::vec2(0.0f, 1.0f));
                if(synchronous_result) {
                    processComputeOutput(std::move(*synchronous_result));
                    statusDisplay_.set("Finished");
                    setProgress(1.0f);
                } else {
                    LERROR("Computation execution failed");
                    statusDisplay_.set("Failed");
                    setProgress(0.0f);
                }
            } else {
                computation_.setInput(std::move(prepareComputeInput()));
                enableRunningState();
                statusDisplay_.set("Running");
                computationStartTime_ = Clock::now();
                computation_.run();
            }
            restartEnqueued_ = false;
        } catch(InvalidInputException& e) {
            std::string msg = e.what();
            std::string header = "Error starting computation:";
            switch(e.severity_) {
                case InvalidInputException::S_ERROR:

                    VoreenApplication::app()->showMessageBox(header, msg, true);
                    LERROR(msg);
                    break;
                case InvalidInputException::S_WARNING:

                    VoreenApplication::app()->showMessageBox(header, msg, false);
                    LWARNING(msg);
                    break;
                case InvalidInputException::S_IGNORE:
                    break;
            }
        }
    }
}

} // namespace voreen

#endif // VRN_ASYNCCOMPUTEPROCESSOR_H
