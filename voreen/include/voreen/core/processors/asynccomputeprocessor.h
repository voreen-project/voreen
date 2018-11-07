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

namespace voreen {

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
     * Note! The processor lock is assumed to be locked when calling this method
     * and will be locked, when returning from this method!
     */
    void stopComputeThread() throw();

    /**
     * Interrupts the computation and sets property states accordingly.
     * In contrast to stopComputeThread, this method should be noly called, if the
     * computation did not exit successfully.
     *
     * Note! The processor lock is assumed to be locked when calling this method
     * and will be locked, when returning from this method!
     */
    void interruptComputation() throw();

    /// Inherited from PortObserver
    void afterConnectionAdded(const Port* source, const Port* connectedPort);
    void beforeConnectionRemoved(const Port* source, const Port*);
    void dataWillChange(const Port* source);
    void dataHasChanged(const Port* source);

    /// Inherited from DataInvalidationObserver
    virtual void dataAboutToInvalidate(const DataInvalidationObservable* source);

private:

    // HACK: this is necessary due to a bug in msvc
    // https://connect.microsoft.com/VisualStudio/feedback/details/858357/
#if (_MSC_VER <= 1800)
    typedef std::chrono::system_clock Clock;
#else
    typedef std::chrono::steady_clock Clock;
#endif
    typedef std::chrono::time_point<Clock> TimePoint;

    const static std::chrono::milliseconds MINIMUM_PROGRESS_UPDATE_INTERVAL;

    /// Used to manage the visibility of properties during compute()
    bool isDisabledDuringComputation(Property* property) const;
    void savePropertyVisibilities();
    void enableRunningState();
    void disableRunningState();

    /**
     * Safely passes progress reports through to the processor in the main thread.
     */
    class ComputeProgressReporter : public ProgressReporter {
    public:
        ComputeProgressReporter(AsyncComputeProcessor<ComputeInput, ComputeOutput>& processor);
        virtual ~ComputeProgressReporter();

        virtual void setProgress(float progress);
        virtual float getProgress() const;

        virtual void setProgressRange(const tgt::vec2& progressRange);
        virtual tgt::vec2 getProgressRange() const;

        virtual void setProgressMessage(const std::string& message);
        virtual std::string getProgressMessage() const;

        virtual void update();
    private:
        AsyncComputeProcessor<ComputeInput, ComputeOutput>& processor_;
        TimePoint startTime_;
        TimePoint lastUpdate_;
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
    };

    BoolProperty continuousUpdate_;
    BoolProperty synchronousComputation_;
    ButtonProperty manualUpdateButton_;
    ButtonProperty stopUpdateButton_;
    ProgressProperty progressDisplay_;
    StringProperty statusDisplay_;

    /// Used to save/restore property writability before/after expensive computation
    std::map<Property*, bool> propertyWritablilityMap_;

    bool updateForced_;
    bool stopForced_;

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
    input_.reset(new I(std::move(input)));
}

template<class I, class O>
std::unique_ptr<O> AsyncComputeProcessor<I, O>::ComputeThread::retrieveOutput() {
    return std::move(output_);
}

template<class I, class O>
void AsyncComputeProcessor<I, O>::ComputeThread::clearOutput() {
    return output_.reset(nullptr);
}

template<class I, class O>
void AsyncComputeProcessor<I, O>::ComputeThread::threadMain() {
    ProcessorBackgroundThread<AsyncComputeProcessor<I, O>>::processor_->lockMutex();
    tgtAssert(input_, "ComputeThread started without input!");
    output_.reset();

    std::unique_ptr<ProgressReporter> progress;
    if(ProcessorBackgroundThread<AsyncComputeProcessor<I,O>>::processor_->synchronousComputation_.get()) {
        progress.reset(new NoopProgressReporter());
    } else {
        progress.reset(new ComputeProgressReporter(*ProcessorBackgroundThread<AsyncComputeProcessor<I, O>>::processor_));
    }
    I* input = input_.release();
    ProcessorBackgroundThread<AsyncComputeProcessor<I, O>>::processor_->unlockMutex();

    O* output = new O(std::move(ProcessorBackgroundThread<AsyncComputeProcessor<I,O>>::processor_->compute(std::move(*input), *progress)));

    ProcessorBackgroundThread<AsyncComputeProcessor<I, O>>::processor_->lockMutex();
    output_.reset(output);
    ProcessorBackgroundThread<AsyncComputeProcessor<I, O>>::processor_->unlockMutex();
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
void AsyncComputeProcessor<I,O>::ComputeProgressReporter::setProgress(float progress) {
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

template<class I, class O>
float AsyncComputeProcessor<I,O>::ComputeProgressReporter::getProgress() const {
    tgtAssert(processor_.getProgress() == processor_.progressDisplay_.getProgress(), "Progress mismatch");
    return processor_.getProgress();
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::ComputeProgressReporter::setProgressRange(const tgt::vec2& progressRange) {
    processor_.setProgressRange(progressRange);
}

template<class I, class O>
tgt::vec2 AsyncComputeProcessor<I,O>::ComputeProgressReporter::getProgressRange() const {
    tgtAssert(processor_.getProgressRange() == processor_.progressDisplay_.getProgressRange(), "ProgressRange mismatch");
    return processor_.getProgressRange();
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::ComputeProgressReporter::setProgressMessage(const std::string& message) {
    processor_.setProgressMessage(message);
}

template<class I, class O>
std::string AsyncComputeProcessor<I,O>::ComputeProgressReporter::getProgressMessage() const {
    tgtAssert(processor_.getProgressMessage() == processor_.progressDisplay_.getProgressMessage(), "ProgressMessage mismatch");
    return processor_.getProgressMessage();
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::ComputeProgressReporter::update() {
    processor_.update();
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
    , continuousUpdate_("continuousUpdate_", "Continuous Update", false, Processor::INVALID_RESULT, Property::LOD_ADVANCED)
    , synchronousComputation_("synchronousComputation", "Wait for Result", false, Processor::INVALID_RESULT, Property::LOD_ADVANCED)
    , manualUpdateButton_("manualUpdateButton_", "Start", Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , stopUpdateButton_("stopUpdateButton", "Stop", Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , progressDisplay_("progressDisplay", "Progress")
    , statusDisplay_("statusDisplay", "Status", "Stopped", Processor::VALID)
    , updateForced_(false)
    , stopForced_(false)
    , computation_(*this)
{
    addProperty(continuousUpdate_);
        continuousUpdate_.setGroupID("ac_processing");
    addProperty(synchronousComputation_);
        synchronousComputation_.setGroupID("ac_processing");

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
        statusDisplay_.setReadOnly(true);
    setPropertyGroupGuiName("ac_processing", "Processing");

    addProgressBar(&progressDisplay_);
}

template<class I, class O>
AsyncComputeProcessor<I,O>::~AsyncComputeProcessor() {
}

template<class I, class O>
bool AsyncComputeProcessor<I,O>::isDisabledDuringComputation(Property* property) const {
    return property != &manualUpdateButton_
        && property != &stopUpdateButton_
        && property != &progressDisplay_
        && property != &statusDisplay_;
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::savePropertyVisibilities() {
    for(Property* property : getProperties()) {
        if(isDisabledDuringComputation(property)) {
            propertyWritablilityMap_[property] = property->isReadOnlyFlagSet();
        }
    }
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::enableRunningState() {
    savePropertyVisibilities();
    for(Property* property : getProperties()) {
        if(isDisabledDuringComputation(property)) {
            property->setReadOnlyFlag(true);
        }
    }

    manualUpdateButton_.setVisibleFlag(false);
    stopUpdateButton_.setVisibleFlag(true);
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::disableRunningState() {
    // Restore property visibilities
    for(Property* property : getProperties()) {
        if(isDisabledDuringComputation(property)) {
            tgtAssert(propertyWritablilityMap_.find(property) != propertyWritablilityMap_.end(), "Invalid property. Did you add properties during computation?");
            property->setReadOnlyFlag(propertyWritablilityMap_[property]);
        }
    }

    manualUpdateButton_.setVisibleFlag(true);
    stopUpdateButton_.setVisibleFlag(false);
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

template<class I, class O>
void AsyncComputeProcessor<I, O>::dataWillChange(const Port* source) {
    lockMutex();
    interruptComputation();
    // Clear the result, if we have one already, as it will not be valid (i.e., corresponding to the input) afterwards.
    computation_.clearOutput();
    // The port will no longer store it's currently contained data.
    // The data itself might continue existing, therefore we need to remove the observer.
    if(source->isDataInvalidationObservable() && source->hasData()) {
        source->removeDataInvalidationObserver(static_cast<DataInvalidationObserver*>(this));
    }
    unlockMutex();
}

template<class I, class O>
void AsyncComputeProcessor<I, O>::dataHasChanged(const Port* source) {
    lockMutex();
    interruptComputation();
    // The port stores new data.
    // The data might change later, so we register ourselves as an observer if this is possible
    if(source->isDataInvalidationObservable() && source->hasData()) {
        source->addDataInvalidationObserver(static_cast<DataInvalidationObserver*>(this));
    }
    unlockMutex();
}

template<class I, class O>
void AsyncComputeProcessor<I, O>::dataAboutToInvalidate(const DataInvalidationObservable* source) {
    lockMutex();
    interruptComputation();
    unlockMutex();
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::initialize() {
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
    savePropertyVisibilities();
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::deinitialize() {
    // We obviously do not want the thread to keep running if the processor is deinitialized,
    // we will never retrieve the result, anyway.
    lockMutex();
    stopComputeThread();

    for(auto port : getCriticalPorts()) {
        // Do not observe this port anymore
        static_cast<const Observable<PortObserver>&>(port.get()).removeObserver(static_cast<PortObserver*>(this));

        if (port.get().hasData() && port.get().isDataInvalidationObservable()) {
            port.get().removeDataInvalidationObserver(static_cast<DataInvalidationObserver*>(this));
        }
    }

    unlockMutex();
    Processor::deinitialize();
}

template<class I, class O>
void AsyncComputeProcessor<I,O>::stopComputeThread() throw() {
    // We need to unlock the processor mutex, before interrupting and joining
    // the computation thread, in order to avoid deadlocks (see backgroundprocess.h)
    unlockMutex();
    computation_.interruptAndJoin();
    lockMutex();

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

    bool abortNecessary = (getInvalidationLevel() >= INVALID_RESULT || updateForced_ || stopForced_) && computation_.isRunning();
    stopForced_ = false;

    if(abortNecessary) {
        interruptComputation();
    }

    std::unique_ptr<O> result = computation_.retrieveOutput();

    bool restartNecessary = (!result && getInvalidationLevel() >= INVALID_RESULT && continuousUpdate_.get()) && !stopForced_ || updateForced_;
    updateForced_ = false;

    //If we need to restart anyways, we don't care about the result
    if(result && !restartNecessary) {
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
            computation_.setInput(std::move(prepareComputeInput()));
            enableRunningState();
            statusDisplay_.set("Running");
            computationStartTime_ = Clock::now();
            computation_.run();

            if(synchronousComputation_.get()) {
                unlockMutex();
                computation_.join();
                lockMutex();

                std::unique_ptr<O> synchronous_result = computation_.retrieveOutput();
                disableRunningState();
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
            }
        } catch(InvalidInputException& e) {
            std::string msg = e.what();
            std::string header = "Error starting computation:";
            switch(e.severity_) {
                case InvalidInputException::S_ERROR:

                    VoreenApplication::app()->showMessageBox(header, msg, true);
                    LERROR(msg);
                    break;
                case InvalidInputException::S_WARNING:

                    VoreenApplication::app()->showMessageBox(header, msg, true);
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
