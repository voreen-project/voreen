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

#include "voreen/core/network/networkevaluator.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/network/processornetwork.h"
#include "voreen/core/interaction/idmanager.h"
#include "voreen/core/network/networkgraph.h"
#include "voreen/core/utils/exception.h"

#include "modules/core/processors/output/canvasrenderer.h" //< core module is always available

#include "tgt/textureunit.h"
#include "tgt/framebufferobject.h"
#include "tgt/glcontextmanager.h"
#include "tgt/immediatemode/immediatemode.h"

#include <vector>

namespace voreen {

const std::string NetworkEvaluator::loggerCat_("voreen.NetworkEvaluator");

NetworkEvaluator::NetworkEvaluator(bool glMode)
    : network_(0)
    , glMode_(glMode)
    //, reuseRenderTargets_(false)
    , renderingOrder_()
    , networkChanged_(false)
    , locked_(false)
    , processPending_(false)
{
#ifdef VRN_DEBUG
    if (glMode_)
        addObserver(new CheckOpenGLStateObserver());
#endif
}

NetworkEvaluator::~NetworkEvaluator() {
#ifdef VRN_DEBUG
    // delete CheckOpenGLStateObserver
    std::vector<NetworkEvaluatorObserver*> observers = getObservers();
    for (size_t i=0; i<observers.size(); ++i) {
        if (dynamic_cast<CheckOpenGLStateObserver*>(observers[i]))
            delete observers[i];
    }
#endif
}

const std::set<Processor*> NetworkEvaluator::getEndProcessors() const {
    tgtAssert(network_, "No network");
    std::set<Processor*> endProcessors;
    for (size_t i = 0; i < network_->getProcessors().size(); ++i) {
        if (network_->getProcessors()[i]->isEndProcessor() == true)
            endProcessors.insert(network_->getProcessors()[i]);
    }
    return endProcessors;
}

bool NetworkEvaluator::initializeNetwork()  {

    if (isLocked()) {
        LDEBUG("initializeNetwork() called on locked evaluator.");
        return false;
    }

    if (!network_) {
        LWARNING("initializeNetwork() called on evaluator without network.");
        return false;
    }

    // prevent parallel execution in multithreaded/event dispatching environments
    lock();

    // Ensure the correct context is active before and after initialization.
    tgt::GLConditionalContextStateGuard guard(glMode_);

    // notify observers
    const std::vector<NetworkEvaluatorObserver*> observers = getObservers();
    for (size_t i = 0; i < observers.size(); ++i)
        observers[i]->beforeNetworkInitialize();
    if (glMode_) LGL_ERROR;

    bool failed = false;
    for (size_t i = 0; i < network_->getProcessors().size(); ++i) {
        Processor* processor = network_->getProcessors()[i];
        if (!processor->isInitialized()) {
            try {
                tgt::GLConditionalContextStateGuard guard(glMode_);
                processor->initialize();
                //update  processor state
                processor->processorState_ = Processor::PROCESSOR_STATE_NOT_READY;
                processor->notifyStateChanged();
                processor->invalidate();
            }
            catch (const tgt::Exception& e) {
                LERROR("Failed to initialize processor '" << processor->getID()
                        << "' (" << processor->getClassName() << "): ");
                LERROR(" - " << e.what());

                // deinitialize processor, in order to make sure that all resources are freed
                LINFO("Deinitializing '" << processor->getID()
                    << "' (" << processor->getClassName() << ") ...");
                //set state as initialized
                processor->processorState_ = Processor::PROCESSOR_STATE_NOT_READY;
                // HACK: deinitializing unintiialized processors can lead to crashes
                //processor->deinitialize();
                processor->processorState_ = Processor::PROCESSOR_STATE_NOT_INITIALIZED;

                // don't break, try to initialize the other processors even if one failed
                failed = true;
                continue;
            }

            // Check for errors on each successfully initialized processor.
            if (glMode_) LGL_ERROR;
        }
    }

    assignRenderTargets();

    unlock();

    // notify observers
    for (size_t i = 0; i < observers.size(); ++i)
        observers[i]->afterNetworkInitialize();
    if (glMode_) LGL_ERROR;

    return !failed;
}

bool NetworkEvaluator::deinitializeNetwork() {
    if (isLocked()) {
        LWARNING("deinitializeNetwork() called on locked evaluator.");
        return false;
    }

    if (!network_) {
        LWARNING("deinitializeNetwork() called on evaluator without network.");
        return false;
    }

    // prevent parallel execution in multithreaded/event dispatching environments
    lock();

    // Ensure the correct context is active before and after deinitialization.
    tgt::GLConditionalContextStateGuard guard(glMode_);

    // notify observers
    const std::vector<NetworkEvaluatorObserver*> observers = getObservers();
    for (size_t i = 0; i < observers.size(); ++i)
        observers[i]->beforeNetworkDeinitialize();
    if (glMode_) LGL_ERROR;
    

    bool failed = false;
    for (size_t i = 0; i < network_->getProcessors().size(); ++i) {
        Processor* processor = network_->getProcessors()[i];
        if (processor->isInitialized()) {
            try {
                tgt::GLConditionalContextStateGuard guard(glMode_);
                processor->deinitialize();
                processor->processorState_ = Processor::PROCESSOR_STATE_NOT_INITIALIZED;
            }
            catch (const tgt::Exception& e) {
                LERROR("Failed to deinitialize processor '" << processor->getID()
                    << "' (" << processor->getClassName() << "): ");
                LERROR(" - " << e.what());
                // don't break, try to deinitialize the other processors even if one failed
                failed = true;
                continue;
            }

            // Check for errors on each successfully deinitialized processor.
            if (glMode_) LGL_ERROR;
        }
    }

    unlock();

    // notify observers
    for (size_t i = 0; i < observers.size(); ++i)
        observers[i]->afterNetworkDeinitialize();
    if (glMode_) LGL_ERROR;

    return !failed;
}

void NetworkEvaluator::onNetworkChange() {
    if (isLocked()) {
        LWARNING("onNetworkChange() called on locked evaluator.");
        return;
    }

    networkChanged_ = false;

    {
        PROFILING_BLOCK("onNetworkChange");

        renderingOrder_.clear();
        loopPortMap_.clear();

        // nothing more to do, if no network is present
        if (!network_)
            return;

        // define rendering order from network
        {
            PROFILING_BLOCK("defineRenderingOrder");
            defineRenderingOrder();
        }

        // initialize processors
        {
            PROFILING_BLOCK("initialize");
            initializeNetwork();
            if (glMode_) LGL_ERROR;
        }
    }
#ifdef VRN_PRINT_PROFILING
    performanceRecord_.getLastSample()->print();
#endif

}

void NetworkEvaluator::process() {

    if (!network_)
        return;

    if (isLocked()) {
        LDEBUG("process() called on locked evaluator. Scheduling.");
        processPending_ = true;
        return;
    }

    // VoreenPainter calls the process function having the related
    // GLCanvas active. Since process needs the main context to be active,
    // we need a guard here.
    tgt::GLConditionalContextStateGuard guard(glMode_);

    // re-analyze and initialize network, if its topology has changed since last process() call
    if (networkChanged_ || checkForInvalidPorts()) {
        onNetworkChange();
    }

    // prevent parallel execution in multithreaded/event dispatching environments
    lock();

    LDEBUG("Starting network evaluation");

    if (renderingOrder_.empty()) {
        LDEBUG("process(): rendering order is not defined!");
    }

    // reset loop ports' iteration counters and progressbars
    for (size_t i = 0; i < renderingOrder_.size(); ++i) {
        Processor* processor = renderingOrder_[i];
        tgtAssert(loopPortMap_.find(processor) != loopPortMap_.end(), "Processor missing in loop port map");
        for (size_t j = 0; j < loopPortMap_[processor].size(); ++j)
            loopPortMap_[processor][j]->setLoopIteration(-1);

        // this may lead to continuous network invalidations/evaluations (infinite loop)
        //if (!processor->isValid())
        //    processor->setProgress(0.f);
    }

    // notify observers
    const std::vector<NetworkEvaluatorObserver*> observers = getObservers();
    for (size_t j = 0; j < observers.size(); ++j)
        observers[j]->beforeNetworkProcess();
    if (glMode_) LGL_ERROR;

    // Iterate over processing in rendering order
    for (size_t i = 0; i < renderingOrder_.size(); ++i) {
        Processor* const currentProcessor = renderingOrder_[i];

        // all processors should have been initialized at this point
        if (!currentProcessor->isInitialized()) {
            LWARNING("process(): Skipping uninitialized processor '" << currentProcessor->getID()
                     << "' (" << currentProcessor->getClassName() << ")");
            continue;
        }

        // trigger property adjustment on on input change
        for (size_t i = 0; i < currentProcessor->inports_.size(); i++) {
            if (currentProcessor->inports_[i]->hasChanged()) {
                currentProcessor->adjustPropertiesToInput();
                break;
            }
        }

        bool needsProcessing = true;
        if (currentProcessor->isValid())
            needsProcessing = false;

        // run the processor, if it needs processing and is ready
        if (needsProcessing) {
            if (currentProcessor->isReady()) {

                // increase iteration counters
                for (size_t j=0; j<loopPortMap_[currentProcessor].size(); ++j) {
                    Port* port = loopPortMap_[currentProcessor][j];
                    // note: modulo is required for nested loops
                    port->setLoopIteration((port->getLoopIteration()+1) % port->getNumLoopIterations());
                }

                // notify observers
                for (size_t j=0; j < observers.size(); ++j)
                    observers[j]->beforeProcess(currentProcessor);
                if (glMode_) LGL_ERROR;

                try {
                    currentProcessor->performanceRecord_.setName(currentProcessor->getID());
                    currentProcessor->lockMutex();

                    {
                        ProfilingBlock block("beforeprocess", currentProcessor->performanceRecord_);
                        currentProcessor->beforeProcess();
                    }
                    if (glMode_) LGL_ERROR;

                    if(currentProcessor->getInvalidationLevel() >= Processor::INVALID_PORTS) {
                        currentProcessor->unlockMutex();
                        unlock();

                        if (glMode_) LGL_ERROR;

                        onNetworkChange();
                        currentProcessor->invalidate();
                        VoreenApplication::app()->scheduleNetworkProcessing();
                        return;
                    }

#ifdef VRN_PRINT_PROFILING
                    currentProcessor->performanceRecord_.getLastSample()->print(0, currentProcessor->getID()+".");
#endif
                    if (!currentProcessor->isValid())
                    {
                        ProfilingBlock block("process", currentProcessor->performanceRecord_);
                        tgt::GLConditionalContextStateGuard guard(glMode_);
                        currentProcessor->process();
                    }
#ifdef VRN_PRINT_PROFILING
                    currentProcessor->performanceRecord_.getLastSample()->print(0, currentProcessor->getID()+".");
#endif
                    if (glMode_) LGL_ERROR;
                    {
                        ProfilingBlock block("afterprocess", currentProcessor->performanceRecord_);
                        currentProcessor->afterProcess();
                    }
#ifdef VRN_PRINT_PROFILING
                    currentProcessor->performanceRecord_.getLastSample()->print(0, currentProcessor->getID()+".");
#endif
                    if (glMode_) LGL_ERROR;

                    currentProcessor->unlockMutex();
                }
                catch (VoreenException& e) {
                    currentProcessor->unlockMutex();
                    LERROR("process(): VoreenException from "
                            << currentProcessor->getClassName()
                            << " (" << currentProcessor->getID() << "): " << e.what());
                }
                catch (std::exception& e) {
                    currentProcessor->unlockMutex();
                    LERROR("process(): Exception from "
                            << currentProcessor->getClassName()
                            << " (" << currentProcessor->getID() << "): " << e.what());
                }

                // notify observers
                for (size_t j = 0; j < observers.size(); ++j)
                    observers[j]->afterProcess(currentProcessor);
                if (glMode_) LGL_ERROR;

                // break loop if network topology has changed (due to changes in loop port configurations)
                if (checkForInvalidPorts()) {
                    unlock();

                    // notify observers
                    for (size_t j = 0; j < observers.size(); ++j)
                        observers[j]->afterNetworkProcess();
                    if (glMode_) LGL_ERROR;

                    onNetworkChange();
                    return;
                }
            }
            else {
                // Processor isn't ready, clear outports:
                currentProcessor->clearOutports();
                // set Processor valid, as it should not be processed while not ready
                // TODO: rename "invalid" to "needsProcessing"
                //currentProcessor->setValid();
            }
        } // needsProcessing

        currentProcessor->firstProcessAfterDeserialization_ = false;

    }   // for (rendering order)

    if (glMode_) LGL_ERROR;

    // notify observers
    for (size_t j = 0; j < observers.size(); ++j)
        observers[j]->afterNetworkProcess();
    if (glMode_) LGL_ERROR;

    LDEBUG("Finished network evaluation");

    unlock();

    if (processPending_ && glMode_) {
        // make sure that canvases are repainted, if their update has been blocked by the locked evaluator
        processPending_ = false;
        updateCanvases();
        LGL_ERROR;
    }

    // TODO: find better solution, very ugly code!
    for (std::vector<Processor*>::const_iterator iter = getProcessorNetwork()->getProcessors().begin(); iter != getProcessorNetwork()->getProcessors().end(); ++iter) {
        if (((*iter)->isReady() && !(*iter)->isValid()) && (((*iter)->getClassName().compare("Canvas") != 0) && ((*iter)->getClassName().compare("StereoCanvas") != 0))){
            tgtAssert(VoreenApplication::app(), "VoreenApplication not instantiated");
            VoreenApplication::app()->scheduleNetworkProcessing();
            break;
        }
    }

    if (glMode_) LGL_ERROR;
}

void NetworkEvaluator::setProcessorNetwork(ProcessorNetwork* network, bool deinitializeCurrent) {
    if (network_) {
        stopObservation(network_);
        if (deinitializeCurrent)
            deinitializeNetwork();
    }

    // assign new network
    ProcessorNetwork* previousNetwork = network_;
    network_ = network;

    // notify observers
    const std::vector<NetworkEvaluatorObserver*> observers = getObservers();
    for (size_t i=0; i<observers.size(); i++)
        observers[i]->networkAssigned(network_, previousNetwork);

    if (network_)
        network_->addObserver(this);

    // If evaluation has been locked by the user, he might not want
    // to trigger an automatic invalidation after loading a workspace.
    // Therefore, schedule evaluation.
    if(isLocked())
        networkChanged_ = true;
    else
        onNetworkChange();
}

void NetworkEvaluator::invalidateProcessors(Processor::InvalidationLevel level) {

    if (isLocked()) {
        LWARNING("invalidateProcessors() called on locked evaluator.");
        return;
    }

    if (!network_)
        return;

    for (size_t i=0; i<network_->getProcessors().size(); ++i) {
        network_->getProcessors()[i]->invalidate(level);
    }
}

void NetworkEvaluator::forceUpdate() {

    if (isLocked()) {
        LWARNING("forceUpdate() called on locked evaluator.");
        return;
    }

    invalidateProcessors();
    process();
}

// protected methods
//

void NetworkEvaluator::defineRenderingOrder() {

    tgtAssert(network_, "No processor network");

    if (isLocked()) {
        LWARNING("defineRenderingOrder() called on locked evaluator.");
        return;
    }

    // construct graph: ignore loop ports for graph construction (second parameter), but unroll them afterwards (third parameter)
    NetworkGraph netGraph(network_->getProcessors(), LoopPortTypeCheck(true), LoopPortTypeCheck());

    // generate subset of processors that are predecessors of the end processors
    std::set<Processor*> endProcessors = getEndProcessors();
    std::set<Processor*> predecessors = netGraph.getPredecessors(endProcessors);

    // sort generated subset topologically
    renderingOrder_ = netGraph.sortTopologically(predecessors);

    // construct mapping from processors to their loop inports
    // (used for incrementing loop ports' iteration counter)
    loopPortMap_.clear();
    for (size_t i=0; i<renderingOrder_.size(); ++i) {
        Processor* processor = renderingOrder_[i];
        if (loopPortMap_.find(processor) == loopPortMap_.end())
            loopPortMap_[processor] = std::vector<Port*>();
        else
            continue;
        for (size_t j=0; j<processor->getInports().size(); ++j) {
            Port* port = processor->getInports()[j];
            if (port->isLoopPort()) {
                if (port->isConnected())
                    loopPortMap_[processor].push_back(port);
                else
                    port->setLoopIteration(0);
            }
        }
    }

    // reduce processors' invalidation level in order to prevent
    // a continuous re-analysis of the network
    for (size_t i=0; i<network_->getProcessors().size(); ++i) {
        if (network_->getProcessors()[i]->getInvalidationLevel() >= Processor::INVALID_PORTS)
            network_->getProcessors()[i]->invalidationLevel_ = Processor::INVALID_PROGRAM;
    }
}

void NetworkEvaluator::assignRenderTargets() {

    // ports to process (connected outports)
    std::vector<RenderPort*> PORTS;

    // direct successors of each port in the network graph (non-transitive)
    std::map< RenderPort*, std::vector<RenderPort*> > SUCC_GRAPH;

    // successors of each port in the rendering order (transitive)
    std::map< RenderPort*, std::vector<RenderPort*> > SUCC_ORDER;

    // predecessors of each port in the rendering order (transitive)
    std::map< RenderPort*, std::vector<RenderPort*> > PRED_ORDER;

    // collect ports to process
    std::vector<Processor*> processors = network_->getProcessors();
    for (size_t p = 0; p < processors.size(); ++p) {
        if (RenderProcessor* rp = dynamic_cast<RenderProcessor*>(processors[p])) {

            //if (!rp->isInitialized())
            //    continue;

            // outports
            std::vector<Port*> outports = rp->getOutports();
            for (size_t o = 0; o < outports.size(); ++o) {
                if (RenderPort* rport = dynamic_cast<RenderPort*>(outports[o])) {
                    if (rport->isConnected())
                        PORTS.push_back(rport);
                }
            }

            // private ports
            const std::vector<RenderPort*>& privatePorts = rp->getPrivateRenderPorts();
            for (size_t o = 0; o < privatePorts.size(); ++o) {
                if (privatePorts[o]) {
                    PORTS.push_back(privatePorts[o]);
                }
            }
        }
    }

    // collect graph successors
    for (size_t portID = 0; portID < PORTS.size(); ++portID) {
        RenderPort* port = PORTS.at(portID);
        std::vector<RenderPort*> successors;

        // iterate over all connected processors and add their render outports to the successors list
        const std::vector<const Port*> connectedPorts = port->getConnected();
        for (size_t connPortID = 0; connPortID < connectedPorts.size(); ++connPortID) {
            tgtAssert(dynamic_cast<const RenderPort*>(connectedPorts.at(connPortID)), "RenderPort connected to Non-RenderPort");
            tgtAssert(connectedPorts.at(connPortID)->getProcessor(), "RenderPort without owner");
            std::vector<Port*> succOutports = connectedPorts.at(connPortID)->getProcessor()->getOutports();
            for (size_t i = 0; i < succOutports.size(); ++i) {
                // add outport to successor list, if it is a RenderPort and element of PORTS
                if (RenderPort* succPort = dynamic_cast<RenderPort*>(succOutports.at(i))) {
                    if (std::find(PORTS.begin(), PORTS.end(), succPort) != PORTS.end())
                        successors.push_back(succPort);
                }
            }
        }

        // add successor list to map
        SUCC_GRAPH.insert(std::pair<RenderPort*, std::vector<RenderPort*> >(port, successors));
    }

}

bool NetworkEvaluator::checkForInvalidPorts() {
    tgtAssert(network_, "no network");

    for (size_t i=0; i<network_->getProcessors().size(); ++i) {
        if (network_->getProcessors()[i]->getInvalidationLevel() >= Processor::INVALID_PORTS) {
            return true;
        }
    }
    return false;
}

const ProcessorNetwork* NetworkEvaluator::getProcessorNetwork() const {
    return network_;
}

std::vector<RenderPort*> NetworkEvaluator::collectRenderPorts() const {

    std::vector<RenderPort*> renderPorts;
    if (isLocked()) {
        LWARNING("collectRenderPorts() called on locked evaluator.");
        return renderPorts;
    }

    if (!network_) {
        //LWARNING("collectRenderPorts() called on evaluator without network.");
        return renderPorts;
    }

    std::vector<Processor*> processors = network_->getProcessors();
    for (size_t p = 0; p < processors.size(); ++p) {
        RenderProcessor* rp = dynamic_cast<RenderProcessor*>(processors[p]);

        if (!rp)
            continue;

        if (!rp->isInitialized())
            continue;

        std::vector<Port*> outports = rp->getOutports();
        for (size_t o = 0; o < outports.size(); ++o) {
            Port* const outport = outports[o];

            RenderPort* rport = dynamic_cast<RenderPort*>(outport);
            if (rport) {
                if (rport->hasRenderTarget())
                    renderPorts.push_back(rport);
            }
        }
        //add private ports:
        const std::vector<RenderPort*>& privatePorts = rp->getPrivateRenderPorts();
        for (size_t o = 0; o < privatePorts.size(); ++o) {
            RenderPort* pport = privatePorts[o];

            if (pport) {
                if (pport->hasRenderTarget())
                    renderPorts.push_back(pport);
            }
        }
    }

    return renderPorts;
}

std::vector<const PerformanceRecord*> NetworkEvaluator::collectPerformanceRecords() const {

    std::vector<const PerformanceRecord*> performanceRecords;
    if (isLocked()) {
        LWARNING("collectPerformanceRecords() called on locked evaluator.");
        return performanceRecords;
    }

    if (!network_) {
        return performanceRecords;
    }

    std::vector<Processor*> processors = network_->getProcessors();
    for (size_t p = 0; p < processors.size(); ++p) {
        Processor* pro = processors[p];

        if (!pro->isInitialized())
            continue;

        performanceRecords.push_back(pro->getPerformanceRecord());
    }

    return performanceRecords;
}

void NetworkEvaluator::clearPerformanceRecords() {

    if (isLocked()) {
        LWARNING("clearPerformanceRecords() called on locked evaluator.");
        return;
    }

    if (!network_) {
        return;
    }

    std::vector<Processor*> processors = network_->getProcessors();
    for (size_t p = 0; p < processors.size(); ++p) {
        Processor* pro = processors[p];

        if (!pro->isInitialized())
            continue;

        pro->resetPerformanceRecord();
    }
}

void NetworkEvaluator::updateCanvases() {

    if (!glMode_) {
        LWARNING("updateCanvases() called in non OpenGL mode!");
        return;
    }

    if (isLocked()) {
        LWARNING("updateCanvases() called on locked evaluator!");
        return;
    }

    if (!network_) {
        LWARNING("updateCanvases() called on evaluator without network!");
        return;
    }

    for (size_t i = 0; i < network_->getProcessors().size(); ++i) {
        CanvasRenderer* canvasRenderer = dynamic_cast<CanvasRenderer*>(network_->getProcessors()[i]);
        if (canvasRenderer)
            canvasRenderer->invalidate();
    }
}

void NetworkEvaluator::lock() {
    locked_ = true;
}

void NetworkEvaluator::unlock() {
    locked_ = false;
}

bool NetworkEvaluator::isLocked() const {
    return locked_;
}

void NetworkEvaluator::networkChanged() {
    networkChanged_ = true;
}

void NetworkEvaluator::processorRenamed(const Processor* /*processor*/, const std::string& /*prevName*/) {
    // ignore
}

void NetworkEvaluator::propertyLinkAdded(const PropertyLink* /*link*/) {
    // ignore
}

void NetworkEvaluator::propertyLinkRemoved(const PropertyLink* /*link*/) {
    // ignore
}

void NetworkEvaluator::processorRemoved(const Processor* processor) {

    if (processor->isInitialized()) {
        try {
            tgt::GLConditionalContextStateGuard guard(glMode_); // Restore correct context afterwards!
            const_cast<Processor*>(processor)->deinitialize();
            const_cast<Processor*>(processor)->processorState_ = Processor::PROCESSOR_STATE_NOT_INITIALIZED;
         }
         catch (const VoreenException& e) {
            LERROR("Failed to deinitialize '" << processor->getID()
                    << "' (" << processor->getClassName() << "): ");
            LERROR(" - " << e.what());
         }
    }

    networkChanged_ = true;
}


// ----------------------------------------------------------------------------

// processWrapper stuff pending...

namespace {

bool checkGL(GLenum pname, bool value) {
    GLboolean b;
    glGetBooleanv(pname, &b);
    return (static_cast<bool>(b) == value);
}

bool checkGL(GLenum pname, GLint value) {
    GLint i;
    glGetIntegerv(pname, &i);
    return (i == value);
}

bool checkGL(GLenum pname, GLfloat value) {
    GLfloat f;
    glGetFloatv(pname, &f);
    return (f == value);
}

bool checkGL(GLenum pname, const tgt::vec4 value) {
    tgt::vec4 v;
    glGetFloatv(pname, reinterpret_cast<float*>(&v.elem));
    return (v == value);
}

} // namespace

void NetworkEvaluator::CheckOpenGLStateObserver::afterProcess(Processor* p) {
    checkState(p);
}

void NetworkEvaluator::CheckOpenGLStateObserver::beforeNetworkProcess() {
    checkState();
}

void NetworkEvaluator::CheckOpenGLStateObserver::checkState(Processor* p) {
    //check for main context
    if (!GLContextMgr.isMainContextActive()) {
        warn(p, "Main GLContext was not current context");
        GLContextMgr.activateMainContext();
    }

    // check OpenGL status
    if (!checkGL(GL_BLEND, false)) {
        glDisable(GL_BLEND);
        warn(p, "GL_BLEND was enabled");
    }

    if (!checkGL(GL_BLEND_SRC, GL_ONE) || !checkGL(GL_BLEND_DST, GL_ZERO)) {
        glBlendFunc(GL_ONE, GL_ZERO);
        warn(p, "Modified BlendFunc");
    }

    if (!checkGL(GL_BLEND_EQUATION_RGB, GL_FUNC_ADD) || !checkGL(GL_BLEND_EQUATION_ALPHA, GL_FUNC_ADD)) {
        glBlendEquation(GL_FUNC_ADD);
        warn(p, "Modified Blend Equation");
    }

    if (!checkGL(GL_DEPTH_TEST, true)) {
        glEnable(GL_DEPTH_TEST);
        warn(p, "GL_DEPTH_TEST was not enabled");
    }

    if (!checkGL(GL_DEPTH_FUNC, GL_LESS)) {
        glDepthFunc(GL_LESS);
        warn(p, "glDepthFunc was not set to GL_LESS");
    }

    if (!checkGL(GL_CULL_FACE, false)) {
        glDisable(GL_CULL_FACE);
        warn(p, "GL_CULL_FACE was enabled");
    }

    if (!checkGL(GL_COLOR_CLEAR_VALUE, tgt::vec4(0.f))) {
        glClearColor(0.f, 0.f, 0.f, 0.f);
        warn(p, "glClearColor() was not set to all zeroes");
    }

    if (!checkGL(GL_DEPTH_CLEAR_VALUE, 1.f)) {
        glClearDepth(1.0);
        warn(p, "glClearDepth() was not set to 1.0");
    }
#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
    if (!checkGL(GL_LIGHTING, false)) {
        glDisable(GL_LIGHTING);
        warn(p, "GL_LIGHTING was enabled");
    }

    if (!checkGL(GL_CURRENT_COLOR, tgt::vec4(1.f))) {
        glColor4f(1.f, 1.f, 1.f, 1.f);
        warn(p, "glColor() was not set to all ones");
    }

    if (!checkGL(GL_LINE_WIDTH, 1.f)) {
        glLineWidth(1.f);
        warn(p, "glLineWidth() was not set to 1.0");
    }

    if (!checkGL(GL_POINT_SIZE, 1.f)) {
        glPointSize(1.f);
        warn(p, "glPointSize() was not set to 1.0");
    }

    if (!checkGL(GL_MATRIX_MODE, GL_MODELVIEW)) {
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        warn(p, "glMatrixMode was not set to GL_MODELVIEW");
    }

    if (!checkGL(GL_ACTIVE_TEXTURE, GL_TEXTURE0)) {
        glActiveTexture(GL_TEXTURE0);
        warn(p, "glActiveTexture was not set to GL_TEXTURE0");
    }

    if (!checkGL(GL_TEXTURE_1D, false)) {
        glDisable(GL_TEXTURE_1D);
        warn(p, "GL_TEXTURE_1D was enabled");
    }

    if (!checkGL(GL_TEXTURE_2D, false)) {
        glDisable(GL_TEXTURE_2D);
        warn(p, "GL_TEXTURE_2D was enabled");
    }

    if (!checkGL(GL_TEXTURE_3D, false)) {
        glDisable(GL_TEXTURE_3D);
        warn(p, "GL_TEXTURE_3D was enabled");
    }
#endif

    if (tgt::Shader::getCurrentProgram() != 0) {
        tgt::Shader::deactivate();
        warn(p, "A shader was active");
    }

    if (tgt::FramebufferObject::getActiveObject() != 0) {
        tgt::FramebufferObject::deactivate();
        warn(p, "A framebuffer object was active (RenderPort::deactivateTarget() missing?)");
    }

    if (!checkGL(GL_CULL_FACE_MODE, GL_BACK)) {
        glCullFace(GL_BACK);
        warn(p, "glCullFace was not set to GL_BACK");
    }

    if (!tgt::TextureUnit::unused()) {
        // Warn only if units that were initialized non-statically are still active (which probably means
        // they were generated with the new operator and not deleted)
        if (tgt::TextureUnit::numLocalActive() > 0)
            warn(p, "tgt::TextureUnit still had active texture units");
        tgt::TextureUnit::cleanup();
    }

    if (MatStack.getMatrixMode() != tgt::MatrixStack::MODELVIEW){
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        warn(p, "MatrixMode is not MODELVIEW.");
    }

    if(MatStack.getStackSize(tgt::MatrixStack::MODELVIEW) != 1){
        while (MatStack.getStackSize(tgt::MatrixStack::MODELVIEW) != 1)
            MatStack.popMatrix();
        warn(p, "ModelView-matrix stack is not empty.");
    }

    if(MatStack.getStackSize(tgt::MatrixStack::PROJECTION) != 1){
        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        while (MatStack.getStackSize(tgt::MatrixStack::PROJECTION) != 1)
            MatStack.popMatrix();
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        warn(p, "Projection-matrix stack is not empty.");
    }

    if(MatStack.getStackSize(tgt::MatrixStack::TEXTURE) != 1){
        MatStack.matrixMode(tgt::MatrixStack::TEXTURE);
        while (MatStack.getStackSize(tgt::MatrixStack::TEXTURE) != 1)
            MatStack.popMatrix();
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        warn(p, "Texture-matrix stack is not empty.");
    }

    if (MatStack.getModelViewMatrix() != tgt::mat4::createIdentity()){
        MatStack.loadIdentity();
        warn(p, "ModelView Matrix is not identity.");
    }

    if (MatStack.getProjectionMatrix() != tgt::mat4::createIdentity()){
        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.loadIdentity();
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        warn(p, "Projection Matrix is not identity.");
    }

    if (MatStack.getTextureMatrix() != tgt::mat4::createIdentity()){
        MatStack.matrixMode(tgt::MatrixStack::TEXTURE);
        MatStack.loadIdentity();
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        warn(p, "Texture Matrix is not identity.");
    }

    if (IMode.textureMode() != tgt::ImmediateMode::TEXNONE) {
        IMode.setTextureMode(tgt::ImmediateMode::TEXNONE);
        warn(p, "Immediate Mode texture mode has not been reset to TEXNONE.");
    }

    if (IMode.getCurrentColor() != tgt::vec4::one) {
        IMode.color(tgt::vec4::one);
        warn(p, "Immediate Mode color was not set to vec4::one.");
    }
}

void NetworkEvaluator::CheckOpenGLStateObserver::warn(Processor* p, const std::string& message) {
    if (p) {
        LWARNING(p->getClassName() << " (" << p->getID()
                 << "): invalid OpenGL state after processing: " << message);
    }
    else {
        LWARNING("Invalid OpenGL state before network processing: " << message);
    }
}

} // namespace
