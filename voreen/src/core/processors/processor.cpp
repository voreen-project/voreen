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

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/property.h"
#include "voreen/core/ports/coprocessorport.h"

#include "tgt/glcanvas.h"
#include "tgt/glmath.h"
#include "tgt/shadermanager.h"
#include "tgt/gpucapabilities.h"
#include "tgt/textureunit.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/voreenmodule.h"
#include "voreen/core/processors/processorwidgetfactory.h"

#include "voreen/core/processors/processorwidget.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/interaction/interactionhandler.h"

#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/io/progressbar.h"

#include <sstream>

using tgt::vec3;
using tgt::vec4;
using tgt::Color;

using std::string;
using std::map;
using std::vector;

namespace voreen {

const std::string Processor::loggerCat_("voreen.Processor");

const static std::string GENERIC_NOT_READY_ERROR_MSG("Not ready.");

Processor::Processor()
    : PropertyOwner("<unknown>", "<unknown>")
    , processorState_(PROCESSOR_STATE_NOT_INITIALIZED)
    , moduleName_("undefined")
    , processorWidget_(0)
    , invalidationVisited_(false)
    , eventVisited_(false)
    , firstProcessAfterDeserialization_(false)
    , locked_(false)
    , notReadyErrorMsg_(GENERIC_NOT_READY_ERROR_MSG)
{
    //setDescriptions(); // not allowed
}

Processor::~Processor() {
    if (isInitialized()) {
        LERROR("~Processor(): '" << getID() << "' has not been deinitialized");
        return;
    }
}

Processor* Processor::clone() const {
    try {
        std::stringstream stream;

        // first serialize
        XmlSerializer s;
        s.serialize("this", this);
        s.write(stream);

        // then deserialize again
        XmlDeserializer d;
        d.read(stream);
        Processor* proc = 0;
        d.deserialize("this", proc);

        proc->setDescriptions();
        return proc;
    }
    catch (std::exception& e) {
        LERROR("Failed to clone processor '" << getID() << "': " << e.what());
        return 0;
    }
}

void Processor::initialize() {
    if(description_.empty())
        setDescriptions();

    if (!VoreenApplication::app()) {
        LERROR("VoreenApplication not instantiated");
        throw VoreenException("VoreenApplication not instantiated");
    }

    if (isInitialized()) {
        LWARNING("initialize(): '" << getID() << "' already initialized");
        return;
    }

    bool glMode = VoreenApplication::app() && VoreenApplication::app()->isInitializedGL();
    if (!glMode)
        LDEBUG("initialize() not in OpenGL mode");

    // create and initialize processor widget
    processorWidget_ = VoreenApplication::app()->createProcessorWidget(this);
    if (processorWidget_) {
        processorWidget_->initialize();
        if (glMode)
            LGL_ERROR;
        // inform the observers about the new widget
        std::vector<ProcessorObserver*> observers = Observable<ProcessorObserver>::getObservers();
        for (size_t i = 0; i < observers.size(); ++i)
            observers[i]->processorWidgetCreated(this);
    }


    // initialize ports
    const std::vector<Port*>& ports = getPorts();
    for (size_t i=0; i < ports.size(); ++i) {
        if (!ports[i]->isInitialized())
            ports[i]->initialize();
        if (glMode)
            LGL_ERROR;
    }

    // initialize properties
    const std::vector<Property*>& properties = getProperties();
    for (size_t i=0; i < properties.size(); ++i) {
        properties[i]->initialize();
        if (glMode)
            LGL_ERROR;
    }
}

void Processor::deinitialize() {
    if (!isInitialized()) {
        LWARNING("deinitialize(): '" << getID() << "' (" << getClassName() << ") not initialized");
        return;
    }

    bool glMode = VoreenApplication::app() && VoreenApplication::app()->isInitializedGL();
    if (!glMode)
        LDEBUG("deinitialize() not in OpenGL mode");

    if (glMode)
        LGL_ERROR;

    // deinitialize ports
    const std::vector<Port*>& ports = getPorts();
    for (size_t i=0; i < ports.size(); ++i)
        ports[i]->deinitialize();
    if (glMode)
        LGL_ERROR;

    // deinitialize properties
    const std::vector<Property*>& properties = getProperties();
    for (size_t i = 0; i < properties.size(); ++i)
        properties[i]->deinitialize();
    if (glMode)
        LGL_ERROR;

    // delete processor widget
    delete processorWidget_;
    processorWidget_ = 0;

    // notify observers about the deleted widget
    std::vector<ProcessorObserver*> observers = Observable<ProcessorObserver>::getObservers();
    for (size_t i = 0; i < observers.size(); ++i)
        observers[i]->processorWidgetDeleted(this);
}

void Processor::beforeProcess() {
}

void Processor::afterProcess() {
    setValid();
    firstProcessAfterDeserialization_ = false;
}

void Processor::adjustPropertiesToInput() {

}

void Processor::clearOutports() {
    const std::vector<Port*>& outports = getOutports();
    for(size_t i=0; i<outports.size(); i++)
        outports[i]->clear();
}

void Processor::lockMutex() {
    mutex_.lock();
    locked_ = true;
}

void Processor::unlockMutex() {
    mutex_.unlock();
    locked_ = false;
}

bool Processor::isLocked() {
    return locked_;
}

bool Processor::isInitialized() const {
    return (processorState_ != PROCESSOR_STATE_NOT_INITIALIZED);
}

void Processor::addPort(Port* port) {
    tgtAssert(port, "Null pointer passed");
    port->setProcessor(this);
    CoProcessorPort* cpp = dynamic_cast<CoProcessorPort*>(port);
    if (port->isOutport()) {
        if (cpp)
            coProcessorOutports_.push_back(cpp);
        else
            outports_.push_back(port);
    }
    else {
        if (cpp)
            coProcessorInports_.push_back(cpp);
        else
            inports_.push_back(port);
    }

    map<std::string, Port*>::const_iterator it = portMap_.find(port->getID());
    if (it == portMap_.end())
        portMap_.insert(std::make_pair(port->getID(), port));
    else {
        LERROR("Port with name " << port->getID() << " has already been inserted!");
        tgtAssert(false, std::string("Port with name " + port->getID() + " has already been inserted").c_str());
    }
    notifyPortsChanged();
}

void Processor::addPort(Port& port) {
    addPort(&port);
}

void Processor::removePort(Port* port) {
    tgtAssert(port, "Null pointer passed");

    if (port->isInitialized()) {
        LWARNING("removePort() Port '" << getID() << "." << port->getID() << "' "
            << "has not been deinitialized");
    }

    //port->setProcessor(nullptr); // As per contract, we are only allowed to set the owner once.
    CoProcessorPort* cpp = dynamic_cast<CoProcessorPort*>(port);
    if (port->isOutport()) {
        if (cpp)
            coProcessorOutports_.erase(std::find(coProcessorOutports_.begin(), coProcessorOutports_.end(), port));
        else
            outports_.erase(std::find(outports_.begin(), outports_.end(), port));
    }
    else {
        if (cpp)
            coProcessorInports_.erase(std::find(coProcessorInports_.begin(), coProcessorInports_.end(), port));
        else
            inports_.erase(std::find(inports_.begin(), inports_.end(), port));
    }

    map<std::string, Port*>::iterator it = portMap_.find(port->getID());
    if (it != portMap_.end())
        portMap_.erase(it);
    else {
        LERROR("Port with name " << port->getID() << " was not found in port map!");
        tgtAssert(false, std::string("Port with name " + port->getID() + " was not found in port map!").c_str());
    }
    notifyPortsChanged();
}

void Processor::removePort(Port& port) {
    removePort(&port);
}

Processor::CodeState Processor::getCodeState() const {
    return Processor::CODE_STATE_EXPERIMENTAL;
}

Processor::ProcessorState Processor::getProcessorState() const {
    return processorState_;
}

void Processor::setGuiName(const std::string& guiName) {
    guiName_ = guiName;
    id_ = guiName;
    if (processorWidget_)
        processorWidget_->processorNameChanged();
}

void Processor::setID(const std::string& id) {
    setGuiName(id);
}

void Processor::setModuleName(const std::string& moduleName) {
    moduleName_ = moduleName;
}

std::string Processor::getModuleName() const {
    return moduleName_;
}

std::string Processor::getModulePath() const {
    return VoreenApplication::app()->getModulePath(getModuleName());
}

bool Processor::isUtility() const {
    return false;
}

const std::vector<Port*>& Processor::getInports() const {
    return inports_;
}

const std::vector<Port*>& Processor::getOutports() const {
    return outports_;
}

const std::vector<CoProcessorPort*>& Processor::getCoProcessorInports() const {
    return coProcessorInports_;
}

const std::vector<CoProcessorPort*>& Processor::getCoProcessorOutports() const {
    return coProcessorOutports_;
}

std::vector<Port*> Processor::getPorts() const {
    std::vector<Port*> result;

    const std::vector<Port*> inports = getInports();
    result.insert(result.begin(), inports.begin(), inports.end());

    const std::vector<Port*> outports = getOutports();
    result.insert(result.end(), outports.begin(), outports.end());

    const std::vector<CoProcessorPort*> coProcessorInports = getCoProcessorInports();
    result.insert(result.begin(), coProcessorInports.begin(), coProcessorInports.end());

    const std::vector<CoProcessorPort*> coProcessorOutports = getCoProcessorOutports();
    result.insert(result.end(), coProcessorOutports.begin(), coProcessorOutports.end());

    return result;
}

Port* Processor::getPort(const std::string& name) const {
    const std::vector<Port*> ports = getPorts();

    for (size_t i=0; i < ports.size(); i++) {
        if (ports[i]->getID() == name)
            return ports[i];
    }

    return 0;
}

const PerformanceRecord* Processor::getPerformanceRecord() const {
    return &performanceRecord_;
}

void Processor::resetPerformanceRecord() {
    performanceRecord_.deleteSamples();
}

void Processor::invalidate(int inv) {
    PropertyOwner::invalidate(inv);

    if (inv == Processor::VALID)
        return;

    if (!isInitialized())
        return;

    if (!invalidationVisited_) {
        invalidationVisited_ = true;

        for (size_t i=0; i<coProcessorOutports_.size(); ++i)
            coProcessorOutports_[i]->invalidatePort();

        invalidationVisited_ = false;

        tgtAssert(VoreenApplication::app(), "VoreenApplication not instantiated");
        VoreenApplication::app()->scheduleNetworkProcessing();
    }
    //check, if is ready (isInitialized)
    bool ready = isReady();
    if(ready && (processorState_ == PROCESSOR_STATE_NOT_READY)){
        processorState_ = PROCESSOR_STATE_READY;
        notifyStateChanged();
    } else if(!ready && (processorState_ == PROCESSOR_STATE_READY)) {
        processorState_ = PROCESSOR_STATE_NOT_READY;
        notifyStateChanged();
    }
}

bool Processor::isReady() const {
    if (!isInitialized())
        return false;

    for(size_t i=0; i<inports_.size(); ++i)
        if (!inports_[i]->isReady()) {
            setNotReadyErrorMessage("Inport(s) not connected!");
            return false;
        }

    for (size_t i=0; i<outports_.size(); ++i)
        if(!outports_[i]->isReady()) {
            setNotReadyErrorMessage("Outport(s) not connected!");
            return false;
        }

    for (size_t i=0; i<coProcessorInports_.size(); ++i)
        if (!coProcessorInports_[i]->isReady()) {
            setNotReadyErrorMessage("Coprocessorinport(s) not connected!");
            return false;
        }

    // Fall back to generic error message in case the child overwrites isReady(),
    // but does not set a better message (shame on them!).
    setNotReadyErrorMessage(GENERIC_NOT_READY_ERROR_MSG);

    return true;
}

bool Processor::isValid() const {
    return PropertyOwner::isValid();
}

void Processor::setValid() {
    PropertyOwner::setValid();

    for (size_t i=0; i<inports_.size(); ++i)
        inports_[i]->setValid();

    for (size_t i=0; i<coProcessorInports_.size(); ++i)
        coProcessorInports_[i]->setValid();
}

bool Processor::isEndProcessor() const {
    return (outports_.empty() && coProcessorOutports_.empty());
}

bool Processor::isSource() const {
    return (inports_.empty() && coProcessorInports_.empty() && coProcessorOutports_.empty());
}

void Processor::serialize(Serializer& s) const {
    // meta data
    metaDataContainer_.serialize(s);

    // misc settings
    s.serialize("name", id_);

    // serialize properties
    PropertyOwner::serialize(s);


    // ---
    // the following entities are static resources (i.e. already existing at this point)
    // that should therefore not be dynamically created by the serializer
    //
    const bool usePointerContentSerialization = s.getUsePointerContentSerialization();
    s.setUsePointerContentSerialization(true);

    // serialize inports using a temporary map
    std::map<std::string, Port*> inportMap;
    for (std::vector<Port*>::const_iterator it = inports_.begin(); it != inports_.end(); ++it)
        inportMap[(*it)->getID()] = *it;
    try {
        s.serialize("Inports", inportMap, "Port", "name");
    }
    catch (SerializationException& e) {
        LWARNING(e.what());
    }

    // serialize outports using a temporary map
    std::map<std::string, Port*> outportMap;
    for (std::vector<Port*>::const_iterator it = outports_.begin(); it != outports_.end(); ++it)
        outportMap[(*it)->getID()] = *it;
    try {
        s.serialize("Outports", outportMap, "Port", "name");
    }
    catch (SerializationException& e) {
        LWARNING(e.what());
    }

     // serialize coprocessorinports using a temporary map
    std::map<std::string, CoProcessorPort*> coProcessorInportMap;
    for (std::vector<CoProcessorPort*>::const_iterator it = coProcessorInports_.begin(); it != coProcessorInports_.end(); ++it)
        coProcessorInportMap[(*it)->getID()] = *it;
    try {
        s.serialize("CoProcessorInports", coProcessorInportMap, "Port", "name");
    }
    catch (SerializationException& e) {
        LWARNING(e.what());
    }

    // serialize coprocessoroutports using a temporary map
    std::map<std::string, CoProcessorPort*> coProcessorOutportMap;
    for (std::vector<CoProcessorPort*>::const_iterator it = coProcessorOutports_.begin(); it != coProcessorOutports_.end(); ++it)
        coProcessorOutportMap[(*it)->getID()] = *it;
    try {
        s.serialize("CoProcessorOutports", coProcessorOutportMap, "Port", "name");
    }
    catch (SerializationException& e) {
        LWARNING(e.what());
    }

    // serialize interaction handlers using a temporary map
    map<string, InteractionHandler*> handlerMap;
    const std::vector<InteractionHandler*>& handlers = getInteractionHandlers();
    for (vector<InteractionHandler*>::const_iterator it = handlers.begin(); it != handlers.end(); ++it)
        handlerMap[(*it)->getID()] = *it;
    try {
        s.serialize("InteractionHandlers", handlerMap, "Handler", "name");
    }
    catch (SerializationException& e) {
        LWARNING(e.what());
    }

    s.setUsePointerContentSerialization(usePointerContentSerialization);
    // --- static resources end ---

}

void Processor::deserialize(Deserializer& s) {
    // meta data
    metaDataContainer_.deserialize(s);

    // misc settings
    s.deserialize("name", id_);
    guiName_ = id_;

    // deserialize properties
    PropertyOwner::deserialize(s);


    // ---
    // the following entities are static resources that should not be dynamically created by the serializer
    //
    const bool usePointerContentSerialization = s.getUsePointerContentSerialization();
    s.setUsePointerContentSerialization(true);

    // deserialize inports using a temporary map
    map<string, Port*> inportMap;
    for (vector<Port*>::const_iterator it = inports_.begin(); it != inports_.end(); ++it)
        inportMap[(*it)->getID()] = *it;
    try {
        s.deserialize("Inports", inportMap, "Port", "name");
    }
    catch (SerializationNoSuchDataException& /*e*/){
        // port key missing => just ignore
        s.removeLastError();
    }

    // deserialize outports using a temporary map
    map<string, Port*> outportMap;
    for (vector<Port*>::const_iterator it = outports_.begin(); it != outports_.end(); ++it)
        outportMap[(*it)->getID()] = *it;
    try {
        s.deserialize("Outports", outportMap, "Port", "name");
    }
    catch (SerializationNoSuchDataException& /*e*/){
        // port key missing => just ignore
        s.removeLastError();
    }

    // serialize coprocessorinports using a temporary map
    std::map<std::string, CoProcessorPort*> coProcessorInportMap;
    for (std::vector<CoProcessorPort*>::const_iterator it = coProcessorInports_.begin(); it != coProcessorInports_.end(); ++it)
        coProcessorInportMap[(*it)->getID()] = *it;
    try {
        s.deserialize("CoProcessorInports", coProcessorInportMap, "Port", "name");
    }
    catch (SerializationNoSuchDataException& /*e*/){
        // port key missing => just ignore
        s.removeLastError();
    }

    // serialize coprocessoroutports using a temporary map
    std::map<std::string, CoProcessorPort*> coProcessorOutportMap;
    for (std::vector<CoProcessorPort*>::const_iterator it = coProcessorOutports_.begin(); it != coProcessorOutports_.end(); ++it)
        coProcessorOutportMap[(*it)->getID()] = *it;
    try {
        s.deserialize("CoProcessorOutports", coProcessorOutportMap, "Port", "name");
    }
    catch (SerializationNoSuchDataException& /*e*/){
        // port key missing => just ignore
        s.removeLastError();
    }


    // deserialize interaction handlers using a temporary map
    map<string, InteractionHandler*> handlerMap;
    const std::vector<InteractionHandler*>& handlers = getInteractionHandlers();
    for (vector<InteractionHandler*>::const_iterator it = handlers.begin(); it != handlers.end(); ++it)
        handlerMap[(*it)->getID()] = *it;
    try {
        s.deserialize("InteractionHandlers", handlerMap, "Handler", "name");
    }
    catch (SerializationNoSuchDataException& /*e*/){
        // interaction handler key missing => just ignore
        s.removeLastError();
    }

    s.setUsePointerContentSerialization(usePointerContentSerialization);
    // --- static resources end ---

    firstProcessAfterDeserialization_ = true;
}

MetaDataContainer& Processor::getMetaDataContainer() const {
    return metaDataContainer_;
}

std::string Processor::getCachePath() const {
    return VoreenApplication::app()->getCachePath() + "/" + getClassName();
}

std::string Processor::getDescription() const {
    if(description_.empty())
        const_cast<Processor*>(this)->setDescriptions(); //make setDescriptions const?

    return description_;
}

void Processor::setDescription(std::string desc) {
    description_ = desc;
}

std::string Processor::getNotReadyErrorMessage() const {
    return notReadyErrorMsg_;
}

void Processor::setNotReadyErrorMessage(const std::string& msg) const {
    notReadyErrorMsg_ = msg;
}

std::string Processor::getPropertyDescription(const std::string& propId) const {
    if(description_.empty())
        const_cast<Processor*>(this)->setDescriptions(); //make setDescriptions const?

    const Property* p = getProperty(propId);
    if(p)
        return p->getDescription();
    else {
        LWARNING("No such property: " << propId);
        return "";
    }
}

std::string Processor::getPortDescription(const std::string& portId) const {
    if(description_.empty())
        const_cast<Processor*>(this)->setDescriptions(); //make setDescriptions const?

    const Port* p = getPort(portId);
    if(p)
        return p->getDescription();
    else {
        LWARNING("No such port: " << portId);
        return "";
    }
}

bool Processor::usesExpensiveComputation() const {
    return false;
}

void Processor::setProgress(float progress) {
    for (size_t i=0; i<progressBars_.size(); i++)
        progressBars_.at(i)->setProgress(progress);
}
float Processor::getProgress() const {
    if (!progressBars_.empty())
        return progressBars_.front()->getProgress();
    else
        return 0.f;
}

void Processor::setProgressRange(const tgt::vec2& range) {
    for (size_t i=0; i<progressBars_.size(); i++)
        progressBars_.at(i)->setProgressRange(range);
}

tgt::vec2 Processor::getProgressRange() const {
    if (!progressBars_.empty())
        return progressBars_.front()->getProgressRange();
    else
        return tgt::vec2(0.f, 1.f);
}

void Processor::setProgressMessage(const std::string& message) {
    for (size_t i=0; i<progressBars_.size(); i++)
        progressBars_.at(i)->setProgressMessage(message);
}
std::string Processor::getProgressMessage() const {
    if (!progressBars_.empty())
        return progressBars_.front()->getProgressMessage();
    else
        return "";
}

void Processor::addProgressBar(ProgressReporter* progressReporter) {
    tgtAssert(progressReporter, "null pointer passed");
    progressBars_.push_back(progressReporter);
}

void Processor::removeProgressBar(ProgressReporter* progressReporter) {
    tgtAssert(progressReporter, "null pointer passed");
    progressBars_.erase(std::remove(progressBars_.begin(), progressBars_.end(), progressReporter), progressBars_.end());
}

ProcessorWidget* Processor::getProcessorWidget() const {
    return processorWidget_;
}

void Processor::addEventProperty(EventPropertyBase* prop) {
    tgtAssert(prop, "Null pointer passed");
    addProperty(prop);
    prop->setOwner(this);
    eventProperties_.push_back(prop);
}

void Processor::addEventProperty(EventPropertyBase& prop) {
    addEventProperty(&prop);
}

void Processor::addInteractionHandler(InteractionHandler* handler) {
    tgtAssert(handler, "Null pointer passed");
    handler->setOwner(this);
    interactionHandlers_.push_back(handler);

    if (handler->getEventProperties().empty())
        LWARNING("Interaction handler '" << handler->getID() << "' has no event property.");
}

void Processor::addInteractionHandler(InteractionHandler& handler) {
    addInteractionHandler(&handler);
}
const std::vector<EventPropertyBase*> Processor::getEventProperties() const {
    return eventProperties_;
}

const std::vector<InteractionHandler*>& Processor::getInteractionHandlers() const {
    return interactionHandlers_;
}

void Processor::onEvent(tgt::Event* e) {
    // cycle prevention
    if (eventVisited_)
        return;

    //if (e->isAccepted())
        //return;

    eventVisited_ = true;
    e->ignore();

    // propagate to interaction handlers (via their event properties)
    for (size_t i=0; i<interactionHandlers_.size() && !e->isAccepted(); ++i) {
        for (size_t j=0; j<interactionHandlers_[i]->getEventProperties().size() && !e->isAccepted(); ++j) {
            interactionHandlers_[i]->getEventProperties().at(j)->execute(e);
        }
    }

    // propagate to event properties
    for (size_t i = 0; (i < eventProperties_.size()) && !e->isAccepted(); ++i)
        eventProperties_[i]->execute(e);

    // pass-over event to processor's default event handling
    if (!e->isAccepted())
        tgt::EventListener::onEvent(e);

    // propagate to coprocessor inports
    for (size_t i=0; i<coProcessorInports_.size() && !e->isAccepted(); ++i)
        coProcessorInports_[i]->distributeEvent(e);

    // propagate to inports
    for (size_t i=0; i<inports_.size() && !e->isAccepted(); ++i)
        inports_[i]->distributeEvent(e);

    eventVisited_ = false;
}

void Processor::onPortEvent(tgt::Event* e, Port* p) {
    onEvent(e);
}

void Processor::deregisterWidget() {
    processorWidget_ = 0;
}

void Processor::initializePort(Port* port) {
    tgtAssert(port, "Null pointer passed");
    if (port->isInitialized()) {
        LWARNING("initializePort() port '" << getID() << "." << port->getID()
            << "' already initialized");
        return;
    }

    port->initialize();
}

void Processor::deinitializePort(Port* port) {
    tgtAssert(port, "Null pointer passed");
    if (!port->isInitialized()) {
        LWARNING("deinitializePort() port '" << getID() << "." << port->getID()
            << "' not initialized");
        return;
    }

    port->deinitialize();
}

void Processor::notifyPortsChanged() const {
    std::vector<ProcessorObserver*> procObservers = Observable<ProcessorObserver>::getObservers();
    for (size_t i = 0; i < procObservers.size(); ++i)
        procObservers[i]->portsChanged(this);
}

void Processor::notifyStateChanged() const {
    std::vector<ProcessorObserver*> procObservers = Observable<ProcessorObserver>::getObservers();
    for (size_t i = 0; i < procObservers.size(); ++i)
        procObservers[i]->stateChanged(this);
}

bool Processor::firstProcessAfterDeserialization() const {
    return firstProcessAfterDeserialization_;
}

void Processor::notifyPropertyValueHasBeenModified(Property* prop) const {
    PropertyOwner::notifyPropertyValueHasBeenModified(prop);
    //notify also the processor observer
    std::vector<ProcessorObserver*> propObservers = Observable<ProcessorObserver>::getObservers();
    for (size_t i = 0; i < propObservers.size(); ++i)
        propObservers[i]->propertyValueHasBeenModified(prop);
}

} // namespace voreen
