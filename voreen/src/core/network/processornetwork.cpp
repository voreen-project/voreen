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

#include "voreen/core/network/processornetwork.h"

#include "tgt/filesystem.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/voreenmodule.h"
#include "voreen/core/network/portconnection.h"
#include "voreen/core/network/networkgraph.h"
#include "voreen/core/network/workspace.h"
#include "voreen/core/ports/coprocessorport.h"
#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/link/propertylink.h"
#include "voreen/core/properties/link/linkevaluatorbase.h"
#include "voreen/core/properties/link/linkevaluatorhelper.h"

#include <vector>
#include <map>
#include <queue>
#include <sstream>

using std::vector;

namespace voreen {

const std::string ProcessorNetwork::loggerCat_("voreen.ProcessorNetwork");

namespace {

// Increase this counter when incompatible changes are introduced into the serialization format
const int NETWORK_VERSION = 19; //< Voreen 5.0

} // namespace

ProcessorNetwork::ProcessorNetwork()
    : version_(NETWORK_VERSION)
    , workspace_(0)
{}

ProcessorNetwork::~ProcessorNetwork() {
    clear();
}

void ProcessorNetwork::addProcessor(Processor* processor, const std::string& name) {
    tgtAssert(!processor->getClassName().empty(), "Processor's class name is empty");
    tgtAssert(!contains(processor), "Processor is already element of the network");

    std::string processorName;
    if (name.empty())
        processorName = processor->getClassName();
    else
        processorName = name;

    if (getProcessor(processorName))
        processorName = generateUniqueProcessorName(processorName);
    tgtAssert(!getProcessor(processorName), "Generated processor name should be unique here");

    // everything fine => rename and add processor
    processor->setID(processorName);
    processors_.push_back(processor);

    static_cast<Observable<ProcessorObserver>* >(processor)->addObserver(this);

    // notify observers
    notifyProcessorAdded(processor);
}

void ProcessorNetwork::removeProcessor(Processor* processor) {
    tgtAssert(processor, "null pointer passed");
    tgtAssert(contains(processor), "trying to remove processor that is not part of the network");

    // clear processor's connections and links
    disconnectPorts(processor);
    removePropertyLinksFromProcessor(processor);

    // remove processor from vector
    std::vector<Processor*>::iterator iter = std::find(processors_.begin(), processors_.end(), processor);
    tgtAssert(iter != processors_.end(), "processor not part of the network");
    processors_.erase(iter);

    notifyProcessorRemoved(processor);
    delete processor;
}

void ProcessorNetwork::clear() {
    while(!empty()) {
        tgtAssert(!processors_.empty(), "processors vector empty");
        removeProcessor(processors_.front());
    }
    errorList_.clear();

    // everything should be cleared now
    tgtAssert(processors_.empty() && propertyLinks_.empty(), "vectors should be empty");
}

MetaDataContainer& ProcessorNetwork::getMetaDataContainer() const {
    return metaDataContainer_;
}

void ProcessorNetwork::replaceProcessor(Processor* oldProc, Processor* newProc) {

    tgtAssert(oldProc && newProc, "Null pointer passed");

    std::vector<Port*> incomingPorts;
    std::vector<Port*> outgoingPorts;
    std::vector<Port*> incomingCoPorts;
    std::vector<Port*> outgoingCoPorts;

    {
        // "normal" ports
        const std::vector<Port*> inports = oldProc->getInports();
        for (size_t i = 0; i < inports.size(); ++i) {
            const std::vector<const Port*> connected = inports[i]->getConnected();
            for (size_t j = 0; j < connected.size(); ++j)
                incomingPorts.push_back(const_cast<Port*>(connected[j]));
        }

        const std::vector<Port*> outports = oldProc->getOutports();
        for (size_t i = 0; i < outports.size(); ++i) {
            const std::vector<const Port*> connected = outports[i]->getConnected();
            for (size_t j = 0; j < connected.size(); ++j)
                outgoingPorts.push_back(const_cast<Port*>(connected[j]));
        }
        // coprocessor ports
        const std::vector<CoProcessorPort*> coInports = oldProc->getCoProcessorInports();
        for (size_t i = 0; i < coInports.size(); ++i) {
            const std::vector<const Port*> connected = coInports[i]->getConnected();
            for (size_t j = 0; j < connected.size(); ++j)
                incomingCoPorts.push_back(const_cast<Port*>(connected[j]));
        }

        const std::vector<CoProcessorPort*> coOutports = oldProc->getCoProcessorOutports();
        for (size_t i = 0; i < coOutports.size(); ++i) {
            const std::vector<const Port*> connected = coOutports[i]->getConnected();
            for (size_t j = 0; j < connected.size(); ++j)
                outgoingCoPorts.push_back(const_cast<Port*>(connected[j]));
        }
    }

    removeProcessor(oldProc);
    addProcessor(newProc);

    {
        // "normal" ports
        const std::vector<Port*> inports = newProc->getInports();
        for (size_t i = 0; i < incomingPorts.size(); ++i) {
            for (size_t j = 0; j < inports.size(); ++j) {
                if (incomingPorts[i]->testConnectivity(inports[j]) != Port::ILLEGAL_CONNECTIVITY) {
                    connectPorts(incomingPorts[i], inports[j]);
                    break;
                }
            }
        }

        const std::vector<Port*> outports = newProc->getOutports();
        for (size_t i = 0; i < outgoingPorts.size(); ++i) {
            for (size_t j = 0; j < outports.size(); ++j) {
                if (outports[j]->testConnectivity(outgoingPorts[i]) != Port::ILLEGAL_CONNECTIVITY) {
                    connectPorts(outports[j], outgoingPorts[i]);
                    break;
                }
            }
        }
        // coprocessor ports
        const std::vector<CoProcessorPort*> coInports = newProc->getCoProcessorInports();
        for (size_t i = 0; i < incomingCoPorts.size(); ++i) {
            for (size_t j = 0; j < coInports.size(); ++j) {
                if (incomingCoPorts[i]->testConnectivity(coInports[j]) != Port::ILLEGAL_CONNECTIVITY) {
                    connectPorts(incomingCoPorts[i], coInports[j]);
                    break;
                }
            }
        }

        const std::vector<CoProcessorPort*> coOutports = newProc->getCoProcessorOutports();
        for (size_t i = 0; i < outgoingCoPorts.size(); ++i) {
            for (size_t j = 0; j < coOutports.size(); ++j) {
                if (coOutports[j]->testConnectivity(outgoingCoPorts[i]) != Port::ILLEGAL_CONNECTIVITY) {
                    connectPorts(coOutports[j], outgoingCoPorts[i]);
                    break;
                }
            }
        }
    }
}

const std::vector<Processor*>& ProcessorNetwork::getProcessors() const {
    return processors_;
}

bool ProcessorNetwork::contains(const Processor* processor) const {
    return (std::find(processors_.begin(), processors_.end(), processor) != processors_.end());
}

void ProcessorNetwork::setProcessorName(Processor* processor, const std::string& name) const {
    if (std::find(processors_.begin(), processors_.end(), processor) == processors_.end())
        throw VoreenException("Passed processor is not part of the network: " + processor->getID());

    if (name.empty())
        throw VoreenException("Passed processor name is empty");

    // processor already has the passed name => do nothing
    if (processor->getID() == name)
        return;

    // check whether passed name is already assigned to a processor in the network
    if (getProcessor(name)) {
        throw VoreenException("Passed processor name is not unique in the network: " + name);
    }

    // everything fine => rename processor
    std::string prevName = processor->getID();
    processor->setID(name);

    // notify observers
    notifyProcessorRenamed(processor, prevName);
}

std::string ProcessorNetwork::generateUniqueProcessorName(const std::string& name) const {

    if (!getProcessor(name))
        return name;

    int num = 2;
    while (getProcessor(name)) {
        std::ostringstream stream;
        stream << name << " " << num;
        std::string newName = stream.str();
        if (!getProcessor(newName))
            return newName;
        else
            num++;
    }

    tgtAssert(false, "should not get here");
    return "";
}

Processor* ProcessorNetwork::getProcessor(const std::string& name) const {
    for (size_t i=0; i<processors_.size(); ++i)
        if (processors_[i]->getID() == name)
            return processors_[i];

    return 0;
}

Processor* ProcessorNetwork::getProcessorByName(const std::string& name) const {
    return getProcessor(name);
}

size_t ProcessorNetwork::numProcessors() const {
    return processors_.size();
}

bool ProcessorNetwork::empty() const {
    return (numProcessors() == 0);
}


ProcessorNetwork* ProcessorNetwork::cloneSubNetwork(std::vector<Processor*> vec){
    // test, if all processors are from this PN
    for (size_t i =0; i < vec.size(); ++i){
        if (!contains(vec[i])){
            LERROR("cloneSubNetwork(): only processors from this ProcessorNetwork can be cloned");
            return 0;
        }
    }

    ProcessorNetwork* result = new ProcessorNetwork();

    // clone processors
    for (size_t i = 0; i < vec.size(); ++i)
        result->addProcessor(vec[i]->clone());

    // set connections
    for(size_t i = 0; i < vec.size(); ++i){ // each processor
        std::vector<Port*> outports = vec[i]->getOutports();
        for (size_t j = 0; j < outports.size(); ++j) { //each outport
            const std::vector<const Port*> inports = outports[j]->getConnected();
            for (size_t k = 0; k < inports.size(); ++k){ // each inport
                for (size_t l = 0; l < vec.size(); ++l){ // each processor
                    if (vec[l] == inports[k]->getProcessor()){
                        result->connectPorts(result->processors_[i]->getOutports()[j],
                            result->processors_[l]->getPort(inports[k]->getID()));
                        break;
                    }
                }
            }
        }
    }

    // set coports
    for(size_t i = 0; i < vec.size(); ++i){ //each processor
        std::vector<CoProcessorPort*> outports = vec[i]->getCoProcessorOutports();
        for (size_t j = 0; j < outports.size(); ++j) { //each outport
            const std::vector<const Port*> inports = outports[j]->getConnected();
            for (size_t k = 0; k < inports.size(); ++k){ //each inport
                for (size_t l = 0; l < vec.size(); ++l){ //each processor
                    if (vec[l] == inports[k]->getProcessor()){
                        result->connectPorts(result->processors_[i]->getCoProcessorOutports()[j],
                            result->processors_[l]->getPort(inports[k]->getID()));
                        break;
                    }
                }
            }
        }
    }

    // set links
    for (size_t i = 0; i < propertyLinks_.size(); ++i){
        Property* src = 0; Property* dest = 0;
        // check, if links are in the subnetwork
        for(size_t j = 0; j < vec.size();++j){
            if (vec[j] == propertyLinks_[i]->getDestinationProperty()->getOwner()){
                dest = result->getProcessors()[j]->getProperty(propertyLinks_[i]->getDestinationProperty()->getID());
            }
            if (vec[j] == propertyLinks_[i]->getSourceProperty()->getOwner()){
                src = result->getProcessors()[j]->getProperty(propertyLinks_[i]->getSourceProperty()->getID());
            }
        }
        if (src && dest) {
            LinkEvaluatorBase* eval = dynamic_cast<LinkEvaluatorBase*>(propertyLinks_[i]->getLinkEvaluator()->create());
            if (eval)
                result->createPropertyLink(src, dest, eval);
            else {
                LERROR(propertyLinks_[i]->getLinkEvaluator()->getClassName() << "::create() " << " did not return a LinkEvaluatorBase");
                delete eval;
            }
        }
    }

    result->notifyNetworkChanged();
    return result;
}

void ProcessorNetwork::mergeSubNetwork(ProcessorNetwork* pnw){
    // add processors
    for (size_t i=0; i < pnw->processors_.size(); ++i){
        addProcessor(pnw->processors_[i]);
    }

    // add connections
    for (size_t i=0; i < pnw->processors_.size(); ++i){
        std::vector<Port*> outports = pnw->processors_[i]->getOutports();
        for (size_t j = 0; j < outports.size(); ++j){
            const std::vector<const Port*> inports = outports[j]->getConnected();
            for (size_t k = 0; k < inports.size(); ++k){
                notifyPortConnectionAdded(outports[j], inports[k]);
            }
        }
    }

    // add coconnections
    for (size_t i=0; i < pnw->processors_.size(); ++i){
        std::vector<CoProcessorPort*> outports = pnw->processors_[i]->getCoProcessorOutports();
        for (size_t j = 0; j < outports.size(); ++j){
            const std::vector<const Port*> inports = outports[j]->getConnected();
            for (size_t k = 0; k < inports.size(); ++k){
                notifyPortConnectionAdded(outports[j], inports[k]);
            }
        }
    }

    for (size_t i=0; i < pnw->propertyLinks_.size(); ++i){
        propertyLinks_.push_back(pnw->propertyLinks_[i]);
        notifyPropertyLinkAdded(pnw->propertyLinks_[i]);
    }

    // delete passed network
    pnw->processors_.clear();
    pnw->propertyLinks_.clear();
    delete pnw;
}

bool ProcessorNetwork::connectPorts(Port* outport, Port* inport) {
    tgtAssert(inport, "passed null pointer");
    tgtAssert(inport->isInport(), "Inport expected");
    tgtAssert(outport, "passed null pointer");
    tgtAssert(outport->isOutport(), "Outport expected");
    tgtAssert(contains(outport->getProcessor()), "port not owned by a processor of this network");
    tgtAssert(contains(inport->getProcessor()), "port not owned by a processor of this network");

    bool result = outport->connect(inport);
    if (result)
        notifyPortConnectionAdded(outport, inport);

    return result;
}

bool ProcessorNetwork::connectPorts(Port& outport, Port& inport) {
    return connectPorts(&outport, &inport);
}

void ProcessorNetwork::disconnectPorts(Port* outport, Port* inport)  {
    tgtAssert(inport->isInport(), "Inport expected");
    tgtAssert(outport->isOutport(), "Outport expected");
    tgtAssert(contains(outport->getProcessor()), "port not owned by a processor of this network");
    tgtAssert(contains(inport->getProcessor()), "port not owned by a processor of this network");

    outport->disconnect(inport);
    notifyPortConnectionRemoved(outport, inport);
}

void ProcessorNetwork::disconnectPorts(Processor* processor) {
    tgtAssert(contains(processor), "processor not part of the network");

    const std::vector<Port*> inports = processor->getInports();
    for (size_t i = 0; i < inports.size(); ++i) {
        const std::vector<const Port*> connected = inports[i]->getConnected();
        for (size_t j = 0; j < connected.size(); ++j) {
            disconnectPorts(const_cast<Port*>(connected[j]), inports[i]);
        }
    }

    const std::vector<CoProcessorPort*> coProcessorInports = processor->getCoProcessorInports();
    for (size_t i = 0; i < coProcessorInports.size(); ++i) {
        const std::vector<const Port*> connected = coProcessorInports[i]->getConnected();
        for (size_t j = 0; j < connected.size(); ++j) {
            disconnectPorts(const_cast<Port*>(connected[j]), coProcessorInports[i]);
        }
    }

    const std::vector<Port*> outports = processor->getOutports();
    for (size_t i = 0; i < outports.size(); ++i) {
        const std::vector<const Port*> connected = outports[i]->getConnected();
        for (size_t j = 0; j < connected.size(); ++j) {
            disconnectPorts(outports[i], const_cast<Port*>(connected[j]));
        }
    }

    const std::vector<CoProcessorPort*> coProcessorOutports = processor->getCoProcessorOutports();
    for (size_t i = 0; i < coProcessorOutports.size(); ++i) {
        const std::vector<const Port*> connected = coProcessorOutports[i]->getConnected();
        for (size_t j = 0; j < connected.size(); ++j) {
            disconnectPorts(coProcessorOutports[i], const_cast<Port*>(connected[j]));
        }
    }
}

void ProcessorNetwork::disconnectPorts(Port& outport, Port& inport) {
    disconnectPorts(&outport, &inport);
}

void ProcessorNetwork::addProcessorInConnection(Port* outport, Port* inport, Processor* processor) {
    tgtAssert(inport, "passed null pointer");
    tgtAssert(inport->isInport(), "no inport");
    tgtAssert(outport, "passed null pointer");
    tgtAssert(outport->isOutport(), "no outport");
    tgtAssert(processor, "passed null pointer");

    addProcessor(processor);
    disconnectPorts(outport, inport);

    const std::vector<Port*> inports = processor->getInports();
    for (size_t i = 0; i < inports.size(); ++i) {
        if (outport->testConnectivity(inports[i]) != Port::ILLEGAL_CONNECTIVITY) {
            connectPorts(outport, inports[i]);
            break;
        }
    }

    const std::vector<Port*> outports = processor->getOutports();
    for (size_t i = 0; i < outports.size(); ++i) {
        if (outports[i]->testConnectivity(inport) != Port::ILLEGAL_CONNECTIVITY) {
            connectPorts(outports[i], inport);
            break;
        }
    }
}

void ProcessorNetwork::takeOverConnections(Port* outport, Processor* processor) {
    tgtAssert(outport, "passed null pointer");
    tgtAssert(outport->isOutport(), "no outport");
    tgtAssert(processor, "passed null pointer");

    addProcessor(processor);

    const std::vector<Port*> connected = outport->connectedPorts_;
    for (size_t i = 0; i < connected.size(); i++) {
        disconnectPorts(outport, connected[i]);
    }

    const std::vector<Port*> inports = processor->getInports();
    for (size_t i = 0; i < inports.size(); ++i) {
        if (outport->testConnectivity(inports[i]) != Port::ILLEGAL_CONNECTIVITY) {
            connectPorts(outport, inports[i]);
            break;
        }
    }

    for (size_t k = 0; k < connected.size(); ++k) {
        const std::vector<Port*> outports = processor->getOutports();
        for (size_t i = 0; i < outports.size(); ++i) {
            if (outports[i]->testConnectivity(connected[k]) != Port::ILLEGAL_CONNECTIVITY) {
                connectPorts(outports[i], connected[k]);
                break;
            }
        }
    }
}

bool ProcessorNetwork::canTakeOverConnections(Port* outport, Processor* processor) {
    tgtAssert(outport, "passed null pointer");
    tgtAssert(outport->isOutport(), "no outport");
    tgtAssert(processor, "passed null pointer");

    const std::vector<Port*> inports = processor->getInports();
    for (size_t i = 0; i < inports.size(); ++i) {
        if (outport->testConnectivity(inports[i]) != Port::ILLEGAL_CONNECTIVITY) {
            return true;
        }
    }

    const std::vector<Port*> connected = outport->connectedPorts_;
    for (size_t k = 0; k < connected.size(); ++k) {
        const std::vector<Port*> outports = processor->getOutports();
        for (size_t i = 0; i < outports.size(); ++i) {
            if (outports[i]->testConnectivity(connected[k]) != Port::ILLEGAL_CONNECTIVITY) {
                return true;
            }
        }
    }

    return false;
}

int ProcessorNetwork::removePropertyLinksFromProperty(Property* property) {
    tgtAssert(property, "null pointer passed");
    tgtAssert(dynamic_cast<Processor*>(property->getOwner()) || dynamic_cast<Port*>(property->getOwner()),
        "passed property is owned neither by a processor nor by a port");
    tgtAssert(contains(getOwningProcessor(property)), "passed property's owner is not part of the network");

    int numRemoved = 0;
    for (size_t i = 0; i < propertyLinks_.size(); ++i) {
        PropertyLink* link = propertyLinks_[i];
        if ((link->getSourceProperty() == property) || (link->getDestinationProperty() == property)) {
            removePropertyLink(link);
            i--; // removePropertyLink also deletes the link object
            numRemoved++;
        }
    }
    return numRemoved;
}

int ProcessorNetwork::removePropertyLinksFromProcessor(Processor* processor) {
    return removePropertyLinksFromProcessor<Property>(processor);
}

PropertyLink* ProcessorNetwork::createPropertyLink(Property* src, Property* dest, LinkEvaluatorBase* linkEvaluator) {
    tgtAssert(src && dest, "null pointer passed");
    tgtAssert(src != dest, "source and destination property are equal");

    tgtAssert(dynamic_cast<Processor*>(src->getOwner()) || dynamic_cast<Port*>(src->getOwner()),
        "passed src property is owned neither by a processor nor by a port");
    tgtAssert(contains(getOwningProcessor(src)), "passed src property's owner is not part of the network");

    tgtAssert(dynamic_cast<Processor*>(dest->getOwner()) || dynamic_cast<Port*>(dest->getOwner()),
        "passed dest property is owned neither by a processor nor by a port");
    tgtAssert(contains(getOwningProcessor(dest)), "passed dest property's owner is not part of the network");

    // check whether the properties are already linked
    if (linkEvaluator && containsPropertyLink(src, dest, linkEvaluator)) {
        LERROR("createPropertyLink(): this link already exists");
        return 0;
    }

    // check for event properties
    if (dynamic_cast<EventPropertyBase*>(src)) {
        LWARNING("Source property " << src->getID() << " is an EventProperty and can therefore not be linked.");
        delete linkEvaluator;
        return 0;
    }
    if (dynamic_cast<EventPropertyBase*>(dest)) {
        LWARNING("Dest property " << dest->getID() << " is an EventProperty and can therefore not be linked.");
        delete linkEvaluator;
        return 0;
    }

    // multiple links (with same direction) between two properties are not allowed
    if (src->isLinkedWith(dest, false)) {
        LWARNING(src->getFullyQualifiedID() << " and " <<
                dest->getFullyQualifiedID() << " are already linked");
        delete linkEvaluator;
        return 0;
    }

    // if no link evaluator passed, iterate over modules' link evaluator factories and use first compatible evaluator
    if (!linkEvaluator) {
        const std::vector<VoreenModule*>& modules = VoreenApplication::app()->getModules();
        for (size_t i=0; i<modules.size() && !linkEvaluator; i++) {
            std::vector<std::pair<std::string, std::string> > compatibleEvals = LinkEvaluatorHelper::getCompatibleLinkEvaluators(src, dest);
            if (!compatibleEvals.empty()) {
                LDEBUG("Creating LinkEvaluator: " << compatibleEvals.front().first);
                linkEvaluator = LinkEvaluatorHelper::createEvaluator(compatibleEvals.front().first);
                if (!linkEvaluator)
                    LWARNING("Failed to create LinkEvaluator: " << compatibleEvals.front().first);
            }
        }

        if (!linkEvaluator) {
            std::string t0 = src->getClassName();
            std::string t1 = dest->getClassName();
            if (t0 == t1)
                LWARNING("Could not create compatible LinkEvaluatorID for Property type: " << t0);
            else
                LWARNING("Could not create compatible LinkEvaluatorID for Property types: " << t0 << ", " << t1);
            return 0;
        }
    }

    tgtAssert(linkEvaluator, "no link evaluator");
    PropertyLink* link = new PropertyLink(src, dest, linkEvaluator);

    // if link is valid
    if (link->testPropertyLink()) {
        src->registerLink(link);
        propertyLinks_.push_back(link);
        notifyPropertyLinkAdded(link);
        return link;
    }
    else {
        LWARNING(src->getFullyQualifiedID() << " is not linkable to " << dest->getFullyQualifiedID());
        delete link;
        return 0;
    }

}

void ProcessorNetwork::removePropertyLink(PropertyLink* propertyLink) {
    tgtAssert(propertyLink, "null pointer passed");

    // remove property link if it is contained by the network
    bool propertyLinkRemoved = false;
    for (size_t i=0; i<propertyLinks_.size(); ++i) {
        if (propertyLink == propertyLinks_.at(i)) {
            propertyLinks_.erase(propertyLinks_.begin() + i);
            propertyLinkRemoved = true;
            break;
        }
    }

    // notify observers and then delete the link
    if (propertyLinkRemoved) {
        notifyPropertyLinkRemoved(propertyLink);
        delete propertyLink;
    }
}

bool ProcessorNetwork::removePropertyLink(Property* src, Property* dst) {
    tgtAssert(src, "null pointer passed");
    tgtAssert(dst, "null pointer passed");

    // remove property link if it is contained by the network
    bool propertyLinkRemoved = false;
    PropertyLink* propertyLink = 0;
    for (size_t i=0; i<propertyLinks_.size(); ++i) {
        if ((src == propertyLinks_.at(i)->getSourceProperty()) && (dst == propertyLinks_.at(i)->getDestinationProperty())) {
            propertyLink = propertyLinks_.at(i);
            propertyLinks_.erase(propertyLinks_.begin() + i);
            propertyLinkRemoved = true;
            break;
        }
    }

    // notify observers and then delete the link
    if (propertyLinkRemoved) {
        notifyPropertyLinkRemoved(propertyLink);
        delete propertyLink;
    }

    return propertyLinkRemoved;
}

const vector<PropertyLink*>& ProcessorNetwork::getPropertyLinks() const {
    return propertyLinks_;
}

int ProcessorNetwork::createPropertyLinksBetweenProperties(const std::vector<Property*>& properties,
        LinkEvaluatorBase* linkEvaluator, bool replace, bool transitive) {

    // filter properties that do not belong to a processor in the network
    std::vector<Property*> linkProperties;
    for (size_t i=0; i<properties.size(); ++i) {
        Processor* owner = getOwningProcessor(properties[i]);
        if (!owner) {
            LWARNING("Property is not (directly or indirectly) owned by a processor: " << properties[i]->getFullyQualifiedID());
            continue;
        }

        if (contains(owner))
            linkProperties.push_back(properties[i]);
        else
            LWARNING("Property's owner is not part of the network: " << properties[i]->getFullyQualifiedID());
    }

    // need at least two properties to link
    if (linkProperties.size() < 2)
        return 0;

    // remove existing links between properties, if in replace mode
    if (replace) {
        for (size_t i=0; i<linkProperties.size()-1; ++i) {
            for (size_t j=i+1; j<linkProperties.size(); ++j) {
                if (PropertyLink* link = linkProperties[i]->getLink(linkProperties[j]))
                    removePropertyLink(link);
                if (PropertyLink* link = linkProperties[j]->getLink(linkProperties[i]))
                    removePropertyLink(link);
            }
        }
    }

    int numLinks = 0;

    if (transitive) {
        // create a cyclic link chain, where each property is linked with its predecessor and successor
        for (size_t i=0; i<linkProperties.size()-1; ++i) {
            if (!linkProperties[i]->isLinkedWith(linkProperties[i+1])) {
                if (createPropertyLink(linkProperties[i], linkProperties[i+1], linkEvaluator ? dynamic_cast<LinkEvaluatorBase*>(linkEvaluator->create()) : 0))
                    numLinks++;
            }
            if (!linkProperties[i+1]->isLinkedWith(linkProperties[i])) {
                if (createPropertyLink(linkProperties[i+1], linkProperties[i], linkEvaluator ? dynamic_cast<LinkEvaluatorBase*>(linkEvaluator->create()) : 0))
                    numLinks++;
            }
        }
        // close cycle by linking last property with first
        if (linkProperties.size() > 2) {
            if (!linkProperties.back()->isLinkedWith(linkProperties.front())) {
                if (createPropertyLink(linkProperties.back(), linkProperties.front(), linkEvaluator ? dynamic_cast<LinkEvaluatorBase*>(linkEvaluator->create()) : 0))
                    numLinks++;
            }
            if (!linkProperties.front()->isLinkedWith(linkProperties.back())) {
                if (createPropertyLink(linkProperties.front(), linkProperties.back(), linkEvaluator ? dynamic_cast<LinkEvaluatorBase*>(linkEvaluator->create()) : 0))
                    numLinks++;
            }
        }
    }
    else {  // not transitive => link all properties with all others
        for (size_t i=0; i<linkProperties.size()-1; ++i) {
            for (size_t j=i+1; j<linkProperties.size(); ++j) {
                if (!linkProperties[i]->isLinkedWith(linkProperties[j])) {
                    if (createPropertyLink(linkProperties[i], linkProperties[j], linkEvaluator ? dynamic_cast<LinkEvaluatorBase*>(linkEvaluator->create()) : 0))
                        numLinks++;
                }
                if (!linkProperties[j]->isLinkedWith(linkProperties[i])) {
                    if (createPropertyLink(linkProperties[j], linkProperties[i], linkEvaluator ? dynamic_cast<LinkEvaluatorBase*>(linkEvaluator->create()) : 0))
                        numLinks++;
                }
            }
        }
    }

    return numLinks;
}

int ProcessorNetwork::removePropertyLinksFromSubNetwork(const std::vector<Processor*>& processors) {
    return removePropertyLinksFromSubNetwork<Property>(processors);
}

PropertyLink* ProcessorNetwork::getRenderSizeLink(RenderPort* source, RenderPort* dest) {
    tgtAssert(source && dest, "null pointer passed");

    Property* sourceProp = 0;
    if (source->getSizeOriginProperty())
        sourceProp = source->getSizeOriginProperty();
    else if (source->getSizeReceiveProperty())
        sourceProp = source->getSizeReceiveProperty();
    else
        return 0;

    Property* destProp = 0;
    if (dest->getSizeOriginProperty())
        destProp = dest->getSizeOriginProperty();
    else if (dest->getSizeReceiveProperty())
        destProp = dest->getSizeReceiveProperty();
    else
        return 0;

    tgtAssert(sourceProp && destProp, "unexpected null pointer");
    PropertyLink* link = sourceProp->getLink(destProp);
    if (link && dynamic_cast<LinkEvaluatorRenderSize*>(link->getLinkEvaluator()))
        return link;
    else
        return 0;
}

PropertyLink* ProcessorNetwork::getRenderSizeLink(RenderPort& source, RenderPort& destination) {
    return getRenderSizeLink(&source, &destination);
}

PropertyLink* ProcessorNetwork::createRenderSizeLink(RenderPort* source, RenderPort* dest) {
    tgtAssert(source && dest, "null pointer passed");

    Property* sourceProp = 0;
    if (source->getSizeOriginProperty())
        sourceProp = source->getSizeOriginProperty();
    else if (source->getSizeReceiveProperty())
        sourceProp = source->getSizeReceiveProperty();
    else {
        LWARNING("createRenderSizeLink(): passed source port is neither a size origin nor a size receiver");
        return 0;
    }

    Property* destProp = 0;
    if (dest->getSizeOriginProperty()) {
        destProp = dest->getSizeOriginProperty();
        if (source->getRenderSizePropagation() == RenderPort::RENDERSIZE_RECEIVER) {
            LWARNING("createRenderSizeLink(): render size link from a receiver to a origin is not possible");
            return 0;
        }
    }
    else if (dest->getSizeReceiveProperty())
        destProp = dest->getSizeReceiveProperty();
    else {
        LWARNING("createRenderSizeLink(): passed dest port is neither a size origin nor a size receiver");
        return 0;
    }

    tgtAssert(sourceProp && destProp, "unexpected null pointer");
    if (getRenderSizeLink(source, dest)) {
        LWARNING("createRenderSizeLink(): size link from '" << source->getQualifiedName() << "' to '" <<
            dest->getQualifiedName() << "' already exists.");
        return 0;
    }
    else {
        return createPropertyLink(sourceProp, destProp, new LinkEvaluatorRenderSize());
    }
}

PropertyLink* ProcessorNetwork::createRenderSizeLink(RenderPort& source, RenderPort& destination) {
    return createRenderSizeLink(&source, &destination);
}

bool ProcessorNetwork::removeRenderSizeLink(RenderPort* source, RenderPort* dest) {
    tgtAssert(source && dest, "null pointer passed");
    PropertyLink* link = getRenderSizeLink(source, dest);
    if (link) {
        removePropertyLink(link);
        return true;
    }
    else
        return false;
}

bool ProcessorNetwork::removeRenderSizeLink(RenderPort& source, RenderPort& destination) {
    return removeRenderSizeLink(&source, &destination);
}

int ProcessorNetwork::removeRenderSizeLinksFromProcessor(Processor* processor) {
    tgtAssert(processor, "null pointer passed");
    tgtAssert(contains(processor), "passed processor is not part of the network");

    // collect render ports of processor
    std::vector<RenderPort*> processorPorts;
    std::vector<Port*> tPorts = processor->getPorts();
    for (size_t i=0; i<tPorts.size(); i++)
        if (RenderPort* renderPort = dynamic_cast<RenderPort*>(tPorts.at(i)))
            processorPorts.push_back(renderPort);

    // iterate over render ports of all processors and remove links to/from processorPorts
    int numRemoved = 0;
    for (size_t i=0; i<processors_.size(); i++) {
        std::vector<Port*> tPorts = processors_.at(i)->getPorts();
        for (size_t j=0; j<tPorts.size(); j++) {
            RenderPort* networkPort = dynamic_cast<RenderPort*>(tPorts.at(j));
            if (networkPort) {
                for (size_t k=0; k<processorPorts.size(); k++) {
                    if (removeRenderSizeLink(networkPort, processorPorts.at(k)))
                        numRemoved++;
                    if (removeRenderSizeLink(processorPorts.at(k), networkPort))
                        numRemoved++;
                }
            }
        }
    }

    return numRemoved;
}

int ProcessorNetwork::createRenderSizeLinksWithinSubNetwork(const std::vector<Processor*>& linkProcessors, bool replaceExisting /*= false*/) {
    // check, if all processors are part of this network
    for (size_t i=0; i < linkProcessors.size(); ++i){
        if (!contains(linkProcessors[i])){
            LERROR("createRenderSizeLinksInSubNetwork(): processor '" << linkProcessors[i]->getID() << "' is not part of this network");
            tgtAssert(false, "processor is not part of this network");
            return 0;
        }
    }

    // create inverse map from render size receivers to linked origins for existing links
    std::map<RenderPort*, std::vector<RenderPort*> > receiverToOriginMapExisting =
        getRenderSizeReceiverToOriginsMap(processors_);

    // determine for each render receiver port (dest) of the passed processors the corresponding render size origins (sources)
    // it should be linked with
    std::map<RenderPort*, std::vector<RenderPort*> > receiverToOriginMap;
    for (std::vector<Processor*>::const_iterator procIt=linkProcessors.begin(); procIt != linkProcessors.end(); ++procIt) {
        Processor* proc = *procIt;
        std::vector<Port*> outports = proc->getOutports();
        for (std::vector<Port*>::iterator outportIt=outports.begin(); outportIt != outports.end(); ++outportIt) {
            RenderPort* rootPort = dynamic_cast<RenderPort*>(*outportIt);
            if (!rootPort || rootPort->getRenderSizePropagation() != RenderPort::RENDERSIZE_RECEIVER)
                continue;
            receiverToOriginMap.insert(std::make_pair(rootPort, getSuccessingRenderSizeOrigins(rootPort, linkProcessors)));
        }
    }

    // create render size links
    int numCreated = 0;
    std::map<RenderPort*, std::vector<RenderPort*> >::iterator receiverIt;
    for (receiverIt = receiverToOriginMap.begin(); receiverIt != receiverToOriginMap.end(); ++receiverIt) {
        RenderPort* receiver = receiverIt->first;
        tgtAssert(receiver->getRenderSizePropagation() == RenderPort::RENDERSIZE_RECEIVER, "not a render size receiver");

        // check if origin has been found
        std::vector<RenderPort*> origins = receiverIt->second;
        if (origins.empty())
            continue;

        // remove duplicates in origins vector
        std::vector<RenderPort*>::iterator uniqueIt = std::unique(origins.begin(), origins.end());
        origins.resize(uniqueIt - origins.begin());

        // check, if receiver already has a size link
        if (receiverToOriginMapExisting.count(receiver) && !receiverToOriginMapExisting[receiver].empty()) {
            std::vector<RenderPort*> existingOrigins = receiverToOriginMapExisting[receiver];
            if (existingOrigins.size() == 1 && origins.size() == 1 && existingOrigins.front() == origins.front()) {
                // single size link to be created for this receiver already exists => do nothing
                continue;
            }
            else if (!replaceExisting) {
                LINFO("Render size receiver '" << receiver->getQualifiedName() << "' already has size link(s). No render size link created.");
                continue;
            }
            else {
                // remove existing size links to this receiver
                LWARNING("Removing existing size links to receiver: " << receiver->getQualifiedName());
                for (size_t i=0; i<existingOrigins.size(); i++)
                    removeRenderSizeLink(existingOrigins.at(i), receiver);
            }
        }

        // create size link, if exactly one size origin has been found for this receiver
        if (origins.size() == 1) {
            RenderPort* origin = origins.front();
            tgtAssert(origin->getRenderSizePropagation() == RenderPort::RENDERSIZE_ORIGIN, "not a render size origin");
            tgtAssert(receiver->getRenderSizePropagation() == RenderPort::RENDERSIZE_RECEIVER, "not a render size receiver");
            createRenderSizeLink(origin, receiver);
            numCreated++;
        }
        else {
            LWARNING("Multiple potential render size origins for render size receiver: " << receiver->getQualifiedName());
        }
    }

    return numCreated;
}

int ProcessorNetwork::createRenderSizeLinksOverConnection(RenderPort* outport, RenderPort* inport, bool replaceExisting) {
    tgtAssert(outport && inport, "null pointer passed");
    tgtAssert(outport->isOutport(), "passed port is not an outport");
    tgtAssert(inport->isInport(), "passed port is not an inport");
    tgtAssert(outport->isConnectedTo(inport), "passed ports are not connected");

    // determine size origin successor of inport
    RenderPort* originPort = 0;
    if (inport->getRenderSizePropagation() == RenderPort::RENDERSIZE_ORIGIN)
        originPort = inport;
    else {
        std::vector<RenderPort*> originSuccessors = getSuccessingRenderSizeOrigins(inport, processors_);
        if (originSuccessors.size() == 1)
            originPort = originSuccessors.front();
        else if (originSuccessors.size() > 1) {
            // multiple origins => no size link over connection possible
            LWARNING("Multiple potential render size origins. No render size link created.");
            return 0;
        }
        else {
            // no origin => no size link over connection possible
            return 0;
        }
    }
    tgtAssert(originPort, "no origin port");

    // create inverse map from render size receivers to linked origins for existing links
    std::map<RenderPort*, std::vector<RenderPort*> > receiverToOriginMapExisting =
        getRenderSizeReceiverToOriginsMap(processors_);

    // determine size receiving predecessors of outport
    std::vector<RenderPort*> receiverPredecessors;
    if (outport->getRenderSizePropagation() == RenderPort::RENDERSIZE_RECEIVER)
        receiverPredecessors.push_back(outport);
    else
        receiverPredecessors = getPredecessingRenderSizeReceivers(outport, processors_);

    // determine for each size predecessor (receiver) its size successors (origins)
    //std::map<RenderPort*, RenderPort*> receiverToOriginMap;
    std::vector<RenderPort*> receiverPortsToLink;
    for (size_t i=0; i<receiverPredecessors.size(); i++) {
        RenderPort* receiverPort = receiverPredecessors.at(i);
        tgtAssert(receiverPort->getRenderSizePropagation() == RenderPort::RENDERSIZE_RECEIVER, "port is not a size receiver");

        // get size origin successors considering the entire network topology
        std::vector<RenderPort*> tOrigins = getSuccessingRenderSizeOrigins(receiverPort, processors_);
        // if receiver has only one potential origin, which equals the connection successor, create size link
        if (tOrigins.size() == 1 && tOrigins.front() == originPort)
            receiverPortsToLink.push_back(receiverPort);
        else if (tOrigins.size() > 1) {
            LWARNING("Multiple potential render size origins for render size receiver: " << receiverPort->getQualifiedName()
                << ". No render size link created.");
        }
    }

    // create render size links
    int numCreated = 0;
    for (size_t i=0; i<receiverPortsToLink.size(); i++) {
        RenderPort* receiverPort = receiverPortsToLink.at(i);
        tgtAssert(receiverPort->getRenderSizePropagation() == RenderPort::RENDERSIZE_RECEIVER, "not a render size receiver");

        // check, if receiver already has a size link
        if (receiverToOriginMapExisting.count(receiverPort) && !receiverToOriginMapExisting[receiverPort].empty()) {
            std::vector<RenderPort*> existingOrigins = receiverToOriginMapExisting[receiverPort];
            if (existingOrigins.size() == 1 && existingOrigins.front() == originPort) {
                // single size link to be created for this receiver already exists => do nothing
                continue;
            }
            else if (!replaceExisting) {
                LINFO("Render size receiver '" << receiverPort->getQualifiedName() << "' already has size link(s). No render size link created.");
                continue;
            }
            else {
                // remove existing size links to this receiver
                LWARNING("Removing existing size links to receiver: " << receiverPort->getQualifiedName());
                for (size_t j=0; j<existingOrigins.size(); j++)
                    removeRenderSizeLink(existingOrigins.at(j), receiverPort);
            }
        }

        // create size link
        tgtAssert(originPort->getRenderSizePropagation() == RenderPort::RENDERSIZE_ORIGIN, "not a render size origin");
        createRenderSizeLink(originPort, receiverPort);
        numCreated++;
    }

    return numCreated;
}

int ProcessorNetwork::createRenderSizeLinksForProcessor(Processor* processor, bool replaceExisting) {
    tgtAssert(processor, "null pointer passed");
    tgtAssert(contains(processor), "passed processor is not part of the network");

    // call createRenderSizeLinksOverConnection() for all incoming and outgoing render port connections
    int numCreated = 0;

    // render inports
    std::vector<Port*> inports = processor->getInports();
    for (size_t i=0; i<inports.size(); i++) {
        RenderPort* renderInport = dynamic_cast<RenderPort*>(inports.at(i));
        if (renderInport) {
            std::vector<const Port*> connectedOutports = renderInport->getConnected();
            for (size_t j=0; j<connectedOutports.size(); j++) {
                const RenderPort* renderOutport = dynamic_cast<const RenderPort*>(connectedOutports.at(j));
                tgtAssert(renderOutport, "RenderPort connected to non-RenderPort");
                numCreated += createRenderSizeLinksOverConnection(const_cast<RenderPort*>(renderOutport), renderInport, replaceExisting);
            }
        }
    }

    // render outport
    std::vector<Port*> outports = processor->getOutports();
    for (size_t i=0; i<outports.size(); i++) {
        RenderPort* renderOutport = dynamic_cast<RenderPort*>(outports.at(i));
        if (renderOutport) {
            std::vector<const Port*> connectedInports = renderOutport->getConnected();
            for (size_t j=0; j<connectedInports.size(); j++) {
                const RenderPort* renderInport = dynamic_cast<const RenderPort*>(connectedInports.at(j));
                tgtAssert(renderInport, "RenderPort connected to non-RenderPort");
                numCreated += createRenderSizeLinksOverConnection(renderOutport, const_cast<RenderPort*>(renderInport), replaceExisting);
            }
        }
    }

    return numCreated;
}

int ProcessorNetwork::removeRenderSizeLinksFromSubNetwork(const std::vector<Processor*>& subNetwork) {
    // check, if all processors are part of this network
    for (size_t i=0; i < subNetwork.size(); ++i){
        if (!contains(subNetwork[i])){
            LERROR("removeRenderSizeLinksFromSubNetwork(): processor '" << subNetwork[i]->getID() << "' is not part of this network");
            tgtAssert(false, "removeRenderSizeLinksFromSubNetwork(): processor is not part of this network");
            return 0;
        }
    }

    int numRemoved = 0;

    // collect render size receivers in passed sub network
    std::set<Property*> sizeReceivers;
    for (size_t i=0; i<subNetwork.size(); i++) {
        Processor* proc = subNetwork.at(i);
        std::vector<Port*> outports = proc->getOutports();
        for (size_t j=0; j<outports.size(); j++) {
            RenderPort* rp = dynamic_cast<RenderPort*>(outports.at(j));
            if (rp && rp->getRenderSizePropagation() == RenderPort::RENDERSIZE_RECEIVER) {
                tgtAssert(rp->getSizeReceiveProperty(), "render size receiver without receive property");
                sizeReceivers.insert(rp->getSizeReceiveProperty());
            }

        }
    }

    // iterate over size origins in sub network and remove all size links to collected size receivers
    for (size_t i=0; i<subNetwork.size(); i++) {
        Processor* proc = subNetwork.at(i);
        std::vector<Port*> inports = proc->getInports();
        for (size_t j=0; j<inports.size(); j++) {
            RenderPort* rp = dynamic_cast<RenderPort*>(inports.at(j));
            if (rp && rp->getRenderSizePropagation() == RenderPort::RENDERSIZE_ORIGIN) {
                tgtAssert(rp->getSizeOriginProperty(), "render size origin without origin property");
                std::vector<PropertyLink*> links = rp->getSizeOriginProperty()->getLinks();
                for (size_t k=0; k<links.size(); k++) {
                    PropertyLink* link = links.at(k);
                    if (dynamic_cast<LinkEvaluatorRenderSize*>(link->getLinkEvaluator()) &&
                            sizeReceivers.count(link->getDestinationProperty())) {
                        removePropertyLink(link);
                        numRemoved++;
                    }
                }
            }
        }
    }

    return numRemoved;
}

int ProcessorNetwork::removeRenderSizeLinksOverConnection(RenderPort* outport, RenderPort* inport) {
    tgtAssert(outport->isOutport(), "passed port is not an outport");
    tgtAssert(inport->isInport(), "passed port is not an inport");

    // determine size origin successor of inport
    RenderPort* originPort = 0;
    if (inport->getRenderSizePropagation() == RenderPort::RENDERSIZE_ORIGIN)
        originPort = inport;
    else {
        std::vector<RenderPort*> originSuccessors = getSuccessingRenderSizeOrigins(inport, processors_);
        if (originSuccessors.size() == 1)
            originPort = originSuccessors.front();
        else if (originSuccessors.size() > 1) {
            // multiple origins => no size link over connection possible
            LWARNING("Multiple potential render size origins. No render size link removed.");
            return 0;
        }
        else {
            // no origin => no size link over connection possible
            return 0;
        }
    }
    tgtAssert(originPort, "no origin port");

    // determine size receiving predecessors of outport
    std::set<RenderPort*> receiverPredecessors;
    if (outport->getRenderSizePropagation() == RenderPort::RENDERSIZE_RECEIVER)
        receiverPredecessors.insert(outport);
    else {
        std::vector<RenderPort*> preds = getPredecessingRenderSizeReceivers(outport, processors_);
        receiverPredecessors.insert(preds.begin(), preds.end());
    }
    if (receiverPredecessors.empty())
        return 0;

    // remove all size links from origins to receivers
    int numRemoved = 0;
    RenderSizeOriginProperty* originProp = originPort->getSizeOriginProperty();
    tgtAssert(originProp, "render size origin without size origin property");
    std::vector<PropertyLink*> links = originProp->getLinks();
    for (size_t i=0; i<links.size(); i++) {
        PropertyLink* link = links.at(i);
        RenderSizeReceiveProperty* receiveProp = dynamic_cast<RenderSizeReceiveProperty*>(link->getDestinationProperty());
        if (receiveProp && dynamic_cast<LinkEvaluatorRenderSize*>(link->getLinkEvaluator())) {
            RenderPort* receiverPort = dynamic_cast<RenderPort*>(receiveProp->getOwner());
            tgtAssert(receiverPort, "RenderSizeReceiveProperty is not owned by a RenderPort");
            if (receiverPredecessors.count(receiverPort)) {
                removePropertyLink(link);
                numRemoved++;
            }
        }
    }

    return numRemoved;
}

int ProcessorNetwork::removeAllPropertyLinks() {
    std::vector<PropertyLink*> links = getPropertyLinks();
    int numRemoved = static_cast<int>(links.size());
    while (!propertyLinks_.empty()) {
        removePropertyLink(propertyLinks_.front());
    }
    return numRemoved;
}

bool ProcessorNetwork::isPropertyLinked(const Property* property) const {
    for (size_t i = 0; i < propertyLinks_.size(); ++i) {
        PropertyLink* link = propertyLinks_[i];
        if ((link->getSourceProperty() == property) || (link->getDestinationProperty() == property))
            return true;
    }
    return false;
}

bool ProcessorNetwork::containsPropertyLink(const Property* src, const Property* dest, LinkEvaluatorBase* evaluator) const {
    for (size_t i = 0; i < propertyLinks_.size(); ++i) {
        PropertyLink* link = propertyLinks_[i];
        bool sourcePropertyMatches = (link->getSourceProperty() == src);
        bool destinationPropertyMatches = (link->getDestinationProperty() == dest);
        bool evaluatorMatches = (typeid(*(link->getLinkEvaluator())) == typeid(*evaluator));

        if (sourcePropertyMatches && destinationPropertyMatches && evaluatorMatches)
            return true;
    }

    return false;
}

void ProcessorNetwork::setVersion(int version) {
    version_ = version;
}

int ProcessorNetwork::getVersion() const {
    return version_;
}

void ProcessorNetwork::serialize(Serializer& s) const {
    const_cast<ProcessorNetwork*>(this)->errorList_.clear();

    // meta data
    metaDataContainer_.serialize(s);

    // misc
    s.serialize("version", version_);

    // Serialize Processors and add them to the network element
    s.serialize("Processors", processors_, "Processor");

    // Serialize the port connections...
    std::vector<PortConnection> portConnections;
    std::vector<PortConnection> coProcessorPortConnections;

    typedef vector<Port*> Ports;
    typedef vector<const Port*> ConPorts;
    typedef vector<CoProcessorPort*> CoProcessorPorts;

    for (vector<Processor*>::const_iterator processorIt = processors_.begin(); processorIt != processors_.end(); ++processorIt) {
        Ports ports = (*processorIt)->getOutports();
        for (Ports::const_iterator portIt = ports.begin(); portIt != ports.end(); ++portIt) {
            const ConPorts connectedPorts = (*portIt)->getConnected();
            for (ConPorts::const_iterator connectedPortIt = connectedPorts.begin(); connectedPortIt != connectedPorts.end(); ++connectedPortIt)
                portConnections.push_back(PortConnection(*portIt, const_cast<Port*>(*connectedPortIt)));
        }

        CoProcessorPorts coProcessorPorts = (*processorIt)->getCoProcessorOutports();
        for (CoProcessorPorts::const_iterator portIt = coProcessorPorts.begin(); portIt != coProcessorPorts.end(); ++portIt) {
            const ConPorts connectedPorts = (*portIt)->getConnected();
            for (ConPorts::const_iterator connectedPortIt = connectedPorts.begin(); connectedPortIt != connectedPorts.end(); ++connectedPortIt)
                coProcessorPortConnections.push_back(PortConnection(*portIt, const_cast<Port*>(*connectedPortIt)));
        }
    }

    s.serialize("Connections", portConnections, "Connection");
    s.serialize("CoProcessorConnections", coProcessorPortConnections, "CoProcessorConnection");

    // Serialize PropertyLinks and add them to the network element
    s.serialize("PropertyLinks", propertyLinks_, "PropertyLink");
}

void ProcessorNetwork::deserialize(Deserializer& s) {
    errorList_.clear();

    // deserialize Processors
    s.deserialize("Processors", processors_, "Processor");
    for (size_t i=0; i<processors_.size(); i++)
        static_cast<Observable<ProcessorObserver>* >(processors_[i])->addObserver(this);

    // deserialize port connections...
    std::vector<PortConnection> portConnections;
    std::vector<PortConnection> coProcessorPortConnections;

    s.deserialize("Connections", portConnections, "Connection");
    s.deserialize("CoProcessorConnections", coProcessorPortConnections, "CoProcessorConnection");

    for (std::vector<PortConnection>::iterator it = portConnections.begin(); it != portConnections.end(); ++it) {
        if (it->getOutport() && it->getInport())
            it->getOutport()->connect(it->getInport());
        else {
            std::string addOn;
            if (it->getOutport()) {
                addOn += "Outport: '";
                if (it->getOutport()->getProcessor())
                    addOn += it->getOutport()->getProcessor()->getID() + "::";
                addOn += it->getOutport()->getID() + "' ";
            }
            if (it->getInport()) {
                addOn += "Inport: '";
                if (it->getInport()->getProcessor())
                    addOn += it->getInport()->getProcessor()->getID() + "::";
                addOn += it->getInport()->getID() + "'";
            }
            s.addError(SerializationException("Port connection could not be established. " + addOn));
        }
    }

    for (std::vector<PortConnection>::iterator it = coProcessorPortConnections.begin(); it != coProcessorPortConnections.end(); ++it) {
        if (it->getOutport() && it->getInport())
            it->getOutport()->connect(it->getInport());
        else
            s.addError(SerializationException("Coprocessor port connection could not be established."));
    }

    // deserialize PropertyLinks
    s.deserialize("PropertyLinks", propertyLinks_, "PropertyLink");
    // remove links without evaluator (may happen due to errors during deserialization)
    for (size_t i=0; i<propertyLinks_.size(); /*no auto-increment*/) {
        if (!propertyLinks_[i]->getLinkEvaluator()) {
            delete propertyLinks_[i];
            propertyLinks_.erase(propertyLinks_.begin() + i);
        }
        else {
            i++;
        }
    }

    // check for duplicate processor names
    vector<vector<Processor*> > duplicates;
    for (size_t i=0; i<processors_.size(); ++i) {
        Processor* processor = processors_[i];
        bool isDuplicate = false;
        for (size_t j=0; j<processors_.size() && !isDuplicate; ++j) {
            if ((i != j) && (processor->getID() == processors_[j]->getID()))
                isDuplicate = true;
        }
        if (isDuplicate) {
            // insert into duplicates vector
            bool newNameGroup = true;
            for (size_t j=0; j<duplicates.size(); ++j) {
                if (duplicates[j][0]->getID() == processor->getID()) {
                    duplicates[j].push_back(processor);
                    newNameGroup = false;
                    break;
                }
            }
            if (newNameGroup) {
                duplicates.push_back(vector<Processor*>());
                duplicates.back().push_back(processor);
            }
        }
    }

    // correct duplicate processor names
    for (size_t i=0; i<duplicates.size(); ++i) {
        LWARNING("Duplicate processor name: \"" << duplicates[i][0]->getID() << "\" (" << duplicates[i].size() << " occurrences)");
        for (size_t j=1; j<duplicates[i].size(); ++j) {
            std::ostringstream stream;
            stream << duplicates[i][j]->getID() << " " << j+1;
            std::string newName = stream.str();
            LINFO("Renaming \"" << duplicates[i][j]->getID() << "\" to \"" << newName << "\"");
            duplicates[i][j]->setID(newName);
        }
    }

    // meta data
    metaDataContainer_.deserialize(s);

    notifyNetworkChanged();
}

std::vector<std::string> ProcessorNetwork::getErrors() const {
    return errorList_;
}

void ProcessorNetwork::setErrors(const std::vector<std::string>& errorList) {
    errorList_ = errorList;
}

Processor* ProcessorNetwork::getOwningProcessor(Property* property) const {
    if (Processor* processor = dynamic_cast<Processor*>(property->getOwner()))
        return processor;
    else if (Port* port = dynamic_cast<Port*>(property->getOwner()))
        return port->getProcessor();
    else
        return 0;
}

std::vector<RenderPort*> ProcessorNetwork::getPredecessingRenderSizeReceivers(RenderPort* rootPort, const std::vector<Processor*>& subNetwork) const {
    tgtAssert(rootPort, "null pointer passed");

    std::vector<RenderPort*> result;

    std::queue<RenderPort*> worklist;
    if (rootPort->isInport())
        worklist.push(rootPort);
    else {
        std::vector<Port*> inports = rootPort->getProcessor()->getInports();
        for (size_t i=0; i<inports.size(); i++) {
            RenderPort* rp = dynamic_cast<RenderPort*>(inports.at(i));
            if (rp)
                worklist.push(rp);
        }
    }

    std::set<RenderPort*> processed; //< for cycle detection
    while (!worklist.empty()) {
        RenderPort* curPort = worklist.front();
        worklist.pop();
        tgtAssert(curPort, "null pointer in worklist");
        tgtAssert(curPort->isInport(), "outport in worklist");
        processed.insert(curPort);
        if (curPort->isLoopPort()) //< do not traverse over loop ports
            continue;
        const std::vector<const Port*> predecessors = curPort->getConnected();
        for (std::vector<const Port*>::const_iterator succIt=predecessors.begin(); succIt != predecessors.end(); ++succIt) {
            const RenderPort* predecessorPort = dynamic_cast<const RenderPort*>(*succIt);
            Processor* predecessorProc = predecessorPort->getProcessor();
            if (!predecessorPort) {
                LWARNING("RenderInport connected to non-RenderOutport");
                continue;
            }
            tgtAssert(predecessorPort->isOutport(), "inport not connected to an outport");
            if (predecessorPort->getRenderSizePropagation() == RenderPort::RENDERSIZE_RECEIVER &&
                std::find(subNetwork.begin(), subNetwork.end(), predecessorProc) != subNetwork.end()  ) {
                    // outport is size receiver and owning processor is to be linked
                    // => save outport as corresponding size receiver and stop traversal
                    result.push_back(const_cast<RenderPort*>(predecessorPort));
            }
            else {
                // outport is not a size receiver port or owning processor is not to be size linked
                // => add owning processor's render inports to worklist, unless they are size origins themselves
                std::vector<Port*> predInports = predecessorProc->getInports();
                for (size_t s=0; s<predInports.size(); s++) {
                    RenderPort* predRenderPort = dynamic_cast<RenderPort*>(predInports.at(s));
                    if (predRenderPort && !processed.count(predRenderPort)) {
                        if (predRenderPort->getRenderSizePropagation() == RenderPort::RENDERSIZE_DEFAULT ||
                                std::find(subNetwork.begin(), subNetwork.end(), predecessorProc) == subNetwork.end())
                            worklist.push(predRenderPort);
                    }
                }
            }
        }
    }

    return result;

}

std::vector<RenderPort*> ProcessorNetwork::getSuccessingRenderSizeOrigins(RenderPort* rootPort, const std::vector<Processor*>& subNetwork) const {
    tgtAssert(rootPort, "null pointer passed");

    std::vector<RenderPort*> result;

    std::queue<RenderPort*> worklist;
    if (rootPort->isOutport())
        worklist.push(rootPort);
    else {
        std::vector<Port*> outports = rootPort->getProcessor()->getOutports();
        for (size_t i=0; i<outports.size(); i++) {
            RenderPort* rp = dynamic_cast<RenderPort*>(outports.at(i));
            if (rp)
                worklist.push(rp);
        }
    }

    std::set<RenderPort*> processed; //< for cycle prevention
    while (!worklist.empty()) {
        RenderPort* curPort = worklist.front();
        worklist.pop();
        tgtAssert(curPort, "null pointer in worklist");
        tgtAssert(curPort->isOutport(), "inport in worklist");
        processed.insert(curPort);
        if (curPort->isLoopPort())  //< do not traverse over loop ports
            continue;
        const std::vector<const Port*> successors = curPort->getConnected();
        for (std::vector<const Port*>::const_iterator succIt=successors.begin(); succIt != successors.end(); ++succIt) {
            const RenderPort* successorPort = dynamic_cast<const RenderPort*>(*succIt);
            Processor* successorProc = successorPort->getProcessor();
            if (!successorPort) {
                LWARNING("RenderOutport connected to non-RenderInport");
                continue;
            }
            tgtAssert(successorPort->isInport(), "outport not connected to an inport");
            if (successorPort->getRenderSizePropagation() == RenderPort::RENDERSIZE_ORIGIN &&
                std::find(subNetwork.begin(), subNetwork.end(), successorProc) != subNetwork.end()  ) {
                    // inport is size origin and owning processor is to be linked
                    // => save inport as corresponding size origin for root port and stop traversal
                    result.push_back(const_cast<RenderPort*>(successorPort));
            }
            else {
                // inport is not a size origin or owning processor is not to be size linked
                // => add owning processor's render outports to worklist, unless they are size receivers themselves
                std::vector<Port*> succOutports = successorProc->getOutports();
                for (size_t s=0; s<succOutports.size(); s++) {
                    RenderPort* succRenderPort = dynamic_cast<RenderPort*>(succOutports.at(s));
                    if (succRenderPort && !processed.count(succRenderPort)) {
                        if (succRenderPort->getRenderSizePropagation() == RenderPort::RENDERSIZE_DEFAULT ||
                                std::find(subNetwork.begin(), subNetwork.end(), successorProc) == subNetwork.end())
                            worklist.push(succRenderPort);
                    }
                }
            }
        }
    }

    return result;
}

std::map<RenderPort*, std::vector<RenderPort*> > ProcessorNetwork::getRenderSizeReceiverToOriginsMap(
    const std::vector<Processor*> subNetwork) const
{
    std::map<RenderPort*, std::vector<RenderPort*> > result;
    for (size_t i=0; i<subNetwork.size(); i++) {
        Processor* proc = subNetwork.at(i);
        std::vector<Port*> inports = proc->getInports();
        for (size_t j=0; j<inports.size(); j++) {
            RenderPort* rp = dynamic_cast<RenderPort*>(inports.at(j));
            if (rp && rp->getRenderSizePropagation() == RenderPort::RENDERSIZE_ORIGIN) {
                tgtAssert(rp->getSizeOriginProperty(), "Render size origin without origin property");
                std::vector<PropertyLink*> links = rp->getSizeOriginProperty()->getLinks();
                for (size_t k=0; k<links.size(); k++) {
                    if (dynamic_cast<LinkEvaluatorRenderSize*>(links.at(k)->getLinkEvaluator())) {
                        Property* destProp = links.at(k)->getDestinationProperty();
                        RenderPort* destPort = dynamic_cast<RenderPort*>(destProp->getOwner());
                        if (destPort && destPort->getRenderSizePropagation() == RenderPort::RENDERSIZE_RECEIVER) {
                            if (!result.count(destPort))
                                result.insert(std::make_pair(destPort, std::vector<RenderPort*>()));
                            result[destPort].push_back(rp);
                        }
                    }
                }
            }
        }
    }

    return result;
}

// observer methods

void ProcessorNetwork::notifyNetworkChanged() const {
    vector<ProcessorNetworkObserver*> observers = getObservers();
    for (size_t i=0; i<observers.size(); ++i)
        observers[i]->networkChanged();
}

void ProcessorNetwork::notifyProcessorAdded(const Processor* processor) const {
    vector<ProcessorNetworkObserver*> observers = getObservers();
    for (size_t i=0; i<observers.size(); ++i)
        observers[i]->processorAdded(processor);
}


void ProcessorNetwork::notifyProcessorRemoved(const Processor* processor) const {
    vector<ProcessorNetworkObserver*> observers = getObservers();
    for (size_t i=0; i<observers.size(); ++i)
        observers[i]->processorRemoved(processor);
}

void ProcessorNetwork::notifyProcessorRenamed(const Processor* processor, const std::string& prevName) const {
    vector<ProcessorNetworkObserver*> observers = getObservers();
    for (size_t i=0; i<observers.size(); ++i)
        observers[i]->processorRenamed(processor, prevName);
}

void ProcessorNetwork::notifyPortConnectionAdded(const Port* outport, const Port* inport) const {
    vector<ProcessorNetworkObserver*> observers = getObservers();
    for (size_t i=0; i<observers.size(); ++i)
        observers[i]->portConnectionAdded(outport, inport);
}

void ProcessorNetwork::notifyPortConnectionRemoved(const Port* outport, const Port* inport) const {
    vector<ProcessorNetworkObserver*> observers = getObservers();
    for (size_t i=0; i<observers.size(); ++i)
        observers[i]->portConnectionRemoved(outport, inport);
}

void ProcessorNetwork::notifyPropertyLinkAdded(const PropertyLink* link) const {
    vector<ProcessorNetworkObserver*> observers = getObservers();
    for (size_t i=0; i<observers.size(); ++i)
        observers[i]->propertyLinkAdded(link);
}

void ProcessorNetwork::notifyPropertyLinkRemoved(const PropertyLink* link) const {
    vector<ProcessorNetworkObserver*> observers = getObservers();
    for (size_t i=0; i<observers.size(); ++i)
        observers[i]->propertyLinkRemoved(link);
}

std::vector<Property*> ProcessorNetwork::getPropertiesByID(const std::string& id) const {
    std::vector<Property*> result;
    for (size_t i = 0; i < processors_.size(); ++i) {
        if (processors_[i]->getProperty(id))
            result.push_back(processors_[i]->getProperty(id));
    }
    return result;
}

void ProcessorNetwork::preparePropertyRemoval(Property* property) {
    removePropertyLinksFromProperty(property);
}

void ProcessorNetwork::propertyValueHasBeenModified(Property* prop) {
    if(workspace_)
        workspace_->setModified(true);
}

void ProcessorNetwork::setWorkspace(Workspace* workspace) {
    workspace_ = workspace;
}

Workspace* ProcessorNetwork::getWorkspace() const {
    return workspace_;
}

} // namespace voreen
