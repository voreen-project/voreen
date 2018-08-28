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

#ifndef VRN_PROCESSORNETWORK_H
#define VRN_PROCESSORNETWORK_H

#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/utils/exception.h"
#include "voreen/core/network/processornetworkobserver.h"
#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/link/propertylink.h"
#include "voreen/core/properties/link/linkevaluatorbase.h"

#include <sstream>

namespace voreen {

class Workspace;
class NetworkSerializer;
class Property;
class LinkEvaluatorBase;
class Port;
class RenderPort;

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API Observable<ProcessorNetworkObserver>;
#endif
/**
 * Manages a network of processors.
 */
class VRN_CORE_API ProcessorNetwork : public Serializable, public ProcessorObserver, public Observable<ProcessorNetworkObserver> {
public:

    /// Current Network version. Used to check for incompatible versions.
    static const int NETWORK_VERSION ;
    /// All imported networks whose version is lower than defined by this constant will emit a warning.
    static const int WARNING_VERSION;


    ProcessorNetwork();
    ~ProcessorNetwork();

    /**
     * Adds a processor to the network and assigns the given name to it.
     * The network takes ownership of the added processor and deletes it on destruction.
     *
     * @param processor the processor to be added to the network.
     *        Must not be already contained by the network.
     * @param name the name to assign to the processor. If no name is specified,
     *        the processor's class name is used.
     *
     * @note If the passed name is already assigned to a processor in the network,
     *       a unique name based on the passed one is generated. The same applies to
     *       the processor's class name, if no name is specified.
     */
    void addProcessor(Processor* processor, const std::string& name = "");

    /**
     * Removes a processor from the network and deletes it. Also
     * removes the processor's port connections and property links.
     *
     * @param processor the processor to be removed from the network.
     *        Must be contained by the network.
     */
    void removeProcessor(Processor* processor);

    /**
     * Removes all processors, port connections and property links
     * from the network and deletes them.
     */
    void clear();

    /**
     * Returns the number of processor of the network.
     */
    size_t numProcessors() const;

    /**
     * Returns true, if the network does not contain any processors.
     */
    bool empty() const;

    /**
     * Returns the network's processors.
     */
    const std::vector<Processor*>& getProcessors() const;

    /**
     * Merges a ProcessorNetwork into this.
     * Adds Processors, Connections and PropertyLinks.
     *
     * @param pnw the network to add. It will be deleted afterwards.
     */
    void mergeSubNetwork(ProcessorNetwork* pnw);

     /**
      * Clones a list of Processors and returns a ProcessorNetwork including these Processors.
      * Connections and Links between these Processors are cloned as well.
      *
      * @param vec vector of processors to bulid a subnetwork. The Processors must be part of this.
      */
    ProcessorNetwork* cloneSubNetwork(std::vector<Processor*> vec);

    /**
     * Returns all processors in the network matching
     * the given template type.
     */
    template<class T>
    std::vector<T*> getProcessorsByType() const;

    /**
     * Returns the processor with the given name, or null if the network does not contain a
     * processor with that name.
     *
     * @note The name is unique among all Processor instances in the network.
     *
     * @see setProcessorName
     */
    Processor* getProcessor(const std::string& name) const;

    /**
     * Returns the processor with the given name, or null if the network does not contain a
     * processor with that name.
     *
     * @deprecated Use getProcessor() instead.
     * @note The name is unique among all Processor instances in the network.
     *
     * @see setProcessorName
     */
    Processor* getProcessorByName(const std::string& name) const;

    /**
     * Returns whether the passed processor is an element of the network.
     */
    bool contains(const Processor* processor) const;

    /**
     * Assigns the passed name to the passed processor.
     *
     *
     * @param processor the processor to rename. Must be part of the network.
     * @param name the name to assign. Must not be empty and must not be
     *      assigned to any processor in the network.
     *
     * @throw VoreenException if the passed processor is not contained by the
     *        network, or if the passed name is empty or already assigned
     *        to a processor in the network.
     *
     * @see getProcessor
     * @see generateUniqueProcessorName
     */
    void setProcessorName(Processor* processor, const std::string& name) const;

    /**
     * Creates a connection between the passed ports.
     *
     * @param inport the inport to be connected.
     *      Must be owned by a processor of the network.
     * @param outport the outport to be connected.
     *      Must be owned by a processor of the network.
     *
     * @return true, if the connection was created,
     *      false if the passed ports are not compatible
     */
    bool connectPorts(Port* outport, Port* inport);

    /// @overload
    bool connectPorts(Port& outport, Port& inport);

    /**
     * Deletes the connection between the passed ports.
     *
     * @param inport the inport to be disconnected.
     *      Must be owned by a processor of the network.
     * @param outport the outport to be disconnected
     *      Must be owned by a processor of the network.
     */
    void disconnectPorts(Port* outport, Port* inport);

    /// @overload
    void disconnectPorts(Port& outport, Port& inport);

    /**
     * Removes all incoming and outgoing port connections
     * from the passed processor.
     *
     * @param processor the processor to disconnect.
     *      Must be part of the network.
     */
    void disconnectPorts(Processor* processor);

    /**
     * Returns all processor properties in the network
     * with the specified id.
     */
    std::vector<Property*> getPropertiesByID(const std::string& id) const;

    /**
     * Returns all processor properties in the network matching
     * the given template type.
     */
    template<class T>
    std::vector<T*> getPropertiesByType() const;

    /**
     * Returns the network's property links.
     */
    const std::vector<PropertyLink*>& getPropertyLinks() const;

    /**
     * Returns true, is the passed property object is
     * the source or destination of any property link.
     */
    bool isPropertyLinked(const Property* property) const;

    /**
     * Tests whether the ProcessorNetwork contains a link with between the two passed Properties with
     * the explicit LinkEvaluatorBase
     *
     * @param src The source property of the link
     * @param dest The destination property of the link
     * @param evaluator The evaluator of the link
     */
    bool containsPropertyLink(const Property* src, const Property* dest, LinkEvaluatorBase* evaluator) const;

    /**
     * Creates a property link, if possible and no link already exists from the source to the destination,
     * and adds it to the network.
     *
     * @note EventProperties cannot be linked.
     *
     * @param src the source property of the link. Must be owned by processor of this network.
     * @param dest the destination property of the link. Must be owned by processor of this network
     *      and must not equal src.
     * @param linkEvaluator the link evaluator defining the type of the link.
     *      If not passed, an identity link is created. Will be deleted by the PropertyLink.
     *
     * @return the created property link or null, if the link could not be created.
     */
    PropertyLink* createPropertyLink(Property* src, Property* dest, LinkEvaluatorBase* linkEvaluator = 0);

    /**
     * Creates incoming and outgoing links from properties of the passed processor
     * to all other properties in the network. Property pairs that are already linked
     * are ignored.
     *
     * @tparam T Type of properties to consider, including subtypes.
     * @param processor Processor whose properties are to be linked.
     * @param evaluator Evaluator to use for creating the links. Is deleted by the function.
     *        If the null pointer is passed, the first compatible evaluator will be used.
     * @return The number of created property links.
     */
    template<class T>
    int createPropertyLinksForProcessor(Processor* processor, LinkEvaluatorBase* evaluator = 0);

    /**
     * Creates links between the passed properties.
     *
     * @param properties the properties to be linked. Must be owned by processors
     *      that are part of the network.
     * @param linkEvaluator the link evaluator defining the type of the link.
     *      If not passed, an identity link is created.
     * @param replace if set to true, already existing links between properties to link
     *      are replaced. Otherwise, the existing links are kept.
     * @param transitive defines whether transitive linking is sufficient or
     *      if each pair of properties has to be linked directly. If set to true,
     *      a cyclic chain of links connecting all of the passed properties is created.
     *      Transitive linking is highly recommended, since the number of links to create is
     *      proportional to the number of properties, whereas pair-wise linking has quadratic complexity!
     *
     * @return the number of links that have been created
     */
    int createPropertyLinksBetweenProperties(const std::vector<Property*>& properties, LinkEvaluatorBase* linkEvaluator = 0,
        bool replace = false, bool transitive = true);

    /**
     * Creates links between the properties of the passed sub network,
     * which are additionally constricted by several filter criteria.
     *
     * @tparam T Type of the properties to be linked, including subtypes.
     *      Use \c Property to embrace all types.
     * @param processors processors whose properties are to be considered.
     *      Must be part of the network. Pass an empty vector for embracing all processors in the network.
     * @param propertyIDs IDs of the properties to link.
     *      Pass an empty vector for embracing all properties.
     * @param linkEvaluator the link evaluator defining the type of the link.
     *      If not passed, an identity link is created.
     * @param replace if set to true, already existing links between properties to link
     *      are replaced. Otherwise, the existing links are kept.
     * @param transitive defines whether transitive linking is sufficient or
     *      if each pair of properties has to be linked directly. If set to true,
     *      a cyclic chain of links connecting all of the matching properties is created.
     *      Transitive linking is highly recommended, since the number of links to create is
     *      proportional to the number of properties, whereas pair-wise linking has quadratic complexity!
     *
     * @return the number of links that have been created
     */
    template<class T>
    int createPropertyLinksWithinSubNetwork(const std::vector<Processor*>& processors, const std::vector<std::string>& propertyIDs,
        LinkEvaluatorBase* linkEvaluator = 0, bool replace = false, bool transitive = true);

    /**
     * Removes a property link from the network and deletes it.
     */
    void removePropertyLink(PropertyLink* propertyLink);

    /**
     * Removes the property link between the passed properties, if such link exists, and deletes it.
     *
     * @return true, if a link between the properties existed and has been removed
     */
    bool removePropertyLink(Property* src, Property* dst);

    /**
     * Removes all property links from the network and deletes them.
     *
     * @return Number of removed links.
     */
    int removeAllPropertyLinks();

    /**
     * Removes all incoming and outgoing links from the passed property and deletes them.
     *
     * @param property the property to unlink. Must be owned by a processor
     *      that is part of the network.
     * @return Number of removed property links.
     */
    int removePropertyLinksFromProperty(Property* property);

    /**
     * Removes all incoming and outgoing property links
     * from the properties of the passed processor and deletes them.
     *
     * @param processor processor to unlink. Must be an element of the network.
     * @return Number of removed property links.
     */
    int removePropertyLinksFromProcessor(Processor* processor);

    /**
     * Removes all incoming and outgoing property links from the properties
     * that are owned by the passed processor and whose type matches the
     * specified template parameter (including subtypes).
     *
     * @tparam T property type (including subtypes) whose links are to be removed.
     * @param processor processor to unlink. Must be an element of the network.
     * @return Number of removed property links.
     */
    template<class T>
    int removePropertyLinksFromProcessor(Processor* processor);

    /**
     * Removes and deletes all links between properties that are owned by
     * the passed processors.
     *
     * @param processors Sub-network whose properties should be unlinked. Pass an empty
     *        empty vector for embracing the entire network.
     * @return Number of removed property links.
     */
    int removePropertyLinksFromSubNetwork(const std::vector<Processor*>& processors);

    /**
     * Removes and deleted all property links between properties that are owned by
     * the passed processors and whose types matches the specified template parameter.
     *
     * @tparam T Property types whose links should be deleted (including subtypes).
     * @param processors Sub-network whose properties should be unlinked. Pass an empty
     *        empty vector for embracing the entire network.
     * @return Number of removed property links.
     */
    template<class T>
    int removePropertyLinksFromSubNetwork(const std::vector<Processor*>& processors);

    /**
     * Creates a render size link from \p source to \p destination, which must both be
     * either a render size origin or a render size receiver.
     * Usually, a size link exists between a render size origin inport and a render size receiving outport.
     * A render size link from a size receiver to a size origin is not possible.
     *
     * @return The created render size link.
     */
    PropertyLink* createRenderSizeLink(RenderPort* source, RenderPort* destination);

    /// @overload
    PropertyLink* createRenderSizeLink(RenderPort& source, RenderPort& destination);

    /**
     * Removes the render size link from \p source to \p dest, if such a link exists.
     *
     * @return true, if a link has been removed.
     */
    bool removeRenderSizeLink(RenderPort* source, RenderPort* destination);

    /// @overload
    bool removeRenderSizeLink(RenderPort& source, RenderPort& destination);

    /**
     * Returns the render size link from \p source to \p destination,
     * or null if no such link exists.
     */
    PropertyLink* getRenderSizeLink(RenderPort* source, RenderPort* destination);

    /// @overload
    PropertyLink* getRenderSizeLink(RenderPort& source, RenderPort& destination);

    /**
     * Creates render size links between the passed processors' RenderPorts
     * according to the current network topology.
     *
     * @param processors The subnetwork to consider
     * @param replaceExisting if true, existing render size links are replaced
     *
     * @return The number of created render size links.
     */
    int createRenderSizeLinksWithinSubNetwork(const std::vector<Processor*>& processors, bool replaceExisting = false);

    /**
     * Creates render size links that span the port connection from \p outport,
     * which must be an outport, to \p inport, which must be an inport.
     *
     * @param replaceExisting if true, existing render size links are replaced
     *
     * @return The number of created render size links.
     */
    int createRenderSizeLinksOverConnection(RenderPort* outport, RenderPort* inport, bool replaceExisting = false);

    /**
     * Creates render size links from/to the passed processor's RenderPorts,
     * by considering the network topology.
     *
     * @param replaceExisting if true, existing render size links are replaced
     *
     * @return The number of created render size links.
     */
    int createRenderSizeLinksForProcessor(Processor* processor, bool replaceExisting = false);

    /**
     * Removes all render size links between the passed processors' RenderPorts.
     *
     * @return The number of removed render size links.
     */
    int removeRenderSizeLinksFromSubNetwork(const std::vector<Processor*>& processors);

    /**
     * Removes all render size links that span the port connection from \p outport,
     * which must be an outport, to \p inport, which must be an inport.
     *
     * @note A actual connection between the passed ports does not necessarily have to exist.
     *
     * @return The number of created render size links.
     */
    int removeRenderSizeLinksOverConnection(RenderPort* outport, RenderPort* inport);

    /**
     * Removes all incoming and outgoing render size links from/to the
     * passed processor's RenderPorts.
     *
     * @return The number of removed render size links.
     */
    int removeRenderSizeLinksFromProcessor(Processor* processor);

    /**
     * Set the version of the network file.
     */
    void setVersion(int version);

    /**
     * Returns the version of the network file.
     */
    int getVersion() const;

    /**
     * Generates a processor name that is unique among all processors in
     * the network by appending a number to the passed name.
     *
     * If the passed name is already unique, it is returned without
     * modifications.
     */
    std::string generateUniqueProcessorName(const std::string& name) const;

    /**
     * Convenience function: Inserts a processor into the network
     * at the connection between two other processors.
     *
     * The connection between those ports (A,B) will be terminated and will be
     * routed through processor C.<br>
     * A->B  will be A->C->B afterwards.<br>
     *
     * @param outport The outport of the predecessing processor
     * @param inport The inport of the successing processor
     * @param processor The new processor which is to be inserted between outport and inport
     */
    void addProcessorInConnection(Port* outport, Port* inport, Processor* processor);

    /**
    * Convenience function: Inserts a processor into the network
    * which takes over the connections of the processor being replaced.
    *
    * The connection between those ports (A,B1,B2,B3,..) will be terminated and will be
    * routed through processor C.<br>
    * A->(B1, B2, B3,..) will be A->C->(B1, B2, B3,..) afterwards.<br>
    *
    * @param outport The outport of the predecessing processor
    * @param processor The new processor which is to be inserted between outport and all connected inports
    * @return true if successful, false otherwise
    */
    void takeOverConnections(Port* outport, Processor* processor);

    /**
    * Convenience function: Tests whether a processor can take over at least one connection
    * of the processor being specified by outport.
    * This even 'allows' to connect a processor directly to a port for convinience reasons.
    *
    * Does not manipulate the network!
    *
    * @see ProcessorNetwork::takeOverConnections
    *
    * @param outport The outport of the predecessing processor
    * @param processor The new processor which is to be inserted between outport and all connected inports
    * @return true if the processor can take over connections of the processor being specified by the port, false otherwise
    */
    bool canTakeOverConnections(Port* outport, Processor* processor);

    /**
     * Convenience function: Replaces a given processor in the network by another.
     *
     * If possible, connections and links from other processors to the old one are redirected
     * to the new processor. Incompatible connections are terminated.
     * If there is no compatible connection, this operation leads to a new isolated processor
     * in the network.
     *
     * @note The replaced processor is deleted after the operation.
     *
     * @param oldProc the old processor which will be replaced, must not be null,
     *  gets deleted.
     * @param newProc the new processor which should be entered, must not be null
     */
    void replaceProcessor(Processor* oldProc, Processor* newProc);

    /**
     * Returns the meta data container of this processor.
     * External objects, such as GUI widgets, can use it
     * to store and retrieve persistent meta data without
     * having to bother with the serialization themselves.
     *
     * @see MetaDataContainer
     */
    MetaDataContainer& getMetaDataContainer() const;

    /**
     * @see Serializable::serialize
     */
    void serialize(Serializer& s) const;

    /**
     * @see Serializable::deserialize
     */
    void deserialize(Deserializer& s);

    /**
     * Serialization helper function: Returns serialization error messages.
     */
    std::vector<std::string> getErrors() const;

    /**
     * Serialization helper function: Sets serialization error messages.
     */
    void setErrors(const std::vector<std::string>& errorList);

    /**
     * This callback is called before a property is removed. It takes care, that all links
     * for this property are removed.
     *
     * \param property The property which shall be prepared for removal
     */
    void preparePropertyRemoval(Property* property);

    /**
     * Sets the modified flag in the workspace (if assigned)
     */
    void propertyValueHasBeenModified(Property* prop);

    /**
     * Get a list of all the referenced files.
     */
    const std::vector<std::string>& getReferencedFiles() const;

    /// Assigns the workspace that owns the network (optional).
    void setWorkspace(Workspace* workspace);

    /// Returns the workspace that owns the network. Returns 0, if workspace has not been set.
    Workspace* getWorkspace() const;

private:
    /**
     * Returns the processor that directly or indirectly owns the passed property.
     * Indirect ownership means that the processor owns another property owner,
     * e.g. a port, that owns the property.
     */
    Processor* getOwningProcessor(Property* property) const;

    /**
     * Returns all RenderPort size receivers that are direct predecessors of \p renderPort,
     * i.e., all size receivers that are connected to \c renderPort through a path, which does
     * not contain any other size receivers. Ports not owned by a processor of \p subNetwork are ignored.
     */
    std::vector<RenderPort*> getPredecessingRenderSizeReceivers(RenderPort* renderPort, const std::vector<Processor*>& subNetwork) const;

    /**
     * Returns all RenderPort size origins that are direct successors of \p renderPort,
     * i.e., all size origins that are connected to \c renderPort through a path, which does
     * not contain any other size origin. Ports not owned by a processor of \p subNetwork are ignored.
     */
    std::vector<RenderPort*> getSuccessingRenderSizeOrigins(RenderPort* renderPort, const std::vector<Processor*>& subNetwork) const;

    /// Returns a map from render inports to the render outports they are linked with, within the passed sub network.
    std::map<RenderPort*, std::vector<RenderPort*> > getRenderSizeReceiverToOriginsMap(const std::vector<Processor*> subNetwork) const;

    /// Calls networkChanged() on the registered observers.
    void notifyNetworkChanged() const;

    /// Calls processorAdded() on the registered observers.
    void notifyProcessorAdded(const Processor* processor) const;

    /// Calls processorRemoved() on the registered observers.
    void notifyProcessorRemoved(const Processor* processor) const;

    /// Calls processorRenamed() on the registered observers.
    void notifyProcessorRenamed(const Processor* processor, const std::string& prevName) const;

    /// Calls portConnectionAdded on the registered observers.
    void notifyPortConnectionAdded(const Port* source, const Port* destination) const;

    /// Calls portConnectionRemoved on the registered observers.
    void notifyPortConnectionRemoved(const Port* source, const Port* destination) const;

    /// Calls propertyLinkAdded() on the registered observers.
    void notifyPropertyLinkAdded(const PropertyLink* link) const;

    /// Calls propertyLinkRemoved() on the registered observers.
    void notifyPropertyLinkRemoved(const PropertyLink* link) const;

    /// List of all \sa Processor contained in the network
    std::vector<Processor*> processors_;

    /// List of all \sa PropertyList contained in the network
    std::vector<PropertyLink*> propertyLinks_;

    /// The network version
    int version_;
    std::vector<std::string> errorList_;

    /// The owning workspace. Is not used internally and may be null.
    Workspace* workspace_;

    /**
     * Contains the associated meta data.
     *
     * We want to return a non-const reference to it from a const member function
     * and since the MetaDataContainer does not affect the processor itself,
     * mutable appears justifiable.
     */
    mutable MetaDataContainer metaDataContainer_;

    static const std::string loggerCat_; ///< category used in logging
};

//---------------------------------------------------------------------------
// template definitions

template<class T>
std::vector<T*> ProcessorNetwork::getProcessorsByType() const {
    std::vector<T*> result;
    for (size_t i = 0; i < processors_.size(); ++i) {
        T* proc = dynamic_cast<T*>(processors_[i]);
        if (proc)
            result.push_back(proc);
    }
    return result;
}

template<class T>
std::vector<T*> ProcessorNetwork::getPropertiesByType() const {
    std::vector<T*> result;
    for (size_t i = 0; i < processors_.size(); ++i) {
        std::vector<T*> processorProps = processors_[i]->getPropertiesByType<T>();
        result.insert(result.end(), processorProps.begin(), processorProps.end());
    }
    return result;
}

template<class T>
int ProcessorNetwork::createPropertyLinksWithinSubNetwork(const std::vector<Processor*>& processors,
    const std::vector<std::string>& propertyIDs,
    LinkEvaluatorBase* linkEvaluator,
    bool replace, bool transitive) {

    // determine processors to be linked
    std::vector<Processor*> linkProcessors = (processors.empty() ? getProcessors() : processors);

    // collect properties to be linked
    std::vector<Property*> linkProps;
    for (size_t i=0; i<linkProcessors.size(); ++i) {
        Processor* curProcessor = linkProcessors[i];
        if (!contains(curProcessor)) {
            LWARNING("Processor is not part of the network: " << curProcessor->getID());
            continue;
        }

        if (propertyIDs.empty()) {
            // no property ids specified => select all properties of matching type
            std::vector<T*> procProps = curProcessor->getPropertiesByType<T>();
            linkProps.insert(linkProps.end(), procProps.begin(), procProps.end());
        }
        else {
            // property ids specified => query each processor for property by id and check type
            for (size_t j=0; j<propertyIDs.size(); ++j) {
                if (dynamic_cast<T*>(curProcessor->getProperty(propertyIDs[j])))
                    linkProps.push_back(curProcessor->getProperty(propertyIDs[j]));
            }
        }
    }

    // passed collected properties to helper function
    return createPropertyLinksBetweenProperties(linkProps, linkEvaluator, replace, transitive);
}

template<class T>
int ProcessorNetwork::createPropertyLinksForProcessor(Processor* processor, LinkEvaluatorBase* evaluator) {
    tgtAssert(processor, "null pointer passed");
    tgtAssert(contains(processor), "passed processor is not part of the network");

    // collect matching properties of passed processor
    std::vector<T*> procProperties = processor->getPropertiesByType<T>();

    // get all matching properties in the network
    std::vector<T*> networkProperties = getPropertiesByType<T>();

    // create links from processor properties to network properties
    int numCreated = 0;
    for (size_t i=0; i<procProperties.size(); i++) {
        T* procProperty = procProperties.at(i);
        tgtAssert(procProperty, "null pointer");
        for (size_t j=0; j<networkProperties.size(); j++) {
            T* networkProperty = networkProperties.at(j);
            tgtAssert(networkProperty, "null pointer");
            if (procProperty == networkProperty)
                continue;
            if (!procProperty->isLinkedWith(networkProperty))
                if (createPropertyLink(procProperty, networkProperty, evaluator ? dynamic_cast<LinkEvaluatorBase*>(evaluator->create()) : 0))
                    numCreated++;
            if (!networkProperty->isLinkedWith(procProperty))
                if (createPropertyLink(networkProperty, procProperty, evaluator ? dynamic_cast<LinkEvaluatorBase*>(evaluator->create()) : 0))
                    numCreated++;
        }
    }

    delete evaluator;
    return numCreated;
}

template<class T>
int ProcessorNetwork::removePropertyLinksFromProcessor(Processor* processor) {
    tgtAssert(processor, "null pointer passed");
    tgtAssert(contains(processor), "passed processor is not part of the network");

    int numRemoved = 0;
    for (size_t i = 0; i < propertyLinks_.size(); ++i) {
        PropertyLink* link = propertyLinks_[i];
        if ( (getOwningProcessor(link->getSourceProperty()) == processor && dynamic_cast<T*>(link->getSourceProperty())) ||
             (getOwningProcessor(link->getDestinationProperty()) == processor && dynamic_cast<T*>(link->getDestinationProperty())) ) {
                removePropertyLink(link);
                i--; // removePropertyLink also deletes the link object
            numRemoved++;
        }
    }
    return numRemoved;
}

template<class T>
int ProcessorNetwork::removePropertyLinksFromSubNetwork(const std::vector<Processor*>& processors) {
    // determine processors whose links are to be removed
    std::vector<Processor*> linkProcessors = (processors.empty() ? getProcessors() : processors);

    // collect properties to unlink
    std::set<Property*> unlinkProps;
    for (size_t i=0; i<linkProcessors.size(); ++i) {
        Processor* curProcessor = linkProcessors[i];
        if (!contains(curProcessor)) {
            LWARNING("Processor is not part of the network: " << curProcessor->getID());
            continue;
        }

        std::vector<T*> procProps = curProcessor->getPropertiesByType<T>();
        unlinkProps.insert(procProps.begin(), procProps.end());
    }

    // delete all links between any pair of unlinkProps
    int numRemoved = 0;
    std::vector<PropertyLink*> links = std::vector<PropertyLink*>(getPropertyLinks());
    for (size_t i=0; i<links.size(); i++) {
        PropertyLink* link = links.at(i);
        if (unlinkProps.count(link->getSourceProperty()) && unlinkProps.count(link->getDestinationProperty())) {
            removePropertyLink(link);
            numRemoved++;
        }
    }
    return numRemoved;
}

} //namespace voreen

#endif
