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

#ifndef VRN_PROCESSORNETWORKOBSERVER_H
#define VRN_PROCESSORNETWORKOBSERVER_H

#include "voreen/core/utils/observer.h"

namespace voreen {

class Port;
class Processor;
class PropertyLink;

/**
 * Interface for network observers. Objects of this type
 * can be registered at a ProcessorNetwork.
 */
class VRN_CORE_API ProcessorNetworkObserver : public Observer {
public:

    /**
     * This method is called on all modifications of the observed network
     * that do not match one of the more specialized modifications below.
     *
     * @note This function is also called by the standard implementations
     *  of the specialized notifier functions below. Therefore, if you do not need
     *  to differentiate between the single modification actions, it is sufficient
     *  to implement just this function.
     */
    virtual void networkChanged() = 0;

    /**
     * This method is called by the observed network after a processor has been added.
     *
     * @note The standard implementation of this function calls \a networkChanged.
     *  Override it for a custom reaction to this particular event.
     *
     * @param processor the processor that has been added
     */
    virtual void processorAdded(const Processor* processor);

    /**
     * This method is called by the observed network immediately \e before
     * a processor is removed from the network and deleted.
     *
     * @note The standard implementation of this function calls \a networkChanged.
     *  Override it for a custom reaction to this particular event.
     *
     * @param processor the processor that is to be removed
     */
    virtual void processorRemoved(const Processor* processor);

    /**
     * This method is called by the observed network after
     * a processor has been renamed.
     *
     * @note The standard implementation of this function calls \a networkChanged.
     *  Override it for a custom reaction to this particular event.
     *
     * @param processor the processor that has been renamed
     * @param prevName the previous name of the processor
     */
    virtual void processorRenamed(const Processor* processor, const std::string& prevName);

    /**
     * This method is called by the observed network when a port connection has been added
     *
     * @note The standard implementation of this function calls \a networkChanged.
     *  Override it for a custom reaction to this particular event.
     *
     */
    virtual void portConnectionAdded(const Port* outport, const Port* inport);

    /**
     * This method is called by the observed network when a port connection has been removed
     *
     * @note The standard implementation of this function calls \a networkChanged.
     *  Override it for a custom reaction to this particular event.
     *
     */
    virtual void portConnectionRemoved(const Port* outport, const Port* inport);

    /**
     * This method is called by the observed network after
     * a PropertyLink has been added.
     *
     * @note The standard implementation of this function calls \a networkChanged.
     *  Override it for a custom reaction to this particular event.
     *
     * @param link the PropertyLink that has been added
     */
    virtual void propertyLinkAdded(const PropertyLink* link);

    /**
     * This method is called by the observed network immediately \e before
     * a PropertyLink is removed from the network and deleted.
     *
     * @note The standard implementation of this function calls \a networkChanged.
     *  Override it for a custom reaction to this particular event.
     *
     * @param link the link that is to be removed
     */
    virtual void propertyLinkRemoved(const PropertyLink* link);
};

} // namespace

#endif // VRN_PROCESSORNETWORKOBSERVER_H
