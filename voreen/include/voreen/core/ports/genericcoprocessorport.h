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

#ifndef VRN_GENERICCOPROCESSORPORT_H
#define VRN_GENERICCOPROCESSORPORT_H

#include "voreen/core/ports/coprocessorport.h"

namespace voreen {

/**
 * @brief
 */
template <typename T>
class GenericCoProcessorPort : public CoProcessorPort  {
public:
    GenericCoProcessorPort(PortDirection direction, const std::string& id, const std::string& guiName = "", bool allowMultipleConnections = false,
                           Processor::InvalidationLevel invalidationLevel = Processor::INVALID_RESULT);

    virtual Port* create(PortDirection direction, const std::string& id, const std::string& guiName = "") const;
    virtual std::string getClassName() const;

    std::vector<T*> getConnectedProcessors() const;

    ///Get the first connected Processor
    T* getConnectedProcessor() const;
protected:
    virtual void setProcessor(Processor* p);
};

// ---------------------------------------- implementation ----------------------------------------

template <typename T>
GenericCoProcessorPort<T>::GenericCoProcessorPort(PortDirection direction, const std::string& id, const std::string& guiName, bool allowMultipleConnections,
                           Processor::InvalidationLevel invalidationLevel)
    : CoProcessorPort(direction, id, guiName, allowMultipleConnections, invalidationLevel)
{
}


template <typename T>
Port* GenericCoProcessorPort<T>::create(PortDirection direction, const std::string& id, const std::string& guiName) const{
    return new GenericCoProcessorPort<T>(direction,id,guiName);
}

template <typename T>
std::string GenericCoProcessorPort<T>::getClassName() const{
    return "GenericCoProcessorPort";
}

template <typename T>
void GenericCoProcessorPort<T>::setProcessor(Processor* p) {
    Port::setProcessor(p);
    if (isOutport()) {
        T* tp = dynamic_cast<T*>(p);
        tgtAssert(tp, "CoProcessor outport attached to processor of wrong type");
        if (!tp)
            LERRORC("voreen.coprocessorport", "CoProcessor outport attached to processor of wrong type" << p->getID() << "|" << this->getID());
    }
}

template <typename T>
std::vector<T*> GenericCoProcessorPort<T>::getConnectedProcessors() const {
    if (isOutport()) {
        tgtAssert(false, "Called CoProcessorPort::getConnectedProcessors on outport!");
        LERRORC("voreen.coprocessorport", "Called CoProcessorPort::getConnectedProcessors on outport!");
        return std::vector<T*>();
    }
    else {
        std::vector<T*> processors;
        for (size_t i = 0; i < connectedPorts_.size(); ++i) {
            T* p = static_cast<T*>(connectedPorts_[i]->getProcessor());

            processors.push_back(p);
        }
        return processors;
    }
}

template <typename T>
///Get the first connected Processor
T* GenericCoProcessorPort<T>::getConnectedProcessor() const {
    if (isOutport()) {
        tgtAssert(false, "Called CoProcessorPort::getConnectedProcessor on outport!");
        LERRORC("voreen.coprocessorport", "Called CoProcessorPort::getConnectedProcessor on outport!");
        return 0;
    }
    else {
        if (isConnected())
            return static_cast<T*>(connectedPorts_[0]->getProcessor());
        else
            return 0;
    }
}

} // namespace

#endif // VRN_GENERICCOPROCESSORPORT_H
