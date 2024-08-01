/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_VORTEXPORT_H
#define VRN_VORTEXPORT_H

#include "voreen/core/ports/genericport.h"
#include "../datastructures/vortex.h"

namespace voreen {
#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API GenericPort<Vortex>;
#endif

class VRN_CORE_API VortexPort : public GenericPort<Vortex> {
public:
    VortexPort( PortDirection direction, const std::string& id, const std::string& guiName = {}, bool allowMultipleConnections = false, Processor::InvalidationLevel invalidationLevel = Processor::INVALID_RESULT );

    virtual Port* create( PortDirection direction, const std::string& id, const std::string& guiName = {} ) const;
    virtual std::string getClassName() const;
    virtual std::string getContentDescriptionHTML() const;
};

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API GenericPort<VortexCollection>;
#endif

class VRN_CORE_API VortexCollectionPort : public GenericPort<VortexCollection> {
public:
    VortexCollectionPort( PortDirection direction, const std::string& id, const std::string& guiName = {}, bool allowMultipleConnections = false, Processor::InvalidationLevel invalidationLevel = Processor::INVALID_RESULT );

    virtual Port* create( PortDirection direction, const std::string& id, const std::string& guiName = {} ) const;
    virtual std::string getClassName() const;
    virtual std::string getContentDescriptionHTML() const;
};

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API GenericPort<std::vector<Vortex>>;
#endif

class VRN_CORE_API VortexListPort : public GenericPort<std::vector<Vortex>> {
public:
    VortexListPort( PortDirection direction, const std::string& id, const std::string& guiName = {}, bool allowMultipleConnections = false, Processor::InvalidationLevel invalidationLevel = Processor::INVALID_RESULT );

    virtual Port* create( PortDirection direction, const std::string& id, const std::string& guiName = {} ) const;
    virtual std::string getClassName() const;
    virtual std::string getContentDescriptionHTML() const;
};

}

#endif // VRN_VORTEXPORT_H
