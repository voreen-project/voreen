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

#ifndef VRN_VESSELGRAPHLISTPORT_H
#define VRN_VESSELGRAPHLISTPORT_H

#include "voreen/core/ports/genericport.h"
#include "../datastructures/vesselgraph.h"
#include <vector>

namespace voreen {

// A generic port for vectors of VesselGraphs. Nothing special to see here.
class VesselGraphListPort : public GenericPort<std::vector<VesselGraph>> {
public:
    VesselGraphListPort(PortDirection direction, const std::string& name, const std::string& guiName = "",
                    bool allowMultipleConnections = false,
                    Processor::InvalidationLevel invalidationLevel = Processor::INVALID_PARAMETERS);
    ~VesselGraphListPort();

    virtual Port* create(PortDirection direction, const std::string& id, const std::string& guiName = "") const;

    virtual std::string getClassName() const {return "VesselGraphListPort";}
    virtual tgt::col3 getColorHint() const;

    virtual std::string getContentDescription() const;

    virtual std::string getContentDescriptionHTML() const;

protected:
    std::string portData_;

    static const std::string loggerCat_; ///< category used in logging
};

} // namespace voreen

#endif // VRN_VESSELGRAPHLISTPORT_H
