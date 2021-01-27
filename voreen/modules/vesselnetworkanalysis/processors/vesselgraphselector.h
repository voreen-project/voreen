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

#ifndef VRN_VESSELGRAPHSELECTOR_H
#define VRN_VESSELGRAPHSELECTOR_H

#include "voreen/core/processors/processor.h"
#include "../ports/vesselgraphport.h"
#include "../ports/vesselgraphlistport.h"
#include "voreen/core/properties/intproperty.h"
#include "../datastructures/vesselgraph.h"

namespace voreen {

class Volume;
class ProcessorWidgetFactory;

/**
 * Selects a single volume out of a input list.
 */
class VRN_CORE_API VesselGraphSelector : public Processor {

public:
    VesselGraphSelector();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "VesselGraphSelector";  }
    virtual std::string getCategory() const  { return "Input";           }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL; }
    virtual bool isUtility() const           { return true; }

    virtual void adjustPropertiesToInput();

protected:
    virtual void setDescriptions() {
        setDescription("Selects a single VesselGraph from the input list. "
                "So far, this processor is mostly intended to be used in conjuction with the debug output of a <b>VesselGraphCreator</b>.");
        graphID_.setDescription("Select the graph with the specified number from the list.");
        debugVolumeID_.setDescription("Select a type of debug volume. This is not actually used in this processor, but used to calculate the resulting debug volume selector id.");
        resultingDebugVesselGraphSelectorID_.setDescription("Link with a VolumeListSelector which is fed by the debug volumes of a <b>VesselGraphCreator</b>.");
    }

    virtual void process();

    void syncResultingDebugVesselGraphSelectorID();

    VesselGraphListPort inport_;
    VesselGraphPort outport_;

    IntProperty graphID_;
    IntProperty debugVolumeID_;
    IntProperty resultingDebugVesselGraphSelectorID_; // Link with VesselGraphSelector

    static const std::string loggerCat_;
};

} // namespace

#endif
