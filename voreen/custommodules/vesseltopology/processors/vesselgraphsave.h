/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2016 University of Muenster, Germany.                        *
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

#ifndef VRN_VESSELGRAPHSAVE_H
#define VRN_VESSELGRAPHSAVE_H

#include "voreen/core/processors/processor.h"

#include "../ports/vesselgraphport.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/boolproperty.h"

namespace voreen {

class VesselGraphSave : public Processor {
public:
    VesselGraphSave();
    virtual ~VesselGraphSave();
    virtual std::string getCategory() const { return "Geometry"; }
    virtual std::string getClassName() const { return "VesselGraphSave"; }
    virtual CodeState getCodeState() const { return Processor::CODE_STATE_EXPERIMENTAL; }
    virtual Processor* create() const { return new VesselGraphSave(); }
    virtual bool isEndProcessor() const       { return true;              }

protected:
    virtual void setDescriptions() {
        setDescription("This processor can be used to save vessel graphs as files that can be later reloaded using VesselGraphSource.");
    }
    virtual void saveCurrentGraph();

    virtual void process();

    VesselGraphPort inport_;

    // properties
    FileDialogProperty graphFilePath_;
    ButtonProperty saveButton_;
    BoolProperty continousSave_;
    BoolProperty prettyJson_;

    static const std::string loggerCat_;
};

} // namespace voreen
#endif // VRN_VESSELGRAPHSAVE_H
