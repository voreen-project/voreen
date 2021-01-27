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
    virtual std::string getCategory() const { return "Output"; }
    virtual std::string getClassName() const { return "VesselGraphSave"; }
    virtual CodeState getCodeState() const { return Processor::CODE_STATE_TESTING; }
    virtual Processor* create() const { return new VesselGraphSave(); }
    virtual bool isEndProcessor() const       { return true;              }

protected:
    virtual void setDescriptions() {
        setDescription("This processor can be used to save vessel graphs as files that can be later reloaded using VesselGraphSource. "
                "Vesselgraphs are serialized in a custom (but simple) json format that is gzip-compressed before writing it to disk.");
        graphFilePath_.setDescription("Location on disk where the serialized version of the graph will be written to.");
        continousSave_.setDescription("Automatically overwrite specified file when the input VesselGraph changes. If not enabled, the export can be triggered manually using the 'Save'-Button.");
        prettyJson_.setDescription("Generate json with human-friendly newlines and identation. If not specified, no superfluous space characters will be written.");
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
