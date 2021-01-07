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

#include "vesselgraphsave.h"

#include "voreen/core/io/serialization/jsonserializer.h"

namespace voreen {

const std::string VesselGraphSave::loggerCat_("voreen.vesseltoplogy.vesselgraphsave");


VesselGraphSave::VesselGraphSave()
    : Processor()
    , inport_(Port::INPORT, "graph.input", "Graph Input", false, Processor::INVALID_RESULT)
#ifdef WIN32
    , graphFilePath_("graphFilePath", "Voreen Vessel Graph File", "Voreen Vessel Graph File", "", "uncompressed (*.vvg)", FileDialogProperty::SAVE_FILE)
#else
    , graphFilePath_("graphFilePath", "Voreen Vessel Graph File", "Voreen Vessel Graph File", "", "compressed (*.vvg.gz);;uncompressed (*.vvg)", FileDialogProperty::SAVE_FILE)
#endif
    , saveButton_("save", "Save")
    , continousSave_("continousSave", "Save on inport change", false)
    , prettyJson_("prettyJson", "Prettify Json", false)
{
    addPort(inport_);

    addProperty(graphFilePath_);
        ON_CHANGE(graphFilePath_, VesselGraphSave, saveCurrentGraph);
    addProperty(saveButton_);
        ON_CHANGE(saveButton_, VesselGraphSave, saveCurrentGraph);
    addProperty(continousSave_);
    addProperty(prettyJson_);
}

VesselGraphSave::~VesselGraphSave() {
}

void VesselGraphSave::saveCurrentGraph() {
    const std::string path = graphFilePath_.get();
    if(path.empty()) {
        return;
    }
    const VesselGraph* input = inport_.getData();
    if(!input) {
        return;
    }

    JsonSerializer serializer;
    try {
        serializer.serialize("graph", *input);
    } catch(SerializationException& s) {
        LERROR("Could not serialize graph: " << s.what());
        return;
    }
    try {
        std::fstream f(path, std::ios::out);
        bool compressed = tgt::FileSystem::fileExtension(path) == "gz";
        serializer.write(f, prettyJson_.get(), compressed);
    } catch(...) {
        LERROR("Could not save graph " << path);
        return;
    }
    LINFO("Saved graph to " << path << ".");
}

void VesselGraphSave::process() {
    if(inport_.hasChanged() && continousSave_.get()) {
        saveCurrentGraph();
    }
}
} // namespace voreen
