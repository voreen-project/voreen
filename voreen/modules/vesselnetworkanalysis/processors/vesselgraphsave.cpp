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

#include "vesselgraphsave.h"

#include "voreen/core/io/serialization/jsonserializer.h"

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#if WIN32 // HACK: Needed due to bug #13164 (described boost bugtracker)
#define BOOST_IOSTREAMS_DYN_LINK
#endif
#include <boost/iostreams/filter/gzip.hpp>

namespace voreen {

const std::string VesselGraphSave::loggerCat_("voreen.vesseltoplogy.vesselgraphsave");


VesselGraphSave::VesselGraphSave()
    : Processor()
    , inport_(Port::INPORT, "graph.input", "Graph Input", false, Processor::INVALID_RESULT)
    , graphFilePath_("graphFilePath", "Voreen Vessel Graph File", "Voreen Vessel Graph File", "",
            "Voreen Vessel Graph File, uncompressed json (*.vvg);;"
#ifndef WIN32 // Compression does not work on windows, see below
            "Voreen Vessel Graph File, compressed json (*.vvg.gz);;"
#endif
            "Wavefront OBJ, centerlines only (*.obj)",
            FileDialogProperty::SAVE_FILE)
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

static size_t writeVertex(std::ostream& s, size_t& counter, const tgt::vec3& v) {
    s << "v " << v.x << " " << v.y << " " << v.z << "\n";
    size_t id = counter;
    ++counter;
    return id;
}

static void saveAsWavefrontObj(std::ostream& s, const VesselGraph& g) {
    size_t vertexCounter = 1;
    for(auto& node : g.getNodes()) {
        s << "#Node " << node.getID().raw() << "\n";
        size_t id = writeVertex(s, vertexCounter, node.pos_);
        s << "p " << id << "\n";
    }

    for(auto& edge : g.getEdges()) {
        s << "#Edge " << edge.getID().raw() << "\n";
        size_t centerlineBegin = vertexCounter;
        for(auto& voxel : edge.getVoxels()) {
            writeVertex(s, vertexCounter, voxel.pos_);
        }
        size_t centerlineEnd = vertexCounter;

        // We use the fact that nodes above are iterated over in id order.
        // However, VesselGraphNode ids start at 0 while OBJ vertex ids start
        // at 1.
        s << "l " << edge.getNodeID1().raw() + 1;
        for(int i = centerlineBegin; i < centerlineEnd; ++i) {
            s << " " << i;
        }
        s << " " << edge.getNodeID2().raw() + 1 << "\n";
    }
}
static void saveAsJson(std::ostream& f, const VesselGraph& g, bool pretty, bool compressed) {
    if(pretty) {
        JsonSerializer serializer;
        serializer.serialize("graph", g);
        serializer.write(f, pretty, compressed);
    } else {
        // Faster and less memory intensive than rapidjson (i.e.,
        // JsonSerializer), but does not support prettifying the json

        using namespace boost::iostreams;
        if(compressed) {
#ifdef WIN32
            // TODO: Using compression produces linker errors in boost.
            // This does, however, not happen in json(de)serializer.cpp.
            // Instead, there the compressed stream is corrupted.
            tgtAssert(false, "Compression currently not supported on windows");
#else
            filtering_ostream compressingStream;
            compressingStream.push(gzip_compressor(
                        gzip_params(
                            gzip::best_compression
                            )
                        ));
            compressingStream.push(f);
            g.serializeToJson(compressingStream);
#endif
        } else {
            g.serializeToJson(f);
        }
    }
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
    try {
        std::fstream f(path, std::ios::out);
        if(endsWith(path, ".vvg.gz")) {
            saveAsJson(f, *input, prettyJson_.get(), true);
        } else if(endsWith(path, ".vvg")) {
            saveAsJson(f, *input, prettyJson_.get(), false);
        } else if(endsWith(path, ".obj")) {
            saveAsWavefrontObj(f, *input);
        }
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
