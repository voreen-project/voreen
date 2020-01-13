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

#include "geometrysave.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include <fstream>

namespace voreen {

const std::string GeometrySave::loggerCat_("voreen.core.GeometrySave");

GeometrySave::GeometrySave()
    : Processor()
    ,  inport_(Port::INPORT, "inport", "Geometry Input")
    ,  fileProp_("file", "Geometry File", "Select Voreen geometry file...", "./",
            "Voreen Geometry files (*.vge)", FileDialogProperty::SAVE_FILE, Processor::INVALID_PATH)
    ,  saveButton_("save", "Save")
    ,  continousSave_("continousSave", "Save continuously", false)
{
    addPort(inport_);

    saveButton_.onChange(MemberFunctionCallback<GeometrySave>(this, &GeometrySave::saveFile));
    addProperty(fileProp_);
    addProperty(saveButton_);
    addProperty(continousSave_);
}

Processor* GeometrySave::create() const {
    return new GeometrySave();
}

void GeometrySave::invalidate(int inv) {
    Processor::invalidate(inv);
    //auto save on path change
    if(inv == Processor::INVALID_PATH && isInitialized())
        saveFile();
}

void GeometrySave::process() {
    if (inport_.hasChanged() && continousSave_.get())
        saveFile();
}

namespace {
    void saveAsVGE(const Geometry* geometry, std::fstream& f) {
        XmlSerializer s;
        try {
            s.serialize("Geometry", geometry);
        }
        catch (VoreenException& e) {
            LERRORC("voreen.core.GeometrySave", "Failed to serialize geometry: " + std::string(e.what()));
            return;
        }
        s.write(f);
        f.close();
    }
    void saveAsObj(const GlMeshGeometryBase* geometry, std::fstream& f) {
        geometry->exportAsObj(f);
    }
    void saveAsStl(const GlMeshGeometryBase* geometry, std::fstream& f) {
        geometry->exportAsStl(f);
    }
}

void GeometrySave::saveFile() {
    const Geometry* geometry = inport_.getData();
    if (!geometry)
        return;

    std::string filename = fileProp_.get();
    if (filename.empty()) {
        LWARNING("Could not save geometry: filename is empty");
        return;
    }

    LINFO("Saving Voreen Geometry to file: " << filename);

    std::fstream stream(filename.c_str(), std::ios::out);
    if (stream.fail()) {
        LERROR("Failed to open file for writing: " << filename);
    }
    else {
        auto extension = tgt::FileSystem::fileExtension(filename, true);
        if(extension == "vge") {
            saveAsVGE(geometry, stream);
        }
        else if(extension == "obj") {
            if(const GlMeshGeometryBase* geom = dynamic_cast<const GlMeshGeometryBase*>(geometry)) {
                saveAsObj(geom, stream);
            } else {
                LERROR("Unsupported output file type '" << extension << "' for input geometry!");
            }
        } else if(extension == "stl") {
            if(const GlMeshGeometryBase* geom = dynamic_cast<const GlMeshGeometryBase*>(geometry)) {
                saveAsStl(geom, stream);
            } else {
                LERROR("Unsupported output file type '" << extension << "' for input geometry!");
            }
        }
        else {
            LERROR("Unsupported output file type '" << extension << "'!");
        }
    }
}

void GeometrySave::adjustPropertiesToInput() {
    const Geometry* geometry = inport_.getData();
    if(const GlMeshGeometryBase* geom = dynamic_cast<const GlMeshGeometryBase*>(geometry)) {
        fileProp_.setFileFilter("Voreen Geometry files (*.vge);; Wavefront Object Files (*.obj);; Stereolithography Files (*.stl)");
    } else {
        fileProp_.setFileFilter("Voreen Geometry files (*.vge)");
    }
}

} // namespace voreen
