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

#include "binarygeometrysave.h"
#include "binarygeometry.h"

#include <fstream>

namespace voreen {

const std::string BinaryGeometrySave::loggerCat_("voreen.sciviscontest2022.BinaryGeometrySave");

BinaryGeometrySave::BinaryGeometrySave()
    : Processor()
    ,  inport_(Port::INPORT, "inport", "Geometry Input")
    ,  fileProp_("file", "Geometry File", "Select Voreen geometry file...", "./",
            "Voreen Binary Geometry files (*.vbge)", FileDialogProperty::SAVE_FILE, Processor::INVALID_PATH)
    ,  saveButton_("save", "Save")
    ,  continousSave_("continousSave", "Save continuously", false)
{
	Processor::addPort(inport_);

    saveButton_.onChange(MemberFunctionCallback<BinaryGeometrySave>(this, &BinaryGeometrySave::saveFile));
	Processor::addProperty(fileProp_);
	Processor::addProperty(saveButton_);
	Processor::addProperty(continousSave_);
}

Processor* BinaryGeometrySave::create() const {
    return new BinaryGeometrySave();
}

void BinaryGeometrySave::invalidate(int inv) {
    Processor::invalidate(inv);
    //auto save on path change
    if(inv == Processor::INVALID_PATH && isInitialized())
        saveFile();
}

void BinaryGeometrySave::process() {
    if (inport_.hasChanged() && continousSave_.get())
        saveFile();
}

void BinaryGeometrySave::saveFile() {
    const Geometry* geometry = inport_.getData();
    if (!geometry)
        return;

    std::string filename = fileProp_.get();
    if (filename.empty()) {
        LWARNING("Could not save geometry: filename is empty");
        return;
    }

    LINFO("Saving Voreen Geometry to file: " << filename);

	auto extension = tgt::FileSystem::fileExtension(filename, true);
	if(extension == "vbge") {
		BinaryFileParser::write_file(geometry, filename.c_str());
	}
	else {
		LERROR("Unsupported output file type '" << extension << "'!");
	}
}

} // namespace voreen
