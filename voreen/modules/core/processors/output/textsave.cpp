/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
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

#include "textsave.h"

#include <fstream>

namespace voreen {

const std::string TextSave::loggerCat_("voreen.core.TextSave");

TextSave::TextSave()
    : Processor()
    ,  inport_(Port::INPORT, "inport", "Text Input")
    ,  fileProp_("file", "Text File", "Select text file...", "./",
            "Text files (*.txt)", FileDialogProperty::SAVE_FILE, Processor::INVALID_PATH)
    , saveButton_("save", "Save")
    , continuousSave_("continousSave", "Save continuously", false)
{
    addPort(inport_);

    saveButton_.onChange(MemberFunctionCallback<TextSave>(this, &TextSave::saveFile));
    addProperty(fileProp_);
    addProperty(saveButton_);
    addProperty(continuousSave_);
}

Processor* TextSave::create() const {
    return new TextSave();
}

void TextSave::invalidate(int inv) {
    Processor::invalidate(inv);
    //auto save on path change
    if(inv == Processor::INVALID_PATH && isInitialized())
        saveFile();
}

void TextSave::process() {
    if (inport_.hasChanged() && continuousSave_.get())
        saveFile();
}

void TextSave::saveFile() {
    if (!inport_.hasData())
        return;
    const std::string text = inport_.getData();

    std::string filename = fileProp_.get();
    if (filename == "") {
        LWARNING("Could not save text file: filename is empty");
        return;
    }

    LINFO("Saving text to file: " << filename);

    std::fstream stream(fileProp_.get().c_str(), std::ios::out);
    if (stream.fail()) {
        LERROR("Failed to open file for writing: " << filename);
    }
    else {
        stream << text;
        stream.close();
    }
}

} // namespace voreen
