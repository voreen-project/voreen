/***********************************************************************************
*                                                                                 *
* Voreen - The Volume Rendering Engine                                            *
*                                                                                 *
* Copyright (C) 2005-2017 University of Muenster, Germany.                        *
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

#include "similaritymatrixsave.h"

#include "voreen/core/voreenapplication.h"

#include "tgt/filesystem.h"

namespace voreen {

const std::string SimilarityMatrixSave::loggerCat_("voreen.ensembleanalysis.SimilarityMatrixSave");

SimilarityMatrixSave::SimilarityMatrixSave()
    : Processor()
    // ports
    , inport_(Port::INPORT, "inport", "Similarity Matrix Input", false)
    // properties
    , filenameProp_("filenameprop", "Save File as", "Select file...", VoreenApplication::app()->getUserDataPath(), "similarity matrix (*.sm)", FileDialogProperty::SAVE_FILE, Processor::INVALID_PATH)
    , saveButton_("saveButton", "Save", Processor::INVALID_PATH)
    // members
    , saveSimilarityMatrix_(true)
{
    addPort(inport_);

    addProperty(filenameProp_);
    addProperty(saveButton_);
}

void SimilarityMatrixSave::invalidate(int inv) {
    Processor::invalidate(inv);

    if (inv == Processor::INVALID_PATH && isInitialized()) {
        saveSimilarityMatrix_ = true;
    }
}

void SimilarityMatrixSave::process() {
    if (saveSimilarityMatrix_){
        saveSimilarityMatrix();
        saveSimilarityMatrix_ = false;
    }
}

void SimilarityMatrixSave::saveSimilarityMatrix() {
    if (!isInitialized())
        return;
    if (!inport_.hasData()) {
        LWARNING("no input field plot");
        return;
    }
    if (filenameProp_.get().empty()) {
        LWARNING("no filename specified");
        return;
    }

    try {
        std::fstream stream(filenameProp_.get(), std::ios::out);
        JsonSerializer json;
        Serializer s(json);
        inport_.getData()->serialize(s);
        json.write(stream, false, true);
        LINFO(filenameProp_.get() << " saved sucessfully!");
    } catch(tgt::FileException& e) {
        LERROR(e.what());
        filenameProp_.set("");
    }
}

}   // namespace
