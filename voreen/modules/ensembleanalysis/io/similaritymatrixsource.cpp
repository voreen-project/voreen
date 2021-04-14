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

#include "similaritymatrixsource.h"

#include "voreen/core/voreenapplication.h"

#include "tgt/filesystem.h"

namespace voreen {

const std::string SimilarityMatrixSource::loggerCat_("voreen.ensembleanalysis.SimilarityMatrixSource");

SimilarityMatrixSource::SimilarityMatrixSource()
    : Processor()
    // ports
    , outport_(Port::OUTPORT, "outport", "Similarity Matrix Output", false)
    // properties
    , filenameProp_("filenameprop", "Load Similarity Matrix File from", "Select file...", VoreenApplication::app()->getUserDataPath(), "Voreen Similarity Matrix (*.vsm)", FileDialogProperty::OPEN_FILE, Processor::INVALID_PATH)
    , loadButton_("loadButton", "Load", INVALID_PATH)
    // members
    , loadSimilarityMatrix_(true)
{
    addPort(outport_);

    addProperty(filenameProp_);
    addProperty(loadButton_);
}

void SimilarityMatrixSource::invalidate(int inv) {
    Processor::invalidate(inv);

    if (inv == Processor::INVALID_PATH && isInitialized()) {
        loadSimilarityMatrix_ = true;
    }
}

void SimilarityMatrixSource::process() {
    if (loadSimilarityMatrix_){
        loadSimilarityMatrix();
        loadSimilarityMatrix_ = false;
    }
}

void SimilarityMatrixSource::loadSimilarityMatrix() {
    if (!isInitialized())
        return;

    outport_.setData(nullptr);

    if (filenameProp_.get().empty()) {
        LWARNING("no filename specified");
        return;
    }

    try {
        std::unique_ptr<SimilarityMatrixList> similarityMatrices(new SimilarityMatrixList());

        std::ifstream stream(filenameProp_.get());
        if(!stream)
            throw tgt::CorruptedFileException("Could not read file: " + filenameProp_.get());
        JsonDeserializer json;
        json.read(stream, false);
        Deserializer s(json);
        s.deserialize("similarity", *similarityMatrices);
        outport_.setData(similarityMatrices.release(), true);
        LINFO(filenameProp_.get() << " loaded sucessfully!");
    } catch(tgt::Exception& e) {
        LERROR(e.what());
        filenameProp_.set("");
    }
}

}   // namespace
