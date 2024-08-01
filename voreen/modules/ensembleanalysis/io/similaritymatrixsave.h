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

#ifndef VRN_SIMILARITYMATRIXSAVE_H
#define VRN_SIMILARITYMATRIXSAVE_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/boolproperty.h"

#include "../ports/similaritymatrixport.h"

namespace voreen {

class VRN_CORE_API SimilarityMatrixSave : public Processor {
public:
    SimilarityMatrixSave();

    virtual Processor* create() const         { return new SimilarityMatrixSave(); }
    virtual std::string getClassName() const  { return "SimilarityMatrixSave"; }
    virtual std::string getCategory() const   { return "Output"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_TESTING; }
    virtual bool isEndProcessor() const       { return true; }

protected:

    virtual void process();
    virtual void invalidate(int inv = 1);

private:

    void saveSimilarityMatrix();

    SimilarityMatrixPort inport_;

    FileDialogProperty filenameProp_;               ///< determines the name of the saved file
    ButtonProperty saveButton_;                     ///< triggers a save

    bool saveSimilarityMatrix_;          ///< used to determine, if process should save or not

    static const std::string loggerCat_;
};

}   //namespace

#endif
