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

#ifndef VRN_TEXTSAVE_H
#define VRN_TEXTSAVE_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/textport.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/filedialogproperty.h"

namespace voreen {

class VRN_CORE_API TextSave : public Processor {
public:
    TextSave();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "TextSave";        }
    virtual std::string getCategory() const  { return "Output";          }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE; }
    virtual bool isEndProcessor() const      { return true;              }

    virtual bool isReady() const { return true; }
    virtual void invalidate(int inv = 1);
protected:
    virtual void setDescriptions() {
        setDescription("Saves the input text to a text file.");
    }

    virtual void process();

    void saveFile();

    TextPort inport_;

    FileDialogProperty fileProp_;
    ButtonProperty saveButton_;
    BoolProperty continuousSave_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif
