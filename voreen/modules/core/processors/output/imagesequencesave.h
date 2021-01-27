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

#ifndef VRN_IMAGESEQUENCESAVE_H
#define VRN_IMAGESEQUENCESAVE_H

#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/ports/genericport.h"

#include <string>

namespace voreen {

class Volume;
class VolumeSerializer;
class VolumeSerializerPopulator;

class VRN_CORE_API ImageSequenceSave : public RenderProcessor {
public:
    ImageSequenceSave();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "ImageSequenceSave";  }
    virtual std::string getCategory() const   { return "Output";             }
    virtual CodeState getCodeState() const    { return CODE_STATE_STABLE;    }
    virtual bool isEndProcessor() const       { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Saves the input image sequence to a user-selected directory. The images will be saved as .png files with the specified basename. If no basename is given, the images' names are used.");
    }

    virtual void process();
    virtual void initialize();

    void saveSequence();
    void saveImage(const std::string& filename, tgt::Texture* image);

private:
    ImageSequencePort inport_;

    FileDialogProperty outputDirectory_;
    StringProperty baseName_;
    ButtonProperty saveButton_;
    BoolProperty continousSave_;

    static const std::string loggerCat_; ///< category used in logging
};

} //namespace

#endif
