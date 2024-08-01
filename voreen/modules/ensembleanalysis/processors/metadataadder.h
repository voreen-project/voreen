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

#ifndef METADATAADDER_H
#define METADATAADDER_H

//base class headers
#include "voreen/core/processors/processor.h"

//Port headers
#include "voreen/core/ports/genericport.h"

//Property headers
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/stringproperty.h"


namespace voreen {

class VRN_CORE_API MetaDataAdder : public Processor {
public:
    MetaDataAdder();
    virtual Processor* create() const;
    virtual std::string getClassName() const         { return "MetaDataAdder"; }
    virtual std::string getCategory() const          { return "Processing"; }
    virtual CodeState getCodeState() const           { return CODE_STATE_EXPERIMENTAL; }

protected:
    virtual void process();
    virtual void setDescriptions() {
        setDescription("This processor can add explicit time steps to the data as well as "
                   "a name. The time information should be given in a file where each line "
                   "contains one number which is used as a timestep.");
    }

private:
    VolumeListPort inport_;
    VolumeListPort outport_;

    BoolProperty addTime_;
    FileDialogProperty timeInformationFile_;
    StringProperty modalityString_;

    std::vector<std::unique_ptr<VolumeBase>> decorators_;
    void clearOutput();
};

} //namespace

#endif // METADATAADDER_H
