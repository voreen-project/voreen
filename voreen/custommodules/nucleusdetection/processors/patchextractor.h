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

#ifndef VRN_PATCHEXTRACTOR_H
#define VRN_PATCHEXTRACTOR_H

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/volumeinfoproperty.h"

namespace voreen {

class VRN_CORE_API PatchExtractor : public VolumeProcessor {
public:
    PatchExtractor();
    //virtual ~VolumeSave();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "PatchExtractor";      }
    virtual std::string getCategory() const   { return "Output";          }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }
    virtual bool isEndProcessor() const       { return true;              }

    virtual bool usesExpensiveComputation() const { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Extracts a 7x7x7 patch around each voxel and stores the (normalized float) values in a binary file.");
    }

    virtual void process();

    virtual void computePatches();

private:

    VolumeListPort inport_;

    FileDialogProperty filename_;

    BoolProperty normalizePatches_;

    ButtonProperty computeButton_;

    ProgressProperty progressProperty_;

    bool startComputation_;

    static const std::string loggerCat_;
};

} // namespace
#endif 
