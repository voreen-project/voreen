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

#ifndef VRN_PATCHLISTREADER_H
#define VRN_PATCHLISTREADER_H

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/volumeinfoproperty.h"
#include "voreen/core/properties/boolproperty.h"

#include "voreen/core/datastructures/volume/volumelist.h"

namespace voreen {

class VRN_CORE_API PatchListReader : public VolumeProcessor {
public:
    PatchListReader();
    //virtual ~VolumeSave();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "PatchListReader";      }
    virtual std::string getCategory() const   { return "Input";          }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }

protected:
    virtual void setDescriptions() {
        setDescription("Reads 7x7x7 (normalized) float patches from an input binary raw file and outputs them as a VolumeList.");
    }

    virtual void deinitialize();

    virtual void process();

    virtual void loadPatches();

    virtual void clearOutputList();


private:

    VolumeListPort outport_;

    VolumeList* outputList_;    ///< list containing all output volumes (memory management is handled by processor)

    FileDialogProperty filename_;
    BoolProperty normalizePatches_; // normalizes the patches with mean 0 and standard deviation 1
    BoolProperty allowAdditionalFeaturesAtEnd_;
    ButtonProperty loadButton_;
    
    BoolProperty autoCompute_;

    static const std::string loggerCat_;
};

} // namespace
#endif 
