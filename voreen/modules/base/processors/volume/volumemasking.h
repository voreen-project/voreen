/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_VOLUMEMASKING_H
#define VRN_VOLUMEMASKING_H

#include <string>
#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/optionproperty.h"

namespace voreen {

class Volume;
/**
 * Masks the incoming volume based on a masking function.
 */
class VRN_CORE_API VolumeMasking : public CachingVolumeProcessor {
public:
    VolumeMasking();
    virtual ~VolumeMasking();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "VolumeMasking";     }
    virtual std::string getCategory() const   { return "Volume Processing"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_STABLE;  }

protected:
    virtual void setDescriptions() {
        setDescription("Masks the incoming volume based on the masking function, such that voxels are discarded, if the masking function assigns zero opacity to their intensities.");
    }

    virtual void process();
    virtual void beforeProcess();

private:
    void forceUpdate();
    void maskVolume();

    VolumePort inport_;
    VolumePort outport_;

    BoolProperty enableProcessingProp_;
    TransFunc1DKeysProperty maskFunction_;
    StringOptionProperty passedVoxelAction_;
    BoolProperty continousUpdate_;
    ButtonProperty updateButton_;

    bool forceUpdate_;
};

}   //namespace

#endif
