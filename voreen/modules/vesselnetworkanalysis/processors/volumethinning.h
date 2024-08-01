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

#ifndef VRN_VOLUMETHINNING_H
#define VRN_VOLUMETHINNING_H

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"

#include "../algorithm/volumemask.h"

namespace voreen {

/**
 * A processor that iteratively thins a (binary) input volume
 * until a voxel skeleton remains.
 */
class VolumeThinning : public VolumeProcessor {
public:
    VolumeThinning();
    virtual ~VolumeThinning();

    virtual std::string getClassName() const         { return "VolumeThinning";      }
    virtual std::string getCategory() const       { return "Volume Processing"; }
    virtual std::string setDescriptions() const       { return "Volume Processing"; }
    virtual VoreenSerializableObject* create() const;
    virtual void setDescriptions() {
        setDescription(
                "Volume processor that reduces a binary volume to a one voxel wide skeleton. "
                "The binary volume input is required and is skeletonized. "
                "Additionally, another binary volume (sample mask) can optionally be used to restrict the skeletonization process to a certain arbitrary subvolume. "
                "Voxels outside the sample mask are considered to be foreground for the skeletonization, but will be written as background to the output volume. "
                );
        binarizationThreshold_.setDescription("Values above this threshold will be considered foreground, others background. If the input volume is not binary already, this property can therefore be used for thresholding.");
    }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;   }
    //virtual bool usesExpensiveComputation() const { return true; }
    virtual bool isReady() const;
    virtual void process();

private:
    VolumePort inport_;
    VolumePort sampleMaskInport_;
    VolumePort outport_;

    BoolProperty enabledProp_;
    OptionProperty<VolumeMask::ThinningAlgorithm> thinningAlgorithm_;
    FloatProperty binarizationThreshold_;
    IntProperty maxSteps_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_VOLUMETHINNING_H
