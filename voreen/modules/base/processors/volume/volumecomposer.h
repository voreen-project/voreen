/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_IMAGESTACKER_H
#define VRN_IMAGESTACKER_H

#include "voreen/core/processors/volumeprocessor.h"

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/boolproperty.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/datastructures/volume/volumelist.h"

namespace voreen {

class ImageSequence;

/**
 * Constructs a volume from a sequence of 2D images.
 */
class VRN_CORE_API VolumeComposer : public VolumeProcessor {
public:
    VolumeComposer();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "VolumeComposer";    }
    virtual std::string getCategory() const   { return "Volume Processing"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_STABLE;   }

    virtual bool usesExpensiveComputation() const { return true; }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Constructs a volume by stacking either the input image sequence or the input volume collection in z-direction. "
                       "Multi-channel input images are either copied channel-wise or converted to grayscale, depending on the respective property setting. "
                       "If the volume is constructed from an image stack, the voxel spacing of the resulting volume has to be specified via property. ");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    void stackImages();
    void stackVolumes();

    ImageSequencePort inportImages_;
    VolumeListPort inportVolumes_;
    VolumePort outport_;

    BoolProperty convertMultiChannelTexToGrayscale_;
    FloatVec3Property voxelSpacing_;

    static const std::string loggerCat_; ///< category used in logging
};

}

#endif //VRN_IMAGESTACKER_H
