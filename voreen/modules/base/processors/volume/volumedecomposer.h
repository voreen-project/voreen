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

#ifndef VRN_VOLUMEDECOMPOSER_H
#define VRN_VOLUMEDECOMPOSER_H

#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"

#include "voreen/core/ports/volumeport.h"

namespace voreen {

class ImageSequence;

/**
 * Decomposes a volume into slices that are put out as image sequence.
 */
class VRN_CORE_API VolumeDecomposer : public RenderProcessor {
public:
    VolumeDecomposer();
    ~VolumeDecomposer();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "VolumeDecomposer";  }
    virtual std::string getCategory() const   { return "Volume Processing"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_STABLE;   }

    virtual bool usesExpensiveComputation() const { return true; }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Decomposes a volume into z-aligned slices and puts them out as image sequence of float textures. "
                       "Multi-channel volumes are either copied channel-wise or converted to single-channel (grayscale), depending on the respective property setting. "
                       "The slice range to be extracted is specified by the \"Start Slice\" and \"End Slice\" properties. " );
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();
    virtual void adjustPropertiesToInput();

    void decomposeVolume();
    void clearSliceSequence();

    void startSliceChanged();
    void endSliceChanged();

    VolumePort inport_;
    ImageSequencePort outport_;

    IntProperty startSlice_;
    IntProperty endSlice_;
    BoolProperty convertMultiChannelToGrayscale_;

    ImageSequence* sliceSequence_;

    static const std::string loggerCat_; ///< category used in logging
};

}

#endif //VRN_VOLUMEDECOMPOSER_H
