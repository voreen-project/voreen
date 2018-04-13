/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_VOLUMESLICELOOPINITIATOR_H
#define VRN_VOLUMESLICELOOPINITIATOR_H

#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/loopport.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"

namespace voreen {

class VolumeSliceLoopInitiator : public RenderProcessor {

public:
    VolumeSliceLoopInitiator();
    virtual Processor* create() const;

    virtual std::string getCategory() const  { return "Utility"; }
    virtual std::string getClassName() const { return "VolumeSliceLoopInitiator"; }
    virtual CodeState getCodeState() const   { return Processor::CODE_STATE_EXPERIMENTAL; }
    virtual bool isUtility() const           { return true; }

    virtual bool isReady() const;

protected:
    /// Determines the current axis-alignment if the displayed slices.
    enum SliceAlignment {
        YZ_PLANE = 0,
        XZ_PLANE = 1,
        XY_PLANE = 2
    };
    friend class OptionProperty<SliceAlignment>;

    virtual void setDescriptions() {
        setDescription("");
    }

    virtual void process();

    virtual void invalidate(int inv = INVALID_RESULT);

    /**
     * Adapts the min/max ranges of the respective properties to the
     * dimensions of the currently connected volume.
     */
    void updateSliceProperties();

    /// Property containing the available alignments: xy (axial), xz (coronal), yz (sagittal)
    OptionProperty<SliceAlignment> sliceAlignment_;
    IntProperty minSliceIndex_;
    IntProperty maxSliceIndex_;

    VolumePort inport_;
    RenderPort outport_;
    LoopPort loopInport_;

    static const std::string loggerCat_; ///< category used in logging
};

}

#endif //VRN_IMAGESEQUENCELOOPINITIATOR_H
