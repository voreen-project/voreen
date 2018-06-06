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

#ifndef VRN_SLICEGENERATOR_H
#define VRN_SLICEGENERATOR_H

#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"
#include "voreen/core/interaction/mwheelnumpropinteractionhandler.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/ports/volumeport.h"

namespace voreen {

class SliceGenerator : public RenderProcessor {
    /// Determines the current axis-alignment if the displayed slices.
    enum SliceAlignment {
        YZ_PLANE = 0,
        XZ_PLANE = 1,
        XY_PLANE = 2
    };
    friend class OptionProperty<SliceAlignment>;

public:
    SliceGenerator();
    virtual ~SliceGenerator();

    virtual std::string getCategory() const { return "Volume"; }
    virtual std::string getClassName() const { return "SliceGenerator"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }
    virtual std::string getProcessorInfo() const;

    virtual Processor* create() const { return new SliceGenerator(); }

    virtual void process();
    virtual void portResized(RenderPort* /*p*/, tgt::ivec2 /*newsize*/) {}

    bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("");
    }

    /**
     * Adapts the min/max ranges of the respective properties to the
     * dimensions of the currently connected volume.
     */
    void updateSliceProperties();

    /// Property containing the available alignments: xy (axial), xz (coronal), yz (sagittal)
    OptionProperty<SliceAlignment> sliceAlignment_;
    IntProperty sliceIndex_;                ///< Property containing the currently selected slice
    BoolProperty applyDatasetTransformationMatrix_;

    MWheelNumPropInteractionHandler<int> mwheelCycleHandler_;

    RenderPort outport_;
    VolumePort inport_;
    GeometryPort geomPort_;

    MeshListGeometry geometry_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif //VRN_SLICEGENERATOR_H
