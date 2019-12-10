/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_ALIGNEDSLICEPROXYGEOMETRY_H
#define VRN_ALIGNEDSLICEPROXYGEOMETRY_H

#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/datastructures/volume/slice/slicehelper.h"
#include "voreen/core/interaction/mwheelnumpropinteractionhandler.h"

#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/buttonproperty.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/textport.h"

namespace voreen {

class AlignedSliceProxyGeometry : public Processor {
    /// Determines the current axis-alignment if the displayed slices.
public:
    AlignedSliceProxyGeometry();
    virtual ~AlignedSliceProxyGeometry();

    virtual std::string getCategory() const { return "Volume"; }
    virtual std::string getClassName() const { return "AlignedSliceProxyGeometry"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_TESTING; }
    virtual std::string getProcessorInfo() const;

    virtual Processor* create() const { return new AlignedSliceProxyGeometry(); }

    virtual void process();
    virtual bool isReady() const;

    virtual void adjustPropertiesToInput();

protected:
    virtual void setDescriptions() {
        setDescription("Generates a proxy-geometry representing slices in a volume. Use in combination with MultiSliceRenderer.");
    }

    /// Adapts the min/max ranges of the respective properties to the dimensions of the currently connected volume.
    void alignCamera();

    void indexChanged();
    void floatIndexChanged();
    void update();

    /// Property containing the available alignments: xy (axial), xz (coronal), yz (sagittal)
    OptionProperty<SliceAlignment> sliceAlignment_;
    /// Property containing the currently selected slice
    IntProperty sliceIndex_;
    FloatProperty floatSliceIndex_;
    BoolProperty restrictToMainVolume_;
    CameraProperty camera_;
    ButtonProperty alignCameraButton_;
    FloatVec3Property plane_;
    FloatProperty planeDist_;

    MWheelNumPropInteractionHandler<int> mwheelCycleHandler_;

    VolumePort inport_;
    VolumePort secondaryVolumePort_;
    GeometryPort geomPort_;
    TextPort textPort_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif
