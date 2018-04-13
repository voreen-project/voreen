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

#ifndef VRN_VOLUMESUBSET_H
#define VRN_VOLUMESUBSET_H

#include "voreen/core/processors/volumeprocessor.h"

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"
#include "voreen/core/ports/geometryport.h"

namespace voreen {

/**
 * Crops the input volume by axis-aligned clipping planes.
 */
class VRN_CORE_API VolumeCrop : public CachingVolumeProcessor {
public:
    VolumeCrop();
    virtual Processor* create() const;

    virtual std::string getClassName() const      { return "VolumeCrop";        }
    virtual std::string getCategory() const       { return "Volume Processing"; }
    virtual CodeState getCodeState() const        { return CODE_STATE_STABLE;   }
    virtual bool usesExpensiveComputation() const { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Crops the input volume by axis-aligned clipping planes. Both the min and max value of each dimension are included in the cropped volume. Optionally, a bounding box (in world coordinates) can be supplied via a Geometry port that will then be used to automatically set dimensions.");
    }

    virtual void process();
    virtual void adjustPropertiesToInput();
    virtual bool isReady() const;

private:
    /// Crops the input volume according to the property values
    /// and writes the result to the outport.
    void crop();
    void inportChanged();
    void updateInfoPropertys();


    VolumePort inport_;
    GeometryPort boundingBoxPort_;
    VolumePort outport_;

    IntBoundingBoxProperty clipRegion_;
    BoolProperty continuousCropping_;   ///< Crop on each change of the input volume.
    ButtonProperty button_;             ///< Perform cropping.

    /// Read-only property displaying the dimensions of the cropped volume.
    IntVec3Property croppedDimensions_;
    /// Read-only property displaying the data size of the cropped volume in MB.
    IntProperty croppedSize_;
    // Read-only property displaying if the clipregion matches the croped region
    BoolProperty clipRegionMatchesCroppedRegion_;

    bool isCropped_;
    tgt::IntBounds croppedRegion_;

    static const std::string loggerCat_; ///< category used in logging
};

}   //namespace

#endif
