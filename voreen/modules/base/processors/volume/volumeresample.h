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

#ifndef VRN_VOLUMERESAMPLE_H
#define VRN_VOLUMERESAMPLE_H

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/optionproperty.h"

namespace voreen {

class Volume;

/**
 * Resizes the input volume to the specified dimensions
 * by using a selectable filtering mode.
 */
class VRN_CORE_API VolumeResample : public CachingVolumeProcessor {
public:
    VolumeResample();
    virtual Processor* create() const;

    virtual std::string getClassName() const      { return "VolumeResample";    }
    virtual std::string getCategory() const       { return "Volume Processing"; }
    virtual CodeState getCodeState() const        { return CODE_STATE_STABLE;   }
    virtual bool usesExpensiveComputation() const { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Resizes the input volume to the specified dimensions by using a selectable filtering mode.");
    }

    virtual void process();

private:
    void resampleVolume();
    void adjustDimensionProperties();

    void forceUpdate();
    void dimensionsChanged(int dim);
    void allowUpsamplingChanged();
    void keepVoxelRatioChanged();

    VolumePort inport_;
    VolumePort outport_;

    BoolProperty enableProcessingProp_;
    BoolProperty allowUpsampling_;
    BoolProperty keepVoxelRatio_;
    IntProperty resampleDimensionX_;
    IntProperty resampleDimensionY_;
    IntProperty resampleDimensionZ_;
    StringOptionProperty filteringMode_;

    /// Read-only properties displaying the data size of the resampled volume in MB/voxels
    IntProperty outputSizeVoxels_;
    IntProperty outputSizeMB_;

    bool forceUpdate_;
    bool blockSignals_;

    static const std::string loggerCat_; ///< category used in logging
};

}   //namespace

#endif
