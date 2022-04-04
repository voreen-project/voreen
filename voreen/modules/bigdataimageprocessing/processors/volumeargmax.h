/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_VOLUMEARGMAX_H
#define VRN_VOLUMEARGMAX_H

#include "voreen/core/processors/asynccomputeprocessor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/temppathproperty.h"

#include "modules/hdf5/io/hdf5filevolume.h"

#include <string>

namespace voreen {

struct VolumeArgMaxInput {
    const VolumeBase& vol0_;
    const VolumeBase* vol1_;
    const VolumeBase* vol2_;
    const VolumeBase* vol3_;
    std::unique_ptr<HDF5FileVolume> outputVolume_;
};

struct VolumeArgMaxOutput {
    std::string outputVolumePath_;
};

class VolumeArgMax : public AsyncComputeProcessor<VolumeArgMaxInput, VolumeArgMaxOutput> {
public:
    VolumeArgMax();
    virtual ~VolumeArgMax();
    virtual Processor* create() const;

    virtual std::string getCategory() const             { return "Volume Processing"; }
    virtual std::string getClassName() const            { return "VolumeArgMax";      }
    virtual Processor::CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL;  }

    virtual bool isReady() const;

    static const std::string loggerCat_; ///< category used in logging

protected:
    virtual void setDescriptions() {
        setDescription("Sets each voxel to the id of the input volume with the highest voxel value");
    }

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);
private:
    VolumePort inportVolume0_;
    VolumePort inportVolume1_;
    VolumePort inportVolume2_;
    VolumePort inportVolume3_;
    VolumePort outportIds_;

    ButtonProperty clearResult_;
    TempPathProperty outputVolumePath_;
};

} //namespace

#endif
