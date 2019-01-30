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

#ifndef VRN_VOLUMEMERGER_H
#define VRN_VOLUMEMERGER_H

#include "voreen/core/processors/asynccomputeprocessor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/temppathproperty.h"

#include "modules/hdf5/io/hdf5filevolume.h"

namespace voreen {

struct VolumeMergerComputeInput {
    const VolumeList* inputVolumes;
    std::unique_ptr<Volume> outputVolume;
};

struct VolumeMergerComputeOutput{
    std::unique_ptr<Volume> outputVolume;
};

/**
 *
 */
class VRN_CORE_API VolumeMerger : public AsyncComputeProcessor<VolumeMergerComputeInput, VolumeMergerComputeOutput>  {
public:
    VolumeMerger();
    virtual ~VolumeMerger();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "VolumeMerger";      }
    virtual std::string getCategory() const   { return "Volume Processing";     }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }

protected:
    virtual void setDescriptions() {
        setDescription("Merges incoming volumes into a single one. No overlap allowed");
    }

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

private:

    VolumeListPort inport_;
    VolumePort outport_;

    BoolProperty allowIntersections_;
};

}   //namespace

#endif
