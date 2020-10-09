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

#ifndef VRN_PARALLELCOORDINATESAXESCREATOR_H
#define VRN_PARALLELCOORDINATESAXESCREATOR_H

#include "voreen/core/processors/asynccomputeprocessor.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/string/stringlistproperty.h"

#include "../ports/ensembledatasetport.h"
#include "../ports/parallelcoordinatesaxesport.h"

namespace voreen {

struct ParallelCoordianesAxesCreatorInput {
    EnsembleDataset ensemble;
    int sampleCount;
    size_t numTimeSteps;
    std::vector<int32_t> validVoxels;
    std::vector<std::string> fieldNames;
    std::vector<std::pair<float, float>> ranges;
};

struct ParallelCoordianesAxesCreatorOutput {
    std::unique_ptr<ParallelCoordinatesAxes> axes;
};

class ParallelCoordinatesAxesCreator : public AsyncComputeProcessor<ParallelCoordianesAxesCreatorInput, ParallelCoordianesAxesCreatorOutput> {
public:
    ParallelCoordinatesAxesCreator();

    Processor* create() const override;
    std::string getClassName() const override;
    std::string getCategory() const override;
    bool isReady() const override;

private:

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

    virtual std::vector<std::reference_wrapper<Port>> getCriticalPorts();

    void updateValidVoxels();

    EnsembleDatasetPort ensembleport_;
    VolumePort volumeport_;
    ParallelCoordinatesAxesPort axesport_;

    StringListProperty propertyMembers_;
    StringListProperty propertyFields_;
    IntProperty propertySampleCount_;
    BoolProperty propertyAggregateMembers_;
    FileDialogProperty propertyFileDialog_;
    ButtonProperty propertySaveButton_;

    std::vector<int32_t> validVoxels_;
};
}

#endif // VRN_PARALLELCOORDINATESAXESCREATOR_H
