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

#ifndef VRN_PARALLELCOORDINATESVOXELSELECTION_H
#define VRN_PARALLELCOORDINATESVOXELSELECTION_H

#include "voreen/core/processors/asynccomputeprocessor.h"
#include "voreen/core/ports/volumeport.h"

#include "../ports/ensembledatasetport.h"
#include "../properties/parallelcoordinatessectionsproperty.h"

namespace voreen {

struct ParallelCoordinatesVoxelSelectionInput {
    PortDataPointer<EnsembleDataset> ensemble;
    tgt::Bounds bounds;
    ParallelCoordinatesSectionsPropertyData sectionData;
    std::vector<std::pair<const VolumeBase*, int>> inputVolumes;
    std::unique_ptr<VolumeRAM_UInt8> outputVolume;
};

struct ParallelCoordinatesVoxelSelectionOutput {
    std::unique_ptr<VolumeBase> volume;
};

/**
 * This processor is to be used with a ParallelCoordinatesViewer (link "sections" property).
 * It creates a binary volume mask according to the selection made in ParallelCoordinatesViewer.
 */
class ParallelCoordinatesVoxelSelection : public AsyncComputeProcessor<ParallelCoordinatesVoxelSelectionInput, ParallelCoordinatesVoxelSelectionOutput> {
public:
    ParallelCoordinatesVoxelSelection();

    virtual Processor* create() const override { return new ParallelCoordinatesVoxelSelection(); }
    virtual std::string getClassName() const override { return "ParallelCoordinatesVoxelSelection"; }
    virtual std::string getCategory() const override { return "ParallelCoordinates"; }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL; }

private:

    virtual void setDescriptions() {
        setDescription("This processor needs to be linked with a <br>ParallelCoordinatesViewer</br> processor. "
                       "It creates a mask where only those voxels are non-zero where the respective value of each field "
                       "is inside all respective ranges selected in the ParallelCoordinatesViewer.");
        propertySections_.setDescription("Link with <br>ParallelCoordinatesViewer</br>");
        outputDimensions_.setDescription("Resolution of the mask");
    }

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

    EnsembleDatasetPort ensembleport_;
    VolumePort volumeport_;

    ParallelCoordinatesSectionsProperty propertySections_;

    IntVec3Property outputDimensions_;

    static const std::string loggerCat_;
};

}

#endif // VRN_PARALLELCOORDINATESVOXELSELECTION_H
