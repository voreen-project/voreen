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

struct ParallelCoordinatesAxesCreatorInput {
    std::unique_ptr<EnsembleDataset> ensemble;
    std::string ensembleHash;
    int temporalSampleCount;
    bool aggregate;
    std::vector<tgt::vec3> seedPoints;
    std::vector<std::pair<std::string, int>> fields;
    std::vector<std::string> axesLabels;
    std::vector<tgt::vec2> ranges;
};

struct ParallelCoordinatesAxesCreatorOutput {
    std::unique_ptr<ParallelCoordinatesAxes> axes;
};

class ParallelCoordinatesAxesCreator : public AsyncComputeProcessor<ParallelCoordinatesAxesCreatorInput, ParallelCoordinatesAxesCreatorOutput> {
public:
    ParallelCoordinatesAxesCreator();

    Processor* create() const override;
    std::string getClassName() const override { return "ParallelCoordinatesAxesCreator"; }
    std::string getCategory() const override { return "ParallelCoordinates"; }
    CodeState getCodeState() const override { return CODE_STATE_TESTING; }
    bool isReady() const override;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

private:

    void setDescriptions() {
        setDescription("Creates parallel coordinates for a given ensemble dataset.<br>"
                       "The ensemble must contain more than a single field.");

        propertyMembers_.setDescription("Use to select considered members from the ensemble.");
        propertyFields_.setDescription("Use to select considered fields from the ensemble.");
        propertySpatialSampleCount_.setDescription("Number of spatial samples");
        propertyTemporalSampleCount_.setDescription("Number of temporal samples");
        propertyAggregateMembers_.setDescription("If enabled, values from selected members will be aggregated");
    }

    ComputeInput prepareComputeInput() override;
    ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const override;
    void processComputeOutput(ComputeOutput output) override;

    std::vector<std::reference_wrapper<Port>> getCriticalPorts() override;

    void adjustToEnsemble();
    void adjustToSelection();

    EnsembleDatasetPort ensembleport_;
    VolumePort seedMask_;
    ParallelCoordinatesAxesPort axesport_;

    StringListProperty propertyMembers_;
    StringListProperty propertyFields_;
    IntProperty propertySpatialSampleCount_;
    IntProperty propertyTemporalSampleCount_;
    StringOptionProperty propertySampleRegion_;
    IntProperty propertySeedTime_;
    BoolProperty propertyAggregateMembers_;

    /// Hash value of last valid data.
    std::string hash_;

    /// Common fields of last member selection.
    std::vector<std::pair<std::string, int>> fields_;

    /// ..and their respective labels.
    std::vector<std::string> axesLabels_;
};
}

#endif // VRN_PARALLELCOORDINATESAXESCREATOR_H
