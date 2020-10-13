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

#ifndef VRN_SIMILARITYMATRIXCREATOR_H
#define VRN_SIMILARITYMATRIXCREATOR_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/optionproperty.h"

#include "../ports/ensembledatasetport.h"
#include "../ports/similaritymatrixport.h"

namespace voreen {

enum SingleChannelSimilarityMeasure {
    MEASURE_ISOCONTOURS,
    MEASURE_GENERALIZED,
    MEASURE_AVG_DIFFERENCE,
};

enum MultiChannelSimilarityMeasure {
    MEASURE_MAGNITUDE,
    MEASURE_ANGLE,
    MEASURE_JIANG,
    MEASURE_CROSSPRODUCT,
    MEASURE_SPLIT_CHANNELS,
    MEASURE_VECTOR_DIFFERENCE_MAGNITUDE
};

struct SimilarityMatrixCreatorInput {
    PortDataPointer<EnsembleDataset> ensemble;
    std::unique_ptr<SimilarityMatrixList> outputMatrices;
    std::vector<tgt::vec3> seedPoints;
    SingleChannelSimilarityMeasure singleChannelSimilarityMeasure;
    float isoValue;
    MultiChannelSimilarityMeasure multiChannelSimilarityMeasure;
    float weight;
    std::string hash;
};

struct SimilarityMatrixCreatorOutput {
    std::unique_ptr<SimilarityMatrixList> outputMatrices;
};

class VRN_CORE_API SimilarityMatrixCreator : public AsyncComputeProcessor<SimilarityMatrixCreatorInput, SimilarityMatrixCreatorOutput> {
public:
    SimilarityMatrixCreator();

    virtual Processor* create() const;
    virtual std::string getClassName() const        { return "SimilarityMatrixCreator"; }
    virtual std::string getCategory() const         { return "Plotting";                }
    virtual CodeState getCodeState() const          { return CODE_STATE_EXPERIMENTAL;   }

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

protected:

    virtual bool isReady() const;
    virtual void setDescriptions();
    virtual void adjustPropertiesToInput();
    virtual std::vector<std::reference_wrapper<Port>> getCriticalPorts();

private:

    std::string calculateHash() const;

    OptionProperty<SingleChannelSimilarityMeasure> singleChannelSimilarityMeasure_;
    FloatProperty isoValue_;
    OptionProperty<MultiChannelSimilarityMeasure> multiChannelSimilarityMeasure_;
    FloatProperty weight_;
    IntProperty numSeedPoints_;
    IntProperty seedTime_;
    ButtonProperty clearCache_;

    /// Inport for the ensemble data structure.
    EnsembleDatasetPort inport_;

    /// Inport for seed mask.
    VolumePort seedMask_;

    /// Plotport used to output eigenvalues.
    SimilarityMatrixPort outport_;

    static const std::string loggerCat_;
};



} // namespace

#endif
