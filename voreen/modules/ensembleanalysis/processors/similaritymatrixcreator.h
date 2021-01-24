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
    MEASURE_EUCLIDEAN_NORM
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

/**
 * This processor Creates a Similarity matrix for each field of the ensemble, containing all pairwise dissimilarities
 * of all time steps of its members.
 */
class VRN_CORE_API SimilarityMatrixCreator : public AsyncComputeProcessor<SimilarityMatrixCreatorInput, SimilarityMatrixCreatorOutput> {
    // Missing values (due to out of bounds sampling, volume could not be loaded) are filled by the following value.
    // TODO: this should either be adjustable by property or handled differently.
    static constexpr float MissingValue = 0.0f;
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
    virtual void adjustPropertiesToInput();
    virtual std::vector<std::reference_wrapper<Port>> getCriticalPorts();

    void setDescriptions() {
        setDescription("Creates a Similarity matrix for each field of the ensemble, containing all pairwise "
                       "dissimilarities of all time steps of its members");
        singleChannelSimilarityMeasure_.setDescription("Specifies how similarity is calculated for single channel fields:<br>"
                                                       "<strong>Iso-Contours</strong>: Uses the isovalue and applies the Jaccard index<br>"
                                                       "<strong>Generalized</strong>: Applies the generalized field similarity by Fofonov et al.<br>"
                                                       "<strong>Avg. Difference</strong>: Calculates the normalized average difference");
        isoValue_.setDescription("Defines the isovalue used for the Iso-Contour similarity measure");
        multiChannelSimilarityMeasure_.setDescription("Specifies how similarity is calculated for multi-channel fields:<br>"
                                                      "<strong>Magnitude</strong>: Uses the vector magnitude and uses the generalized field similarity by Fofonov et al.<br>"
                                                      "<strong>Angle</strong>: Calculates the average normalized angle between vectors<br>"
                                                      "<strong>Jiang</strong>: Uses the vector similarity introduced by Jiang et al.<br>"
                                                      "<strong>Crossproduct Magnitude</strong>: Calculates the cross product magnitude of two vectors<br>"
                                                      "<strong>Split Channels</strong>: Treats each channel as separate scalar field and applies the generalized field similarity by Fofonov et al.<br>"
                                                      "<strong>Euclidean Norm</strong>: Uses the magnitude of the difference between two vectors");
        weight_.setDescription("Defines the weights between magnitude and direction used by similarity measure of Jiang et al.");
        numSeedPoints_.setDescription("Number of seeds to define a Monte Carlo sampling");
        seedTime_.setDescription("Allows to use a fixed seed");
        seedMask_.setDescription("Optional seed mask. Only voxels different from zero are seeding candidates");
    }

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
