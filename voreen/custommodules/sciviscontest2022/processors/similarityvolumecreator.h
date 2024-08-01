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

#ifndef VRN_SIMILARITYVOLUMECREATOR_H
#define VRN_SIMILARITYVOLUMECREATOR_H

#include "voreen/core/processors/asynccomputeprocessor.h"
#include "modules/ensembleanalysis/ports/ensembledatasetport.h"
#include "voreen/core/ports/volumeport.h"
#include "modules/plotting/ports/plotport.h"
#include <Eigen/Eigenvalues>

#include "voreen/core/properties/numeric/intervalproperty.h"

namespace voreen {

    enum MDSAlgorithms {
        CLASSICAL_MDS,
        LANDMARK_MDS
    };

    struct SimilarityVolumeCreatorInput {
        PortDataPointer<EnsembleDataset> ensemble;
        std::unique_ptr<VolumeRAM> outputVolume;
        tgt::Bounds bounds;
        std::string field;
        MDSAlgorithms algorithm;
        int k;
        int iterations;
    };

    struct SimilarityVolumeCreatorOutput {
        std::unique_ptr<VolumeBase> volume;
        std::vector<float> eigenvalues;
    };

/**
 * Creates an RGB volume with a similarity image as defined in the EuroVis 2021 Paper
 */
    class VRN_CORE_API SimilarityVolumeCreator : public AsyncComputeProcessor<SimilarityVolumeCreatorInput, SimilarityVolumeCreatorOutput> {
    public:
        SimilarityVolumeCreator();
        virtual Processor* create() const;

        virtual std::string getClassName() const    { return "SimilarityVolumeCreator";      }
        virtual std::string getCategory() const     { return "Data Processing";  }
        virtual CodeState getCodeState() const      { return CODE_STATE_EXPERIMENTAL;  }

        virtual ComputeInput prepareComputeInput();
        virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
        virtual void processComputeOutput(ComputeOutput output);

    protected:

    private:
        /**
         * Adapts properties to input data
         */
        void adjustToEnsemble();
        /**
         * Performs a classical MDS to 3D with iterative eigenvalue determination
         * based on the correlation matrix and stores the result as CIELAB colors (transformed to RGB) in volume
         * @param correlations Correlations among all pairwise samples
         * @param iterations Number of iterations for eigenvalue decomposition
         * @param output Output similarity volume
         */
        std::vector<float> mds(Eigen::MatrixXf  &distances, int numIterations, VolumeRAM &output) const;
        /**
         * Approximates MDS
         * @param correlations Correlations among all pairwise samples
         * @param iterations Number of iterations for eigenvalue decomposition
         * @param output Output similarity volume
         */
        std::vector<float> landmarkMDS(Eigen::MatrixXf  &distances, int numIterations, int k, VolumeRAM &output) const;
        /**
         * Converts the resulting projection to CIELAB colors and stores them in volume
         * @param result MDS result
         * @param output Output volume to store the colors
         */
        void processResult(Eigen::MatrixXf result, VolumeRAM &output) const;
        void mds(Eigen::MatrixXf &distances, int numIterations, Eigen::MatrixXf &result,
                 std::vector<float> &eigenvalues, Eigen::MatrixXf &eigenvecs) const;
        tgt::vec3 cielabToRGB(tgt::vec3 labVector, float min, float max) const;
        // Ports
        EnsembleDatasetPort inport_;
        VolumePort outport_;
        /// Plotport used to output eigenvalues.
        PlotPort eigenValueOutport_;

        // Properties
        FloatIntervalProperty xRange_;
        FloatIntervalProperty yRange_;
        FloatIntervalProperty zRange_;
        // TODO: Change to Stringlistproperty to allow multiple fields (maybe, one day)
        OptionProperty<int> member_;
        OptionProperty<std::string> field_;
        OptionProperty<MDSAlgorithms> algorithm_;
        IntProperty landmarkNumber_;
        IntProperty iterations_;
    };

}


#endif //VRN_SIMILARITYVOLUMECREATOR_H
