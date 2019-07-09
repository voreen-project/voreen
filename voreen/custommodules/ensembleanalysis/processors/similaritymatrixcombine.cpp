/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
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

#include "similaritymatrixcombine.h"

#include "voreen/core/voreenapplication.h"

namespace voreen {

const std::string SimilarityMatrixCombine::loggerCat_("voreen.ensembleanalysis.SimilarityMatrixCombine");

SimilarityMatrixCombine::SimilarityMatrixCombine()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , inport_(Port::INPORT, "inport", "Ensemble Datastructure Input", true)
    , outport_(Port::OUTPORT, "outport", "Similarity Matrix Output", false)
    , similarityCombinationMethod_("similarityCombinationMethod", "Combination Method")
{
    // Ports
    addPort(inport_);
    addPort(outport_);

    // Calculation
    addProperty(similarityCombinationMethod_);
    similarityCombinationMethod_.addOption("method-max", "Max.", METHOD_MAX);
    similarityCombinationMethod_.addOption("method-avg", "Avg.", METHOD_AVG);
    similarityCombinationMethod_.addOption("method-std", "Std.", METHOD_STD);
    similarityCombinationMethod_.selectByValue(METHOD_MAX);
}

Processor* SimilarityMatrixCombine::create() const {
    return new SimilarityMatrixCombine();
}


bool SimilarityMatrixCombine::isReady() const {
    if (!isInitialized()) {
        setNotReadyErrorMessage("Not initialized.");
        return false;
    }
    if (!inport_.isReady()) {
        setNotReadyErrorMessage("Inport not ready.");
        return false;
    }

    // Note: Seed Mask is optional!

    return true;
}

SimilarityMatrixCombineInput SimilarityMatrixCombine::prepareComputeInput() {
    const auto& inputData = inport_.getThreadSafeAllData();
    if (inputData.empty())
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);

    const SimilarityMatrixList* referenceMatrices = inputData.front();
    std::vector<const SimilarityMatrixList*> inputMatrixLists;
    for(const auto& matrix : inputData) {
        if(referenceMatrices->getHash() != matrix->getHash()) {
            throw InvalidInputException("Hashes do not match", InvalidInputException::S_ERROR);
        }

        inputMatrixLists.push_back(matrix);
    }

    std::vector<Statistics> statisticsMatrix(referenceMatrices->getSize()*referenceMatrices->getSize(), false);
    std::unique_ptr<SimilarityMatrixList> outputMatrices(new SimilarityMatrixList(*referenceMatrices));

    return SimilarityMatrixCombineInput{
            inputMatrixLists,
            std::move(statisticsMatrix),
            std::move(outputMatrices),
            similarityCombinationMethod_.getValue()
    };
}

SimilarityMatrixCombineOutput SimilarityMatrixCombine::compute(SimilarityMatrixCombineInput input, ProgressReporter& progress) const {

    std::vector<const SimilarityMatrixList*> inputMatrixLists = input.inputMatrixLists;
    std::vector<Statistics> statisticsMatrix = std::move(input.statisticsMatrix);
    std::unique_ptr<SimilarityMatrixList> outputMatrices = std::move(input.outputMatrices);

    std::vector<std::string> fieldNames = outputMatrices->getFieldNames();
    for(size_t ch = 0; ch < fieldNames.size(); ch++) {
        const std::string& fieldName = fieldNames[ch];

        size_t size = outputMatrices->getSimilarityMatrix(fieldName).getSize();
#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for shared(outputMatrices), shared(statisticsMatrix), shared(inputMatrixLists)
#endif
        for(long i=0; i<static_cast<long>(size); i++) {
            for(long j=0; j<=i; j++) {

                Statistics& statistics = statisticsMatrix[i * size + j];

                for (const SimilarityMatrixList* inputMatrixList : inputMatrixLists) {
                    const SimilarityMatrix& inputDistanceMatrix = inputMatrixList->getSimilarityMatrix(fieldName);
                    statistics.addSample(inputDistanceMatrix(i, j));
                }

                SimilarityMatrix& outputDistanceMatrix = outputMatrices->getSimilarityMatrix(fieldName);
                switch (input.method) {
                case METHOD_MAX:
                    outputDistanceMatrix(i, j) = statistics.getMax();
                    break;
                case METHOD_AVG:
                    outputDistanceMatrix(i, j) = statistics.getMean();
                    break;
                case METHOD_STD:
                    outputDistanceMatrix(i, j) = statistics.getStdDev();
                    break;
                }
            }
        }

        progress.setProgress(1.0f * ch / fieldNames.size());
    }

    progress.setProgress(1.0f);

    return SimilarityMatrixCombineOutput{
            std::move(outputMatrices)
    };
}

void SimilarityMatrixCombine::processComputeOutput(ComputeOutput output) {
    outport_.setData(output.outputMatrices.release(), true);
}


} // namespace
