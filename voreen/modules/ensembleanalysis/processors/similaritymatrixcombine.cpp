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

#include "similaritymatrixcombine.h"

#include "voreen/core/voreenapplication.h"

namespace voreen {

const std::string SimilarityMatrixCombine::loggerCat_("voreen.ensembleanalysis.SimilarityMatrixCombine");

SimilarityMatrixCombine::SimilarityMatrixCombine()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , inport_(Port::INPORT, "inport", "Ensemble Datastructure Input", true)
    , outport_(Port::OUTPORT, "outport", "Similarity Matrix Output", false)
    , similarityCombinationMethod_("similarityCombinationMethod", "Combination Method")
    , ignoreHash_("ignoreHash", "Ignore Hash", false)
{
    // Ports
    addPort(inport_);
    addPort(outport_);

    // Calculation
    addProperty(similarityCombinationMethod_);
    similarityCombinationMethod_.addOption("method-max", "Max.", METHOD_MAX);
    similarityCombinationMethod_.addOption("method-avg", "Avg.", METHOD_AVG);
    //similarityCombinationMethod_.addOption("method-std", "Std.", METHOD_STD); // No proper use-case yet-
    similarityCombinationMethod_.addOption("method-multiply", "Multiply", METHOD_MULTIPLY);
    similarityCombinationMethod_.selectByValue(METHOD_MAX);

    addProperty(ignoreHash_);
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
        if(!ignoreHash_.get() && referenceMatrices->getHash() != matrix->getHash()) {
            throw InvalidInputException("Hashes do not match", InvalidInputException::S_ERROR);
        }

        if (referenceMatrices->getFieldNames().size() != matrix->getFieldNames().size()) {
            throw InvalidInputException("Number of contained fields does not match", InvalidInputException::S_ERROR);
        }

        if (referenceMatrices->getSize() != matrix->getSize()) {
            throw InvalidInputException("Matrix size does not match", InvalidInputException::S_ERROR);
        }

        inputMatrixLists.push_back(matrix);
    }

    std::unique_ptr<SimilarityMatrixList> outputMatrices(new SimilarityMatrixList(*referenceMatrices));

    return SimilarityMatrixCombineInput{
            inputMatrixLists,
            std::move(outputMatrices),
            similarityCombinationMethod_.getValue()
    };
}

SimilarityMatrixCombineOutput SimilarityMatrixCombine::compute(SimilarityMatrixCombineInput input, ProgressReporter& progress) const {

    std::vector<const SimilarityMatrixList*> inputMatrixLists = input.inputMatrixLists;
    std::unique_ptr<SimilarityMatrixList> outputMatrices = std::move(input.outputMatrices);

    std::vector<std::string> fieldNames = outputMatrices->getFieldNames();
    for(size_t fi = 0; fi < fieldNames.size(); fi++) {
        const std::string& fieldName = fieldNames[fi];

        SimilarityMatrix& outputDistanceMatrix = outputMatrices->getSimilarityMatrix(fieldName);
        long size = static_cast<long>(outputDistanceMatrix.getSize());

        SubtaskProgressReporter flagsProgress(progress,tgt::vec2(fi, fi+1) / tgt::vec2(fieldNames.size()));
        ThreadedTaskProgressReporter threadedProgress(flagsProgress, size);
        bool aborted = false;

#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for(long i=0; i<size; i++) {
#ifdef VRN_MODULE_OPENMP
            if(aborted) {
                continue;
            }
#else
            if(aborted) {
                break;
            }
#endif

            for(long j=0; j<=i; j++) {

                float d = 1.0f;
                Statistics statistics(false);

                for (const SimilarityMatrixList* inputMatrixList : inputMatrixLists) {
                    std::string inputFieldName = inputMatrixList->getFieldNames()[fi]; // Never ever use reference here!
                    const SimilarityMatrix& inputDistanceMatrix = inputMatrixList->getSimilarityMatrix(inputFieldName);
                    // We just need one or the other measure, but simply do both, since the calculation is cheap.
                    float dissimilarity = inputDistanceMatrix(i, j);
                    statistics.addSample(dissimilarity);
                    d *= (1.0f - dissimilarity); // Transform to similarity
                }

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
                case METHOD_MULTIPLY:
                    outputDistanceMatrix(i, j) = 1.0f - d; // Transform to dissimilarity.
                    break;
                }
            }

            if (threadedProgress.reportStepDone()) {
#ifdef VRN_MODULE_OPENMP
#pragma omp critical
                aborted = true;
#else
                aborted = true;
                break;
#endif
            }
        }

        if (aborted) {
            throw boost::thread_interrupted();
        }
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
