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

#ifndef VRN_SIMILARITYMATRIXCOMBINE_H
#define VRN_SIMILARITYMATRIXCOMBINE_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/utils/statistics.h"

#include "../ports/ensembledatasetport.h"
#include "../ports/similaritymatrixport.h"

namespace voreen {

enum SimilarityCombinationMethod {
    METHOD_MAX,
    METHOD_AVG,
    METHOD_STD,
    METHOD_MULTIPLY,
};

struct SimilarityMatrixCombineInput {
    std::vector<const SimilarityMatrixList*> inputMatrixLists;
    std::unique_ptr<SimilarityMatrixList> outputMatrices;
    SimilarityCombinationMethod method;
};

struct SimilarityMatrixCombineOutput {
    std::unique_ptr<SimilarityMatrixList> outputMatrices;
};

/**
 * This processor combines two equally shaped similarity matrices as described in "Interactive Visual Similarity Analysis of Measured and Simulated Multi-field Tubular Flow Ensembles" by Leistikow et al.
 */
class VRN_CORE_API SimilarityMatrixCombine : public AsyncComputeProcessor<SimilarityMatrixCombineInput, SimilarityMatrixCombineOutput> {
public:
    SimilarityMatrixCombine();

    virtual Processor* create() const;
    virtual std::string getClassName() const        { return "SimilarityMatrixCombine"; }
    virtual std::string getCategory() const         { return "Plotting";                }
    virtual CodeState getCodeState() const          { return CODE_STATE_TESTING;        }

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

protected:

    virtual bool isReady() const;
    virtual void setDescriptions() {
        setDescription("Combines two equally shaped similarity matrices as described in \"Interactive Visual Similarity Analysis of Measured and Simulated Multi-field Tubular Flow Ensembles\" by Leistikow et al. ");
        similarityCombinationMethod_.setDescription("Defines which approach is to be used for combining the input matrices");
        ignoreHash_.setDescription("It is possible to combine two matrices of the same shape for different ensembles. If it's ensured that combining these different ensembles is semantically useful, tick this box to ignore the respective test.");
    }

private:

    OptionProperty<SimilarityCombinationMethod> similarityCombinationMethod_;
    BoolProperty ignoreHash_;

    /// Inport for the ensemble data structure.
    SimilarityMatrixPort inport_;

    /// Plotport used to output eigenvalues.
    SimilarityMatrixPort outport_;

    static const std::string loggerCat_;
};


} // namespace

#endif
