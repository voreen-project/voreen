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

#ifndef VRN_ENSEMBLEVARIANCEANALYSIS_H
#define VRN_ENSEMBLEVARIANCEANALYSIS_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/ports/volumeport.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/optionproperty.h"

#include "../ports/ensembledatasetport.h"

namespace voreen {


struct EnsembleVarianceAnalysisInput {
    PortDataPointer<EnsembleDataset> ensemble;
    const VolumeBase* meanVolume;
    std::unique_ptr<VolumeRAM_Float> outputVolume;
    std::string field;
    float vectorMagnitudeThreshold;
    int vectorComponent;
    float time;
};

struct EnsembleVarianceAnalysisOutput {
    std::unique_ptr<VolumeBase> volume;
};

/**
 * Calculates ensemble variance for scalar and vector fields as described in \"Interactive Visual
 * Similarity Analysis of Measured and Simulated Multi-field Tubular Flow Ensembles\" by Leistikow et al.
 */
class VRN_CORE_API EnsembleVarianceAnalysis : public AsyncComputeProcessor<EnsembleVarianceAnalysisInput, EnsembleVarianceAnalysisOutput>  {
public:
    EnsembleVarianceAnalysis();
    virtual ~EnsembleVarianceAnalysis();
    virtual Processor* create() const;

    virtual std::string getClassName() const      { return "EnsembleVarianceAnalysis"; }
    virtual std::string getCategory() const       { return "Ensemble Processing";      }
    virtual CodeState getCodeState() const        { return CODE_STATE_TESTING;         }

protected:

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

    void adjustToEnsemble();

protected:

    virtual void setDescriptions() {
        setDescription("Calculates ensemble variance for scalar and vector fields as described in \"Interactive Visual "
                       "Similarity Analysis of Measured and Simulated Multi-field Tubular Flow Ensembles\" by Leistikow et al. ");
        ensembleMeanPort_.setDescription("Mean volume as calculated by <br>EnsembleMeanCreator</br>");
        selectedField_.setDescription("Link with <br>EnsembleMeanCreator</br>");
        time_.setDescription("Link with <br>EnsembleMeanCreator</br>");
        vectorComponent_.setDescription("Defines if magnitude, direction or both should be used for calculating vector variance");
        vectorMagnitudeThreshold_.setDescription("Only calculate vector variance based on direction if magnitude "
                                                 "surpasses this threshold (direction calculation of very small vectors "
                                                 "is unstable.");
    }

    enum VectorComponent {
        BOTH = 0,
        MAGNITUDE = 1,
        DIRECTION = 2
    };

    EnsembleDatasetPort ensembleInport_;
    VolumePort ensembleMeanPort_;
    VolumePort outport_;

    StringOptionProperty selectedField_;
    FloatProperty vectorMagnitudeThreshold_;
    IntOptionProperty vectorComponent_;
    FloatProperty time_;

    static const std::string loggerCat_;
};

} // namespace

#endif
