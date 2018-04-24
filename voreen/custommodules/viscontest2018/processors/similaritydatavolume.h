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

#ifndef VRN_SIMILARITYDATAVOLUME_H
#define VRN_SIMILARITYDATAVOLUME_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/ports/volumeport.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/optionproperty.h"

#include "../ports/ensembledatasetport.h"
#include "../properties/stringlistproperty.h"

namespace voreen {


struct SimilarityDataVolumeCreatorInput {
    const EnsembleDataset& dataset;
    VolumeRAM_Float* volumeData;
    float resampleFactor;
    float time;
    std::string runGroup1;
    std::string runGroup2;
    std::string channel;
};

struct SimilarityDataVolumeCreatorOutput {
    VolumeBase* volume;
};

/**
 *
 */
class VRN_CORE_API SimilartyDataVolume : public AsyncComputeProcessor<SimilarityDataVolumeCreatorInput, SimilarityDataVolumeCreatorOutput>  {
public:
    SimilartyDataVolume();
    virtual ~SimilartyDataVolume();
    virtual Processor* create() const;

    virtual std::string getClassName() const      { return "SimilartyDataVolume"; }
    virtual std::string getCategory() const       { return "Plotting";                }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;   }

protected:

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

    void adjustToEnsemble();
    float calculateVariance(const std::vector<float>& voxelData) const;
    float calculateMinMaxDiff(const std::vector<float>& voxelData) const;
    const std::vector<float> applyGroupLogic(const std::vector<float>& rawVoxelData) const;

protected:

    virtual void setDescriptions() {
        setDescription("");
    }

    EnsembleDatasetPort inport_;
    VolumePort outport_;

    FloatProperty resampleFactor_;
    FloatProperty time_;

    StringOptionProperty similarityMethod_;
    StringOptionProperty selectedChannel_;

    StringListProperty group1_;
    StringListProperty group2_;

    static const std::string loggerCat_;
};

} // namespace

#endif
