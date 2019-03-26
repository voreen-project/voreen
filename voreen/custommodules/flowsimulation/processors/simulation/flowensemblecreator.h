/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_FLOWENSEMBLECREATOR_H
#define VRN_FLOWENSEMBLECREATOR_H

#include "voreen/core/processors/asynccomputeprocessor.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/vectorproperty.h"

#include "custommodules/bigdataimageprocessing/properties/interactivelistproperty.h"
#include "custommodules/bigdataimageprocessing/volumefiltering/slicereader.h"

namespace voreen {

class SliceReader;
class VolumeFilter;
class FlowFeature;

struct FlowEnsembleCreatorInput {
    std::string simulationResultPath;
    std::string ensembleOutputPath;
    const VolumeList* measuredData;
};

struct FlowEnsembleCreatorOutput {
};

/**
 * Creates an ensemble dataset with disk representations (HDF5) from raw simulated data (local or cluster)
 * and adds measured data. The created ensemble can be used by the ensemble analysis module.
 */
class VRN_CORE_API FlowEnsembleCreator : public AsyncComputeProcessor<FlowEnsembleCreatorInput, FlowEnsembleCreatorOutput> {
public:
    FlowEnsembleCreator();
    virtual ~FlowEnsembleCreator();

    Processor* create() const;

    std::string getClassName() const { return "FlowEnsembleCreator"; }

    std::string getCategory() const { return "Volume Processing"; }

    CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }

    virtual bool isReady() const;

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

protected:
    void setDescriptions() {
        setDescription(
                "Creates an ensemble dataset with disk representations (HDF5) from raw simulated data (local or cluster)"
                "and adds measured data. The created ensemble can be used by the ensemble analysis module."
        );
    }

    virtual void adjustPropertiesToInput();

private:

    void addFeature(FlowFeature* filterProperties);

    /// Filter list
    std::vector<std::unique_ptr<FlowFeature>> flowFeatures_;

    VolumeListPort inport_;

    // Input properties.
    IntProperty spatialResolution_;
    FloatVec3Property spacing_;
    FloatVec3Property offset_;

    // Output properties.
    FileDialogProperty simulationResultPath_;
    FileDialogProperty ensembleOutputPath_;
    IntProperty outputVolumeDeflateLevel_;
    BoolProperty deleteOriginalData_;

    InteractiveListProperty featureList_;

    static const std::string loggerCat_; ///< category used in logging
};

}

#endif // VRN_FLOWENSEMBLECREATOR_H
