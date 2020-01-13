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

#ifndef VRN_PATCHTRAININGDATAEXTRACTOR_H
#define VRN_PATCHTRAININGDATAEXTRACTOR_H

#include "patchfeatureextractor.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/volumeinfoproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/stringproperty.h"

namespace voreen {

class VRN_CORE_API PatchTrainingDataExtractor : public PatchFeatureExtractor {
public:
    PatchTrainingDataExtractor();

    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "PatchTrainingDataExtractor";      }
    virtual std::string getCategory() const   { return "Classification";          }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }
    virtual bool isEndProcessor() const       { return true;              }

    virtual bool usesExpensiveComputation() const { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Exports Patch Training data for Caffe in HDF5 format.");
    }

    virtual void process();

    virtual void computePatches();

    virtual void deactivateProperties();

    virtual void reactivateProperties();

private:

    // data input
    VolumeListPort dataInport_;

    // input of labels
    VolumeListPort labelInport_;
    
    // for scaling intensities
    VolumePort intensityReferenceVolume_;

    IntProperty numSamplesPerFile_;

    ButtonProperty computeButton_;

    StringProperty timeEstimation_;

    bool startComputation_;

    static const std::string loggerCat_;

    FileDialogProperty filename_;

    IntProperty foregroundModulus_;
    BoolProperty capBackground_;

    BoolProperty autoCompute_;
};

} // namespace
#endif 
