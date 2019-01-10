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

#ifndef VRN_TILECAFFECLASSIFIER_H
#define VRN_TILECAFFECLASSIFIER_H

#include "voreen/core/processors/volumeprocessor.h"

#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/volumeinfoproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/stringproperty.h"

namespace voreen {

class VRN_CORE_API TileCaffeClassifier : public VolumeProcessor {
public:
    TileCaffeClassifier();

    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "TileCaffeClassifier";      }
    virtual std::string getCategory() const   { return "Classification";          }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }
    //virtual bool isEndProcessor() const       { return true;              }

    virtual bool usesExpensiveComputation() const { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Extracts 2D tiles(x-y) plus multiple channels in z-direction fomr a volume and uses a trained Caffe network (CNN) to segment these patches (or do any kind of voxel-wise classification) and assemble an output segmentation as a float volume.");
    }

    virtual void process();

    virtual void classify();

    virtual void deactivateProperties();

    virtual void reactivateProperties();

private:

    VolumePort dataInport_;

    VolumePort classifierOutput_;

    IntProperty tileSize_;
    IntProperty borderSize_;
    IntProperty channelWidth_;
    BoolProperty useMirroring_;


    ButtonProperty computeButton_;
    
    StringProperty timeEstimation_;

    bool startComputation_;

    static const std::string loggerCat_;

    FileDialogProperty protoFilename_;
    FileDialogProperty modelFilename_;

    BoolProperty autoCompute_;
};

} // namespace
#endif 