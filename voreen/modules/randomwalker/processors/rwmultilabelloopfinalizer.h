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

#ifndef VRN_RWMULTILABELLOOPFINALIZER_H
#define VRN_RWMULTILABELLOOPFINALIZER_H

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumelist.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "voreen/core/ports/loopport.h"
#include "voreen/core/properties/optionproperty.h"


#include <string>
#include "tgt/vector.h"


namespace voreen {

class Volume;

class RWMultiLabelLoopFinalizer : public VolumeProcessor {

public:
    RWMultiLabelLoopFinalizer();
    virtual ~RWMultiLabelLoopFinalizer();
    virtual Processor* create() const;

    virtual std::string getCategory() const { return "Volume Processing"; }
    virtual std::string getClassName() const { return "RWMultiLabelLoopFinalizer"; }
    virtual CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Loop-wrapper for the RandomWalker processor for multi-label segmentations (experimental). <br/>"
                       "See: RWMultiLabelLoopInitializer");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

private:
    VolumePort probabilitiesInport_;
    VolumePort segmentationInport_;
    VolumeListPort probabilitiesCollectionOutport_;
    VolumeListPort segmentationCollectionOutport_;
    VolumePort segmentationOutport_;
    VolumePort segmentationProbOutport_;

    LoopPort loopOutport_;

    VolumeList* segmentationCollection_;
    VolumeList* probabilitiesCollection_;
    VolumeRAM_UInt8* segmentationVolume_;
    VolumeRAM_UInt16* segmentationProbVolume_;

    clock_t runtime_;

    static const std::string loggerCat_; ///< category used in logging
};

} //namespace

#endif
