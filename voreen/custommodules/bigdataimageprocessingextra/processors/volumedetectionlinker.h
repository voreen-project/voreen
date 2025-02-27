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

#ifndef VRN_VOLUMEDETECTIONLINKER_H
#define VRN_VOLUMEDETECTIONLINKER_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/buttonproperty.h"

namespace voreen {

/**
 * Allows the quantification of the volume of a single (binary) segmentation or the quantification of two segmentationi volumes (including the density of one within the other).
 */
class VRN_CORE_API VolumeDetectionLinker : public Processor {

public:
    VolumeDetectionLinker();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "VolumeDetectionLinker";     }
    virtual std::string getCategory() const  { return "Quantification";         }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL;  }

    virtual bool isEndProcessor() const {return true;   }

protected:
    virtual void setDescriptions() {
        setDescription("Processor that simply changes a property state when a number of volumes are connected to the inport");
    }

    virtual void process();

    BoolProperty volumesDetected_;
    IntProperty numVolumesRequiredForDetection_;
    IntProperty numDetectionsRequired_;
    IntProperty numDetections_;
    ButtonProperty reset_;

    VolumePort volumePort_;

    static const std::string loggerCat_;

private:
    void checkCondition();

};

} // namespace

#endif
