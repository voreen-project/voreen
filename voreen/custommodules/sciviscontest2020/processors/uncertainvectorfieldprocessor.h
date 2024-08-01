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

#ifndef VRN_UNCERTAINVECTORFIELDPROCESSOR_H
#define VRN_UNCERTAINVECTORFIELDPROCESSOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/string/stringlistproperty.h"
#include "voreen/core/ports/volumeport.h"

#include "modules/ensembleanalysis/ports/ensembledatasetport.h"

namespace voreen {

class UncertainVectorFieldProcessor : public Processor {
public:
    UncertainVectorFieldProcessor();

    Processor* create() const override
    {
        return new UncertainVectorFieldProcessor();
    }
    std::string getClassName() const override
    {
        return "UncertainVectorFieldProcessor";
    }
    std::string getCategory() const override
    {
        return "Volume Processing";
    }

private:
    void process() override {}
    void update();

    EnsembleDatasetPort _inportEnsemble;
    VolumePort _inportMask, _outportQ, _outportL;

    StringListProperty _propertySelectedMembers;
    IntProperty _propertyTimestep, _propertySampleCount;
    FloatProperty _propertyThresholdQ, _propertyThresholdL;
    ButtonProperty _propertyUpdate;
};

}

#endif // VRN_UNCERTAINVECTORFIELDPROCESSOR_H