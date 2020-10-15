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

#ifndef VRN_VORTEXPROCESSOR_H
#define VRN_VORTEXPROCESSOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

class VortexProcessor : public Processor {
public:
    VortexProcessor();

    Processor* create() const override
    {
        return new VortexProcessor();
    }
    std::string getClassName() const override
    {
        return "VortexProcessor";
    }
    std::string getCategory() const override
    {
        return "Vortex Extraction";
    }
    bool isReady() const override
    {
        return _inportVelocity.isReady() && _inportMask.isReady();
    }

    static void Process( const VolumeRAM& velocity, const VolumeRAM& mask, VolumeRAM_Mat3Float& outJacobi, VolumeRAM_Float* outDelta, VolumeRAM_Float* outLambda2, VolumeRAM_Float* outQ );

private:
    void process() override;

    VolumePort _inportVelocity, _inportMask;
    VolumePort _outportJacobi, _outportDelta, _outportLambda2, _outportQ;
};

}

#endif // VRN_VORTEXPROCESSOR_H