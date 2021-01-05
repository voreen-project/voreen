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

#ifndef VRN_HELICITYDENSITY_H
#define VRN_HELICITYDENSITY_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/volumeport.h"

namespace voreen {

/**
 * This processor calculates the helicity density for a velocity vector field and its vorticity.
 */
class VRN_CORE_API HelicityDensity : public Processor {
public:
    HelicityDensity();
    virtual ~HelicityDensity();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "HelicityDensity"; }
    virtual std::string getCategory() const   { return "Volume Processing"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_TESTING; }

protected:
    virtual void setDescriptions() {
        setDescription("Calculates the helicity density for a velocity vector field and its vorticity. "
                       "Helicity density is defined by the voxel-wise dot-product of velocity and vorticity.");
        velocityInport_.setDescription("Velocity input volume (3D vector field)");
        vorticityInport_.setDescription("Vorticity input volume (3D vector field)");
    }

    virtual bool isReady() const;
    virtual void process();

private:

    VolumePort velocityInport_;
    VolumePort vorticityInport_;
    VolumePort helicityDensityOutport_;
    BoolProperty normalize_;
};

}   //namespace

#endif
