/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2015 University of Muenster, Germany.                        *
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

#ifndef VRN_COSMOLOGYVOLUMECONVERTER_H
#define VRN_COSMOLOGYVOLUMECONVERTER_H

#include "voreen/core/processors/volumeprocessor.h"


#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/ports/volumeport.h"

#include "../ports/cmparticleport.h"

namespace voreen{
class CosmologyVolumeConverter : public CachingVolumeProcessor{
public:
    CosmologyVolumeConverter();
    virtual ~CosmologyVolumeConverter();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "CosmologyVolumeConverter";            }
    virtual std::string getCategory() const  { return "Viscontest201x";                      }
    virtual void        setDescriptions()    { setDescription("Cosmology Volume Converter"); }
    virtual CodeState   getCodeState() const { return CODE_STATE_EXPERIMENTAL;               }

    enum ParticleType{
        BARYON                        = 0,
        DARK_MATTER                   = 1,
        WIND                          = 2,
        STAR                          = 3,
        GAS                           = 4,
        AGN                           = 5,
        ALL                           = 6,
    };

    enum ParticleProperty{
        MASS              = 0,
        PHI               = 1,
        VELOCITY          = 2,
        SMOOTHING_LENGTH  = 3,
        VELOCITY_X        = 4,
        VELOCITY_Y        = 5,
        VELOCITY_Z        = 6,
        UU                = 7,
        TYPE              = 8,
        TEMPERATURE       = 9,
        SPHDensity        = 10,
        Velocities       = 11,

    };

    FloatProperty                    timeStep_;

protected:
    virtual void initialize();
    virtual void deinitialize();

    virtual void process();

    void changedVolumeDimensions();
    void changedParticleProperty();
    void changedSpreadMode();

private:

    enum SpreadMode{
        BACKWARD_LINEAR_INTERPOLATION = 0,
        WEIGHT_ADDITION               = 1,
        NEAREST_VOXEL                 = 2,
        SPH                           = 3,
        AMOUNT                        = 4,
        VELOCITYFIELD                 = 5,
    };

    VolumePort                       outport_;
    CMParticlePort                   inport_;

    FloatProperty                    volumeColor_;
    IntVec3Property                  volumeDimensions_;
    OptionProperty<SpreadMode>       spreadMode_;
    OptionProperty<ParticleProperty> particleProperty_;
    OptionProperty<ParticleType>     particleType_;
    StringProperty                   propertyUnit_;
    
    Volume *                         volume_;

    static const std::string loggerCat_;
};
}

#endif
