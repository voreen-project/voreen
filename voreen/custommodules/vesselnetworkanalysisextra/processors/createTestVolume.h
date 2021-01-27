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

#ifndef VRN_createTestVolume_H
#define VRN_createTestVolume_H

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/optionproperty.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "tgt/vector.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include <random>

namespace voreen {

typedef std::mt19937 random_engine;

class createTestVolume : public VolumeProcessor {
public:
    createTestVolume();
    virtual ~createTestVolume();

    virtual std::string getClassName() const         { return "createTestVolume";      }
    virtual std::string getCategory() const       { return "Volume Processing"; }
    virtual std::string setDescriptions() const       { return "Volume Processing"; }
    virtual VoreenSerializableObject* create() const;
    virtual void setDescriptions() { setDescription( "Creates test Volume with a given distribution"); }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;   }
    virtual void process();
    virtual void setRandomPoints(VolumeRAM*, voreen::random_engine&);
    virtual void setMean(random_engine &randomEngine,std::uniform_int_distribution<uint32_t> dist,std::uniform_int_distribution<uint32_t> distYZ,std::vector<int>& meanX,std::vector<int>& meanY,std::vector<int>& meanZ , int iteration);
    virtual void densityFunction(int , int , int , int, int, int, int , float&);
    virtual void qMeanDistance(int&, int, int&);



private:

    VolumePort outport_;

    IntProperty volX_;
    IntProperty volY_;
    IntProperty volZ_;

    IntProperty stdevx_;

    IntProperty amount_;
    IntProperty predeterminedSeed_;

    std::random_device randomDevice;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_VOLUMEFLOODFILL_H
