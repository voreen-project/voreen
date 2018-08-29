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

#ifndef VRN_VOLUMEEALWORLDMAPPING_H
#define VRN_VOLUMEEALWORLDMAPPING_H

#include <string>
#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/optionproperty.h"

namespace voreen {

class Volume;

class VRN_CORE_API VolumeRealWorldMapping : public VolumeProcessor {
public:
    VolumeRealWorldMapping();
    virtual Processor* create() const;

    virtual std::string getCategory() const   { return "Volume Processing"; }
    virtual std::string getClassName() const  { return "VolumeRealWorldMapping";     }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL;   }

protected:
    virtual void setDescriptions() {
        setDescription("Modifies the volume's real world mapping");
    }

    virtual void initialize();
    virtual void process();


private:
    enum Mode{
        SCALE_OFFSET,
        INTERVAL,
        NORMALIZED,
        DENORMALIZED
    };

    void reset();
    void updatePropertyVisibility();
    void updatePropertyValues();
    RealWorldMapping getRealWorldMapping();
    VolumePort inport_;
    VolumePort outport_;

    BoolProperty enableProcessing_;

    OptionProperty<Mode> mode_;
    FloatIntervalProperty realWorldRange_;
    FloatProperty offset_;
    FloatProperty scale_;


    static const std::string loggerCat_; ///< category used in logging
};

}   //namespace

#endif // VRN_VOLUMERVMNORMALIZATION_H
