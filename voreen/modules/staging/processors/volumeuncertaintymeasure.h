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

#ifndef VRN_VOLUMEUNCERTAINTYMEASURE_H
#define VRN_VOLUMEUNCERTAINTYMEASURE_H

#include <string>
#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/boolproperty.h"

namespace voreen {

class Volume;

/**
 * Transforms an incoming uncertainty volume by the use of a VolumeOperator.
 * Currently only handles floats but in theory could handle every datatype,
 * the VolumeOperator is capable of.
 */
class VRN_CORE_API VolumeUncertaintyMeasure : public CachingVolumeProcessor {
public:
    VolumeUncertaintyMeasure();
    virtual ~VolumeUncertaintyMeasure();
    virtual Processor* create() const;

    virtual std::string getCategory() const   { return "Volume Processing"; }
    virtual std::string getClassName() const  { return "VolumeUncertaintyMeasure"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_TESTING; }

    virtual bool usesExpensiveComputation() const { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Calculates the uncertainty within the incoming volume. Each voxel v gets transformed by the function f(v) = 1 - |2v - 1|");
        inport_.setDescription("Requires float volume whose values are within range [0, 1].");
        enableProcessing_.setDescription("Toggles calculation. If enabled, a new VolumeRAM will be generated. The input volume will be passed through otherwise.");
    }

    virtual void process();

private:

    void forceUpdate();
    void measureVolumeUncertainty();

    VolumePort inport_;
    VolumePort outport_;

    BoolProperty enableProcessing_;

    bool forceUpdate_;

    static const std::string loggerCat_; ///< category used in logging
};

}   //namespace

#endif // VRN_VOLUMEUNCERTAINTYMEASURE_H
