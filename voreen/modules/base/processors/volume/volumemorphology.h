/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
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

#ifndef VRN_VOLUMEMORPHOLOGY_H
#define VRN_VOLUMEMORPHOLOGY_H

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"

namespace voreen {

/**
 * Provides the basic morphologic operators dilation and erosion.
 */
class VRN_CORE_API VolumeMorphology : public CachingVolumeProcessor {
public:
    VolumeMorphology();
    virtual ~VolumeMorphology();
    virtual Processor* create() const;

    virtual std::string getClassName() const      { return "VolumeMorphology";  }
    virtual std::string getCategory() const       { return "Volume Processing"; }
    virtual CodeState getCodeState() const        { return CODE_STATE_STABLE;   }
    virtual bool usesExpensiveComputation() const { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Provides the basic morphologic operators dilation and erosion for volumetric data.\n\n"
                "Available shapes of the operators are \"Cube\" or \"Sphere\"."
                "Be aware that choosing \"Sphere\" may result in a higher runtime as in that case the kernel is not separable.");
    }

    virtual void process();

private:

    enum Shape {
        CUBE,
        SPHERE,
    };

    void forceUpdate();
    void applyOperator();

    VolumePort inport_;
    VolumePort outport_;

    BoolProperty enableProcessing_;
    StringOptionProperty morphologicOperator_;
    OptionProperty<Shape> morphologicOperatorShape_;
    IntOptionProperty kernelSize_;

    bool forceUpdate_;

    static const std::string loggerCat_;
};

}   //namespace

#endif // VRN_VOLUMEMORPHOLOGY_H
