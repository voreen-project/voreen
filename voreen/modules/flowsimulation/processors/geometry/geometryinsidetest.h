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

#ifndef VRN_GEOMETRYINSIDETEST_H
#define VRN_GEOMETRYINSIDETEST_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/temppathproperty.h"

namespace voreen {

class VRN_CORE_API GeometryInsideTest : public Processor {
public:
    GeometryInsideTest();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "GeometryInsideTest"; }
    virtual std::string getCategory() const  { return "Geometry";               }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL;  }

protected:

    virtual void setDescriptions() {
        setDescription("Creates a binary volume for the incoming geometry."
                       "Value 1 determines inside, value 0 means outside."
                       "Ensure that the geometry has been converted to mm!");
        dimensions_.setDescription("Dimension of the longest side of the output volume.");
    }

    virtual bool isReady() const;
    virtual void process();

private:

    GeometryPort inport_;
    VolumePort outport_;

    IntProperty dimensions_;

    static const std::string loggerCat_;
};

}   //namespace

#endif
