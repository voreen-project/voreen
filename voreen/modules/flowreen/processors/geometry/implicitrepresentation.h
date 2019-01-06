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

#ifndef VRN_IMPLICITREPRESENTATION_H
#define VRN_IMPLICITREPRESENTATION_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/temppathproperty.h"

namespace voreen {

class VRN_CORE_API ImplicitRepresentation : public Processor {
public:
    ImplicitRepresentation();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "ImplicitRepresentation"; }
    virtual std::string getCategory() const  { return "Geometry";               }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL;  }
    virtual bool usesExpensiveComputation() const { return false; }

protected:

    virtual void setDescriptions() {
        setDescription("Creates an implicit representation of the incoming geometry."
                       "Ensure that the geometry has been converted to mm!");
        method_.setDescription("Method to be used for inside/outside test."
                               "Fast requires the geometry to be a close mesh."
                               "Accucate, however, can take a while to calculate.");
        dimensions_.setDescription("Dimensions (cubed) of the output volume.");
        path_.setDescription("Path where the input geometry will be stored."
                             "It is also possible to open an existing file for the"
                             "calculation. In this case, the input geometry will be ignored.");
    }

    virtual bool isReady() const;
    virtual void process();

private:

    GeometryPort inport_;
    VolumePort outport_;

    IntOptionProperty method_;
    IntProperty dimensions_;
    TempPathProperty path_;

    static const std::string loggerCat_;
};

}   //namespace

#endif // VRN_IMPLICITREPRESENTATION_H
