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

#ifndef VRN_ZEBRACOORDINATECONVERTER_H
#define VRN_ZEBRACOORDINATECONVERTER_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/geometryport.h"
#include "voreen/core/ports/volumeport.h"

#include "modules/plotting/ports/plotport.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/matrixproperty.h"

#include "voreen/core/datastructures/geometry/pointlistgeometry.h"

namespace voreen {

class VRN_CORE_API CoordinateConverter : public Processor {

public:
    CoordinateConverter();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "CoordinateConverter";     }
    virtual std::string getCategory() const  { return "Utility";         }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL;  }

protected:
    virtual void setDescriptions() {
        setDescription("To Do");
    }

    virtual void process();

    /// Inport for the plot data.
    PlotPort inport_;

    /// Inport for the volume
    VolumePort volumeInport_;


    /// Outports for lines and point
    GeometryPort spineOutport_;
    GeometryPort gutOutport_;
    GeometryPort yolkOnsetOutport_;

    // "output" of plane that connects spine and gut
    FloatVec3Property planeNormal_;        
    FloatProperty planePosition_;         

    /// Can be used to discard tracks with less than the minimum length
    //IntProperty minimumLength_;

    // select a range of time steps
    //IntIntervalProperty columnRange_;


    static const std::string loggerCat_;

    FloatMat4Property transformMatrix_;
};

} // namespace

#endif
