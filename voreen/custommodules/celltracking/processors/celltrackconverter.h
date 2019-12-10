/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_CELLTRACKCONVERTER_H
#define VRN_CELLTRACKCONVERTER_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/geometryport.h"
#include "modules/plotting/ports/plotport.h"

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"

#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

namespace voreen {

/**
 * Gets plot data where each row specifies a track, consisting of coordinates x-y-z for consecutive time steps.
 * The processor converts this data into a PointSegmentListVec3 for rendering.
 */
class VRN_CORE_API CellTrackConverter : public Processor {

public:
    CellTrackConverter();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "CellTrackConverter";     }
    virtual std::string getCategory() const  { return "Utility";         }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL;  }

protected:
    virtual void setDescriptions() {
        setDescription("Converts plot data containing cell tracks into a PointSegmentListVec3. Each row has to resemble one track and x,y,z coordinates are written consecutively for each time step. \
                        A tuple (-1,-1,-1) is interpreted as not being a valid occurence of the object and is thus discarded.");
    }

    virtual void process();

    // adjust minimum range max value, adjust range of time step max value
    virtual void adjustPropertiesToInput();

    /// Inport for the plot data.
    PlotPort inport_;

    /// Outport for the point segment list
    GeometryPort outport_;

    /// Can be used to discard tracks with less than the minimum length
    IntProperty minimumLength_;

    // select a range of time steps
    IntIntervalProperty columnRange_;


    static const std::string loggerCat_;

};

} // namespace

#endif
