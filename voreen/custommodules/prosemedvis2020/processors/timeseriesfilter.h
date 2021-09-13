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

#ifndef VRN_TIMESERIESFILTER_H
#define VRN_TIMESERIESFILTER_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/string/stringlistproperty.h"

#include <string>
#include <vector>

#include "../modules/plotting/ports/plotport.h"
#include "../modules/plotting/datastructures/plotdata.h"

namespace voreen {
    

class VRN_CORE_API TimeseriesFilter : public Processor {
public:

    TimeseriesFilter();
    virtual Processor* create() const;
    virtual std::string getClassName() const;
    virtual std::string getCategory() const;

protected:

    virtual void setDescriptions();
    virtual void process();

private:

    // Ports
    PlotPort inport_;        ///< inport used to pass the plot-data
    PlotPort outport_;       ///< outport used to pass the plot-data

    // Propertys
    StringProperty regionProp_;  // Filter-Prefix
    StringListProperty tracerProp_;  // Tracer

    std::vector<std::string> tracers_;

    /**
     *	Returns a vector containing the names of all tracers using the column names
     *	from the plotDataInport_
     **/
    std::vector<std::string> getTracerNames() const;

};

} // namespace

#endif // VRN_TIMESERIESFILTER_H
