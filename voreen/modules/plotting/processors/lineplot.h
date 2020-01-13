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

#ifndef VRN_LINEPLOT_H
#define VRN_LINEPLOT_H

#include "plotprocessor.h"

#include "voreen/core/properties/floatproperty.h"

namespace voreen {

class VRN_CORE_API LinePlot : public PlotProcessor {

public:
    LinePlot();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "LinePlot";        }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE; }

protected:
    virtual void setDescriptions() {
        setDescription("This processor is able to plot lines and candlesticks.");
    }

private:
    // inherited methods
    virtual void render();
    virtual void renderData();
    virtual void renderAxes();
    virtual void setPlotStatus();
    virtual void readFromInport();
    virtual void calcDomains();
    virtual void toggleProperties();
    virtual void createPlotLabels();

    /// create line labels
    void createLineLabels();

    // properties
    FloatProperty lineWidth_;
    FloatProperty pointSize_;
    BoolProperty logXAxis_;
    BoolProperty logYAxis_;
    BoolProperty renderLineLabel_;

    static const std::string loggerCat_;
};

}   //namespace

#endif // VRN_LINEPLOT_H
