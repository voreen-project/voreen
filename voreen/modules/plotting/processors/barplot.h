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

#ifndef VRN_BARPlOT_H
#define VRN_BARPlOT_H

#include "plotprocessor.h"

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/vectorproperty.h"

namespace voreen {

class VRN_CORE_API BarPlot : public PlotProcessor {

public:
    BarPlot();

    virtual Processor* create() const;

    virtual std::string getClassName() const { return "BarPlot"; }
    virtual CodeState getCodeState() const   { return Processor::CODE_STATE_STABLE; }

protected:

    virtual void setDescriptions() {
        setDescription("This processor is able to plot bars. These bars can be grouped, merged or stacked.");
    }

    // inherited methods
    virtual void render();
    virtual void renderData();
    virtual void renderAxes();
    virtual void setPlotStatus();
    virtual void readFromInport();
    virtual void calcDomains();
    virtual void toggleProperties();
    virtual void createPlotLabels();
    virtual void select(tgt::MouseEvent* e, bool highlight, bool label, bool zoom, bool additive);

    IntProperty barWidth_;                                 ///< bar width in percent of maximal width
    FloatVec2Property shear_;                              ///< shearfactor for fake 3d effect
    FloatProperty squeezeFactor_;                          ///< squeeze factor for merged bars
    OptionProperty<PlotLibrary::BarGroupingMode> barMode_; ///< select bar mode (grouped, stacked, merged)

    static const std::string loggerCat_;
};

}   //namespace

#endif // VRN_BARPLOT_H
