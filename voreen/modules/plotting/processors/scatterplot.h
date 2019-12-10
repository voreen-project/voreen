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

#ifndef VRN_SCATTERPLOT_H
#define VRN_SCATTERPLOT_H

#include "plotprocessor.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"

namespace voreen {

/**
 * This processor is able to create 2D and 3D scatter plots.
 */
class VRN_CORE_API ScatterPlot : public PlotProcessor {
public:
    ScatterPlot();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "ScatterPlot";     }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE; }

private:
    virtual void setDescriptions() {
        setDescription("This processor is able to create 2D and 3D Scatterplots.");
    }

    // inherited methods
    virtual void render();
    virtual void renderData();
    virtual void renderAxes();
    virtual void setPlotStatus();
    virtual void readFromInport();
    virtual void calcDomains();
    virtual void toggleProperties();

    /// Loads (and creates) needed textures.
    void loadTextures(std::string path);

    BoolProperty lighting_;

    static const std::string loggerCat_;
};

}   //namespace

#endif // VRN_SCATTERPLOT_H
