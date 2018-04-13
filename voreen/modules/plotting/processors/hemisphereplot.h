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

#ifndef VRN_HEMISPHEREPLOT_H
#define VRN_HEMISPHEREPLOT_H

#include "surfaceplot.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"

#include "../ext/triangle/include/del_interface.hpp"

#include <vector>
#include <list>

namespace voreen {

/**
 * This processor is able to create 3D surface plots.
 */
class VRN_CORE_API HemispherePlot : public SurfacePlot {
public:
    HemispherePlot();

    virtual Processor* create() const;
    virtual std::string getClassName() const { return "HemispherePlot"; }
    virtual CodeState getCodeState() const { return Processor::CODE_STATE_TESTING; }

protected:
    virtual void setDescriptions() {
        setDescription("Shows the Hemisphere of the Plotdata.");
    }

    virtual void initialize();
    virtual void deinitialize();

    virtual void readFromInport();
    virtual void generateDelaunay();
    virtual void calcDomains();
    virtual void setPlotStatus();

    //virtual void selectDataFromFunction();
    static const std::string loggerCat_;
};

}   //namespace

#endif // VRN_HEMISPHEREPLOT_H
