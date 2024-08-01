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

#ifndef VRN_LP_PlOT_H
#define VRN_LP_PLOT_H

#include "../../plotting/processors/plotprocessor.h"

#include "voreen/core/ports/textport.h"

#include "voreen/core/properties/floatproperty.h"


namespace voreen {

class LP_Plot : public PlotProcessor {

public:
    LP_Plot();
    ~LP_Plot();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "LP_Plot"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL; }

    virtual bool isReady() const;
    //virtual void onEvent(tgt::Event* e);
protected:
    virtual void setDescriptions() {
        setDescription("");
    }

     virtual void initialize();
     virtual void deinitialize();

private:
    //ports
    TextPort textInport_;

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

    GLuint dataList_;         ///< display list of the data
    GLuint pickingList_;      ///< display list of the picking data


};

}   //namespace

#endif // VRN_LINEPLOT_H
