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

#ifndef VRN_FLOWINDICATORSELECTION_H
#define VRN_FLOWINDICATORSELECTION_H

#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/matrixproperty.h"

#include "../../datastructures/flowparameters.h"
#include "../../ports/flowparametrizationport.h"

namespace voreen {

/**
 * This processor is being used to select in and out flow.
 */
class VRN_CORE_API FlowIndicatorSelection : public RenderProcessor {
public:
    FlowIndicatorSelection();
    virtual Processor* create() const         { return new FlowIndicatorSelection();    }

    virtual std::string getClassName() const  { return "FlowIndicatorSelection";        }
    virtual std::string getCategory() const   { return "Simulation";                    }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL;         }

    virtual bool isReady() const;
    virtual void process();

protected:
    virtual void setDescriptions() {
        setDescription("This processor is being used to select in and out flow.");
    }

private:

    void selectRegion(tgt::MouseEvent* e);

    RenderPort renderInport_;
    RenderPort renderOutport_;

    FlowParametrizationPort flowParametrizationPort_;

    StringProperty ensembleName_;
    FloatProperty simulationTime_;
    FloatProperty temporalResolution_;

    FloatMat4Property pickingMatrix_;       ///< Picking matrix from SliceViewer

    tgt::plane plane_;
    tgt::vec3 xVec_;
    tgt::vec3 yVec_;
    tgt::vec3 origin_;
    float samplingRate_;
    tgt::ivec2 resolution_;

    std::vector<FlowIndicator> flowIndicators_;

    bool rebuildOutput_;

    static const std::string loggerCat_;
};

}   //namespace

#endif
