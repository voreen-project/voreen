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

#ifndef VRN_FIELDPARALLELPLOTHISTOGRAM_H
#define VRN_FIELDPARALLELPLOTHISTOGRAM_H

#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/datastructures/volume/histogram.h"
#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "modules/plotting/ports/plotport.h"
#include "modules/plotting/datastructures/plotdata.h"
#include "modules/plotting/datastructures/plotcell.h"

namespace voreen {

class PlotLibrary;

/**
 *
 */
class VRN_CORE_API FieldParallelPlotHistogram : public RenderProcessor {

public:
    FieldParallelPlotHistogram();
    virtual ~FieldParallelPlotHistogram();
    virtual Processor* create() const;

    virtual std::string getClassName() const      { return "FieldParallelPlotHistogram"; }
    virtual std::string getCategory() const       { return "Plotting";                   }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;      }

protected:
    virtual void process();
    virtual void initialize();
    virtual void onEvent(tgt::Event* e);
    virtual bool isReady() const;

    static const std::string loggerCat_;

protected:

    void initHistogram();
    void updateSelectedValues();
    void mouseEvent(tgt::MouseEvent* e);

    VolumePort inportVolume_;
    RenderPort outport_;

    FloatIntervalProperty valueRange_;

    std::unique_ptr<VolumeHistogramIntensity> volumeHistogramIntensity_;
    std::unique_ptr<PlotData> data_;

    std::unique_ptr<PlotLibrary> plotLib_;
    std::pair<tgt::vec2, tgt::vec2> margins_;
    tgt::vec2 selectedRange_;
    bool isSelectionMode_;
};

} // namespace

#endif // VRN_FIELDPARALLELPLOTHISTOGRAM_H
