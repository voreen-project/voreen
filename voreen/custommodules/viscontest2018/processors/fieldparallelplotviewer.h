/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
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

#ifndef VRN_FIELDPARALLELPLOTVIEWER_H
#define VRN_FIELDPARALLELPLOTVIEWER_H

#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/io/volumeserializerpopulator.h"

#include "voreen/core/ports/volumeport.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"

#include "../datastructures/fieldplotdata.h"
#include "../ports/ensembledatasetport.h"
#include "../ports/fieldplotdataport.h"
#include "../properties/stringlistproperty.h"

namespace voreen {

class PlotLibrary;

/**
 *
 */
class VRN_CORE_API FieldParallelPlotViewer : public RenderProcessor {

public:
    FieldParallelPlotViewer();
    virtual ~FieldParallelPlotViewer();
    virtual Processor* create() const;

    virtual std::string getClassName() const      { return "FieldParallelPlotViewer"; }
    virtual std::string getCategory() const       { return "Plotting";                }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;   }

protected:

    virtual void onEvent(tgt::Event* e);

    virtual void initialize();
    virtual void deinitialize();

    virtual void beforeProcess();
    virtual void process();

    virtual bool isReady() const;

    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);
    virtual bool rebuildShader();

protected:

    void renderPlot();
    void renderAxes();

    void switchChannel();
    void loadPlotData();
    void updatePlotData();
    bool checkConsistency();

    void applyThreshold(VolumeRAM_Float* volume);

    void mouseEvent(tgt::MouseEvent* e);
    void updateSelection();

    virtual void setDescriptions() {
        setDescription("");
    }

    FieldPlotDataPort plotDataInport_;
    EnsembleDatasetPort ensembleInport_;
    RenderPort outport_;
    VolumePort volumeOutport_;
    RenderPort privatePort_;

    TransFunc1DKeysProperty transferFunc_;
    OptionProperty<std::string> renderedChannel_;
    StringListProperty renderedRuns_;

    TransFunc1DKeysProperty volumeTransferFunc_;
    FloatIntervalProperty valueRange_;
    FloatIntervalProperty timeInterval_;
    StringListProperty selectedRuns_;
    BoolProperty hasLogarithmicDensity_;
    StringProperty xUnit_;
    StringProperty yUnit_;

    ShaderProperty plotShader_;

    VolumeBase* plotData_;
    VolumeBase* channelSlices_;
    std::unique_ptr<Volume> currentPlot_;
    std::unique_ptr<tgt::Texture> plotTexture_;

    float timeStepPosition_;
    float viewPortWidth_;

    bool requirePlotDataUpdate_;
    bool requireConsistencyCheck_;
    bool consistent_;

    // UI
    std::unique_ptr<PlotLibrary> plotLib_;
    tgt::vec2 selectionStart_;
    tgt::vec2 selectionEnd_;
    bool isSelectionMode_;

    static const std::string loggerCat_;
};

} // namespace

#endif
