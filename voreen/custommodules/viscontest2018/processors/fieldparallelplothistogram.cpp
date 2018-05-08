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

#include "fieldparallelplothistogram.h"

#include "tgt/immediatemode/immediatemode.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/volume/histogram.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/utils/glsl.h"

namespace voreen {

static const tgt::vec2 NO_SELECTION(-1.0f, 1.0f);

const std::string FieldParallelPlotHistogram::loggerCat_("voreen.viscontest2018.FieldParallelPlotHistogram");

FieldParallelPlotHistogram::FieldParallelPlotHistogram()
        : RenderProcessor()
        , inportVolume_(Port::INPORT, "fpp.histogram.data.volume", "Field Plot Volume Input")
        , outport_(Port::OUTPORT, "outport", "Histogram Image Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
        , outportPlot_(Port::OUTPORT, "plotOutport", "Plot Outport")
        , valueRange_("valueRange", "Value Range", tgt::vec2(0.0f, 0.0f), 0.0f, 0.0f)
        , selectedRange_(NO_SELECTION)
        , volumeHistogramIntensity_(nullptr)
        , isSelectionMode_(false)
{
    addPort(inportVolume_);
    addPort(outport_);
    addPort(outportPlot_);

    addProperty(valueRange_);

    ON_CHANGE(inportVolume_, FieldParallelPlotHistogram, initHistogram);
    initHistogram();
}

FieldParallelPlotHistogram::~FieldParallelPlotHistogram() {
}

Processor* FieldParallelPlotHistogram::create() const {
    return new FieldParallelPlotHistogram();
}

void FieldParallelPlotHistogram::initHistogram() {
    if(!inportVolume_.hasData()) return;

    const VolumeRAM_Float* volume = dynamic_cast<const VolumeRAM_Float*>(inportVolume_.getData()->getRepresentation<VolumeRAM>());

    valueRange_.setMinValue(volume->min());
    valueRange_.setMaxValue(volume->max());
    valueRange_.set(tgt::vec2(volume->min(), volume->max()));

    volumeHistogramIntensity_.reset(static_cast<VolumeHistogramIntensity*>(VolumeHistogramIntensity().createFrom(inportVolume_.getData())));


    int bucketCount = volumeHistogramIntensity_->getBucketCount();

    PlotData* plotPortData = new PlotData(1,1);
    plotPortData->setColumnLabel(0,"Index");
    plotPortData->setColumnLabel(1,"Frequency");
    for (int i = 0; i<bucketCount; i++) {
        float value = volumeHistogramIntensity_->getLogNormalized(i);
        std::vector<PlotCellValue> cells(0);
        cells.push_back((PlotCellValue(i)));
        cells.push_back((PlotCellValue(value)));
        plotPortData->insert(cells);
    }
    outportPlot_.setData(plotPortData, true);

    invalidate();
}

bool FieldParallelPlotHistogram::isReady() const {
    if(!inportVolume_.isReady()) {
        setNotReadyErrorMessage("Inport volume is not set.");
        return false;
    }
    if(!outport_.isReady() && !outportPlot_.isReady()) {
        setNotReadyErrorMessage("No outport connected.");
        return false;
    }
    return true;
}

void FieldParallelPlotHistogram::process() {
    outport_.activateTarget();
    outport_.clearTarget();

    drawHistogram();

    outport_.deactivateTarget();
}

void FieldParallelPlotHistogram::drawHistogram() {
    if(!inportVolume_.hasData()) return;

    float y_offset = -1.0f;
    float globalWidth = 2.0f;
    float globalHeight = 2.0f;

    float height_factor = globalHeight;
    float width = globalWidth / static_cast<float>(volumeHistogramIntensity_->getBucketCount());

    // draw buckets
    IMode.begin(IMode.QUADS);
    IMode.color(0.0f, 0.6f, 1.0f, 1.0f);
    for (int i = 0; i<volumeHistogramIntensity_->getBucketCount(); i++){

        float value = volumeHistogramIntensity_->getLogNormalized(i);
        float x = (i * width) - 1;
        float y_max = y_offset + value * height_factor; 

        IMode.vertex(x, y_offset);
        IMode.vertex(x+width, y_offset);
        IMode.vertex(x+width, y_max);
        IMode.vertex(x, y_max);
    }
    IMode.end();

    // draw selection
    if(selectedRange_ != NO_SELECTION && selectedRange_.x != selectedRange_.y) {
        glDepthFunc(GL_ALWAYS);
        IMode.color(1.0f, 1.0f, 1.0f, 0.5f);
        glLineWidth(1.0f);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_BLEND);

        IMode.begin(tgt::ImmediateMode::QUADS);
            float xPos = selectedRange_.x;
            float xPosMax = selectedRange_.y;
            float yPos = -1.0f;
            IMode.vertex(xPos, 1.0f);
            IMode.vertex(xPosMax, 1.0f);
            IMode.vertex(xPosMax, -1.0f);
            IMode.vertex(xPos, -1.0f);
        IMode.end();

        glLineWidth(1.0f);
        glDisable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ZERO);
        glDepthFunc(GL_LESS);
    }

    IMode.color(tgt::vec4::one);
}

void FieldParallelPlotHistogram::updateSelectedValues() {

    const tgt::vec2& valueRange = tgt::vec2(valueRange_.getMinValue(), valueRange_.getMaxValue());

    float xValue = ((selectedRange_.x + 1.0f) / 2.0f) * (valueRange.y - valueRange.x) + valueRange.x;
    float yValue = ((selectedRange_.y + 1.0f) / 2.0f) * (valueRange.y - valueRange.x) + valueRange.x;

    float lower = std::min(xValue, yValue);
    float upper = std::max(xValue, yValue);

    valueRange_.set(tgt::vec2(lower, upper));
}

void FieldParallelPlotHistogram::mouseEvent(tgt::MouseEvent* e) {

    viewPortWidth_ = static_cast<float>(e->viewport().x);

    int x = tgt::clamp(e->x(), 0, e->viewport().x);
    float mousePosition = (x / viewPortWidth_) * 2.0f - 1.0f;

    switch (e->action()) {
    case tgt::MouseEvent::PRESSED:
        if (selectedRange_ == NO_SELECTION) {
            isSelectionMode_ = true;
            selectedRange_.x = mousePosition;
            selectedRange_.y = mousePosition;
        }
        else {
            selectedRange_ = NO_SELECTION;
        }
        break;
    case tgt::MouseEvent::MOTION:
        if (isSelectionMode_) {
            selectedRange_.y = mousePosition;
        }
        break;
    case tgt::MouseEvent::EXIT:
        selectedRange_.y = tgt::clamp(selectedRange_.y, -1.0f, 1.0f);
        break;
    case tgt::MouseEvent::RELEASED:
        if (selectedRange_.x == selectedRange_.y) {
            selectedRange_ = NO_SELECTION;
        }
        isSelectionMode_ = false;
        updateSelectedValues();
        break;
    default:
        e->ignore();
        return;
    }

    e->accept();

    // Trigger redraw.
    invalidate();
}

void FieldParallelPlotHistogram::onEvent(tgt::Event* e) {
    tgt::MouseEvent* event = dynamic_cast<tgt::MouseEvent*>(e);

    if (!event)
        RenderProcessor::onEvent(e);
    else
        mouseEvent(event);
}


} // namespace voreen
