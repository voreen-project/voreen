/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "fieldparallelplothistogram.h"

#include "tgt/immediatemode/immediatemode.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/volume/histogram.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/utils/glsl.h"

#include "modules/plotting/utils/plotlibrary/plotlibrary.h"
#include "modules/plotting/utils/plotlibrary/plotlibraryopengl.h"

#include "modules/ensembleanalysis/utils/utils.h"

namespace voreen {

static const tgt::ivec2 MARGINS(50);
static const tgt::vec2 NO_SELECTION(-1.0f, 1.0f);

const std::string FieldParallelPlotHistogram::loggerCat_("voreen.ensembleanalysis.FieldParallelPlotHistogram");

FieldParallelPlotHistogram::FieldParallelPlotHistogram()
        : RenderProcessor()
        , inportVolume_(Port::INPORT, "fpp.histogram.data.volume", "Field Plot Volume Input")
        , outport_(Port::OUTPORT, "outport", "Histogram Image Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
        , valueRange_("valueRange", "Value Range", tgt::vec2(0.0f, 0.0f), 0.0f, 0.0f)
        , selectedRange_(NO_SELECTION)
        , volumeHistogramIntensity_(nullptr)
        , data_(nullptr)
        , plotLib_(new PlotLibraryOpenGl())
        , isSelectionMode_(false)
{
    addPort(inportVolume_);
    addPort(outport_);

    addProperty(valueRange_);

    ON_CHANGE(inportVolume_, FieldParallelPlotHistogram, initHistogram);
    initHistogram();
}

FieldParallelPlotHistogram::~FieldParallelPlotHistogram() {
}

Processor* FieldParallelPlotHistogram::create() const {
    return new FieldParallelPlotHistogram();
}

void FieldParallelPlotHistogram::initialize() {
    RenderProcessor::initialize();

    plotLib_->setDimension(PlotLibrary::TWO);
    plotLib_->setBarWidth(1.0);
    plotLib_->setBarGroupingMode(PlotLibrary::GROUPED);
    plotLib_->setAxesWidth(1.0f);
    plotLib_->setDrawingColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
    plotLib_->setLineWidth(1.0f);
    plotLib_->setMaxGlyphSize(1.0f);
    plotLib_->setMarginBottom(MARGINS.y);
    plotLib_->setMarginTop(MARGINS.y);
    plotLib_->setMarginLeft(MARGINS.x);
    plotLib_->setMarginRight(MARGINS.x);
    plotLib_->setMinimumScaleStep(32, PlotLibrary::X_AXIS);
    plotLib_->setMinimumScaleStep(32, PlotLibrary::Y_AXIS);
    plotLib_->setLogarithmicAxis(true, PlotLibrary::Y_AXIS);
}

void FieldParallelPlotHistogram::initHistogram() {
    if(!inportVolume_.hasData()) return;

    const VolumeRAM_Float* volume = dynamic_cast<const VolumeRAM_Float*>(inportVolume_.getData()->getRepresentation<VolumeRAM>());

    valueRange_.setMinValue(volume->min());
    valueRange_.setMaxValue(volume->max());
    valueRange_.set(tgt::vec2(volume->min(), volume->max()));

    volumeHistogramIntensity_.reset(dynamic_cast<VolumeHistogramIntensity*>(VolumeHistogramIntensity().createFrom(inportVolume_.getData())));

    data_.reset(new PlotData(1, 1));
    data_->setColumnLabel(0,"Index");
    data_->setColumnLabel(1,"Frequency");
    for (int i = 0; i<volumeHistogramIntensity_->getBucketCount(); i++) {
        float value = volumeHistogramIntensity_->getLogNormalized(i);
        std::vector<PlotCellValue> cells;
        cells.push_back(PlotCellValue(i));
        cells.push_back(PlotCellValue(value));
        data_->insert(cells);
    }

    invalidate();
}

bool FieldParallelPlotHistogram::isReady() const {
    if(!inportVolume_.isReady()) {
        setNotReadyErrorMessage("Inport volume is not set.");
        return false;
    }
    if(!outport_.isReady()) {
        setNotReadyErrorMessage("No outport connected.");
        return false;
    }
    return true;
}

void FieldParallelPlotHistogram::process() {
    outport_.activateTarget();
    outport_.clearTarget();

    // Recalculate margins.
    margins_.first  = mapRange(tgt::vec2(MARGINS), tgt::vec2::zero, tgt::vec2(outport_.getSize()), -tgt::vec2::one, tgt::vec2::one);
    margins_.second = mapRange(tgt::vec2(outport_.getSize()-MARGINS), tgt::vec2::zero, tgt::vec2(outport_.getSize()), -tgt::vec2::one, tgt::vec2::one);

    // Set Plot status.
    plotLib_->setWindowSize(outport_.getSize());
    plotLib_->setDomain(Interval<plot_t>(valueRange_.getMinValue(), valueRange_.getMaxValue()), PlotLibrary::X_AXIS);
    plotLib_->setDomain(Interval<plot_t>(volumeHistogramIntensity_->getHistogram().getMinValue(), volumeHistogramIntensity_->getHistogram().getMaxValue()), PlotLibrary::Y_AXIS);

    if (plotLib_->setRenderStatus()) {
        plotLib_->setDrawingColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
        plotLib_->renderAxes();
        plotLib_->setDrawingColor(tgt::Color(0, 0, 0, .5f));
        plotLib_->setFontSize(10);
        plotLib_->setFontColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
        plotLib_->renderAxisScales(PlotLibrary::X_AXIS, false);
        plotLib_->renderAxisScales(PlotLibrary::Y_AXIS, false);
        plotLib_->setFontSize(12);
        plotLib_->renderAxisLabel(PlotLibrary::X_AXIS, "Index");
        plotLib_->renderAxisLabel(PlotLibrary::Y_AXIS, "log. Frequency");
    }
    plotLib_->resetRenderStatus();

    // draw buckets
    IMode.begin(IMode.QUADS);
    IMode.color(0.0f, 0.6f, 1.0f, 1.0f);
    float width = (margins_.second.x - margins_.first.x) / volumeHistogramIntensity_->getBucketCount();
    for (int i = 0; i<volumeHistogramIntensity_->getBucketCount(); i++){
        float x = margins_.first.x + i * width;
        float y = margins_.first.y;

        float value = volumeHistogramIntensity_->getLogNormalized(i);
        float height = (margins_.second.y - margins_.first.y) * value;

        IMode.vertex(x, y);
        IMode.vertex(x+width, y);
        IMode.vertex(x+width, y+height);
        IMode.vertex(x, y+height);
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
        IMode.vertex(selectedRange_.x, margins_.second.y);
        IMode.vertex(selectedRange_.y, margins_.second.y);
        IMode.vertex(selectedRange_.y, margins_.first.y);
        IMode.vertex(selectedRange_.x, margins_.first.y);
        IMode.end();

        glLineWidth(1.0f);
        glDisable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ZERO);
        glDepthFunc(GL_LESS);
    }

    IMode.color(tgt::vec4::one);
    outport_.deactivateTarget();
}

void FieldParallelPlotHistogram::updateSelectedValues() {

    if(selectedRange_ == NO_SELECTION) {
        valueRange_.set(tgt::vec2(valueRange_.getMinValue(), valueRange_.getMaxValue()));
    }
    else {
        tgt::vec2 value = mapRange(selectedRange_,
                                   tgt::vec2(margins_.first.x),
                                   tgt::vec2(margins_.second.x),
                                   tgt::vec2(valueRange_.getMinValue()),
                                   tgt::vec2(valueRange_.getMaxValue())
        );

        float lower = std::min(value.x, value.y);
        float upper = std::max(value.x, value.y);

        valueRange_.set(tgt::vec2(lower, upper));
    }
}

void FieldParallelPlotHistogram::mouseEvent(tgt::MouseEvent* e) {

    int ex = tgt::clamp(e->x(), MARGINS.x, e->viewport().x-MARGINS.x);
    int ey = tgt::clamp(e->y(), MARGINS.y, e->viewport().y-MARGINS.y);

    float mx = mapRange(ex, 0, e->viewport().x, -1.0f,  1.0f);

    switch (e->action()) {
    case tgt::MouseEvent::PRESSED:
        // TODO: test
        if(ex != e->x() || ey != e->y())
            break;

        if (selectedRange_ == NO_SELECTION) {
            isSelectionMode_ = true;
            selectedRange_.x = mx;
            selectedRange_.y = mx;
        }
        else {
            selectedRange_ = NO_SELECTION;
            updateSelectedValues();
        }
        break;
    case tgt::MouseEvent::MOTION:
        if (isSelectionMode_) {
            selectedRange_.y = mx;
        }
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
