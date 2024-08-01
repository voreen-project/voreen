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

#include "timeseriesplot.h"

#include "modules/ensembleanalysis/utils/utils.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/interaction/camerainteractionhandler.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"

#include "modules/plotting/datastructures/plotcell.h"
#include "modules/plotting/datastructures/plotdata.h"
#include "modules/plotting/utils/plotlibrary/plotlibrary.h"
#include "modules/plotting/utils/plotlibrary/plotlibraryopengl.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"

#include "tgt/immediatemode/immediatemode.h"
#include "tgt/texture.h"

#include <random>
#include <chrono>

#include <Eigen/Eigenvalues>

namespace voreen {

static const int MAX_NUM_DIMENSIONS = 3;
static const int GRAB_ENABLE_DISTANCE = 0;
static const tgt::ivec2 MARGINS(75, 50);
static const tgt::vec3 FIRST_TIME_STEP_COLOR(1.0f, 0.0f, 0.0f);
static const tgt::vec3 LAST_TIME_STEP_COLOR = tgt::vec3::one;
static const tgt::vec3 MIN_DURATION_COLOR(1.0f, 0.0f, 0.0f);
static const tgt::vec3 MAX_DURATION_COLOR(0.0f, 0.0f, 1.0f);
static const tgt::vec3 FADE_OUT_COLOR = tgt::vec3::one;
static const float SELECTED_LINE_WIDTH = 6.0f;
static const float UNSELECTED_LINE_WIDTH = 3.0f;
static const tgt::vec3 SMALLER_THRESHOLD_COLOR(1.0f, 0.0f, 0.0f);
static const tgt::vec3 BETWEEN_THRESHOLD_COLOR(1.0f, 1.0f, 0.0f);
static const tgt::vec3 LARGER_THRESHOLD_COLOR(0.0f, 1.0f, 0.0f);
static const tgt::vec3 NOTHING_COLOR(0.5f, 0.5f, 0.5f);
static const tgt::vec3 SLAB_COLOR(0.0f, 0.0f, 1.0f);
static const tgt::vec3 PLUME_COLOR(1.0f, 0.0f, 0.0f);
const float PI = 3.14159265358979f;

const std::string TimeseriesPlot::fontName_("Vera.ttf");
const std::string TimeseriesPlot::loggerCat_("voreen.scivis21.TimeseriesPlot");

void TimeseriesPlot::Embedding::serialize(Serializer& s) const {
    s.serialize("nVectors", nVectors_);
    s.serialize("eigenvalues", eigenvalues_);
    if(!names_.empty())
        s.serialize("names", names_);
}

void TimeseriesPlot::Embedding::deserialize(Deserializer& s) {
    s.deserialize("nVectors", nVectors_);
    s.deserialize("eigenvalues", eigenvalues_);
    s.optionalDeserialize("names", names_, decltype(names_)());
}

TimeseriesPlot::TimeseriesPlot()
    : RenderProcessor()
    , timeseriesInport_(Port::INPORT, "timeseriesInport", "Timeseries Input", false)
    //, similarityMatrixInport_(Port::INPORT, "similaritymatrixin", "Similarity Matrix Input", false)
    , outport_(Port::OUTPORT, "outport", "Outport", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , privatePort_(Port::OUTPORT, "image.tmp", "image.tmp", false)
    , pickingBuffer_(Port::OUTPORT, "picking", "Picking", false)
    , eigenValueOutport_(Port::OUTPORT, "eigenvalueOutport", "Eigenvalues", false)
    //, pointListOutport_(Port::OUTPORT, "pointListOutport", "Selected Points Geometry")
    , calculateButton_("calculate", "Create Embeddings")
    , autoCalculate_("autoCalculate", "Auto Calculate", true)
    , progressBar_("progressBar", "Progress")
    , numIterations_("numIterations", "Number of Iterations", 1000, 1, 10000)
    , numEigenvalues_("numEigenvalues", "Number of Eigenvalues", MAX_NUM_DIMENSIONS, MAX_NUM_DIMENSIONS, MAX_NUM_DIMENSIONS)
    , numDimensions_("numDimensions", "Number of Dimensions", MAX_NUM_DIMENSIONS, 1, MAX_NUM_DIMENSIONS)
    , principleComponent_("principalComponent", "Principle Component", 1, 1, MAX_NUM_DIMENSIONS)
    //, scaleToMagnitude_("scaleToMagnitude", "Scale to Magnitude of Principle Component")
    , sphereRadius_("sphereRadius", "Sphere Radius", 0.01, 0.0f, 0.1f)
    , fontSize_("fontSize", "Font Size", 10, 1, 30)
    , showTooltip_("showTooltip", "Show Tooltip", true)
    , renderTimeSelection_("renderTimeSelection", "Render Time Selection", false, Processor::INVALID_RESULT, Property::LOD_ADVANCED)
    , colorCoding_("colorCoding", "Color Coding")
    , angleRange_("angleRange", "Angle Range", 90.0, 0.0, 180.0)
    , tempAnomalyRange_("tempAnomalyRange", "Temp. Anomaly Range", 0.0, -1000, 1000)
    , transferFunc_("transferFunction", "Transfer Function", Processor::INVALID_RESULT)
    //, renderedField_("renderedChannel", "Field")
    , renderedMembers_("renderedMembers", "Rendered Members")
    , firstSelectedMember_("selectedMembers", "First selected Members")
    , firstSelectedTimeInterval_("selectedTimeSteps", "First selected Time Interval", tgt::vec2(0.0f, 0.0f), 0.0f, 0.0f)
    , secondSelectedMember_("referenceMember", "Second selected Member")
    , secondSelectedTimeInterval_("referenceTimeStep", "Second selected Time Interval", tgt::vec2(0.0f, 0.0f), 0.0f, 0.0f)
    , saveFileDialog_("saveFileDialog", "Export Embedding", "Select file...", VoreenApplication::app()->getUserDataPath(),
                      "Voreen MDS Embedding (*.vmds)", FileDialogProperty::SAVE_FILE, Processor::INVALID_PATH, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
    , saveButton_("saveButton", "Save")
    , loadFileDialog_("loadFile", "Import Embedding", "Select file...", VoreenApplication::app()->getUserDataPath(),
                      "Voreen MDS Embedding (*.vmds)", FileDialogProperty::OPEN_FILE, Processor::INVALID_PATH, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
    , loadButton_("loadButton", "Load")
    , volumeDimension_("volume", "Volume", 100, 1, 2048)
    , position_("position", "Position", tgt::ivec3::zero, tgt::ivec3::zero, tgt::ivec3::one*2048)
    , camera_("camera", "Camera", tgt::Camera(tgt::vec3(0.0f, 0.0f, 3.5f), tgt::vec3(0.0f, 0.0f, 0.0f), tgt::vec3(0.0f, 1.0f, 0.0f)))
    , cameraHandler_(nullptr)
    , plotLib_(new PlotLibraryOpenGl())
    , lastHit_(boost::none)
{
    // Ports
    //NOTE: don't use adjustPropertiesToInput() as callback, this seems to be triggered each frame for RenderProcessors.
    addPort(timeseriesInport_);
    ON_CHANGE(timeseriesInport_, TimeseriesPlot, adjustToTimeserieslist);
    //addPort(similarityMatrixInport_);
    //ON_CHANGE(similarityMatrixInport_, TimeseriesPlot, adjustToTimeserieslist);
    addPort(outport_);
    addPort(eigenValueOutport_);
    //addPort(pointListOutport_);
    addPrivateRenderPort(privatePort_);
    addPrivateRenderPort(pickingBuffer_);

    // Calculation
    addProperty(calculateButton_);
        calculateButton_.setGroupID("calculation");
        calculateButton_.setReadOnlyFlag(true);
        ON_CHANGE(calculateButton_, TimeseriesPlot, createEmbeddings);
    addProperty(autoCalculate_);
        autoCalculate_.setGroupID("calculation");
    addProperty(progressBar_);
        progressBar_.setGroupID("calculation");
    addProgressBar(&progressBar_);

    addProperty(numIterations_);
        numIterations_.setGroupID("calculation");
    addProperty(numEigenvalues_);
        numEigenvalues_.setGroupID("calculation");
    setPropertyGroupGuiName("calculation", "Calculation");

    // Rendering
    addProperty(numDimensions_);
        numDimensions_.setGroupID("rendering");
        ON_CHANGE_LAMBDA(numDimensions_, [this] {
            // We need at least as many eigenvalues as dimensions to be rendered.
            numEigenvalues_.setMinValue(numDimensions_.get());

            // Eigenvalue index only useful when rendering 1D.
            principleComponent_.setVisibleFlag(numDimensions_.get() == 1);

            // Scaling only useful when rendering 3D
            //scaleToMagnitude_.setVisibleFlag(numDimensions_.get() == 3);

            // Adjust camera correctly when using 3D visualization.
            if (numDimensions_.get() == 3) {
                camera_.adaptInteractionToScene(tgt::Bounds(-tgt::vec3::one, tgt::vec3::one), 0.1f, true);
            }

            // We currently only allow for manual time range selection in 1D mode.
            firstSelectedTimeInterval_.setVisibleFlag(numDimensions_.get() == 1);
        });
    addProperty(principleComponent_);
        principleComponent_.setVisibleFlag(numDimensions_.get() == 1);
        principleComponent_.setGroupID("rendering");
//    addProperty(scaleToMagnitude_);
//        scaleToMagnitude_.setVisibleFlag(numDimensions_.get() == 3);
//        scaleToMagnitude_.setGroupID("rendering");
    addProperty(sphereRadius_);
        sphereRadius_.setGroupID("rendering");
    addProperty(fontSize_);
        fontSize_.setGroupID("rendering");
    addProperty(showTooltip_);
        showTooltip_.setGroupID("rendering");
    addProperty(renderTimeSelection_);
        renderTimeSelection_.setGroupID("rendering");
    addProperty(colorCoding_);
        //colorCoding_.addOption("member", "Only Member", COLOR_MEMBER);
        //colorCoding_.addOption("timeStep", "Only Time Step", COLOR_TIMESTEP);
        //colorCoding_.addOption("memberAndTimeStep", "Member and Time Step", COLOR_MEMBER_AND_TIMESTEP);
        //colorCoding_.addOption("duration", "Time Step Duration", COLOR_DURATION);
        //colorCoding_.selectByValue(COLOR_MEMBER_AND_TIMESTEP);
        colorCoding_.setGroupID("rendering");
       ON_CHANGE_LAMBDA(colorCoding_, [this] {
            if(colorCoding_.getValue() == COLOR_FIELD) {
                int index = (colorCoding_.getSelectedIndex()-3)/2;
                transferFunc_.get()->setDomain(getRange(index));
            }
        });
//    addProperty(renderedField_);
//        ON_CHANGE(renderedField_, TimeseriesPlot, outputEigenValues);
//        renderedField_.setGroupID("rendering");
    addProperty(angleRange_);
        angleRange_.setGroupID("rendering");
    addProperty(tempAnomalyRange_);
        tempAnomalyRange_.setGroupID("rendering");
    addProperty(transferFunc_);
        transferFunc_.setGroupID("rendering");
    addProperty(renderedMembers_);
        ON_CHANGE(renderedMembers_, TimeseriesPlot, renderedMembersChanged);
        renderedMembers_.setGroupID("rendering");
    setPropertyGroupGuiName("rendering", "Rendering");

    // Selection (Linking)
    addProperty(firstSelectedMember_);
        firstSelectedMember_.setGroupID("selection");
    addProperty(firstSelectedTimeInterval_);
        firstSelectedTimeInterval_.setGroupID("selection");
    addProperty(secondSelectedMember_);
        secondSelectedMember_.setGroupID("selection");
    addProperty(secondSelectedTimeInterval_);
        secondSelectedTimeInterval_.setGroupID("selection");
    setPropertyGroupGuiName("selection", "Selection");
    setPropertyGroupVisible("selection", false);

    // IO
    addProperty(saveFileDialog_);
        ON_CHANGE(saveFileDialog_, TimeseriesPlot, saveEmbeddings);
        saveFileDialog_.setGroupID("io");
    addProperty(saveButton_);
        ON_CHANGE(saveButton_, TimeseriesPlot, saveEmbeddings);
        saveButton_.setGroupID("io");
    addProperty(loadFileDialog_);
        ON_CHANGE(loadFileDialog_, TimeseriesPlot, loadEmbeddings);
        loadFileDialog_.setGroupID("io");
    addProperty(loadButton_);
        ON_CHANGE(loadButton_, TimeseriesPlot, loadEmbeddings);
        loadButton_.setGroupID("io");
    setPropertyGroupGuiName("io", "Save/Load MDS Plot");

    addProperty(volumeDimension_);
    ON_CHANGE_LAMBDA(volumeDimension_, [this] {
        position_.setMaxValue(tgt::ivec3::one*volumeDimension_.get());
    });
    volumeDimension_.setVisibleFlag(false);
    addProperty(position_);
    position_.setVisibleFlag(false);

    // Camera
    addProperty(camera_);
    cameraHandler_ = new CameraInteractionHandler("cameraHandler", "Camera", &camera_);
    cameraHandler_->setEnabled(true);
    addInteractionHandler(cameraHandler_);
}

TimeseriesPlot::~TimeseriesPlot() {
    delete cameraHandler_;
}

Processor* TimeseriesPlot::create() const {
    return new TimeseriesPlot();
}

void TimeseriesPlot::initialize() {
    RenderProcessor::initialize();
    sphere_.setSphereGeometry(1.0f, tgt::vec3::zero, tgt::vec4::one, 16);
}

void TimeseriesPlot::deinitialize() {
    sphere_.clear();
    RenderProcessor::deinitialize();
}

void TimeseriesPlot::process() {

    // Resize frame buffers accordingly.
    if(pickingBuffer_.getSize() != outport_.getSize())
        pickingBuffer_.resize(outport_.getSize());

    if(privatePort_.getSize() != outport_.getSize())
        privatePort_.resize(outport_.getSize());

    // Change depth func to apply rendering order properly.
    glDepthFunc(GL_LEQUAL);

    // Render picking pass.
    pickingBuffer_.activateTarget();
    pickingBuffer_.clearTarget();
    glLineWidth(7.0f);
    //renderingPass(true);
    pickingBuffer_.deactivateTarget();

    // Draw the members in new parameter space.
    bool threeDimensional = numDimensions_.get() == 3;
    if(threeDimensional) {
        outport_.activateTarget();
        outport_.clearTarget();

        renderingPass(false);

        // Draw tooltip.
        if(showTooltip_.get()) {
            renderTooltip();
        }

        outport_.deactivateTarget();
    }
    else {
        privatePort_.activateTarget();
        privatePort_.clearTarget();

        renderingPass(false);

        // Draw tooltip.
        if(showTooltip_.get()) {
            renderTooltip();
        }

        privatePort_.deactivateTarget();

        // Draw axes.
        renderAxes();
    }

    // Restore state.
    glDepthFunc(GL_LESS);
    glLineWidth(1.0f);
    IMode.color(tgt::col4::one);
}

void TimeseriesPlot::renderAxes() {

    outport_.activateTarget();
    outport_.clearTarget();

    // Set Plot status.
    plotLib_->setWindowSize(outport_.getSize());
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
    plotLib_->setMinimumScaleStep(32, PlotLibrary::Z_AXIS);
    if (numDimensions_.get() == 1) {
        plotLib_->setDomain(Interval<plot_t>(timeseriesInport_.getData()->getOverallStartTime(),
                                             timeseriesInport_.getData()->getOverallEndTime()), PlotLibrary::X_AXIS);
    }
    else {
        plotLib_->setDomain(Interval<plot_t>(-1.0, 1.0), PlotLibrary::X_AXIS);
    }
    plotLib_->setDomain(Interval<plot_t>(-1.0, 1.0), PlotLibrary::Y_AXIS);
    if (numDimensions_.get() == 3) {
        plotLib_->setDimension(PlotLibrary::THREE);
        plotLib_->setDomain(Interval<plot_t>(-1.0, 1.0), PlotLibrary::Z_AXIS);
    } else
        plotLib_->setDimension(PlotLibrary::TWO);

    if (plotLib_->setRenderStatus()) {
        plotLib_->setDrawingColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
        plotLib_->renderAxes();
        plotLib_->setDrawingColor(tgt::Color(0, 0, 0, .5f));
        plotLib_->setFontSize(fontSize_.get() + 2);

        switch (numDimensions_.get()) {
            case 1:
                plotLib_->renderAxisLabel(PlotLibrary::X_AXIS, "t [s]");
                if (principleComponent_.get() == 1)
                    plotLib_->renderAxisLabel(PlotLibrary::Y_AXIS, "1st PC");
                else if (principleComponent_.get() == 2)
                    plotLib_->renderAxisLabel(PlotLibrary::Y_AXIS, "2nd PC");
                else
                    plotLib_->renderAxisLabel(PlotLibrary::Y_AXIS, "3rd PC");

                plotLib_->setDrawingColor(tgt::Color(0, 0, 0, .5f));
                plotLib_->setFontSize(fontSize_.get());
                plotLib_->setFontColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
                plotLib_->renderAxisScales(PlotLibrary::X_AXIS, false);
                plotLib_->setFontSize(fontSize_.get() + 2);
                break;
            case 3:
                plotLib_->renderAxisLabel(PlotLibrary::Z_AXIS, "3rd PC");
                // Fallthrough
            case 2:
                plotLib_->renderAxisLabel(PlotLibrary::X_AXIS, "1st PC");
                plotLib_->renderAxisLabel(PlotLibrary::Y_AXIS, "2nd PC");
                break;
        }
    }
    plotLib_->resetRenderStatus();

    // Plot data
    float xMinMarginNDC = mapRange(MARGINS.x, 0, outport_.getSize().x, -1.0f, 1.0f);
    float xMaxMarginNDC = mapRange(outport_.getSize().x-MARGINS.x, 0, outport_.getSize().x, -1.0f, 1.0f);
    float yMinMarginNDC = mapRange(MARGINS.y, 0, outport_.getSize().y, -1.0f, 1.0f);
    float yMaxMarginNDC = mapRange(outport_.getSize().y-MARGINS.y, 0, outport_.getSize().y, -1.0f, 1.0f);

    tgt::TextureUnit::setZeroUnit();
    tgt::Texture* texture = privatePort_.getColorTexture();
    if(texture) {
        IMode.setTextureMode(tgt::ImmediateMode::TEX2D);
        texture->bind();

        IMode.begin(tgt::ImmediateMode::QUADS);
            IMode.texcoord(0.0f, 0.0f); IMode.vertex(xMinMarginNDC, yMinMarginNDC);
            IMode.texcoord(1.0f, 0.0f); IMode.vertex(xMaxMarginNDC, yMinMarginNDC);
            IMode.texcoord(1.0f, 1.0f); IMode.vertex(xMaxMarginNDC, yMaxMarginNDC);
            IMode.texcoord(0.0f, 1.0f); IMode.vertex(xMinMarginNDC, yMaxMarginNDC);
        IMode.end();

        IMode.setTextureMode(tgt::ImmediateMode::TEXNONE);
    }

    // Reset state.
    IMode.color(1.0f, 1.0f, 1.0f, 1.0f);
    outport_.deactivateTarget();
    LGL_ERROR;
}

void TimeseriesPlot::renderTooltip() const {

    const TimeSeriesList* timeserieslist = timeseriesInport_.getData();
    if(!timeserieslist || !lastHit_) {
        return;
    }

    glDisable(GL_DEPTH_TEST);

    tgt::Font font(VoreenApplication::app()->getFontPath(fontName_));
    font.setFontSize(fontSize_.get());

    const int legendOffsetX = 10.f;
    const int legendOffsetY = 10.f;

    const tgt::ivec2& screensize = outport_.getSize();

    const TimeSeries& series = timeserieslist->getTimeSeries(lastHit_->memberIdx);
    const TimeSeriesStep& timeStep = series.getTimeSeriesStep(lastHit_->timeStepIdx);

    std::string tooltip = "Member: ";
    tooltip += std::to_string(series.getPosition().x) + " " + std::to_string(series.getPosition().y) + " " + std::to_string(series.getPosition().z);//.getName();
    tooltip += "\n";
    tooltip += "Time:";
    tooltip += std::to_string(timeStep.time_);

    tgt::vec3 pos(legendOffsetX + lastHit_->x, screensize.y - lastHit_->y - 1, 0);
    tgt::vec2 size = font.getSize(pos, tooltip, screensize);
    pos.y -= size.y + legendOffsetY;

    tgt::vec2 ll = mapRange(pos.xy(), tgt::vec2::zero, tgt::vec2(screensize), -tgt::vec2::one, tgt::vec2::one);
    tgt::vec2 ur = mapRange(pos.xy() + size, tgt::vec2::zero,  tgt::vec2(screensize), -tgt::vec2::one, tgt::vec2::one);

    // Draw background.
    IMode.color(tgt::vec4(0.2f, 0.2f, 0.2f, 0.5f));
    IMode.begin(tgt::ImmediateMode::QUADS);
    IMode.vertex(ll.x, ll.y);
    IMode.vertex(ur.x, ll.y);
    IMode.vertex(ur.x, ur.y);
    IMode.vertex(ll.x, ur.y);
    IMode.end();
    IMode.color(tgt::vec4::one);

    // Draw font.
    font.render(pos, tooltip, screensize);

    glEnable(GL_DEPTH_TEST);
}

bool TimeseriesPlot::isReady() const {

    if(!timeseriesInport_.isReady()) {
        setNotReadyErrorMessage("No time series connected");
        return false;
    }

    // Note: Similarity Matrix is optional.

    return true;
}

void TimeseriesPlot::renderingPass(bool picking) {

    switch(numDimensions_.get()) {
    case 1:
        renderEmbedding1D(picking);
        break;
    case 2:
        renderEmbedding2D(picking);
        break;
    case 3:
        renderEmbedding3D(picking);
        break;
    default:
        // No visualization available
        break;
    }
}

void TimeseriesPlot::renderEmbedding1D(bool picking) {

    // Retrieve dataset.
    const TimeSeriesList* timeserieslist = timeseriesInport_.getData();

    // Retrieve selected embedding.
    const Embedding& embedding = embedding_;

    tgt::vec2 timeRange = tgt::vec2(timeserieslist->getOverallStartTime(), timeserieslist->getOverallEndTime());
    if(renderTimeSelection_.get()) {
        timeRange = firstSelectedTimeInterval_.get();
    }

    for(int seriesIdx : renderingOrder_) {

        glLineWidth((subSelection_.count(seriesIdx) != 0) ? SELECTED_LINE_WIDTH : UNSELECTED_LINE_WIDTH);

        const TimeSeries& series = timeserieslist->getTimeSeries(seriesIdx);
        size_t numTimeSteps = series.getNumberOfTimeSteps();
        int eigenValueIdx = principleComponent_.get() - 1;
        const auto& vertices = embedding.nVectors_.at(seriesIdx);

        // In case we have a single time step, we draw it across the whole range,
        // since it doesn't change. This could (and should!) be improved, however,
        // such that it becomes clear at which t the time step is recorded.
        if(numTimeSteps == 1) {
            IMode.begin(tgt::ImmediateMode::FAKE_LINES);
            IMode.color(getColor(seriesIdx, 0, picking));
            const int segments = 80;
            for(int i=0; i<segments; i+=2) {
                float x0 = mapRange(i+0, 0, segments-1, -1.0f, 1.0f);
                float x1 = mapRange(i+1, 0, segments-1, -1.0f, 1.0f);
                IMode.vertex(tgt::vec2(x0, vertices[0][eigenValueIdx]));
                IMode.vertex(tgt::vec2(x1, vertices[0][eigenValueIdx]));
            }
            IMode.end();
        }
        else {
            IMode.begin(tgt::ImmediateMode::FAKE_LINE_STRIP);
            for (size_t j = 0; j < numTimeSteps; j++) {
                float colorSaturation = 1.0f;
                if(!picking && (series.getTimeSeriesStep(j).time_ < timeRange.x || series.getTimeSeriesStep(j).time_ > timeRange.y)) {
                    colorSaturation = 0.25f;
                }
                float t = mapRange(series.getTimeSeriesStep(j).time_, timeserieslist->getOverallStartTime(), timeserieslist->getOverallEndTime(), -1.0f, 1.0f);
                IMode.color(getColor(seriesIdx, j, picking) * colorSaturation + tgt::vec3(1.0f - colorSaturation));
                IMode.vertex(tgt::vec2(t, vertices[j][eigenValueIdx]));
            }
            IMode.end();
        }
    }

    if(!picking && renderTimeSelection_.get()) {
        tgt::vec2 mappedTimeRange = mapRange(firstSelectedTimeInterval_.get(), tgt::vec2(timeserieslist->getOverallStartTime()), tgt::vec2(timeserieslist->getOverallEndTime()), -tgt::vec2::one, tgt::vec2::one);

        glLineWidth(3.0f);
        IMode.color(tgt::vec3::zero);
        IMode.begin(tgt::ImmediateMode::LINES);
        IMode.vertex(tgt::vec2(mappedTimeRange.x, -1.0f));
        IMode.vertex(tgt::vec2(mappedTimeRange.x,  1.0f));
        IMode.vertex(tgt::vec2(mappedTimeRange.y, -1.0f));
        IMode.vertex(tgt::vec2(mappedTimeRange.y,  1.0f));
        IMode.end();
    }
}

void TimeseriesPlot::renderEmbedding2D(bool picking) {

    // Retrieve dataset.
    const TimeSeriesList* list = timeseriesInport_.getData();

    for (int seriesIdx : renderingOrder_) {

        glLineWidth((subSelection_.count(seriesIdx) != 0) ? SELECTED_LINE_WIDTH : UNSELECTED_LINE_WIDTH);

        size_t numTimeSteps = list->getTimeSeries(seriesIdx).getNumberOfTimeSteps();
        const auto& vertices = embedding_.nVectors_.at(seriesIdx);

        IMode.begin(tgt::ImmediateMode::FAKE_LINE_STRIP);
        for(size_t j=0; j<numTimeSteps; j++) {
            IMode.color(getColor(seriesIdx, j, picking));
            IMode.vertex(tgt::vec2(vertices[j][0], vertices[j][1]));
        }
        IMode.end();
//        IMode.begin(tgt::ImmediateMode::POINTS);
//        for(size_t j=0; j<numTimeSteps; j++) {
//            IMode.color(getColor(seriesIdx, j, picking));
//            IMode.vertex(tgt::vec2(vertices[j][0], vertices[j][1]));
//        }
//        IMode.end();

        if((!picking && renderTimeSelection_.get()) || numTimeSteps == 1) {
            size_t selectedTimeStep = list->getTimeSeries(seriesIdx).getTimeStepIndex(firstSelectedTimeInterval_.get().x);
            tgt::vec3 position(vertices[selectedTimeStep][0], vertices[selectedTimeStep][1], 0.0f);
            tgt::vec3 color = (numTimeSteps == 1) ? getColor(seriesIdx, selectedTimeStep, picking) : tgt::vec3::one;
            renderTimeStepSelection(seriesIdx, selectedTimeStep, position, color);
        }
    }
}

void TimeseriesPlot::renderEmbedding3D(bool picking) {

    // Retrieve dataset.
    const TimeSeriesList* list = timeseriesInport_.getData();

    // If using 3D visualisation, use camera interaction.
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadMatrix(camera_.get().getProjectionMatrix(outport_.getSize()));
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadMatrix(camera_.get().getViewMatrix());

    tgt::vec3 scale = tgt::vec3::one;
//    if(scaleToMagnitude_.get()) {
//        // Scale each axis to it's eigenvalues size.
//        scale = tgt::vec3::fromPointer(&embedding_.eigenvalues_[0]);
//        scale /= tgt::vec3(embedding_.eigenvalues_[0]);
//    }

    for (int seriesIdx : renderingOrder_) {

        glLineWidth((subSelection_.count(seriesIdx) != 0) ? SELECTED_LINE_WIDTH : UNSELECTED_LINE_WIDTH);

        size_t numTimeSteps = list->getTimeSeries(seriesIdx).getNumberOfTimeSteps();
        const auto& vertices = embedding_.nVectors_.at(seriesIdx);

        IMode.begin(tgt::ImmediateMode::FAKE_LINE_STRIP);
        for(size_t j=0; j<numTimeSteps; j++) {
            IMode.color(getColor(seriesIdx, j, picking));
            IMode.vertex(tgt::vec3::fromPointer(&vertices[j][0]) * scale);
        }
        IMode.end();

        if((!picking && renderTimeSelection_.get()) || numTimeSteps == 1) {
            size_t selectedTimeStep = list->getTimeSeries(seriesIdx).getTimeStepIndex(firstSelectedTimeInterval_.get().x);
            tgt::vec3 position = tgt::vec3::fromPointer(&vertices[selectedTimeStep][0])*scale;
            tgt::vec3 color = (numTimeSteps == 1) ? getColor(seriesIdx, selectedTimeStep, picking) : tgt::vec3::one;
            renderTimeStepSelection(seriesIdx, selectedTimeStep, position, color);
        }
    }

    // restore matrices
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();
}

void TimeseriesPlot::renderTimeStepSelection(size_t seriesIdx, size_t timeStepIdx, const tgt::vec3& position, const tgt::vec3& color) const {

    // Skip rendering, if not visible anyways.
    if(sphereRadius_.get() <= std::numeric_limits<float>::epsilon())
        return;

    MatStack.pushMatrix();
    MatStack.translate(position);
    MatStack.scale(tgt::vec3(sphereRadius_.get()));

    const TimeSeriesList* list = timeseriesInport_.getData();
    size_t numTimeSteps = list->getTimeSeries(seriesIdx).getNumberOfTimeSteps();

    bool timeStepAvailable = timeStepIdx < numTimeSteps;
    if(!timeStepAvailable) {
        IMode.color(tgt::vec4(color, 0.5f));
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_BLEND);
    }
    else {
        IMode.color(tgt::vec4(color, 1.0f));
    }

    sphere_.render(GL_TRIANGLES);

    if(!timeStepAvailable) {
        glDisable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ZERO);
    }

    MatStack.popMatrix();
}

tgt::vec3 TimeseriesPlot::getColor(size_t seriesIdx, size_t timeStepIdx, bool picking) const {

    const TimeSeriesList* list = timeseriesInport_.getData();
    const TimeSeries& series = list->getTimeSeries(seriesIdx);

    float ts = static_cast<float>(timeStepIdx) / series.getNumberOfTimeSteps();
    if (picking)
        return tgt::vec3(static_cast<float>(seriesIdx) / list->getNumberOfTimeSeries(), ts, 1.0f);

    tgt::vec3 colorPosition = series.getPosition();// use position as color
    colorPosition.x = ((float) colorPosition.x)/list->getVolumeDimensions().x;
    colorPosition.y = ((float) colorPosition.y)/list->getVolumeDimensions().y;
    colorPosition.z = ((float) colorPosition.z)/list->getVolumeDimensions().z;

    int tempAnomalyIndex = 0;
    int angleIndex = 0;
    int numFields = series.getTimeSeriesStep(timeStepIdx).fieldValues_.size();
    for(int i = 0; i < numFields; i++) {
        if(list->getFieldName(i)=="angle") {
            angleIndex = i;
        }
        if(list->getFieldName(i)=="temperature anomaly") {
            tempAnomalyIndex = i;
        }
    }
    float angleVal = series.getTimeSeriesStep(timeStepIdx).fieldValues_[angleIndex];
    tgt::vec2 angleRange = tgt::vec2(0,0);
    if(list->getFieldName(angleIndex)=="angle") {
        angleRange = list->getRange("angle");
    }
    angleVal = angleVal * (angleRange.y-angleRange.x) + angleRange.x;
    float tempAnomalyVal = series.getTimeSeriesStep(timeStepIdx).fieldValues_[tempAnomalyIndex];
    tgt::vec2 tempAnomalyRange = tgt::vec2(0,0);
    if(list->getFieldName(tempAnomalyIndex)=="temperature anomaly") {
        tempAnomalyRange = list->getRange("temperature anomaly");
    }
    tempAnomalyVal = tempAnomalyVal * (tempAnomalyRange.y-tempAnomalyRange.x) + tempAnomalyRange.x;

    switch (colorCoding_.getValue()) {
    case COLOR_MEMBER:
        return colorPosition;
    case COLOR_TIMESTEP:
        return  (1.0f - ts) * FIRST_TIME_STEP_COLOR + ts * LAST_TIME_STEP_COLOR;
    case COLOR_MEMBER_AND_TIMESTEP:
        return (1.0f - ts) * colorPosition + ts * FADE_OUT_COLOR;
//    case COLOR_DURATION:
//    {
//        const Statistics& stats = member.getTimeStepDurationStats();
//        if(std::abs(dataset->getMembers()[seriesIdx].getTimeSteps()[timeStepIdx].getDuration() - stats.getMean()) > stats.getStdDev()) {
//            ts = mapRange(dataset->getMembers()[seriesIdx].getTimeSteps()[timeStepIdx].getDuration(), stats.getMin(), stats.getMax(), 0.0f, 1.0f);
//            return (1.0f - ts) * MIN_DURATION_COLOR + ts * MAX_DURATION_COLOR;
//        }
//        return MIN_DURATION_COLOR;
//    }
    case COLOR_FIELD:
    {
        int index = (colorCoding_.getSelectedIndex()-3);
        float val = series.getTimeSeriesStep(timeStepIdx).fieldValues_[index];
        tgt::vec3 color = getColorForValue(val);
        return color;
        //return  color;//(1.0f - val) * FIRST_TIME_STEP_COLOR + val * LAST_TIME_STEP_COLOR;
    }
//    case COLOR_THRESHOLD:
//    {
//        int index = (colorCoding_.getSelectedIndex()-3)/2;
//        float val = series.getTimeSeriesStep(timeStepIdx).fieldValues_[index];
//        tgt::vec2 range = parameterRange_.get();
//        if(val < range.x)
//            return SMALLER_THRESHOLD_COLOR;
//        else if (val > range.y)
//            return LARGER_THRESHOLD_COLOR;
//        else
//            return BETWEEN_THRESHOLD_COLOR;
//    }
    case COLOR_SLAB:
    {
        if(angleVal > angleRange_.get().y &&  tempAnomalyVal < tempAnomalyRange_.get().x) {
            return SLAB_COLOR;
        }
        return NOTHING_COLOR;
    }
    case COLOR_PLUME:
    {
        if(angleVal < angleRange_.get().x && tempAnomalyVal > tempAnomalyRange_.get().y) {
            return PLUME_COLOR;
        }
        return NOTHING_COLOR;
    }
    case COLOR_SLAB_AND_PLUME:
    {
        if(angleVal > angleRange_.get().y && tempAnomalyVal < tempAnomalyRange_.get().x) {
            return SLAB_COLOR;
        }
        if(angleVal < angleRange_.get().x && tempAnomalyVal > tempAnomalyRange_.get().y) {
            return PLUME_COLOR;
        }
        return NOTHING_COLOR;
    }
    default:
        return tgt::vec3::one;
    }
}

void TimeseriesPlot::mouseEvent(tgt::MouseEvent* e) {

    // We always accept the event and force a redraw.
    e->accept();
    invalidate();

    // Right click resets subselection to rendered members.
    if(e->button() == tgt::MouseEvent::MOUSE_BUTTON_MIDDLE && e->action() == tgt::MouseEvent::PRESSED
    && e->modifiers() == tgt::MouseEvent::MODIFIER_NONE) {
        //subSelection_.clear(); // Old behavior.
        subSelection_ = std::set<int>(renderedMembers_.get().begin(), renderedMembers_.get().end());
    }
    // Otherwise, we look for a hit.
    else {

        // Reset last hit.
        lastHit_ = boost::none;

        int x = e->x();
        int y = e->y();

        // In 1D, and 2D case, we use margins to show axis and labels.
        if (numDimensions_.get() < 3) {
            x = tgt::clamp(x, MARGINS.x, e->viewport().x - MARGINS.x);
            y = tgt::clamp(y, MARGINS.y, e->viewport().y - MARGINS.y);
        }

        // Not inside margins.
        if (x != e->x() || (y != e->y() && numDimensions_.get() > 1)) {
            return;
        }

        const TimeSeriesList* list = timeseriesInport_.getData();
        if(y != e->y() && numDimensions_.get() == 1) {
            float t0 = firstSelectedTimeInterval_.get().x;
            float t1 = firstSelectedTimeInterval_.get().y;

            // Select a time interval outside of the axis
            if(e->action() == tgt::MouseEvent::PRESSED) {
                t0 = mapRange(x, MARGINS.x, e->viewport().x - MARGINS.x, list->getOverallStartTime(), list->getOverallEndTime());
            }
            else if(e->action() == tgt::MouseEvent::RELEASED) {
                t1 = mapRange(x, MARGINS.x, e->viewport().x - MARGINS.x, list->getOverallStartTime(), list->getOverallEndTime());
            } // We don't support MOTION atm., since it would invalidate too often and would trigger linked filters.
            else {
                return;
            }

            // Fix range.
            if(t0 > t1) {
                std::swap(t0, t1);
            }

            // Add small offset, if necessary to include at least a single time step in the selected range.
//            if(t0 < t1 + dataset->getMinTimeStepDuration()) {
//                t1 += dataset->getMinTimeStepDuration();
//            }

            // Apply selection for both ranges.
            firstSelectedTimeInterval_.set(tgt::vec2(t0, t1));
            secondSelectedTimeInterval_.set(tgt::vec2(t0, t1));
            return;
        }

        // Handle margins.
        RenderTarget* target = pickingBuffer_.getRenderTarget();
        tgt::ivec2 pixel;
        if (numDimensions_.get() < 3) {
            pixel = tgt::ivec2(
                    mapRange(tgt::vec2(x, y), tgt::vec2(MARGINS), tgt::vec2(outport_.getSize() - MARGINS),
                             tgt::vec2::zero, tgt::vec2(target->getSize())));
        }
        else {
            pixel = mapRange(tgt::ivec2(x, y), tgt::ivec2::zero, outport_.getSize(), tgt::ivec2::zero,
                             target->getSize());
        }

        tgt::vec4 texel = target->getColorAtPos(tgt::ivec2(pixel.x, target->getSize().y - pixel.y - 1));
        if (texel.z != 1.0f) {
            return; // No hit - preserve last selection.
        }

        // Calculate member index.
        //const std::vector<TimeSeriesList>& members = dataset->getMembers();
        size_t numMembers = list->getNumberOfTimeSeries();
        int r = tgt::clamp<int>(std::round(texel.r * numMembers), 0, numMembers - 1);

        // Calculate time step index.
        size_t numTimeSteps = list->getTimeSeries(r).getNumberOfTimeSteps();
        int t = tgt::clamp<int>(std::round(texel.g * numTimeSteps), 0, numTimeSteps - 1);

        // Update last hit.
        lastHit_ = Hit{x, y, r, t};

        // Handle simple selection, for both main selection (left button) and reference selection (right button).
        if((e->button() == tgt::MouseEvent::MOUSE_BUTTON_LEFT || e->button() == tgt::MouseEvent::MOUSE_BUTTON_RIGHT)
        && e->action() == tgt::MouseEvent::PRESSED && e->modifiers() == tgt::MouseEvent::MODIFIER_NONE) {

            // Reset subselection.
            //subSelection_.clear();

            std::vector<int> memberIndices;
            memberIndices.push_back(r);
            subSelection_.insert(r);

            // add point for member r to list
            tgt::vec3 point = list->getTimeSeries(r).getPosition();
            // Set point coordinates to spherical coordinates
            //transformedPoints.push_back(point);
            float radius = volumeDimension_.get()/2;
            float theta = point.x * PI/180;
            float phi = point.z * PI/180 - PI;
            point.x = radius * std::sin(theta) * std::cos(phi) + radius;
            point.y = radius * std::sin(theta) * std::sin(phi) + radius;
            point.z = radius * std::cos(theta) + radius;
            position_.set(point);
            //std::cout << point << std::endl;
            //point.x = radius*std::sin()
//            PointListGeometryVec3* transformedPositions = new PointListGeometryVec3();
//            transformedPositions->setData(transformedPoints);
//            pointListOutport_.setData(transformedPositions, true);

            const TimeSeriesStep& timeStep = list->getTimeSeries(r).getTimeSeriesStep(t);
            //TODO: Duration is workaround
            float duration = 1;
            if(t < list->getTimeSeries(r).getNumberOfTimeSteps())
                duration = list->getTimeSeries(r).getTimeSeriesStep(t+1).time_ - timeStep.time_;
            float lower = std::floor(timeStep.time_ * 100.0f) / 100.0f;
            float upper = std::ceil((timeStep.time_ + duration) * 100.0f) / 100.0f;
            if (e->button() == tgt::MouseEvent::MOUSE_BUTTON_LEFT) {
                firstSelectedMember_.setSelectedRowIndices(memberIndices);
                firstSelectedTimeInterval_.set(tgt::vec2(lower, upper));
            }
            else if (e->button() == tgt::MouseEvent::MOUSE_BUTTON_RIGHT) {
                secondSelectedMember_.setSelectedRowIndices(memberIndices);
                secondSelectedTimeInterval_.set(tgt::vec2(lower, upper));
            }

            return;
        }

        // Add (CTRL) or remove (ALT) members from subselection.
        if(e->button() == tgt::MouseEvent::MOUSE_BUTTON_NONE) {
            auto iter = subSelection_.find(r);
            if(e->modifiers() == tgt::MouseEvent::CTRL && iter == subSelection_.end()) {
                subSelection_.insert(r);
            }
            else if(e->modifiers() == tgt::MouseEvent::ALT && iter != subSelection_.end()) {
                subSelection_.erase(iter);
            }
            return;
        }

        // Modify rendering order as simple workaround for occlusion.
        if (e->button() == tgt::MouseEvent::MOUSE_BUTTON_MIDDLE && e->action() == tgt::MouseEvent::PRESSED
        && e->modifiers() == tgt::MouseEvent::SHIFT) {
            // Push selected member to the front of the rendering order.
            renderingOrder_.erase(std::find(renderingOrder_.begin(), renderingOrder_.end(), r));
            renderingOrder_.push_front(r);

            return;
        }
    }
}

void TimeseriesPlot::onEvent(tgt::Event* e) {
    tgt::MouseEvent* event = dynamic_cast<tgt::MouseEvent*>(e);
    if (event) {
        mouseEvent(event);
    }

    RenderProcessor::onEvent(e);
}

void TimeseriesPlot::adjustToTimeserieslist() {

    subSelection_.clear();
    renderedMembers_.reset();
    firstSelectedMember_.reset();
    secondSelectedMember_.reset();
    colorCoding_.setOptions(std::deque<Option<ColorCoding>>());
    //colorCoding_.reset();
    calculateButton_.setReadOnlyFlag(true);
    embedding_.eigenvalues_.clear();
    embedding_.nVectors_.clear();
    embedding_.names_.clear();

    // Check for data set.
    const TimeSeriesList* list = timeseriesInport_.getData();
    if (!list)
        return;

    // Check for similarity matrices.
    //const SimilarityMatrixList* similarityMatrices = similarityMatrixInport_.getData();
    //if (!similarityMatrices)
    //    return;

    numEigenvalues_.setMinValue(MAX_NUM_DIMENSIONS);
    numEigenvalues_.setMaxValue(5*MAX_NUM_DIMENSIONS);

    std::vector<int> seriesIndices;
    for (const TimeSeries& series : list->getAllTimeSeries()) {
        std::string seriesPosition = std::to_string(series.getPosition().x) + " " + std::to_string(series.getPosition().y) + " " + std::to_string(series.getPosition().z);
        tgt::vec3 positionColor = series.getPosition()/list->getVolumeDimensions();
        renderedMembers_.addRow(seriesPosition, positionColor);
        firstSelectedMember_.addRow(seriesPosition, positionColor);
        secondSelectedMember_.addRow(seriesPosition, positionColor);
        subSelection_.insert(static_cast<int>(seriesIndices.size()));
        seriesIndices.push_back(static_cast<int>(seriesIndices.size()));
    }
    renderedMembers_.setSelectedRowIndices(seriesIndices);
    firstSelectedMember_.setSelectedRowIndices(seriesIndices);
    secondSelectedMember_.setSelectedRowIndices(seriesIndices);

    firstSelectedTimeInterval_.setMinValue(list->getOverallStartTime());
    firstSelectedTimeInterval_.setMaxValue(list->getOverallEndTime());
    firstSelectedTimeInterval_.set(tgt::vec2(list->getOverallStartTime(), list->getOverallEndTime()));

    secondSelectedTimeInterval_.setMinValue(list->getOverallStartTime());
    secondSelectedTimeInterval_.setMaxValue(list->getOverallEndTime());
    secondSelectedTimeInterval_.set(tgt::vec2(list->getOverallStartTime(), list->getOverallEndTime()));

    // Color Coding options
    colorCoding_.addOption("member", "Only Member", COLOR_MEMBER);
    colorCoding_.addOption("timeStep", "Only Time Step", COLOR_TIMESTEP);
    colorCoding_.addOption("memberAndTimeStep", "Member and Time Step", COLOR_MEMBER_AND_TIMESTEP);
    for(auto fieldname : list->getFieldNames()) {
        colorCoding_.addOption(fieldname, fieldname, COLOR_FIELD);
        //colorCoding_.addOption(fieldname+"t", fieldname+" Threshold", COLOR_THRESHOLD);
    }
    colorCoding_.addOption("slab", "Slabs", COLOR_SLAB);
    colorCoding_.addOption("plume", "Plumes", COLOR_PLUME);
    colorCoding_.addOption("slabsAndPlume", "Slabs and Plumes", COLOR_SLAB_AND_PLUME);

    // Try to load plot data, if already set.
    if(!loadFileDialog_.get().empty()) {
        loadEmbeddings();
    }

    calculateButton_.setReadOnlyFlag(false);

    if (autoCalculate_.get()) {
        createEmbeddings();
    }
}

void TimeseriesPlot::createEmbeddings() {

    if(subSelection_.empty()) {
        LERROR("No member selected");
        return;
    }

    size_t numSeries = timeseriesInport_.getData()->getNumberOfTimeSeries();
    if(subSelection_.size() == numSeries) {
        LINFO("Calculating for whole data set..");
    }
    else {
        LINFO("Calculating for subset..");
    }

    std::vector<std::string> names;
    for(const auto& series : timeseriesInport_.getData()->getAllTimeSeries()) {
        std::string seriesPosition = std::to_string(series.getPosition().x) + " " + std::to_string(series.getPosition().y) + " " + std::to_string(series.getPosition().z);
        names.push_back(seriesPosition);
    }

    calculateButton_.setReadOnlyFlag(true);

    //const SimilarityMatrixList* matrices = similarityMatrixInport_.getData();
    std::vector<std::string> fieldNames = timeseriesInport_.getData()->getFieldNames();

    setProgress(0.0f);
    //for (size_t i=0; i<fieldNames.size(); i++) {
    SubtaskProgressReporter progressReporter(*this, tgt::vec2(0, 1));

    // Get distance matrix for current field.
    //const SimilarityMatrix& distanceMatrix = matrices->getSimilarityMatrix(fieldNames[i]);

    // Compute Principal components and corresponding eigenvectors.
    embedding_ = createEmbedding(timeseriesInport_.getData(), progressReporter);

    // Add member name.
    embedding_.names_ = names;
    //}

    // Finally, output eigen values to allow an eigen value analysis.
    outputEigenValues();

    // Update rendering order.
    renderedMembers_.setSelectedRowIndices(std::vector<int>(subSelection_.rbegin(), subSelection_.rend()));
    //renderingOrder_.assign(subSelection_.rbegin(), subSelection_.rend());


    // Done.
    calculateButton_.setReadOnlyFlag(false);
    setProgress(1.0f);
    invalidate();
}

TimeseriesPlot::Embedding TimeseriesPlot::createEmbedding(const TimeSeriesList* timeseries, ProgressReporter& progressReporter, float epsilon) const {
    using namespace Eigen;

    auto start = std::chrono::steady_clock::now();

    progressReporter.setProgress(0.1f);

    const std::vector<TimeSeries>& series = timeseriesInport_.getData()->getAllTimeSeries();

    const size_t numDimensions = numEigenvalues_.get();
    const int numIterations = numIterations_.get();
    const size_t numMembers = series.size();

    Embedding embedding;

    size_t numPoints = 0;
    for(int seriesIdx : subSelection_) {
        size_t numTimeSteps = series[seriesIdx].getNumberOfTimeSteps();
        numPoints += numTimeSteps;
        embedding.nVectors_[seriesIdx] = std::vector<std::vector<float>>(numTimeSteps, std::vector<float>(numDimensions, 0.0f));
    }

    MatrixXf result = MatrixXf::Zero(numPoints, numDimensions);
    MatrixXf EigSq = MatrixXf::Zero(numDimensions, numDimensions);
    MatrixXf PMatrix(numPoints, numPoints);

    // Init datastructures.
    size_t offsetA = 0;
    size_t positionA = 0;

    auto end = std::chrono::steady_clock::now();
    std::cout << "Time for preparations: " << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << "ms" << std::endl;

    // Iterate each member, call it A.
    for(size_t memberIdxA=0; memberIdxA < numMembers; memberIdxA++) {
        size_t offsetB = 0;
        size_t positionB = 0;
        size_t numTimeStepsA = series[memberIdxA].getNumberOfTimeSteps();
        if(subSelection_.count(memberIdxA) != 0) { // Do we consider member A?
            // Again iterate each member, call it B. Now looking at pairs of member A and B.
            for(size_t memberIdxB=0; memberIdxB<=memberIdxA; memberIdxB++) {
                size_t numTimeStepsB = series[memberIdxB].getNumberOfTimeSteps();
                if (subSelection_.count(memberIdxB) != 0) { // Do we consider member B?
                    // Iterate time steps of member A.
                    for (size_t i = 0; i < numTimeStepsA; i++) {
                        // Iterate time steps of member B.
                        for (size_t j = 0; j < numTimeStepsB; j++) {
                            std::vector<float> pointA = timeseries->getTimeSeries(memberIdxA).getTimeSeriesStep(i).fieldValues_;
                            std::vector<float> pointB = timeseries->getTimeSeries(memberIdxB).getTimeSeriesStep(j).fieldValues_;
                            // squared distance
                            float v = 0;//distanceMatrix(i + offsetA, j + offsetB);
                            for(size_t k = 0; k < timeseries->getNumberOfComponentsPerField(); k++) {
                                v += (pointB[k]-pointA[k])*(pointB[k]-pointA[k]);
                            }
                            PMatrix(i + positionA, j + positionB) = PMatrix(j + positionB, i + positionA) = v;
                        }
                    }
                    positionB+=numTimeStepsB;
                }
                offsetB+=numTimeStepsB;
            }
            positionA+=numTimeStepsA;
        }
        offsetA+=numTimeStepsA;
    }

    auto end2 = std::chrono::steady_clock::now();
    std::cout << "Time for filling matrices: " << std::chrono::duration_cast<std::chrono::milliseconds>(end2-end).count() << "ms" << std::endl;

    MatrixXf JMatrix = MatrixXf::Identity (numPoints, numPoints) -
                       (1.0f / numPoints) * MatrixXf::Ones (numPoints, numPoints);

    MatrixXf BMatrix = -0.5f*JMatrix*PMatrix*JMatrix;

    VectorXf* eigenVectors = new VectorXf[numDimensions];

    end = std::chrono::steady_clock::now();
    std::cout << "Time for some matrix calculations: " << std::chrono::duration_cast<std::chrono::milliseconds>(end-end2).count() << "ms" << std::endl;

    progressReporter.setProgress(0.1f);
    //for(size_t i=0; i < numDimensions; i++) {

//        VectorXf& eigenVector = eigenVectors[i];
//        eigenVector = VectorXf::Ones(numPoints);
//        float eigenValue = 0.0;
//
//        VectorXf PrVector = VectorXf::Zero(numPoints);
//        VectorXf TVector = VectorXf::Ones(numPoints);
//
//        MatrixXf EMVector(numPoints, 1);

        //auto startIter = std::chrono::steady_clock::now();

//        for (int iter=0; iter < numIterations; iter++) {
//
//            EMVector.col(0) = eigenVector;
//            VectorXf TempVector = BMatrix * EMVector;
//            eigenVector = TempVector;
//            for(size_t j=0; j < i; j++)
//                eigenVector -= eigenVectors[j] * (eigenVectors[j].dot(TempVector));
//            eigenValue = eigenVector.norm();
//            eigenVector.normalize();
//            TVector = eigenVector - PrVector;
//            PrVector = eigenVector;
//
//            if(TVector.norm() <= epsilon) {
//                break;
//            }
//
//            //auto endIter = std::chrono::steady_clock::now();
//            //std::cout << "Time for iteration: " << std::chrono::duration_cast<std::chrono::milliseconds>(endIter-startIter).count() << "ms" << std::endl;
//        }
//
//        // Don't continue calculating when eigenvalues get too small in relation to biggest eigenvalue.
//        if(i >= static_cast<size_t>(numDimensions_.get()) && // We do have calculated at least as many PCs as we want to display.
//        (eigenValue <= std::numeric_limits<float>::epsilon() || // Eigenvalue is de-facto zero.
//        (i > 0 && eigenValue < embedding.eigenvalues_[i - 1]))) {// EV must not be smaller than predecessor.
//            std::cout << "Eigenvalue small" << std::endl;
//            std::cout << eigenValue << " " << embedding.eigenvalues_[i - 1] << std::endl;
//            std::cout << "<" << (eigenValue < embedding.eigenvalues_[i - 1]) << std::endl;
//            std::cout << i << std::endl;
//            break;
//        }
//
//        EMVector.col(0) = eigenVector;
//        result.col(i) = eigenVector;
//        EigSq(i, i) = std::sqrt(eigenValue);
//        embedding.eigenvalues_.push_back(eigenValue);
//
//        progressReporter.setProgress(0.1f + 0.8f * i / numDimensions);
//    }
//    delete [] eigenVectors;
    SelfAdjointEigenSolver<MatrixXf> es(BMatrix);
    for(size_t i=0; i < numDimensions; i++) {
        int numEigenvalues = es.eigenvalues().size();
        float eigenvalue = es.eigenvalues()[numEigenvalues-i-1];
        embedding.eigenvalues_.push_back(eigenvalue);
        result.col(i) = es.eigenvectors().real().col(numEigenvalues-1-i);
        EigSq(i, i) = std::sqrt(eigenvalue);
    }
    end2 = std::chrono::steady_clock::now();
    std::cout << "Time for eigenvalue calculations: " << std::chrono::duration_cast<std::chrono::milliseconds>(end2-end).count() << "ms" << std::endl;

    // Now get the resulting matrix.
    result = result * EigSq;

    float maxValue = result(0, 0);
    float minValue = result(0, 0);
    for (size_t i=0; i < embedding.eigenvalues_.size(); i++) {
        for (size_t j = 1; j < numPoints; j++) {
            if (maxValue < result(j, i))
                maxValue = result(j, i);
            else if (minValue > result(j, i))
                minValue = result(j, i);
        }
    }

    for (size_t i=0; i < embedding.eigenvalues_.size(); i++) {
        // Interpret as member and time steps.
        size_t j=0;
        for (int memberIdx : subSelection_) {
            for(size_t t=0; t<series[memberIdx].getNumberOfTimeSteps(); t++) {
                float value = mapRange(result(j, i), minValue, maxValue, -1.0f, 1.0f);
                embedding.nVectors_[memberIdx][t][i] = value;
                j++;
            }
        }
    }

    progressReporter.setProgress(1.0f);

    return embedding;
}

void TimeseriesPlot::outputEigenValues() {

    if (embedding_.eigenvalues_.empty()) {
        eigenValueOutport_.setData(nullptr);
        return;
    }

    PlotData* data = new PlotData(1, 1);

    data->setColumnLabel(0, "Index");
    data->setColumnLabel(1, "Eigenvalue");

    for(size_t i = 0; i < embedding_.eigenvalues_.size(); i++) {
        std::vector<PlotCellValue> values;
        values.push_back(PlotCellValue(i+1));
        values.push_back(PlotCellValue(static_cast<plot_t>(embedding_.eigenvalues_[i])));
        data->insert(values);
    }
    eigenValueOutport_.setData(data, true);
}

void TimeseriesPlot::renderedMembersChanged() {
    renderingOrder_.clear();

    if(embedding_.eigenvalues_.empty())
        return;

    for(int memberIdx : renderedMembers_.getSelectedRowIndices()) {
        // The first field is representative for all fields here.
        // We just want to check if the embedding was already calculated for the member.
        if(embedding_.nVectors_.count(memberIdx)) {
            renderingOrder_.push_back(memberIdx);
        }
    }
}

void TimeseriesPlot::saveEmbeddings() {
    if (saveFileDialog_.get().empty()) {
        LWARNING("no filename specified");
        return;
    }
    if(embedding_.eigenvalues_.empty()) {
        LWARNING("No embeddings calculated.");
        return;
    }

    // write file
    try {
        std::ofstream outFile;
        outFile.open(saveFileDialog_.get().c_str());
        LINFO("Writing embeddings to file " << saveFileDialog_.get());

        JsonSerializer s;
        Serializer serializer(s);
        serializer.serialize("embeddings", embedding_);
        s.write(outFile, true);

        outFile.close();
        LINFO("Saving " << saveFileDialog_.get() << " was successful.");
    }
    catch(tgt::Exception& e) {
        VoreenApplication::app()->showMessageBox("saving MDS Plot failed", e.what(), true);
        LERROR(e.what());
        saveFileDialog_.set("");
    }
}

void TimeseriesPlot::loadEmbeddings() {
    if (loadFileDialog_.get().empty()) {
        LWARNING("no filename specified");
        return;
    }
    if (!timeseriesInport_.isReady()) {
        LWARNING("no time series connected");
        return;
    }

    std::string absPath = loadFileDialog_.get();
    try {
        std::ifstream inFile;
        inFile.open(absPath.c_str());

        JsonDeserializer d;
        d.read(inFile);
        Deserializer deserializer(d);
        deserializer.deserialize("embeddings", embedding_);

        inFile.close();
        outputEigenValues();
        LINFO("Reading embeddings from file " << absPath << " was successful");
    }
    catch(tgt::Exception& e) {
        VoreenApplication::app()->showMessageBox("loading embeddings failed", e.what(), true);
        LERROR(e.what());
        loadFileDialog_.set("");
    }
}

tgt::vec3 TimeseriesPlot::getColorForValue(float val) const {
    auto keys = transferFunc_.get()->getKeys();
    if(keys.size() == 0)
        return tgt::vec3::one;
    int index = 0;
    tgt::vec2 range = transferFunc_.get()->getDomain();
    float diff = range.y-range.x;
    while(index < keys.size() && val > keys[index]->getIntensity()*diff+range.x) {
        index++;
    }
    if(index == 0) {
        return keys[index]->getColorR().xyz();
    }
    else if(index == keys.size()) {
        return keys[index-1]->getColorR().xyz();
    }
    else {
        float distanceUp = keys[index]->getIntensity()*diff+range.x - val;
        float distanceDown = val - keys[index - 1]->getIntensity()*diff+range.x;
        float distance = distanceUp + distanceDown;
        tgt::vec3 color = distanceUp / distance * (tgt::vec3) keys[index - 1]->getColorR().xyz() +
                          distanceDown / distance * (tgt::vec3) keys[index]->getColorL().xyz();
        return color/255.f;
    }
}

tgt::vec2 TimeseriesPlot::getRange(size_t fieldIdx) {
    tgt::vec2 range(1.f,0.f);
    for(auto series : timeseriesInport_.getData()->getAllTimeSeries()) {
        for(auto step : series.getTimeSeriesSteps()) {
            range.x = std::min(range.x, step.fieldValues_[fieldIdx]);
            range.y = std::max(range.y, step.fieldValues_[fieldIdx]);
        }
    }
    return range;
}

} // namespace
