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

#include "similarityplot.h"

#include "../utils/colorpool.h"
#include "../utils/ensemblehash.h"
#include "../utils/utils.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/interaction/camerainteractionhandler.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"

#include "modules/plotting/datastructures/plotcell.h"
#include "modules/plotting/datastructures/plotdata.h"
#include "modules/plotting/utils/plotlibrary/plotlibrary.h"
#include "modules/plotting/utils/plotlibrary/plotlibraryopengl.h"

#include "tgt/immediatemode/immediatemode.h"
#include "tgt/texture.h"

#include <random>

#include <Eigen/Eigenvalues>

namespace voreen {

static const int MAX_NUM_DIMENSIONS = 3;
static const int GRAB_ENABLE_DISTANCE = 0;
static const tgt::ivec2 MARGINS(75, 50);
static const float EIGENVALUE_RELATIVE_THRESHOLD = 0.01f;
static const tgt::vec3 FIRST_TIME_STEP_COLOR = tgt::vec3(1.0f, 0.0f, 0.0f);
static const tgt::vec3 LAST_TIME_STEP_COLOR = tgt::vec3::one;
static const tgt::vec3 MIN_DURATION_COLOR = tgt::vec3(1.0f, 0.0f, 0.0f);
static const tgt::vec3 MAX_DURATION_COLOR = tgt::vec3(0.0f, 0.0f, 1.0f);
static const tgt::vec3 FADE_OUT_COLOR = tgt::vec3::one;

const std::string SimilarityPlot::loggerCat_("voreen.ensembleanalysis.SimilarityPlot");

void SimilarityPlot::MDSData::serialize(Serializer& s) const {
    s.serialize("nVectors", nVectors_);
    s.serialize("eigenvalues", eigenvalues_);
}

void SimilarityPlot::MDSData::deserialize(Deserializer& s) {
    s.deserialize("nVectors", nVectors_);
    s.deserialize("eigenvalues", eigenvalues_);
}

SimilarityPlot::SimilarityPlot()
    : RenderProcessor()
    , ensembleInport_(Port::INPORT, "ensembledatastructurein", "Ensemble Datastructure Input", false)
    , seedMaskInport_(Port::INPORT, "seedmask", "Seed Mask Input (optional)")
    , outport_(Port::OUTPORT, "outport", "Outport", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , privatePort_(Port::OUTPORT, "image.tmp", "image.tmp", false)
    , pickingBuffer_(Port::OUTPORT, "picking", "Picking", false)
    , eigenValueOutport_(Port::OUTPORT, "eigenvalueOutport", "Eigenvalues", false)
    , calculateButton_("calculate", "Calculate")
    , progressBar_("progressBar", "Progress")
    , fieldSimilarityMeasure_("fieldSimilarityMeasure", "Field Similarity Measure")
    , isoValue_("isovalue", "Iso-Value", 0.5f, 0.0f, 1.0f)
    , numSeedPoints_("numSeedPoints", "Number of Seed Points", 0, 0, 0)
    , numIterations_("numIterations", "Number of Iterations", 1000, 1, 10000)
    , numEigenvalues_("numEigenvalues", "Number of Eigenvalues", MAX_NUM_DIMENSIONS, MAX_NUM_DIMENSIONS, MAX_NUM_DIMENSIONS)
    , seedTime_("seedTime", "Current Random Seed", static_cast<int>(time(0)), std::numeric_limits<int>::min(), std::numeric_limits<int>::max())
    , numDimensions_("numDimensions", "Number of Dimensions", MAX_NUM_DIMENSIONS, 1, MAX_NUM_DIMENSIONS)
    , principalComponent_("principalComponent", "Principal Component", 1, 1, MAX_NUM_DIMENSIONS)
    , scaleToMagnitude_("scaleToMagnitude", "Scale to Magnitude of Principle Component")
    , sphereRadius_("sphereRadius", "Sphere Radius", 0.01, 0.0f, 0.1f)
    , toggleAxes_("toggleAxes", "Render Axes", true, Processor::INVALID_RESULT, Property::LOD_ADVANCED)
    , fontSize_("fontSize", "Font Size", 10, 1, 30)
    , colorCoding_("colorCoding", "Color Coding")
    , renderedChannel_("renderedChannel", "Channel")
    , renderedRuns_("renderedRuns", "Rendered Runs")
    , selectedTimeSteps_("selectedTimeSteps", "Selected Time Interval", tgt::vec2(0.0f, 0.0f), 0.0f, 0.0f)
    , selectedRuns_("selectedRuns", "Selected Runs")
    , saveFileDialog_("saveFileDialog", "Export MDS Plot", "Select file...", VoreenApplication::app()->getUserDataPath(),
                      "MDS Plot data (*.mds)", FileDialogProperty::SAVE_FILE, Processor::INVALID_PATH, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
    , loadFileDialog_("loadFile", "Import MDS Plot", "Select file...", VoreenApplication::app()->getUserDataPath(),
                      "MDS Plot data (*.mds)", FileDialogProperty::OPEN_FILE, Processor::INVALID_PATH, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
    , camera_("camera", "Camera", tgt::Camera(tgt::vec3(0.0f, 0.0f, 3.5f), tgt::vec3(0.0f, 0.0f, 0.0f), tgt::vec3(0.0f, 1.0f, 0.0f)))
    , cameraHandler_(nullptr)
    , plotLib_(new PlotLibraryOpenGl())
    /*
    , grabAnchorPosition_(tgt::ivec2::zero)
    , interaction_(false)
    , clicked_(false)
    */
{
    // Ports
    ON_CHANGE(ensembleInport_, SimilarityPlot, adjustToEnsemble);
    addPort(ensembleInport_);
    addPort(seedMaskInport_);
    seedMaskInport_.addCondition(new PortConditionVolumeTypeUInt8());
    addPort(outport_);
    addPort(eigenValueOutport_);
    addPrivateRenderPort(privatePort_);
    addPrivateRenderPort(pickingBuffer_);

    // Calculation
    addProperty(calculateButton_);
        calculateButton_.setGroupID("calculation");
        calculateButton_.setReadOnlyFlag(true);
        ON_CHANGE(calculateButton_, SimilarityPlot, calculate);
    addProperty(progressBar_);
        progressBar_.setGroupID("calculation");
    addProgressBar(&progressBar_);

    addProperty(fieldSimilarityMeasure_);
        fieldSimilarityMeasure_.setGroupID("calculation");
        fieldSimilarityMeasure_.addOption("isovalue", "Iso-Surface", MEASURE_ISOSURFACE);
        fieldSimilarityMeasure_.addOption("multifield", "Multi-Field", MEASURE_MULTIFIELD);
        fieldSimilarityMeasure_.set("multifield");
        ON_CHANGE_LAMBDA(fieldSimilarityMeasure_, [this] {
            isoValue_.setVisibleFlag(fieldSimilarityMeasure_.getValue() == MEASURE_ISOSURFACE);
        });
    addProperty(isoValue_);
        isoValue_.setGroupID("calculation");
        isoValue_.setVisibleFlag(false);

    addProperty(numSeedPoints_);
        numSeedPoints_.setGroupID("calculation");
    addProperty(numIterations_);
        numIterations_.setGroupID("calculation");
    addProperty(numEigenvalues_);
        numEigenvalues_.setGroupID("calculation");
    addProperty(seedTime_);
        seedTime_.setGroupID("calculation");
    setPropertyGroupGuiName("calculation", "Calculation");

    // Rendering
    addProperty(numDimensions_);
        numDimensions_.setGroupID("rendering");
        ON_CHANGE_LAMBDA(numDimensions_, [this] {
            // We need at least as many eigenvalues as dimensions to be rendered.
            numEigenvalues_.setMinValue(numDimensions_.get());

            // Eigenvalue index only useful when rendering 1D.
            principalComponent_.setVisibleFlag(numDimensions_.get() == 1);

            // Scaling only useful when rendering 3D
            scaleToMagnitude_.setVisibleFlag(numDimensions_.get() == 3);

            // Adjust camera correctly when using 3D visualization.
            if (numDimensions_.get() == 3)
                camera_.adaptInteractionToScene(tgt::Bounds(-tgt::vec3::one, tgt::vec3::one), 0.1f, true);
        });
    addProperty(principalComponent_);
        principalComponent_.setVisibleFlag(numDimensions_.get() == 1);
        principalComponent_.setGroupID("rendering");
    addProperty(scaleToMagnitude_);
        scaleToMagnitude_.setVisibleFlag(numDimensions_.get() == 3);
        scaleToMagnitude_.setGroupID("rendering");
    addProperty(sphereRadius_);
        sphereRadius_.setGroupID("rendering");
    addProperty(fontSize_);
        fontSize_.setGroupID("rendering");
    addProperty(toggleAxes_);
        toggleAxes_.setGroupID("rendering");
    addProperty(colorCoding_);
        colorCoding_.addOption("run", "Only Run", COLOR_RUN);
        colorCoding_.addOption("timeStep", "Only Time Step", COLOR_TIMESTEP);
        colorCoding_.addOption("runAndTimeStep", "Run and Time Step", COLOR_RUN_AND_TIMESTEP);
        colorCoding_.addOption("duration", "Time Step Duration", COLOR_DURATION);
        colorCoding_.selectByValue(COLOR_RUN_AND_TIMESTEP);
        colorCoding_.setGroupID("rendering");
    addProperty(renderedChannel_);
        ON_CHANGE(renderedChannel_, SimilarityPlot, outputEigenValues);
        renderedChannel_.setGroupID("rendering");
    addProperty(renderedRuns_);
        ON_CHANGE(renderedRuns_, SimilarityPlot, renderedChannelsChanged);
        renderedRuns_.setGroupID("rendering");
    setPropertyGroupGuiName("rendering", "Rendering");

    // Selection
    addProperty(selectedTimeSteps_);
        selectedTimeSteps_.setGroupID("selection");
    addProperty(selectedRuns_);
        selectedRuns_.setGroupID("selection");
        selectedRuns_.setReadOnlyFlag(true);
    setPropertyGroupGuiName("selection", "Selection");

    // IO
    addProperty(saveFileDialog_);
        ON_CHANGE(saveFileDialog_, SimilarityPlot, save);
        saveFileDialog_.setGroupID("io");
    addProperty(loadFileDialog_);
        ON_CHANGE(loadFileDialog_, SimilarityPlot, load);
        loadFileDialog_.setGroupID("io");
    setPropertyGroupGuiName("io", "Save/Load MDS Plot");

    // Camera
    addProperty(camera_);
    cameraHandler_ = new CameraInteractionHandler("cameraHandler", "Camera", &camera_);
    cameraHandler_->setEnabled(true);
    addInteractionHandler(cameraHandler_);
}

SimilarityPlot::~SimilarityPlot() {
    delete cameraHandler_;
}

Processor* SimilarityPlot::create() const {
    return new SimilarityPlot();
}

void SimilarityPlot::initialize() {
    RenderProcessor::initialize();
    sphere_.setSphereGeometry(1.0f, tgt::vec3::zero, tgt::vec4::one, 16);
}

void SimilarityPlot::deinitialize() {
    sphere_.clear();
    RenderProcessor::deinitialize();
}

void SimilarityPlot::process() {

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
    renderingPass(true);
    pickingBuffer_.deactivateTarget();

    // Draw the runs in new parameter space.
    bool threeDimensional = numDimensions_.get() == 3;
    if(threeDimensional) {
        outport_.activateTarget();
        outport_.clearTarget();

        glLineWidth(6.0f);
        renderingPass(false);

        outport_.deactivateTarget();
    }
    else {
        privatePort_.activateTarget();
        privatePort_.clearTarget(); // TODO: may change clear color to..?

        glLineWidth(6.0f);
        renderingPass(false);

        privatePort_.deactivateTarget();

        // Draw axes.
        renderAxes();
    }

    // Restore state.
    glDepthFunc(GL_LESS);
    glLineWidth(1.0f);
    IMode.color(tgt::col4::one);
}

void SimilarityPlot::renderAxes() {

    outport_.activateTarget();
    outport_.clearTarget();

    if(toggleAxes_.get()) {
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
        if (numDimensions_.get() == 1)
            plotLib_->setDomain(Interval<plot_t>(ensembleInport_.getData()->getStartTime(),
                                                 ensembleInport_.getData()->getEndTime()), PlotLibrary::X_AXIS);
        else
            plotLib_->setDomain(Interval<plot_t>(-1.0, 1.0), PlotLibrary::X_AXIS);
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
                    if (principalComponent_.get() == 1)
                        plotLib_->renderAxisLabel(PlotLibrary::Y_AXIS, "1st PC");
                    else if (principalComponent_.get() == 2)
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

    }

    // Plot data
    float xMinMarginNDC = mapRange(MARGINS.x, 0, outport_.getSize().x, -1.0f, 1.0f);
    float xMaxMarginNDC = mapRange(outport_.getSize().x-MARGINS.x, 0, outport_.getSize().x, -1.0f, 1.0f);
    float yMinMarginNDC = mapRange(MARGINS.y, 0, outport_.getSize().y, -1.0f, 1.0f);
    float yMaxMarginNDC = mapRange(outport_.getSize().y-MARGINS.y, 0, outport_.getSize().y, -1.0f, 1.0f);

    tgt::TextureUnit::setZeroUnit();
    tgt::Texture* texture = privatePort_.getColorTexture();
    if(texture) {
        texture->enable();
        texture->bind();
    }

    IMode.begin(tgt::ImmediateMode::QUADS);
        IMode.texcoord(0.0f, 0.0f); IMode.vertex(xMinMarginNDC, yMinMarginNDC);
        IMode.texcoord(1.0f, 0.0f); IMode.vertex(xMaxMarginNDC, yMinMarginNDC);
        IMode.texcoord(1.0f, 1.0f); IMode.vertex(xMaxMarginNDC, yMaxMarginNDC);
        IMode.texcoord(0.0f, 1.0f); IMode.vertex(xMinMarginNDC, yMaxMarginNDC);
    IMode.end();

    if(texture)
        texture->disable();

    // Reset state.
    IMode.color(1.0f, 1.0f, 1.0f, 1.0f);
    outport_.deactivateTarget();
    LGL_ERROR;
}

bool SimilarityPlot::isReady() const {

    if(!ensembleInport_.isReady()) {
        setNotReadyErrorMessage("No Ensemble connected");
        return false;
    }

    if(mdsData_.empty()) {
        setNotReadyErrorMessage("No Plot data available");
        return false;
    }

    return true;
}

void SimilarityPlot::renderingPass(bool picking) {

    // Retrieve dataset.
    const EnsembleDataset* dataset = ensembleInport_.getData();

    // Retrieve selected mds data.
    const MDSData& mdsData = mdsData_[renderedChannel_.getSelectedIndex()];

    switch(numDimensions_.get()) {
    case 1:
    {
        for(int runIdx : renderingOrder_) {

            size_t runOffset = 0;
            for (int i = 0; i < runIdx; i++)
                runOffset += dataset->getRuns()[i].timeSteps_.size();

            const EnsembleDataset::Run& run = dataset->getRuns()[runIdx];
            size_t numTimeSteps = run.timeSteps_.size();
            int eigenValueIdx = principalComponent_.get() - 1;

            IMode.begin(tgt::ImmediateMode::LINE_STRIP);
            for(size_t j=0; j<numTimeSteps; j++) {
                float t = mapRange(run.timeSteps_[j].time_, dataset->getStartTime(), dataset->getEndTime(), -1.0f, 1.0f);
                IMode.color(getColor(runIdx, j, picking));
                IMode.vertex(tgt::vec2(t, mdsData.nVectors_[j+runOffset][eigenValueIdx]));
            }
            IMode.end();

            if(!picking && numTimeSteps > 0) {
                size_t selectedTimeStep = dataset->pickTimeStep(runIdx, selectedTimeSteps_.get().x);
                float x = mapRange(run.timeSteps_[selectedTimeStep].time_, dataset->getStartTime(), dataset->getEndTime(), -1.0f, 1.0f);
                tgt::vec3 position(x, mdsData.nVectors_[runOffset+selectedTimeStep][eigenValueIdx], 0.0f);
                drawTimeStepSelection(runIdx, selectedTimeStep, position);
            }
        }
        break;
    }
    case 2:
    {
        for (int runIdx : renderingOrder_) {

            size_t runOffset = 0;
            for (int i = 0; i < runIdx; i++)
                runOffset += dataset->getRuns()[i].timeSteps_.size();

            size_t numTimeSteps = dataset->getRuns()[runIdx].timeSteps_.size();

            IMode.begin(tgt::ImmediateMode::LINE_STRIP);
            for(size_t j=0; j<numTimeSteps; j++) {
                IMode.color(getColor(runIdx, j, picking));
                IMode.vertex(tgt::vec2(mdsData.nVectors_[j+runOffset][0], mdsData.nVectors_[j+runOffset][1]));
            }
            IMode.end();

            if(!picking && numTimeSteps > 0) {
                size_t selectedTimeStep = dataset->pickTimeStep(runIdx, selectedTimeSteps_.get().x);
                tgt::vec3 position(mdsData.nVectors_[runOffset+selectedTimeStep][0], mdsData.nVectors_[runOffset+selectedTimeStep][1], 0.0f);
                drawTimeStepSelection(runIdx, selectedTimeStep, position);
            }
        }
        break;
    }
    case 3:
    {
        // If using 3D visualisation, use camera interaction.
        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.loadMatrix(camera_.get().getProjectionMatrix(outport_.getSize()));
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.loadMatrix(camera_.get().getViewMatrix());

        tgt::vec3 scale = tgt::vec3::one;
        if(scaleToMagnitude_.get()) {
            // Scale each axis to it's eigenvalues size.
            scale = tgt::vec3::fromPointer(&mdsData.eigenvalues_[0]);
            scale /= tgt::vec3(mdsData.eigenvalues_[0]);
        }

        for (int runIdx : renderingOrder_) {

            size_t runOffset = 0;
            for (int i = 0; i < runIdx; i++)
                runOffset += dataset->getRuns()[i].timeSteps_.size();

            size_t numTimeSteps = dataset->getRuns()[runIdx].timeSteps_.size();

            IMode.begin(tgt::ImmediateMode::LINE_STRIP);
            for(size_t j=0; j<numTimeSteps; j++) {
                IMode.color(getColor(runIdx, j, picking));
                IMode.vertex(tgt::vec3::fromPointer(&mdsData.nVectors_[j+runOffset][0]) * scale);
            }
            IMode.end();

            if(!picking && numTimeSteps > 0) {
                size_t selectedTimeStep = dataset->pickTimeStep(runIdx, selectedTimeSteps_.get().x);
                drawTimeStepSelection(runIdx, selectedTimeStep, tgt::vec3::fromPointer(&mdsData.nVectors_[runOffset+selectedTimeStep][0])*scale);
            }
        }

        // restore matrices
        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.loadIdentity();
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.loadIdentity();
        break;
    }
    default:
        // No visualization available
        break;
    }
}

void SimilarityPlot::drawTimeStepSelection(size_t runIdx, size_t timeStepIdx, const tgt::vec3& position) const {

    // Skip rendering, if not visible anyways.
    if(sphereRadius_.get() <= std::numeric_limits<float>::epsilon())
        return;

    MatStack.pushMatrix();
    MatStack.translate(position);
    MatStack.scale(tgt::vec3(sphereRadius_.get()));

    size_t numTimeSteps = ensembleInport_.getData()->getRuns()[runIdx].timeSteps_.size();
    bool timeStepAvailable = timeStepIdx < numTimeSteps;
    if(timeStepAvailable)
        IMode.color(tgt::vec4::one);
    else {
        IMode.color(tgt::vec4(1.0f, 1.0f, 1.0f, 0.5f));
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_BLEND);
    }

    sphere_.render(GL_TRIANGLES);

    if(!timeStepAvailable) {
        glDisable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ZERO);
    }

    MatStack.popMatrix();
}

tgt::vec3 SimilarityPlot::getColor(size_t runIdx, size_t timeStepIdx, bool picking) const {

    const EnsembleDataset* dataset = ensembleInport_.getData();

    float ts = static_cast<float>(timeStepIdx) / dataset->getRuns()[runIdx].timeSteps_.size();
    if (picking)
        return tgt::vec3(static_cast<float>(runIdx) / dataset->getRuns().size(), ts, 0.0f);
    
    switch (colorCoding_.getValue()) {
    case COLOR_RUN:
        return dataset->getColor(runIdx);
    case COLOR_TIMESTEP:
        return  (1.0f - ts) * FIRST_TIME_STEP_COLOR + ts * LAST_TIME_STEP_COLOR;
    case COLOR_RUN_AND_TIMESTEP:
        return (1.0f - ts) * dataset->getColor(runIdx) + ts * FADE_OUT_COLOR;
    case COLOR_DURATION:
    {
        const Statistics& stats = dataset->getTimeStepDurationStats(runIdx);
        if(std::abs(dataset->getRuns()[runIdx].timeSteps_[timeStepIdx].duration_ - stats.getMean()) > stats.getStdDev()) {
            ts = mapRange(dataset->getRuns()[runIdx].timeSteps_[timeStepIdx].duration_, stats.getMin(), stats.getMax(), 0.0f, 1.0f);
            return (1.0f - ts) * MIN_DURATION_COLOR + ts * MAX_DURATION_COLOR;
        }
        return MIN_DURATION_COLOR;
    }
    default:
        return tgt::vec3::one;
    }
}

void SimilarityPlot::mouseClickEvent(tgt::MouseEvent* e) {

    RenderTarget* target = pickingBuffer_.getRenderTarget();

    int x = tgt::clamp(e->x(), MARGINS.x, e->viewport().x-MARGINS.x);
    int y = tgt::clamp(e->y(), MARGINS.y, e->viewport().y-MARGINS.y);

    // Not inside margins.
    if(x != e->x() || y != e->y())
        return;

    // Handle margins.
    tgt::ivec2 pixel;
    if(numDimensions_.get() < 3)
        pixel = tgt::ivec2(mapRange(tgt::vec2(x, y), tgt::vec2(MARGINS), tgt::vec2(outport_.getSize()-MARGINS), tgt::vec2::zero, tgt::vec2(target->getSize())));
    else
        pixel = mapRange(tgt::ivec2(x, y), tgt::ivec2::zero, outport_.getSize(), tgt::ivec2::zero, target->getSize());

    tgt::vec4 texel = target->getColorAtPos(tgt::ivec2(pixel.x, target->getSize().y - pixel.y - 1));
    if(texel.xyz() == tgt::vec3::zero) {
        e->accept();
        return; // No hit - preserve last selection.
    }

    // Calculate run index.
    const std::vector<EnsembleDataset::Run>& runs = ensembleInport_.getData()->getRuns();
    size_t numRuns = runs.size();
    int r = static_cast<int>(std::round(texel.r * numRuns));

    // Right-click selection changes rendering order.
    if (e->button() & tgt::MouseEvent::MOUSE_BUTTON_RIGHT) {
        renderingOrder_.erase(std::find(renderingOrder_.begin(), renderingOrder_.end(), r));
        renderingOrder_.push_front(r);
        invalidate();
    }
    // Left-click performs time step selection.
    else {
        // Calculate time step index.
        size_t numTimeSteps = runs[r].timeSteps_.size();
        size_t t = static_cast<size_t>(std::round(texel.g * numTimeSteps));

        // Update selection.
        std::vector<int> run;
        run.push_back(r);
        selectedRuns_.setSelectedRowIndices(run);

        const EnsembleDataset::TimeStep& timeStep = runs[r].timeSteps_[tgt::clamp<size_t>(t, 0, numTimeSteps - 1)];
        float lower = std::floor(timeStep.time_ * 100.0f) / 100.0f;
        float upper = std::ceil((timeStep.time_+timeStep.duration_) * 100.0f) / 100.0f;
        selectedTimeSteps_.set(tgt::vec2(lower, upper));
    }

    e->accept();
}

void SimilarityPlot::onEvent(tgt::Event* e) {
    tgt::MouseEvent* event = dynamic_cast<tgt::MouseEvent*>(e);

    //*
    if(event && event->getEventType() == tgt::MouseEvent::PRESSED && (event->modifiers() & tgt::MouseEvent::SHIFT))
        mouseClickEvent(event);
    else
        RenderProcessor::onEvent(e);
    /*/
    // TODO: does not work as expected..
    if (event) {
        tgt::ivec2 position = tgt::ivec2(event->x(), event->y());

        if(e->getEventType() == tgt::MouseEvent::PRESSED && !interaction_) {
            if(!clicked_) {
                clicked_ = true;
                grabAnchorPosition_ = position;
            }
            else if(tgt::distanceSq(grabAnchorPosition_, position) > GRAB_ENABLE_DISTANCE*GRAB_ENABLE_DISTANCE) {
                interaction_ = true;
            }
        }
        else if(e->getEventType() == tgt::MouseEvent::RELEASED) {
            clicked_ = false;
            if(!interaction_) {
                mouseClickEvent(event);
                return;
            } else
                interaction_ = false;
        }

        if(interaction_) {
            for (size_t j=0; j<cameraHandler_->getEventProperties().size() && !e->isAccepted(); ++j) {
                cameraHandler_->getEventProperties().at(j)->execute(event);
            }
            return;
        }

        event->accept();
    } else
        RenderProcessor::onEvent(e);
    //*/
}

void SimilarityPlot::adjustToEnsemble() {

    ensembleHash_.clear();
    mdsData_.clear();
    renderedChannel_.setOptions(std::deque<Option<std::string>>());
    renderedRuns_.reset();
    selectedRuns_.reset();
    calculateButton_.setReadOnlyFlag(true);

    if (!ensembleInport_.isReady())
        return;

    const EnsembleDataset* dataset = ensembleInport_.getData();

    numSeedPoints_.setMinValue(1);
    numSeedPoints_.setMaxValue(131072);
    numSeedPoints_.set(32768);

    numEigenvalues_.setMinValue(MAX_NUM_DIMENSIONS);
    numEigenvalues_.setMaxValue(static_cast<int>(dataset->getTotalNumTimeSteps()));

    renderedChannel_.blockCallbacks(true);
    for (const std::string& channel : dataset->getCommonChannels())
        renderedChannel_.addOption(channel, channel, channel);
    renderedChannel_.blockCallbacks(false);

    std::vector<int> renderedRunsIndices;
    for (const EnsembleDataset::Run& run : dataset->getRuns()) {
        renderedRuns_.addRow(run.name_, dataset->getColor(renderedRunsIndices.size()));
        selectedRuns_.addRow(run.name_, dataset->getColor(renderedRunsIndices.size()));
        renderedRunsIndices.push_back(static_cast<int>(renderedRunsIndices.size()));
    }
    renderedRuns_.setSelectedRowIndices(renderedRunsIndices);

    selectedTimeSteps_.setMinValue(dataset->getStartTime());
    selectedTimeSteps_.setMaxValue(dataset->getEndTime());
    selectedTimeSteps_.set(tgt::vec2(dataset->getStartTime(), dataset->getEndTime()));

    // Try to load plot data, if already set.
    loadFileDialog_.invalidate();

    calculateButton_.setReadOnlyFlag(false);
}

void SimilarityPlot::calculate() {

    calculateButton_.setReadOnlyFlag(true);
    ensembleHash_.clear();
    mdsData_.clear();

    setProgress(0.0f);
    for (const std::string& channel : ensembleInport_.getData()->getCommonChannels()) {

        // Calculate distance matrix.
        DMMatrix distanceMatrix = calculateDistanceMatrixFromField(channel);

        // Compute Principal components and corresponding eigenvectors.
        MDSData mdsData = computeFromDM(distanceMatrix);

        // Add the result.
        mdsData_.push_back(mdsData);
    }

    outputEigenValues();
    setProgress(1.0f);

    ensembleHash_ = EnsembleHash(*ensembleInport_.getData()).getHash();
    calculateButton_.setReadOnlyFlag(false);
    invalidate();
}

SimilarityPlot::DMMatrix SimilarityPlot::calculateDistanceMatrixFromField(const std::string& channel) {

    const float initialProgress = getProgress();
    const float progressPerChannel = 1.0f; // actually less, but good enough.

    const EnsembleDataset* dataset = ensembleInport_.getData();

    const size_t totalPoints = static_cast<size_t>(numSeedPoints_.get());

    size_t DistanceMatrixSize = dataset->getTotalNumTimeSteps();
    DMMatrix DistanceMatrix(DistanceMatrixSize);
    for(size_t i=0; i<DistanceMatrix.size(); i++)
        DistanceMatrix[i].resize(i+1);

    ///////////////////////////////////////
    //Calculate flags for all the samples//
    ///////////////////////////////////////

    DMMatrix Flags(DistanceMatrixSize);
    for (size_t i=0; i<DistanceMatrixSize; i++)
        Flags[i].resize(totalPoints);

    // Calculate random seed points.
    std::function<float()> rnd(std::bind(std::uniform_real_distribution<float>(0.0f, 1.0f), std::mt19937(seedTime_.get())));
    const tgt::Bounds& roi = dataset->getRoi();

    const VolumeBase* seedMask = seedMaskInport_.getData();
    tgt::Bounds seedMaskBounds;
    tgt::mat4 seedMaskPhysicalToVoxelMatrix;
    if(seedMask) {
        seedMaskBounds = seedMask->getBoundingBox(false).getBoundingBox();
        seedMaskPhysicalToVoxelMatrix = seedMask->getPhysicalToVoxelMatrix();
        LINFO("Restricting seed points to volume mask");
    }

    std::vector<tgt::vec3> seedPoints(totalPoints);
    for (size_t k=0; k<totalPoints; k++) {
        seedPoints[k] = tgt::vec3(rnd(), rnd(), rnd());
        seedPoints[k] = tgt::vec3(roi.getLLF()) + seedPoints[k] * tgt::vec3(roi.diagonal());

        // TODO: very rough and dirty restriction, implement something more intelligent.
        if (seedMask && (!seedMaskBounds.containsPoint(seedPoints[k]) ||
                         seedMask->getRepresentation<VolumeRAM>()->getVoxelNormalized(seedMaskPhysicalToVoxelMatrix*seedPoints[k]) == 0.0f)) {
            k--; // reseed this point.
            continue;
        }
    }

    setProgress(initialProgress + 0.01f);
    float progressIncrement = (progressPerChannel * 0.85f / DistanceMatrixSize) / dataset->getCommonChannels().size();

    size_t i = 0;
    for (const EnsembleDataset::Run& run : dataset->getRuns()) {
        for (const EnsembleDataset::TimeStep& timeStep : run.timeSteps_) {

            const tgt::vec2& valueRange = dataset->getValueRange(channel);
            const VolumeBase* volume = timeStep.channels_.at(channel);
            const VolumeRAM* data = volume->getRepresentation<VolumeRAM>();

            for (size_t k=0; k<totalPoints; k++) {
                float value = data->getVoxelNormalizedLinear(volume->getPhysicalToVoxelMatrix() * seedPoints[k]);

                switch(fieldSimilarityMeasure_.getValue()) {
                case MEASURE_ISOSURFACE:
                {
                    bool inside = value > (isoValue_.get() * (valueRange.y - valueRange.x) - valueRange.x);
                    Flags[i][k] = inside ? 1.0f : 0.0f;
                    break;
                }
                case MEASURE_MULTIFIELD:
                    Flags[i][k] = (value - valueRange.x) / (valueRange.y - valueRange.x);
                    break;
                default:
                    break;
                }
            }
            i++;

            // Update progress.
            setProgress(getProgress() + progressIncrement);
        }
    }
    setProgress(initialProgress + 0.85f / dataset->getCommonChannels().size());

    ////////////////////////////////////////////////////////////
    //Calculate distances for up-right corner and reflect them//
    ////////////////////////////////////////////////////////////

#ifdef VRN_MODULE_OPENMP
    #pragma omp parallel for shared(Flags)
#endif
    for (long i=0; i<static_cast<long>(DistanceMatrixSize); i++) {
        for (long j=0; j<i+1; j++) {
            float ScaleSum = 0.0f;
            float resValue = 0.0f;

            DistanceMatrix[i][j] = 0.0f;

            for (long k=0; k<static_cast<long>(totalPoints); k++) {

                float a = Flags[i][k];
                float b = Flags[j][k];

                ScaleSum += (1.0f - (a < b ? a : b));
                resValue += (1.0f - (a > b ? a : b));
            }

            if (ScaleSum > 0.0f)
                DistanceMatrix[i][j] = (ScaleSum-resValue) / ScaleSum;
            else
                DistanceMatrix[i][j] = 1.0f;
        }
    }
    setProgress(initialProgress + progressPerChannel / dataset->getCommonChannels().size());

    return DistanceMatrix;
}

SimilarityPlot::MDSData SimilarityPlot::computeFromDM(const DMMatrix& DistanceMatrix, float epsilon) {
    using namespace Eigen;

    const size_t dimNum = numEigenvalues_.get();
    const int iterNum = numIterations_.get();
    const size_t PointsNumber = DistanceMatrix.size();

    MDSData result;
    result.nVectors_.resize(PointsNumber);
    result.eigenvalues_.clear();

    MatrixXf Result(PointsNumber, dimNum);
    MatrixXf EigSq = MatrixXf::Zero(dimNum, dimNum);
    MatrixXf PMatrix(PointsNumber, PointsNumber);

    // Init datastructures.
    for (size_t i=0; i<PointsNumber; i++) {
        result.nVectors_[i].resize(dimNum);
        for (size_t j=0; j<=i; j++) {
            float v = DistanceMatrix[i][j];
            PMatrix(i, j) = v*v;
            PMatrix(j, i) = PMatrix(i, j);
        }
    }

    MatrixXf JMatrix = MatrixXf::Identity (PointsNumber, PointsNumber) -
              (1.0f / PointsNumber) * MatrixXf::Ones (PointsNumber, PointsNumber);

    MatrixXf BMatrix = -0.5f*JMatrix*PMatrix*JMatrix;

    VectorXf TempVector(PointsNumber);
    VectorXf* EVectors = new VectorXf[dimNum];

    for(size_t i=0; i<dimNum; i++) {

        VectorXf& EVector = EVectors[i];
        EVector = VectorXf::Ones(PointsNumber);
        float EValue = 0.0;

        VectorXf PrVector = VectorXf::Zero(PointsNumber);
        VectorXf TVector = VectorXf::Ones(PointsNumber);

        int iter = 0;

        MatrixXf EMVector(PointsNumber, 1);

        while (iter < iterNum && TVector.norm() > epsilon) {

            EMVector.col(0) = EVector;
            TempVector = BMatrix * EMVector;
            EVector = TempVector;
            for(size_t j=0; j < i; j++)
                EVector -= EVectors[j] * (EVectors[j].dot(TempVector));
            EValue = EVector.norm();
            EVector.normalize();
            TVector = EVector - PrVector;
            PrVector = EVector;
            iter++;
        }

        EMVector.col(0) = EVector;

        // Don't continue calculating when eigenvalues get too small in relation to biggest eigenvalue.
        if(i >= static_cast<size_t>(numDimensions_.get()) && EValue < result.eigenvalues_[0] * EIGENVALUE_RELATIVE_THRESHOLD)
            break;

        Result.col(i) = EVector;
        EigSq(i, i) = std::sqrt(EValue);
        result.eigenvalues_.push_back(EValue);
    }
    delete [] EVectors;

    //Now get the resulting matrix
    Result = Result*EigSq;

    for (size_t i=0; i<result.eigenvalues_.size(); i++) {

        float MaxValue = Result(0, i);
        float MinValue = Result(0, i);

        for (size_t j=1; j<PointsNumber; j++) {
            if (MaxValue < Result(j, i))
                MaxValue = Result(j, i);
            else if (MinValue > Result(j, i))
                MinValue = Result(j, i);
        }

        for (size_t j = 0; j < PointsNumber; j++) {
            float value = (Result(j, i) - MinValue) / (MaxValue - MinValue) * 2.0f - 1.0f;
            result.nVectors_[j][i] = value;
        }
    }

    return result;
}

void SimilarityPlot::outputEigenValues() {

    if (mdsData_.empty()) {
        eigenValueOutport_.setData(nullptr);
        return;
    }

    PlotData* data = new PlotData(1, 1);

    data->setColumnLabel(0, "Index");
    data->setColumnLabel(1, "Eigenvalue");

    const MDSData& mdsData = mdsData_[renderedChannel_.getSelectedIndex()];
    for(size_t i = 0; i < mdsData.eigenvalues_.size(); i++) {
        std::vector<PlotCellValue> values;
        values.push_back(PlotCellValue(i+1));
        values.push_back(PlotCellValue(static_cast<plot_t>(mdsData.eigenvalues_[i])));
        data->insert(values);
    }
    eigenValueOutport_.setData(data, true);
}

void SimilarityPlot::save() {
    if (saveFileDialog_.get().empty()) {
        LWARNING("no filename specified");
        return;
    }
    if(mdsData_.empty()) {
        LWARNING("No vectors calculated.");
        return;
    }

    // write file
    try {
        std::ofstream outFile;
        outFile.open(saveFileDialog_.get().c_str());
        LINFO("Writing mds to file " << saveFileDialog_.get());

        XmlSerializer s;
        Serializer serializer(s);
        serializer.serialize("hash", ensembleHash_);
        serializer.serialize("mds_data", mdsData_);
        s.write(outFile);

        outFile.close();
        LINFO("Saving " << saveFileDialog_.get() << " was successful.");
    }
    catch(tgt::Exception& e) {
        VoreenApplication::app()->showMessageBox("saving MDS Plot failed", e.what(), true);
        LERROR(e.what());
        saveFileDialog_.set("");
    }
}

void SimilarityPlot::renderedChannelsChanged() {
    renderingOrder_ = std::deque<int>(renderedRuns_.getSelectedRowIndices().begin(), renderedRuns_.getSelectedRowIndices().end());
}

void SimilarityPlot::load() {
    if (loadFileDialog_.get().empty()) {
        LWARNING("no filename specified");
        return;
    }
    if (!ensembleInport_.isReady()) {
        LWARNING("no ensemble connected");
        return;
    }

    std::string absPath = loadFileDialog_.get();
    try {
        std::ifstream inFile;
        inFile.open(absPath.c_str());

        XmlDeserializer d;
        d.read(inFile);
        Deserializer deserializer(d);
        deserializer.deserialize("hash", ensembleHash_);
        if (ensembleHash_ != EnsembleHash(*ensembleInport_.getData()).getHash())
            throw VoreenException("The plot was probably not generated by the currently loaded ensemble dataset");
        deserializer.deserialize("mds_data", mdsData_);

        inFile.close();
        outputEigenValues();
        LINFO("Reading mds plot vectors from file " << absPath << " was successful");
    }
    catch(tgt::Exception& e) {
        VoreenApplication::app()->showMessageBox("loading MDS Plot failed", e.what(), true);
        LERROR(e.what());
        loadFileDialog_.set("");
    }
}

} // namespace
