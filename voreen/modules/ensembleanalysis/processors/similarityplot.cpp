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

#include "similarityplot.h"

#include "../utils/ensemblehash.h"
#include "../utils/utils.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/datastructures/volume/volume.h"
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
static const tgt::vec3 FIRST_TIME_STEP_COLOR(1.0f, 0.0f, 0.0f);
static const tgt::vec3 LAST_TIME_STEP_COLOR = tgt::vec3::one;
static const tgt::vec3 MIN_DURATION_COLOR(1.0f, 0.0f, 0.0f);
static const tgt::vec3 MAX_DURATION_COLOR(0.0f, 0.0f, 1.0f);
static const tgt::vec3 FADE_OUT_COLOR = tgt::vec3::one;
static const float SELECTED_LINE_WIDTH = 6.0f;
static const float UNSELECTED_LINE_WIDTH = 3.0f;
static const float PICKING_LINE_WIDTH = SELECTED_LINE_WIDTH;

const std::string SimilarityPlot::fontName_("Vera.ttf");
const std::string SimilarityPlot::loggerCat_("voreen.ensembleanalysis.SimilarityPlot");

void SimilarityPlot::Embedding::serialize(Serializer& s) const {
    s.serialize("nVectors", nVectors_);
    s.serialize("eigenvalues", eigenvalues_);
    if(!names_.empty())
        s.serialize("names", names_);
}

void SimilarityPlot::Embedding::deserialize(Deserializer& s) {
    s.deserialize("nVectors", nVectors_);
    s.deserialize("eigenvalues", eigenvalues_);
    s.optionalDeserialize("names", names_, decltype(names_)());
}

SimilarityPlot::SimilarityPlot()
    : RenderProcessor()
    , ensembleInport_(Port::INPORT, "ensembledatastructurein", "Ensemble Datastructure Input", false)
    , similarityMatrixInport_(Port::INPORT, "similaritymatrixin", "Similarity Matrix Input", false)
    , outport_(Port::OUTPORT, "outport", "Outport", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , privatePort_(Port::OUTPORT, "image.tmp", "image.tmp", false)
    , pickingBuffer_(Port::OUTPORT, "picking", "Picking", false)
    , eigenValueOutport_(Port::OUTPORT, "eigenvalueOutport", "Eigenvalues", false)
    , calculateButton_("calculate", "Create Embeddings")
    , autoCalculate_("autoCalculate", "Auto Calculate", true)
    , progressBar_("progressBar", "Progress")
    , numIterations_("numIterations", "Number of Iterations", 1000, 1, 10000)
    , numEigenvalues_("numEigenvalues", "Number of Eigenvalues", MAX_NUM_DIMENSIONS, MAX_NUM_DIMENSIONS, MAX_NUM_DIMENSIONS)
    , numDimensions_("numDimensions", "Number of Dimensions", MAX_NUM_DIMENSIONS, 1, MAX_NUM_DIMENSIONS)
    , principleComponent_("principalComponent", "Principle Component", 1, 1, MAX_NUM_DIMENSIONS)
    , scaleToMagnitude_("scaleToMagnitude", "Scale to Magnitude of Principle Component")
    , sphereRadius_("sphereRadius", "Sphere Radius", 0.01, 0.0f, 0.1f)
    , fontSize_("fontSize", "Font Size", 10, 1, 30)
    , showTooltip_("showTooltip", "Show Tooltip", true)
    , renderTimeSelection_("renderTimeSelection", "Render Time Selection", false, Processor::INVALID_RESULT, Property::LOD_ADVANCED)
    , colorCoding_("colorCoding", "Color Coding")
    , renderedField_("renderedChannel", "Field")
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
    , camera_("camera", "Camera", tgt::Camera(tgt::vec3(0.0f, 0.0f, 3.5f), tgt::vec3(0.0f, 0.0f, 0.0f), tgt::vec3(0.0f, 1.0f, 0.0f)))
    , cameraHandler_(nullptr)
    , plotLib_(new PlotLibraryOpenGl())
    , lastHit_(boost::none)
{
    // Ports
    //NOTE: don't use adjustPropertiesToInput() as callback, this seems to be triggered each frame for RenderProcessors.
    addPort(ensembleInport_);
    ON_CHANGE(ensembleInport_, SimilarityPlot, adjustToEnsemble);
    addPort(similarityMatrixInport_);
    ON_CHANGE(similarityMatrixInport_, SimilarityPlot, adjustToEnsemble);
    addPort(outport_);
    addPort(eigenValueOutport_);
    addPrivateRenderPort(privatePort_);
    addPrivateRenderPort(pickingBuffer_);

    // Calculation
    addProperty(calculateButton_);
        calculateButton_.setGroupID("calculation");
        calculateButton_.setReadOnlyFlag(true);
        ON_CHANGE(calculateButton_, SimilarityPlot, createEmbeddings);
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
            scaleToMagnitude_.setVisibleFlag(numDimensions_.get() == 3);

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
    addProperty(scaleToMagnitude_);
        scaleToMagnitude_.setVisibleFlag(numDimensions_.get() == 3);
        scaleToMagnitude_.setGroupID("rendering");
    addProperty(sphereRadius_);
        sphereRadius_.setGroupID("rendering");
    addProperty(fontSize_);
        fontSize_.setGroupID("rendering");
    addProperty(showTooltip_);
        showTooltip_.setGroupID("rendering");
    addProperty(renderTimeSelection_);
        renderTimeSelection_.setGroupID("rendering");
    addProperty(colorCoding_);
        colorCoding_.addOption("member", "Only Member", COLOR_MEMBER);
        colorCoding_.addOption("timeStep", "Only Time Step", COLOR_TIMESTEP);
        colorCoding_.addOption("memberAndTimeStep", "Member and Time Step", COLOR_MEMBER_AND_TIMESTEP);
        colorCoding_.addOption("duration", "Time Step Duration", COLOR_DURATION);
        colorCoding_.selectByValue(COLOR_MEMBER_AND_TIMESTEP);
        colorCoding_.setGroupID("rendering");
    addProperty(renderedField_);
        ON_CHANGE(renderedField_, SimilarityPlot, outputEigenValues);
        renderedField_.setGroupID("rendering");
    addProperty(renderedMembers_);
        ON_CHANGE(renderedMembers_, SimilarityPlot, renderedMembersChanged);
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
        ON_CHANGE(saveFileDialog_, SimilarityPlot, saveEmbeddings);
        saveFileDialog_.setGroupID("io");
    addProperty(saveButton_);
        ON_CHANGE(saveButton_, SimilarityPlot, saveEmbeddings);
        saveButton_.setGroupID("io");
    addProperty(loadFileDialog_);
        ON_CHANGE(loadFileDialog_, SimilarityPlot, loadEmbeddings);
        loadFileDialog_.setGroupID("io");
    addProperty(loadButton_);
        ON_CHANGE(loadButton_, SimilarityPlot, loadEmbeddings);
        loadButton_.setGroupID("io");
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
    glLineWidth(PICKING_LINE_WIDTH);
    renderingPass(true);
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

void SimilarityPlot::renderAxes() {

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
        plotLib_->setDomain(Interval<plot_t>(ensembleInport_.getData()->getStartTime(),
                                             ensembleInport_.getData()->getEndTime()), PlotLibrary::X_AXIS);
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

void SimilarityPlot::renderTooltip() const {

    const EnsembleDataset* ensemble = ensembleInport_.getData();
    if(!ensemble || !lastHit_) {
        return;
    }

    glDisable(GL_DEPTH_TEST);

    tgt::Font font(VoreenApplication::app()->getFontPath(fontName_));
    font.setFontSize(fontSize_.get());

    const int legendOffsetX = 10.f;
    const int legendOffsetY = 10.f;

    const tgt::ivec2& screensize = outport_.getSize();

    const EnsembleMember& member = ensemble->getMembers()[lastHit_->memberIdx];
    const TimeStep& timeStep = member.getTimeSteps()[lastHit_->timeStepIdx];

    std::string tooltip = "Member: ";
    tooltip += member.getName();
    tooltip += "\n";
    tooltip += "Time:";
    tooltip += std::to_string(timeStep.getTime());

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

bool SimilarityPlot::isReady() const {

    if(!ensembleInport_.isReady()) {
        setNotReadyErrorMessage("No Ensemble connected");
        return false;
    }

    if(ensembleInport_.getData()->getCommonFieldNames().empty()) {
        setNotReadyErrorMessage("No common fields available");
        return false;
    }

    // Note: Similarity Matrix is optional.

    if(embeddings_.empty()) {
        setNotReadyErrorMessage("No embedding available");
        return false;
    }

    return true;
}

void SimilarityPlot::renderingPass(bool picking) {

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

void SimilarityPlot::renderEmbedding1D(bool picking) {

    // Retrieve dataset.
    const EnsembleDataset* dataset = ensembleInport_.getData();

    // Retrieve selected embedding.
    const Embedding& embedding = embeddings_[renderedField_.getSelectedIndex()];

    tgt::vec2 timeRange = tgt::vec2(dataset->getStartTime(), dataset->getEndTime());
    if(renderTimeSelection_.get()) {
        timeRange = firstSelectedTimeInterval_.get();
    }

    for(int memberIdx : renderingOrder_) {

        glLineWidth((subSelection_.count(memberIdx) != 0) ? SELECTED_LINE_WIDTH : UNSELECTED_LINE_WIDTH);

        const EnsembleMember& member = dataset->getMembers()[memberIdx];
        size_t numTimeSteps = member.getTimeSteps().size();
        int eigenValueIdx = principleComponent_.get() - 1;
        const auto& vertices = embedding.nVectors_.at(memberIdx);

        // In case we have a single time step, we draw it across the whole range,
        // since it doesn't change. This could (and should!) be improved, however,
        // such that it becomes clear at which t the time step is recorded.
        if(numTimeSteps == 1) {
            IMode.begin(tgt::ImmediateMode::FAKE_LINES);
            IMode.color(getColor(memberIdx, 0, picking));
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
                if(!picking && (member.getTimeSteps()[j].getTime() < timeRange.x || member.getTimeSteps()[j].getTime() > timeRange.y)) {
                    colorSaturation = 0.25f;
                }
                float t = mapRange(member.getTimeSteps()[j].getTime(), dataset->getStartTime(), dataset->getEndTime(), -1.0f, 1.0f);
                IMode.color(getColor(memberIdx, j, picking) * colorSaturation + tgt::vec3(1.0f - colorSaturation));
                IMode.vertex(tgt::vec2(t, vertices[j][eigenValueIdx]));
            }
            IMode.end();
        }
    }

    if(!picking && renderTimeSelection_.get()) {
        tgt::vec2 mappedTimeRange = mapRange(firstSelectedTimeInterval_.get(), tgt::vec2(dataset->getStartTime()), tgt::vec2(dataset->getEndTime()), -tgt::vec2::one, tgt::vec2::one);

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

void SimilarityPlot::renderEmbedding2D(bool picking) {

    // Retrieve dataset.
    const EnsembleDataset* dataset = ensembleInport_.getData();

    // Retrieve selected embedding.
    const Embedding& embedding = embeddings_[renderedField_.getSelectedIndex()];

    for (int memberIdx : renderingOrder_) {

        bool selected = (subSelection_.count(memberIdx) != 0);
        glLineWidth(selected ? SELECTED_LINE_WIDTH : UNSELECTED_LINE_WIDTH);

        size_t numTimeSteps = dataset->getMembers()[memberIdx].getTimeSteps().size();
        const auto& vertices = embedding.nVectors_.at(memberIdx);

        IMode.begin(tgt::ImmediateMode::FAKE_LINE_STRIP);
        for(size_t j=0; j<numTimeSteps; j++) {
            IMode.color(getColor(memberIdx, j, picking));
            IMode.vertex(tgt::vec2(vertices[j][0], vertices[j][1]));
        }
        IMode.end();

        if((!picking && renderTimeSelection_.get()) || numTimeSteps == 1) {
            size_t selectedTimeStep = dataset->getMembers()[memberIdx].getTimeStep(firstSelectedTimeInterval_.get().x);
            tgt::vec3 position(vertices[selectedTimeStep][0], vertices[selectedTimeStep][1], 0.0f);
            tgt::vec3 color = (numTimeSteps == 1) ? getColor(memberIdx, selectedTimeStep, picking) : tgt::vec3::one;
            renderSphere(position, color, selected);
        }
    }
}

void SimilarityPlot::renderEmbedding3D(bool picking) {

    // Retrieve dataset.
    const EnsembleDataset* dataset = ensembleInport_.getData();

    // Retrieve selected embedding.
    const Embedding& embedding = embeddings_[renderedField_.getSelectedIndex()];

    // If using 3D visualisation, use camera interaction.
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadMatrix(camera_.get().getProjectionMatrix(outport_.getSize()));
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadMatrix(camera_.get().getViewMatrix());

    tgt::vec3 scale = tgt::vec3::one;
    if(scaleToMagnitude_.get()) {
        // Scale each axis to it's eigenvalues size.
        scale = tgt::vec3::fromPointer(&embedding.eigenvalues_[0]);
        scale /= tgt::vec3(embedding.eigenvalues_[0]);
    }

    for (int memberIdx : renderingOrder_) {

        bool selected = (subSelection_.count(memberIdx) != 0);
        glLineWidth(selected ? SELECTED_LINE_WIDTH : UNSELECTED_LINE_WIDTH);

        size_t numTimeSteps = dataset->getMembers()[memberIdx].getTimeSteps().size();
        const auto& vertices = embedding.nVectors_.at(memberIdx);

        IMode.begin(tgt::ImmediateMode::FAKE_LINE_STRIP);
        for(size_t j=0; j<numTimeSteps; j++) {
            IMode.color(getColor(memberIdx, j, picking));
            IMode.vertex(tgt::vec3::fromPointer(&vertices[j][0]) * scale);
        }
        IMode.end();

        if((!picking && renderTimeSelection_.get()) || numTimeSteps == 1) {
            size_t selectedTimeStep = dataset->getMembers()[memberIdx].getTimeStep(firstSelectedTimeInterval_.get().x);
            tgt::vec3 position = tgt::vec3::fromPointer(&vertices[selectedTimeStep][0])*scale;
            tgt::vec3 color = (numTimeSteps == 1) ? getColor(memberIdx, selectedTimeStep, picking) : tgt::vec3::one;
            renderSphere(position, color, selected);
        }
    }

    // restore matrices
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();
}

void SimilarityPlot::renderSphere(const tgt::vec3& position, const tgt::vec3& color, bool drawBorder) const {

    // Skip rendering, if not visible anyways.
    if(sphereRadius_.get() <= std::numeric_limits<float>::epsilon())
        return;

    // TODO: Draw using fragment shader!

    MatStack.pushMatrix();
    MatStack.translate(position);
    MatStack.scale(tgt::vec3(sphereRadius_.get()));
    /*
    // TODO: currently not working.
    if(drawBorder) {
        IMode.color(tgt::vec4(0.0f, 0.0f, 0.0f, 1.0f));
        sphere_.render(GL_TRIANGLES);
        MatStack.scale(tgt::vec3(0.9f));
    }
    */
    IMode.color(tgt::vec4(color, 1.0f));
    sphere_.render(GL_TRIANGLES);
    MatStack.popMatrix();
}

tgt::vec3 SimilarityPlot::getColor(size_t memberIdx, size_t timeStepIdx, bool picking) const {

    const EnsembleDataset* dataset = ensembleInport_.getData();
    const EnsembleMember& member = dataset->getMembers()[memberIdx];

    float ts = static_cast<float>(timeStepIdx) / member.getTimeSteps().size();
    if (picking)
        return tgt::vec3(static_cast<float>(memberIdx) / dataset->getMembers().size(), ts, 1.0f);

    switch (colorCoding_.getValue()) {
    case COLOR_MEMBER:
        return member.getColor();
    case COLOR_TIMESTEP:
        return  (1.0f - ts) * FIRST_TIME_STEP_COLOR + ts * LAST_TIME_STEP_COLOR;
    case COLOR_MEMBER_AND_TIMESTEP:
        return (1.0f - ts) * member.getColor() + ts * FADE_OUT_COLOR;
    case COLOR_DURATION:
    {
        const Statistics& stats = member.getTimeStepDurationStats();
        if(std::abs(dataset->getMembers()[memberIdx].getTimeSteps()[timeStepIdx].getDuration() - stats.getMean()) > stats.getStdDev()) {
            ts = mapRange(dataset->getMembers()[memberIdx].getTimeSteps()[timeStepIdx].getDuration(), stats.getMin(), stats.getMax(), 0.0f, 1.0f);
            return (1.0f - ts) * MIN_DURATION_COLOR + ts * MAX_DURATION_COLOR;
        }
        return MIN_DURATION_COLOR;
    }
    default:
        return tgt::vec3::one;
    }
}

void SimilarityPlot::mouseEvent(tgt::MouseEvent* e) {

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

        const EnsembleDataset* dataset = ensembleInport_.getData();
        if(y != e->y() && numDimensions_.get() == 1) {
            float t0 = firstSelectedTimeInterval_.get().x;
            float t1 = firstSelectedTimeInterval_.get().y;

            // Select a time interval outside of the axis
            if(e->action() == tgt::MouseEvent::PRESSED) {
                t0 = mapRange(x, MARGINS.x, e->viewport().x - MARGINS.x, dataset->getStartTime(), dataset->getEndTime());
            }
            else if(e->action() == tgt::MouseEvent::RELEASED) {
                t1 = mapRange(x, MARGINS.x, e->viewport().x - MARGINS.x, dataset->getStartTime(), dataset->getEndTime());
            } // We don't support MOTION atm., since it would invalidate too often and would trigger linked filters.
            else {
                return;
            }

            // Fix range.
            if(t0 > t1) {
                std::swap(t0, t1);
            }

            // Add small offset, if necessary to include at least a single time step in the selected range.
            if(t0 < t1 + dataset->getMinTimeStepDuration()) {
                t1 += dataset->getMinTimeStepDuration();
            }

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
        const std::vector<EnsembleMember>& members = dataset->getMembers();
        size_t numMembers = members.size();
        int r = tgt::clamp<int>(std::round(texel.r * numMembers), 0, numMembers - 1);

        // Calculate time step index.
        size_t numTimeSteps = members[r].getTimeSteps().size();
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

            const TimeStep& timeStep = members[r].getTimeSteps()[t];
            float lower = std::floor(timeStep.getTime() * 100.0f) / 100.0f;
            float upper = std::ceil((timeStep.getTime() + timeStep.getDuration()) * 100.0f) / 100.0f;
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

void SimilarityPlot::onEvent(tgt::Event* e) {

    // Events will be triggered for RenderProcessors, even if the processor is not ready it seems.
    if (isReady()) {
        tgt::MouseEvent* event = dynamic_cast<tgt::MouseEvent*>(e);
        if(event) {
            mouseEvent(event);
        }
    }

    RenderProcessor::onEvent(e);
}

void SimilarityPlot::adjustToEnsemble() {

    ensembleHash_.clear();
    embeddings_.clear();
    subSelection_.clear();
    renderedField_.setOptions(std::deque<Option<std::string>>());
    renderedMembers_.reset();
    firstSelectedMember_.reset();
    secondSelectedMember_.reset();
    calculateButton_.setReadOnlyFlag(true);

    // Check for data set.
    const EnsembleDataset* dataset = ensembleInport_.getData();
    if (!dataset)
        return;

    // Check for similarity matrices.
    const SimilarityMatrixList* similarityMatrices = similarityMatrixInport_.getData();
    if (!similarityMatrices)
        return;

    // Check if both match.
    if(EnsembleHash(*dataset).getHash() != similarityMatrices->getHash())
        return;

    numEigenvalues_.setMinValue(MAX_NUM_DIMENSIONS);
    numEigenvalues_.setMaxValue(static_cast<int>(dataset->getTotalNumTimeSteps()));

    renderedField_.blockCallbacks(true);
    for (const std::string& fieldName : dataset->getCommonFieldNames())
        renderedField_.addOption(fieldName, fieldName, fieldName);
    renderedField_.blockCallbacks(false);

    std::vector<int> memberIndices;
    for (const EnsembleMember& member : dataset->getMembers()) {
        renderedMembers_.addRow(member.getName(), member.getColor());
        firstSelectedMember_.addRow(member.getName(), member.getColor());
        secondSelectedMember_.addRow(member.getName(), member.getColor());
        subSelection_.insert(static_cast<int>(memberIndices.size()));
        memberIndices.push_back(static_cast<int>(memberIndices.size()));
    }
    renderedMembers_.setSelectedRowIndices(memberIndices);
    firstSelectedMember_.setSelectedRowIndices(memberIndices);
    secondSelectedMember_.setSelectedRowIndices(memberIndices);

    firstSelectedTimeInterval_.setMinValue(dataset->getStartTime());
    firstSelectedTimeInterval_.setMaxValue(dataset->getEndTime());
    firstSelectedTimeInterval_.set(tgt::vec2(dataset->getStartTime(), dataset->getEndTime()));

    secondSelectedTimeInterval_.setMinValue(dataset->getStartTime());
    secondSelectedTimeInterval_.setMaxValue(dataset->getEndTime());
    secondSelectedTimeInterval_.set(tgt::vec2(dataset->getStartTime(), dataset->getEndTime()));

    // Try to load plot data, if already set.
    if(!loadFileDialog_.get().empty()) {
        loadEmbeddings();
    }

    calculateButton_.setReadOnlyFlag(false);

    if (embeddings_.empty() && autoCalculate_.get()) {
        createEmbeddings();
    }
}

void SimilarityPlot::createEmbeddings() {

    if(subSelection_.empty()) {
        LERROR("No member selected");
        return;
    }

    size_t numMembers = ensembleInport_.getData()->getMembers().size();
    if(subSelection_.size() == numMembers) {
        LINFO("Calculating for whole data set..");
    }
    else {
        LINFO("Calculating for subset..");
    }

    std::vector<std::string> names;
    for(const auto& member : ensembleInport_.getData()->getMembers()) {
        names.push_back(member.getName());
    }

    const SimilarityMatrixList* matrices = similarityMatrixInport_.getData();
    std::vector<std::string> fieldNames = matrices->getFieldNames();

    calculateButton_.setReadOnlyFlag(true);
    ensembleHash_.clear();
    embeddings_.resize(fieldNames.size());

    setProgress(0.0f);
    ThreadedTaskProgressReporter progressReporter(*this, fieldNames.size());

#ifdef VRN_MODULE_OPENMP
    #pragma omp parallel for
#endif
    for (long i=0; i<static_cast<long>(fieldNames.size()); i++) {

        // Get distance matrix for current field.
        const SimilarityMatrix& distanceMatrix = matrices->getSimilarityMatrix(fieldNames[i]);

        // Compute Principal components and corresponding eigenvectors.
        Embedding embedding = createEmbedding(distanceMatrix);

        // Add member name.
        embedding.names_ = names;

        // Add the result.
        embeddings_[i] = std::move(embedding);

        progressReporter.reportStepDone();
    }

    // Finally, output eigen values to allow an eigen value analysis.
    outputEigenValues();

    // Update rendering order.
    renderedMembers_.setSelectedRowIndices(std::vector<int>(subSelection_.rbegin(), subSelection_.rend()));
    //renderingOrder_.assign(subSelection_.rbegin(), subSelection_.rend());

    // Update hash.
    ensembleHash_ = EnsembleHash(*ensembleInport_.getData()).getHash();

    // Done.
    calculateButton_.setReadOnlyFlag(false);
    setProgress(1.0f);
    invalidate();
}

SimilarityPlot::Embedding SimilarityPlot::createEmbedding(const SimilarityMatrix& distanceMatrix, float epsilon) const {
    using namespace Eigen;

    const std::vector<EnsembleMember>& members = ensembleInport_.getData()->getMembers();

    const size_t numDimensions = numEigenvalues_.get();
    const int numIterations = numIterations_.get();
    const size_t numMembers = members.size();
    const std::set<int> selection(subSelection_); // Copy selection for thread safety.

    Embedding embedding;

    size_t numPoints = 0;
    for(int memberIdx : selection) {
        size_t numTimeSteps = members[memberIdx].getTimeSteps().size();
        numPoints += numTimeSteps;
        embedding.nVectors_[memberIdx] = std::vector<std::vector<float>>(numTimeSteps, std::vector<float>(numDimensions, 0.0f));
    }

    MatrixXf result(numPoints, numDimensions);
    MatrixXf EigSq = MatrixXf::Zero(numDimensions, numDimensions);
    MatrixXf PMatrix(numPoints, numPoints);

    // Init datastructures.
    size_t offsetA = 0;
    size_t positionA = 0;

    // Iterate each member, call it A.
    for(size_t memberIdxA=0; memberIdxA < numMembers; memberIdxA++) {
        size_t offsetB = 0;
        size_t positionB = 0;
        size_t numTimeStepsA = members[memberIdxA].getTimeSteps().size();
        if(selection.count(memberIdxA) != 0) { // Do we consider member A?
            // Again iterate each member, call it B. Now looking at pairs of member A and B.
            for(size_t memberIdxB=0; memberIdxB<=memberIdxA; memberIdxB++) {
                size_t numTimeStepsB = members[memberIdxB].getTimeSteps().size();
                if (selection.count(memberIdxB) != 0) { // Do we consider member B?
                    // Iterate time steps of member A.
                    for (size_t i = 0; i < numTimeStepsA; i++) {
                        // Iterate time steps of member B.
                        for (size_t j = 0; j < numTimeStepsB; j++) {
                            float v = distanceMatrix(i + offsetA, j + offsetB);
                            PMatrix(i + positionA, j + positionB) = PMatrix(j + positionB, i + positionA) = v * v;
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

    MatrixXf JMatrix = MatrixXf::Identity (numPoints, numPoints) -
                       (1.0f / numPoints) * MatrixXf::Ones (numPoints, numPoints);

    MatrixXf BMatrix = -0.5f*JMatrix*PMatrix*JMatrix;

    VectorXf* eigenVectors = new VectorXf[numDimensions];

    for(size_t i=0; i < numDimensions; i++) {

        VectorXf& eigenVector = eigenVectors[i];
        eigenVector = VectorXf::Ones(numPoints);
        float eigenValue = 0.0;

        VectorXf PrVector = VectorXf::Zero(numPoints);
        VectorXf TVector = VectorXf::Ones(numPoints);

        MatrixXf EMVector(numPoints, 1);

        for (int iter=0; iter < numIterations; iter++) {

            EMVector.col(0) = eigenVector;
            VectorXf TempVector = BMatrix * EMVector;
            eigenVector = TempVector;
            for(size_t j=0; j < i; j++)
                eigenVector -= eigenVectors[j] * (eigenVectors[j].dot(TempVector));
            eigenValue = eigenVector.norm();
            eigenVector.normalize();
            TVector = eigenVector - PrVector;
            PrVector = eigenVector;

            if(TVector.norm() <= epsilon) {
                break;
            }
        }

        // Don't continue calculating when eigenvalues get too small in relation to biggest eigenvalue.
        if(i >= static_cast<size_t>(numDimensions_.get()) && // We do have calculated at least as many PCs as we want to display.
        (eigenValue <= std::numeric_limits<float>::epsilon() || // Eigenvalue is de-facto zero.
        (i > 0 && eigenValue > embedding.eigenvalues_[i - 1]))) // EV must be smaller than predecessor.
            break;

        EMVector.col(0) = eigenVector;
        result.col(i) = eigenVector;
        EigSq(i, i) = std::sqrt(eigenValue);
        embedding.eigenvalues_.push_back(eigenValue);
    }
    delete [] eigenVectors;

    // Now get the resulting matrix.
    result = result * EigSq;

    for (size_t i=0; i < embedding.eigenvalues_.size(); i++) {

        float maxValue = result(0, i);
        float minValue = result(0, i);

        for (size_t j=1; j < numPoints; j++) {
            if (maxValue < result(j, i))
                maxValue = result(j, i);
            else if (minValue > result(j, i))
                minValue = result(j, i);
        }

        // Interpret as member and time steps.
        size_t j=0;
        for (int memberIdx : selection) {
            for(size_t t=0; t<members[memberIdx].getTimeSteps().size(); t++) {
                float value = mapRange(result(j, i), minValue, maxValue, -1.0f, 1.0f);
                embedding.nVectors_[memberIdx][t][i] = value;
                j++;
            }
        }
    }

    return embedding;
}

void SimilarityPlot::outputEigenValues() {

    if (embeddings_.empty()) {
        eigenValueOutport_.setData(nullptr);
        return;
    }

    PlotData* data = new PlotData(1, 1);

    data->setColumnLabel(0, "Index");
    data->setColumnLabel(1, "Eigenvalue");

    const Embedding& embedding = embeddings_[renderedField_.getSelectedIndex()];
    for(size_t i = 0; i < embedding.eigenvalues_.size(); i++) {
        std::vector<PlotCellValue> values;
        values.push_back(PlotCellValue(i+1));
        values.push_back(PlotCellValue(static_cast<plot_t>(embedding.eigenvalues_[i])));
        data->insert(values);
    }
    eigenValueOutport_.setData(data, true);
}

void SimilarityPlot::renderedMembersChanged() {
    renderingOrder_.clear();

    if(embeddings_.empty())
        return;

    for(int memberIdx : renderedMembers_.getSelectedRowIndices()) {
        // The first field is representative for all fields here.
        // We just want to check if the embedding was already calculated for the member.
        if(embeddings_.front().nVectors_.count(memberIdx)) {
            renderingOrder_.push_back(memberIdx);
        }
    }
}

void SimilarityPlot::saveEmbeddings() {
    if (saveFileDialog_.get().empty()) {
        LWARNING("no filename specified");
        return;
    }
    if(embeddings_.empty()) {
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
        serializer.serialize("hash", ensembleHash_);
        serializer.serialize("embeddings", embeddings_);
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

void SimilarityPlot::loadEmbeddings() {
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

        JsonDeserializer d;
        d.read(inFile);
        Deserializer deserializer(d);
        deserializer.deserialize("hash", ensembleHash_);
        if (ensembleHash_ != EnsembleHash(*ensembleInport_.getData()).getHash())
            throw VoreenException("The plot was probably not generated by the currently loaded ensemble dataset");
        deserializer.deserialize("embeddings", embeddings_);

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

} // namespace
