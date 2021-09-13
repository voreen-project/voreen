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

#include "distancetooriginplot.h"

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

#include "../../../modules/ensembleanalysis/utils/utils.h"


namespace voreen {

static const int MAX_NUM_DIMENSIONS = 3;
static const tgt::ivec2 MIN_MARGINS(75, 50);
static const tgt::vec3 FIRST_TIME_STEP_COLOR(1.0f, 0.0f, 0.0f);
static const tgt::vec3 LAST_TIME_STEP_COLOR = tgt::vec3(0.9f, 0.9f, 0.9f);
static const tgt::vec3 SELECTED_TIME_STEP_COLOR = tgt::vec3::zero;
static const tgt::vec3 MIN_DURATION_COLOR(1.0f, 0.0f, 0.0f);
static const tgt::vec3 MAX_DURATION_COLOR(0.0f, 0.0f, 1.0f);
static const tgt::vec3 FADE_OUT_COLOR = LAST_TIME_STEP_COLOR;
static const float SELECTED_LINE_WIDTH = 6.0f;
static const float UNSELECTED_LINE_WIDTH = 3.0f;
static const float PICKING_LINE_WIDTH = SELECTED_LINE_WIDTH;

const std::string DistanceToOriginPlot::fontName_("Vera.ttf");
const std::string DistanceToOriginPlot::loggerCat_("voreen.ensembleanalysis.DistanceToOriginPlot");

DistanceToOriginPlot::DistanceToOriginPlot()
    : RenderProcessor()
    , volumeListPort_(Port::INPORT, "distancetooriginVolumeListPort", "Feature VolumeList", false)
    , originalInport_(Port::INPORT, "ensembledatastructurein", "Ensemble Datastructure Input", false)
    , outport_(Port::OUTPORT, "outport", "Outport", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , privatePort_(Port::OUTPORT, "image.tmp", "image.tmp", false)
    , pickingBuffer_(Port::OUTPORT, "picking", "Picking", false)
    , eigenValueOutport_(Port::OUTPORT, "eigenvalueOutport", "Eigenvalues", false)
    , calculateButton_("calculate", "Calculate distances")
    , progressBar_("progressBar", "Progress")
    , sphereRadius_("sphereRadius", "Sphere Radius", 0.01, 0.0f, 0.1f)
    , fontSize_("fontSize", "Font Size", 10, 1, 30)
    , showTooltip_("showTooltip", "Show Tooltip", true)
    , renderTimeSelection_("renderTimeSelection", "Render Time Selection", false, Processor::INVALID_RESULT, Property::LOD_ADVANCED)
    , colorCoding_("colorCoding", "Color Coding")
    , distanceMethod_("distanceMethod", "Distance Method")
    , firstSelectedMember_("selectedMembers", "First selected Members")
    , selectedFeatureIdProperty_("selectedFeatureIdProperty", "Selected Feature IDs")

    , firstSelectedTimeInterval_("selectedTimeSteps", "First selected Time Interval", tgt::vec2(0.0f, 0.0f), 0.0f, 300.0f)
    , firstSelectedRadiusInterval_("firstSelectedRadiusInterval_", "First selected Radius Interval", tgt::vec2(0.0f, 0.0f), 0.0f, 230.0f)
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
    , largestFeatureVoxelCount_(0)
    , selectionStart_(tgt::vec2(0,0))
    , selectionEnd_(tgt::vec2(0,0))
    , isSelecting_(false)
    , temperatureRange_(tgt::vec2(9999999, -9999999))
    , vectormagnitudeRange_(tgt::vec2(9999999, -9999999))
    , condactivityanamolyRange_(tgt::vec2(9999999, -9999999))
    , expansivityanamolyRange_(tgt::vec2(9999999, -9999999))
    , temperatureanamolyRange_(tgt::vec2(9999999, -9999999))
    , densityanamolyRange_(tgt::vec2(9999999, -9999999))
    , transferFunc_("transferFunction", "Transfer Function", Processor::INVALID_RESULT)
    , selectField_("selectField", "Select field")
{
    // Ports
    addPort(volumeListPort_);
    addPort(originalInport_);

    //NOTE: don't use adjustPropertiesToInput() as callback, this seems to be triggered each frame for RenderProcessors.
    addPort(outport_);
    addPrivateRenderPort(privatePort_);
    addPrivateRenderPort(pickingBuffer_);

    // Calculation
    addProperty(calculateButton_);
        calculateButton_.setGroupID("calculation");
        calculateButton_.setReadOnlyFlag(false);
        ON_CHANGE(calculateButton_, DistanceToOriginPlot, calculateDistances);

    addProperty(progressBar_);
        progressBar_.setGroupID("calculation");
    addProgressBar(&progressBar_);

    setPropertyGroupGuiName("calculation", "Calculation");

    // Rendering
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

    addProperty(distanceMethod_);
    distanceMethod_.addOption("max", "Max. Distance to origin", "max");
    distanceMethod_.addOption("min", "Min. Distance to origin", "min");
    distanceMethod_.addOption("size", "Size", "size");
    distanceMethod_.selectByValue("max");
    distanceMethod_.setGroupID("rendering");

    addProperty(selectedFeatureIdProperty_);
    selectedFeatureIdProperty_.setGroupID("rendering");

    int lastId = 10000;
    for(int id = 0; id <= lastId; id++) {
        selectedFeatureIdProperty_.addRow(std::to_string(id));
    }

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
    //setPropertyGroupVisible("selection", false);

    addProperty(transferFunc_);
    addProperty(selectField_);
    selectField_.addOption("temperature", "Temperature", "temperature");
    selectField_.addOption("magnitude", "Vector Magnitude", "vectormagnitude");
    selectField_.addOption("thermal conductivity", "Conductivity Anamoly", "conductivityanamoly");
    selectField_.addOption("thermal expansivity", "Expansivity Anamoly", "expansivityanamoly");
    selectField_.addOption("temperature anomaly", "Temperature Anomaly", "temperatureanomaly");
    selectField_.addOption("spin transition-induced density anomaly", "Density Anamoly", "densityanamoly");
    selectField_.addOption("id", "ID", "id");
    selectField_.addOption("distancemax", "Distance Max","distancemax");
    ON_CHANGE_LAMBDA(selectField_, [this] {
        if(!transferFunc_.get()) {
            return;
        }

        if (selectField_.getValue() == "temperature") {
            transferFunc_.get()->setDomain(temperatureRange_);
        }
        if (selectField_.getValue() == "vectormagnitude") {
            transferFunc_.get()->setDomain(vectormagnitudeRange_);
        }
        if (selectField_.getValue() == "conductivityanamoly") {
            transferFunc_.get()->setDomain(condactivityanamolyRange_);
        }
        if (selectField_.getValue() == "expansivityanamoly") {
            transferFunc_.get()->setDomain(expansivityanamolyRange_);
        }
        if (selectField_.getValue() == "temperatureanomaly") {
            transferFunc_.get()->setDomain(temperatureanamolyRange_);
        }
        if (selectField_.getValue() == "densityanamoly") {
            transferFunc_.get()->setDomain(densityanamolyRange_);
        }
        if (selectField_.getValue() == "distancemax") {
            transferFunc_.get()->setDomain(tgt::vec2(0, 201));
        }
    });

    // IO
    addProperty(saveFileDialog_);
        ON_CHANGE(saveFileDialog_, DistanceToOriginPlot, save);
        saveFileDialog_.setGroupID("io");
    addProperty(saveButton_);
        ON_CHANGE(saveButton_, DistanceToOriginPlot, save);
        saveButton_.setGroupID("io");
    addProperty(loadFileDialog_);
        ON_CHANGE(loadFileDialog_, DistanceToOriginPlot, load);
        loadFileDialog_.setGroupID("io");
    addProperty(loadButton_);
        ON_CHANGE(loadButton_, DistanceToOriginPlot, load);
        loadButton_.setGroupID("io");
    setPropertyGroupGuiName("io", "Save/Load MDS Plot");

    // Camera
    addProperty(camera_);
    cameraHandler_ = new CameraInteractionHandler("cameraHandler", "Camera", &camera_);
    cameraHandler_->setEnabled(true);
    addInteractionHandler(cameraHandler_);
}

DistanceToOriginPlot::~DistanceToOriginPlot() {
    delete cameraHandler_;
}

Processor* DistanceToOriginPlot::create() const {
    return new DistanceToOriginPlot();
}

void DistanceToOriginPlot::initialize() {
    RenderProcessor::initialize();
    sphere_.setSphereGeometry(1.0f, tgt::vec3::zero, tgt::vec4::one, 16);
}

void DistanceToOriginPlot::deinitialize() {
    sphere_.clear();
    RenderProcessor::deinitialize();
}

void DistanceToOriginPlot::process() {

    // Resize frame buffers accordingly.
    if(pickingBuffer_.getSize() != outport_.getSize())
        pickingBuffer_.resize(outport_.getSize());

    if(privatePort_.getSize() != outport_.getSize())
        privatePort_.resize(outport_.getSize());

    // Change depth func to apply rendering order properly.
    glDepthFunc(GL_LEQUAL);
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

    // Render picking pass.
//    pickingBuffer_.activateTarget();
//    pickingBuffer_.clearTarget();
//    glLineWidth(PICKING_LINE_WIDTH);
//    renderingPass(true);
//    pickingBuffer_.deactivateTarget();

    // Draw the members in new parameter space.
    privatePort_.activateTarget();
    privatePort_.clearTarget();

    renderingPass(false);

    privatePort_.deactivateTarget();

    // Draw axes.
    renderAxes();

    // Draw tooltip.
    if(showTooltip_.get()) {
        outport_.activateTarget();
        renderTooltip();
        outport_.deactivateTarget();
    }

    // Restore state.
    glDepthFunc(GL_LESS);
    glLineWidth(1.0f);
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    IMode.color(tgt::col4::one);
}

void DistanceToOriginPlot::renderAxes() {

    outport_.activateTarget();
    outport_.clearTarget();

    tgt::ivec2 size = outport_.getSize();
    tgt::ivec2 margins = getMargins();

    // Set Plot status.
    plotLib_->setWindowSize(size);
    plotLib_->setAxesWidth(2.0f);
    plotLib_->setDrawingColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
    plotLib_->setLineWidth(1.0f);
    plotLib_->setMaxGlyphSize(1.0f);
    plotLib_->setMarginBottom(margins.y);
    plotLib_->setMarginTop(margins.y);
    plotLib_->setMarginLeft(margins.x);
    plotLib_->setMarginRight(margins.x);
    plotLib_->setMinimumScaleStep(32, PlotLibrary::X_AXIS);
    plotLib_->setMinimumScaleStep(32, PlotLibrary::Y_AXIS);
    plotLib_->setMinimumScaleStep(32, PlotLibrary::Z_AXIS);

    plotLib_->setDomain(Interval<plot_t>(0, volumeListPort_.getData()->size() - 1), PlotLibrary::X_AXIS);

    plotLib_->setDomain(Interval<plot_t>(3485, 6371), PlotLibrary::Y_AXIS);
//    plotLib_->setDomain(Interval<plot_t>(3485, 6371), PlotLibrary::Y_AXIS);
    plotLib_->setDimension(PlotLibrary::TWO);

    if (plotLib_->setRenderStatus()) {
        plotLib_->setDrawingColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
        plotLib_->renderAxes();
        plotLib_->setDrawingColor(tgt::Color(0, 0, 0, .5f));
        plotLib_->setFontSize(fontSize_.get() + 2);

        plotLib_->renderAxisLabel(PlotLibrary::X_AXIS, "t [2 MYr]");
        plotLib_->renderAxisLabel(PlotLibrary::Y_AXIS, "Distance to center");

        plotLib_->setDrawingColor(tgt::Color(0, 0, 0, .5f));
        plotLib_->setFontSize(fontSize_.get());
        plotLib_->setFontColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
        plotLib_->renderAxisScales(PlotLibrary::X_AXIS, false);
        plotLib_->setFontSize(fontSize_.get() + 2);
        plotLib_->renderAxisScales(PlotLibrary::Y_AXIS, true);
    }
    plotLib_->resetRenderStatus();

    // Plot data.
    float x0 = mapRange(margins.x, 0, size.x, -1.0f, 1.0f);
    float x1 = mapRange(size.x - margins.x, 0, size.x, -1.0f, 1.0f);
    float y0 = mapRange(margins.y, 0, size.y, -1.0f, 1.0f);
    float y1 = mapRange(size.y - margins.y, 0, size.y, -1.0f, 1.0f);

    tgt::TextureUnit::setZeroUnit();
    tgt::Texture* texture = privatePort_.getColorTexture();
    if(texture) {
        IMode.setTextureMode(tgt::ImmediateMode::TEX2D);
        texture->bind();

        IMode.begin(tgt::ImmediateMode::QUADS);
        IMode.texcoord(0.0f, 0.0f); IMode.vertex(x0, y0);
        IMode.texcoord(1.0f, 0.0f); IMode.vertex(x1, y0);
        IMode.texcoord(1.0f, 1.0f); IMode.vertex(x1, y1);
        IMode.texcoord(0.0f, 1.0f); IMode.vertex(x0, y1);
        IMode.end();

        IMode.setTextureMode(tgt::ImmediateMode::TEXNONE);
    }

    // Reset state.
    IMode.color(1.0f, 1.0f, 1.0f, 1.0f);
    outport_.deactivateTarget();
    LGL_ERROR;
}

void DistanceToOriginPlot::renderTooltip() const {
    return;
    // todo: do it
    /*
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

    */
}

bool DistanceToOriginPlot::isReady() const {

    if(!originalInport_.isReady()) {
        setNotReadyErrorMessage("No Ensemble connected");
        return false;
    }

    if(!volumeListPort_.isReady()) {
        return false;
    }

    if(volumeListPort_.getData()->size() == 0) {
        setNotReadyErrorMessage("Volume list is empty!");
        return false;
    }
    if(volumeListPort_.getData()->size() != originalInport_.getData()->getMaxNumTimeSteps()) {
        setNotReadyErrorMessage("Timestep count of feature do not match timestep count of ensemble");
    }
    return true;
}

void DistanceToOriginPlot::renderingPass(bool picking) {
    render(picking);
}

void DistanceToOriginPlot::render(bool picking) {
    if(!isReady()) return;
    VolumeRAMRepresentationLock volume(volumeListPort_.getData()->at(0));
    size_t maxRadius = volume->getDimensions().y;

    firstSelectedTimeInterval_.setMaxValue(300);
    firstSelectedRadiusInterval_.setMaxValue(300);

    // DRAW 660 KM LINE
    glLineWidth(SELECTED_LINE_WIDTH);
    IMode.begin(tgt::ImmediateMode::LINES);
    IMode.color(tgt::vec3(1., 0, 0));
    float l0 = mapRange(static_cast<size_t>(0), static_cast<size_t>(0), volumeListPort_.getData()->size() - 1, -1.0f, 1.0f);
    float l1 = mapRange(static_cast<size_t>(volumeListPort_.getData()->size() - 1), static_cast<size_t>(0), volumeListPort_.getData()->size() - 1, -1.0f, 1.0f);
    float l2 = mapRange(static_cast<float>(maxRadius - 45.966735967f), static_cast<float>(0), static_cast<float>(maxRadius), -1.0f, 1.0f);
    IMode.vertex(tgt::vec2(l0, l2));
    IMode.vertex(tgt::vec2(l1, l2));
    IMode.end();

    // DRAW 1600 KM LINE
    glLineWidth(SELECTED_LINE_WIDTH);
    IMode.begin(tgt::ImmediateMode::LINES);
    IMode.color(tgt::vec3(1., 0, 0));
    l0 = mapRange(static_cast<size_t>(0), static_cast<size_t>(0), volumeListPort_.getData()->size() - 1, -1.0f, 1.0f);
    l1 = mapRange(static_cast<size_t>(volumeListPort_.getData()->size() - 1), static_cast<size_t>(0), volumeListPort_.getData()->size() - 1, -1.0f, 1.0f);
    l2 = mapRange(static_cast<float>(maxRadius - 111.434511436f), static_cast<float>(0), static_cast<float>(maxRadius), -1.0f, 1.0f);
    IMode.vertex(tgt::vec2(l0, l2));
    IMode.vertex(tgt::vec2(l1, l2));
    IMode.end();

    const auto& selectedFeatures = selectedFeatureIdProperty_.getSelectedRowIndices();

    for(auto& feature : distances_) {

        unsigned int id = feature.first;

        std::vector<tgt::vec2> temptest = temperature_.find(id)->second;
        std::vector<tgt::vec2> magnitest = vectormagnitude_.find(id)->second;
        std::vector<tgt::vec2> condatest = condactivityanamoly_.find(id)->second;
        std::vector<tgt::vec2> expatest = expansivityanamoly_.find(id)->second;
        std::vector<tgt::vec2> tempantest = temperatureanomaly_.find(id)->second;
        std::vector<tgt::vec2> denstest = densityanamoly_.find(id)->second;


        std::vector<tgt::ivec3> timesteps = feature.second;
//        if ( std::find(selectedFeatures.begin(), selectedFeatures.end(), id-1) != selectedFeatures.end() )
//
//            glLineWidth(SELECTED_LINE_WIDTH);
//        else
//            glLineWidth(UNSELECTED_LINE_WIDTH);

        float saturation = 1.0;
        if(selectedFeatures.size() > 0) {
            if ( std::find(selectedFeatures.begin(), selectedFeatures.end(), id-1) != selectedFeatures.end() )
                saturation = 1.0;
            else
                saturation = 0.3;
        }

        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        size_t numTimeSteps = timesteps.size();
        for (size_t t = 0; t < numTimeSteps - 1; t++) {
            float val0 = -1;
            float val1 = -1;

            if(distanceMethod_.get() == "max") {
                val0 = timesteps[t].x;
                val1 = timesteps[t+1].x;
            }
            else if(distanceMethod_.get() == "min") {
                val0 = timesteps[t].y;
                val1 = timesteps[t+1].y;
            }
            else {
                val0 = timesteps[t].z;
                val1 = timesteps[t+1].z;
            }
            if (val0 < 0 || val0 > 999998) continue;
            float x0 = mapRange(t, static_cast<size_t>(0), volumeListPort_.getData()->size() - 1, -1.0f, 1.0f);
            float y0 = mapRange(static_cast<float>(val0), static_cast<float>(0), static_cast<float>(maxRadius), -1.0f, 1.0f);

            /*
            // THIS IS HOW WE ACCESS THE MEAN VALUES:
            // TODO: FROM HERE
            std::cout << "mean temp: " << temptest[t].x / temptest[t].y << std::endl;
            std::cout << "mean magnitude: " << magnitest[t].x / magnitest[t].y << std::endl;
            std::cout << "mean cond: " << condatest[t].x / condatest[t].y << std::endl;
            std::cout << "mean expansion: " << expatest[t].x / expatest[t].y << std::endl;
            std::cout << "mean temperature anamoly: " << tempantest[t].x / tempantest[t].y << std::endl;
            std::cout << "mean density: " << denstest[t].x / denstest[t].y << std::endl;
             */

            tgt::vec3 color = tgt::vec3::zero;
            tgt::vec3 color2 = tgt::vec3::zero;
            if(selectField_.get() == "id") {
                color = colors_.at(id);
                color2 = colors_.at(id);
            }
            else {
                float val = 0.0f;
                float val2 = 0.0f;
                if (selectField_.getValue() == "temperature") {
                    val = temptest[t].x / temptest[t].y;
                    val2 = temptest[t+1].x / temptest[t+1].y;
                }
                else if (selectField_.getValue() == "vectormagnitude") {
                    val = magnitest[t].x / magnitest[t].y;
                    val2 = magnitest[t+1].x / magnitest[t+1].y;
                }
                else if (selectField_.getValue() == "conductivityanamoly") {
                    val = condatest[t].x / condatest[t].y;
                    val2 = condatest[t+1].x / condatest[t+1].y;
                }
                else if (selectField_.getValue() == "expansivityanamoly") {
                    val = expatest[t].x / expatest[t].y;
                    val2 = expatest[t+1].x / expatest[t+1].y;
                }
                else if (selectField_.getValue() == "temperatureanomaly") {
                    val = tempantest[t].x / tempantest[t].y;
                    val2 = tempantest[t+1].x / tempantest[t+1].y;
                }
                else if (selectField_.getValue() == "densityanamoly") {
                    val = denstest[t].x / denstest[t].y;
                    val2 = denstest[t+1].x / denstest[t+1].y;
                }
                else if (selectField_.getValue() == "distancemax") {
                    val = y0;
                    val2 = y0;
                }
                color = getColorForValue(val);
                color2 = getColorForValue(val2);
            }

            IMode.begin(tgt::ImmediateMode::POINTS);
            glPointSize(3.f);
            IMode.color(tgt::vec4(color, saturation));
            IMode.vertex(tgt::vec2(x0, y0));
            IMode.end();
            glPointSize(1.);

            if(val1 >= 0 && val1 < 999998) {
                IMode.begin(tgt::ImmediateMode::LINES);
                IMode.color(tgt::vec4(color, saturation));

                float x1 = mapRange(t+1, static_cast<size_t>(0), volumeListPort_.getData()->size() - 1, -1.0f, 1.0f);
                float y1 = mapRange(static_cast<float>(val1), static_cast<float>(0), static_cast<float>(maxRadius), -1.0f, 1.0f);
//
//                if(std::abs(y1 - y0) > 100) {
//                    continue;
//                }

                if (selectField_.getValue() == "distancemax") {
                    color2 = getColorForValue(y1);
                }

                IMode.vertex(tgt::vec2(x0, y0));
                IMode.color(tgt::vec4(color2, saturation));
                IMode.vertex(tgt::vec2(x1, y1));
                IMode.end();
            }
            else {
//                IMode.begin(tgt::ImmediateMode::POINTS);
//                glPointSize(5.);
//                IMode.color(colors_.at(id));
//                IMode.vertex(tgt::vec2(x0, y0));
//                IMode.end();
//                glPointSize(1.);
                // todo draw point
            }

            // CHECK IF ADDITIONAL VALUES WORK:

        }


        glDisable(GL_BLEND);

        float x0 = selectionStart_.x;
        float x1 = selectionEnd_.x;
        float y0 = selectionStart_.y;
        float y1 = selectionEnd_.y;
        if(x0 > x1) std::swap(x0, x1);
        if(y0 > y1) std::swap(y0, y1);

        glDepthFunc(GL_ALWAYS);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_BLEND);
        IMode.color(0.0f, 0.0f, 0.0f, 0.3f);
        IMode.begin(tgt::ImmediateMode::LINES);

        IMode.vertex(x0, y0);
        IMode.vertex(x0, y1);

        IMode.vertex(x0, y1);
        IMode.vertex(x1, y1);

        IMode.vertex(x1, y1);
        IMode.vertex(x1, y0);

        IMode.vertex(x1, y0);
        IMode.vertex(x0, y0);

        IMode.end();
        glDisable(GL_BLEND);
        glDepthFunc(GL_LESS);
        glBlendFunc(GL_ONE, GL_ZERO);


    }
}

void DistanceToOriginPlot::renderCircle(tgt::vec2 position, float size, tgt::vec3 color, tgt::vec2 canvasSize) {
    tgt::vec2 aspect = (canvasSize.x > canvasSize.y) ? tgt::vec2(canvasSize.y / canvasSize.x, 1) : tgt::vec2(1, canvasSize.x / canvasSize.y);

    //render Cirlce
    size_t circleResolution = 64;
    IMode.begin(tgt::ImmediateMode::POLYGON);
    IMode.color(color);
    for (size_t i = 0; i < circleResolution; ++i) {
        float angle = (static_cast<float>(i) / static_cast<float>(circleResolution)) * (2 * tgt::PIf);
        tgt::vec2 direction(cosf(angle), sinf(angle));
        direction *= size;

        tgt::vec2 pos = position;
        IMode.vertex(tgt::vec3(pos.x, pos.y, 0));
    }
    IMode.end();

    //render black Line around
    IMode.begin(tgt::ImmediateMode::LINE_LOOP);
    glLineWidth(1.f);
    IMode.color(tgt::vec3(0, 0, 0));
    for (size_t i = 0; i < circleResolution; ++i) {
        float angle = (static_cast<float>(i) / static_cast<float>(circleResolution)) * (2 * tgt::PIf);
        tgt::vec2 direction(cosf(angle), sinf(angle));
        direction *= size;

        tgt::vec2 pos = position;
        IMode.vertex(tgt::vec3(pos.x, pos.y, 0));
    }
    IMode.end();
}

void DistanceToOriginPlot::renderSphere(const tgt::vec3& position, const tgt::vec3& color, bool drawBorder) const {

    // Skip rendering, if not visible anyways.
    if(sphereRadius_.get() <= std::numeric_limits<float>::epsilon())
        return;

    // TODO: Draw using fragment shader!

    MatStack.pushMatrix();
    MatStack.translate(position);
    MatStack.scale(tgt::vec3(sphereRadius_.get()));
    IMode.color(tgt::vec4(color, 1.0f));
    sphere_.render(GL_TRIANGLES);
    MatStack.popMatrix();
}

tgt::vec3 DistanceToOriginPlot::getColor(size_t memberIdx, size_t timeStepIdx, bool picking) const {
    return tgt::vec3(1,0,0);
}

void DistanceToOriginPlot::mouseEvent(tgt::MouseEvent* e) {
    // We always accept the event and force a redraw.
    e->accept();
    invalidate();

    // Reset last hit.
    lastHit_ = boost::none;

    int x = e->x();
    int y = e->y();

    // In 1D, and 2D case, we use margins to show axis and labels.
    tgt::ivec2 margins = getMargins();
    x = tgt::clamp(x, margins.x, e->viewport().x - margins.x);
    y = tgt::clamp(y, margins.y, e->viewport().y - margins.y);

    float t0 = -1;
    float t1 = -1;
    float r0 = -1;
    float r1 = -1;

    // Select a time interval outside of the axis
    if(e->action() == tgt::MouseEvent::PRESSED) {

        selectionStart_.x =mapRange(x, margins.x, e->viewport().x - margins.x, -1.f, 1.f);
        selectionStart_.y =  mapRange(y, margins.y, e->viewport().y - margins.y, 1.f, -1.f);
        selectionEnd_.x = selectionStart_.x;
        selectionEnd_.y = selectionStart_.y;

        t0 = mapRange(static_cast<float>(x), static_cast<float>(margins.x), static_cast<float>(e->viewport().x - margins.x), 0, static_cast<int>(volumeListPort_.getData()->size() - 1));
        firstSelectedTimeInterval_.set(tgt::vec2(t0, t0 + 1));

        VolumeRAMRepresentationLock volume(volumeListPort_.getData()->at(0));
        size_t maxRadius = volume->getDimensions().y;
        r0 = mapRange(y, margins.y, e->viewport().y - margins.y, 0, static_cast<int>(maxRadius - 1));
        firstSelectedRadiusInterval_.set(tgt::vec2(r0, r0 + 1));

        isSelecting_ = true;
    }
    else if(e->action() == tgt::MouseEvent::MOTION) {
        if(isSelecting_) {
            selectionEnd_.x =mapRange(x, margins.x, e->viewport().x - margins.x, -1.f, 1.f);
            selectionEnd_.y =  mapRange(y, margins.y, e->viewport().y - margins.y, 1.f, -1.f);
        }
    }
    else if(e->action() == tgt::MouseEvent::RELEASED) {
        selectionEnd_.x =mapRange(x, margins.x, e->viewport().x - margins.x, -1.f, 1.f);
        selectionEnd_.y =  mapRange(y, margins.y, e->viewport().y - margins.y, 1.f, -1.f);

        t1 = mapRange(x, margins.x, e->viewport().x - margins.x, 0, static_cast<int>(volumeListPort_.getData()->size() - 1));
        if(t1 > firstSelectedTimeInterval_.get().x) {
            firstSelectedTimeInterval_.set(tgt::vec2(firstSelectedTimeInterval_.get().x, t1));
        }
        else {
            firstSelectedTimeInterval_.set(tgt::vec2(t1, firstSelectedTimeInterval_.get().x));
        }

        VolumeRAMRepresentationLock volume(volumeListPort_.getData()->at(0));
        size_t maxRadius = volume->getDimensions().y;
        r1 = mapRange(y, margins.y, e->viewport().y - margins.y, 0, static_cast<int>(maxRadius - 1));
        if(r1 > firstSelectedRadiusInterval_.get().x) {
            firstSelectedRadiusInterval_.set(tgt::vec2(firstSelectedRadiusInterval_.get().x, r1));
        }
        else {
            firstSelectedRadiusInterval_.set(tgt::vec2(r1, firstSelectedRadiusInterval_.get().x));
        }
        isSelecting_ = false;


        t0 = firstSelectedTimeInterval_.get().x;
        t1 = firstSelectedTimeInterval_.get().y;
        r1 = maxRadius - firstSelectedRadiusInterval_.get().x;
        r0 = maxRadius - firstSelectedRadiusInterval_.get().y;

        std::vector<int> selectedRowIndices;
        for(auto& feature: distances_) {
            for(int i = t0; i <= t1; i++) {
                if(feature.second[i].x >= r0 && feature.second[i].x <= r1) {
                    selectedRowIndices.emplace_back(feature.first-1);
                    LINFO("Selected Feature: " << feature.first);
                    break;
                }
            }
        }
        selectedFeatureIdProperty_.setSelectedRowIndices(selectedRowIndices);
    }


    // Apply selection for both ranges.
    return;
}

void DistanceToOriginPlot::onEvent(tgt::Event* e) {

    // Events will be triggered for RenderProcessors, even if the processor is not ready it seems.
    if (isReady()) {
        tgt::MouseEvent* event = dynamic_cast<tgt::MouseEvent*>(e);
        if(event) {
            mouseEvent(event);
        }
    }

    RenderProcessor::onEvent(e);
}

void DistanceToOriginPlot::adjustToEnsemble() {

    subSelection_.clear();
    firstSelectedMember_.reset();
    secondSelectedMember_.reset();
    calculateButton_.setReadOnlyFlag(false);
}

void DistanceToOriginPlot::renderedMembersChanged() {
    renderingOrder_.clear();
}

tgt::ivec2 DistanceToOriginPlot::getMargins() const {
    return MIN_MARGINS;
}


void DistanceToOriginPlot::calculateDistances() {
    auto timesteps = volumeListPort_.getData();
    float progressStepSize = 1.f / timesteps->size();
    float currProgress = 0.f;
    setProgress(currProgress);

    for (size_t t = 0; t < timesteps->size(); t++) {

        VolumeRAMRepresentationLock volume(timesteps->at(t));
        
        // Carefull! Only works if ensemble is connected!
        
        const EnsembleDataset& ensemble = *originalInport_.getData();
        for (int i = 0; i < ensemble.getMembers()[0].getTimeSteps()[t].getFieldNames().size(); i++)
        {
            std::cout << ensemble.getMembers()[0].getTimeSteps()[t].getFieldNames()[i] << std::endl;
        }
        
        /* Fields in HDF5 Data:
            Acceleration Magnitude 0
            Rotated Velocity 1
            Velocity 2
            Vorticity Magnitude 3
            angle 4
            magnitude 5
            radius 6
            spin transition - induced density anomaly 7
            temperature 8
            temperature anomaly 9
            thermal conductivity 10
            thermal expansivity 11
        */

        LINFO(ensemble.getMembers()[0].getTimeSteps()[t].getFieldNames()[5]);
        LINFO(ensemble.getValueRange("magnitude"));

        //VolumeRAMRepresentationLock anglevol(ensemble.getMembers()[0].getTimeSteps()[t].getVolume(ensemble.getMembers()[0].getTimeSteps()[t].getFieldNames()[4]));
        VolumeRAMRepresentationLock magnitudevol(ensemble.getMembers()[0].getTimeSteps()[t].getVolume(ensemble.getMembers()[0].getTimeSteps()[t].getFieldNames()[5]));
        VolumeRAMRepresentationLock radiusvol(ensemble.getMembers()[0].getTimeSteps()[t].getVolume(ensemble.getMembers()[0].getTimeSteps()[t].getFieldNames()[6]));
        VolumeRAMRepresentationLock densvol(ensemble.getMembers()[0].getTimeSteps()[t].getVolume(ensemble.getMembers()[0].getTimeSteps()[t].getFieldNames()[7]));
        VolumeRAMRepresentationLock tempvol(ensemble.getMembers()[0].getTimeSteps()[t].getVolume(ensemble.getMembers()[0].getTimeSteps()[t].getFieldNames()[8]));
        VolumeRAMRepresentationLock tempanvol(ensemble.getMembers()[0].getTimeSteps()[t].getVolume(ensemble.getMembers()[0].getTimeSteps()[t].getFieldNames()[9]));
        VolumeRAMRepresentationLock condvol(ensemble.getMembers()[0].getTimeSteps()[t].getVolume(ensemble.getMembers()[0].getTimeSteps()[t].getFieldNames()[10]));
        VolumeRAMRepresentationLock expvol(ensemble.getMembers()[0].getTimeSteps()[t].getVolume(ensemble.getMembers()[0].getTimeSteps()[t].getFieldNames()[11]));

//#pragma omp parallel for
        for (long r = volume->getDimensions().y - 1; r >= 0; r-- ) {
            for (size_t phi = 0; phi < volume->getDimensions().x; phi++) {
                for(size_t theta = 0; theta < volume->getDimensions().z; theta++) {
                    uint32_t id = std::round(timesteps->at(t)->getRealWorldMapping().normalizedToRealWorld(volume->getVoxelNormalized(phi, r, theta)));
                    
                    if(id == 0) continue;

                    if(distances_.count(id) == 0) {
                        // Here we calculate a bunch of properties:
                        // Radius
                        std::vector<tgt::ivec3> radii(timesteps->size(), tgt::ivec3(-1, 999999, 0));

                        // Mean Temperature:
                        std::vector<tgt::vec2> meantemp(timesteps->size(), tgt::vec2(0.0, 0.0));

                        // Mean Velocity Magnitude
                        std::vector<tgt::vec2> meanvelocitymagnitude(timesteps->size(), tgt::vec2( 0.0, 0.0));

                        // Mean Conducitvity Anamoly
                        std::vector<tgt::vec2> meanvcondan(timesteps->size(), tgt::vec2(0.0, 0.0));

                        // Mean Expansivitiy Anamoly
                        std::vector<tgt::vec2> meanexpan(timesteps->size(), tgt::vec2(0.0, 0.0));

                        // Mean Temperature Anamoly
                        std::vector<tgt::vec2> meantempan(timesteps->size(), tgt::vec2(0.0, 0.0));

                        // Mean Density Anamoly
                        std::vector<tgt::vec2> meandensan(timesteps->size(), tgt::vec2(0.0, 0.0));


                        distances_.insert(std::pair<uint32_t , std::vector<tgt::ivec3>>(id, radii));


                        temperature_.insert(std::pair<uint32_t, std::vector<tgt::vec2>>(id, meantemp));
                        vectormagnitude_.insert(std::pair<uint32_t, std::vector<tgt::vec2>>(id, meanvelocitymagnitude));
                        condactivityanamoly_.insert(std::pair<uint32_t, std::vector<tgt::vec2>>(id, meanvcondan));
                        expansivityanamoly_.insert(std::pair<uint32_t, std::vector<tgt::vec2>>(id, meanexpan));
                        temperatureanomaly_.insert(std::pair<uint32_t, std::vector<tgt::vec2>>(id, meantempan));
                        densityanamoly_.insert(std::pair<uint32_t, std::vector<tgt::vec2>>(id, meandensan));


                        colors_.insert(std::pair<uint32_t , tgt::vec3>(id, tgt::vec3(float(rand()), float(rand()), float(rand())) / (float)RAND_MAX));
                    }
                    tgtAssert(t < distances_.find(id)->second.size(), "ERRÖR" );
                    
                    // Insert information:

                    // Distance
                    auto& iter = distances_.find(id)->second.at(t);
                    iter.x = std::max<int>(iter.x, r);
                    iter.y = std::min<int>(iter.y, r);
                    iter.z++;

                    // Mean Temperature:
                    auto& seciter = temperature_.find(id)->second.at(t);
                    seciter.x += tempvol->getVoxelNormalized(phi, r, theta);
                    seciter.y++;

                    // Mean Velocity Magnitude
                    auto& seciter2 = vectormagnitude_.find(id)->second.at(t);
                    seciter2.x += magnitudevol->getVoxelNormalized(phi, r, theta);
                    if(seciter2.x < 0) {
                        LINFO(seciter2.y);
                        LINFO(seciter2.x);
                        LINFO(magnitudevol->getVoxelNormalized(phi, r, theta));
                    }
                    seciter2.y++;

                    // Mean Conducitvity Anamoly
                    auto& seciter3 = condactivityanamoly_.find(id)->second.at(t);
                    seciter3.x += condvol->getVoxelNormalized(phi, r, theta);
                    seciter3.y++;

                    // Mean Expansivitiy Anamoly
                    auto& seciter4 = expansivityanamoly_.find(id)->second.at(t);
                    seciter4.x += expvol->getVoxelNormalized(phi, r, theta);
                    seciter4.y++;

                    // Mean Temperature Anamoly
                    auto& seciter5 = temperatureanomaly_.find(id)->second.at(t);
                    seciter5.x += tempanvol->getVoxelNormalized(phi, r, theta);
                    seciter5.y++;

                    // Mean Density Anamoly
                    auto& seciter6 = densityanamoly_.find(id)->second.at(t);
                    seciter6.x += densvol->getVoxelNormalized(phi, r, theta);
                    seciter6.y++;

                    largestFeatureVoxelCount_ = std::max(largestFeatureVoxelCount_, iter.z);
                }
            }
        }

        currProgress += progressStepSize;
        setProgress(currProgress);

        // find ranges
        for(auto it = distances_.begin(); it != distances_.end(); it++) {
            float val0 = -1;

            if(distanceMethod_.get() == "max") {
                val0 = it->second.at(t).x;
            }
            else if(distanceMethod_.get() == "min") {
                val0 = it->second.at(t).y;
            }
            else {
                val0 = it->second.at(t).z;
            }

            // if valid feature
            if(val0 >= 0 && val0 < 999998) {
                if(temperature_[it->first].at(t).y > 0) {
                    temperatureRange_.x = std::min(temperatureRange_.x, temperature_[it->first].at(t).x/temperature_[it->first].at(t).y);
                    temperatureRange_.y = std::max(temperatureRange_.y, temperature_[it->first].at(t).x/temperature_[it->first].at(t).y);
                    vectormagnitudeRange_.x = std::min(vectormagnitudeRange_.x, vectormagnitude_[it->first].at(t).x/vectormagnitude_[it->first].at(t).y);
                    vectormagnitudeRange_.y = std::max(vectormagnitudeRange_.y, vectormagnitude_[it->first].at(t).x/vectormagnitude_[it->first].at(t).y);
                    condactivityanamolyRange_.x = std::min(condactivityanamolyRange_.x, condactivityanamoly_[it->first].at(t).x/condactivityanamoly_[it->first].at(t).y);
                    condactivityanamolyRange_.y = std::max(condactivityanamolyRange_.y, condactivityanamoly_[it->first].at(t).x/condactivityanamoly_[it->first].at(t).y);
                    expansivityanamolyRange_.x = std::min(expansivityanamolyRange_.x, expansivityanamoly_[it->first].at(t).x/expansivityanamoly_[it->first].at(t).y);
                    expansivityanamolyRange_.y = std::max(expansivityanamolyRange_.y, expansivityanamoly_[it->first].at(t).x/expansivityanamoly_[it->first].at(t).y);
                    temperatureanamolyRange_.x = std::min(temperatureanamolyRange_.x, temperatureanomaly_[it->first].at(t).x/temperatureanomaly_[it->first].at(t).y);
                    temperatureanamolyRange_.y = std::max(temperatureanamolyRange_.y, temperatureanomaly_[it->first].at(t).x/temperatureanomaly_[it->first].at(t).y);
                    densityanamolyRange_.x = std::min(densityanamolyRange_.x, densityanamoly_[it->first].at(t).x/densityanamoly_[it->first].at(t).y);
                    densityanamolyRange_.y = std::max(densityanamolyRange_.y, densityanamoly_[it->first].at(t).x/densityanamoly_[it->first].at(t).y);
                }
            }
        }


//        for(auto it = temperature_.begin(); it != temperature_.end(); it++) {
//            if(it->second.at(t).y > 0) {
//                temperatureRange_.x = std::min(temperatureRange_.x, it->second.at(t).x/it->second.at(t).y);
//                temperatureRange_.y = std::max(temperatureRange_.y, it->second.at(t).x/it->second.at(t).y);
//            }
//        }
//
//        for(auto it = vectormagnitude_.begin(); it != vectormagnitude_.end(); it++) {
//            if(it->second.at(t).y > 0) {
//                vectormagnitudeRange_.x = std::min(vectormagnitudeRange_.x, it->second.at(t).x/it->second.at(t).y);
//                vectormagnitudeRange_.y = std::max(vectormagnitudeRange_.y, it->second.at(t).x/it->second.at(t).y);
//            }
//        }
//
//        for(auto it = condactivityanamoly_.begin(); it != condactivityanamoly_.end(); it++) {
//            if(it->second.at(t).y > 0) {
//                condactivityanamolyRange_.x = std::min(condactivityanamolyRange_.x, it->second.at(t).x/it->second.at(t).y);
//                condactivityanamolyRange_.y = std::max(condactivityanamolyRange_.y, it->second.at(t).x/it->second.at(t).y);
//            }
//        }
//
//        for(auto it = expansivityanamoly_.begin(); it != expansivityanamoly_.end(); it++) {
//            if(it->second.at(t).y > 0) {
//                expansivityanamolyRange_.x = std::min(expansivityanamolyRange_.x, it->second.at(t).x/it->second.at(t).y);
//                expansivityanamolyRange_.y = std::max(expansivityanamolyRange_.y, it->second.at(t).x/it->second.at(t).y);
//            }
//        }
//
//        for(auto it = temperatureanomaly_.begin(); it != temperatureanomaly_.end(); it++) {
//            if(it->second.at(t).y > 0) {
//                temperatureanamolyRange_.x = std::min(temperatureanamolyRange_.x, it->second.at(t).x/it->second.at(t).y);
//                temperatureanamolyRange_.y = std::max(temperatureanamolyRange_.y, it->second.at(t).x/it->second.at(t).y);
//            }
//        }
//
//        for(auto it = densityanamoly_.begin(); it != densityanamoly_.end(); it++) {
//            if(it->second.at(t).y > 0) {
//                densityanamolyRange_.x = std::min(densityanamolyRange_.x, it->second.at(t).x/it->second.at(t).y);
//                densityanamolyRange_.y = std::max(densityanamolyRange_.y, it->second.at(t).x/it->second.at(t).y);
//            }
//        }
    }
    LINFO("Temp. range: " << temperatureRange_.x << " " << temperatureRange_.y);
    LINFO("Vector magn. range: " << vectormagnitudeRange_.x << " " << vectormagnitudeRange_.y);
    LINFO("Therm. Cond. range: " << condactivityanamolyRange_.x << " " << condactivityanamolyRange_.y);
    LINFO("Therm. Exp. range: " << expansivityanamolyRange_.x << " " << expansivityanamolyRange_.y);
    LINFO("Temp. Anom. range: " << temperatureanamolyRange_.x << " " << temperatureanamolyRange_.y);
    LINFO("Density Anom. range: " << densityanamolyRange_.x << " " << densityanamolyRange_.y);

}

void DistanceToOriginPlot::save() {
    if (saveFileDialog_.get().empty()) {
        LWARNING("no filename specified");
        return;
    }
    if(distances_.empty()) {
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
        serializer.serialize("distances", distances_);
        serializer.serialize("temperature", temperature_);
        serializer.serialize("vectormagnitude", vectormagnitude_);
        serializer.serialize("condactivity", condactivityanamoly_);
        serializer.serialize("expansivity", expansivityanamoly_);
        serializer.serialize("temperatureanamoly", temperatureanomaly_);
        serializer.serialize("density", densityanamoly_);
        serializer.serialize("colors", colors_);
        serializer.serialize("largestFeatureVoxelCount", largestFeatureVoxelCount_);
        serializer.serialize("temperatureRange", temperatureRange_);
        serializer.serialize("vectormagnitudeRange", vectormagnitudeRange_);
        serializer.serialize("condactivityanamolyRange", condactivityanamolyRange_);
        serializer.serialize("expansivityanamolyRange", expansivityanamolyRange_);
        serializer.serialize("temperatureanamolyRange", temperatureanamolyRange_);
        serializer.serialize("densityanamolyRange", densityanamolyRange_);
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

void DistanceToOriginPlot::load() {
    if (loadFileDialog_.get().empty()) {
        LWARNING("no filename specified");
        return;
    }

    std::string absPath = loadFileDialog_.get();
    try {
        std::ifstream inFile;
        inFile.open(absPath.c_str());

        JsonDeserializer d;
        d.read(inFile);
        Deserializer deserializer(d);
        deserializer.deserialize("distances", distances_);
        deserializer.deserialize("temperature", temperature_);
        deserializer.deserialize("vectormagnitude", vectormagnitude_);
        deserializer.deserialize("condactivity", condactivityanamoly_);
        deserializer.deserialize("expansivity", expansivityanamoly_);
        deserializer.deserialize("temperatureanamoly", temperatureanomaly_);
        deserializer.deserialize("density", densityanamoly_);
        deserializer.deserialize("colors", colors_);
        deserializer.deserialize("largestFeatureVoxelCount", largestFeatureVoxelCount_);
        deserializer.deserialize("temperatureRange", temperatureRange_);
        deserializer.deserialize("vectormagnitudeRange", vectormagnitudeRange_);
        deserializer.deserialize("condactivityanamolyRange", condactivityanamolyRange_);
        deserializer.deserialize("expansivityanamolyRange", expansivityanamolyRange_);
        deserializer.deserialize("temperatureanamolyRange", temperatureanamolyRange_);
        deserializer.deserialize("densityanamolyRange", densityanamolyRange_);

        inFile.close();
        LINFO("Reading embeddings from file " << absPath << " was successful");
    }
    catch(tgt::Exception& e) {
        VoreenApplication::app()->showMessageBox("loading embeddings failed", e.what(), true);
        LERROR(e.what());
        loadFileDialog_.set("");
    }
}

tgt::vec3 DistanceToOriginPlot::getColorForValue(float val) const {
    auto keys = transferFunc_.get()->getKeys();
    if(keys.size() == 0)
        return tgt::vec3::one;
    int index = 0;
    tgt::vec2 range = transferFunc_.get()->getDomain();
    //LINFO("Domain " << range.x << " " << range.y << " " << val);
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

} // namespace
