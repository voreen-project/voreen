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

#include "starcoordinates.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include "modules/ensembleanalysis/utils/utils.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/interaction/camerainteractionhandler.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "voreen/core/utils/multisampler.h"

#include "modules/plotting/datastructures/plotcell.h"
#include "modules/plotting/datastructures/plotdata.h"
#include "modules/plotting/utils/plotlibrary/plotlibrary.h"
#include "modules/plotting/utils/plotlibrary/plotlibraryopengl.h"

#include "tgt/immediatemode/immediatemode.h"
#include "tgt/texture.h"

#include <random>



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
    static const tgt::vec3 NOTHING_COLOR(0.5f, 0.5f, 0.5f);
    static const tgt::vec3 SLAB_COLOR(0.0f, 0.0f, 1.0f);
    static const tgt::vec3 PLUME_COLOR(1.0f, 0.0f, 0.0f);

    const float PI = 3.14159265358979f;

    static const float BACKGROUND_ID = 0.f;
    static const float AXIS_ID = 0.5f;
    static const float DATA_ID = 1.f;

    const std::string StarCoordinates::fontName_("Vera.ttf");
    const std::string StarCoordinates::loggerCat_("voreen.scivis21.StarCoordinates");

    StarCoordinates::StarCoordinateAxis::StarCoordinateAxis()
        : label_(""), color_(tgt::vec3(0)), direction_(tgt::vec2(0)), length_(0.f), id_(0), enabled_(true) { }
    StarCoordinates::StarCoordinateAxis::StarCoordinateAxis(std::string label, tgt::vec3 color, tgt::vec2 direction, float length, size_t id)
    : label_(label), color_(color), direction_(direction), length_(length), id_(id), enabled_(true) { }

    void StarCoordinates::StarCoordinateAxis::serialize(Serializer& s) const {
        s.serialize("color", color_);
        s.serialize("direction", direction_);
        s.serialize("length", length_);
        s.serialize("axisid", id_);
        s.serialize("label", label_);
        s.serialize("enabled", enabled_);
    }
    void StarCoordinates::StarCoordinateAxis::deserialize(Deserializer& s) {
        s.deserialize("color", color_);
        s.deserialize("direction", direction_);
        s.deserialize("length", length_);
        s.deserialize("axisid", id_);
        s.deserialize("label", label_);
        s.deserialize("enabled", enabled_);
    }

    StarCoordinates::StarCoordinatePoint::StarCoordinatePoint()
        : seriesID_(0), timeStepID_(0) { }
    StarCoordinates::StarCoordinatePoint::StarCoordinatePoint(size_t seriesID, size_t timeStepID)
        : seriesID_(seriesID), timeStepID_(timeStepID) { }

    void StarCoordinates::StarCoordinatePoint::serialize(Serializer& s) const {
        s.serialize("seriesid", seriesID_);
        s.serialize("timestepid", timeStepID_);
    }
    void StarCoordinates::StarCoordinatePoint::deserialize(Deserializer& s) {
        s.deserialize("seriesid", seriesID_);
        s.deserialize("timestepid", timeStepID_);
    }

    StarCoordinates::StarCoordinates()
        : RenderProcessor()
        , timeseriesInport_(Port::INPORT, "timeseriesInport", "Timeseries Input", false)
        , outport_(Port::OUTPORT, "outport", "Outport", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
        , privatePort_(Port::OUTPORT, "image.tmp", "image.tmp", false)
        , pickingBuffer_(Port::OUTPORT, "picking", "Picking", false)
        , camera_("camera", "Camera", tgt::Camera(tgt::vec3(0.0f, 0.0f, 3.5f), tgt::vec3(0.0f, 0.0f, 0.0f), tgt::vec3(0.0f, 1.0f, 0.0f)))
        , cameraHandler_(nullptr)
        , calculateButton_("calculate", "Create Embeddings")
        , renderedMembers_("renderedMembers", "Rendered Members")
        , selectedMembers_("selectedMembers", "Selected Members")
        , plotLib_(new PlotLibraryOpenGl())
        , sphereRadius_("sphere.radius", "Sphere Radius", 0.02f, 1e-6f, 10.f)
        , zoom_(0.6f)
        , mouseWheelZoomAcuteness_(4.f)
        , fontScaling_("font.scaling", "Font Scaling", 18.f, 1e-6f, 100.f)
        , fontSize_("fontSize", "Font Size", 10, 1, 30)
        , colorCoding_("colorCoding", "Color Coding")
        , labelOffset_("labelOffset", "Label Offset", 0.01f, 1e-6f, 100.f)
        , pointSize_("point.size", "Point Size", 4.f, 1.f, 100.f)
        //, angleRange_("angleRange", "Angle Range", 90.0, 0.0, 180.0)
        , tempAnomalyRange_("tempAnomalyRange", "Temp. Anomaly Range", 0.0, -1000, 1000)
        , transferFunc_("transferFunction", "Transfer Function", Processor::INVALID_RESULT)
        , resetAxesButton_("reset.axes", "Reset axes")
        , resetSelectionButton_("resetSelectionButton", "Reset selection")
        , showTooltip_("showTooltip", "Show Tooltip", true)
        , volumeDimension_("volume", "Volume", 100, 1, 2048)
        , position_("position", "Position", tgt::ivec3::zero, tgt::ivec3::zero, tgt::ivec3::one*2048)
    {
        // Ports
        //NOTE: don't use adjustPropertiesToInput() as callback, this seems to be triggered each frame for RenderProcessors.
        addPort(timeseriesInport_);
        ON_CHANGE(timeseriesInport_, StarCoordinates, adjustToTimeserieslist);
        addPort(outport_);
        addPrivateRenderPort(privatePort_);
        addPrivateRenderPort(pickingBuffer_);

        // Camera
        addProperty(camera_);
        cameraHandler_ = new CameraInteractionHandler("cameraHandler", "Camera", &camera_);
        cameraHandler_->setEnabled(true);
        addInteractionHandler(cameraHandler_);
        addProperty(sphereRadius_);

        // Calculation
        addProperty(calculateButton_);
        calculateButton_.setGroupID("calculation");
        calculateButton_.setReadOnlyFlag(true);
        //ON_CHANGE(calculateButton_, TimeseriesPlot, createEmbeddings);
        addProperty(renderedMembers_);
        //ON_CHANGE(renderedMembers_, StarCoordinates, renderedMembersChanged);
        renderedMembers_.setGroupID("rendering");
        addProperty(selectedMembers_);
        selectedMembers_.setGroupID("rendering");
        addProperty(fontScaling_);

        addProperty(colorCoding_);
        //colorCoding_.addOption("member", "Only Member", COLOR_MEMBER);
        //colorCoding_.addOption("timeStep", "Only Time Step", COLOR_TIMESTEP);
        //colorCoding_.addOption("memberAndTimeStep", "Member and Time Step", COLOR_MEMBER_AND_TIMESTEP);
        //colorCoding_.selectByValue(COLOR_MEMBER_AND_TIMESTEP);
        //colorCoding_.setGroupID("rendering");
        ON_CHANGE_LAMBDA(colorCoding_, [this] {
            if (colorCoding_.getValue() == COLOR_FIELD) {
                int index = (colorCoding_.getSelectedIndex() - 3) / 2;
                transferFunc_.get()->setDomain(getRange(index));
            }
        });

        addProperty(labelOffset_);
        addProperty(pointSize_);
        addProperty(tempAnomalyRange_);
        addProperty(transferFunc_);
        addProperty(resetAxesButton_);
        ON_CHANGE(resetAxesButton_, StarCoordinates, resetAxes);
        addProperty(resetSelectionButton_);
        ON_CHANGE(resetSelectionButton_, StarCoordinates, resetSelection);

        addProperty(showTooltip_);
        showTooltip_.setGroupID("rendering");

        addProperty(fontSize_);
        fontSize_.setGroupID("rendering");

        //Linking
        addProperty(volumeDimension_);
        ON_CHANGE_LAMBDA(volumeDimension_, [this] {
            position_.setMaxValue(tgt::ivec3::one*volumeDimension_.get());
        });
        volumeDimension_.setVisibleFlag(false);
        addProperty(position_);
        position_.setVisibleFlag(false);

        updatedProjection = true;
    }
    tgt::vec2 StarCoordinates::getRange(size_t fieldIdx) {
        tgt::vec2 range(1.f, 0.f);
        for (auto series : timeseriesInport_.getData()->getAllTimeSeries()) {
            for (auto step : series.getTimeSeriesSteps()) {
                range.x = std::min(range.x, step.fieldValues_[fieldIdx]);
                range.y = std::max(range.y, step.fieldValues_[fieldIdx]);
            }
        }
        return range;
    }

    StarCoordinates::~StarCoordinates() {
        delete cameraHandler_;
    }

    Processor* StarCoordinates::create() const {
        return new StarCoordinates();
    }

    void StarCoordinates::serialize(Serializer& s) const {
        RenderProcessor::serialize(s);
        s.serialize("starcoordinateaxes", starCoordinateAxes_);
        s.serialize("starcoordinatepoints", starCoordinatePoints_);
        s.serialize("zoom", zoom_);
        s.serialize("selectedaxis", selectedAxis_);
        s.serialize("updateaxes", updateAxes);
    }
    void StarCoordinates::deserialize(Deserializer& s) {
        RenderProcessor::deserialize(s);
        s.deserialize("starcoordinateaxes", starCoordinateAxes_);
        s.deserialize("starcoordinatepoints", starCoordinatePoints_);
        s.deserialize("zoom", zoom_);
        s.deserialize("selectedaxis", selectedAxis_);
        s.deserialize("updateaxes", updateAxes);
    }

    void StarCoordinates::initialize() {
        RenderProcessor::initialize();

        bypassProgram_ = ShdrMgr.loadSeparate("passthrough.vert", "copyimage.frag", generateHeader(), false);

        sphere_.setSphereGeometry(1.0f, tgt::vec3::zero, tgt::vec4::one, 16);

        // initialize selection
        selectedAxis_ = -1;
        updateAxes = false;
    }

    void StarCoordinates::deinitialize() {
        ShdrMgr.dispose(bypassProgram_);
        bypassProgram_ = 0;
        LGL_ERROR;

        starCoordinateAxes_.clear();

        sphere_.clear();

        RenderProcessor::deinitialize();
    }

    void StarCoordinates::process() {

        // Resize frame buffers accordingly.
        if (pickingBuffer_.getSize() != outport_.getSize())
            pickingBuffer_.resize(outport_.getSize());

        if (privatePort_.getSize() != outport_.getSize())
            privatePort_.resize(outport_.getSize());

        // Set state
        glDepthFunc(GL_LEQUAL);
        glLineWidth(4.0f);

        // update projection in case axis has moved
        if (updatedProjection) {
            createProjection();
            projectPoints();
            updatedProjection = false;
        }

        // Render picking pass.
        pickingBuffer_.activateTarget();
        renderingPass(true);
        pickingBuffer_.deactivateTarget();

        // Draw the members in new parameter space.
        {
            //Multisampler m(outport_);
            outport_.activateTarget();
            renderingPass(false);
            outport_.deactivateTarget();
        }

        // Restore state.
        glDepthFunc(GL_LESS);
        glLineWidth(1.0f);
        IMode.color(tgt::col4::one);

        updateAxes = false;
    }

    void StarCoordinates::renderingPass(bool picking) {
        setPlotStatus();
        if (plotLib_->setRenderStatus()) {
            if (picking) {
                glClearColor(0.f, 0.f, 0.f, 0.f);
                pickingBuffer_.clearTarget();
                renderSelection();
                renderProjection(picking);
            }
            else {
                renderAxes();
                renderProjection(picking);
            }
        }
        plotLib_->resetRenderStatus();
    }

    void StarCoordinates::setPlotStatus() {
        plotLib_->setWindowSize(outport_.getSize());
        plotLib_->setFontSize(fontScaling_.get());
        plotLib_->setLineWidth(4.f);
        plotLib_->setMaxGlyphSize(0);
        plotLib_->setMarginBottom(outport_.getSize().y / 2.f);
        plotLib_->setMarginTop(0);
        plotLib_->setMarginLeft(outport_.getSize().x / 2.f);
        plotLib_->setMarginRight(0);
    }

    void StarCoordinates::renderAxes() {
        // extract window size and aspect ratio
        tgt::vec2 size = outport_.getSize();
        tgt::vec2 aspect = (size.x > size.y) ? tgt::vec2(size.y / size.x, 1) : tgt::vec2(1, size.x / size.y);
        
        // draw circle that indicates 1 length
        size_t circleResolution = 64;
        tgt::vec3 circleColor(0);
        glLineWidth(4.f);
        IMode.begin(tgt::ImmediateMode::LINE_LOOP);
        for (size_t i = 0; i < circleResolution; ++i) {
            IMode.color(0, 0, 0);

            float angle = (static_cast<float>(i) / static_cast<float>(circleResolution)) * (2 * tgt::PIf);
            tgt::vec2 direction(cosf(angle), sinf(angle));

            tgt::vec2 pos = direction * aspect * zoom_;
            IMode.vertex(tgt::vec3(pos.x, pos.y, 0));
        }
        IMode.end();

        // draw axes
        IMode.begin(tgt::ImmediateMode::LINES);
        int axisIdx = 0;
        for (auto& axis : starCoordinateAxes_) {
            tgt::vec3 color = (axis.enabled_) ? tgt::vec3(0) : tgt::vec3(0.4);
            IMode.color(color.r, color.g, color.b);
            IMode.vertex(tgt::vec3(0, 0, 0));

            tgt::vec2 pos = axis.direction_ * axis.length_ * aspect * zoom_;
            IMode.vertex(tgt::vec3(pos.x, pos.y, 0));

            pos = (axis.direction_ * axis.length_ + axis.direction_ * labelOffset_.get()) * aspect * zoom_;
            tgt::vec2 absolutePos = size / 2.f + pos * size / 2.f;
            auto alignment = (axis.direction_.x >= 0) 
                ? SmartLabel::Alignment::MIDDLELEFT 
                : SmartLabel::Alignment::MIDDLERIGHT;
            plotLib_->renderLabel(absolutePos, alignment, axis.label_);

            ++axisIdx;
        }
        IMode.end();

        // draw line handles
        glPointSize(10.f);
        IMode.begin(tgt::ImmediateMode::POINTS);
        for (auto& axis : starCoordinateAxes_) {
            tgt::vec3 color = (axis.enabled_) ? tgt::vec3(1, 0, 0) : tgt::vec3(0.1);
            IMode.color(color.r, color.g, color.b);
            tgt::vec2 pos = axis.direction_ * axis.length_ * aspect * zoom_;
            IMode.vertex(tgt::vec3(pos.x, pos.y, 0));
        }
        IMode.end();
    }

    void StarCoordinates::renderProjection(bool picking) {
        tgt::vec2 size = outport_.getSize();
        tgt::vec2 aspect = (size.x > size.y) ? tgt::vec2(size.y / size.x, 1) : tgt::vec2(1, size.x / size.y);

        // Check for data set.
        const TimeSeriesList* list = timeseriesInport_.getData();
        if (!list)
            return;
        size_t numberOfTimeSeries = list->getNumberOfTimeSeries();
        size_t numberOfTimeSteps = list->getTimeSeries(0).getNumberOfTimeSteps();

        // draw lines
        glLineWidth(pointSize_.get()*0.5f);
        IMode.begin(tgt::ImmediateMode::LINES);
        std::vector<int> selection = selectedMembers_.get();
        float alpha = 1.0;
        for (size_t y = 0; y < numberOfTimeSeries; ++y) {
            if(selection.size() > 0 && picking == false) {
                alpha = 0.03;
            }
            if(std::find(selection.begin(), selection.end(), y)!=selection.end())
                alpha = 1.0;
            size_t off = y * numberOfTimeSteps;
            for (size_t x = 1; x < numberOfTimeSteps; ++x) {
                size_t idx1 = off + x - 1;
                size_t idx2 = off + x;

                StarCoordinatePoint point = starCoordinatePoints_.at(idx1);
                tgt::vec3 color = getColor(point.seriesID_, point.timeStepID_, picking);
                IMode.color(tgt::vec4(color, alpha));
                float px = projectedPoints_(0, idx1);
                float py = projectedPoints_(1, idx1);
                tgt::vec2 pos = tgt::vec2(px, py) * aspect * zoom_;
                IMode.vertex(tgt::vec3(pos.x, pos.y, 0));

                point = starCoordinatePoints_.at(idx2);
                color = getColor(point.seriesID_, point.timeStepID_, picking);
                IMode.color(tgt::vec4(color, alpha));
                px = projectedPoints_(0, idx2);
                py = projectedPoints_(1, idx2);
                pos = tgt::vec2(px, py) * aspect * zoom_;
                IMode.vertex(tgt::vec3(pos.x, pos.y, 0));
            }
        }
        IMode.end();

        // iterate over matrix columnwise
        glPointSize(pointSize_.get());
        IMode.begin(tgt::ImmediateMode::POINTS);
        for (size_t c = 0; c < projectedPoints_.cols(); ++c) {
            StarCoordinatePoint point = starCoordinatePoints_.at(c);
            if(selection.size() > 0 && picking == false) {
                alpha = 0.03;
            }
            if(std::find(selection.begin(), selection.end(), point.seriesID_)!=selection.end())
                alpha = 1.0;
            tgt::vec3 color = getColor(point.seriesID_, point.timeStepID_, picking);
            IMode.color(tgt::vec4(color, alpha));
            float x = projectedPoints_(0, c);
            float y = projectedPoints_(1, c);
            tgt::vec2 pos = tgt::vec2(x, y) * aspect * zoom_;
            IMode.vertex(tgt::vec3(pos.x, pos.y, 0));
        }
        IMode.end();

        // Draw tooltip.
        if (showTooltip_.get()) {
            renderTooltip();
        }
    }

    void StarCoordinates::renderSelection() {
        tgt::vec2 size = outport_.getSize();
        tgt::vec2 aspect = (size.x > size.y) ? tgt::vec2(size.y / size.x, 1) : tgt::vec2(1, size.x / size.y);

        for (auto& axis : starCoordinateAxes_) {
            tgt::vec2 position = axis.direction_ * axis.length_ * aspect * zoom_;

            size_t axisCount = starCoordinateAxes_.size();
            float axisIdx = (axisCount > 1) ? static_cast<float>(axis.id_) / (axisCount - 1) : 0.f;
            tgt::vec3 color(axisIdx, 0, AXIS_ID);

            MatStack.pushMatrix();
            MatStack.translate(tgt::vec3(position, 0));
            MatStack.scale(tgt::vec3(sphereRadius_.get()));

            IMode.color(tgt::vec4(color, 1.0f));
            sphere_.render(GL_TRIANGLES);

            MatStack.popMatrix();
        }
    }

    void StarCoordinates::renderTooltip() const {

        const TimeSeriesList* timeserieslist = timeseriesInport_.getData();
        if (!timeserieslist || !lastHit_) {
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
        tgt::vec2 ur = mapRange(pos.xy() + size, tgt::vec2::zero, tgt::vec2(screensize), -tgt::vec2::one, tgt::vec2::one);

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

    tgt::vec3 StarCoordinates::getColor(size_t seriesIdx, size_t timeStepIdx, bool picking) const {
        const TimeSeriesList* list = timeseriesInport_.getData();
        const TimeSeries& series = list->getTimeSeries(seriesIdx);

        float ts = static_cast<float>(timeStepIdx) / series.getNumberOfTimeSteps();
        if (picking)
            return tgt::vec3(static_cast<float>(seriesIdx) / list->getNumberOfTimeSeries(), ts, DATA_ID);

        tgt::vec3 colorPosition = series.getPosition();// use position as color
        colorPosition.x = ((float)colorPosition.x) / list->getVolumeDimensions().x;
        colorPosition.y = ((float)colorPosition.y) / list->getVolumeDimensions().y;
        colorPosition.z = ((float)colorPosition.z) / list->getVolumeDimensions().z;

        int tempAnomalyIndex = 0;
        //int angleIndex = 0;
        list->getNumberOfFields();
        //int numFields = series.getTimeSeriesStep(timeStepIdx).fieldValues_.size();
        int numFields = list->getNumberOfFields();
        for (int i = 0; i < numFields; i++) {
//            if (list->getFieldName(i) == "angle") {
//                //angleIndex = i;
//            }
            if (list->getFieldName(i) == "temperature anomaly") {
                tempAnomalyIndex = i;
            }
        }
//        float angleVal = series.getTimeSeriesStep(timeStepIdx).fieldValues_[angleIndex];
//        tgt::vec2 angleRange = list->getRange("angle");
//        angleVal = angleVal * (angleRange.y - angleRange.x) + angleRange.x;
        float tempAnomalyVal = series.getTimeSeriesStep(timeStepIdx).fieldValues_[tempAnomalyIndex];
        tgt::vec2 tempAnomalyRange = list->getRange("temperature anomaly");
        tempAnomalyVal = tempAnomalyVal * (tempAnomalyRange.y - tempAnomalyRange.x) + tempAnomalyRange.x;

        switch (colorCoding_.getValue()) {
        case COLOR_MEMBER:
            return colorPosition;
        case COLOR_TIMESTEP:
            return  (1.0f - ts) * FIRST_TIME_STEP_COLOR + ts * LAST_TIME_STEP_COLOR;
        case COLOR_MEMBER_AND_TIMESTEP:
            return (1.0f - ts) * colorPosition + ts * FADE_OUT_COLOR;
        case COLOR_FIELD:
        {
            int index = (colorCoding_.getSelectedIndex() - 3);
            float val = series.getTimeSeriesStep(timeStepIdx).fieldValues_[index];
            tgt::vec3 color = getColorForValue(val);
            return color;
        }
        case COLOR_SLAB:
        {
            if (tempAnomalyVal < tempAnomalyRange_.get().x) {
                return SLAB_COLOR;
            }
            return NOTHING_COLOR;
        }
        case COLOR_PLUME:
        {
            if (tempAnomalyVal > tempAnomalyRange_.get().y) {
                return PLUME_COLOR;
            }
            return NOTHING_COLOR;
        }
        case COLOR_SLAB_AND_PLUME:
        {
            if (tempAnomalyVal < tempAnomalyRange_.get().x) {
                return SLAB_COLOR;
            }
            if (tempAnomalyVal > tempAnomalyRange_.get().y) {
                return PLUME_COLOR;
            }
            return NOTHING_COLOR;
        }
        default:
            return tgt::vec3::one;
        }
    }

    tgt::vec3 StarCoordinates::getColorForValue(float val) const {
        auto keys = transferFunc_.get()->getKeys();
        if (keys.size() == 0)
            return tgt::vec3::one;
        int index = 0;
        tgt::vec2 range = transferFunc_.get()->getDomain();
        float diff = range.y - range.x;
        while (index < keys.size() && val > keys[index]->getIntensity()* diff + range.x) {
            index++;
        }
        if (index == 0) {
            return keys[index]->getColorR().xyz();
        }
        else if (index == keys.size()) {
            return keys[index - 1]->getColorR().xyz();
        }
        else {
            float distanceUp = keys[index]->getIntensity() * diff + range.x - val;
            float distanceDown = val - keys[index - 1]->getIntensity() * diff + range.x;
            float distance = distanceUp + distanceDown;
            tgt::vec3 color = distanceUp / distance * (tgt::vec3) keys[index - 1]->getColorR().xyz() +
                distanceDown / distance * (tgt::vec3) keys[index]->getColorL().xyz();
            return color / 255.f;
        }
    }

    void StarCoordinates::adjustToTimeserieslist() {

        subSelection_.clear();
        renderedMembers_.reset();
        selectedMembers_.reset();
        calculateButton_.setReadOnlyFlag(true);
        colorCoding_.setOptions(std::deque<Option<ColorCoding>>());

        // Check for data set.
        const TimeSeriesList* list = timeseriesInport_.getData();
        if (!list)
            return;

        // initialize selection
        selectedAxis_ = -1;

        // check if the new axis match the previous axes
        size_t axisCount = list->getNumberOfFields();
        auto fieldNames = list->getFieldNames();
        bool fieldsChanged = false;
        if (starCoordinateAxes_.size() > 0 && starCoordinateAxes_.size() == axisCount) {
            for (size_t i = 0; i < axisCount; ++i) {
                std::string label = fieldNames.at(i);
                if (starCoordinateAxes_.at(i).label_.compare(label) != 0) {
                    fieldsChanged == true;
                }
            }
        }
        else {
            fieldsChanged = true;
        }
        
        // update axes if the fields changed
        if (fieldsChanged) {
            starCoordinateAxes_.clear();

            float radius = 0.8f;
            for (size_t i = 0; i < axisCount; ++i) {
                float alpha = i * 2.f * M_PI / axisCount;
                float x = cosf(alpha);
                float y = sinf(alpha);

                tgt::vec2 direction(x, y);
                tgt::vec3 color(0.5 * (x + 1), 0.5 * (y + 1), 0);

                std::string label = fieldNames.at(i);

                starCoordinateAxes_.push_back(StarCoordinateAxis{ label, color, direction, 1.f, i });
            }
        }

        // update points
        extractPoints(list);

        std::vector<int> seriesIndices;
        for (const TimeSeries& series : list->getAllTimeSeries()) {
            std::string seriesPosition = std::to_string(series.getPosition().x) + " " + std::to_string(series.getPosition().y) + " " + std::to_string(series.getPosition().z);
            tgt::vec3 positionColor = series.getPosition() / list->getVolumeDimensions();
            renderedMembers_.addRow(seriesPosition, positionColor);
            selectedMembers_.addRow(seriesPosition, positionColor);
            subSelection_.insert(static_cast<int>(seriesIndices.size()));
            seriesIndices.push_back(static_cast<int>(seriesIndices.size()));
        }
        renderedMembers_.setSelectedRowIndices(seriesIndices);

        // Color Coding options
        colorCoding_.addOption("member", "Only Member", COLOR_MEMBER);
        colorCoding_.addOption("timeStep", "Only Time Step", COLOR_TIMESTEP);
        colorCoding_.addOption("memberAndTimeStep", "Member and Time Step", COLOR_MEMBER_AND_TIMESTEP);
        for (auto fieldname : list->getFieldNames()) {
            colorCoding_.addOption(fieldname, fieldname, COLOR_FIELD);
        }
        colorCoding_.addOption("slab", "Slabs", COLOR_SLAB);
        colorCoding_.addOption("plume", "Plumes", COLOR_PLUME);
        colorCoding_.addOption("slabsAndPlume", "Slabs and Plumes", COLOR_SLAB_AND_PLUME);
    }

    void StarCoordinates::createProjection() {
        // create the projection matrix based on the position of the axes
        projection_ = Eigen::MatrixXf(2, starCoordinateAxes_.size());
        for (size_t i = 0; i < starCoordinateAxes_.size(); ++i) {
            auto& axis = starCoordinateAxes_.at(i);
            float length = (axis.enabled_) ? axis.length_ : 0.f;
            tgt::vec2 coordinates = axis.direction_ * length;
            projection_(0, i) = coordinates.x;
            projection_(1, i) = coordinates.y;
        }
    }

    void StarCoordinates::extractPoints(const TimeSeriesList* timeseries) {
        // we assume that we have an equal amount of time steps for each field
        size_t numberOfPoints = timeseries->getNumberOfTimeSeries() * timeseries->getTimeSeries(0).getNumberOfTimeSteps();
        nDpoints_ = Eigen::MatrixXf(starCoordinateAxes_.size(), numberOfPoints);

        starCoordinatePoints_.clear();

        size_t j = 0;
        size_t seriesID = 0;
        for (auto& timeSeriesMember : timeseries->getAllTimeSeries()) {
            size_t timeStepID = 0;
            for (auto& timeSeriesStep : timeSeriesMember.getTimeSeriesSteps()) {
                size_t off = 0;
                for (size_t i = 0; i < starCoordinateAxes_.size(); ++i) {
                    int components = timeseries->getComponentForField(i);

                    // extract scalar value or calculate magnitude
                    float value = 0;
                    if (components == 1) {
                        value = timeSeriesStep.fieldValues_.at(off);
                    }
                    else {
                        // calculate magnitude
                        for (size_t k = 0; k < timeseries->getComponentForField(i); ++k) {
                            float v = timeSeriesStep.fieldValues_.at(off + k);
                            value += v * v;
                        }
                        value = sqrt(value);
                    }

                    // normalize value
                    tgt::vec2 range = timeseries->getRange(timeseries->getFieldName(i));
                    //value = ((range.y - range.x) < 1e-6) ? 0 : (value - range.x) / (range.y - range.x);
                    /*if (j == 0) {
                        LINFO("Range " + timeseries->getFieldName(i) + " " + std::to_string(range.x) + ", " + std::to_string(range.y));
                    }*/

                    nDpoints_(i, j) = value;

                    off += components; 
                }
                // create star coordinates point for the color and tool tip
                starCoordinatePoints_.push_back(StarCoordinatePoint(seriesID, timeStepID));
                /*if (j == 0) {
                    LINFO("IN " + std::to_string(seriesID) + "," + std::to_string(timeStepID));
                }*/

                ++j;
                ++timeStepID;
            }
            ++seriesID;
        }
    }

    void StarCoordinates::projectPoints() {
        projectedPoints_ = projection_ * nDpoints_;
    }

    void StarCoordinates::resetAxes() {
        size_t axisCount = starCoordinateAxes_.size();
        float radius = 0.8f;

        for (size_t i = 0; i < axisCount; ++i) {
            float alpha = i * 2.f * M_PI / axisCount;
            float x = cosf(alpha);
            float y = sinf(alpha);

            tgt::vec2 direction(x, y);

            starCoordinateAxes_.at(i).direction_ = direction;
            starCoordinateAxes_.at(i).length_ = 1.f;
        }
    }

    void StarCoordinates::resetSelection() {
        selectedMembers_.setSelectedRowIndices(std::vector<int>());
    }

    void StarCoordinates::mouseEvent(tgt::MouseEvent* e) {

        // We always accept the event and force a redraw.
        e->accept();
        invalidate();

        if (updateAxes)
            return;

        int x = e->x();
        int y = e->y();

        // Right click resets subselection to rendered members.
        int doubleClickAxis = -1;
        tgt::vec4 texel;
        bool selectedData = false;
        bool addToSelection = false;
        bool removeFromSelection = false;
        if (e->button() == tgt::MouseEvent::MOUSE_BUTTON_LEFT && e->action() == tgt::MouseEvent::DOUBLECLICK) {
            // handle margins
            RenderTarget* target = pickingBuffer_.getRenderTarget();
            tgt::ivec2 pixel;
            pixel = tgt::ivec2(
                mapRange(tgt::vec2(x, y), tgt::vec2::zero, tgt::vec2(outport_.getSize()),
                    tgt::vec2::zero, tgt::vec2(outport_.getSize())));
            tgt::ivec2 cursor = tgt::ivec2(pixel.x, outport_.getSize().y - pixel.y - 1);

            // select axis
            tgt::vec4 texel = target->getColorAtPos(cursor);
            if (almostEqual(texel.z, AXIS_ID, 1e-4)) {
                doubleClickAxis = std::round(texel.x * (starCoordinateAxes_.size() - 1));
            }
        }
        else if (e->button() == tgt::MouseEvent::MOUSE_BUTTON_LEFT && e->action() == tgt::MouseEvent::PRESSED
            && e->modifiers() == tgt::MouseEvent::MODIFIER_NONE) {
            
            // Reset last hit.
            lastHit_ = boost::none;
            
            // handle margins
            RenderTarget* target = pickingBuffer_.getRenderTarget();
            tgt::ivec2 pixel;
            pixel = tgt::ivec2(
                mapRange(tgt::vec2(x, y), tgt::vec2::zero, tgt::vec2(outport_.getSize()),
                    tgt::vec2::zero, tgt::vec2(outport_.getSize())));
            tgt::ivec2 cursor = tgt::ivec2(pixel.x, outport_.getSize().y - pixel.y - 1);

            // select axis
            texel = target->getColorAtPos(cursor);
            LINFO("texel " + std::to_string(texel.x) + ',' + std::to_string(texel.y) + ',' + std::to_string(texel.z));
            if (almostEqual(texel.z, AXIS_ID, 1e-4)) {
                selectedAxis_ = std::round(texel.x * (starCoordinateAxes_.size() - 1));
            }
            else if (almostEqual(texel.z, DATA_ID, 1e-4)) {
                selectedData = true;
            }
            else {
                selectedAxis_ = -1;
            }
        }
        else if (e->button() == tgt::MouseEvent::MOUSE_BUTTON_LEFT && e->action() == tgt::MouseEvent::RELEASED) {
            selectedAxis_ = -1;
        }
        else if (e->button() == tgt::MouseEvent::MOUSE_BUTTON_LEFT && e->action() == tgt::MouseEvent::PRESSED) {

            // Reset last hit.
            lastHit_ = boost::none;

            // handle margins
            RenderTarget* target = pickingBuffer_.getRenderTarget();
            tgt::ivec2 pixel;
            pixel = tgt::ivec2(
                    mapRange(tgt::vec2(x, y), tgt::vec2::zero, tgt::vec2(outport_.getSize()),
                             tgt::vec2::zero, tgt::vec2(outport_.getSize())));
            tgt::ivec2 cursor = tgt::ivec2(pixel.x, outport_.getSize().y - pixel.y - 1);

            // select axis
            texel = target->getColorAtPos(cursor);
            LINFO("texel " + std::to_string(texel.x) + ',' + std::to_string(texel.y) + ',' + std::to_string(texel.z));
            if (almostEqual(texel.z, DATA_ID, 1e-4)) {
                if(e->modifiers() == tgt::MouseEvent::CTRL)
                    addToSelection = true;
                else if(e->modifiers() == tgt::MouseEvent::ALT)
                    removeFromSelection = true;
                selectedData = true;
            }
        }

        // adapt axis position to mouse cursor position
        if (selectedAxis_ >= 0) {
            tgt::vec2 size = outport_.getSize();
            tgt::vec2 offset = size / 2.f;
            tgt::vec2 aspect = (size.x > size.y) ? tgt::vec2(size.x / size.y, 1) : tgt::vec2(1, size.y / size.x);

            tgt::vec2 rescaledDirection = ((2.f * tgt::vec2(x, y) / size) - tgt::vec2(1)) * aspect;
            rescaledDirection.y = -rescaledDirection.y;

            // get the selected axis
            auto& axis = starCoordinateAxes_.at(selectedAxis_);

            // split axis into direction and length
            axis.direction_ = tgt::normalize(rescaledDirection);
            axis.length_ = tgt::length(rescaledDirection / zoom_);


            // Calculate member index.  
            const TimeSeriesList* list = timeseriesInport_.getData();
            size_t numMembers = list->getNumberOfTimeSeries();
            int r = tgt::clamp<int>(std::round(texel.r * numMembers), 0, numMembers - 1);

            // Calculate time step index.
            size_t numTimeSteps = list->getTimeSeries(r).getNumberOfTimeSteps();
            int t = tgt::clamp<int>(std::round(texel.g * numTimeSteps), 0, numTimeSteps - 1);
            // Update last hit.
            lastHit_ = Hit{ x, y, r, t };

            updatedProjection = true;
        }

        if (selectedData) {
            // Calculate member index.
            const TimeSeriesList* list = timeseriesInport_.getData();
            size_t numMembers = list->getNumberOfTimeSeries();
            int r = tgt::clamp<int>(std::round(texel.r * numMembers), 0, numMembers - 1);
            if(addToSelection) {
                std::vector<int> selection = selectedMembers_.get();
                if(std::find(selection.begin(), selection.end(), r) == selection.end()) {
                    selection.push_back(r);
                    selectedMembers_.setSelectedRowIndices(selection);
                }
            }
            else if(removeFromSelection) {
                std::vector<int> selection = selectedMembers_.get();
                auto it = std::find(selection.begin(), selection.end(), r);
                if(it != selection.end()) {
                    selection.erase(it);
                    selectedMembers_.setSelectedRowIndices(selection);
                }
            }
            else {

                // Calculate time step index.
                size_t numTimeSteps = list->getTimeSeries(r).getNumberOfTimeSteps();
                int t = tgt::clamp<int>(std::round(texel.g * numTimeSteps), 0, numTimeSteps - 1);
                // Update last hit.
                lastHit_ = Hit{ x, y, r, t };

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
            }
        }

        // handle double click
        if (doubleClickAxis >= 0) {
            LINFO("collapse axis " + std::to_string(doubleClickAxis));
            auto& axis = starCoordinateAxes_.at(doubleClickAxis);

            // split axis into direction and length
            axis.enabled_ = !axis.enabled_;

            selectedAxis_ = -1;
            updateAxes = true;

            updatedProjection = true;
        }
    }

    void StarCoordinates::wheelEvent(tgt::MouseEvent* e) {
        float dz = getZoomFactor(mouseWheelZoomAcuteness_, ((e->button() & tgt::MouseEvent::MOUSE_WHEEL_UP)));
        zoom_ = zoom_ * dz;
        zoom_ = std::max(0.01f, zoom_);
        e->ignore();
    }

    void StarCoordinates::onEvent(tgt::Event* e) {
        tgt::MouseEvent* event = dynamic_cast<tgt::MouseEvent*>(e);
        if (event) {
            mouseEvent(event);

            if (event->action() == tgt::Event::WHEELEVENT && (event->modifiers() & event->CTRL)) {
                wheelEvent(event);
            }
        }

        RenderProcessor::onEvent(e);
    }

    float StarCoordinates::getZoomFactor(const float& acuteness, const bool& zoomIn) const {
        if (zoomIn) {
            return 1.f + 1.f / acuteness;
        }
        else {
            return 1.f - 1.f / acuteness;
        }
    }

} // namespace
