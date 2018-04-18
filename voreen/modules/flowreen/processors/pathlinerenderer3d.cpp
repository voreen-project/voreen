/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
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

#include "pathlinerenderer3d.h"

#include "voreen/core/interaction/camerainteractionhandler.h"
#include "floworthogonalslicerenderer.h"
#include "voreen/core/datastructures/volume/volumelist.h"
#include "modules/flowreen/utils/flowmath.h"
#include "modules/flowreen/datastructures/deprecated/volumeoperatorintensitymask.h"

#include <limits>

namespace voreen {

const std::string PathlineRenderer3D::loggerCat_("Processor.PathlineRenderer3D");

PathlineRenderer3D::PathlineRenderer3D()
    : RenderProcessor(),
    FlowreenProcessor(),
    seedingStrategyProp_(0),
    numSeedpointsProp_("numSeedpointsProp_", "# seedpoints:", 100, 1, 10000),
    seedOnXYSliceProp_("seedOnXYSliceProp_", "seed on XY slice:", true),
    seedOnXZSliceProp_("seedOnXZSliceProp_", "seed on XZ slice:", true),
    seedOnYZSliceProp_("seedOnYZSliceProp_", "seed on YZ slice:", true),
    gridSpacingProp_("gridSpacingProp_", "grid spacing: ", 8, 1, 1024),
    objectSizeProp_("objectSizeProp_", "size: ", 4, 1, 100),
    lineStyleProp_(0),
    thresholdingProp_(0),
    lineLengthProp_("lineLengthProp_", "min./max. line length:", tgt::vec2(0.0, 150.0f),
        tgt::vec2(0.0f), tgt::vec2(500.0f)),
    intensityProp_("intensityProp_", "min./max. intensity:", tgt::vec2(0, 100000.0f),
        tgt::vec2(0.0f), tgt::vec2(100000.0f)),
    lineColorProp_("lineColorProp_", "line color: ", tgt::vec4(0.875f, 0.0f, 0.0f, 1.0f)),
    segmentLengthProp_("lineSegmentLengthProp_", "max. segment length: ", 5, 1, 1000),
    integrationStepProp_("integrationStepProp", "integration stepwidth:", 0.25f, 0.1f, 1.0f),
    timestepProp_("timeStepProp", "timestep:", 0.0f, 0.0f, std::numeric_limits<float>::max()),
    camProp_( "camera", "Camera", tgt::Camera(tgt::vec3(0.0f, 0.0f, 3.5f),
        tgt::vec3(0.0f, 0.0f, 0.0f), tgt::vec3(0.0f, 1.0f, 0.0f)) ),
    cameraHandler_(0),
    flows_(),
    contextCollection_(0),
    flowCollection_(0),
    intensityMasks_(),
    shader_(0),
    currentStyle_(STYLE_LINES),
    seedingStrategy_(SEED_SLICES_GRID),
    thresholding_(THRESHOLDING_LINELENGTH),
    numPathlines_(0),
    pathlines_(0),
    currentTimestep_(0),
    previousTimestep_(0),
    slicePositions_(1, 1, 1),
    coInport_(Port::INPORT, "coprocessor.slicepositions"),
    imgOutport_(Port::OUTPORT, "image.pathlines", "Pathlines Image Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER),
    inportContext_(Port::INPORT, "volumecollection.context", "VolumeCollection Context"),
    inportFlows_(Port::INPORT, "volumecollection.flow", "VolumeCollection Flow")
{
    lineStyleProp_ = new OptionProperty<LineStyle>("lineStyleProp", "pathline style:");
    lineStyleProp_->addOption("points", "points", STYLE_POINTS);
    lineStyleProp_->addOption("line segments", "line segments", STYLE_SEGMENTS);
    lineStyleProp_->addOption("lines", "lines", STYLE_LINES);
    lineStyleProp_->addOption("tubes", "tubes", STYLE_TUBES);
    lineStyleProp_->addOption("arrows", "arrows", STYLE_ARROWS);
    lineStyleProp_->onChange(MemberFunctionCallback<PathlineRenderer3D>(this,
        &PathlineRenderer3D::onLineStyleChange));

    MemberFunctionCallback<PathlineRenderer3D> cmaSeeding(this,
        &PathlineRenderer3D::onSeedingStrategyChange);
    MemberFunctionCallback<PathlineRenderer3D> cmaTimestep(this,
        &PathlineRenderer3D::onTimestepChange);

    seedingStrategyProp_ = new OptionProperty<SeedingStrategy>("seedingStrategyProp_",
        "seeding strategy:");
    seedingStrategyProp_->addOption("random position", "random position", SEED_RANDOM);
    seedingStrategyProp_->addOption("regular grid", "regular grid", SEED_GRID);
    seedingStrategyProp_->addOption("regular grid on slices", "regular grid on slices", SEED_SLICES_GRID);
    seedingStrategyProp_->onChange(cmaSeeding);

    integrationStepProp_.onChange(cmaSeeding);
    integrationStepProp_.onChange(cmaTimestep);
    integrationStepProp_.onChange(MemberFunctionCallback<PathlineRenderer3D>(this,
        &PathlineRenderer3D::adjustTimestepProperty));
    timestepProp_.onChange(cmaTimestep);

    thresholdingProp_ = new OptionProperty<Thresholding>("thresholdingProp_", "thresholding:");
    thresholdingProp_->addOption("none", "none", THRESHOLDING_NONE);
    thresholdingProp_->addOption("line length", "line length", THRESHOLDING_LINELENGTH);
    thresholdingProp_->addOption("intensity", "intensity", THRESHOLDING_INTENSITY);
    thresholdingProp_->addOption("line length OR intensity", "line length OR intensity", THRESHOLDING_OR);
    thresholdingProp_->addOption("line length AND intensity", "line length AND intensity", THRESHOLDING_AND);
    thresholdingProp_->onChange(MemberFunctionCallback<PathlineRenderer3D>(this,
        &PathlineRenderer3D::onThresholdingChange));

    MemberFunctionCallback<PathlineRenderer3D> cmaClearLines(this, &PathlineRenderer3D::clearPathlines);
    lineLengthProp_.onChange(cmaClearLines);
    intensityProp_.onChange(cmaClearLines);
    intensityProp_.onChange(MemberFunctionCallback<PathlineRenderer3D>(this,
        &PathlineRenderer3D::onIntensityChange));
    integrationStepProp_.setStepping(0.1f);

    addProperty(useBoundingBoxProp_);
    addProperty(boundingBoxColorProp_);
    addProperty(useCoordinateAxisProp_);
    addProperty(thresholdingProp_);
    addProperty(lineLengthProp_);
    addProperty(intensityProp_);
    addProperty(lineStyleProp_);
    addProperty(segmentLengthProp_);
    addProperty(integrationStepProp_);
    addProperty(lineColorProp_);
    addProperty(objectSizeProp_);
    addProperty(seedingStrategyProp_);
    addProperty(numSeedpointsProp_);
    addProperty(gridSpacingProp_);
    addProperty(seedOnXYSliceProp_);
    addProperty(seedOnXZSliceProp_);
    addProperty(seedOnYZSliceProp_);
    addProperty(timestepProp_);
    addProperty(camProp_);

    cameraHandler_ = new CameraInteractionHandler("cameraHandler", "Camera Handler", &camProp_);
    addInteractionHandler(cameraHandler_);

    addPort(inportFlows_);
    addPort(inportContext_);
    addPort(imgOutport_);
    addPort(coInport_);

    onLineStyleChange();
    onSeedingStrategyChange();
    onThresholdingChange();
}

PathlineRenderer3D::~PathlineRenderer3D() {
    for (size_t i = 0; i < intensityMasks_.size(); ++i)
        delete intensityMasks_[i];
    intensityMasks_.clear();
    flows_.clear();
    delete [] pathlines_;
    delete lineStyleProp_;
    delete seedingStrategyProp_;
    delete thresholdingProp_;
    if (shader_ != 0) {
        shader_->deactivate();
        ShdrMgr.dispose(shader_);
    }
    delete cameraHandler_;
}

void PathlineRenderer3D::process() {
    imgOutport_.activateTarget();

    // Check the input VolumeSet: it must contain a VolumeSeries of modality 'flow'
    // and may of course not be NULL.
    //
    bool inputChanged = (inportContext_.hasChanged() || inportFlows_.hasChanged());
    contextCollection_ = inportContext_.getData();
    flowCollection_ = inportFlows_.getData();

    // Check whether the number of flows has changed and whether are all of the same
    // dimensions and eventually enforce an update of the flow data.
    //
    if (checkFlows(inputChanged) < 1) {
        LERROR("process(): no valid flows for input! Cannot continue.");
        return;
    }

    const size_t gridSpacing = static_cast<size_t>(gridSpacingProp_.get());
    switch (seedingStrategy_) {
        case SEED_RANDOM:
            initPathlines(static_cast<size_t>(numSeedpointsProp_.get()));
            break;

        case SEED_GRID:
            initPathlinesGrid(gridSpacing);
            break;

        case SEED_SLICES_GRID:
            {
                FlowOrthogonalSliceRenderer* osr = coInport_.getConnectedProcessor();
                if (osr != 0) {
                    tgt::ivec3 positions = tgt::clamp(osr->getSlicePositions(),
                        tgt::ivec3(-1), tgt::ivec3(flowDimensions_)) - tgt::ivec3(1);
                    if (positions != slicePositions_) {
                        slicePositions_ = positions;
                        delete [] pathlines_;
                        pathlines_ = 0;
                    }
                    initPathlinesSliceGrid(gridSpacing);
                } else {
                    LINFO("There seems to be no co-processor attached providing slices.");
                    LINFO("Switching seeding to regular grid...");
                    slicePositions_ = tgt::ivec3(-1);
                    seedingStrategyProp_->selectByValue(SEED_SLICES_GRID);
                    onSeedingStrategyChange();
                    seedingStrategy_ = SEED_GRID;
                }
            }
            break;
        default:
            break;
    }   // switch

    if (pathlines_ == 0)
        return;

    // important: save current camera state before using the processor's camera or
    // successive processors will use those settings!
    //
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();

    camProp_.look(imgOutport_.getSize());

    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

    if (useCoordinateAxisProp_.get() == true)
        renderCoordinateAxis();
    if (useBoundingBoxProp_.get() == true)
        renderBoundingBox(boundingBoxColorProp_.get());

    switch (currentStyle_) {
        case STYLE_POINTS:
            renderAsPoints(lineColorProp_.get());
            break;
        case STYLE_ARROWS:
            renderAsArrows(lineColorProp_.get());
            break;
        case STYLE_TUBES:
            renderAsTubes(lineColorProp_.get());
            break;
        case STYLE_LINES:
            renderAsLines(lineColorProp_.get());
            break;
        case STYLE_SEGMENTS:
            renderAsLineSegments(lineColorProp_.get());
            break;
    }

    if (shader_ != 0)
        shader_->deactivate();
    glActiveTexture(GL_TEXTURE0);
    imgOutport_.deactivateTarget();

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();
}

void PathlineRenderer3D::initialize() {
    FlowreenProcessor::init();
    RenderProcessor::initialize();
    //adjustTimestepProperty();
}

void PathlineRenderer3D::invalidate(int inv) {
    Processor::invalidate(inv);
    adjustTimestepProperty();
}

// private methods
//

bool PathlineRenderer3D::applyThresholds(std::vector<tgt::vec3>& pathline, const float& length) const {
    if (pathline.empty() == true)
        return false;

    const tgt::vec2& lineLength = lineLengthProp_.get();

    switch (thresholding_) {
        case THRESHOLDING_NONE:
            return true;

        case THRESHOLDING_LINELENGTH:
            if ((length < lineLength.x) || (length > lineLength.y)) {
                pathline.clear();
                return false;
            }
            break;

        case THRESHOLDING_AND:
            if ((length < lineLength.x) || (length > lineLength.y)) {
                pathline.clear();
                return false;
            }
            // no break here, go and check intensity, becaus the length check is passed.
        case THRESHOLDING_INTENSITY:
            if (intensityMasks_.empty() == false) {
                tgt::ivec3 r0 = static_cast<tgt::ivec3>(pathline[0]);
                bool keep = false;
                for (size_t i = 0; i < intensityMasks_.size(); ++i) {
                    const bool& value = static_cast<const VolumeAtomic<bool>*>(intensityMasks_[i]->getRepresentation<VolumeRAM>())->voxel(r0);
                    if (value == true) {
                        keep = true;
                        break;
                    }
                }

                if (keep == false)
                    pathline.clear();
                return keep;
            }
            break;

        case THRESHOLDING_OR:
            if ((length < lineLength.x) || (length > lineLength.y)) {
                // First check (line length) has failed, so give a second
                // chance on intensity check.
                //
                bool keep = false;
                if (intensityMasks_.empty() == false) {
                    tgt::ivec3 r0 = static_cast<tgt::ivec3>(pathline[0]);
                    for (size_t i = 0; i < intensityMasks_.size(); ++i) {
                        const bool& value = static_cast<const VolumeAtomic<bool>*>(intensityMasks_[i]->getRepresentation<VolumeRAM>())->voxel(r0);
                        if (value == true) {
                            keep = true;
                            break;
                        }
                    }   // for
                }

                if (keep == false)
                    pathline.clear();
            }
            break;
    }   // switch
    return (! pathline.empty());
}

size_t PathlineRenderer3D::markContextData(const tgt::vec2& thresholds) {
    for (size_t i = 0; i < intensityMasks_.size(); ++i)
        delete intensityMasks_[i];
    intensityMasks_.clear();

    if ((contextCollection_ == 0) || (flowDimensions_ == tgt::svec3::zero))
        return 0;

    for (size_t i = 0; i < contextCollection_->size(); ++i) {
        if (contextCollection_->at(i) == 0)
            continue;

        const VolumeBase* volume = contextCollection_->at(i);
        if ((volume != 0) && (volume->getDimensions() == flowDimensions_))
            intensityMasks_.push_back(VolumeOperatorIntensityMask::APPLY_OP(volume, thresholds));
    }

    return intensityMasks_.size();
}

size_t PathlineRenderer3D::checkFlows(const bool forceUpdate) {
    if ((forceUpdate == false) && (flows_.empty() == false))
        return flows_.size();

    flowDimensions_ = tgt::svec3::zero;
    flows_.clear();
    if ((flowCollection_ == 0) || (flowCollection_->size() < 1))
        return 0;

    for (size_t i = 0; i < flowCollection_->size(); ++i) {
        if (flowCollection_->at(i) != 0) {
            const VolumeFlow3D* flowVolume = dynamic_cast<const VolumeFlow3D*>(flowCollection_->at(i)->getRepresentation<VolumeRAM>());
            if (flowVolume == 0)
                continue;

            if (flowDimensions_ == tgt::svec3::zero) {
                flowDimensions_ = flowVolume->getDimensions();
                flows_.push_back(&(flowVolume->getFlow3D()));
            } else if (flowDimensions_ == flowVolume->getDimensions())
                flows_.push_back(&(flowVolume->getFlow3D()));
            else
                LWARNING("checkFlows(): dimensions of flow volume #" << i << " mismatch!");
        }
    }

    if (flows_.empty() == true)
        flowDimensions_ = tgt::ivec3::zero;

    adjustTimestepProperty();
    return flows_.size();
}

void PathlineRenderer3D::initPathlines(const size_t numPoints)
{
    if ((numPoints == numPathlines_) && (pathlines_ != 0))
        return;

    delete [] pathlines_;
    pathlines_ = 0;

    numPathlines_ = numPoints;
    if (numPathlines_ == 0)
        return;

    const float deltaT = integrationStepProp_.get();
    tgt::vec3 dim = static_cast<tgt::vec3>(flowDimensions_);
    pathlines_ = new std::vector<tgt::vec3>[numPathlines_];
    float length = 0.0f;
    for (size_t i = 0; i < numPathlines_; ++i) {
        tgt::vec3 r0 = FlowMath::uniformRandomVec3() * dim;
        pathlines_[i] = FlowMath::computePathline(flows_, r0, deltaT, &length);
        applyThresholds(pathlines_[i], length);
    }
}

void PathlineRenderer3D::initPathlinesGrid(const size_t spacing)
{
    tgt::svec3 grid = ((flowDimensions_ - tgt::svec3(1)) / spacing) + tgt::svec3(1);
    size_t numPathlines = (grid.x * grid.y * grid.z);
    if ((numPathlines == numPathlines_) && (pathlines_ != 0))
        return;

    delete [] pathlines_;
    pathlines_ = 0;
    numPathlines_ = numPathlines;
    if (numPathlines_ == 0)
        return;

    const float deltaT = integrationStepProp_.get();
    pathlines_ = new std::vector<tgt::vec3>[numPathlines_];
    for (size_t z = 0; z < grid.z; ++z) {
        float fz = static_cast<float>(z * spacing);
        for (size_t y = 0; y < grid.y; ++y) {
            float fy = static_cast<float>(y * spacing);
            for (size_t x = 0; x < grid.x; ++x) {
                float fx = static_cast<float>(x * spacing);
                tgt::vec3 pos(fx, fy, fz);
                size_t n = z * (grid.x * grid.y) + y * (grid.x) + x;
                float length = 0.0f;
                pathlines_[n] = FlowMath::computePathline(flows_, pos, deltaT, &length);
                applyThresholds(pathlines_[n], length);
            }   // for (x
        }   // for (y
    }   // for (z
}

void PathlineRenderer3D::initPathlinesSliceGrid(const size_t spacing)
{
    tgt::svec3 grid = ((flowDimensions_ - tgt::svec3(1)) / spacing) + tgt::svec3(1);
    size_t numPathlines = 0;
    if ((slicePositions_.z >= 0) && (seedOnXYSliceProp_.get() == true))
        numPathlines += (grid.x * grid.y);
    if ((slicePositions_.y >= 0) && (seedOnXZSliceProp_.get() == true))
        numPathlines += (grid.x * grid.z);
    if ((slicePositions_.x >= 0) && (seedOnYZSliceProp_.get() == true))
        numPathlines += (grid.y * grid.z);

    if ((numPathlines == numPathlines_) && (pathlines_ != 0))
        return;

    delete [] pathlines_;
    pathlines_ = 0;
    numPathlines_ = numPathlines;
    if (numPathlines_ == 0)
        return;

    pathlines_ = new std::vector<tgt::vec3>[numPathlines_];

    const float deltaT = integrationStepProp_.get();
    float fx = 0.0f, fy = 0.0f, fz = 0.0f;
    size_t n = 0;
    if ((slicePositions_.x >= 0) && (seedOnYZSliceProp_.get() == true)) {
        fx = static_cast<float>(slicePositions_.x);
        for (size_t z = 0; z < grid.z; ++z) {
            fz = static_cast<float>(z * spacing);
            for (size_t y = 0; y < grid.y; ++y, ++n) {
                fy = static_cast<float>(y * spacing);
                tgt::vec3 pos(fx, fy, fz);
                float length = 0.0f;
                pathlines_[n] = FlowMath::computePathline(flows_, pos, deltaT, &length);
                applyThresholds(pathlines_[n], length);
            }   // for (y
        }   // for (z
    }

    if ((slicePositions_.y >= 0) && (seedOnXZSliceProp_.get() == true)) {
        fy = static_cast<float>(slicePositions_.y);
        for (size_t z = 0; z < grid.z; ++z) {
            fz = static_cast<float>(z * spacing);
            for (size_t x = 0; x < grid.x; ++x, ++n) {
                fx = static_cast<float>(x * spacing);
                tgt::vec3 pos(fx, fy, fz);
                float length = 0.0f;
                pathlines_[n] = FlowMath::computePathline(flows_, pos, deltaT, &length);
                applyThresholds(pathlines_[n], length);
            }   // for (x
        }   // for (z
    }

    if ((slicePositions_.z >= 0) && (seedOnXYSliceProp_.get() == true)) {
        fz = static_cast<float>(slicePositions_.z);
        for (size_t y = 0; y < grid.y; ++y) {
            fy = static_cast<float>(y * spacing);
            for (size_t x = 0; x < grid.x; ++x, ++n) {
                fx = static_cast<float>(x * spacing);
                tgt::vec3 pos(fx, fy, fz);
                float length = 0.0f;
                pathlines_[n] = FlowMath::computePathline(flows_, pos, deltaT, &length);
                applyThresholds(pathlines_[n], length);
            }   // for (x
        }   // for (y
    }
}

void PathlineRenderer3D::adjustTimestepProperty() {
    const VolumeList* collection = inportFlows_.getData();
    if (flows_.empty() == false) {
        float stepping = integrationStepProp_.get();
        float maxTimestep = static_cast<float>(flows_.size()) - stepping;
        if (timestepProp_.get() > maxTimestep)
            timestepProp_.set(maxTimestep);
        timestepProp_.setMaxValue(maxTimestep);
        timestepProp_.setStepping(stepping);
        timestepProp_.setReadOnlyFlag(false);
    } else if (collection != 0) {
        // Only adjust the property's value, if the flows are empty,
        // VolumeCollection containing flows is present and the
        // collection is of size 1 or less, because when deserializing
        // the timestep property, all these conditions are met, but
        // the correct value for the property is already set and it
        // will become valid, as soon the the entire network is ready.
        //
        /*if (collection->size() <= 1)
            timestepProp_.set(0.0f);
        timestepProp_.setMaxValue(0.0f); */
        timestepProp_.setReadOnlyFlag(true);
    }
    timestepProp_.updateWidgets();
}

void PathlineRenderer3D::clearPathlines() {
    delete [] pathlines_;
    pathlines_ = 0;
    numPathlines_ = 0;
}

void PathlineRenderer3D::onIntensityChange() {
    markContextData(intensityProp_.get());
}

void PathlineRenderer3D::onLineStyleChange() {
    currentStyle_ = lineStyleProp_->getValue();
    if (currentStyle_ != STYLE_SEGMENTS)
        segmentLengthProp_.setVisibleFlag(false);
    else
        segmentLengthProp_.setVisibleFlag(true);
}

void PathlineRenderer3D::onSeedingStrategyChange() {
    seedingStrategy_ = seedingStrategyProp_->getValue();

    if (seedingStrategy_ == SEED_RANDOM) {
        gridSpacingProp_.setVisibleFlag(false);
        numSeedpointsProp_.setVisibleFlag(true);
    } else {
        gridSpacingProp_.setVisibleFlag(true);
        numSeedpointsProp_.setVisibleFlag(false);
    }

    if (seedingStrategy_ == SEED_SLICES_GRID) {
        seedOnXYSliceProp_.setVisibleFlag(true);
        seedOnXZSliceProp_.setVisibleFlag(true);
        seedOnYZSliceProp_.setVisibleFlag(true);
    } else {
        seedOnXYSliceProp_.setVisibleFlag(false);
        seedOnXZSliceProp_.setVisibleFlag(false);
        seedOnYZSliceProp_.setVisibleFlag(false);
    }

    clearPathlines();
}

void PathlineRenderer3D::onTimestepChange() {
    previousTimestep_ = currentTimestep_;
    currentTimestep_ = static_cast<size_t>(timestepProp_.get() / integrationStepProp_.get());
}

void PathlineRenderer3D::onThresholdingChange() {
    thresholding_ = thresholdingProp_->getValue();
    switch (thresholding_) {
        default:
        case THRESHOLDING_NONE:
            lineLengthProp_.setVisibleFlag(false);
            intensityProp_.setVisibleFlag(false);
            break;
        case THRESHOLDING_LINELENGTH:
            lineLengthProp_.setVisibleFlag(true);
            intensityProp_.setVisibleFlag(false);
            break;
        case THRESHOLDING_INTENSITY:
            lineLengthProp_.setVisibleFlag(false);
            intensityProp_.setVisibleFlag(true);
            break;
        case THRESHOLDING_OR:
        case THRESHOLDING_AND:
            lineLengthProp_.setVisibleFlag(true);
            intensityProp_.setVisibleFlag(true);
            break;
    }
    clearPathlines();
}

void PathlineRenderer3D::renderAsArrows(const tgt::vec4& lineColor) {
    if (setupShader() == true) {
        shader_->setUniform("camPos_", camProp_.get().getPosition());
        shader_->setUniform("lightDir_", camProp_.get().getLook());
    }

    const float arrowSize = static_cast<float>(objectSizeProp_.get());

    glEnable(GL_DEPTH_TEST);
    glColor4fv(lineColor.elem);
    for (size_t i = 0; i < numPathlines_; ++i) {
        const std::vector<tgt::vec3>& pathline = pathlines_[i];
        if (pathline.empty() == true)
            continue;

        float length = 0.0f;
        size_t a = (currentTimestep_ < pathline.size()) ? currentTimestep_ : (pathline.size() - 1);
        for (size_t j = 0; j <= a; ++j) {
            if (j > 0)
                length += tgt::length(pathline[j] - pathline[j - 1]);
        }
        length /= arrowSize;
        length += 0.5f;

        tgt::mat4 trafo = FlowMath::getTransformationMatrix(pathline, a,
            ((length < 1.0f) ? (arrowSize * length) : arrowSize));
        if (previousTimestep_ > currentTimestep_) {
            trafo.t00 *= -1.0f; trafo.t02 *= -1.0f;
            trafo.t10 *= -1.0f; trafo.t12 *= -1.0f;
            trafo.t20 *= -1.0f; trafo.t22 *= -1.0f;
        }
        renderArrow(trafo, length);
    }
    glColor4f(0.0f, 0.0f, 0.0f, 0.0f);
}

void PathlineRenderer3D::renderAsLines(const tgt::vec4& lineColor) {
    if (setupShader() == true) {
        shader_->setUniform("camPos_", camProp_.get().getPosition());
        shader_->setUniform("lightDir_", camProp_.get().getLook());
    }

    const float deltaT = integrationStepProp_.get();
    const float lineWidth = static_cast<float>(objectSizeProp_.get());
    glEnable(GL_DEPTH_TEST);
    glLineWidth(lineWidth);
    glColor4fv(lineColor.elem);
    for (size_t i = 0; i < numPathlines_; ++i) {
        const std::vector<tgt::vec3>& pathline = pathlines_[i];
        if (pathline.size() < 1)
            continue;

        size_t b = (currentTimestep_ < pathline.size()) ? currentTimestep_ : (pathline.size() - 1);
        glBegin(GL_LINE_STRIP);
        for (size_t j = 0; j <= b; ++j) {
            tgt::vec3 r0 = mapToFlowBoundingBox(pathline[j]);
            tgt::vec3 tangent;
            if (j < b) {
                tgt::vec3 r1 = mapToFlowBoundingBox(pathline[j + 1]);
                tangent = FlowMath::normalize((r1 - r0) / deltaT);
            } else if (j > 0) {
                tgt::vec3 r1 = mapToFlowBoundingBox(pathline[j - 1]);
                tangent = FlowMath::normalize((r0 - r1) / deltaT);
            } else
                tangent = FlowMath::normalize(r0);

            glNormal3fv(tangent.elem);
            glVertex3fv(r0.elem);
        }
        glEnd();
    }   // for (i
    glLineWidth(1.0f);
}

void PathlineRenderer3D::renderAsLineSegments(const tgt::vec4& lineColor) {
    if (setupShader() == true) {
        shader_->setUniform("camPos_", camProp_.get().getPosition());
        shader_->setUniform("lightDir_", camProp_.get().getLook());
    }

    const float deltaT = integrationStepProp_.get();
    const float lineWidth = static_cast<float>(objectSizeProp_.get());
    const int segmentLength = segmentLengthProp_.get();

    glEnable(GL_DEPTH_TEST);
    glLineWidth(lineWidth);
    for (size_t i = 0; i < numPathlines_; ++i) {
        const std::vector<tgt::vec3>& pathline = pathlines_[i];
        if (pathline.size() < 1)
            continue;

        int size = static_cast<int>(pathline.size() - 1);
        int b = static_cast<int>(currentTimestep_) + 1;
        if (b > size)
            b = size;

        int a = (b - segmentLength);
        if (a < 0)
            a = 0;

        glBegin(GL_LINE_STRIP);
        for (int j = a; j <= b; ++j) {
            tgt::vec3 r0 = mapToFlowBoundingBox(pathline[j]);
            tgt::vec3 tangent;
            if (j < b) {
                tgt::vec3 r1 = mapToFlowBoundingBox(pathline[j + 1]);
                tangent = FlowMath::normalize((r1 - r0) / deltaT);
            } else if (j > 0) {
                tgt::vec3 r1 = mapToFlowBoundingBox(pathline[j - 1]);
                tangent = FlowMath::normalize((r0 - r1) / deltaT);
            } else
                tangent = FlowMath::normalize(r0);

            float alpha = 1.0f;
            if (b > a) {
                alpha = (j - a) / static_cast<float>(b - a);
                if (previousTimestep_ > currentTimestep_)
                    alpha = 1.0f - alpha;
            }

            glColor4f(lineColor.r, lineColor.g, lineColor.b, alpha);
            glNormal3fv(tangent.elem);
            glVertex3fv(r0.elem);
        }
        glEnd();
    }   // for (i
    glLineWidth(1.0f);
}

void PathlineRenderer3D::renderAsPoints(const tgt::vec4& lineColor) {
    const float pointSize = static_cast<float>(objectSizeProp_.get());

    if (shader_ != 0)
        shader_->deactivate();

    glEnable(GL_DEPTH_TEST);
    glPointSize(pointSize);
    glBegin(GL_POINTS);

    for (size_t i = 0; i < numPathlines_; ++i) {
        const std::vector<tgt::vec3>& pathline = pathlines_[i];
        if (pathline.size() < 1)
            continue;

        if (previousTimestep_ <= currentTimestep_) {
            size_t b = (currentTimestep_ < pathline.size()) ? currentTimestep_ : (pathline.size() - 1);
            for (size_t j = 0; j <= b; ++j) {
                float alpha = j / static_cast<float>(b);
                glColor4f(lineColor.r, lineColor.g, lineColor.b, alpha);
                glVertex3fv(mapToFlowBoundingBox(pathline[j]).elem);
            }
        } else {
            int a = static_cast<int>(pathline.size() - 1);
            for (int j = a; j >= static_cast<int>(currentTimestep_); --j) {
                float alpha = (a - j) / static_cast<float>(a - currentTimestep_);
                glColor4f(lineColor.r, lineColor.g, lineColor.b, alpha);
                glVertex3fv(mapToFlowBoundingBox(pathline[j]).elem);
            }
        }
    }
    glEnd();
    glColor4f(0.0f, 0.0f, 0.0f, 0.0f);
    glPointSize(1.0f);
}

void PathlineRenderer3D::renderAsTubes(const tgt::vec4& lineColor) {
    if (setupShader() == true) {
        shader_->setUniform("camPos_", camProp_.get().getPosition());
        shader_->setUniform("lightDir_", camProp_.get().getLook());
    }

    const size_t stepwidth = 1;
    const float tubeRadius = (objectSizeProp_.get() / 10.0f);

    glEnable(GL_DEPTH_TEST);
    glColor4fv(lineColor.elem);
    for (size_t i = 0; i < numPathlines_; ++i) {
        const std::vector<tgt::vec3>& pathline = pathlines_[i];
        if (pathline.size() < 1)
            continue;

        bool forwardDirection = true;
        switch ((currentTimestep_ / pathline.size()) % 2) {
            case 0:
                break;
            case 1:
                forwardDirection = false;
                break;
        }

        tgt::mat4 trans = FlowMath::getTransformationMatrix(pathline, 0, tubeRadius);
        std::vector<tgt::vec3> circleInit = getTransformedCircle(trans);
        std::vector<tgt::vec3>& circle1 = circleInit;

        size_t mod = (currentTimestep_ % pathline.size());
        size_t b = 0;

        if (forwardDirection == true) {
            b = mod;
            if (b >= pathline.size())
                b = (pathline.size() - 1);
        } else {
            int aux = (static_cast<int>(pathline.size()) - 1) - static_cast<int>(mod);
            b = (aux < 0) ? 0 : static_cast<size_t>(aux);
        }

        for (size_t n = 1; n <= b; n += stepwidth) {
            trans = FlowMath::getTransformationMatrix(pathline, n, tubeRadius);
            std::vector<tgt::vec3> circle2 = getTransformedCircle(trans);

            glBegin(GL_QUAD_STRIP);
            for (size_t k = 0; k <= circle1.size(); ++k) {
                const tgt::vec3& r1 = circle1[k % circle1.size()];
                const tgt::vec3& r2 = circle2[k % circle2.size()];
                const tgt::vec3& r4 = circle1[(k + 1) % circle1.size()];

                tgt::vec3 normal = FlowMath::normalize(tgt::cross((r2 - r1), (r4 - r1)));
                glNormal3fv(normal.elem);
                glVertex3fv(mapToFlowBoundingBox(r1).elem);
                glVertex3fv(mapToFlowBoundingBox(r2).elem);
            }   // for (k
            glEnd();

            circle1.clear();
            circle1 = circle2;
        }   // for (n
        circle1.clear();
    }   // for (i
    glColor4f(0.f, 0.f, 0.f, 0.f);
}

bool PathlineRenderer3D::setupShader() {
    if (shader_ == 0) {
        try {
            shader_ = ShdrMgr.load("phong", generateHeader(), false);
        } catch(tgt::Exception) {
            LGL_ERROR;
            return false;
        }
    }

    // activate the shader and set the needed uniforms
    if (shader_->isActivated() == false)
        shader_->activate();

    shader_->setUniform("usePhongShading_", true);
    shader_->setUniform("ka_", 0.3f);
    shader_->setUniform("kd_", 0.7f);
    shader_->setUniform("ks_", 1.0f);
    shader_->setUniform("shininess_", 32.0f);

    if ((currentStyle_ == STYLE_LINES) || (currentStyle_ == STYLE_SEGMENTS))
        shader_->setUniform("useForLines_", true);
    else
        shader_->setUniform("useForLines_", false);

    return true;
}

}   // namespace
