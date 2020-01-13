/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "explosionproxygeometry.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"
#include "voreen/core/datastructures/geometry/geometrysequence.h"
#include "tgt/glmath.h"
#include <vector>

namespace voreen {

using tgt::vec3;

ExplosionProxyGeometry::ExplosionProxyGeometry()
    : RenderProcessor()
    , canvasSize_("canvasSize", "Canvas Size", tgt::ivec2(256), tgt::ivec2(64), tgt::ivec2(10000), Processor::VALID)
    , brickingMode_("BrickingMode", "Bricking mode", Processor::INVALID_PROGRAM)
    , showBoundingBoxes_("showBoundingBoxes", "Show bounding boxes", Processor::INVALID_PROGRAM)
    , brickColor_("brickColor", "Brick color", tgt::vec4(0.0f, 0.0f, 0.0f, 1.0f))
    , selectedBrickColor_("selectedBrickColor", "Selected brick color", tgt::vec4(1.0f, 1.0f, 0.0f, 1.0f))
    , hideNonSelectedBricks_("hideNonSelectedBricks", "Hide non-selected bricks")
    , hideSelectedBricks_("hideSelectedBricks", "Hide selected bricks")
    , showAllBricks_("showAllBricks", "Show all bricks")
    , numBricksXp_("numberBricksX", "Number of bricks on X", 1, 1, 10)
    , numBricksYp_("numberBricksY", "Number of bricks on Y", 1, 1, 10)
    , numBricksZp_("numberBricksZ", "Number of bricks on Z", 1, 1, 10)
    , regularExplosionGap_("regularBrickExplosionGap", "Explosion-gap between bricks", 0.f, 0.f, 10.f)
    , manualExplosionGap_("manualBrickExplosionGap", "Explosion-gap between bricks", 0.f, 0.f, 10.f)
    , XbrickingPlane_("XbrickingPlane", "X-brickingplane", 0, 0, 100000)
    , addXbrickingPlane_("addXbrickingPlane", "Add X-brickingplane")
    , deleteLastXbrickingPlane_("deleteLastXbrickingPlane", "Delete last X-brickingplane")
    , clearXbrickingList_("clearXbrickingList", "Clear X-bricking list")
    , YbrickingPlane_("YbrickingPlane", "Y-brickingplane", 0, 0, 100000)
    , addYbrickingPlane_("addYbrickingPlane", "Add Y-brickingplane")
    , deleteLastYbrickingPlane_("deleteLastYbrickingPlane", "Delete last Y-brickingplane")
    , clearYbrickingList_("clearYbrickingList", "Clear Y-bricking list")
    , ZbrickingPlane_("ZbrickingPlane", "Z-brickingplane", 0, 0, 100000)
    , addZbrickingPlane_("addZbrickingPlane", "Add Z-brickingplane")
    , deleteLastZbrickingPlane_("deleteLastZbrickingPlane", "Delete last Z-brickingplane")
    , clearZbrickingList_("clearZbrickingList", "Clear Z-bricking list")
    , resetMode_("resetMode", "Reset mode")
    , undo_("Undo", "Undo")
    , takeForCustom_("takeForCustomMode", "Take this for custom mode")
    , translation_("translation", "Translation", tgt::vec3(0.f, 0.f, 0.f), tgt::vec3(-10.f), tgt::vec3(10.f))
    , camera_("camera", "Camera", tgt::Camera(), Processor::VALID)
    , selectingOneBrickEvent_("mouseEvent.selectingOneBrick", "Select one brick", this,
                               &ExplosionProxyGeometry::onSelectingOneBrickEvent,
                               tgt::MouseEvent::MOUSE_BUTTON_LEFT, tgt::MouseEvent::PRESSED, tgt::Event::CTRL)
    , selectingSeveralBricksEvent_("mouseEvent.selectingSeveralBricks", "Select several bricks", this,
                                    &ExplosionProxyGeometry::onSelectingSeveralBricksEvent,
                                    tgt::MouseEvent::MOUSE_BUTTON_LEFT, tgt::MouseEvent::PRESSED, tgt::Event::ALT)
    , translateSelectedBricksEvent_("mouseEvent.translateSelectedBricks", "Translate selected bricks", this,
                                     &ExplosionProxyGeometry::onTranslateSelectedBricksEvent,
                                     tgt::MouseEvent::MOUSE_BUTTON_LEFT,
                                     tgt::MouseEvent::PRESSED | tgt::MouseEvent::MOTION | tgt::MouseEvent::RELEASED,
                                     tgt::Event::Modifier(tgt::Event::CTRL | tgt::Event::ALT))
    , inportVolume_(Port::INPORT, "volumehandle.in", "Volume Input")
    , outportProxyGeometry_(Port::OUTPORT, "geometry.proxy", "Proxy Geometry")
    , outportRenderGeometry_(Port::OUTPORT, "geometry.render", "Render Geometry")
    , pickingBuffer_(Port::OUTPORT, "pickingBuffer")
{
    addPort(inportVolume_);
    addPort(outportProxyGeometry_);
    addPort(outportRenderGeometry_);
    addPrivateRenderPort(pickingBuffer_);

    canvasSize_.onChange(MemberFunctionCallback<ExplosionProxyGeometry>(this, &ExplosionProxyGeometry::canvasSizeChanged));
    canvasSize_.setReadOnlyFlag(true);
    canvasSize_.setVisibleFlag(false);
    addProperty(canvasSize_);

    brickingMode_.addOption("regularGrid", "Regular Grid");
    brickingMode_.addOption("manualGrid", "Manual Grid");
    brickingMode_.addOption("custom", "Custom");
    addProperty(brickingMode_);
    brickingMode_.onChange(MemberFunctionCallback<ExplosionProxyGeometry>(this, &ExplosionProxyGeometry::onBrickingModeChange));
    showBoundingBoxes_.addOption("all", "All");
    showBoundingBoxes_.addOption("none", "None");
    showBoundingBoxes_.addOption("selected", "Selected");
    addProperty(showBoundingBoxes_);
    addProperty(brickColor_);
    addProperty(selectedBrickColor_);
    addProperty(hideNonSelectedBricks_);
    hideNonSelectedBricks_.onChange(MemberFunctionOneParameterCallback<ExplosionProxyGeometry, bool>
                                    (this, &ExplosionProxyGeometry::onHideBricksChange, false));
    addProperty(hideSelectedBricks_);
    hideSelectedBricks_.onChange(MemberFunctionOneParameterCallback<ExplosionProxyGeometry, bool>
                                 (this, &ExplosionProxyGeometry::onHideBricksChange, true));
    addProperty(showAllBricks_);
    showAllBricks_.onChange(MemberFunctionCallback<ExplosionProxyGeometry>
                            (this, &ExplosionProxyGeometry::onShowAllBricksChange));


    //Properties only for regularGrid bricking mode
    addProperty(numBricksXp_);
    numBricksXp_.onChange(MemberFunctionCallback<ExplosionProxyGeometry>(this, &ExplosionProxyGeometry::doBricking));
    addProperty(numBricksYp_);
    numBricksYp_.onChange(MemberFunctionCallback<ExplosionProxyGeometry>(this, &ExplosionProxyGeometry::doBricking));
    addProperty(numBricksZp_);
    numBricksZp_.onChange(MemberFunctionCallback<ExplosionProxyGeometry>(this, &ExplosionProxyGeometry::doBricking));
    addProperty(regularExplosionGap_);
    regularExplosionGap_.onChange(MemberFunctionCallback<ExplosionProxyGeometry>(this, &ExplosionProxyGeometry::onExplosionGapChange));

    //Properties for manual and custom bricking mode
    numBricksXm_ = 1;
    numBricksYm_ = 1;
    numBricksZm_ = 1;
    addProperty(XbrickingPlane_);
    addProperty(addXbrickingPlane_);
    addXbrickingPlane_.onChange(MemberFunctionOneParameterCallback<ExplosionProxyGeometry, char>
                                (this, &ExplosionProxyGeometry::onAddBrickingPlaneChange, 'x'));
    addProperty(deleteLastXbrickingPlane_);
    deleteLastXbrickingPlane_.onChange(MemberFunctionOneParameterCallback<ExplosionProxyGeometry, char>
                                       (this, &ExplosionProxyGeometry::onDeleteLastBrickingPlaneChange, 'x'));
    addProperty(clearXbrickingList_);
    clearXbrickingList_.onChange(MemberFunctionOneParameterCallback<ExplosionProxyGeometry, char>
                                 (this, &ExplosionProxyGeometry::onClearBrickingListChange, 'x'));
    addProperty(YbrickingPlane_);
    addProperty(addYbrickingPlane_);
    addYbrickingPlane_.onChange(MemberFunctionOneParameterCallback<ExplosionProxyGeometry, char>
                                (this, &ExplosionProxyGeometry::onAddBrickingPlaneChange, 'y'));
    addProperty(deleteLastYbrickingPlane_);
    deleteLastYbrickingPlane_.onChange(MemberFunctionOneParameterCallback<ExplosionProxyGeometry, char>
                                       (this, &ExplosionProxyGeometry::onDeleteLastBrickingPlaneChange, 'y'));
    addProperty(clearYbrickingList_);
    clearYbrickingList_.onChange(MemberFunctionOneParameterCallback<ExplosionProxyGeometry, char>
                                 (this, &ExplosionProxyGeometry::onClearBrickingListChange, 'y'));
    addProperty(ZbrickingPlane_);
    addProperty(addZbrickingPlane_);
    addZbrickingPlane_.onChange(MemberFunctionOneParameterCallback<ExplosionProxyGeometry, char>
                                (this, &ExplosionProxyGeometry::onAddBrickingPlaneChange, 'z'));
    addProperty(deleteLastZbrickingPlane_);
    deleteLastZbrickingPlane_.onChange(MemberFunctionOneParameterCallback<ExplosionProxyGeometry, char>
                                       (this, &ExplosionProxyGeometry::onDeleteLastBrickingPlaneChange, 'z'));
    addProperty(clearZbrickingList_);
    clearZbrickingList_.onChange(MemberFunctionOneParameterCallback<ExplosionProxyGeometry, char>
                                 (this, &ExplosionProxyGeometry::onClearBrickingListChange, 'z'));
    addProperty(manualExplosionGap_);
    manualExplosionGap_.onChange(MemberFunctionCallback<ExplosionProxyGeometry>(this, &ExplosionProxyGeometry::onExplosionGapChange));
    addProperty(translation_);
    translation_.onChange(MemberFunctionCallback<ExplosionProxyGeometry>(this, &ExplosionProxyGeometry::onTranslationChange));


    addProperty(takeForCustom_);
    takeForCustom_.onChange(MemberFunctionCallback<ExplosionProxyGeometry>(this, &ExplosionProxyGeometry::takeForCustom));

    addProperty(undo_);
    undo_.onChange(MemberFunctionCallback<ExplosionProxyGeometry>(this, &ExplosionProxyGeometry::undo));

    addProperty(resetMode_);
    resetMode_.onChange(MemberFunctionCallback<ExplosionProxyGeometry>(this, &ExplosionProxyGeometry::resetMode));
    addProperty(camera_);

    addEventProperty(selectingOneBrickEvent_);
    addEventProperty(selectingSeveralBricksEvent_);
    addEventProperty(translateSelectedBricksEvent_);
    isClicked_ = false;

    onBrickingModeChange();

}

ExplosionProxyGeometry::~ExplosionProxyGeometry() {
}

Processor* ExplosionProxyGeometry::create() const {
    return new ExplosionProxyGeometry();
}

void ExplosionProxyGeometry::initialize() {
    RenderProcessor::initialize();

    idManager_.setRenderTarget(pickingBuffer_.getRenderTarget());
    idManager_.initializeTarget();

    onBrickingModeChange();
}

void ExplosionProxyGeometry::deinitialize() {
    outportProxyGeometry_.setData(0);
    outportRenderGeometry_.setData(0);

    RenderProcessor::deinitialize();
}

bool ExplosionProxyGeometry::isReady() const {
    return (inportVolume_.isReady() && (outportProxyGeometry_.isConnected() || outportRenderGeometry_.isConnected()));
}

void ExplosionProxyGeometry::beforeProcess() {
    RenderProcessor::beforeProcess();

    if (inportVolume_.hasChanged()) {
        if (brickingMode_.isSelected("regularGrid")) {
            regularBricklist_.clear();
            updatePlaneWidgets();
        }
        else {
            resetMode();
        }
    }
}

void ExplosionProxyGeometry::process() {
    tgtAssert(inportVolume_.getData(), "no volume");

    /////////////////////////////Bricking process beginning/////////////////////////////////////
    GeometrySequence* proxyGeometry = new GeometrySequence();
    MeshListGeometry* renderGeometry = new MeshListGeometry();

    tgt::mat4 tToW = inportVolume_.getData()->getTextureToWorldMatrix();

    //regularGrid bricking
    if (brickingMode_.get() == "regularGrid") {
        if (regularBricklist_.empty())
            createRegularOrManualBricklist();

        for(size_t i=0; i<regularBricklist_.size(); i++) {
            TriangleMeshGeometryColorNormal* tmg = TriangleMeshGeometryColorNormal::createCube(regularBricklist_.at(i).coordLlf, regularBricklist_.at(i).coordUrb, regularBricklist_.at(i).texLlf, regularBricklist_.at(i).texUrb, 1.0f);
            tmg->setTransformationMatrix(tToW * tgt::mat4::createTranslation(regularBricklist_.at(i).translation));
            proxyGeometry->addGeometry(tmg);
        }

        if (showBoundingBoxes_.get() == "all" || showBoundingBoxes_.get() == "selected") {
            for(size_t i=0; i<regularBricklist_.size(); i++) {
                MeshGeometry mg = MeshGeometry::createCube(
                    regularBricklist_.at(i).coordLlf,
                    regularBricklist_.at(i).coordUrb,
                    regularBricklist_.at(i).texLlf,
                    regularBricklist_.at(i).texUrb,
                    brickColor_.get().xyz(),
                    brickColor_.get().xyz());

                mg.setTransformationMatrix(tToW * tgt::mat4::createTranslation(regularBricklist_.at(i).translation));
                renderGeometry->addMesh(mg);
            }
        }
    }
    //manualGrid bricking
    else if (brickingMode_.get() == "manualGrid") {
        if (manualBricklist_.empty())
            doBricking();

        for(size_t i=0; i<manualBricklist_.size(); i++) {
            TriangleMeshGeometryColorNormal* tmg = TriangleMeshGeometryColorNormal::createCube(manualBricklist_.at(i).coordLlf, manualBricklist_.at(i).coordUrb, manualBricklist_.at(i).texLlf, manualBricklist_.at(i).texUrb, 1.0f);
            tmg->transform(tToW);
            tmg->transform(tgt::mat4::createTranslation(manualBricklist_.at(i).translation));
            proxyGeometry->addGeometry(tmg);
        }

        if (showBoundingBoxes_.get() == "all" || showBoundingBoxes_.get() == "selected") {
            for(size_t i=0; i<manualBricklist_.size(); i++) {
                MeshGeometry mg = MeshGeometry::createCube(
                            manualBricklist_.at(i).coordLlf,
                            manualBricklist_.at(i).coordUrb,
                            manualBricklist_.at(i).texLlf,
                            manualBricklist_.at(i).texUrb,
                            brickColor_.get().xyz(),
                            brickColor_.get().xyz());
                mg.setTransformationMatrix(tToW * tgt::mat4::createTranslation(regularBricklist_.at(i).translation));
                renderGeometry->addMesh(mg);
            }
        }
    }
    //custom bricking
    else{
        if(customBricklist_.empty())
            doBricking();

        registerForSelecting(customBricklist_);

        for(size_t i=0; i<customBricklist_.size(); i++) {
            if (customBricklist_.at(i).selected) {
                if (showBoundingBoxes_.get() == "all" || showBoundingBoxes_.get() == "selected") {
                    MeshGeometry mg = MeshGeometry::createCube(
                                customBricklist_.at(i).coordLlf,
                                customBricklist_.at(i).coordUrb,
                                customBricklist_.at(i).texLlf,
                                customBricklist_.at(i).texUrb,
                                selectedBrickColor_.get().xyz(),
                                selectedBrickColor_.get().xyz());

                    mg.setTransformationMatrix(tToW * tgt::mat4::createTranslation(regularBricklist_.at(i).translation));
                    renderGeometry->addMesh(mg);
                }

                TriangleMeshGeometryColorNormal* tmg = TriangleMeshGeometryColorNormal::createCube(customBricklist_.at(i).coordLlf, customBricklist_.at(i).coordUrb, customBricklist_.at(i).texLlf, customBricklist_.at(i).texUrb, 1.0f);
                tmg->transform(tToW);
                tmg->transform(tgt::mat4::createTranslation(customBricklist_.at(i).translation));
                proxyGeometry->addGeometry(tmg);
            }
            else {
                if (!customBricklist_.at(i).hidden) {
                    if (showBoundingBoxes_.get() == "all") {
                        MeshGeometry mg = MeshGeometry::createCube(
                                    customBricklist_.at(i).coordLlf,
                                    customBricklist_.at(i).coordUrb,
                                    customBricklist_.at(i).texLlf,
                                    customBricklist_.at(i).texUrb,
                                    brickColor_.get().xyz(),
                                    brickColor_.get().xyz());

                        mg.setTransformationMatrix(tToW * tgt::mat4::createTranslation(regularBricklist_.at(i).translation));
                        renderGeometry->addMesh(mg);
                    }
                    TriangleMeshGeometryColorNormal* tmg = TriangleMeshGeometryColorNormal::createCube(customBricklist_.at(i).coordLlf, customBricklist_.at(i).coordUrb, customBricklist_.at(i).texLlf, customBricklist_.at(i).texUrb, 1.0f);
                    tmg->transform(tToW);
                    tmg->transform(tgt::mat4::createTranslation(customBricklist_.at(i).translation));
                    proxyGeometry->addGeometry(tmg);
                }
            }
        }

    }

    outportProxyGeometry_.setData(proxyGeometry);
    outportRenderGeometry_.setData(renderGeometry);
    LGL_ERROR;
}

void ExplosionProxyGeometry::onBrickingModeChange() {

    if (brickingMode_.get() == "regularGrid") {
        //showBoundingBoxes_.set("all");
        selectedBrickColor_.setVisibleFlag(false);
        hideNonSelectedBricks_.setVisibleFlag(false);
        hideSelectedBricks_.setVisibleFlag(false);
        showAllBricks_.setVisibleFlag(false);
        numBricksXp_.setVisibleFlag(true);
        numBricksYp_.setVisibleFlag(true);
        numBricksZp_.setVisibleFlag(true);
        XbrickingPlane_.setVisibleFlag(false);
        addXbrickingPlane_.setVisibleFlag(false);
        deleteLastXbrickingPlane_.setVisibleFlag(false);
        clearXbrickingList_.setVisibleFlag(false);
        YbrickingPlane_.setVisibleFlag(false);
        addYbrickingPlane_.setVisibleFlag(false);
        deleteLastYbrickingPlane_.setVisibleFlag(false);
        clearYbrickingList_.setVisibleFlag(false);
        ZbrickingPlane_.setVisibleFlag(false);
        addZbrickingPlane_.setVisibleFlag(false);
        deleteLastZbrickingPlane_.setVisibleFlag(false);
        clearZbrickingList_.setVisibleFlag(false);
        regularExplosionGap_.setVisibleFlag(true);
        manualExplosionGap_.setVisibleFlag(false);
        translation_.setVisibleFlag(false);
        takeForCustom_.setVisibleFlag(true);
        undo_.setVisibleFlag(false);
        resetMode_.setVisibleFlag(true);
        if (inportVolume_.hasData())
            updatePlaneWidgets();
    }
    else if (brickingMode_.get() == "manualGrid") {
        doBricking();
        showBoundingBoxes_.set("all");
        selectedBrickColor_.setVisibleFlag(false);
        hideNonSelectedBricks_.setVisibleFlag(false);
        hideSelectedBricks_.setVisibleFlag(false);
        showAllBricks_.setVisibleFlag(false);
        numBricksXp_.setVisibleFlag(false);
        numBricksYp_.setVisibleFlag(false);
        numBricksZp_.setVisibleFlag(false);
        XbrickingPlane_.setVisibleFlag(true);
        addXbrickingPlane_.setVisibleFlag(true);
        deleteLastXbrickingPlane_.setVisibleFlag(true);
        clearXbrickingList_.setVisibleFlag(true);
        YbrickingPlane_.setVisibleFlag(true);
        addYbrickingPlane_.setVisibleFlag(true);
        deleteLastYbrickingPlane_.setVisibleFlag(true);
        clearYbrickingList_.setVisibleFlag(true);
        ZbrickingPlane_.setVisibleFlag(true);
        addZbrickingPlane_.setVisibleFlag(true);
        deleteLastZbrickingPlane_.setVisibleFlag(true);
        clearZbrickingList_.setVisibleFlag(true);
        regularExplosionGap_.setVisibleFlag(false);
        manualExplosionGap_.setVisibleFlag(true);
        translation_.setVisibleFlag(false);
        if (selectedBricks_.size() == 1)
            translation_.set(selectedBricks_.at(0).translation);
        else if (selectedBricks_.size() > 1)
            translation_.set(findMiddle(selectedBricks_));
        else
            translation_.set(tgt::vec3(0.f, 0.f, 0.f));

        takeForCustom_.setVisibleFlag(true);
        undo_.setVisibleFlag(false);
        resetMode_.setVisibleFlag(true);
        if (inportVolume_.hasData())
            updatePlaneWidgets();
    }
    else{
        doBricking();
        showBoundingBoxes_.set("all");
        selectedBrickColor_.setVisibleFlag(true);
        hideNonSelectedBricks_.setVisibleFlag(true);
        hideSelectedBricks_.setVisibleFlag(true);
        showAllBricks_.setVisibleFlag(true);
        numBricksXp_.setVisibleFlag(false);
        numBricksYp_.setVisibleFlag(false);
        numBricksZp_.setVisibleFlag(false);
        XbrickingPlane_.setVisibleFlag(true);
        addXbrickingPlane_.setVisibleFlag(true);
        deleteLastXbrickingPlane_.setVisibleFlag(false);
        clearXbrickingList_.setVisibleFlag(false);
        YbrickingPlane_.setVisibleFlag(true);
        addYbrickingPlane_.setVisibleFlag(true);
        deleteLastYbrickingPlane_.setVisibleFlag(false);
        clearYbrickingList_.setVisibleFlag(false);
        ZbrickingPlane_.setVisibleFlag(true);
        addZbrickingPlane_.setVisibleFlag(true);
        deleteLastZbrickingPlane_.setVisibleFlag(false);
        clearZbrickingList_.setVisibleFlag(false);
        regularExplosionGap_.setVisibleFlag(false);
        manualExplosionGap_.setVisibleFlag(false);
        translation_.setVisibleFlag(true);
        takeForCustom_.setVisibleFlag(false);
        undo_.setVisibleFlag(true);
        resetMode_.setVisibleFlag(true);
        if (inportVolume_.hasData())
            updatePlaneWidgets();
    }
}

void ExplosionProxyGeometry::onHideBricksChange(bool selected) {
    size_t i;

    if (selected)
        selectedBricks_.clear();

    for (i = 0; i < customBricklist_.size(); ++i) {
        if (customBricklist_.at(i).selected) {
            if (selected) {
                customBricklist_.at(i).selected = false;
                customBricklist_.at(i).hidden = true;
            }
        }
        else if (!customBricklist_.at(i).hidden) {
            if (selected) {
                customBricklist_.at(i).selected = true;
                selectedBricks_.insert(selectedBricks_.begin(), customBricklist_.at(i));
            }
            else{
                customBricklist_.at(i).hidden = true;
            }
        }
    }

    updatePlaneWidgets();
    invalidate();
}

void ExplosionProxyGeometry::onShowAllBricksChange() {
    size_t i;

    for (i = 0; i < customBricklist_.size(); ++i) {
        customBricklist_.at(i).hidden = false;
        if (selectedBricks_.empty())
            customBricklist_.at(i).selected = true;
    }

    if (selectedBricks_.empty())
        selectedBricks_ = customBricklist_;

    updatePlaneWidgets();
    invalidate();
}

void ExplosionProxyGeometry::doBricking() {
    if (!inportVolume_.getData())
        return;

    const VolumeBase* volume = inportVolume_.getData();

    //tgt::vec3 volumeSize = volume->getCubeSize();;
    tgt::vec3 volumeSize(1.0f);
    tgt::vec3 offsetLlf;
    tgt::vec3 offsetUrb;
    tgt::vec3 texLlf;
    tgt::vec3 texUrb;

    if (brickingMode_.isSelected("regularGrid")) {
        createRegularOrManualBricklist();
    }
    else if (brickingMode_.isSelected("manualGrid")) {
        offsetLlf = tgt::vec3(0.f, 0.f, 0.f);
        offsetUrb = tgt::vec3(volumeSize.x, volumeSize.y, -volumeSize.z);
        texLlf = tgt::vec3(0.f, 0.f, 1.f);
        texUrb = tgt::vec3(1.f, 1.f, 0.f);

        if (xPlanescoord_.empty()) {
            xPlanescoord_.insert(xPlanescoord_.begin(), offsetLlf.x);
            xPlanescoord_.insert(xPlanescoord_.end(), offsetUrb.x);
            xPlanestex_.insert(xPlanestex_.begin(), texLlf.x);
            xPlanestex_.insert(xPlanestex_.end(), texUrb.x);
        }

        if (yPlanescoord_.empty()) {
            yPlanescoord_.insert(yPlanescoord_.begin(), offsetLlf.y);
            yPlanescoord_.insert(yPlanescoord_.end(), offsetUrb.y);
            yPlanestex_.insert(yPlanestex_.begin(), texLlf.y);
            yPlanestex_.insert(yPlanestex_.end(), texUrb.y);
        }

        if (zPlanescoord_.empty()) {
            zPlanescoord_.insert(zPlanescoord_.begin(), offsetLlf.z);
            zPlanescoord_.insert(zPlanescoord_.end(), offsetUrb.z);
            zPlanestex_.insert(zPlanestex_.begin(), texLlf.z);
            zPlanestex_.insert(zPlanestex_.end(), texUrb.z);
        }

        createRegularOrManualBricklist();
    }
    else{
        if (customBricklist_.empty())
            createCustomBricklist(' ', 0, 0);
    }
}

void ExplosionProxyGeometry::onAddBrickingPlaneChange(char axis) {
    if (axis == 'x')
        addBrickingPlane('x', numBricksXm_, xPlanescoord_, xPlanestex_, xCoords_);

    else if (axis == 'y')
        addBrickingPlane('y', numBricksYm_, yPlanescoord_, yPlanestex_, yCoords_);

    else
        addBrickingPlane('z', numBricksZm_, zPlanescoord_, zPlanestex_, zCoords_);
}

void ExplosionProxyGeometry::addBrickingPlane(char axis, size_t& numBricks, std::vector<float>& planescoord,
                                             std::vector<float>& planestex, std::stack<float>& coords) {

    const VolumeBase* volume = inportVolume_.getData();
    //tgt::vec3 volumeSize = volume->getCubeSize();
    tgt::vec3 volumeSize(1.0f);
    tgt::ivec3 numSlices = volume->getDimensions();
    tgt::vec3 texLlf(0, 0, 1);
    tgt::vec3 texUrb(1, 1, 0);
    size_t i;

    float tex;
    float offset;
    float offsetLlf;
    float offsetUrb;

    if (axis == 'x') {
        tex = XbrickingPlane_.get() / static_cast<tgt::vec3::ElemType>(numSlices.x);
        offset = -volumeSize.x * (texLlf.x - tex);
        if (brickingMode_.get() == "manualGrid") {
            offsetLlf = 0.f;
            offsetUrb = volumeSize.x;
        }
        else if (!selectedBricks_.empty()) {
            createCustomBricklist('x', offset, tex);
        }
    }
    else if (axis == 'y') {
        tex = YbrickingPlane_.get() / static_cast<tgt::vec3::ElemType>(numSlices.y);
        offset = -volumeSize.y * (texLlf.y - tex);
        if (brickingMode_.get() == "manualGrid") {
            offsetLlf = 0.f;
            offsetUrb = volumeSize.y;
        }
        else if (!selectedBricks_.empty()) {
            createCustomBricklist('y', offset, tex);
        }
    }
    else if (axis == 'z') {
        tex = ZbrickingPlane_.get() / static_cast<tgt::vec3::ElemType>(numSlices.z);
        offset = -volumeSize.z * (texLlf.z - tex);
        if (brickingMode_.get() == "manualGrid") {
            offsetLlf = 0.f;
            offsetUrb = -volumeSize.z;
        }
        else if (!selectedBricks_.empty()) {
            createCustomBricklist('z', offset, tex);
        }
    }

    if (brickingMode_.get() == "manualGrid") {
        if (planescoord.empty()) {
            planescoord.insert(planescoord.begin(), offset);
            planescoord.insert(planescoord.begin(), offsetLlf);
            planescoord.insert(planescoord.end(), offsetUrb);
            planestex.insert(planestex.begin(), tex);
            planestex.insert(planestex.begin(), texLlf.x);
            planestex.insert(planestex.end(), texUrb.x);
        }
        else{
            i=1;
            while (numBricks > (planescoord.size() - 2) && i < numBricks) {
                if (axis == 'x' || axis == 'y') {
                    if (planescoord.at(i) > offset) {
                        planescoord.insert(planescoord.begin() + i, offset);
                        planestex.insert(planestex.begin() + i, tex);
                    }
                    else{
                        ++i;
                    }
                }
                else if (axis == 'z') {
                    if (planescoord.at(i) < offset) {
                        planescoord.insert(planescoord.begin() + i, offset);
                        planestex.insert(planestex.begin() + i, tex);
                    }
                    else{
                        ++i;
                    }
                }
            }
            if (numBricks > (planescoord.size() - 2)) {
                planescoord.insert(planescoord.end() - 1, offset);
                planestex.insert(planestex.end() - 1, tex);
            }
        }
        coords.push(offset);
        numBricks = numBricks + 1;
        doBricking();
    }
    invalidate();
}

void ExplosionProxyGeometry::onDeleteLastBrickingPlaneChange(char axis) {
    if (axis == 'x')
        deleteLastBrickingPlane(numBricksXm_, xPlanescoord_, xPlanestex_, xCoords_);

    else if (axis == 'y')
        deleteLastBrickingPlane(numBricksYm_, yPlanescoord_, yPlanestex_, yCoords_);

    else if (axis == 'z')
        deleteLastBrickingPlane(numBricksZm_, zPlanescoord_, zPlanestex_, zCoords_);

    invalidate();
}

void ExplosionProxyGeometry::deleteLastBrickingPlane(size_t& numBricks, std::vector<float>& planescoord,
                                                    std::vector<float>& planestex, std::stack<float>& coords) {

    if (!coords.empty()) {
        float del = coords.top();
        coords.pop();
        bool deleted = false;
        int i = 0;

        while (!deleted) {
            if (planescoord.at(i) == del) {
                planescoord.erase(planescoord.begin() + i);
                planestex.erase(planestex.begin() + i);
                deleted = true;
            }
            else{
                ++i;
            }
        }

        numBricks = numBricks - 1;
        doBricking();
    }

    invalidate();
}

void ExplosionProxyGeometry::onClearBrickingListChange(char axis) {
    clearBrickingList(axis);
    doBricking();
}

void ExplosionProxyGeometry::clearBrickingList(char axis) {

    if (axis == 'x') {
        numBricksXm_ = 1;
        xPlanescoord_.clear();
        xPlanestex_.clear();
        while (!xCoords_.empty())
            xCoords_.pop();
    }
    else if (axis == 'y') {
        numBricksYm_ = 1;
        yPlanescoord_.clear();
        yPlanestex_.clear();
        while (!yCoords_.empty())
            yCoords_.pop();
    }
    else{
        numBricksZm_ = 1;
        zPlanescoord_.clear();
        zPlanestex_.clear();
        while (!zCoords_.empty())
            zCoords_.pop();
    }

    invalidate();
}

void ExplosionProxyGeometry::onExplosionGapChange() {
    if (inportVolume_.hasData()) {
        createRegularOrManualBricklist();
    }
}

void ExplosionProxyGeometry::onTranslationChange() {
    size_t i;

    if (selectedBricks_.size() == 1) {
        customBricklist_.at(selectedBricks_.at(0).brickId).translation = translation_.get();
        selectedBricks_.at(0).translation = translation_.get();
    }
    else if (selectedBricks_.size() > 1) {
        tgt::vec3 middle = findMiddle(selectedBricks_);
        tgt::vec3 distance = translation_.get() - middle;
        for (i = 0; i < selectedBricks_.size(); ++i) {
            customBricklist_.at(selectedBricks_.at(i).brickId).translation += distance;
            selectedBricks_.at(i).translation += distance;
        }
    }
    invalidate();
}

void ExplosionProxyGeometry::takeForCustom() {
    size_t i;

    if (brickingMode_.get() == "regularGrid") {
        brickingMode_.set("custom");
        resetMode();
        customBricklist_ = regularBricklist_;
    }
    else if (brickingMode_.get() == "manualGrid") {
        brickingMode_.set("custom");
        resetMode();
        customBricklist_ = manualBricklist_;
    }

    for (i = 0; i < customBricklist_.size(); ++i) {
        customBricklist_.at(i).brickId = static_cast<int>(i);

        if (customBricklist_.at(i).selected)
            selectedBricks_.insert(selectedBricks_.begin(), customBricklist_.at(i));
    }

    registerForSelecting(customBricklist_);
    onBrickingModeChange();

}

void ExplosionProxyGeometry::undo() {
    if (!customBricklistStack_.empty()) {
        customBricklist_.clear();
        customBricklist_ = customBricklistStack_.top();
        customBricklistStack_.pop();
        selectedBricks_.clear();

        size_t i;
        for (i = 0; i < customBricklist_.size(); ++i) {
            if (customBricklist_.at(i).selected)
                selectedBricks_.insert(selectedBricks_.begin(), customBricklist_.at(i));
        }

        if (selectedBricks_.empty())
            translation_.set(tgt::vec3(0.f, 0.f, 0.f));

        else if (selectedBricks_.size() == 1)
            translation_.set(selectedBricks_.at(0).translation);

        else
            translation_.set(findMiddle(selectedBricks_));

        registerForSelecting(customBricklist_);
        updatePlaneWidgets();
        invalidate();
    }
}

void ExplosionProxyGeometry::resetMode() {

    if (brickingMode_.isSelected("regularGrid")) {
        numBricksXp_.set(1);
        numBricksYp_.set(1);
        numBricksZp_.set(1);
        regularExplosionGap_.set(0);
        regularBricklist_.clear();
        updatePlaneWidgets();
    }
    else if (brickingMode_.isSelected("manualGrid")) {
        clearBrickingList('x');
        clearBrickingList('y');
        clearBrickingList('z');
        doBricking();
        updatePlaneWidgets();
        manualExplosionGap_.set(0);
    }
    else if (brickingMode_.isSelected("custom")) {
        while (!customBricklistStack_.empty())
            customBricklistStack_.pop();

        customBricklist_.clear();
        selectedBricks_.clear();
        translation_.set(tgt::vec3(0.f, 0.f, 0.f));
        updatePlaneWidgets();
        doBricking();
    }
    else {
        LERROR("Unknown bricking mode: " << brickingMode_.get());
    }

}

void ExplosionProxyGeometry::createRegularOrManualBricklist() {
    if (!inportVolume_.hasData()) {
        LWARNING("No volume!");
        return;
    }

    const VolumeBase* volume = inportVolume_.getData();
    //tgt::vec3 volumeSize = volume->getCubeSize();
    tgt::vec3 volumeSize(1.0f, 1.0f, 1.0f);
    tgt::ivec3 numSlices = volume->getDimensions();
    tgt::vec3 sliceWidth = getSliceWidth(volume);
    //tgt::vec3 coordLlf = volume->getLLF();
    //tgt::vec3 coordUrb = volume->getURB();;
    tgt::vec3 coordLlf(0, 0, 0);
    tgt::vec3 coordUrb(1, 1, 1);
    //std::swap(coordLlf.z, coordUrb.z);

    Brick current;
    int xOffset;
    int yOffset;
    int zOffset;

    if (brickingMode_.get() == "regularGrid") {
        regularBricklist_.clear();

        vec3 bricks = vec3(static_cast<float>(numBricksXp_.get()), static_cast<float>(numBricksYp_.get()), static_cast<float>(numBricksZp_.get()));
        vec3 g = vec3(-1.0f) * (bricks - vec3(1.0f));
        vec3 h = volumeSize / bricks;

        for (int x = 0; x < numBricksXp_.get(); ++x) {
            for (int y = 0; y < numBricksYp_.get(); ++y) {
                for (int z = 0; z < numBricksZp_.get(); ++z) {
                    vec3 coord(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z));

                    current.coordLlf = coordLlf + (coord * h);
                    current.coordUrb = coordLlf + ((coord + vec3(1.0f)) * h);

                    current.texLlf = coord / bricks;
                    current.texUrb = (coord + vec3(1.0f)) / bricks;

                    current.sliceIndexLlf = tgt::ivec3((coord * h) / sliceWidth);
                    current.sliceIndexUrb = tgt::ivec3(static_cast<int>((current.coordUrb.x - coordLlf.x) / sliceWidth.x), static_cast<int>((current.coordUrb.y - coordLlf.y) / sliceWidth.y), numSlices.z + static_cast<int>((current.coordUrb.z - coordLlf.z) / sliceWidth.z));

                    //current.translation = tgt::vec3(0.0f);
                    current.translation = (g * regularExplosionGap_.get() * 0.5f) + coord * regularExplosionGap_.get();
                    current.brickId = -1;
                    current.selected = true;
                    current.hidden = false;
                    regularBricklist_.push_back(current);
                }
            }
        }
    }
    else if (brickingMode_.get() == "manualGrid") {
        manualBricklist_.clear();
        tgt::vec3 startExplosionGap = -1.f * tgt::vec3((numBricksXm_ - 1) * manualExplosionGap_.get() / 2,
                                              (numBricksYm_ - 1) * manualExplosionGap_.get() / 2,
                                              -1.f * (numBricksZm_ - 1) * manualExplosionGap_.get() / 2);

        for (int it = 0; it < static_cast<int>(numBricksXm_ * numBricksYm_ * numBricksZm_); ++it) {
            xOffset = it / static_cast<int>(numBricksYm_ * numBricksZm_);
            yOffset = (it / static_cast<int>(numBricksZm_)) % numBricksYm_;
            zOffset = it % numBricksZm_;

            current.coordLlf = coordLlf + tgt::vec3(xPlanescoord_.at(xOffset),yPlanescoord_.at(yOffset),(zPlanescoord_.at(zOffset)));
            current.coordUrb = coordLlf + tgt::vec3(xPlanescoord_.at(xOffset + 1),yPlanescoord_.at(yOffset + 1), (zPlanescoord_.at(zOffset + 1)));
            current.texLlf = tgt::vec3(xPlanestex_.at(xOffset),yPlanestex_.at(yOffset),zPlanestex_.at(zOffset));
            current.texUrb = tgt::vec3(xPlanestex_.at(xOffset + 1),yPlanestex_.at(yOffset + 1),zPlanestex_.at(zOffset + 1));
            current.sliceIndexLlf = tgt::ivec3(static_cast<int>(xPlanescoord_.at(xOffset) / sliceWidth.x),
                                               static_cast<int>(yPlanescoord_.at(yOffset) / sliceWidth.y),
                                               numSlices.z + static_cast<int>(zPlanescoord_.at(zOffset) / sliceWidth.z));
            current.sliceIndexUrb = tgt::ivec3(static_cast<int>(xPlanescoord_.at(xOffset + 1) / sliceWidth.x),
                                               static_cast<int>(yPlanescoord_.at(yOffset + 1) / sliceWidth.y),
                                               numSlices.z+static_cast<int>(zPlanescoord_.at(zOffset + 1) / sliceWidth.z));
            current.translation = startExplosionGap +  tgt::vec3(static_cast<float>(xOffset) * manualExplosionGap_.get(),
                                                        static_cast<float>(yOffset) * manualExplosionGap_.get(),
                                                        -static_cast<float>(zOffset) * manualExplosionGap_.get());
            current.brickId = -1;
            current.selected = true;
            current.hidden = false;
            manualBricklist_.insert(manualBricklist_.end(), current);
        }
    }
}

void ExplosionProxyGeometry::createCustomBricklist(char addPlaneAxis, float addPlaneCoordOffset, float addPlaneTex) {
    if (!inportVolume_.hasData()) {
        LWARNING("No volume!");
        return;
    }

    const VolumeBase* volume = inportVolume_.getData();
    tgt::ivec3 numSlices = volume->getDimensions();
    tgt::vec3 sliceWidth = getSliceWidth(volume);
    //tgt::vec3 coordLlf = volume->getLLF();
    //tgt::vec3 coordUrb = volume->getURB();
    tgt::vec3 coordLlf(0, 0, 0);
    tgt::vec3 coordUrb(1, 1, 1);

    const tgt::vec3 texLlf(0, 0, 0);
    const tgt::vec3 texUrb(1, 1, 1);

    Brick current;
    std::vector<Brick> newBricks;
    size_t it;
    size_t i;

    if (customBricklist_.empty()) {
        current.coordLlf = coordLlf;
        current.coordUrb = coordUrb;
        current.texLlf = texLlf;
        current.texUrb = texUrb;
        current.sliceIndexLlf = tgt::ivec3(0, 0, numSlices.z);
        current.sliceIndexUrb = tgt::ivec3(numSlices.x, numSlices.y, 0);
        current.translation = tgt::vec3(0,0,0);
        current.brickId = 0;
        current.selected = true;
        current.hidden = false;
        customBricklist_.insert(customBricklist_.begin(), current);
    }
    else{
        for (i = 0; i < selectedBricks_.size(); ++i) {
            newBricks.clear();
            if (addPlaneAxis == 'x') {
                if (selectedBricks_.at(i).coordLlf.x < (coordLlf.x + addPlaneCoordOffset) &&
                    selectedBricks_.at(i).coordUrb.x > (coordLlf.x + addPlaneCoordOffset)) {

                    customBricklist_.erase(customBricklist_.begin() + selectedBricks_.at(i).brickId);

                    current.coordLlf = selectedBricks_.at(i).coordLlf;
                    current.coordUrb = tgt::vec3(coordLlf.x + addPlaneCoordOffset, selectedBricks_.at(i).coordUrb.y,
                                                 selectedBricks_.at(i).coordUrb.z);
                    current.texLlf = selectedBricks_.at(i).texLlf;
                    current.texUrb = tgt::vec3(addPlaneTex, selectedBricks_.at(i).texUrb.y, selectedBricks_.at(i).texUrb.z);
                    current.sliceIndexLlf = selectedBricks_.at(i).sliceIndexLlf;
                    current.sliceIndexUrb = tgt::ivec3(static_cast<int>(addPlaneCoordOffset / sliceWidth.x),
                                                       selectedBricks_.at(i).sliceIndexUrb.y,
                                                       selectedBricks_.at(i).sliceIndexUrb.z);
                    current.translation = selectedBricks_.at(i).translation;
                    current.brickId = -1;
                    current.selected = true;
                    current.hidden = false;
                    newBricks.insert(newBricks.end(), current);

                    current.coordLlf = tgt::vec3(coordLlf.x + addPlaneCoordOffset, selectedBricks_.at(i).coordLlf.y,
                                                 selectedBricks_.at(i).coordLlf.z);
                    current.coordUrb = selectedBricks_.at(i).coordUrb;
                    current.texLlf = tgt::vec3(addPlaneTex, selectedBricks_.at(i).texLlf.y, selectedBricks_.at(i).texLlf.z);
                    current.texUrb = selectedBricks_.at(i).texUrb;
                    current.sliceIndexLlf = tgt::ivec3(static_cast<int>(addPlaneCoordOffset / sliceWidth.x),
                                                       selectedBricks_.at(i).sliceIndexLlf.y,
                                                       selectedBricks_.at(i).sliceIndexLlf.z);
                    current.sliceIndexUrb = selectedBricks_.at(i).sliceIndexUrb;
                    current.translation = selectedBricks_.at(i).translation;
                    current.brickId = -1;
                    current.selected = true;
                    current.hidden = false;
                    newBricks.insert(newBricks.end(), current);

                    for (it = 0; it < newBricks.size(); ++it) {
                        customBricklist_.insert(customBricklist_.begin() + selectedBricks_.at(i).brickId, newBricks.at(it));
                    }
                }

            }
            else if (addPlaneAxis == 'y') {
                if (selectedBricks_.at(i).coordLlf.y < (coordLlf.y + addPlaneCoordOffset) &&
                    selectedBricks_.at(i).coordUrb.y > (coordLlf.y + addPlaneCoordOffset)) {

                    customBricklist_.erase(customBricklist_.begin() + selectedBricks_.at(i).brickId);

                    current.coordLlf = selectedBricks_.at(i).coordLlf;
                    current.coordUrb = tgt::vec3(selectedBricks_.at(i).coordUrb.x, coordLlf.y + addPlaneCoordOffset,
                                                 selectedBricks_.at(i).coordUrb.z);
                    current.texLlf = selectedBricks_.at(i).texLlf;
                    current.texUrb = tgt::vec3(selectedBricks_.at(i).texUrb.x, addPlaneTex, selectedBricks_.at(i).texUrb.z);
                    current.sliceIndexLlf = selectedBricks_.at(i).sliceIndexLlf;
                    current.sliceIndexUrb = tgt::ivec3(selectedBricks_.at(i).sliceIndexUrb.x,
                                                       static_cast<int>(addPlaneCoordOffset / sliceWidth.y),
                                                       selectedBricks_.at(i).sliceIndexUrb.z);
                    current.translation = selectedBricks_.at(i).translation;
                    current.brickId = -1;
                    current.selected = true;
                    current.hidden = false;
                    newBricks.insert(newBricks.end(), current);

                    current.coordLlf = tgt::vec3(selectedBricks_.at(i).coordLlf.x, coordLlf.y + addPlaneCoordOffset,
                                                 selectedBricks_.at(i).coordLlf.z);
                    current.coordUrb = selectedBricks_.at(i).coordUrb;
                    current.texLlf = tgt::vec3(selectedBricks_.at(i).texLlf.x, addPlaneTex, selectedBricks_.at(i).texLlf.z);
                    current.texUrb = selectedBricks_.at(i).texUrb;
                    current.sliceIndexLlf = tgt::ivec3(selectedBricks_.at(i).sliceIndexLlf.x,
                                                       static_cast<int>(addPlaneCoordOffset / sliceWidth.y),
                                                       selectedBricks_.at(i).sliceIndexLlf.z);
                    current.sliceIndexUrb = selectedBricks_.at(i).sliceIndexUrb;
                    current.translation = selectedBricks_.at(i).translation;
                    current.brickId = -1;
                    current.selected = true;
                    current.hidden = false;
                    newBricks.insert(newBricks.end(), current);

                    for (it = 0; it < newBricks.size(); ++it) {
                        customBricklist_.insert(customBricklist_.begin() + selectedBricks_.at(i).brickId, newBricks.at(it));
                    }
                }
            }
            else{
                if (selectedBricks_.at(i).coordLlf.z > (coordLlf.z + addPlaneCoordOffset) &&
                    selectedBricks_.at(i).coordUrb.z < (coordLlf.z + addPlaneCoordOffset)) {

                    customBricklist_.erase(customBricklist_.begin() + selectedBricks_.at(i).brickId);

                    current.coordLlf = selectedBricks_.at(i).coordLlf;
                    current.coordUrb = tgt::vec3(selectedBricks_.at(i).coordUrb.x, selectedBricks_.at(i).coordUrb.y,
                                                 coordLlf.z + addPlaneCoordOffset);
                    current.texLlf = selectedBricks_.at(i).texLlf;
                    current.texUrb = tgt::vec3(selectedBricks_.at(i).texUrb.x, selectedBricks_.at(i).texUrb.y,  addPlaneTex);
                    current.sliceIndexLlf = selectedBricks_.at(i).sliceIndexLlf;
                    current.sliceIndexUrb = tgt::ivec3(selectedBricks_.at(i).sliceIndexUrb.x, selectedBricks_.at(i).sliceIndexUrb.y,
                                                       numSlices.z + static_cast<int>(addPlaneCoordOffset / sliceWidth.z));
                    current.translation = selectedBricks_.at(i).translation;
                    current.brickId = -1;
                    current.selected = true;
                    current.hidden = false;
                    newBricks.insert(newBricks.end(), current);

                    current.coordLlf = tgt::vec3(selectedBricks_.at(i).coordLlf.x, selectedBricks_.at(i).coordLlf.y,
                                                 coordLlf.z + addPlaneCoordOffset);
                    current.coordUrb = selectedBricks_.at(i).coordUrb;
                    current.texLlf = tgt::vec3(selectedBricks_.at(i).texLlf.x, selectedBricks_.at(i).texLlf.y, addPlaneTex);
                    current.texUrb = selectedBricks_.at(i).texUrb;
                    current.sliceIndexLlf = tgt::ivec3(selectedBricks_.at(i).sliceIndexLlf.x, selectedBricks_.at(i).sliceIndexLlf.y,
                                                       numSlices.z + static_cast<int>(addPlaneCoordOffset / sliceWidth.y));
                    current.sliceIndexUrb = selectedBricks_.at(i).sliceIndexUrb;
                    current.translation = selectedBricks_.at(i).translation;
                    current.brickId = -1;
                    current.selected = true;
                    current.hidden = false;
                    newBricks.insert(newBricks.end(), current);

                    for (it = 0; it < newBricks.size(); ++it) {
                        customBricklist_.insert(customBricklist_.begin() + selectedBricks_.at(i).brickId, newBricks.at(it));
                    }
                }
            }
        }
    }

    registerForSelecting(customBricklist_);
    selectedBricks_.clear();

    for (i = 0; i < customBricklist_.size(); ++i) {
        customBricklist_.at(i).brickId=static_cast<int>(i);

        if (customBricklist_.at(i).selected)
            selectedBricks_.insert(selectedBricks_.begin(), customBricklist_.at(i));
    }

    customBricklistStack_.push(customBricklist_);
}

tgt::vec3 ExplosionProxyGeometry::findMiddle(std::vector<Brick> bricklist) {
    float xMiddle;
    float yMiddle;
    float zMiddle;
    float xMin = bricklist.at(0).coordLlf.x + bricklist.at(0).translation.x;
    float xMax = bricklist.at(0).coordUrb.x + bricklist.at(0).translation.x;
    float yMin = bricklist.at(0).coordLlf.y + bricklist.at(0).translation.y;
    float yMax = bricklist.at(0).coordUrb.y + bricklist.at(0).translation.y;
    float zMin = bricklist.at(0).coordLlf.z + bricklist.at(0).translation.z;
    float zMax = bricklist.at(0).coordUrb.z + bricklist.at(0).translation.z;
    size_t i;


    for (i = 1; i < bricklist.size(); ++i) {
        if (xMin > bricklist.at(i).coordLlf.x + bricklist.at(i).translation.x)
            xMin = bricklist.at(i).coordLlf.x + bricklist.at(i).translation.x;

        if (xMax < bricklist.at(i).coordUrb.x + bricklist.at(i).translation.x)
            xMax = bricklist.at(i).coordUrb.x + bricklist.at(i).translation.x;

        if (yMin > bricklist.at(i).coordLlf.y + bricklist.at(i).translation.y)
            yMin = bricklist.at(i).coordLlf.y + bricklist.at(i).translation.y;

        if (yMax < bricklist.at(i).coordUrb.y + bricklist.at(i).translation.y)
            yMax = bricklist.at(i).coordUrb.y + bricklist.at(i).translation.y;

        if (zMin < bricklist.at(i).coordLlf.z + bricklist.at(i).translation.z)
            zMin = bricklist.at(i).coordLlf.z + bricklist.at(i).translation.z;

        if (zMax > bricklist.at(i).coordUrb.z + bricklist.at(i).translation.z)
            zMax = bricklist.at(i).coordUrb.z + bricklist.at(i).translation.z;
    }

    xMiddle = xMax - ((xMax - xMin) / 2.f);
    yMiddle = yMax - ((yMax - yMin) / 2.f);
    zMiddle = zMin - ((zMin - zMax) / 2.f);

    return tgt::vec3(xMiddle, yMiddle, zMiddle);
}

float ExplosionProxyGeometry::getCameraDistance(tgt::vec3 cameraposition, tgt::vec3 brickmiddle) {
    float distance = sqrt((cameraposition.x - brickmiddle.x) * (cameraposition.x - brickmiddle.x)
                            + (cameraposition.y - brickmiddle.y) * (cameraposition.y - brickmiddle.y)
                            +(cameraposition.z - brickmiddle.z) * (cameraposition.z - brickmiddle.z));
    return distance;
}

tgt::vec3 ExplosionProxyGeometry::getSliceWidth(const VolumeBase* volume) {
    tgt::vec3 sliceWidth;
    //sliceWidth = volume->getCubeSize() / tgt::vec3(volume->getDimensions());
    sliceWidth = tgt::vec3(1.0f) / tgt::vec3(volume->getDimensions());
    return sliceWidth;
}

void ExplosionProxyGeometry::updatePlaneWidgets() {
    if (brickingMode_.get() == "regularGrid" || brickingMode_.get() == "manualGrid" || selectedBricks_.empty()) {
        tgt::ivec3 numSlices = inportVolume_.getData()->getRepresentation<VolumeRAM>()->getDimensions();

        XbrickingPlane_.setMinValue(0);
        XbrickingPlane_.set(0);
        XbrickingPlane_.setMaxValue(numSlices.x);
        XbrickingPlane_.updateWidgets();

        YbrickingPlane_.setMinValue(0);
        YbrickingPlane_.set(0);
        YbrickingPlane_.setMaxValue(numSlices.y);
        YbrickingPlane_.updateWidgets();

        ZbrickingPlane_.setMinValue(0);
        ZbrickingPlane_.set(0);
        ZbrickingPlane_.setMaxValue(numSlices.z);
        ZbrickingPlane_.updateWidgets();
    }
    else{
        XbrickingPlane_.setMinValue(findMinSlice(selectedBricks_, 'x'));
        XbrickingPlane_.set(XbrickingPlane_.getMinValue());
        XbrickingPlane_.setMaxValue(findMaxSlice(selectedBricks_, 'x'));
        XbrickingPlane_.updateWidgets();
        YbrickingPlane_.setMinValue(findMinSlice(selectedBricks_, 'y'));
        YbrickingPlane_.set(YbrickingPlane_.getMinValue());
        YbrickingPlane_.setMaxValue(findMaxSlice(selectedBricks_, 'y'));
        YbrickingPlane_.updateWidgets();
        ZbrickingPlane_.setMinValue(findMaxSlice(selectedBricks_, 'z'));
        ZbrickingPlane_.set(ZbrickingPlane_.getMinValue());
        ZbrickingPlane_.setMaxValue(findMinSlice(selectedBricks_, 'z'));
        ZbrickingPlane_.updateWidgets();
    }
}

int ExplosionProxyGeometry::findMinSlice(std::vector<Brick> bricklist, char axis) {
    tgt::ivec3 sliceIndex = bricklist.at(0).sliceIndexLlf;
    size_t i;

    if (axis == 'x') {
        for (i = 1; i < bricklist.size(); ++i) {
            if (sliceIndex.x > bricklist.at(i).sliceIndexLlf.x)
                    sliceIndex.x = bricklist.at(i).sliceIndexLlf.x;
        }
        return sliceIndex.x;
    }
    if (axis == 'y') {
        for (i = 1; i < bricklist.size(); ++i) {
            if (sliceIndex.y > bricklist.at(i).sliceIndexLlf.y)
                    sliceIndex.y = bricklist.at(i).sliceIndexLlf.y;
        }
        return sliceIndex.y;
    }
    else{
        for (i = 1; i < bricklist.size(); ++i) {
            if (sliceIndex.z < bricklist.at(i).sliceIndexLlf.z)
                    sliceIndex.z = bricklist.at(i).sliceIndexLlf.z;
        }
        return sliceIndex.z;
    }
}

int ExplosionProxyGeometry::findMaxSlice(std::vector<Brick> bricklist, char axis) {
    tgt::ivec3 sliceIndex = bricklist.at(0).sliceIndexUrb;
    size_t i;

    if (axis == 'x') {
        for (i = 1; i < bricklist.size(); ++i) {
            if (sliceIndex.x < bricklist.at(i).sliceIndexUrb.x)
                    sliceIndex.x = bricklist.at(i).sliceIndexUrb.x;
        }
        return sliceIndex.x;
    }
    if (axis == 'y') {
        for (i = 1; i < bricklist.size(); ++i) {
            if (sliceIndex.y < bricklist.at(i).sliceIndexUrb.y)
                    sliceIndex.y = bricklist.at(i).sliceIndexUrb.y;
        }
        return sliceIndex.y;
    }
    else{
        for (i = 1; i < bricklist.size(); ++i) {
            if (sliceIndex.z > bricklist.at(i).sliceIndexUrb.z)
                    sliceIndex.z = bricklist.at(i).sliceIndexUrb.z;
        }
        return sliceIndex.z;
    }
}

void ExplosionProxyGeometry::registerForSelecting(std::vector<Brick> bricklist) {
    size_t i;

    idManager_.clearRegisteredObjects();
    idManager_.activateTarget(getID());
    idManager_.clearTarget();

    for (i = 0; i < bricklist.size(); ++i) {
        idManager_.registerObject((void *)(i+1));
    }

    // set modelview and projection matrices
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadMatrix(camera_.get().getProjectionMatrix(pickingBuffer_.getSize()));
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadMatrix(camera_.get().getViewMatrix());
    LGL_ERROR;

    for (i = 0; i < bricklist.size(); ++i) {
        tgt::col4 color = idManager_.getColorFromObject((void *)(i + 1));
        MeshGeometry selectingBrick = MeshGeometry::createCube(
                                                        bricklist.at(i).coordLlf + bricklist.at(i).translation,
                                                        bricklist.at(i).coordUrb + bricklist.at(i).translation,
                                                        bricklist.at(i).texLlf,
                                                        bricklist.at(i).texUrb,
                                                        tgt::vec3(color.xyz()) / 255.f,
                                                        tgt::vec3(color.xyz()) / 255.f);
        selectingBrick.render();
    }

    // restore matrices
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();
    LGL_ERROR;
}

void ExplosionProxyGeometry::onSelectingOneBrickEvent(tgt::MouseEvent * e) {
    if (!pickingBuffer_.hasRenderTarget()) {
        return;
    }
    registerForSelecting(customBricklist_);

    e->accept();
    size_t i;
    tgt::ivec2 pos = tgt::ivec2(e->coord().x, canvasSize_.get().y - e->coord().y);
    const void* brick = idManager_.getObjectAtPos(pos);
    if (brick != 0) {
        int brickIndex = static_cast<int>(reinterpret_cast<size_t>(brick)) - 1;
        if  (!customBricklist_.at(brickIndex).hidden) {
            customBricklistStack_.push(customBricklist_);

            for (i = 0; i < selectedBricks_.size(); ++i) {
                customBricklist_.at(selectedBricks_.at(i).brickId).selected = false;
            }

            customBricklist_.at(brickIndex).selected = true;
            selectedBricks_.clear();
            selectedBricks_.insert(selectedBricks_.begin(), customBricklist_.at(brickIndex));
            translation_.set(customBricklist_.at(brickIndex).translation);
            updatePlaneWidgets();
            invalidate();
        }
    }
    else{
        customBricklistStack_.push(customBricklist_);
        for (i = 0; i < selectedBricks_.size(); ++i) {
            customBricklist_.at(selectedBricks_.at(i).brickId).selected = false;
        }
        selectedBricks_.clear();
        translation_.set(tgt::vec3(0.f, 0.f, 0.f));
        updatePlaneWidgets();
        invalidate();
    }
}

void ExplosionProxyGeometry::onSelectingSeveralBricksEvent(tgt::MouseEvent * e) {
    if (!pickingBuffer_.hasRenderTarget()) {
        return;
    }
    registerForSelecting(customBricklist_);

    size_t i;
    e->accept();
    tgt::ivec2 pos = tgt::ivec2(e->coord().x, canvasSize_.get().y - e->coord().y);
    const void* brick = idManager_.getObjectAtPos(pos);
    if (brick != 0) {
        int brickIndex = static_cast<int>(reinterpret_cast<size_t>(brick)) - 1;
        if (!customBricklist_.at(brickIndex).hidden) {
            customBricklistStack_.push(customBricklist_);

            if (!customBricklist_.at(brickIndex).selected) {
                customBricklist_.at(brickIndex).selected = true;
                i = 0;
                bool inserted = false;
                while(i < selectedBricks_.size() && !inserted) {
                    if (customBricklist_.at(brickIndex).brickId > selectedBricks_.at(i).brickId) {
                        selectedBricks_.insert(selectedBricks_.begin() + i, customBricklist_.at(brickIndex));
                        inserted = true;
                    }
                    else{
                    ++i;
                    }
                }
                if (!inserted) {
                    selectedBricks_.insert(selectedBricks_.end(), customBricklist_.at(brickIndex));
                }
            }
            else{
                if (selectedBricks_.size() == 1) {
                    for (i = 0; i < customBricklist_.size(); ++i) {
                        customBricklist_.at(i).selected = true;
                    }
                    selectedBricks_.clear();
                    selectedBricks_ = customBricklist_;
                }
                else{
                    customBricklist_.at(brickIndex).selected = false;
                    i = 0;
                    bool del = false;
                    while(!del) {
                        if (selectedBricks_.at(i).brickId == brickIndex) {
                            selectedBricks_.erase(selectedBricks_.begin() + i);
                            del = true;
                        }
                        else{
                        ++i;
                        }
                    }
                }
            }
            if (selectedBricks_.size() == 1)
                translation_.set(selectedBricks_.at(0).translation);

            else
                translation_.set(findMiddle(selectedBricks_));

            updatePlaneWidgets();
            invalidate();
        }
    }
    else{
        customBricklistStack_.push(customBricklist_);
        for (i = 0; i < selectedBricks_.size(); ++i) {
            customBricklist_.at(selectedBricks_.at(i).brickId).selected = false;
        }
        selectedBricks_.clear();
        translation_.set(tgt::vec3(0.f, 0.f, 0.f));
        updatePlaneWidgets();
        invalidate();
    }
}

void ExplosionProxyGeometry::onEvent(tgt::Event* e) {
    tgt::MouseEvent* me = dynamic_cast<tgt::MouseEvent*>(e);
    if(me)
        pickingBuffer_.resize(me->viewport());

    RenderProcessor::onEvent(e);
}

void ExplosionProxyGeometry::onTranslateSelectedBricksEvent(tgt::MouseEvent * e) {
    if (!pickingBuffer_.hasRenderTarget()) {
        return;
    }

    tgt::ivec2 pos = tgt::ivec2(e->coord().x, canvasSize_.get().y - e->coord().y);

    if  (e->action() & tgt::MouseEvent::PRESSED) {
        const void* brick = idManager_.getObjectAtPos(pos);
            if (brick != 0) {
                customBricklistStack_.push(customBricklist_);
                int brickIndex = static_cast<int>(reinterpret_cast<size_t>(brick)) - 1;

                if (customBricklist_.at(brickIndex).selected) {
                    e->accept();
                    isClicked_ = true;
                    startCoord_.x = e->coord().x;
                    startCoord_.y = e->coord().y;
                    startPos_ = translation_.get();
                }
            }
        return;
    }

    if (e->action() & tgt::MouseEvent::MOTION) {
        if (isClicked_) {
            e->accept();

            LGL_ERROR;
            GLint deltaX, deltaY;

            GLint viewport[4];
            GLdouble modelview[16];
            GLdouble projection[16];
            GLdouble winX, winY, winZ;
            GLdouble posX, posY, posZ;

            deltaX = e->coord().x - startCoord_.x;
            deltaY = startCoord_.y - e->coord().y;

            tgt::mat4 projection_tgt = camera_.get().getProjectionMatrix(idManager_.getRenderTarget()->getSize());
            tgt::mat4 modelview_tgt = camera_.get().getViewMatrix();
            for (int i = 0; i < 4; ++i) {
                modelview[i+0]   = modelview_tgt[i].x;
                modelview[i+4]   = modelview_tgt[i].y;
                modelview[i+8]   = modelview_tgt[i].z;
                modelview[i+12]  = modelview_tgt[i].w;
                projection[i+0]  = projection_tgt[i].x;
                projection[i+4]  = projection_tgt[i].y;
                projection[i+8]  = projection_tgt[i].z;
                projection[i+12] = projection_tgt[i].w;
            }

            viewport[0] = 0;
            viewport[1] = 0;
            viewport[2] = static_cast<GLint>(idManager_.getRenderTarget()->getSize().x);
            viewport[3] = static_cast<GLint>(idManager_.getRenderTarget()->getSize().y);

            posX = startPos_.x;
            posY = startPos_.y;
            posZ = startPos_.z;

            LGL_ERROR;
            tgt::gluProject(posX, posY,posZ,modelview,projection, viewport,&winX, &winY, &winZ);

            winX = winX + deltaX;
            winY = winY + deltaY;

            LGL_ERROR;
            tgt::gluUnProject( winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ);
            LGL_ERROR;

            translation_.set(tgt::vec3(static_cast<float>(posX), static_cast<float>(posY), static_cast<float>(posZ)));

            LGL_ERROR;
            invalidate();
            LGL_ERROR;
        }
        return;
    }

    if (e->action() & tgt::MouseEvent::RELEASED) {
        if (isClicked_) {
            e->accept();
            isClicked_ = false;
            registerForSelecting(customBricklist_);
            invalidate();
        }
        return;
    }
}

void ExplosionProxyGeometry::canvasSizeChanged() {
    //RenderProcessor::portResized(0, canvasSize_.get());
}


} // namespace
