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

//header file
#include "cosmologyhalorenderer.h"
//needed headers (used in process())
#include "tgt/textureunit.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include <assert.h>
#include <random>

//we are in namespace voreen
namespace voreen {

CosmologyHaloRenderer::CosmologyHaloRenderer()
    : RenderProcessor()
    , outport_(Port::OUTPORT, "outport","Modified Image", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , inport_(Port::INPORT, "halohandle.output", "Halo Data Output")
	, enable_("enabled", "Enable Halo Rendering", true)
    , radiusDivProp_("radiusDivProp", "Radius divider" ,100.0f ,1.0f, 100000.0f)
    , camera_("camera", "Camera")
    , handleCameraMovement_("handleCameraMovement_", "Handle camera movement", true)
    , focusMode_("focusMode", "Focus of rendering")
    , initialZoomDistanceProp_("initialZoomDistanceProp", "Distance to Halo in Focus mode", 25.0f, 1.0f, 100.0f)
    , animationDurationProp_("animationDurationProp", "Duration of camera movement in ms", 1000, 1, 10000, Processor::INVALID_RESULT, NumericProperty<int>::STATIC, Property::LOD_DEBUG)
    , shaderProp_("shaderProp", "ASDSAF", "halorenderer.frag", "halorenderer.vert", "halorenderer.geom" )
    , shaderPropSelection_("shaderPropSelection", "ASDSAF", "haloselectionrenderer.frag", "haloselectionrenderer.vert", "haloselectionrenderer.geom" )
    , shaderPropDVC_("shaderPropDVC", "DVCASDSAF", "halorenderer_direct_velocity_color.frag", "halorenderer.vert", "halorenderer.geom" )
    , shaderPropSatelliteLinks_("shaderPropSatelliteLinks", "SatLinkShader", "satellitelinks.frag", "satellitelinks.vert", "satellitelinks.geom" )
    , rotationAxisProp_("rotationAxis", "RotationsAxis", "rotation_axis_color.frag", "rotationaxis.vert", "rotationaxis.geom" )
    , rotationOrbitProp_("rotationOrbit", "RotationsOrbit", "rotation_orbit_color.frag", "rotationorbit.vert", "rotationorbit.geom" )
    , transFunc_("transFunc", "Misc TF (body)")
    , axisTransFunc_("axisTransFunc", "Velocity TF (axis)")
    , orbitTransFunc_("orbitTransFunc", "Spin TF (orbit)")
    , renderMode_("renderMode", "Mode for rendering")
    , useAlpha_("useAlpha", "Use alpha", false)
    , alphaFactor_("alphaFactor", "Halo alpha factor", 1.0f, 0.0f, 1.0f)
    , previewHaloAlphaFactor_("previewHaloAlphaFactor", "Preview Halo Alpha factor", 0.5f, 0.0f, 1.0f)
    , selectedHaloIDProp_("selectedHaloIDProp", "ID of focused halo", 0, -1, 1000000)
    , mouseOverHaloIDProp_("mouseOverHaloIDProp", "ID of hovered over halo", CMMergerTree::NO_HALO_ID, -1, 1000000, Processor::INVALID_RESULT, NumericProperty<int>::STATIC, Property::LOD_DEBUG)
    , timeStep_("timeStep", "Time Step", 0.0f, 0.0f, 624.0f)
    , overviewFocusProp_("overviewFocusProp_", "Focus for overView mode", tgt::vec3(30, 30, 30), tgt::vec3(0,0,0), tgt::vec3(100, 100, 100))
    , overviewZoomProp_("overviewZoomProp", "Zoom for overview mode", 100.0f, 0.0f, 1000.0f)
    , satelliteLinkWidthProp_("satelliteLinkWidthProp", "Width of satellite links", 0.01f, 0.0f, 1.0f)
    , satelliteLinkDistanceProp_("satelliteLinkDistanceProp", "Distance for satellite links", 1.5f, 1.0f, 10.0f)
    , transferFunctionVisibleProp_("transferFunctionVisibleProp", "Main Transfer function visible", true)
    , bodyUnitProp_("bodyUnitProp", "Unit of body data displayed", "")
    , cameraAnimator_(60)
    , timeAnimator_(60)
    , vbo_(0)
    , vao_(0)
    , selectionVao_(0)
    , linkIbo_(0)
    , velocityValues_(nullptr)
    , bodyValues_(nullptr)
    , spinValues_(nullptr)
{

    selectionManager_ = new CMSelectionManager("selectionHandler", "Selection", tgt::ivec2(256, 256), &selectedHaloIDProp_, &mouseOverHaloIDProp_, CMMergerTree::NO_HALO_ID);
    addInteractionHandler(selectionManager_);
    cameraHandler_ = new CameraInteractionHandler("cameraHandler", "Camera", &camera_);
    addInteractionHandler(cameraHandler_);
    addPort(outport_);
    addPort(inport_);
    //register properties

    addProperty(enable_);
    addProperty(renderMode_);
    addProperty(radiusDivProp_);
    addProperty(transFunc_);
    addProperty(axisTransFunc_);
    addProperty(orbitTransFunc_);

    addProperty(useAlpha_);
    addProperty(alphaFactor_);
    addProperty(previewHaloAlphaFactor_);

    addProperty(selectedHaloIDProp_);
    addProperty(mouseOverHaloIDProp_);
    addProperty(timeStep_);

    addProperty(camera_);
    addProperty(handleCameraMovement_);
    addProperty(focusMode_);
    addProperty(initialZoomDistanceProp_);
    addProperty(animationDurationProp_);
    addProperty(shaderProp_);
    addProperty(shaderPropSelection_);
    // addProperty(shaderPropRotation_);
    addProperty(shaderPropDVC_);
    addProperty(shaderPropSatelliteLinks_);
    addProperty(rotationAxisProp_);
    addProperty(rotationOrbitProp_);
    addProperty(overviewFocusProp_);
    addProperty(overviewZoomProp_);
    addProperty(satelliteLinkWidthProp_);
    addProperty(satelliteLinkDistanceProp_);
    addProperty(transferFunctionVisibleProp_);
    addProperty(bodyUnitProp_);


    renderMode_.addOption("mass", "Halo mass", MASS_TRANSFER_FUNC);
    renderMode_.addOption("directVelocityColor",  "Direct velocity Color" , DIRECT_VELOCITY_COLOR);

    focusMode_.addOption("halo", "Focus Halo", HALO);
    focusMode_.addOption("overview",  "Overview" , OVERVIEW);


    ON_CHANGE(renderMode_, CosmologyHaloRenderer, changedRenderMode);
    ON_CHANGE(useAlpha_, CosmologyHaloRenderer, changeUseAlpha);
    inport_.onNewData(MemberFunctionCallback<CosmologyHaloRenderer>(this, &CosmologyHaloRenderer::inputDataChanged));
    ON_CHANGE(timeStep_, CosmologyHaloRenderer, timeStepChanged);
    ON_CHANGE(selectedHaloIDProp_, CosmologyHaloRenderer, selectedHaloChanged);
    ON_CHANGE(mouseOverHaloIDProp_, CosmologyHaloRenderer, mouseOverHaloChanged);
    ON_CHANGE(focusMode_, CosmologyHaloRenderer, focusModeChanged);
    ON_CHANGE(radiusDivProp_, CosmologyHaloRenderer, setupBuffers);

    changedRenderMode();
    changeUseAlpha();
}

CosmologyHaloRenderer::~CosmologyHaloRenderer(){
    delete cameraHandler_;
    delete selectionManager_;
}

void CosmologyHaloRenderer::initialize() {
    // call superclass function first
    RenderProcessor::initialize();

    //move to constructor?
    camera_.setTrackballCenterBehaviour(CameraProperty::CAMSHIFT);

    // load fragment shader 'sample.frag'
    shaderProp_.rebuild();
    shaderPropSelection_.rebuild();
    //shaderPropRotation_.rebuild();
    shaderPropDVC_.rebuild();
    shaderPropSatelliteLinks_.rebuild();
    rotationAxisProp_.rebuild();
    rotationOrbitProp_.rebuild();

    selectionManager_->initialize();

    //Fill texture with zeros
    selectionManager_->resize(outport_.getSize());
    selectionManager_->activate();
    selectionManager_->clear();
    selectionManager_->deactivate();

	glGenVertexArrays(1, &vao_);
	glGenBuffers(1, &vbo_);
	glGenVertexArrays(1, &selectionVao_);
	glGenBuffers(1, &linkIbo_);

	glBindVertexArray(vao_);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_);

	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(HaloVertexLayout), (void*)offsetof(HaloVertexLayout, pos));
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(HaloVertexLayout), (void*)offsetof(HaloVertexLayout, vel));
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(HaloVertexLayout), (void*)offsetof(HaloVertexLayout, angularMomenta));
	glEnableVertexAttribArray(3);
	glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(HaloVertexLayout), (void*)offsetof(HaloVertexLayout, spinParameter));
	glEnableVertexAttribArray(4);
	glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, sizeof(HaloVertexLayout), (void*)offsetof(HaloVertexLayout, radius));
	glEnableVertexAttribArray(5);
	glVertexAttribPointer(5, 1, GL_FLOAT, GL_FALSE, sizeof(HaloVertexLayout), (void*)offsetof(HaloVertexLayout, bodyValue));


	glBindVertexArray(selectionVao_);

	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(HaloVertexLayout), (void*)offsetof(HaloVertexLayout, pos));
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(HaloVertexLayout), (void*)offsetof(HaloVertexLayout, radius));
	glEnableVertexAttribArray(2);
	glVertexAttribIPointer(2, 1, GL_UNSIGNED_INT, sizeof(HaloVertexLayout), (void*)offsetof(HaloVertexLayout, haloID));

	glBindVertexArray(0);

}

void CosmologyHaloRenderer::deinitialize() {
    delete(velocityValues_);
    delete(bodyValues_);
    delete(spinValues_);
    selectionManager_->deinitialize();

	glDeleteVertexArrays(1, &vao_);
	glDeleteBuffers(1, &vbo_);
	glDeleteVertexArrays(1, &selectionVao_);
	glDeleteBuffers(1, &linkIbo_);

    // call superclass function last
    RenderProcessor::deinitialize();
}


void CosmologyHaloRenderer::selectedHaloChanged() {
    if(focusMode_.getValue() == FocusMode::HALO && handleCameraMovement_.get()) {
        const CMMergerTree* tree = inport_.getData();
        if(!tree) {
            return;
        }
        const CMHalo* selectedHalo = tree->haloByID(selectedHaloIDProp_.get());
        if(!selectedHalo) {
            return;
        }
        std::vector<CMPropertyInterpolator*> props;
        props.push_back(new CMCameraPropertyInterpolator(&camera_, selectedHalo->pos, initialZoomDistanceProp_.get()*selectedHalo->radius/radiusDivProp_.get()));
        cameraAnimator_.animate(animationDurationProp_.get(), props);
    }
}
void CosmologyHaloRenderer::mouseOverHaloChanged() {
    if(!vbo_) {
        return;
    }
    const CMMergerTree* tree = inport_.getData();
    if(!tree) {
        return;
    }
    const CMHalo* mouseOverHalo = tree->haloByID(mouseOverHaloIDProp_.get());
    if(!mouseOverHalo) {
        return;
    }
    HaloVertexLayout l;
    l.pos = mouseOverHalo->pos;
    l.vel = mouseOverHalo->velocity;
    l.angularMomenta = tgt::normalize(mouseOverHalo->angularMomenta);
    l.spinParameter = mouseOverHalo->spinParameter;
    l.radius = mouseOverHalo->radius/radiusDivProp_.get();
    l.haloID = mouseOverHalo->ID;
    l.bodyValue = getSelectedBodyValue(mouseOverHalo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferSubData(GL_ARRAY_BUFFER, sizeof(HaloVertexLayout)*particleCount_, sizeof(HaloVertexLayout), &l);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

}
void CosmologyHaloRenderer::timeStepChanged() {
    setupBuffers();
}
void CosmologyHaloRenderer::focusModeChanged() {
    FocusMode mode = focusMode_.getValue();
    if(handleCameraMovement_.get()) {
       if(mode == FocusMode::HALO) {
           const CMMergerTree* tree = inport_.getData();
           if(!tree) {
               return;
           }
           const CMHalo* selectedHalo = tree->haloByID(selectedHaloIDProp_.get());
           if(!selectedHalo) {
               return;
           }
           std::vector<CMPropertyInterpolator*> props;
           props.push_back(new CMCameraPropertyInterpolator(&camera_, selectedHalo->pos, 10*selectedHalo->radius/radiusDivProp_.get()));
           cameraAnimator_.animate(animationDurationProp_.get(), props);
       } else if(mode == FocusMode::OVERVIEW) {
           std::vector<CMPropertyInterpolator*> props;
           props.push_back(new CMCameraPropertyInterpolator(&camera_, overviewFocusProp_.get(), overviewZoomProp_.get()));
           cameraAnimator_.animate(animationDurationProp_.get(), props);
       }
    }
}
void CosmologyHaloRenderer::inputDataChanged() {
    const CMMergerTree* tree = inport_.getData();
    if(!tree) {
        return;
    }
    renderMode_.removeOption("blackHoleMass");
    renderMode_.removeOption("blackHoleSpin");
    renderMode_.removeOption("spheroidRadius");
    renderMode_.removeOption("spheroidMassGas");
    renderMode_.removeOption("spheroidVelocity");
    renderMode_.removeOption("diskRadius");
    renderMode_.removeOption("diskMassGas");
    renderMode_.removeOption("diskVelocity");
    if(tree->containsGalacticusData()) {
        renderMode_.addOption("blackHoleMass", "Black hole mass", BLACK_HOLE_MASS_TRANSFER_FUNC);
        renderMode_.addOption("blackHoleSpin", "Black hole spin", BLACK_HOLE_SPIN_TRANSFER_FUNC);
        renderMode_.addOption("spheroidRadius", "Spheroid radius", SPHEROID_RADIUS_TRANSFER_FUNC);
        renderMode_.addOption("spheroidMassGas", "Spheroid mass gas", SPHEROID_MASS_GAS_TRANSFER_FUNC);
        renderMode_.addOption("spheroidVelocity", "Spheroid velocity", SPHEROID_VELOCITY_TRANSFER_FUNC);
        renderMode_.addOption("diskRadius", "Disk radius", DISK_RADIUS_TRANSFER_FUNC);
        renderMode_.addOption("diskMassGas", "Disk mass gas", DISK_MASS_GAS_TRANSFER_FUNC);
        renderMode_.addOption("diskVelocity", "Disk velocity", DISK_VELOCITY_TRANSFER_FUNC);
    }
    setupBuffers();
}
float CosmologyHaloRenderer::getSelectedBodyValue(const CMHalo* h) const {
    float value = -1;
    switch(renderMode_.getValue()) {
        case MASS_TRANSFER_FUNC:
            value = h->mass;
            break;
        case BLACK_HOLE_MASS_TRANSFER_FUNC:
            value = h->blackHoleMass;
            break;
        case BLACK_HOLE_SPIN_TRANSFER_FUNC:
            value = h->blackHoleSpin;
            break;
        case SPHEROID_RADIUS_TRANSFER_FUNC:
            value = h->spheroidRadius;
            break;
        case SPHEROID_MASS_GAS_TRANSFER_FUNC:
            value = h->spheroidMassGas;
            break;
        case SPHEROID_VELOCITY_TRANSFER_FUNC:
            value = h->spheroidVelocity;
            break;
        case DISK_RADIUS_TRANSFER_FUNC:
            value = h->diskRadius;
            break;
        case DISK_MASS_GAS_TRANSFER_FUNC:
            value = h->diskMassGas;
            break;
        case DISK_VELOCITY_TRANSFER_FUNC:
            value = h->diskVelocity;
            break;
    }
    return value;
}
void CosmologyHaloRenderer::changedRenderMode(){
    RenderMode rendermode = renderMode_.getValue();

    if (rendermode != DIRECT_VELOCITY_COLOR){
        transFunc_.setVisibleFlag(true);
        transferFunctionVisibleProp_.set(true);
    }else{
        transFunc_.setVisibleFlag(false);
        transferFunctionVisibleProp_.set(false);
    }
    switch(renderMode_.getValue()) {
        //Is this even correct?
        case MASS_TRANSFER_FUNC:
        case SPHEROID_MASS_GAS_TRANSFER_FUNC:
        case BLACK_HOLE_MASS_TRANSFER_FUNC:
        case DISK_MASS_GAS_TRANSFER_FUNC:
            bodyUnitProp_.set("Msun/h");
            break;
        case BLACK_HOLE_SPIN_TRANSFER_FUNC:
            bodyUnitProp_.set("Spin");
            break;
        case SPHEROID_RADIUS_TRANSFER_FUNC:
        case DISK_RADIUS_TRANSFER_FUNC:
            bodyUnitProp_.set("Mpc/h");
            break;
        case SPHEROID_VELOCITY_TRANSFER_FUNC:
        case DISK_VELOCITY_TRANSFER_FUNC:
            bodyUnitProp_.set("km/s");
            break;
    }
    //Read currently need values from mergerTree
    setupBuffers();
}

void CosmologyHaloRenderer::changeUseAlpha(){
    if (useAlpha_.get()){
        alphaFactor_.setVisibleFlag(true);
    }else{
        alphaFactor_.setVisibleFlag(false);
    }
}
void CosmologyHaloRenderer::collectHalosAndLinks(const CMHalo* hostHalo, int hostVertID, std::vector<HaloVertexLayout>& verts, std::vector<GLushort>& satelliteLinks) const {
    HaloVertexLayout l;
    l.pos = hostHalo->pos;
    l.vel = hostHalo->velocity;
    l.angularMomenta = tgt::normalize(hostHalo->angularMomenta);
    l.spinParameter = hostHalo->spinParameter;
    l.radius = hostHalo->radius/radiusDivProp_.get();
    l.haloID = hostHalo->ID;
    l.bodyValue = getSelectedBodyValue(hostHalo);

    GLushort ownVertID = verts.size();
    verts.push_back(l);

    if(hostVertID != -1) {
        satelliteLinks.push_back((GLushort)hostVertID);
        satelliteLinks.push_back(ownVertID);
    }

    const CMHalo* currentSatellite = hostHalo->satellite();
    while(currentSatellite) {
        if(currentSatellite->scale != hostHalo->scale) {
            LWARNING("satellite halo's (OID=" + std::to_string(currentSatellite->origID) + ")scale ("
                    + std::to_string(currentSatellite->scale) + ") differs from host's (OID="
                    + std::to_string(hostHalo->origID) + ") (" + std::to_string(hostHalo->scale) + ").");
            break;
        } else {
            collectHalosAndLinks(currentSatellite, ownVertID, verts, satelliteLinks);
            currentSatellite = currentSatellite->siblingSatellite();
        }
    }
}

void CosmologyHaloRenderer::setupBuffers(){

    if(!inport_.getData()) {
        vao_ = 0;
        vbo_ = 0;
        selectionVao_ = 0;
        linkIbo_ = 0;
        return;
    }
	float timeStepScaleFactor = (float)(((1.0l - (1.0l / 201.0l)) / 625.0l) * (timeStep_.get() + 1.0l)) + (1.0l / 201.0l);
	const CMTimeStepHaloList* timeSlicePointer = inport_.getData()->haloDataAt(timeStepScaleFactor);
    if(!timeSlicePointer) {
        vao_ = 0;
        vbo_ = 0;
        selectionVao_ = 0;
        linkIbo_ = 0;
        return;
    }

    tgt::svec3 dataVolumeDim = tgt::svec3(timeSlicePointer->size(), 1, 1);
    VolumeRAM_Float* velocityvram = new VolumeRAM_Float(dataVolumeDim);
    VolumeRAM_Float* bodyvram     = new VolumeRAM_Float(dataVolumeDim);
    VolumeRAM_Float* spinvram     = new VolumeRAM_Float(dataVolumeDim);

    std::vector<HaloVertexLayout> verts;
    std::vector<GLushort> satelliteLinks;
    size_t i=0;
    for(auto h = timeSlicePointer->halosBegin(); h != timeSlicePointer->halosEnd(); ++h, ++i){
        velocityvram->setVoxelNormalized(tgt::length(h->velocity), i, 0, 0);
        spinvram->setVoxelNormalized(h->spinParameter, i, 0, 0);
        bodyvram->setVoxelNormalized(getSelectedBodyValue(h), i, 0, 0);
        //Halo is a root host
        if(h->rootHostID==CMMergerTree::NO_HALO_ID) {
            collectHalosAndLinks(h, -1, verts, satelliteLinks);
        }
    }
    delete(velocityValues_);
    delete(bodyValues_);
    delete(spinValues_);
    velocityValues_ = new Volume(velocityvram, tgt::vec3::zero, tgt::vec3::zero);
    bodyValues_ = new Volume(bodyvram, tgt::vec3::zero, tgt::vec3::zero);
    spinValues_ = new Volume(spinvram, tgt::vec3::zero, tgt::vec3::zero);
    transFunc_.setVolume(bodyValues_);
    axisTransFunc_.setVolume(velocityValues_);
    orbitTransFunc_.setVolume(spinValues_);

    const CMHalo* mouseOverHalo = inport_.getData()->haloByID(mouseOverHaloIDProp_.get());
    HaloVertexLayout l;
    if(mouseOverHalo) {
        l.pos = mouseOverHalo->pos;
        l.vel = mouseOverHalo->velocity;
        l.angularMomenta = tgt::normalize(mouseOverHalo->angularMomenta);
        l.spinParameter = mouseOverHalo->spinParameter;
        l.radius = mouseOverHalo->radius/radiusDivProp_.get();
        l.haloID = mouseOverHalo->ID;
        l.bodyValue = getSelectedBodyValue(mouseOverHalo);
    }
    verts.push_back(l);

    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(HaloVertexLayout)*verts.size(), verts.data(), GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);


    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, linkIbo_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLushort)*satelliteLinks.size(), satelliteLinks.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    particleCount_ = verts.size() - 1; //The last halo in verts is the one hovered over
    linkCount_ = satelliteLinks.size();
}

void CosmologyHaloRenderer::process() {
    if (!vao_ || !vbo_ || !selectionVao_ || !particleCount_) return;

    outport_.activateTarget();
    outport_.clearTarget();

	if (timeStep_.get() < 460.51f || !enable_.get()) {
		outport_.deactivateTarget();
		return;
	}

    tgt::Camera cam = camera_.get();
    const tgt::mat4 projectionMatrix = cam.getProjectionMatrix(outport_.getSize());
    const tgt::mat4 viewMatrix = cam.getViewMatrix();
    bool useAlpha = useAlpha_.get();
    float alphaFactor = alphaFactor_.get();

    if(useAlpha) {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    } else {
        glDisable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ZERO);
    }
    if(alphaFactor>0) {

        processRotationAxis(projectionMatrix, viewMatrix, alphaFactor, 0, particleCount_);
        processRotationOrbit(projectionMatrix, viewMatrix, alphaFactor, 0, particleCount_);
        processHalo(projectionMatrix, viewMatrix, alphaFactor, 0, particleCount_);

        processSatelliteLinks(projectionMatrix, viewMatrix, alphaFactor*0.5f);
    }

    if(inport_.getData() && mouseOverHaloIDProp_.get() != CMMergerTree::NO_HALO_ID
            /*&& inport_.getData()->haloByID(mouseOverHaloIDProp_.get())->scale != timeStep_.get()*/) {
        float previewHaloAlphaFactor = previewHaloAlphaFactor_.get();

        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        processRotationAxis(projectionMatrix, viewMatrix, previewHaloAlphaFactor, particleCount_, 1);
        processRotationOrbit(projectionMatrix, viewMatrix, previewHaloAlphaFactor, particleCount_, 1);
        processHalo(projectionMatrix, viewMatrix, previewHaloAlphaFactor, particleCount_, 1);
    }
    glDisable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ZERO);


    // cleanup
    outport_.deactivateTarget();

    selectionManager_->resize(outport_.getSize());
    selectionManager_->activate();
    selectionManager_->clear();
    processHaloSelection(projectionMatrix, viewMatrix);
    selectionManager_->deactivate();

    tgt::TextureUnit::setZeroUnit();
    LGL_ERROR;
}
void CosmologyHaloRenderer::processHaloSelection(const tgt::mat4& projectionMatrix, const tgt::mat4& viewMatrix) {
    tgt::Shader* shader;
    shader = shaderPropSelection_.getShader();
    shader->activate();

    shader->setUniform("projectionMatrix", projectionMatrix);
    shader->setUniform("viewMatrix", viewMatrix);

    glBindVertexArray(selectionVao_);
    glDisable(GL_BLEND);

    glDrawArrays(GL_POINTS, 0, particleCount_);

    glBindVertexArray(0);

    shader->deactivate();
    LGL_ERROR;
    assert(glGetError() == GL_NO_ERROR);
}
void CosmologyHaloRenderer::processHalo(const tgt::mat4& projectionMatrix, const tgt::mat4& viewMatrix, float alphaFactor, GLint first, GLsizei count) {
    RenderMode rendermode = renderMode_.getValue();

    tgt::Shader* shader;
    if (rendermode != DIRECT_VELOCITY_COLOR){
        shader = shaderProp_.getShader();
    }else{
        shader = shaderPropDVC_.getShader();
    }
    shader->activate();

    shader->setUniform("projectionMatrix", projectionMatrix);
    shader->setUniform("viewMatrix", viewMatrix);

    shader->setUniform("alphaFactor_", alphaFactor);

    if (rendermode != DIRECT_VELOCITY_COLOR){
        //tgt::vec2 domain = transFunc_.get()->getDomain();
        //shader->setUniform("scale_", domain.y-domain.x);
        //shader->setUniform("offset_", domain.x);

        tgt::TextureUnit transferUnit;
        if (transFunc_.get()) {
            transferUnit.activate();
            transFunc_.get()->getTexture()->bind();
            transFunc_.get()->setUniform(shader, "transferFunc_", "transferFuncTex_", transferUnit.getUnitNumber());
        }
    }

    glBindVertexArray(vao_);
    glDrawArrays(GL_POINTS, first, count);
    glBindVertexArray(0);

    shader->deactivate();
    LGL_ERROR;
    assert(glGetError() == GL_NO_ERROR);
}

void CosmologyHaloRenderer::processRotationAxis(const tgt::mat4& projectionMatrix, const tgt::mat4& viewMatrix, float alphaFactor, GLint first, GLsizei count) {
    RenderMode rendermode = renderMode_.getValue();

    glLineWidth(2);
    tgt::Shader* shader;

    shader = rotationAxisProp_.getShader();
    shader->activate();

    shader->setUniform("projectionMatrix", projectionMatrix);
    shader->setUniform("viewMatrix", viewMatrix);
    shader->setUniform("alphaFactor_", alphaFactor);

    tgt::TextureUnit transferUnit;
    if (axisTransFunc_.get()) {
        transferUnit.activate();
        axisTransFunc_.get()->getTexture()->bind();
        axisTransFunc_.get()->setUniform(shader, "axisTransferFunc_", "axisTransferFuncTex_", transferUnit.getUnitNumber());
    }

    glBindVertexArray(vao_);
    glDrawArrays(GL_POINTS, first, count);
    glBindVertexArray(0);

    glLineWidth(1);
    shader->deactivate();
    assert(glGetError() == GL_NO_ERROR);
}

void CosmologyHaloRenderer::processRotationOrbit(const tgt::mat4& projectionMatrix, const tgt::mat4& viewMatrix, float alphaFactor, GLint first, GLsizei count) {

    RenderMode rendermode = renderMode_.getValue();


    outport_.activateTarget();

    glLineWidth(2);

    tgt::Shader* shader;

    shader = rotationOrbitProp_.getShader();
    shader->activate();

    shader->setUniform("projectionMatrix", projectionMatrix);
    shader->setUniform("viewMatrix", viewMatrix);
    shader->setUniform("alphaFactor_", alphaFactor);

    tgt::TextureUnit transferUnit;
    if (orbitTransFunc_.get()) {
        transferUnit.activate();
        orbitTransFunc_.get()->getTexture()->bind();
        orbitTransFunc_.get()->setUniform(shader, "orbitTransferFunc_", "orbitTransferFuncTex_", transferUnit.getUnitNumber());
    }

    glBindVertexArray(vao_);
    glDrawArrays(GL_POINTS, first, count);
    glBindVertexArray(0);

    glLineWidth(1);

    shader->deactivate();
    assert(glGetError() == GL_NO_ERROR);
}

void CosmologyHaloRenderer::processSatelliteLinks(const tgt::mat4& projectionMatrix, const tgt::mat4& viewMatrix, float alphaFactor) {
    tgt::Shader* linkshader = shaderPropSatelliteLinks_.getShader();
    linkshader->activate();

    linkshader->setUniform("projectionMatrix", projectionMatrix);
    linkshader->setUniform("viewMatrix", viewMatrix);
    linkshader->setUniform("baseWidth", satelliteLinkWidthProp_.get());
    linkshader->setUniform("spacing", satelliteLinkDistanceProp_.get());
    linkshader->setUniform("alphaFactor_", alphaFactor);

    glBindVertexArray(vao_);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, linkIbo_);
    glDrawElements(GL_LINES, linkCount_, GL_UNSIGNED_SHORT, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    linkshader->deactivate();
    LGL_ERROR;
}


} // namespace
