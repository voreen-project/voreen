/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2014 University of Muenster, Germany.                        *
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

//header file
#include "cmhalosubstructurerenderer.h"
//needed headers (used in process())
#include "tgt/textureunit.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include <assert.h>
#include <random>

//we are in namespace voreen
namespace voreen {

CMHaloSubstructureRenderer::CMHaloSubstructureRenderer()
    : RenderProcessor()
    , outport_(Port::OUTPORT, "outport","Modified Image", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , inport_(Port::INPORT, "halohandle.output", "Halo Data Output")
    , radiusDivProp_("radiusDivProp", "Radius divider" ,1000.0f ,1.0f, 100000.f)
    , hostHaloRingThicknessProp_("hostHaloRingThicknessProp_", "Relative size of host halo ring" , 0.1f, 0.0f, 1.0f)
    , camera_("camera", "Camera")
    , handleCameraMovement_("handleCameraMovement_", "Handle camera movement", true)
    , initialZoomDistanceProp_("initialZoomDistanceProp", "Distance to Halo in Focus mode", 2.50f, 1.0f, 100.0f)
    , animationDurationProp_("animationDurationProp", "Duration of camera movement in ms", 1000, 1, 10000, Processor::INVALID_RESULT, NumericProperty<int>::STATIC, Property::LOD_DEBUG)
    , shaderProp_("shaderProp", "HaloShader", "cmhalosubstructurerenderer.frag", "cmhalosubstructurerenderer.vert", "cmhalosubstructurerenderer.geom" )
    , shaderPropDVC_("shaderPropDVC", "HaloShaderDVC", "cmhalosubstructurerenderer_direct_velocity_color.frag", "cmhalosubstructurerenderer.vert", "cmhalosubstructurerenderer.geom" )
    , shaderPropSelection_("shaderPropSelection", "HaloSelectionShader", "cmhalosubstructureselectionrenderer.frag", "cmhalosubstructureselectionrenderer.vert", "cmhalosubstructureselectionrenderer.geom" )
    , transFunc_("transFunc", "Transfer function")
    , renderMode_("renderMode", "Mode for rendering")
    , useAlpha_("useAlpha", "Use alpha", false)
    , alphaFactor_("alphaFactor", "Alpha Factor", 1.0f, 0.0f, 1.0f)
    , selectedHaloIDProp_("selectedHaloIDProp", "ID of focused halo", 0, -1, 1000000)
    , mouseOverHaloIDProp_("mouseOverHaloIDProp", "ID of hovered over halo", CMMergerTree::NO_HALO_ID, -1, 1000000, Processor::INVALID_RESULT, NumericProperty<int>::STATIC, Property::LOD_DEBUG)
    , timeStep_("timeStep", "Time Step", 0.0f, 0.0f, 624.0f)
    , displayTimeStep_("displayTimeStep", "currently displayed time step", 0.0f, 0.0f, 624.0f)
    , timeAnimator_(60)
    , cameraAnimator_(60)
    , vbo_(0)
    , vao_(0)
    , selectionVao_(0)
    , bodyValues_(nullptr)
{

    selectionManager_ = new CMSelectionManager("selectionHandler", "Selection", tgt::ivec2(256, 256), &selectedHaloIDProp_, &mouseOverHaloIDProp_, CMMergerTree::NO_HALO_ID);
    addInteractionHandler(selectionManager_);
    cameraHandler_ = new CameraInteractionHandler("cameraHandler", "Camera", &camera_);
    addInteractionHandler(cameraHandler_);
    addPort(outport_);
    addPort(inport_);
    //register properties

    addProperty(renderMode_);
    addProperty(radiusDivProp_);
    addProperty(hostHaloRingThicknessProp_);
    addProperty(transFunc_);

    addProperty(useAlpha_);
    addProperty(alphaFactor_);

    addProperty(selectedHaloIDProp_);
    addProperty(mouseOverHaloIDProp_);
    addProperty(timeStep_);
    addProperty(displayTimeStep_);

    addProperty(camera_);
    addProperty(handleCameraMovement_);
    addProperty(initialZoomDistanceProp_);
    addProperty(shaderProp_);
    addProperty(shaderPropSelection_);
    addProperty(shaderPropDVC_);


    renderMode_.addOption("mass", "Halo mass", MASS_TRANSFER_FUNC);
    renderMode_.addOption("directVelocityColor",  "Direct velocity Color" , DIRECT_VELOCITY_COLOR);

    ON_CHANGE(renderMode_, CMHaloSubstructureRenderer, changedRenderMode);
    ON_CHANGE(useAlpha_, CMHaloSubstructureRenderer, changeUseAlpha);
    inport_.onNewData(MemberFunctionCallback<CMHaloSubstructureRenderer>(this, &CMHaloSubstructureRenderer::inputDataChanged));
    ON_CHANGE(timeStep_, CMHaloSubstructureRenderer, timeStepChanged);
    ON_CHANGE(displayTimeStep_, CMHaloSubstructureRenderer, setupBuffers);
    ON_CHANGE(selectedHaloIDProp_, CMHaloSubstructureRenderer, selectedHaloChanged);
    ON_CHANGE(radiusDivProp_, CMHaloSubstructureRenderer, setupBuffers);

    changedRenderMode();
    changeUseAlpha();
}

CMHaloSubstructureRenderer::~CMHaloSubstructureRenderer(){
    delete cameraHandler_;
    delete selectionManager_;
}

void CMHaloSubstructureRenderer::initialize() {
    // call superclass function first
    RenderProcessor::initialize();

    //move to constructor?
    camera_.setTrackballCenterBehaviour(CameraProperty::CAMSHIFT);

    // load fragment shader 'sample.frag'
    shaderProp_.rebuild();
    shaderPropSelection_.rebuild();
    shaderPropDVC_.rebuild();

    selectionManager_->initialize();
    //Fill texture with zeros
    selectionManager_->resize(outport_.getSize());
    selectionManager_->activate();
    selectionManager_->clear();
    selectionManager_->deactivate();
}

void CMHaloSubstructureRenderer::deinitialize() {


    // call superclass function last
    RenderProcessor::deinitialize();

    selectionManager_->deinitialize();
}


void CMHaloSubstructureRenderer::selectedHaloChanged() {
    if(handleCameraMovement_.get()) {
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
    setupBuffers();
}
void CMHaloSubstructureRenderer::timeStepChanged() {
    float timeStepDiff = fabs(timeStep_.get() - displayTimeStep_.get());
    if(timeStepDiff!=0) {
        std::vector<CMPropertyInterpolator*> props;
        props.push_back(new CMNumericPropertyInterpolator<float>(CMPropertyInterpolator::LINEAR, &displayTimeStep_, timeStep_.get()));
        timeAnimator_.animate(timeStepDiff*animationDurationProp_.get(), props);
    }
}
void CMHaloSubstructureRenderer::inputDataChanged() {
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

float CMHaloSubstructureRenderer::getSelectedBodyValue(const CMHalo* h) const {
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
void CMHaloSubstructureRenderer::changedRenderMode(){
    RenderMode rendermode = renderMode_.getValue();

    if (rendermode != DIRECT_VELOCITY_COLOR){
        transFunc_.setVisibleFlag(true);
    }else{
        transFunc_.setVisibleFlag(false);
    }
    setupBuffers();
}
void CMHaloSubstructureRenderer::changeUseAlpha(){
    if (useAlpha_.get()){
        alphaFactor_.setVisibleFlag(true);
    }else{
        alphaFactor_.setVisibleFlag(false);
    }
}
void CMHaloSubstructureRenderer::collectHalos(const CMHalo* hostHalo, std::vector<HaloVertexLayout>& verts) const{
    HaloVertexLayout l;
    l.pos = hostHalo->pos;
    l.vel = hostHalo->velocity;
    l.angularMomenta = tgt::normalize(hostHalo->angularMomenta);
    l.spinParameter = hostHalo->spinParameter;
    l.radius = hostHalo->radius/radiusDivProp_.get();
    //After a click on the ring we want to focus the parent to go "back"
    if(hostHalo->hostID != CMMergerTree::NO_HALO_ID) {
       l.haloID = hostHalo->hostID;
    } else {
       l.haloID = hostHalo->ID;
    }
    l.bodyValue = getSelectedBodyValue(hostHalo);

    verts.push_back(l);

    const CMHalo* currentSatellite = hostHalo->satellite();
    while(currentSatellite) {
        if(currentSatellite->scale != hostHalo->scale) {
            LWARNING("satellite halo's (OID=" + std::to_string(currentSatellite->origID) + ")scale ("
                    + std::to_string(currentSatellite->scale) + ") differs from host's (OID="
                    + std::to_string(hostHalo->origID) + ") (" + std::to_string(hostHalo->scale) + ").");
            break;
        } else {
            HaloVertexLayout h;
            h.pos = currentSatellite->pos;
            h.vel = currentSatellite->velocity;
            h.angularMomenta = tgt::normalize(currentSatellite->angularMomenta);
            h.spinParameter = currentSatellite->spinParameter;
            h.radius = currentSatellite->radius/radiusDivProp_.get();
            h.haloID = currentSatellite->ID;
            h.bodyValue = getSelectedBodyValue(currentSatellite);

            verts.push_back(h);
            currentSatellite = currentSatellite->siblingSatellite();
        }
    }
}

void CMHaloSubstructureRenderer::setupBuffers(){
    if (vao_) glDeleteVertexArrays(1, &vao_);
    if (vbo_) glDeleteBuffers(1, &vbo_);
    if (selectionVao_) glDeleteVertexArrays(1, &selectionVao_);


    if(!inport_.getData()) {
        vao_ = 0;
        vbo_ = 0;
        selectionVao_ = 0;
        return;
    }
    const CMTimeStepHaloList* timeSlicePointer = inport_.getData()->haloDataAt(displayTimeStep_.get());
    if(!timeSlicePointer) {
        vao_ = 0;
        vbo_ = 0;
        selectionVao_ = 0;
        return;
    }
    tgt::svec3 dataVolumeDim = tgt::svec3(timeSlicePointer->size(), 1, 1);

    VolumeRAM_Float* bodyvram     = new VolumeRAM_Float(dataVolumeDim);

    std::vector<HaloVertexLayout> verts;
    std::vector<GLushort> satelliteLinks;
    size_t i=0;
    for(auto h = timeSlicePointer->halosBegin(); h != timeSlicePointer->halosEnd(); ++h, ++i){
        bodyvram->setVoxelNormalized(getSelectedBodyValue(h), i, 0, 0);
        //Halo is a root host
        if(h->ID==selectedHaloIDProp_.get()) {
           collectHalos(h, verts);
        }
    }
    delete(bodyValues_);
    bodyValues_ = new Volume(bodyvram, tgt::vec3::zero, tgt::vec3::zero);
    transFunc_.setVolume(bodyValues_);

    glGenVertexArrays(1, &vao_);
    glGenBuffers(1, &vbo_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(HaloVertexLayout)*verts.size(), verts.data(), GL_STATIC_DRAW);

    glBindVertexArray(vao_);

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

    glGenVertexArrays(1, &selectionVao_);
    glBindVertexArray(selectionVao_);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(HaloVertexLayout), (void*)offsetof(HaloVertexLayout, pos));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(HaloVertexLayout), (void*)offsetof(HaloVertexLayout, radius));
    glEnableVertexAttribArray(2);
    glVertexAttribIPointer(2, 1, GL_UNSIGNED_INT, sizeof(HaloVertexLayout), (void*)offsetof(HaloVertexLayout, haloID));

    glBindVertexArray(0);

    particleCount_ = verts.size();
}

void CMHaloSubstructureRenderer::process() {
    if (!vao_ || !vbo_ || !selectionVao_) return;

    outport_.activateTarget();
    outport_.clearTarget();

    tgt::Camera cam = camera_.get();
    tgt::mat4 projectionMatrix = cam.getProjectionMatrix(outport_.getSize());
    tgt::mat4 viewMatrix = cam.getViewMatrix();

    processHalo(projectionMatrix, viewMatrix);


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
void CMHaloSubstructureRenderer::processHaloSelection(const tgt::mat4& projectionMatrix, const tgt::mat4& viewMatrix) {
    tgt::Shader* shader;
    shader = shaderPropSelection_.getShader();
    shader->activate();

    shader->setUniform("projectionMatrix", projectionMatrix);
    shader->setUniform("viewMatrix", viewMatrix);
    shader->setUniform("hostHaloRingThickness_", hostHaloRingThicknessProp_.get());

    glBindVertexArray(selectionVao_);
    glDisable(GL_BLEND);

    glDrawArrays(GL_POINTS, 0, particleCount_);

    glBindVertexArray(0);

    shader->deactivate();
    LGL_ERROR;
    assert(glGetError() == GL_NO_ERROR);
}
void CMHaloSubstructureRenderer::processHalo(const tgt::mat4& projectionMatrix, const tgt::mat4& viewMatrix) {
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
    shader->setUniform("alphaFactor_", alphaFactor_.get());
    shader->setUniform("hostHaloRingThickness_", hostHaloRingThicknessProp_.get());

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
    glDrawArrays(GL_POINTS, 0, particleCount_);
    bool useAlpha = useAlpha_.get();
    if (useAlpha){
        glEnable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ONE);
    }else{
        glDisable(GL_BLEND);
    }

    glBindVertexArray(0);


    glDisable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ZERO);
    shader->deactivate();
    LGL_ERROR;
    assert(glGetError() == GL_NO_ERROR);
}

} // namespace
