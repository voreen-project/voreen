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

#include "illuminationlineraycaster.h"

#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"
#include "tgt/textureunit.h"
#include "tgt/glmath.h"

using tgt::vec3;
using tgt::TextureUnit;

struct Vertex2D{
    float x, y;
    float s, t;
};

#define BUFFER_OFFSET(i) ((char *)NULL + (i))

namespace voreen {

IlluminationLineRaycaster::IlluminationLineRaycaster()
    : VolumeRaycaster()
    //, volumeInport_(Port::INPORT, "volumehandle.volumehandle")
    , geometryPort_(Port::INPORT, "proxgeometry.geometry")
    //, entryPort_(Port::INPORT, "image.entrypoints")
    //, exitPort_(Port::INPORT, "image.exitpoints")
#ifdef ERT_DOWNSAMPLE_SCHEME
    , ertPort_(Port::INPORT, "image.ertPort")
#endif
    //, outport1_(Port::OUTPORT, "image.output", "image.output", true)
    , lightportOne_(Port::OUTPORT, "image.lightportOne", "image.lightportOne", true, Processor::INVALID_PROGRAM, RenderPort::RENDERSIZE_DEFAULT, GL_RGBA32F_ARB)
    , lightportTwo_(Port::OUTPORT, "image.lightportTwo", "image.lightportTwo", true, Processor::INVALID_PROGRAM, RenderPort::RENDERSIZE_DEFAULT, GL_RGBA32F_ARB)
    , lightportOnePos_(Port::OUTPORT, "image.lightportOnePos", "image.lightportOnePos", true, Processor::INVALID_PROGRAM, RenderPort::RENDERSIZE_DEFAULT, GL_RGBA32F_ARB)
    , lightportTwoPos_(Port::OUTPORT, "image.lightportTwoPos", "image.lightportTwoPos", true, Processor::INVALID_PROGRAM, RenderPort::RENDERSIZE_DEFAULT, GL_RGBA32F_ARB)
    , propagationZonesPort_(Port::OUTPORT, "image.propagationZones", "image.propagationZones", true, Processor::INVALID_PROGRAM, RenderPort::RENDERSIZE_DEFAULT, GL_RGBA32F_ARB)
    , propagationErtPort_(Port::OUTPORT, "image.propagationErtPort", "image.propagationErtPort", true, Processor::INVALID_PROGRAM, RenderPort::RENDERSIZE_DEFAULT, GL_RGBA32F_ARB)
    , renderQuad_("renderQuad", "Render quad, thus disabling lighting (using GLSL 3.30)")
    , autoLightCacheSize_("autoLightCacheSize", "Use canvas size as light cache size.", true)
    , lightCacheWidth_("lightCacheWidth", "Width of Light Cache", 512, 1, 4096)
    , lightCacheHeight_("lightCacheHeight", "Height of Light Cache", 512, 1, 4096)
    , minLightDirLength_("minLightDirLength", "Minimum Allowed Length of Light Vector to Render Shadows with", 0.f, 0.f, 1.f)
    , diffCamLightDistance_("diffCamLightDistance", "The difference in light distance to focus point compared to camera", 0.f, -10.f, 10.f)
    , adaptiveNearPlane_("adaptiveNearPlane", "Adapt near plane i.e. light cache for optimal use.")
    , ertMod_("ertMod", "Perform line sweeping on ERT generated texture")
#ifdef ERT_DOWNSAMPLE_SCHEME
    , ertIncrement_("ertIncrement", "Specify increment of ert sweeping", 1, 1, 10)
    , ertTextureSizeWidth_("ertTextureSizeWidth", "Specify width of ert texture", 128, 1, 1024)
    , ertTextureSizeHeight_("ertTextureSizeHeight", "Specify height of ert texture", 128, 1, 1024)
#endif
    , twoPassRenderingOn_("twoPassRenderingOn", "Rendering using two pass when inside volume")
    , zoneRenderingOn_("zoneRenderingOn", "Rendering using propagation zones.")
    , zoneOnOff_("zoneOnOff", "Define which zone to be rendered or all of them.")
    , zoneEdgeBlendingFactor_("zoneEdgeBlendingFactor", "Factor that defines size of edge blending area.", 2, 1, 10)
    , accumulatedRayVis_("accumulatedRayVis", "Different lighting method saving only the resulting ray.")
    , colorContribution_("colorContribution", "Incorporate voxel color in shadow apperance.")
    , shadowBrightness_("shadowBrightness", "Brightness of shadow when color is incorporated.", 1.f, 0.f, 1.f)
    , interpolation_("interpolation", "Type of shadow interpolation.")
    , piecewise_("piecewise_", "Piecewise Integration On/Off", true)
    , adjustToProxyGeometry_("adjustToProxyGeometry", "Adjust light map based on proxy geometry size seen from camera.")
#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
    , vboOn_("vboOn", "Use VBO/IBO when drawing lines.")
    , drawTrianglesAsLines_("drawTrianglesAsLines", "Draw a triangle strip to represent a line")
    , linesOnlyAboveProxyGeometry_("linesOnlyAboveProxyGeometry", "Draw lines that just cover an enclosed area around the proxy geometry.", true)
    , visualizeProxySize_("visualizeProxySize", "Visualize the screen size of the proxy geometry.")
    , pingPong_("pingPong", "Turn ping-pong buffering on/off.", true)
    , outputLineToDebugRT_("outputLineToDebugRT", "Visibility from certain processed line. 0 = OFF", 512, 0, 2048)
#endif
    , lineWidth_("lineWidth", "Line width", 1, 1, 10)
    , glslVersion_("glslVersion", "GLSL Version")
    , raycastPrg_(0)
    , copyImagePrg_(0)
    //, blendPrg_(0)
    , propZonePrg_(0)
    , propErtPrg_(0)
    , transferFunc_("transferFunction", "Transfer Function", Processor::INVALID_RESULT)
    , albedo_("albedo", "Albedo", Processor::INVALID_RESULT)
    , surfaceness_("surfaceness", "Surfaceness", Processor::INVALID_RESULT)
    , ambientIntensity_("ambientIntensity", "Ambient intensity", 0.1f, 0.0f, 1.0f)
    , lightCamera_(tgt::Camera(tgt::vec3(0.f, -3.5f, -3.5f), tgt::vec3(0.f, 0.f, 0.f), tgt::vec3(0.f, 1.f, 0.f)))
    , camera_("camera", "Camera", tgt::Camera(vec3(0.f, 0.f, 3.5f), vec3(0.f, 0.f, 0.f), vec3(0.f, 1.f, 0.f)))
{

    addProperty(renderQuad_);

    addProperty(autoLightCacheSize_);
    addProperty(lightCacheWidth_);
    addProperty(lightCacheHeight_);
    addProperty(minLightDirLength_);
    addProperty(diffCamLightDistance_);

    addProperty(adaptiveNearPlane_);
    addProperty(ertMod_);
#ifdef ERT_DOWNSAMPLE_SCHEME
    addProperty(ertIncrement_);
    addProperty(ertTextureSizeWidth_);
    addProperty(ertTextureSizeHeight_);
#endif
    addProperty(twoPassRenderingOn_);
    addProperty(zoneRenderingOn_);

    zoneOnOff_.addOption("A", "ALL_ZONES", 0);
    zoneOnOff_.addOption("R", "ZONE_1(RED)", 1);
    zoneOnOff_.addOption("G", "ZONE_2(GREEN)", 2);
    zoneOnOff_.addOption("B", "ZONE_3(BLUE)", 3);
    zoneOnOff_.addOption("Y", "ZONE_4(YELLOW)", 4);
    zoneOnOff_.select("A");

    addProperty(zoneOnOff_);
    addProperty(zoneEdgeBlendingFactor_);

    addProperty(accumulatedRayVis_);
    addProperty(colorContribution_);
    addProperty(shadowBrightness_);

    interpolation_.addOption("NN", "NEAREST_NEIGHBOR", 0);
    interpolation_.addOption("BI", "BILINEAR_INTERPOLATION", 1);
    interpolation_.select("NN");

    addProperty(interpolation_);
    addProperty(piecewise_);
    addProperty(adjustToProxyGeometry_);
#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
    addProperty(vboOn_);
    addProperty(drawTrianglesAsLines_);
    addProperty(linesOnlyAboveProxyGeometry_);
    addProperty(visualizeProxySize_);
    addProperty(pingPong_);
    addProperty(outputLineToDebugRT_);
#endif
    addProperty(lineWidth_);

    glslVersion_.addOption("400", "GLSL 4.00", tgt::GpuCapabilities::GlVersion::SHADER_VERSION_400);
    glslVersion_.addOption("330", "GLSL 3.30", tgt::GpuCapabilities::GlVersion::SHADER_VERSION_330);
    glslVersion_.addOption("150", "GLSL 1.50", tgt::GpuCapabilities::GlVersion::SHADER_VERSION_150);
    glslVersion_.addOption("140", "GLSL 1.40", tgt::GpuCapabilities::GlVersion::SHADER_VERSION_140);
    glslVersion_.addOption("130", "GLSL 1.30", tgt::GpuCapabilities::GlVersion::SHADER_VERSION_130);
    glslVersion_.select("130");

    addProperty(glslVersion_);
    addProperty(lightPosition_);
    addProperty(lightDiffuse_);

    //addPort(volumeInport_);
    addPort(geometryPort_);
    //addPort(entryPort_);
    //addPort(exitPort_);
#ifdef ERT_DOWNSAMPLE_SCHEME
    addPort(ertPort_);
#endif
    //addPort(outport1_);

    addPrivateRenderPort(lightportOne_);
    addPrivateRenderPort(lightportTwo_);
    addPrivateRenderPort(lightportOnePos_);
    addPrivateRenderPort(lightportTwoPos_);
    addPrivateRenderPort(propagationZonesPort_);
    addPrivateRenderPort(propagationErtPort_);

    addProperty(transferFunc_);
    addProperty(camera_);
    camera_.setVisibleFlag(false);

    addProperty(ambientIntensity_);
    addProperty(albedo_);
    addProperty(surfaceness_);
}

IlluminationLineRaycaster::~IlluminationLineRaycaster() {
}

Processor* IlluminationLineRaycaster::create() const {
    return new IlluminationLineRaycaster();
}

void IlluminationLineRaycaster::initialize() {
    VolumeRaycaster::initialize();

    float lineWidth[2];
    glGetFloatv(GL_LINE_WIDTH_RANGE, lineWidth);
    lineWidth_.setMaxValue((int)lineWidth[1]);

    // load shaders
    raycastPrg_ = ShdrMgr.loadSeparate("passthrough.vert", "rc_image_based_illum.frag", generateHeader(), false);
    propZonePrg_ = ShdrMgr.loadSeparate("passthrough.vert", "define_propagation_zones.frag", generateHeader(), false);
    propErtPrg_ = ShdrMgr.loadSeparate("passthrough.vert", "define_ert_propagation_image.frag", generateHeader(), false);
    copyImagePrg_ = ShdrMgr.loadSeparate("passthrough.vert", "copyimage.frag",
        generateHeader() + "#define NO_DEPTH_TEX\n", false);

    quadList_ = 0;
    drawQuadIntoList();

    vboID_[0] = vboID_[1] = 0;
    iboID_ = 0;
    currentVBO_ = 0;

    lightPosition_.onChange(MemberFunctionCallback<IlluminationLineRaycaster>(this, &IlluminationLineRaycaster::updatePropagation));
    camera_.onChange(MemberFunctionCallback<IlluminationLineRaycaster>(this, &IlluminationLineRaycaster::updatePropagation));
    minLightDirLength_.onChange(MemberFunctionCallback<IlluminationLineRaycaster>(this, &IlluminationLineRaycaster::updatePropagation));
    diffCamLightDistance_.onChange(MemberFunctionCallback<IlluminationLineRaycaster>(this, &IlluminationLineRaycaster::updatePropagation));
    renderQuad_.onChange(MemberFunctionCallback<IlluminationLineRaycaster>(this, &IlluminationLineRaycaster::updateHeader));
    zoneRenderingOn_.onChange(MemberFunctionCallback<IlluminationLineRaycaster>(this, &IlluminationLineRaycaster::updateHeader));
    accumulatedRayVis_.onChange(MemberFunctionCallback<IlluminationLineRaycaster>(this, &IlluminationLineRaycaster::updateHeader));
    ertMod_.onChange(MemberFunctionCallback<IlluminationLineRaycaster>(this, &IlluminationLineRaycaster::updateHeader));
    ertMod_.onChange(MemberFunctionCallback<IlluminationLineRaycaster>(this, &IlluminationLineRaycaster::showAndHideProperties));
    colorContribution_.onChange(MemberFunctionCallback<IlluminationLineRaycaster>(this, &IlluminationLineRaycaster::updateHeader));
    colorContribution_.onChange(MemberFunctionCallback<IlluminationLineRaycaster>(this, &IlluminationLineRaycaster::showAndHideProperties));
    zoneRenderingOn_.onChange(MemberFunctionCallback<IlluminationLineRaycaster>(this, &IlluminationLineRaycaster::showAndHideProperties));
    autoLightCacheSize_.onChange(MemberFunctionCallback<IlluminationLineRaycaster>(this, &IlluminationLineRaycaster::showAndHideProperties));
    interpolation_.onChange(MemberFunctionCallback<IlluminationLineRaycaster>(this, &IlluminationLineRaycaster::updateHeader));
    piecewise_.onChange(MemberFunctionCallback<IlluminationLineRaycaster>(this, &IlluminationLineRaycaster::updateHeader));
    glslVersion_.onChange(MemberFunctionCallback<IlluminationLineRaycaster>(this, &IlluminationLineRaycaster::updateHeader));

    screenSizeChanged();

#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
    blendPrg_ = ShdrMgr.loadSeparate("passthrough.vert", "blendwithimage.frag",
        generateHeader(), false);

    adjustToProxyGeometry_.onChange(MemberFunctionCallback<IlluminationLineRaycaster>(this, &IlluminationLineRaycaster::showAndHideProperties));
    lineWidth_.onChange(MemberFunctionCallback<IlluminationLineRaycaster>(this, &IlluminationLineRaycaster::screenSizeChanged));
    outputLineToDebugRT_.set(outputLineToDebugRT_.getMaxValue());
#endif

    readWriteValues_ = tgt::ivec2(0, 1);

    plfFunc_ = &IlluminationLineRaycaster::postProcessLine;

    showAndHideProperties();

    updatePropagation();
}

void IlluminationLineRaycaster::deinitialize() {
    ShdrMgr.dispose(raycastPrg_);
    raycastPrg_ = 0;

    ShdrMgr.dispose(propZonePrg_);
    propZonePrg_ = 0;

    ShdrMgr.dispose(propErtPrg_);
    propErtPrg_ = 0;

    ShdrMgr.dispose(copyImagePrg_);
    copyImagePrg_ = 0;

    if(quadList_)
        glDeleteLists(quadList_, 1);

    if(!vboID_[0])
        glDeleteBuffers(2, &vboID_[0]);

    if(!iboID_)
        glDeleteBuffers(1, &iboID_);

#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
    ShdrMgr.dispose(blendPrg_);
    blendPrg_ = 0;
#endif

    LGL_ERROR;

    VolumeRaycaster::deinitialize();
}

void IlluminationLineRaycaster::beforeProcess() {
    VolumeRaycaster::beforeProcess();

    tgtAssert(volumeInport_.getData()->getRepresentation<VolumeRAM>(), "no volume");

    // assign volume to transfer function
    transferFunc_.setVolume(volumeInport_.getData());
    LGL_ERROR;
}

void IlluminationLineRaycaster::initBufferObjects(){
    //Generate Vertex Buffer Objects
    if(!vboID_[0])
        glGenBuffers(2, &vboID_[0]);

    int size_x = screenSize_.x;
    int size_y = screenSize_.y;
    int max_size = std::max(size_x, size_y);

    //Fill vbo with data for drawing vertical lines
    glBindBuffer(GL_ARRAY_BUFFER, vboID_[0]);
    glBufferData(GL_ARRAY_BUFFER, 2*size_x*sizeof(Vertex2D), NULL, GL_STATIC_DRAW);

    Vertex2D* verticesX = new Vertex2D[2*size_x];

    float x = 0.f;
    for (int i=0; i<2*size_x; i += 2) {
        verticesX[i].x = verticesX[i+1].x = ((x/(float)size_x)*2.f)-1.f;
        verticesX[i].y = -1.f;
        verticesX[i+1].y = 1.f;
        verticesX[i].s = verticesX[i+1].s = ((x/(float)size_x));
        verticesX[i].t = 0.f;
        verticesX[i+1].t = 1.f;
        x++;
    }

    glBufferSubData(GL_ARRAY_BUFFER, 0, 2*size_x*sizeof(Vertex2D), &verticesX[0].x);

    //Fill vbo with data for drawing horizontal lines
    glBindBuffer(GL_ARRAY_BUFFER, vboID_[1]);
    glBufferData(GL_ARRAY_BUFFER, 2*size_y*sizeof(Vertex2D), NULL, GL_STATIC_DRAW);

    Vertex2D* verticesY = new Vertex2D[2*size_y];

    float y = 0.f;
    for (int i=0; i<2*size_y; i += 2) {
        verticesY[i].x = -1.f;
        verticesY[i+1].x = 1.f;
        verticesY[i].y = verticesY[i+1].y = ((y/(float)size_y)*2.f)-1.f;
        verticesY[i].s = verticesY[i+1].s = ((y/(float)size_y));
        verticesY[i].t = 0.f;
        verticesY[i+1].t = 1.f;
        y++;
    }

    glBufferSubData(GL_ARRAY_BUFFER, 0, 2*size_y*sizeof(Vertex2D), &verticesY[0].x);

    //Generate Index Buffer Objects
    if(!iboID_)
        glGenBuffers(1, &iboID_);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, iboID_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 2*max_size*sizeof(unsigned short), NULL, GL_STATIC_DRAW);

    //Fill IBOs
    unsigned short* indices = new unsigned short[2*max_size];

    for(unsigned short i=0; i<2*(unsigned short)max_size; i++)
        indices[(int)i]=i;

    glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, 2*max_size*sizeof(unsigned short), indices);

    delete[] verticesX;
    delete[] verticesY;
    delete[] indices;
}

//Determine propogation axis
void IlluminationLineRaycaster::updatePropagation() {
    //Project focus point and light position onto screen.
    tgt::mat4 pointToScreen = camera_.get().getProjectionMatrix(outport1_.getSize()) * camera_.get().getViewMatrix();
    projLightPos_ = pointToScreen * tgt::vec4(lightPosition_.get().xyz(), 1.0);
    tgt::vec4 projFocusPos = pointToScreen * tgt::vec4(camera_.get().getFocus(), 1.0);

    //Division by the homogeneous coordinate
    projLightPos_ = projLightPos_/std::abs((float)projLightPos_._w);
    projFocusPos = projFocusPos/std::abs((float)projFocusPos._w);

    //Determine closest major axis (positive or negative x,y) in camera eye space.
    //Also, calculate length in light direction to cover the spacing between lines.
    projLightDirCam_ = projFocusPos - projLightPos_;

    //Check if close to middle of screen, i.e. if camera and light vector are almost parallell.
    if(minLightDirLength_.get() > 0.f || noLinesOverride_){
        float lightVecLength = tgt::length(projLightDirCam_);
        bool tmp_override = (lightVecLength < minLightDirLength_.get());

        if(tmp_override != noLinesOverride_){
            noLinesOverride_ = tmp_override;
            updateHeader();
        }

        if(noLinesOverride_)
            return;
    }

    projLightDirCam_ = tgt::normalize(projLightDirCam_);

    propigationAxis_ = tgt::ivec2(0, 0);
    verticalLines_ = false;
    negativeDir_ = false;
    if(abs(projLightDirCam_.x) > abs(projLightDirCam_.y)){
        verticalLines_ = true;
        projLightDirCam_ *= ((float)lineWidth_.get())/abs(projLightDirCam_.x);
        currentVBO_ = 0;
        if(projLightDirCam_.x < 0){
            negativeDir_ = true;
            propigationAxis_.x = -1;
        }
        else
            propigationAxis_.x = 1;
    }
    else{
        projLightDirCam_ *= ((float)lineWidth_.get())/abs(projLightDirCam_.y);
        currentVBO_ = 1;
        if(projLightDirCam_.y < 0){
            negativeDir_ = true;
            propigationAxis_.y = -1;
        }
        else
            propigationAxis_.y = 1;
    }

    //Use regular light position when inside volume
    if(zoneRenderingOn_.get()){
        lightCamera_.setPosition(lightPosition_.get().xyz());
    }
    else{
        //Calculate distance to camera
        float dist_to_cam = tgt::distance(camera_.get().getFocus(), camera_.get().getPosition());

        //Define light vector
        tgt::vec3 lightVec = tgt::normalize(lightPosition_.get().xyz() - camera_.get().getFocus());

        //Set light camera
        lightCamera_.setPosition(camera_.get().getFocus() + (diffCamLightDistance_.get()+dist_to_cam)*lightVec);
        lightCamera_.setFocus(camera_.get().getFocus());

        //Define proxy geometry size seen from light source, to define new near plane
        if(adaptiveNearPlane_.get() && !geometry_.empty() && lightCamera_.getProjectionMode() == tgt::Camera::PERSPECTIVE){
            tgt::mat4 objToLightMat = lightCamera_.getProjectionMatrix(outport1_.getSize()) * lightCamera_.getViewMatrix() * volumeInport_.getData()->getPhysicalToWorldMatrix();

            //Get proxy size which will be used to define the near plane
            tgt::vec4 resultProxyScreenSize = calculateProxyGeometryScreenSize(&geometry_, &objToLightMat);

            if(resultProxyScreenSize.x > 0.0 && resultProxyScreenSize.y > 0.0 && resultProxyScreenSize.z < 1.0 && resultProxyScreenSize.y < 1.0){
                //Calculate max half height and width scaling to cover the proxy geometry
                float Hscale = std::max(abs(resultProxyScreenSize.y)+0.5f, abs(resultProxyScreenSize.w)-0.5f);

                //Calculate height and width of near plane
                float Hnear = Hscale*lightCamera_.getFrustTop();

                //Calculate and set new field of view
                float fov = tgt::rad2deg(2.f*atanf(Hnear/lightCamera_.getNearDist()));
                lightCamera_.setFovy(fov);
            }
        }
    }
}

void IlluminationLineRaycaster::clearLightMaps(){
    bool outportActive = outport1_.isActive();
    bool propErtPortActive = propagationErtPort_.isActive();

    if(outportActive)
        outport1_.deactivateTarget();

    if(propErtPortActive)
        propagationErtPort_.deactivateTarget();

    glClearColor(lightDiffuse_.get().x, lightDiffuse_.get().y, lightDiffuse_.get().z, 1.f);

    lightportOne_.activateTarget();
    lightportOne_.clearTarget();
    lightportOne_.deactivateTarget();

    lightportTwo_.activateTarget();
    lightportTwo_.clearTarget();
    lightportTwo_.deactivateTarget();

    glClearColor(0.f, 0.f, 0.f, 0.f);

    lightportOnePos_.activateTarget();
    lightportOnePos_.clearTarget();
    lightportOnePos_.deactivateTarget();

    lightportTwoPos_.activateTarget();
    lightportTwoPos_.clearTarget();
    lightportTwoPos_.deactivateTarget();

    if(outportActive)
        outport1_.activateTarget();

    if(propErtPortActive)
        propagationErtPort_.activateTarget();
}

void IlluminationLineRaycaster::drawQuadIntoList(){
    if(quadList_)
        glDeleteLists(quadList_, 1);

    quadList_ = glGenLists(1);

    glNewList(quadList_, GL_COMPILE);
    renderQuad();
    glEndList();
}

#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
bool IlluminationLineRaycaster::postProcessLine(int *processed_lines){
    glFinish();

    if(processed_lines)
        return true;
    else
        return false;
}
#else
void IlluminationLineRaycaster::postProcessLine(){
    glFinish();
}
#endif

#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
bool IlluminationLineRaycaster::postProcessLineGL4(int *processed_lines){
#else
void IlluminationLineRaycaster::postProcessLineGL4(){
#endif
#ifdef GL_EXT_shader_image_load_store
        glMemoryBarrierEXT(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT_EXT);
#endif

#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
    if(outputLineToDebugRT_.get() > 0 && outputLineToDebugRT_.get() <= (*processed_lines)){
        //Stop rendering.
        return false;
    }

    if(pingPong_.get()){
#endif
        //Ping-pong buffering
        int tmp = readWriteValues_.x;
        readWriteValues_.x = readWriteValues_.y;
        readWriteValues_.y = tmp;
        raycastPrg_->setUniform("lightMapReadXWriteY", readWriteValues_);
#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
    }
    return true;
#endif
}

#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
bool IlluminationLineRaycaster::drawSynchronizedLine(tgt::vec2 *startVertex, tgt::vec2 *endVertex, tgt::vec2 *startTexCoord, tgt::vec2 *endTexCoord, int *numLines){
    (*numLines)++;
#else
void IlluminationLineRaycaster::drawSynchronizedLine(tgt::vec2 *startVertex, tgt::vec2 *endVertex, tgt::vec2 *startTexCoord, tgt::vec2 *endTexCoord){
#endif

#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
    if(drawTrianglesAsLines_.get()){
        // specify line between points
        tgt::vec2 lv = (*endVertex) - (*startVertex);
        tgt::vec2 lt = (*endTexCoord) - (*startTexCoord);

        // specify perpendicular vector
        tgt::vec2 perpv = tgt::normalize(tgt::vec2(-lv.y, lv.x));
        tgt::vec2 perpt = tgt::normalize(tgt::vec2(-lt.y, lt.x));
        perpv *= triangleThickness_;
        perpt *= triangleThickness_;

        // specify new vertices
        tgt::vec2 Av = (*startVertex) + perpv;
        tgt::vec2 Bv = (*startVertex) - perpv;
        tgt::vec2 Cv = (*endVertex) + perpv;
        tgt::vec2 Dv = (*endVertex) - perpv;

        // specify new texcoords
        tgt::vec2 At = (*startTexCoord) + perpt;
        tgt::vec2 Bt = (*startTexCoord) - perpt;
        tgt::vec2 Ct = (*endTexCoord) + perpt;
        tgt::vec2 Dt = (*endTexCoord) - perpt;

        // draw triangle strip
        glBegin(GL_TRIANGLE_STRIP);
            glTexCoord2f(At.x, At.y);
            glVertex2f(Av.x, Av.y);
            glTexCoord2f(Bt.x, Bt.y);
            glVertex2f(Bv.x, Bv.y);
            glTexCoord2f(Ct.x, Ct.y);
            glVertex2f(Cv.x, Cv.y);
            glTexCoord2f(Dt.x, Dt.y);
            glVertex2f(Dv.x, Dv.y);
        glEnd();
    }
    else{
#endif
        glBegin(GL_LINES);
        glTexCoord2f(startTexCoord->x, startTexCoord->y);
        glVertex2f(startVertex->x, startVertex->y);
        glVertex2f(endVertex->x, endVertex->y);
        glTexCoord2f(endTexCoord->x, endTexCoord->y);
        glEnd();
#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
    }
    return (this->*plfFunc_)(numLines);
#else
    (this->*plfFunc_)();
#endif
}

//Lines will be drawn in four different direction: x+, x-, y+, y-, according to the projected light source direction.
//Coords have following alignment: x = x_min, y = y_min, z = x_max, w = y_max
bool IlluminationLineRaycaster::drawLines(tgt::ivec4 coords, tgt::ivec2 light_coords, int lineWidth){
    tgt::ivec2 lastLine = screenSize_ - 1;
    bool abort = false;
    bool negDir = !negativeDir_;
    int increment = lineWidth;
    int start_offset = std::max(1, (int)floor(((float)increment)/2.f)+1);
    float offset = (lineWidth % 2) ? 0.f : 0.5f;

#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
    triangleThickness_ = (1.f)/abs((float)tgt::dot(screenSize_, propigationAxis_));
#endif

    //Perform one or two pass dependent on light source position in relation to camera
    int num_of_passes = 1;
    tgt::ivec4 divided_coords[2] = { coords, coords };
    if(light_coords.x > 0.f){
        negDir = false;
        num_of_passes = 2;
        divided_coords[0].z = light_coords.x;
        divided_coords[0].w = light_coords.y;
        divided_coords[1].x = light_coords.x;
        divided_coords[1].y = light_coords.y;
    }

    tgt::vec2 start_vert, end_vert;
    tgt::vec2 start_tex, end_tex;

    start_tex.y = 0.f;
    end_tex.y = 1.f;

    //Draw in one or two passes
    for(int p=0; p<num_of_passes; p++){
        //Clear light maps
        clearLightMaps();

        negDir = !negDir;

        #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
        int proccessed_lines = 0;
        #endif

        if (verticalLines_) {
            start_vert.y = -1.f;
            end_vert.y = 1.f;

            if(negDir) {
                int firstLine = divided_coords[p].z - 1;
                raycastPrg_->setUniform("firstLine_", firstLine);
                raycastPrg_->setUniform("lineWidth_", 1);
                start_tex.x = end_tex.x = (((float)firstLine)-offset)/(float)lastLine.x;
                start_vert.x = end_vert.x = (start_tex.x*2.f)-1.f;
                #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                if(!drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex, &proccessed_lines)){
                        abort = true;
                #else
                    drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex);
                #endif
                #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                }
                else{
                #endif

                #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                triangleThickness_ *= (float)lineWidth;
                if(!drawTrianglesAsLines_.get())
                #endif
                glLineWidth(((float)lineWidth));
                raycastPrg_->setUniform("lineWidth_", lineWidth_.get());

                for (int x=firstLine-start_offset; x>=divided_coords[p].x; x-=increment) {

                    start_tex.x = end_tex.x = (((float)x)-offset)/(float)lastLine.x;
                    start_vert.x = end_vert.x = (start_tex.x*2.f)-1.f;

                #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                    if(!drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex, &proccessed_lines)){
                        abort = true;
                        break;
                    }
                #else
                    drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex);
                #endif

                }

                #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                }
                #endif
            }
            else{
                int firstLine = divided_coords[p].x;
                raycastPrg_->setUniform("firstLine_", firstLine);
                raycastPrg_->setUniform("lineWidth_", 1);
                start_tex.x = end_tex.x = (((float)firstLine)+offset)/(float)lastLine.x;
                start_vert.x = end_vert.x = (start_tex.x*2.f)-1.f;
                #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                if(!drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex, &proccessed_lines)){
                        abort = true;
                #else
                    drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex);
                #endif
                #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                }
                else{
                #endif

                #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                triangleThickness_ *= (float)lineWidth;
                if(!drawTrianglesAsLines_.get())
                #endif
                glLineWidth(((float)lineWidth));
                raycastPrg_->setUniform("lineWidth_", lineWidth_.get());

                for (int x=firstLine+start_offset; x<=divided_coords[p].z; x+=increment) {

                    start_tex.x = end_tex.x = (((float)x)+offset)/(float)lastLine.x;
                    start_vert.x = end_vert.x = (start_tex.x*2.f)-1.f;

                #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                    if(!drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex, &proccessed_lines)){
                        abort = true;
                        break;
                    }
                #else
                    drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex);
                #endif
                }

                #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                }
                #endif
            }
        } else {
            start_vert.x = -1.f;
            end_vert.x = 1.f;

            if(negDir) {
                int firstLine = divided_coords[p].w-1;
                raycastPrg_->setUniform("firstLine_", firstLine);
                raycastPrg_->setUniform("lineWidth_", 1);
                start_tex.x = end_tex.x = (((float)firstLine)-offset)/(float)lastLine.y;
                start_vert.y = end_vert.y = (start_tex.x*2.f)-1.f;
                #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                if(!drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex, &proccessed_lines)){
                        abort = true;
                #else
                    drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex);
                #endif
                #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                }
                else{
                #endif

                #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                triangleThickness_ *= (float)lineWidth;
                if(!drawTrianglesAsLines_.get())
                #endif
                glLineWidth(((float)lineWidth));
                raycastPrg_->setUniform("lineWidth_", lineWidth_.get());

                for (int y=firstLine-start_offset; y>=divided_coords[p].y; y-=increment) {

                    start_tex.x = end_tex.x = (((float)y)-offset)/(float)lastLine.y;
                    start_vert.y = end_vert.y = (start_tex.x*2.f)-1.f;

                #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                    if(!drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex, &proccessed_lines)){
                        abort = true;
                        break;
                    }
                #else
                    drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex);
                #endif
                }

                #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                }
                #endif
            }
            else{
                int firstLine = divided_coords[p].y;
                raycastPrg_->setUniform("firstLine_", firstLine);
                raycastPrg_->setUniform("lineWidth_", 1);
                start_tex.x = end_tex.x = (((float)firstLine)+offset)/(float)lastLine.y;
                start_vert.y = end_vert.y = (start_tex.x*2.f)-1.f;
                #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                if(!drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex, &proccessed_lines)){
                        abort = true;
                #else
                    drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex);
                #endif
                #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                }
                else{
                #endif

                #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                triangleThickness_ *= (float)lineWidth;
                if(!drawTrianglesAsLines_.get())
                #endif
                glLineWidth(((float)lineWidth));
                raycastPrg_->setUniform("lineWidth_", lineWidth_.get());

                for (int y=firstLine+start_offset; y<=divided_coords[p].w; y+=increment) {

                    start_tex.x = end_tex.x = (((float)y)+offset)/(float)lastLine.y;
                    start_vert.y = end_vert.y = (start_tex.x*2.f)-1.f;

                #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                    if(!drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex, &proccessed_lines)){
                        abort = true;
                        break;
                    }
                #else
                    drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex);
                #endif
                }

                #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                }
                #endif
            }
        }
    }

    #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
        if(!drawTrianglesAsLines_.get())
    #endif
    glLineWidth(1.f);

    return abort;
}

bool IlluminationLineRaycaster::drawLinesVBO(tgt::ivec4 coords, tgt::ivec2 light_coords, int lineWidth){
    //Define variables
    int increment = lineWidth;
    int amount_of_lines = abs(tgt::dot(propigationAxis_, screenSize_));
    bool abort = false;
    bool negDir = !negativeDir_;

    //Perform one or two pass dependent on light source position in relation to camera
    int num_of_passes = 1;
    tgt::ivec4 divided_coords[2] = { coords, coords };
    if(light_coords.x > 0.f){
        negDir = false;
        num_of_passes = 2;
        divided_coords[0].z = light_coords.x;
        divided_coords[0].w = light_coords.y;
        divided_coords[1].x = light_coords.x;
        divided_coords[1].y = light_coords.y;
    }

    //Enable State
    glEnableClientState(GL_VERTEX_ARRAY);

    //Bind VBO
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, vboID_[currentVBO_]);

    //Set Pointers, Interleaved Array
    glVertexPointer(2, GL_FLOAT, 16, BUFFER_OFFSET(0));
    glTexCoordPointer(2, GL_FLOAT, 16, BUFFER_OFFSET(8));

    //Bind IBO
    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, iboID_);

    //Draw in one or two passes
    for(int p=0; p<num_of_passes; p++){
        //Clear light maps
        clearLightMaps();

        negDir = !negDir;

        //Calculate order
        int start_offset = abs(tgt::dot(propigationAxis_, divided_coords[p].xy()))*2;
        int end_offset = abs(tgt::dot(propigationAxis_, screenSize_-divided_coords[p].zw()))*2;

        //Draw
        #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
        int proccessed_lines = 0;
        #endif
        if(negDir){
            for(int i=(2*amount_of_lines)-1-end_offset; i>start_offset; i-=2*increment){
                glDrawRangeElements(GL_LINES, i, i+1, 2*increment, GL_UNSIGNED_SHORT, BUFFER_OFFSET(i*2));

        #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                if(!(this->*plfFunc_)(&proccessed_lines)){ abort = true; break; }
        #else
                (this->*plfFunc_)();
        #endif
            }
        }
        else{
            for(int i=start_offset; i<(2*amount_of_lines)-end_offset; i+=2*increment){
                glDrawRangeElements(GL_LINES, i, i+1, 2*increment, GL_UNSIGNED_SHORT, BUFFER_OFFSET(i*2));

        #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                if(!(this->*plfFunc_)(&proccessed_lines)){ abort = true; break; }
        #else
                (this->*plfFunc_)();
        #endif
            }
        }
    }

    //Deativate VBO
    glDisableClientState(GL_VERTEX_ARRAY);

    //Bind with 0 -> switch back to normal pointer operation
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, 0);

    return abort;
}

//Draw lines in four passes based on the propagation zones
bool IlluminationLineRaycaster::drawZones(tgt::ivec4 proxy_coords, tgt::ivec2 screen_space_light_pos, int lineWidth){
    #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
    triangleThickness_ = ((float)lineWidth)/abs((float)tgt::dot(screenSize_, propigationAxis_));
        if(!drawTrianglesAsLines_.get())
    #endif
    glLineWidth(((float)lineWidth));
    bool abort = false;
    int increment = lineWidth;
    tgt::vec2 lastLine = tgt::vec2(screenSize_ - 1);

    //Define the four vectors in screen space with origin at the light source
    tgt::ivec2 vec12 = tgt::ivec2(-1, -1);
    tgt::ivec2 vec32 = tgt::ivec2(1, -1);
    tgt::ivec2 vec34 = tgt::ivec2(1, 1);
    tgt::ivec2 vec14 = tgt::ivec2(-1, 1);

    //Check which zones to render
    bool renderZone1 = (screen_space_light_pos.x >= proxy_coords.x) ? (zoneOnOff_.getValue() == 1 || zoneOnOff_.getValue() == 0) : false;
    bool renderZone2 = (screen_space_light_pos.y >= proxy_coords.y) ? (zoneOnOff_.getValue() == 2 || zoneOnOff_.getValue() == 0) : false;
    bool renderZone3 = (screen_space_light_pos.x <= proxy_coords.z) ? (zoneOnOff_.getValue() == 3 || zoneOnOff_.getValue() == 0) : false;
    bool renderZone4 = (screen_space_light_pos.y <= proxy_coords.w) ? (zoneOnOff_.getValue() == 4 || zoneOnOff_.getValue() == 0) : false;

    //Render Algorithm
    // - Start from start pos
    // - While position doesn't lie beyond the proxy geometry, i.e. further away from light source then proxy geometry
    //        - Move a step using the two vectors forming the zone with origin in the light source
    //        - Clamp values to the proxy geometry for the coordinate that is not static for each line.
    //        - If middle of line inside proxy geometry, then render (to avoid line-line intersection check)
    //        - Draw synchronized line

    tgt::ivec2 start, end, start_clamp, end_clamp;
    tgt::vec2 start_vert, end_vert, start_tex, end_tex;
    #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
    int    proccessed_lines;
    #endif

    //Set factor for edge blending
    raycastPrg_->setUniform("edgeFactor_", zoneEdgeBlendingFactor_.get());
    tgt::ivec2 factor_x = tgt::ivec2(zoneEdgeBlendingFactor_.get(), 1);
    tgt::ivec2 factor_y = tgt::ivec2(1, zoneEdgeBlendingFactor_.get());

    tgt::vec3 lightVec = camera_.get().getFocus() - lightPosition_.get().xyz();

    //Calculate new field of view related to overlap.
    float fov = lightCamera_.getFovy();
    float new_fov = tgt::rad2deg(acosf(tgt::dot(
        tgt::normalize(tgt::vec2(-(float)zoneEdgeBlendingFactor_.get(), 1.f)),
        tgt::normalize(tgt::vec2((float)zoneEdgeBlendingFactor_.get(), 1.f)))))/2.f;
    lightCamera_.setFovy(new_fov);

    //Draw Zone 1: Vertical lines, processed from right to left, angle 135-225.
    if(renderZone1){

        //Clear light maps
        clearLightMaps();

        //Set light camera
        lightCamera_.setFocus(lightPosition_.get().xyz() + (-camera_.get().getStrafe() + lightVec));
        raycastPrg_->setUniform("objToLightMat_", lightCamera_.getViewMatrix() * volumeInport_.getData()->getPhysicalToWorldMatrix());

        //Set Vectors that form area
        raycastPrg_->setUniform("posDirZoneVec_", vec14);
        raycastPrg_->setUniform("negDirZoneVec_", vec12);

        //Set Correct Propagation Axis
        raycastPrg_->setUniform("propagationDirection_", tgt::ivec2(-1, 0));

        //Set current zone
        raycastPrg_->setUniform("currentRenderedZone_", 1);

        #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
        proccessed_lines = 0;
        #endif
        start = screen_space_light_pos;// + (vec14*increment); //start = top
        end = screen_space_light_pos;// + (vec12*increment);   //end = bottom

        while(start.x > proxy_coords.x){
            if(start.x < proxy_coords.z && (start.y > proxy_coords.y || end.y < proxy_coords.w)){
                start_clamp.y = std::min(start.y, proxy_coords.w);
                end_clamp.y = std::max(end.y, proxy_coords.y);

                start_tex.x = ((float)start.x)/lastLine.x;
                start_tex.y = ((float)start_clamp.y)/lastLine.y;
                end_tex.x = ((float)end.x)/lastLine.x;
                end_tex.y = ((float)end_clamp.y)/lastLine.y;

                start_vert = (start_tex*2.f)-1.f;
                end_vert = (end_tex*2.f)-1.f;

#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                if(!drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex, &proccessed_lines)){
                    abort = true;
                    break;
                }
#else
                drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex);
#endif
            }

            start += vec14*factor_y*increment;
            end += vec12*factor_y*increment;
        }
    }

    //Draw Zone 2: Horizontal lines, processed from top to bottom, angle 225-315.
    if(renderZone2){

        //Clear light maps
        clearLightMaps();

        //Set light camera
        lightCamera_.setFocus(lightPosition_.get().xyz() + (-camera_.get().getUpVector() + lightVec));
        raycastPrg_->setUniform("objToLightMat_", lightCamera_.getViewMatrix() * volumeInport_.getData()->getPhysicalToWorldMatrix());

        //Set Vectors that form area
        raycastPrg_->setUniform("posDirZoneVec_", vec32);
        raycastPrg_->setUniform("negDirZoneVec_", vec12);

        //Set Correct Propagation Axis
        raycastPrg_->setUniform("propagationDirection_", tgt::ivec2(0, -1));

        //Set current zone
        raycastPrg_->setUniform("currentRenderedZone_", 2);

        #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
        proccessed_lines = 0;
        #endif
        start = screen_space_light_pos;// + (vec12*increment); //start = left
        end = screen_space_light_pos;// + (vec32*increment);   //end = right

        while(start.y > proxy_coords.y){
            if(start.y < proxy_coords.w && (start.x < proxy_coords.z || end.x > proxy_coords.x)){
                start_clamp.x = std::max(start.x, proxy_coords.x);
                end_clamp.x = std::min(end.x, proxy_coords.z);

                start_tex.x = ((float)start_clamp.x)/lastLine.x;
                start_tex.y = ((float)start.y)/lastLine.y;
                end_tex.x = ((float)end_clamp.x)/lastLine.x;
                end_tex.y = ((float)end.y)/lastLine.y;
                start_vert = (start_tex*2.f)-1.f;
                end_vert = (end_tex*2.f)-1.f;

#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                if(!drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex, &proccessed_lines)){
                    abort = true;
                    break;
                }
#else
                drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex);
#endif
            }

            start += vec12*factor_x*increment;
            end += vec32*factor_x*increment;
        }
    }

    //Draw Zone 3: Vertical lines, processed from left to right, angle 315-45.
    if(renderZone3){

        //Clear light maps
        clearLightMaps();

        //Set light camera
        lightCamera_.setFocus(lightPosition_.get().xyz() + (camera_.get().getStrafe() + lightVec));
        raycastPrg_->setUniform("objToLightMat_", lightCamera_.getViewMatrix() * volumeInport_.getData()->getPhysicalToWorldMatrix());

        //Set Vectors that form area
        raycastPrg_->setUniform("posDirZoneVec_", vec34);
        raycastPrg_->setUniform("negDirZoneVec_", vec32);

        //Set Correct Propagation Axis
        raycastPrg_->setUniform("propagationDirection_", tgt::ivec2(1, 0));

        //Set current zone
        raycastPrg_->setUniform("currentRenderedZone_", 3);

        #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
        proccessed_lines = 0;
        #endif
        start = screen_space_light_pos;// + (vec34*increment); //start = top
        end = screen_space_light_pos;// + (vec32*increment);   //end = bottom

        while(start.x < proxy_coords.z){
            if(start.x > proxy_coords.x && (start.y > proxy_coords.x || end.y < proxy_coords.w)){
                start_clamp.y = std::min(start.y, proxy_coords.w);
                end_clamp.y = std::max(end.y, proxy_coords.y);

                start_tex.x = ((float)start.x)/lastLine.x;
                start_tex.y = ((float)start_clamp.y)/lastLine.y;
                end_tex.x = ((float)end.x)/lastLine.x;
                end_tex.y = ((float)end_clamp.y)/lastLine.y;

                start_vert = (start_tex*2.f)-1.f;
                end_vert = (end_tex*2.f)-1.f;

#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                if(!drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex, &proccessed_lines)){
                    abort = true;
                    break;
                }
#else
                drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex);
#endif
            }

            start += vec34*factor_y*increment;
            end += vec32*factor_y*increment;
        }
    }

    //Draw Zone 4: Horizontal lines, processed from bottom to top, angle 45-135.
    if(renderZone4){

        //Clear light maps
        clearLightMaps();

        //Set light camera
        lightCamera_.setFocus(lightPosition_.get().xyz() + (camera_.get().getUpVector() + lightVec));
        raycastPrg_->setUniform("objToLightMat_", lightCamera_.getViewMatrix() * volumeInport_.getData()->getPhysicalToWorldMatrix());

        //Set Vectors that form area
        raycastPrg_->setUniform("posDirZoneVec_", vec34);
        raycastPrg_->setUniform("negDirZoneVec_", vec14);

        //Set Correct Propagation Axis
        raycastPrg_->setUniform("propagationDirection_", tgt::ivec2(0, 1));

        //Set current zone
        raycastPrg_->setUniform("currentRenderedZone_", 4);

        #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
        proccessed_lines = 0;
        #endif
        start = screen_space_light_pos;// + (vec14*increment); //start = left
        end = screen_space_light_pos;// + (vec34*increment);   //end = right

        while(start.y < proxy_coords.w){
            if(start.y > proxy_coords.y && (start.x < proxy_coords.w || end.x > proxy_coords.x)){
                start_clamp.x = std::max(start.x, proxy_coords.x);
                end_clamp.x = std::min(end.x, proxy_coords.w);

                start_tex.x = ((float)start_clamp.x)/lastLine.x;
                start_tex.y = ((float)start.y)/lastLine.y;
                end_tex.x = ((float)end_clamp.x)/lastLine.x;
                end_tex.y = ((float)end.y)/lastLine.y;

                start_vert = (start_tex*2.f)-1.f;
                end_vert = (end_tex*2.f)-1.f;

#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
                if(!drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex, &proccessed_lines)){
                    abort = true;
                    break;
                }
#else
                drawSynchronizedLine(&start_vert, &end_vert, &start_tex, &end_tex);
#endif
            }

            start += vec14*factor_x*increment;
            end += vec34*factor_x*increment;
        }
    }

    lightCamera_.setFovy(fov);
    lightCamera_.setFocus(camera_.get().getFocus());

    #ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
        if(!drawTrianglesAsLines_.get())
    #endif
    glLineWidth(1.f);
    return abort;
}

//Calculate a proxy geometry boundary
tgt::vec4 IlluminationLineRaycaster::calculateProxyGeometryScreenSize(MeshListGeometry *geom, const tgt::mat4 *transform, const tgt::vec2 pos, bool *isInside){
    tgt::vec4 limits = tgt::vec4(1.f, 1.f, 0.f, 0.f);

    if(geom){
        tgt::vec4 p;
        for (MeshListGeometry::const_iterator itMesh = geom->begin(); itMesh != geom->end(); ++itMesh){
            for (MeshGeometry::const_iterator itFace = (*itMesh).begin(); itFace != (*itMesh).end(); ++itFace){
                for (FaceGeometry::const_iterator itVertex = (*itFace).begin(); itVertex != (*itFace).end(); ++itVertex){
                    p = tgt::vec4((*itVertex).getCoords(), 1.f);

                    if(transform){
                        p = (*transform) * p;
                        p /= p.w;
                    }

                    p = p*0.5f + 0.5f;
                    limits.x = (p.x < limits.x) ? p.x : limits.x;
                    limits.y = (p.y < limits.y) ? p.y : limits.y;
                    limits.z = (p.x > limits.z) ? p.x : limits.z;
                    limits.w = (p.y > limits.w) ? p.y : limits.w;
                }
            }
        }
    }

    if(isInside){
        (*isInside) = (pos.x > limits.x && pos.x < limits.z &&
                       pos.y > limits.y && pos.y < limits.w)? true : false;
    }

    return limits;
}

void IlluminationLineRaycaster::calculateERTmap(const GLint entryUnit, const GLint exitUnit){
#ifdef ERT_DOWNSAMPLE_SCHEME
    //propagationErtPort_.resize(ertTextureSizeWidth_.get(), ertTextureSizeHeight_.get());
#endif

    // activate and clear output render target
    propagationErtPort_.activateTarget();
    propagationErtPort_.clearTarget();

    TextureUnit texUnit;
    propagationErtPort_.bindColorTexture(texUnit);

#ifdef ERT_DOWNSAMPLE_SCHEME
    //TextureUnit ertPropTexUnit;
    //ertPort_.bindColorTexture(ertPropTexUnit);
#endif

    propErtPrg_->activate();
    tgt::Camera cam = camera_.get();
    setGlobalShaderParameters(propErtPrg_, &cam);

    // pass texture units to the shader
    propErtPrg_->setUniform("entryPoints_", entryUnit);
    entryPort_.setTextureParameters(propErtPrg_, "entryParameters_");

    propErtPrg_->setUniform("exitPoints_", exitUnit);
    exitPort_.setTextureParameters(propErtPrg_, "exitParameters_");

#ifdef ERT_DOWNSAMPLE_SCHEME
    //propErtPrg_->setUniform("ertExitPoints_", ertPropTexUnit.getUnitNumber());
    //ertPort_.setTextureParameters(propErtPrg_, "ertExitParameters_");
#endif

    propErtPrg_->setUniform("projectedLightDir_", projLightDirCam_.xy());

    propErtPrg_->setUniform("ertTexture_", texUnit.getUnitNumber());
    propagationErtPort_.setTextureParameters(propErtPrg_, "ertTextureParams_");

    //Don't use image store extension functionality
    plfFunc_ = &IlluminationLineRaycaster::postProcessLine;

    //Reverse processing, and override some variables
    int tmp_neg = negativeDir_;
    tgt::ivec2 tmp_size = screenSize_;
    bool tmp_glsl = glsl4On_;

    negativeDir_ = !negativeDir_;
#ifdef ERT_DOWNSAMPLE_SCHEME
    screenSize_ = tgt::ivec2(ertTextureSizeWidth_.get(), ertTextureSizeHeight_.get());
#else
    screenSize_ = propagationErtPort_.getSize();
#endif
    glsl4On_ = false;

#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
    if(vboOn_.get())
        drawLinesVBO(tgt::ivec4(tgt::ivec2(0), screenSize_ - 1), tgt::ivec2(-1, -1));
    else
#endif
        drawLines(tgt::ivec4(tgt::ivec2(0), screenSize_ - 1), tgt::ivec2(-1, -1));

    // clean up
    negativeDir_ = tmp_neg;
    screenSize_ = tmp_size;
    glsl4On_ = tmp_glsl;

    propErtPrg_->deactivate();
    propagationErtPort_.deactivateTarget();
}

//Calculate a texture which defines zones with different propagation direction.
//4 zones -> x-(1), y-(2), x+(3), y+(4)
//Visualized for now as RGBY, suffices with 2D one component later
void IlluminationLineRaycaster::calculatePropagationZones(const tgt::vec2 light_pos, const GLint entryUnit, const GLint exitUnit){
    // activate and clear output render target
    propagationZonesPort_.activateTarget();
    propagationZonesPort_.clearTarget();

    propZonePrg_->activate();
    tgt::Camera cam = camera_.get();
    setGlobalShaderParameters(propZonePrg_, &cam);

    // pass texture units to the shader
    propZonePrg_->setUniform("entryPoints_", entryUnit);
    entryPort_.setTextureParameters(propZonePrg_, "entryParameters_");

    propZonePrg_->setUniform("exitPoints_", exitUnit);
    exitPort_.setTextureParameters(propZonePrg_, "exitParameters_");

    propZonePrg_->setUniform("projectedLight_", light_pos);

    glCallList(quadList_);

    propZonePrg_->deactivate();
    propagationZonesPort_.deactivateTarget();
}

#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
void IlluminationLineRaycaster::copyRenderPort(RenderPort *src, RenderPort *dst){
    // activate and clear output render target
    dst->activateTarget();
    dst->clearTarget();

    // activate shader and set common uniforms
    copyImagePrg_->activate();
    tgt::Camera cam = camera_.get();
    setGlobalShaderParameters(copyImagePrg_, &cam);

    TextureUnit texUnit;
    src->bindColorTexture(texUnit);

    // pass texture units to the shader
    copyImagePrg_->setUniform("colorTex_", texUnit.getUnitNumber());
    src->setTextureParameters(copyImagePrg_, "texParams_");

    // render screen aligned quad
    glDepthFunc(GL_ALWAYS);

    MatStack.pushMatrix();

    glCallList(quadList_);

    MatStack.popMatrix();

    glDepthFunc(GL_LESS);

    //clean up
    copyImagePrg_->deactivate();
    dst->deactivateTarget();
}

void IlluminationLineRaycaster::visualizeProxyGeometrySize(RenderPort *dst, RenderPort *src, GLint src_unit, tgt::vec4 size){
    copyRenderPort(src, dst);

    dst->activateTarget();

    blendPrg_->activate();

    tgt::Camera cam = camera_.get();
    setGlobalShaderParameters(blendPrg_, &cam);
    blendPrg_->setUniform("colorTex_", src_unit);
    src->setTextureParameters(blendPrg_, "textureParameters_");

    glPushAttrib(GL_TRANSFORM_BIT | GL_VIEWPORT_BIT);

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);

    MatStack.loadIdentity();

    glOrtho(0, dst->getColorTexture()->getWidth(), 0, dst->getColorTexture()->getHeight(), -1, 1);

    MatStack.pushMatrix();

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);

    MatStack.loadIdentity();

    glDisable(GL_DEPTH_TEST);

    MatStack.pushMatrix();

    glColor4f(0.f, 1.f, 0.f, 1.f);

    drawSquareWithLines(tgt::ivec2((int)size.x, (int)size.y),
        tgt::ivec2((int)size.z, (int)size.y),
        tgt::ivec2((int)size.z, (int)size.w),
        tgt::ivec2((int)size.x, (int)size.w));

    glFlush();

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);

    MatStack.popMatrix();

    glEnable(GL_DEPTH_TEST);

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);

    MatStack.popMatrix();

    glPopAttrib();

    blendPrg_->deactivate();

    dst->deactivateTarget();
}

//Draw square, default between -1 and 1
void IlluminationLineRaycaster::drawSquareWithLines(tgt::ivec2 v1, tgt::ivec2 v2, tgt::ivec2 v3, tgt::ivec2 v4){
    glLineWidth(4.f);
    glBegin(GL_LINE_STRIP);
    glVertex2i(v1.x, v1.y);
    glVertex2i(v2.x, v2.y);
    glVertex2i(v3.x, v3.y);
    glVertex2i(v4.x, v4.y);
    glVertex2i(v1.x, v1.y);
    glEnd();
    glLineWidth(1.f);
}
#endif

void IlluminationLineRaycaster::showAndHideProperties(){
#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
    linesOnlyAboveProxyGeometry_.setVisibleFlag(adjustToProxyGeometry_.get());
    visualizeProxySize_.setVisibleFlag(adjustToProxyGeometry_.get());
#endif
    shadowBrightness_.setVisibleFlag(colorContribution_.get());
    zoneOnOff_.setVisibleFlag(zoneRenderingOn_.get());
    zoneEdgeBlendingFactor_.setVisibleFlag(zoneRenderingOn_.get());
    lightCacheWidth_.setVisibleFlag(!autoLightCacheSize_.get());
    lightCacheHeight_.setVisibleFlag(!autoLightCacheSize_.get());
#ifdef ERT_DOWNSAMPLE_SCHEME
    ertTextureSizeWidth_.setVisibleFlag(ertIncrement_.get());
    ertTextureSizeWidth_.setVisibleFlag(ertMod_.get());
    ertTextureSizeHeight_.setVisibleFlag(ertMod_.get());
#endif
}

void IlluminationLineRaycaster::screenSizeChanged() {
    screenSize_ = outport1_.getSize();
#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
    float ratio = (float)outputLineToDebugRT_.get()/(float)outputLineToDebugRT_.getMaxValue();
    outputLineToDebugRT_.setMaxValue((verticalLines_) ? (screenSize_.x/lineWidth_.get())+2 : (screenSize_.y/lineWidth_.get())+2);
    outputLineToDebugRT_.set((int)((float)outputLineToDebugRT_.getMaxValue()*ratio));
    initBufferObjects();
#endif
}

void IlluminationLineRaycaster::process() {
    if(screenSize_ != outport1_.getSize()){
        screenSizeChanged();
    }

    glsl4On_ = (!renderQuad_.get() && !noLinesOverride_ && glslVersion_.getValue() >= tgt::GpuCapabilities::GlVersion::SHADER_VERSION_400
                    && GpuCaps.getShaderVersion() >= tgt::GpuCapabilities::GlVersion::SHADER_VERSION_400);

    tgt::Camera cam = camera_.get();

    // set lightport sizes
    if(!autoLightCacheSize_.get()){
        const tgt::ivec2 light_size = tgt::ivec2(lightCacheWidth_.get(), lightCacheHeight_.get());
        lightportOne_.resize(light_size);
        lightportTwo_.resize(light_size);
    }

    // setup light units
    TextureUnit lightUnitOne, lightUnitTwo;
    lightportOne_.bindColorTexture(lightUnitOne);
    lightportTwo_.bindColorTexture(lightUnitTwo);
    lightMapUnits_[0] = lightUnitOne.getUnitNumber();
    lightMapUnits_[1] = lightUnitTwo.getUnitNumber();

    TextureUnit lightUnitPosOne, lightUnitPosTwo;
    lightportOnePos_.bindColorTexture(lightUnitPosOne);
    lightportTwoPos_.bindColorTexture(lightUnitPosTwo);
    lightMapPosUnits_[0] = lightUnitPosOne.getUnitNumber();
    lightMapPosUnits_[1] = lightUnitPosTwo.getUnitNumber();

    // bind entry and exit params and pass texture units to the shader
    TextureUnit entryUnit, entryDepthUnit, exitUnit, exitDepthUnit;
    entryPort_.bindTextures(entryUnit, entryDepthUnit);
    exitPort_.bindTextures(exitUnit, exitDepthUnit);

    // computations of proxy size, propagation zones etc
    tgt::vec4 proxyMinMaxScreenCoords = tgt::vec4(0.f, 0.f, (float)screenSize_.x - 1, (float)screenSize_.y - 1);
    tgt::ivec2 dividePassCoords = tgt::ivec2(-1, -1);

    // convert projected light position from -1:1 -> 0:1
    tgt::vec2 projLightPos = (projLightPos_.xy()*0.5f) + 0.5f;
#if defined(GL_EXT_direct_state_access) && defined(GL_EXT_shader_image_load_store)
    tgt::vec2 lightMapSize = tgt::vec2((float)lightportOne_.getColorTexture()->getDimensions().x - 1, (float)lightportOne_.getColorTexture()->getDimensions().y - 1);

    if(glsl4On_){
        bool light_inside_volume = false;

        if(adjustToProxyGeometry_.get()){
            // retrieve input geometry
            const MeshListGeometry* meshListGeometry = dynamic_cast<const MeshListGeometry*>(geometryPort_.getData());
            const MeshGeometry* meshGeometry = dynamic_cast<const MeshGeometry*>(geometryPort_.getData());
            if (meshListGeometry) {
                geometry_ = *meshListGeometry;
            }
            else if (meshGeometry) {
                geometry_.clear();
                geometry_.addMesh(*meshGeometry);
            }
            else {
                geometry_.clear();
                geometry_.addMesh(MeshGeometry::createCube());
            }

            if(geometryPort_.hasChanged()){
                updatePropagation();
                projLightPos = (projLightPos_.xy()*0.5f) + 0.5f;
            }

            // calculate proxy geometry screen size
            tgt::mat4 objToCamMat = camera_.get().getProjectionMatrix(outport1_.getSize()) * camera_.get().getViewMatrix() * volumeInport_.getData()->getPhysicalToWorldMatrix();
            tgt::vec4 resultProxyScreenSize = calculateProxyGeometryScreenSize(&geometry_, &objToCamMat, projLightPos, &light_inside_volume);

            // clamp and scale coordinates
            proxyMinMaxScreenCoords.x = std::floor(resultProxyScreenSize.x*lightMapSize.x);
            proxyMinMaxScreenCoords.y = std::floor(resultProxyScreenSize.y*lightMapSize.y);
            proxyMinMaxScreenCoords.z = std::ceil(resultProxyScreenSize.z*lightMapSize.x);
            proxyMinMaxScreenCoords.w = std::ceil(resultProxyScreenSize.w*lightMapSize.y);

            // calculate light map size
            lightMapSize.x = std::min(proxyMinMaxScreenCoords.z - proxyMinMaxScreenCoords.x, lightMapSize.x);
            lightMapSize.y = std::min(proxyMinMaxScreenCoords.w - proxyMinMaxScreenCoords.y, lightMapSize.y);

            // maintain ratio
            /*if(lightMapSize.x > lightMapSize.y){
                lightMapSize.y = ((float)lightportOne_.getColorTexture()->getDimensions().y/(float)lightportOne_.getColorTexture()->getDimensions().x)*lightMapSize.x;
            }
            else if(lightMapSize.x < lightMapSize.y){
                lightMapSize.x = ((float)lightportOne_.getColorTexture()->getDimensions().x/(float)lightportOne_.getColorTexture()->getDimensions().y)*lightMapSize.y;
            }*/

            // clamp and scale coordinates
            proxyMinMaxScreenCoords.x = std::max(std::floor(resultProxyScreenSize.x*(float)screenSize_.x), 0.f);
            proxyMinMaxScreenCoords.y = std::max(std::floor(resultProxyScreenSize.y*(float)screenSize_.y), 0.f);
            proxyMinMaxScreenCoords.z = std::min(std::ceil(resultProxyScreenSize.z*(float)screenSize_.x), (float)screenSize_.x);
            proxyMinMaxScreenCoords.w = std::min(std::ceil(resultProxyScreenSize.w*(float)screenSize_.y), (float)screenSize_.y);
        }

        // calculate propagation zone texture
        if(zoneRenderingOn_.get())
            calculatePropagationZones(projLightPos, entryUnit.getUnitNumber(), exitUnit.getUnitNumber());

        // determine if light source lies in the view of the proxy geometry seen from the camera
        if(twoPassRenderingOn_.get() || zoneRenderingOn_.get()){
            projLightPos *= tgt::vec2(screenSize_);
            dividePassCoords = (light_inside_volume) ? tgt::ivec2(projLightPos) : dividePassCoords;
        }
    }
#endif

    TextureUnit propErtUnit;

    //Modify incoming ERT texture
    if(ertMod_.get()){
        calculateERTmap(entryUnit.getUnitNumber(), exitUnit.getUnitNumber());

        //Bind new exit texture to the unit
        propagationErtPort_.bindColorTexture(propErtUnit);
    }

    // activate and clear output render target
    outport1_.activateTarget();
    outport1_.clearTarget();

    // activate shader and set common uniforms
    raycastPrg_->activate();
    setGlobalShaderParameters(raycastPrg_, &cam);

    // pass texture units to the shader
    raycastPrg_->setUniform("entryPoints_", entryUnit.getUnitNumber());
    raycastPrg_->setUniform("entryPointsDepth_", entryDepthUnit.getUnitNumber());
    entryPort_.setTextureParameters(raycastPrg_, "entryParameters_");

    raycastPrg_->setUniform("exitPoints_", exitUnit.getUnitNumber());
    raycastPrg_->setUniform("exitPointsDepth_", exitDepthUnit.getUnitNumber());
    exitPort_.setTextureParameters(raycastPrg_, "exitParameters_");

    if(ertMod_.get()){
        raycastPrg_->setUniform("ertPoints_", propErtUnit.getUnitNumber());
        propagationErtPort_.setTextureParameters(raycastPrg_, "ertParameters_");
    }

    // bind volume texture and pass it to the shader
    std::vector<VolumeStruct> volumeTextures;
    TextureUnit volUnit;
     volumeTextures.push_back(VolumeStruct(
        volumeInport_.getData(),
        &volUnit,
        "volume_","volumeStruct_",
        true)
    );
    bindVolumes(raycastPrg_, volumeTextures, &cam, lightPosition_.get());

    // bind transfer function and pass it to the shader
    TextureUnit transferUnit;
    if (transferFunc_.get()) {
        transferUnit.activate();
        transferFunc_.get()->getTexture()->bind();
        transferFunc_.get()->setUniform(raycastPrg_, "transferFunc_", "transferFuncTex_", transferUnit.getUnitNumber());
    }

    TextureUnit albedoUnit, surfacenessUnit;
    if (albedo_.get()) {
        albedoUnit.activate();
        albedo_.get()->getTexture()->bind();
        raycastPrg_->setUniform("albedoFunc_", albedoUnit.getUnitNumber());
    }
    if (surfaceness_.get()) {
        surfacenessUnit.activate();
        surfaceness_.get()->getTexture()->bind();
        raycastPrg_->setUniform("surfacenessFunc_", surfacenessUnit.getUnitNumber());
    }
    raycastPrg_->setUniform("ambientIntensity_", ambientIntensity_.get());

#if defined(GL_EXT_direct_state_access) && defined(GL_EXT_shader_image_load_store)
    if(glsl4On_){
        raycastPrg_->setIgnoreUniformLocationError(true);

        // bind light map unit to write
        glBindImageTextureEXT(lightMapUnits_[0],
            lightportOne_.getColorTexture()->getId(), 0, false, 0,  GL_READ_WRITE, GL_RGBA32F_ARB);
        glBindImageTextureEXT(lightMapUnits_[1],
            lightportTwo_.getColorTexture()->getId(), 0, false, 0,  GL_READ_WRITE, GL_RGBA32F_ARB);

        glBindImageTextureEXT(lightMapPosUnits_[0],
            lightportOnePos_.getColorTexture()->getId(), 0, false, 0,  GL_READ_WRITE, GL_RGBA32F_ARB);
        glBindImageTextureEXT(lightMapPosUnits_[1],
            lightportTwoPos_.getColorTexture()->getId(), 0, false, 0,  GL_READ_WRITE, GL_RGBA32F_ARB);

        glProgramUniform1ivEXT(raycastPrg_->getCurrentProgram(),
            raycastPrg_->getUniformLocation("lightMaps_"), 2, lightMapUnits_);

        glProgramUniform1ivEXT(raycastPrg_->getCurrentProgram(),
            raycastPrg_->getUniformLocation("lightPosMaps_"), 2, lightMapPosUnits_);

#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
        readWriteValues_ = (pingPong_.get()) ? tgt::ivec2(0, 1) : tgt::ivec2(0, 0);
#endif

        raycastPrg_->setUniform("lightMapReadXWriteY", readWriteValues_);

        lightportOne_.setTextureParameters(raycastPrg_, "lightMapParameters_");

        raycastPrg_->setUniform("lineWidth_", lineWidth_.get());
        raycastPrg_->setUniform("shadowBrightness_", shadowBrightness_.get());
        raycastPrg_->setUniform("propagationDirection_", propigationAxis_);
        raycastPrg_->setUniform("objToLightMat_", lightCamera_.getProjectionMatrix(outport1_.getSize()) * lightCamera_.getViewMatrix() * volumeInport_.getData()->getPhysicalToWorldMatrix());

        raycastPrg_->setUniform("lightMapSize_", lightMapSize);

        //bind propagation zone port
        if(zoneRenderingOn_.get()){
            TextureUnit zoneUnit, renderTUnit;
            propagationZonesPort_.bindColorTexture(zoneUnit);
            raycastPrg_->setUniform("zoneMap_", zoneUnit.getUnitNumber());
            propagationZonesPort_.setTextureParameters(raycastPrg_, "zoneMapParameters_");

            outport1_.bindColorTexture(renderTUnit);
            raycastPrg_->setUniform("renderTarget_", renderTUnit.getUnitNumber());
            outport1_.setTextureParameters(raycastPrg_, "renderTargetParameters_");

            raycastPrg_->setUniform("projLightPos_", tgt::ivec2(projLightPos));
        }
    }
    else
        raycastPrg_->setIgnoreUniformLocationError(false);
#endif

    bool aborted = false;

    plfFunc_ = &IlluminationLineRaycaster::postProcessLineGL4;

    // render screen aligned quad
    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);

    MatStack.pushMatrix();

    if(!glsl4On_ || renderQuad_.get())
        glCallList(quadList_);
    else if(zoneRenderingOn_.get())
        aborted = drawZones(tgt::ivec4(proxyMinMaxScreenCoords), tgt::ivec2(projLightPos), lineWidth_.get());
#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
    else if(vboOn_.get()){
        if(linesOnlyAboveProxyGeometry_.get()){
            aborted = drawLinesVBO(tgt::ivec4(proxyMinMaxScreenCoords), dividePassCoords, lineWidth_.get());
        }
        else{
            aborted = drawLinesVBO(tgt::ivec4(tgt::ivec2(0), lightportOne_.getSize() - 1), dividePassCoords, lineWidth_.get());
        }
    }
    else{
        if(linesOnlyAboveProxyGeometry_.get()){
            aborted = drawLines(tgt::ivec4(proxyMinMaxScreenCoords), dividePassCoords, lineWidth_.get());
        }
        else{
            aborted = drawLines(tgt::ivec4(tgt::ivec2(0), lightportOne_.getSize() - 1), dividePassCoords, lineWidth_.get());
        }
    }
#else
    else
        aborted = drawLines(tgt::ivec4(proxyMinMaxScreenCoords), dividePassCoords, lineWidth_.get());
#endif

    MatStack.popMatrix();

    glDepthMask(GL_TRUE);
    glEnable(GL_DEPTH_TEST);

    raycastPrg_->deactivate();

    // clean up
    outport1_.deactivateTarget();

#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
    if(aborted && pingPong_.get() && readWriteValues_.y == 0)
        copyRenderPort(&lightportTwo_, &lightportOne_);

    // visualize proxy geometry size by drawing a screen aligned lined square that covers this size
    if(adjustToProxyGeometry_.get() && visualizeProxySize_.get())
        visualizeProxyGeometrySize(&outport1_, &entryPort_, entryUnit.getUnitNumber(), proxyMinMaxScreenCoords);
#endif
    if(aborted) {} // prevent unused warning

    TextureUnit::setZeroUnit();
    LGL_ERROR;
}

std::string IlluminationLineRaycaster::generateHeader(const tgt::GpuCapabilities::GlVersion*) {
    std::string header;

    if((noLinesOverride_ && glslVersion_.getValue() == tgt::GpuCapabilities::GlVersion::SHADER_VERSION_400)
        || (renderQuad_.get() && glslVersion_.getValue() == tgt::GpuCapabilities::GlVersion::SHADER_VERSION_400)){
        const tgt::GpuCapabilities::GlVersion version = tgt::GpuCapabilities::GlVersion::SHADER_VERSION_330;
        header = VolumeRaycaster::generateHeader(&version);
    }
    else {
        const tgt::GpuCapabilities::GlVersion version = glslVersion_.getValue();
        header = VolumeRaycaster::generateHeader(&version);
    }

    header += transferFunc_.get()->getShaderDefines();

    if(zoneRenderingOn_.get())
        header += "#define ZONE_RENDERING\n";

    if(accumulatedRayVis_.get())
        header += "#define VIS_ACCUMULATED_RAY\n";

    if(colorContribution_.get())
        header += "#define COLOR_CONTRIBUTION\n";

    if(piecewise_.get())
        header += "#define PIECEWISE_INTEGRATION\n";

    if(interpolation_.getValue() == 1)
        header += "#define BILINEAR_INTERPOLATION\n";

    if(ertMod_.get())
        header += "#define USE_EARLY_RAY_TERMINATION_SWEEP\n";

    return header;
}

void IlluminationLineRaycaster::updateHeader() {
    raycastPrg_->setHeaders(generateHeader());
    raycastPrg_->rebuild();
}

} // namespace voreen
