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

#include "voreen/core/ports/renderport.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/datastructures/rendertarget/rendertarget.h"
#include "voreen/core/processors/renderprocessor.h"
#include "tgt/filesystem.h"

#ifdef VRN_MODULE_DEVIL
#include <IL/il.h>
#include <IL/ilu.h>
#include "modules/devil/devilmodule.h"
#endif

namespace voreen {

using namespace tgt;

const std::string RENDERSIZE_ORIGIN_PROPERTY_ID = "renderSizeOrigin";
const std::string RENDERSIZE_RECEIVE_PROPERTY_ID = "renderSizeReceive";

const std::string RenderPort::loggerCat_("voreen.RenderPort");

RenderSizeOriginProperty::RenderSizeOriginProperty(const std::string& id, const std::string& guiText,
      const tgt::ivec2& value, int invalidationLevel/*=Processor::INVALID_RESULT*/)
    : IntVec2Property(id, guiText, value, tgt::ivec2(2), tgt::ivec2(8192), invalidationLevel)
{}
RenderSizeOriginProperty::RenderSizeOriginProperty()
    : IntVec2Property("", "", tgt::ivec2(0), tgt::ivec2(2), tgt::ivec2(8192), Processor::INVALID_RESULT)
{}

RenderSizeReceiveProperty::RenderSizeReceiveProperty(const std::string& id, const std::string& guiText, const tgt::ivec2& value,
      int invalidationLevel/*=Processor::INVALID_RESULT*/)
    : IntVec2Property(id, guiText, value, tgt::ivec2(2), tgt::ivec2(8192), invalidationLevel)
{}
RenderSizeReceiveProperty::RenderSizeReceiveProperty()
    : IntVec2Property("", "", tgt::ivec2(0), tgt::ivec2(2), tgt::ivec2(8192), Processor::INVALID_RESULT)
{}

bool LinkEvaluatorRenderSize::arePropertiesLinkable(const Property* p1, const Property* p2) const  {
    return ( (dynamic_cast<const RenderSizeOriginProperty*>(p1) && dynamic_cast<const RenderSizeOriginProperty*>(p2)) ||
             (dynamic_cast<const RenderSizeReceiveProperty*>(p1) && dynamic_cast<const RenderSizeReceiveProperty*>(p2)) ||
             (dynamic_cast<const RenderSizeOriginProperty*>(p1) && dynamic_cast<const RenderSizeReceiveProperty*>(p2)) );
}

void LinkEvaluatorRenderSize::eval(Property* src, Property* dst) {
    tgtAssert(arePropertiesLinkable(src,dst),"RenderSizeLink between wrong propery types");

    IntVec2Property* src_ = dynamic_cast<IntVec2Property*>(src);
    IntVec2Property* dst_ = dynamic_cast<IntVec2Property*>(dst);

    dst_->set(src_->get());
}

//-----------------------------------------------------------------------------

RenderPort::RenderPort(PortDirection direction, const std::string& id, const std::string& guiName,
                       bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel,
                       RenderSizePropagation renderSizePropagation,
                       GLint internalColorFormat, GLint internalDepthFormat)
    : Port(direction, id, guiName, allowMultipleConnections, invalidationLevel)
    , renderTarget_(0)
    , validResult_(false)
    , cleared_(false)
    , size_(128,128)
    , renderSizePropagation_(renderSizePropagation)
    , internalColorFormat_(internalColorFormat)
    , internalDepthFormat_(internalDepthFormat)
    , renderTargetSharing_(false)
    , deinitializeOnDisconnect_(true)
{
    if (renderSizePropagation_ == RENDERSIZE_ORIGIN) {
        if (direction == INPORT) {
            RenderSizeOriginProperty* originProp = new RenderSizeOriginProperty(id + "." + RENDERSIZE_ORIGIN_PROPERTY_ID,"Render Size Origin");
            originProp->setGroupID(id);
            addProperty(originProp);
            setPropertyGroupGuiName(id, (isInport() ? "Inport: " : "Outport: ") + guiName);
            originProp->setVisibleFlag(false);
        }
        else {
            LERROR("Render size propagation mode 'RENDERSIZE_ORIGIN' only allowed for inports (" << getQualifiedName() << ")");
        }
    }

    if (renderSizePropagation_ == RENDERSIZE_RECEIVER) {
        if (direction == INPORT) {
            LERROR("Render size propagation mode 'RENDERSIZE_RECEIVER' only allowed for outports (" << getQualifiedName() << ")");
        }
        else {
            RenderSizeReceiveProperty* receiveProp = new RenderSizeReceiveProperty(id + "." + RENDERSIZE_RECEIVE_PROPERTY_ID,"Render Size Receive");
            receiveProp->setGroupID(id);
            addProperty(receiveProp);
            setPropertyGroupGuiName(id, (isInport() ? "Inport: " : "Outport: ") + guiName);
            receiveProp->setVisibleFlag(false);
        }
    }

}

RenderPort::~RenderPort() {
    if (renderTarget_ && !renderTargetSharing_)
        LERROR("~RenderPort(): '" << getID()
                << "' has not been deinitialized before destruction");

    RenderSizeOriginProperty* originProp = dynamic_cast<RenderSizeOriginProperty*>(getProperty(getID() + "." + RENDERSIZE_ORIGIN_PROPERTY_ID));
    if (originProp) {
        removeProperty(originProp);
        delete originProp;
    }

    RenderSizeReceiveProperty* receiveProp = dynamic_cast<RenderSizeReceiveProperty*>(getProperty(getID() + "." + RENDERSIZE_RECEIVE_PROPERTY_ID));
    if (receiveProp) {
        removeProperty(receiveProp);
        delete receiveProp;
    }
}

std::string RenderPort::getContentDescription() const {
    std::stringstream strstr;
    strstr  << Port::getContentDescription();
    if(hasData())
        strstr << std::endl << "Size: " << getSize().x << " x " << getSize().y;
    return strstr.str();
}

std::string RenderPort::getContentDescriptionHTML() const {
    std::stringstream strstr;
    strstr  << Port::getContentDescriptionHTML();
    if(hasData())
        strstr << "<br>" << "Size: " << getSize().x << " x " << getSize().y;
    return strstr.str();
}

void RenderPort::setProcessor(Processor* p) {
    Port::setProcessor(p);

    RenderProcessor* rp = dynamic_cast<RenderProcessor*>(p);
    tgtAssert(rp, "RenderPort attached to processor of wrong type (RenderProcessor expected)");
    if (!rp)
        LERROR("RenderPort attached to processor of wrong type (RenderProcessor): " << p->getID() << "." << getID());
}

void RenderPort::initialize() {

    Port::initialize();

    if (!isOutport())
        return;

    if (renderTargetSharing_)
        return;

    renderTarget_ = new RenderTarget();
    renderTarget_->initialize(internalColorFormat_, internalDepthFormat_);

    tgtAssert(processor_, "Not attached to processor!");
    renderTarget_->setDebugLabel(processor_->getID()+ "::" + getID());
    renderTarget_->resize(size_);
    clear(); // Clear the potentially uninitialized texture(s)
    validResult_ = false;
    LGL_ERROR;
}

void RenderPort::deinitialize() {
    if (isOutport() && renderTarget_ && !renderTargetSharing_) {
        renderTarget_->deinitialize();
        delete renderTarget_;
        renderTarget_ = 0;
    }
    LGL_ERROR;

    Port::deinitialize();
}

void RenderPort::activateTarget(const std::string& debugLabel) {
    if (isOutport()) {
        if (renderTarget_) {
            renderTarget_->activateTarget(processor_->getID()+ ":" + getID()
                + (debugLabel.empty() ? "" : ": " + debugLabel));
            validateResult();
        }
        else
            LERROR("Trying to activate RenderPort without RenderTarget (" <<
            processor_->getID() << ":" << getID() << ")");
    }
    else {
        if (getRenderTarget()) {
            getRenderTarget()->activateTarget(processor_->getID()+ ":" + getID()
                + (debugLabel.empty() ? "" : ": " + debugLabel));
            //validateResult();
        }
        else
            LERROR("Trying to activate RenderPort without RenderTarget (" <<
            processor_->getID() << ":" << getID() << ")");
        //LERROR("activateTarget() called on inport (" <<
            //processor_->getID() << ":" << getID() << ")");
    }

    cleared_ = false;
}

void RenderPort::deactivateTarget() {
    if (isOutport()) {
        if (renderTarget_){
            renderTarget_->deactivateTarget();
            invalidatePort();
        }
        else
            LERROR("Trying to activate RenderPort without RenderTarget");
    }
    else {
        if(getRenderTarget())
            getRenderTarget()->deactivateTarget();
        else
            LERROR("deactivateTarget() called on inport without RenderTarget");
    }

    cleared_ = false;
}

bool RenderPort::isActive() const {
    return (renderTarget_ && renderTarget_->isActive());
}

void RenderPort::clearTarget() {
    if (!isOutport())
        LERROR("clearTarget() called on inport");
    else if (!isActive())
        LERROR("clearTarget() called on inactive outport");
    else
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void RenderPort::clear() {
    if (!isOutport())
        LERROR("clear() called on inport");
    else {
        tgtAssert(!isActive(), "clear() called on active render target");
        if (hasRenderTarget() && !cleared_) {
            activateTarget();
            clearTarget();
            invalidateResult();
            deactivateTarget();
            cleared_ = true;
        }
    }
}

void RenderPort::changeFormat(GLint internalColorFormat, GLint internalDepthFormat) {
    if (!isOutport()) {
        LERROR("changeFormat() called on inport");
        return;
    }

    tgt::ivec2 s = getSize();
    if (renderTarget_) {
        renderTarget_->deinitialize();
        delete renderTarget_;
    }
    renderTarget_ = new RenderTarget();
    renderTarget_->initialize(internalColorFormat, internalDepthFormat);
    renderTarget_->resize(s);
    invalidatePort();

    internalColorFormat_ = internalColorFormat;
    internalDepthFormat_ = internalDepthFormat;
}

bool RenderPort::hasValidResult() const {
    if (isOutport())
        return renderTarget_ && validResult_;
    else { // inport
        if (!isConnected())
            return false;

        // first connected port is authoritative
        if (RenderPort* p = dynamic_cast<RenderPort*>(connectedPorts_[0]))
            return p->hasValidResult();
        else {
            LERROR("RenderPort is connected to Non-RenderPort");
            return false;
        }
    }
}

void RenderPort::validateResult() {
    if (isOutport()) {
        if (renderTarget_) {
            validResult_ = true;
            renderTarget_->increaseNumUpdates();
        }
        else
            LERROR("validateResult(): no RenderTarget");
    }
    else {
        LERROR("validateResult() called on inport");
    }
}

tgt::ivec2 RenderPort::getSize() const {
    if (hasRenderTarget())
        return getRenderTarget()->getSize();
    else
        return tgt::ivec2(0);
}

void RenderPort::invalidateResult() {
    if (isOutport())
        validResult_ = false;
    else
        LERROR("invalidateResult() called on inport");
}

bool RenderPort::isReady() const {
    bool validInport = isInport() && hasValidResult();
    bool validOutport = isOutport() && hasRenderTarget();
    return (isConnected() && (validInport || validOutport));
}

void RenderPort::setTextureParameters(tgt::Shader* shader, const std::string& uniform) {
    tgtAssert(shader, "Null pointer passed");
    if (hasRenderTarget()) {
        bool oldIgnoreError = shader->getIgnoreUniformLocationError();
        shader->setIgnoreUniformLocationError(true);
        shader->setUniform(uniform + ".dimensions_", tgt::vec2(getSize()));
        shader->setUniform(uniform + ".dimensionsRCP_", tgt::vec2(1.f) / tgt::vec2(getSize()));
        shader->setUniform(uniform + ".matrix_", tgt::mat4::identity);
        shader->setIgnoreUniformLocationError(oldIgnoreError);
    } else {
        LWARNING("setTextureParameter called on RenderPort without RenderTarget. No parameters have been set.");
    }
}

bool RenderPort::connect(Port* inport) {
    if (Port::connect(inport)) {
        RenderPort* rp = static_cast<RenderPort*>(inport);
        return true;
    }
    return false;
}

void RenderPort::disconnect(Port* other) {
    Port::disconnect(other);

    RenderPort* rp = static_cast<RenderPort*>(other);
    if (isOutport()) {
        other->invalidatePort();
    }
}

void RenderPort::resize(const tgt::ivec2& newsize) {
    if (isOutport()) {
        if (size_ == newsize)
            return;
        if (newsize == tgt::ivec2(0)) {
            LWARNING("resize(): invalid size " << newsize);
            return;
        }
        if (renderTarget_) {
            renderTarget_->resize(newsize);
            activateTarget();
            clearTarget();
            deactivateTarget();
        }
        validResult_ = false;
        size_ = newsize;
    }
    else {
        LERROR("resize() called on render inport: " << getQualifiedName());
        //size_ = newsize;
    }
}

void RenderPort::resize(int x, int y) {
    resize(tgt::ivec2(x, y));
}

void RenderPort::bindColorTexture() {
    if (getRenderTarget())
        getRenderTarget()->bindColorTexture();
}

void RenderPort::bindColorTexture(GLint texUnit, GLint filterMode/* = GL_LINEAR*/, GLint wrapMode /*= GL_CLAMP_TO_EDGE*/, tgt::vec4 borderColor /*= tgt::vec4(0.f)*/) {
    if (getRenderTarget())
        getRenderTarget()->bindColorTexture(texUnit, filterMode, wrapMode, borderColor);
}

void RenderPort::bindColorTexture(tgt::TextureUnit& texUnit, GLint filterMode /*= GL_LINEAR*/, GLint wrapMode /*= GL_CLAMP_TO_EDGE*/, tgt::vec4 borderColor /*= tgt::vec4(0.f)*/) {
    bindColorTexture(texUnit.getEnum(), filterMode, wrapMode, borderColor);
}

void RenderPort::bindDepthTexture() {
    if (getRenderTarget())
        getRenderTarget()->bindDepthTexture();
}

void RenderPort::bindDepthTexture(GLint texUnit, GLint filterMode /*= GL_LINEAR*/, GLint wrapMode /*= GL_CLAMP_TO_EDGE*/, tgt::vec4 borderColor /*= tgt::vec4(0.f)*/) {
    if (getRenderTarget())
        getRenderTarget()->bindDepthTexture(texUnit, filterMode, wrapMode, borderColor);
}

void RenderPort::bindDepthTexture(tgt::TextureUnit& texUnit, GLint filterMode /*= GL_LINEAR*/, GLint wrapMode /*= GL_CLAMP_TO_EDGE*/, tgt::vec4 borderColor /*= tgt::vec4(0.f)*/) {
    bindDepthTexture(texUnit.getEnum(), filterMode, wrapMode, borderColor);
}

void RenderPort::bindTextures(GLint colorUnit, GLint depthUnit, GLint filterMode /*= GL_LINEAR*/, GLint wrapMode /*= GL_CLAMP_TO_EDGE*/, tgt::vec4 borderColor /*= tgt::vec4(0.f)*/) {
    bindColorTexture(colorUnit, filterMode, wrapMode, borderColor);
    bindDepthTexture(depthUnit, filterMode, wrapMode, borderColor);
}

void RenderPort::bindTextures(tgt::TextureUnit& colorUnit, tgt::TextureUnit& depthUnit, GLint filterMode /*= GL_LINEAR*/, GLint wrapMode /*= GL_CLAMP_TO_EDGE*/, tgt::vec4 borderColor /*= tgt::vec4(0.f)*/) {
    bindColorTexture(colorUnit, filterMode, wrapMode, borderColor);
    bindDepthTexture(depthUnit, filterMode, wrapMode, borderColor);
}

const tgt::Texture* RenderPort::getColorTexture() const {
    if (hasRenderTarget())
        return getRenderTarget()->getColorTexture();
    else
        return 0;
}

tgt::Texture* RenderPort::getColorTexture() {
    if (hasRenderTarget())
        return getRenderTarget()->getColorTexture();
    else
        return 0;
}

const tgt::Texture* RenderPort::getDepthTexture() const {
    if (hasRenderTarget())
        return getRenderTarget()->getDepthTexture();
    else
        return 0;
}

tgt::Texture* RenderPort::getDepthTexture() {
    if (hasRenderTarget())
        return getRenderTarget()->getDepthTexture();
    else
        return 0;
}

#ifdef VRN_MODULE_DEVIL

void RenderPort::saveToImage(const std::string& filename) {

    if (filename.empty())
        throw VoreenException("filename is empty");
    else if (FileSystem::fileExtension(filename).empty())
        throw VoreenException("filename has no extension");

    // get color buffer content
    tgt::Vector4<uint16_t>* colorBuffer = readColorBuffer<uint16_t>();
    tgt::ivec2 size = getSize();

    // create Devil image from image data and write it to file
    ILuint img;
    ilGenImages(1, &img);
    ilBindImage(img);
    // put pixels into IL-Image
    ilTexImage(size.x, size.y, 1, 4, IL_RGBA, IL_UNSIGNED_SHORT, colorBuffer);
    ilEnable(IL_FILE_OVERWRITE);
    ilResetWrite();
    ILboolean success = ilSaveImage(const_cast<char*>(filename.c_str()));
    ilDeleteImages(1, &img);

    delete[] colorBuffer;

    if (!success) {
        throw VoreenException(DevILModule::getDevILError());
    }
}

#else

void RenderPort::saveToImage(const std::string& /*filename*/) {
    throw VoreenException("Unable to write rendering to file: Voreen was compiled without Devil module.");
}

#endif // VRN_MODULE_DEVIL

void RenderPort::setSharedRenderTarget(RenderTarget* renderTarget) {
    if (isOutport()) {
        if(renderTargetSharing_)
            renderTarget_ = renderTarget;
        else {
            setRenderTargetSharing(true);
            renderTarget_ = renderTarget;
        }

        invalidatePort();
    }
    else {
        LERROR("setRenderTarget() called on inport");
    }
}

const RenderTarget* RenderPort::getRenderTarget() const {
    if (isOutport())
        return renderTarget_;
    else {
        const std::vector<const Port*> connectedPorts = getConnected();
        // first connected port is authoritative
        for (size_t i = 0; i < connectedPorts.size(); ++i) {
            if (!connectedPorts[i]->isOutport())
                continue;
            else if (const RenderPort* p = dynamic_cast<const RenderPort*>(connectedPorts[i]))
                return p->getRenderTarget();
        }
    }
    return 0;
}

RenderTarget* RenderPort::getRenderTarget() {
    // call the const version of getRenderTarget and remove the const from the result
    return const_cast<RenderTarget*>(
        static_cast<const RenderPort*>(this)->getRenderTarget());
}

bool RenderPort::hasRenderTarget() const {
    return (getRenderTarget() != 0);
}

void RenderPort::setRenderTargetSharing(bool sharing) {
    if(sharing == renderTargetSharing_)
        return;
    if(sharing && renderTarget_){
        renderTarget_->deinitialize();
        delete renderTarget_;
        renderTarget_ = 0;
    } else {
        renderTarget_ = 0;
    }
    renderTargetSharing_ = sharing;
}

bool RenderPort::getRenderTargetSharing() const {
    return renderTargetSharing_;
}

void RenderPort::setDeinitializeOnDisconnect(bool deinitializeOnDisconnect) {
    if(isInport())
        LWARNING("Called setDeinitializeOnDisconnect() on inport!");

    deinitializeOnDisconnect_ = deinitializeOnDisconnect;
}

bool RenderPort::getDeinitializeOnDisconnect() const {
    if(isInport())
        LWARNING("Called getDeinitializeOnDisconnect() on inport!");

    return deinitializeOnDisconnect_;
}

bool RenderPort::hasData() const {
    return hasValidResult();
}

tgt::col3 RenderPort::getColorHint() const {
    return tgt::col3(0, 0, 255);
}

RenderPort::RenderSizePropagation RenderPort::getRenderSizePropagation() const {
    return renderSizePropagation_;
}

RenderSizeOriginProperty* RenderPort::getSizeOriginProperty() const {
    return dynamic_cast<RenderSizeOriginProperty*>(getProperty(getID() + "." + RENDERSIZE_ORIGIN_PROPERTY_ID));
}

RenderSizeReceiveProperty* RenderPort::getSizeReceiveProperty() const {
    return dynamic_cast<RenderSizeReceiveProperty*>(getProperty(getID() + "." + RENDERSIZE_RECEIVE_PROPERTY_ID));
}

void RenderPort::requestSize(const tgt::ivec2& size) {
    if (!isInport()) {
        LERROR("requestSize() called on outport: " << getQualifiedName());
        return;
    }
    if (getRenderSizePropagation() != RENDERSIZE_ORIGIN) {
        LERROR("requestSize() called on inport that is not a port size origin: " << getQualifiedName());
        return;
    }

    tgtAssert(getSizeOriginProperty(), "render size origin without SizeOriginProperty");
    getSizeOriginProperty()->set(size);
}

tgt::ivec2 RenderPort::getReceivedSize() const {
    if (!isOutport()) {
        LERROR("getReceivedSize() called on inport: " << getQualifiedName());
        return tgt::ivec2(0);
    }
    if (getRenderSizePropagation() != RENDERSIZE_RECEIVER) {
        LERROR("getReceivedSize() called on outport that is no port size receiver: " << getQualifiedName());
        return tgt::ivec2(0);
    }

    tgtAssert(getSizeReceiveProperty(), "port size receiver has no SizeReceiveProperty");
    return getSizeReceiveProperty()->get();
}

//-------------------------------------------------------------------------------
// PortGroup

PortGroup::PortGroup(bool ignoreConnectivity) : fbo_(0),
    ignoreConnectivity_(ignoreConnectivity)
{

}

PortGroup::~PortGroup() {
}

void PortGroup::addPort(RenderPort* rp) {
    ports_.push_back(rp);
    reattachTargets();
}

void PortGroup::addPort(RenderPort& rp) {
    addPort(&rp);
}

void PortGroup::initialize() {
    if (fbo_)
        return;

    fbo_ = new tgt::FramebufferObject();
}

void PortGroup::deinitialize() {
    delete fbo_;
    fbo_ = 0;
}

void PortGroup::removePort(RenderPort* rp) {
    ports_.erase(std::find(ports_.begin(), ports_.end(), rp));
    reattachTargets();
}

void PortGroup::removePort(RenderPort& rp) {
    removePort(&rp);
}

bool PortGroup::containsPort(RenderPort* rp) {
    return (std::find(ports_.begin(), ports_.end(), rp) != ports_.end());
}

bool PortGroup::containsPort(RenderPort& rp) {
    return containsPort(&rp);
}

void PortGroup::activateTargets(const std::string& debugLabel) {
    fbo_->activate();
    if (ports_.empty())
        return;

    GLenum* buffers = new GLenum[ports_.size()];

    int count=0;
    for (size_t i=0; i < ports_.size();++i) {
        if (ignoreConnectivity_ || ports_[i]->isConnected()) {
            buffers[count] = static_cast<GLenum>(GL_COLOR_ATTACHMENT0+i);
            ports_[i]->validateResult();
            ports_[i]->getRenderTarget()->setDebugLabel(ports_[i]->getProcessor()->getID() + "::"
                                                + ports_[i]->getID() + (debugLabel.empty() ? "" : ": " + debugLabel));
            count++;
        }
    }

    glDrawBuffers(count, buffers);
    glViewport(0, 0, ports_[0]->getSize().x, ports_[0]->getSize().y);
    delete[] buffers;
}

void PortGroup::deactivateTargets() {
    fbo_->deactivate();
    for (size_t i=0; i < ports_.size();++i) {
        if (ignoreConnectivity_ || ports_[i]->isConnected()) {
            ports_[i]->invalidatePort();
        }
    }
}

void PortGroup::clearTargets() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void PortGroup::reattachTargets() {
    if (!fbo_)
        return;

    fbo_->activate();
    fbo_->detachAll();

    if (ports_.empty())
        return;

    bool hasDepth_ = false;

    for (size_t i=0; i<ports_.size(); ++i) {
        RenderPort* p = ports_[i];

        if (!ignoreConnectivity_ && !p->isConnected())
            continue;

        if (p->getColorTexture())
            fbo_->attachTexture(p->getColorTexture(), static_cast<GLenum>(GL_COLOR_ATTACHMENT0+i));
        if (!hasDepth_ && p->getDepthTexture()) {
            hasDepth_ = true;
            fbo_->attachTexture(p->getDepthTexture(), GL_DEPTH_ATTACHMENT);
        }
    }

    LGL_ERROR;

    if (hasDepth_)
        fbo_->isComplete();
    fbo_->deactivate();
}

void PortGroup::resize(const tgt::ivec2& newsize) {
    for (size_t i=0; i<ports_.size(); ++i) {
        ports_[i]->resize(newsize);
    }
    LGL_ERROR;
}

std::string PortGroup::generateHeader(tgt::Shader* shader) {
    // map ports to render targets
    std::string headerSource;
    int targetidx = 0;
    for (size_t i=0; i<ports_.size(); ++i) {
        std::ostringstream op, num, out;
        op << i;
        num << targetidx;
        out << "FragData" << targetidx;
        if (ignoreConnectivity_ || ports_[i]->isConnected()) {
            headerSource += "#define OP" + op.str() + " " + num.str() + "\n";
            if (tgt::GpuCapabilities::isInited() && GpuCaps.getShaderVersion() >= tgt::GpuCapabilities::GlVersion::SHADER_VERSION_130) {
                if (targetidx > 0)
                    headerSource += "out vec4 " + out.str() + ";\n";
                if (shader)
                    shader->bindFragDataLocation(static_cast<GLuint>(targetidx), out.str());
                LGL_ERROR;
            } else {
                if (targetidx > 0)
                    headerSource += "#define FragData" + num.str() + " gl_FragData[" + num.str() + "]\n";
            }
        }
        targetidx++;
    }

    return headerSource;
}

} // namespace
