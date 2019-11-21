/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_RENDERPORT_H
#define VRN_RENDERPORT_H

#include "voreen/core/ports/port.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"
#include "voreen/core/properties/link/linkevaluatorbase.h"

#include "tgt/textureunit.h"
#include "tgt/tgt_gl.h"
#include "tgt/gpucapabilities.h"
#include "tgt/texture.h"
#include "tgt/vector.h"

namespace tgt {
    class FramebufferObject;
    class Shader;
}

namespace voreen {

class RenderTarget;

/**
 * Used by a render inport for requesting certain input dimensions of an predecessing render outport via a property link.
 * Serves as source of a render size link, i.e., it is expected to be linked to a render outport's RenderSizeReceiveProperty.
 */
class VRN_CORE_API RenderSizeOriginProperty : public IntVec2Property {
public:
    RenderSizeOriginProperty();
    RenderSizeOriginProperty(const std::string& id, const std::string& guiText, const tgt::ivec2& value = tgt::ivec2(128),
        int invalidationLevel=1);
    virtual ~RenderSizeOriginProperty(){};
    virtual Property* create() const               { return new RenderSizeOriginProperty(); }
    virtual std::string getClassName() const       { return "RenderSizeOriginProperty";     }
    virtual std::string getTypeDescription() const { return "IntVector2";                   }
};

/**
 * Used by a render outport for receiving its rendering dimensions via a property link from an render inport.
 * Serves as destination of a render size link, i.e., it is expected to be linked to a render inport's RenderSizeOriginProperty.
 */
class VRN_CORE_API RenderSizeReceiveProperty : public IntVec2Property {
public:
    RenderSizeReceiveProperty();
    RenderSizeReceiveProperty(const std::string& id, const std::string& guiText, const tgt::ivec2& value = tgt::ivec2(128),
        int invalidationLevel=1);
    virtual ~RenderSizeReceiveProperty(){};
    virtual Property* create() const               { return new RenderSizeReceiveProperty(); }
    virtual std::string getClassName() const       { return "RenderSizeReceiveProperty";     }
    virtual std::string getTypeDescription() const { return "IntVector2";                    }
};

/**
 * Uni-directionally propagates a render size (tgt::ivec2) from a RenderSizeOriginProperty to a RenderSizeReceiveProperty.
 */
class VRN_CORE_API LinkEvaluatorRenderSize : public LinkEvaluatorBase {
public:
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorRenderSize(); }
    virtual std::string getClassName() const  { return "LinkEvaluatorRenderSize";     }
    virtual std::string getGuiName() const    { return "Render Size";                 }

    virtual bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
    virtual void eval(Property* src, Property* dst);
};

//-----------------------------------------------------------------------------

class VRN_CORE_API RenderPort : public Port {

    friend class RenderProcessor;
    friend class StereoCanvasRenderer;
    friend class NetworkEvaluator;

public:
    /**
     * Determines a port's rendering size propagation behaviour.
     */
    enum RenderSizePropagation {
        RENDERSIZE_DEFAULT,     ///< Port is neither a size origin nor does it receive an output size.
                                ///  If the port is an outport, then its rendering size is automatically set to the processor's input dimensions.
        RENDERSIZE_STATIC,      ///< Port is neither a size origin nor does it receive an output size,
                                ///  and its rendering dimensions are also not propagated from the processor's inports.
        RENDERSIZE_ORIGIN,      ///< Port requests an input size via a RenderSizeOriginProperty (inports only).
        RENDERSIZE_RECEIVER     ///< Port receives its output dimensions via a RenderSizeReceiveProperty (outports only).
    };

    RenderPort(PortDirection direction, const std::string& id, const std::string& guiName = "", bool allowMultipleConnections = false,
               Processor::InvalidationLevel invalidationLevel = Processor::INVALID_RESULT,
               RenderSizePropagation renderSizePropagation = RENDERSIZE_DEFAULT,
               GLint internalColorFormat = GL_RGBA16, GLint internalDepthFormat = GL_DEPTH_COMPONENT24);
    virtual ~RenderPort();

    virtual std::string getClassName() const { return "RenderPort"; }
    virtual Port* create(PortDirection direction, const std::string& id, const std::string& guiName = "") const {return new RenderPort(direction,id,guiName);}
    virtual std::string getContentDescription() const;
    virtual std::string getContentDescriptionHTML() const;

    /**
     * @brief Activates the outport's RenderTarget, so that all subsequent rendering operations
     * are performed on this target. The target has to be deactivated after all rendering operations
     * have been finished or before any other RenderPort is activated.
     *
     * @note Must not be called on an inport.
     *
     * @see deactivateTarget
     *
     * @param debugLabel Additional description to be presented in the GUI.
     */
    void activateTarget(const std::string& debugLabel = "");

    /**
     * Deactivates the outport's RenderTarget after rendering
     * has been finished.
     *
     * @note Must not be called on an inport.
     */
    void deactivateTarget();

    /**
     * Clears the contents of an activated outport's RenderTarget,
     * by calling glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT).
     *
     * @note May only be called on an activated outport.
     */
    void clearTarget();

    /**
     * Do not use this method while the RenderPort is activated (use clearTarget() for rendering instead).
     */
    virtual void clear();

    /**
     * Returns true, if the port is connected and (if it is an inport) has a valid rendering.
     */
    virtual bool isReady() const;

    /**
     * Returns true, if the port is an outport and stores a valid rendering,
     * or if is an inport, and its (first) connected outport stores a valid rendering.
     */
    bool hasValidResult() const;

    /**
     * Returns the result of hasValidResult().
     */
    virtual bool hasData() const;

    /**
     * Returns the dimensions of the associated RenderTarget.
     *
     * If no RenderTarget is associated, ivec2(0) is returned.
     */
    tgt::ivec2 getSize() const;

    /**
     * Returns true, if the associated RenderTarget is currently active.
     *
     * @see activateTarget, deactivateTarget
     */
    bool isActive() const;

    /**
     * Returns the port's RenderTarget, may be null.
     *
     * If this function is called on an outport, its own RenderTarget is returned.
     * If it is called on an inport, the RenderTarget of the (first) connected outport is returned.
     *
     * @see hasRenderTarget
     */
    RenderTarget* getRenderTarget();

    /// @overload
    const RenderTarget* getRenderTarget() const;

    /**
     * Returns the associated color texture.
     *
     * @return The color texture or null, if none is present.
     */
    const tgt::Texture* getColorTexture() const;

    /// @overload
    tgt::Texture* getColorTexture();

    /**
     * Returns the associated depth texture.
     *
     * @return The depth texture or null, if none is present.
     */
    const tgt::Texture* getDepthTexture() const;

    /// @overload
    tgt::Texture* getDepthTexture();

    /**
     * Returns true, if a RenderTarget is assigned to the RenderPort.
     *
     * If this function is called on an inport, it checks whether
     * its (first) connected outport has a RenderTarget.
     */
    bool hasRenderTarget() const;

    /**
     * Resizes the associated RenderTarget to the passed dimensions.
     *
     * If this function is called on an inport, the call is propagated to
     * all connected outports.
     */
    void resize(const tgt::ivec2& newsize);

    /**
     * Resizes the associated RenderTarget to the passed dimensions.
     *
     * If this function is called on an inport, the call is propagated to
     * all connected outports.
     */
    void resize(int x, int y);

    /**
     * Re-create and initialize the outport's RenderTarget with the given format.
     * The RenderTarget's resolution is preserved.
     *
     * @note This function must not be called on an inport.
     */
    void changeFormat(GLint internalColorFormat, GLint internalDepthFormat=GL_DEPTH_COMPONENT24);

    /**
     * Flags the outport's rendering as valid.
     *
     * @note It is normally not necessary to call this function,
     *  since it is called implicitly by activateTarget/deactivateTarget.
     *
     * @note Must not be called on an inport.
     */
    void validateResult();

    /**
     * Flags the outport's rendering as invalid.
     *
     * @note It is normally not necessary to call this function,
     *  since it is called implicitly.
     *
     * @note Must not be called on an inport.
     */
    void invalidateResult();

    /**
     * @brief Passes the assigned textures' TextureParameters struct to the passed shader.
     *  Needs to be called for each inport whose textures are to be accessed by the shader.
     *
     * @param uniform The name of the struct used in the shader.
     */
    void setTextureParameters(tgt::Shader* shader, const std::string& uniform);

    /**
     * Bind the color texture of this port's RenderTarget to the currently active texture unit.
     */
    void bindColorTexture();

    /**
     * Bind the color texture of this port's RenderTarget to a specific texture unit, using the specified filter settings.
     *
     * @param texUnit The texture unit to activate before binding.
     */
    void bindColorTexture(GLint texUnit, GLint filterMode = GL_LINEAR, GLint wrapMode = GL_CLAMP_TO_EDGE, tgt::vec4 borderColor = tgt::vec4(0.f));

    /// @overload
    void bindColorTexture(tgt::TextureUnit& texUnit, GLint filterMode = GL_LINEAR, GLint wrapMode = GL_CLAMP_TO_EDGE, tgt::vec4 borderColor = tgt::vec4(0.f));

    /**
     * Bind the depth texture of this port's RenderTarget to the currently active texture unit.
     */
    void bindDepthTexture();

    /**
     * Bind the depth texture of this port's RenderTarget to a specific texture unit, using the specified filter settings.
     *
     * @param texUnit The texture unit to activate before binding.
     */
    void bindDepthTexture(GLint texUnit, GLint filterMode = GL_LINEAR, GLint wrapMode = GL_CLAMP_TO_EDGE, tgt::vec4 borderColor = tgt::vec4(0.f));

    /// @overload
    void bindDepthTexture(tgt::TextureUnit& texUnit, GLint filterMode = GL_LINEAR, GLint wrapMode = GL_CLAMP_TO_EDGE, tgt::vec4 borderColor = tgt::vec4(0.f));

    /**
     * Convenience function, calling bindColorTexture(colorUnit) and bindDepthTexture(depthUnit), passing the specified filter settings.
     */
    void bindTextures(GLint colorUnit, GLint depthUnit, GLint filterMode = GL_LINEAR, GLint wrapMode = GL_CLAMP_TO_EDGE, tgt::vec4 borderColor = tgt::vec4(0.f));

    /// @overload
    void bindTextures(tgt::TextureUnit& colorUnit, tgt::TextureUnit& depthUnit, GLint filterMode = GL_LINEAR, GLint wrapMode = GL_CLAMP_TO_EDGE, tgt::vec4 borderColor = tgt::vec4(0.f));

    /// Returns the rendering size propagation behaviour of the port.
    RenderSizePropagation getRenderSizePropagation() const;

    /**
     * Returns the port's size origin property, if it has one,
     * or 0 otherwise.
     */
    RenderSizeOriginProperty* getSizeOriginProperty() const;

    /**
     * Returns the port's size receiver property, if it has one,
     * or 0 otherwise.
     */
    RenderSizeReceiveProperty* getSizeReceiveProperty() const;

    /**
     * Requests the passed input size by assigning it to the port's RenderSizeOriginProperty.
     *
     * @note Only allowed for an inport that is a RENDERSIZE_ORIGIN.
     */
    void requestSize(const tgt::ivec2& size);

    /**
     * Returns the output size the port has received via its RenderSizeReceiveProperty.
     *
     * @note Only allowed for an outport that is a RENDERSIZE_RECEIVER.
     */
    tgt::ivec2 getReceivedSize() const;

    /**
     * Registers a function that will be called when the value of the port's
     * RenderSizeReceiveProperty changes.
     *
     * @note Only allowed for an outport that is a RENDERSIZE_RECEIVER.
     *
     * @param target The class whose member function is to be called.
     * @param function Function pointer of the function to be called.
     *
     */
    template<typename T>
    void onSizeReceiveChange(T* target, void (T::*function)());

    /**
     * Writes the currently stored rendering to an image file.
     *
     * @note This function requires Voreen to be built with DevIL support.
     *
     * @param filename the filename of the output file. Must have an
     *      extension known by the DevIL library. *.jpg and *.png
     *      should work fine.
     *
     * @throw VoreenException if the image saving failed
     */
    void saveToImage(const std::string &filename);

    /**
     * \brief Returns the current content of the color buffer, converted
     * to RGBA format. The data type is specified as template parameter.
     *
     * @tparam T Specifies the data type of the returned image data.
     *         Supported types: uint8_t, uint16_t, float
     *
     * @note Releasing the allocated memory is up to the caller!
     *
     * @throw VoreenException if the color buffer content could not be read
     *        or converted
     */
    template<typename T>
    tgt::Vector4<T>* readColorBuffer();

    /// @see getDeinitializeOnDisconnect
    void setDeinitializeOnDisconnect(bool deinitializeOnDisconnect);
    /// Deinitialize (delete) the rendertarget stored in this port upon disconnect? This is used to conserve GPU memory, default is true.
    bool getDeinitializeOnDisconnect() const;

protected:
    /**
     * Assigns the passed RenderTarget to this RenderPort
     * and sets the shared render target flag.
     */
    virtual void setSharedRenderTarget(RenderTarget* renderTarget);

    /**
     * Initializes the RenderTarget, if the port is an outport.
     *
     * @see Port::initialize
     */
    virtual void initialize();

    /**
     * Deinitializes the RenderTarget, if the port is an outport.
     *
     * @see Port::deinitialize
     */
    virtual void deinitialize();

    /**
     * Additionally propagates the connected port's size origin.
     *
     * @see Port::connect
     */
    virtual bool connect(Port* inport);

    /**
     * Additionally propagates itself as size origin.
     *
     * @see Port::disconnect
     */
    virtual void disconnect(Port* other);

    /**
     * Additionally checks, if the passed processor is of type RenderProcessor.
     *
     * @see Port::setProcessor
     */
    virtual void setProcessor(Processor* p);

    void setRenderTargetSharing(bool sharing);

    bool getRenderTargetSharing() const;

    virtual tgt::col3 getColorHint() const;

private:
    RenderTarget* renderTarget_;

    bool validResult_;
    bool cleared_;    //< false, if the render buffer is empty
    tgt::ivec2 size_; //neccessary for inports

    RenderSizePropagation renderSizePropagation_;

    GLint internalColorFormat_;
    GLint internalDepthFormat_;

    bool renderTargetSharing_;
    bool deinitializeOnDisconnect_;

    static const std::string loggerCat_; ///< category used in logging
};

template<typename T>
void voreen::RenderPort::onSizeReceiveChange(T* target, void (T::*function)()) {
    if (!isOutport()) {
        LERROR("onSizeReceiveChange() may only be called on an outport: " << getQualifiedName());
        return;
    }
    if (getRenderSizePropagation() != RENDERSIZE_RECEIVER) {
        LERROR("onSizeReceiveChange() may only be called on an outport that is a port size receiver: " << getQualifiedName());
        return;
    }

    tgtAssert(getSizeReceiveProperty(), "port size receiver has no SizeReceiveProperty");
    getSizeReceiveProperty()->onChange(MemberFunctionCallback<T>(target, function));
}

template<typename T>
tgt::Vector4<T>* voreen::RenderPort::readColorBuffer() {
    if (!getColorTexture()) {
        throw VoreenException("RenderPort::readColorBuffer() called on an empty render port");
    }

    // determine OpenGL data type from template parameter
    GLenum dataType;
    if (typeid(T) == typeid(uint8_t))
        dataType = GL_UNSIGNED_BYTE;
    else if (typeid(T) == typeid(uint16_t))
        dataType = GL_UNSIGNED_SHORT;
    else if (typeid(T) == typeid(float))
        dataType = GL_FLOAT;
    else
        throw VoreenException("RenderPort::readColorBuffer(): unsupported data type. "
            "Expected: uint8_t, uint16_t, float");

    GLubyte* pixels = 0;
    try {
        pixels = getColorTexture()->downloadTextureToBuffer(GL_RGBA, dataType);
    }
    catch (std::bad_alloc&) {
        throw VoreenException("RenderPort::readColorBuffer(): bad allocation");
    }
    LGL_ERROR;

    if (pixels)
        return reinterpret_cast<tgt::Vector4<T>*>(pixels);
    else
        throw VoreenException("RenderPort::readColorBuffer(): failed to download texture");
}


//------------------------------------------------------------------------------------------------------


/**
 * @brief This class groups RenderPorts to allow easy rendering to multiple targets.
 */
class VRN_CORE_API PortGroup {
public:
    /**
     * @param ignoreConnectivity If this is true all ports are attached to the FBO, otherwise only the connected ones.
     */
    PortGroup(bool ignoreConnectivity = false);
    ~PortGroup();

    /// Initializes the FBO
    void initialize();

    /// Deinitializes the FBO
    void deinitialize();

    ///Add or remove a port to/from the group.
    void addPort(RenderPort* rp);
    void addPort(RenderPort& rp);
    void removePort(RenderPort* rp);
    void removePort(RenderPort& rp);
    // Check if the group contains the specified port.
    bool containsPort(RenderPort* rp);
    bool containsPort(RenderPort& rp);

    /**
     * @brief Start rendering to all attached ports.
     *
     * @param debugLabel @see RenderPort::activateTarget
     */
    void activateTargets(const std::string& debugLabel = "");

    void deactivateTargets();

    /**
     * Clears the contents of the associated RenderPorts,
     * by calling glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT).
     *
     * @note May only be called on an activated port group.
     */
    void clearTargets();

    ///Resize all attached ports.
    void resize(const tgt::ivec2& newsize);

    /**
     * @brief Defines OP0, OP1, ... OPn to address targets in shader.
     *
     * Example:
     * CPP:
     * addPort(p1); addPort(p2); addPort(p3);
     *
     * glsl:
     * FragData0 = vec4(1.0); //write to p1
     * FragData1 = vec4(1.0); //write to p2
     * FragData2 = vec4(1.0); //write to p3
     *
     * The defines can also be used to test if a port is connected:
     * #ifdef OP0
     *    FragData0 = result;
     * #endif
     *
     * @return A string containing the defines.
     */
    std::string generateHeader(tgt::Shader* shader);

    /**
     * Re-attach all rendertargets to the FBO.
     * If ignoreConnectivity is true only connected ports are attached.
     */
    void reattachTargets();

protected:

    std::vector<RenderPort*> ports_;
    tgt::FramebufferObject* fbo_;
    bool ignoreConnectivity_;
};


} // namespace

#endif // VRN_RENDERPORT_H
