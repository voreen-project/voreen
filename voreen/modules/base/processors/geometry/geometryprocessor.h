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

#ifndef VRN_GEOMETRYPROCESSOR_H
#define VRN_GEOMETRYPROCESSOR_H

#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/processors/geometryrendererbase.h"
#include "voreen/core/interaction/idmanager.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/buttonproperty.h"


namespace voreen {

class CameraInteractionHandler;

/**
 * Renders GeometryRenderer objects.
 *
 * This processor renders all to it connected coprocessors. The projection and modelview matrices
 * are set according to the current viewing position.
 */
class VRN_CORE_API GeometryProcessor : public RenderProcessor {
public:
    GeometryProcessor();
    ~GeometryProcessor();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "GeometryProcessor"; }
    virtual std::string getCategory() const { return "Geometry"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_STABLE; }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Manages the connected GeometryRenderer objects. Holds a vector of GeometryRenderer Objects and renders all of them on <i>render()</i>. The current projection and modelview matrices are propagated to the GeometryRenderer objects.");
        maxFragmentsPerPixel_.setDescription(
                "Controls the size of global/local buffers used for rendering using order independent transparency. "
                "The local buffer size is a hard limit for the number of samples blended during rendering. "
                "The global buffer size is a soft limit in so far as it poses a maximum value for the average number of samples taken from the whole image "
                "and can usually be lower than the local if the rendered objects do not occupy the whole viewport, "
                "as samples from all rendered transparent objects share the same global buffer, which is viewport size depended. "
                "Rendering artifacts can be caused by setting either the global (white regions image regions) or local (other artifacts) buffer size too low. "
                "Too large buffer sizes can impair performance."
                );
    }

    virtual void beforeProcess();
    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);
    virtual void adjustOITProperties();

    // OIT methods
    void blendOITBuffer();
    void renderIntoOITBuffer(RenderPort& p);
    void clearOITBuffers();
    void setupOITBuffers();
    void setupOITShaders();
    void compose(RenderPort& r1, RenderPort& r2);

    void maxFragmentsPerPixelChanged();

private:
    /// Requests the received output render size from the render inport.
    void passThroughSizeRequest();

    /// Updates enclosing Bounding Box and adapts to scene.
    void adaptToScene(bool forced = false);

    // Shaders for OIT
    ShaderProperty oirAddImageShader_;
    ShaderProperty oirBlendShader_;

    // Buffers for OIT
    GLuint headPointerTexture_;
    GLuint headPointerInitializer_;
    GLuint atomicFragmentCounterBuffer_;
    GLuint fragmentStorageBuffer_;
    GLuint fragmentStorageTexture_;

    // Configuration properties for OIT
    IntIntervalProperty maxFragmentsPerPixel_;
    BoolProperty applyOrderIndependentTransparency_;

    IDManager idManager_;

    BoolProperty renderGeometries_;
    BoolProperty adaptToScene_;
    ButtonProperty forceAdaptToScene_;
    CameraProperty camera_;
    CameraInteractionHandler* cameraHandler_;
    ShaderProperty composeShader_;
    BoolProperty useFloatRenderTargets_;

    RenderPort inport_;
    RenderPort outport_;
    RenderPort tempPort_; ///< holds image data if we cannot directly render into outport
    RenderPort oitImageAddTargetPort_; ///< another temp port only used by renderIntoOITBuffer
    RenderPort pickingPort_;
    GenericCoProcessorPort<GeometryRendererBase> cpPort_;

    bool updateSceneAdaptation_;
};

} // namespace voreen

#endif // VRN_GEOMETRYPROCESSOR_H
