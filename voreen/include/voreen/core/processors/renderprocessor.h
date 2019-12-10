/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_RENDERPROCESSOR_H
#define VRN_RENDERPROCESSOR_H

#include <vector>

#include "tgt/shadermanager.h"
#include "tgt/camera.h"

#include "voreen/core/datastructures/rendertarget/rendertarget.h"
#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/renderport.h"
#include "voreen/core/ports/genericport.h"
#include "voreen/core/ports/geometryport.h"


namespace voreen {

/**
 * The base class for all processor classes that render to RenderPorts.
 */
class VRN_CORE_API RenderProcessor : public Processor {
public:

    RenderProcessor();

    /**
     * In addition to the invalidations performed by the Processor base class,
     * the results of all RenderPort outports are invalidated.
     *
     * @see Processor::invalidate
     */
    virtual void invalidate(int inv = INVALID_RESULT);

    /**
     * Returns the registered private render ports of this processor.
     *
     * \sa addPrivateRenderPort
     */
    const std::vector<RenderPort*>& getPrivateRenderPorts() const;

    /// Renders a screen-aligned quad with depth func GL_ALWAYS.
    static void renderQuad();

protected:
    /**
     * Initializes all private RenderPorts.
     */
    virtual void initialize();

    /**
     * Deinitializes all private RenderPorts.
     */
    virtual void deinitialize();

    /**
     * Calls adjustRenderOutportSizes().
     *
     * @see Processor::beforeProcess
     */
    virtual void beforeProcess();

    /**
     * Adjusts the dimensions of all render outports
     * TODO: expain strategy
     */
    virtual void adjustRenderOutportSizes();

    /**
     * Allocate or deallocate RenderTargets in outports.
     */
    void manageRenderTargets();

    /// @todo documentation
    void addPrivateRenderPort(RenderPort* port);
    void addPrivateRenderPort(RenderPort& port);

    /**
     * This generates the header that will be used at the beginning of the shaders. It includes the necessary #defines that
     * are to be used with the shaders.
     * \note If you overwrite this function in a subclass, you have to the call the superclass' function first and
     * append your additions to its result!
     */
    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);

    /**
     * Sets some uniforms potentially needed by every shader.
     * @note This function should be called for every shader before every rendering pass!
     *
     * @param shader the shader to set up
     * @param camera camera whose position is passed to uniform cameraPosition_, also needed for passing matrices
     * @param screenDim dimensions of the render target's viewport. Is no parameter is passed,
     *      the dimensions of the first outport is chosen.
     */
    virtual void setGlobalShaderParameters(tgt::Shader* shader, const tgt::Camera* camera = 0, tgt::ivec2 screenDim = tgt::ivec2(-1));

    static const std::string loggerCat_; ///< category used in logging

private:
    /// The private render ports this processor has. Private ports
    /// are mapped to rendertargets no other processor has access to.
    std::vector<RenderPort*> privateRenderPorts_;

    static GLuint fullscreenVBO__;
    static tgt::Shader* fullscreenShader_;
    static size_t fullScreenShaderInUse_;
};

} // namespace voreen

#endif // VRN_RENDERPROCESSOR_H
