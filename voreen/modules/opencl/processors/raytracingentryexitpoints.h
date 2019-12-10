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

#ifndef VRN_RAYTRACINGENTRYEXITPOINTS_H
#define VRN_RAYTRACINGENTRYEXITPOINTS_H

#include "voreen/core/processors/volumerenderer.h"
#include "voreen/core/interaction/camerainteractionhandler.h"

#include "voreen/core/datastructures/geometry/meshlistgeometry.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"

#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/boolproperty.h"

#include "tgt/shadermanager.h"
#include "modules/opencl/utils/clwrapper.h"
#include "modules/opencl/processors/openclprocessor.h"

namespace voreen {

/**
 * Calculates the entry and exit points for GPU raycasting.
 */
class VRN_CORE_API RaytracingEntryExitPoints : public cl::OpenCLProcessor<RenderProcessor> {
public:
    RaytracingEntryExitPoints();
    virtual ~RaytracingEntryExitPoints();

    virtual std::string getClassName() const    { return "RaytracingEntryExitPoints"; }
    virtual std::string getCategory() const     { return "EntryExitPoints"; }
    virtual CodeState getCodeState() const      { return CODE_STATE_EXPERIMENTAL; }

    virtual Processor* create() const;
    virtual bool isReady() const;

    virtual void beforeProcess();

    virtual void initializeCL();
    virtual void deinitializeCL();

protected:
    virtual void setDescriptions() {
        setDescription("Raytraces the Proxygeometry to create Entry/Exit points. The generated image color-codes the ray parameters for a subsequent VolumeRaycaster. See CubeMeshGeometry");
    }

    virtual void process();

    // processor properties
    BoolProperty useFloatRenderTargets_;
    CameraProperty camera_;  ///< camera used for rendering the proxy geometry

    // interaction handlers
    CameraInteractionHandler* cameraHandler_;

    cl::OpenCL* opencl_;
    cl::Context* context_;
    cl::Device device_;
    cl::CommandQueue* queue_;
    cl::Program* prog_;

    // ports
    RenderPort entryPort_;
    RenderPort exitPort_;
    GeometryPort inport_;

    /// Category used for logging.
    static const std::string loggerCat_;
};

} // namespace voreen

#endif //VRN_RAYTRACINGENTRYEXITPOINTS_H
