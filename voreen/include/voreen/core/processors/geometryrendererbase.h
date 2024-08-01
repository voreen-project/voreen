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

#ifndef VRN_GEOMETRYRENDERERBASE_H
#define VRN_GEOMETRYRENDERERBASE_H

#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/interaction/idmanager.h"
#include "voreen/core/ports/genericcoprocessorport.h"

namespace voreen {

/**
 * @brief Abstract base class for rendering Geometry onto images.
 * Derived processors are expected to be connected to a GeometryProcessor
 * through a co-processor connection.
 */
class VRN_CORE_API GeometryRendererBase : public RenderProcessor {

public:
    GeometryRendererBase();

    virtual bool usesTransparency() const;

    virtual void renderTransparent();
    virtual void render() = 0;
    virtual tgt::Bounds getBoundingBox() const;

    virtual void renderPicking();

    virtual void setCamera(const tgt::Camera& c) { camera_ = c; }
    virtual void setViewport(tgt::ivec2 viewport) { viewport_ = viewport; }
    virtual void setIDManager(IDManager* idm) { idManager_ = idm; }

    virtual const IDManager* getIDManager() const { return idManager_; }

protected:

    virtual void process();

    ///Get 3D vector from screen position (unproject)
    tgt::vec3 getOGLPos(int x, int y, float z) const;

    /**
     * Get screen position from 3D vector (project),
     * using the camera property's view and projection matrices.
     *
     * @param pos the point to project
     * @param modelview the modelview matrix to apply in addition to the camera's matrix.
     *      The camera's matrix is multiplied by the passed one, i.e., the passed one is applied first.
     * @param projection the projection matrix to apply in addition to the camera's matrix.
     *      The camera's matrix is multiplied by the passed one, i.e., the passed one is applied first.
     */
    tgt::vec3 getWindowPos(tgt::vec3 pos, tgt::mat4 modelview = tgt::mat4::identity,
        tgt::mat4 projection = tgt::mat4::identity) const;

    tgt::Camera camera_;
    tgt::ivec2 viewport_;
    GenericCoProcessorPort<GeometryRendererBase> outPort_;
    IDManager* idManager_;
};

} // namespace voreen

#endif // VRN_GEOMETRYRENDERERBASE_H
