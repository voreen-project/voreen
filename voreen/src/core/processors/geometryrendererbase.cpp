/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
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

#include "voreen/core/processors/geometryrendererbase.h"

#include "tgt/glmath.h"
#include "tgt/vector.h"
#include "tgt/matrix.h"
#include "tgt/matrixstack.h"

using tgt::vec4;
using tgt::vec3;

namespace voreen {

GeometryRendererBase::GeometryRendererBase()
    : RenderProcessor()
    , camera_()
    , outPort_(Port::OUTPORT, "coprocessor.geometryrenderer", "GeometryProcessor")
    , idManager_(0)
{
    //addProperty(camera_);
    addPort(outPort_);
}

tgt::vec3 GeometryRendererBase::getOGLPos(int x, int y,float z) const {
    if (!idManager_)
        return tgt::vec3(0.0f);
    // taken from NEHE article 13
    // http://nehe.gamedev.net/data/articles/article.asp?article=13
    GLint viewport[4];
    GLdouble modelview[16];
    GLdouble projection[16];
    GLdouble winX, winY, winZ;
    GLdouble posXh, posYh, posZh;

#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
    tgt::mat4 projection_tgt = tgt::getTransposeProjectionMatrix();
    tgt::mat4 modelview_tgt = tgt::getTransposeModelViewMatrix();
#else
    tgt::mat4 projection_tgt = transpose(MatStack.getProjectionMatrix());
    tgt::mat4 modelview_tgt = transpose(MatStack.getModelViewMatrix());
#endif
    for (int i = 0; i < 4; ++i) {
        modelview[4*i]    = modelview_tgt[i].x;
        modelview[4*i+1]  = modelview_tgt[i].y;
        modelview[4*i+2]  = modelview_tgt[i].z;
        modelview[4*i+3]  = modelview_tgt[i].w;
        projection[4*i]   = projection_tgt[i].x;
        projection[4*i+1] = projection_tgt[i].y;
        projection[4*i+2] = projection_tgt[i].z;
        projection[4*i+3] = projection_tgt[i].w;
    }
    viewport[0] = 0;
    viewport[1] = 0;
    viewport[2] = idManager_->getRenderTarget()->getSize().x;
    viewport[3] = idManager_->getRenderTarget()->getSize().y;

    winX = static_cast<GLdouble>(x);
    winY = static_cast<GLdouble>(viewport[3]) - static_cast<GLint>(y);
    winZ = static_cast<GLdouble>(z);

    tgt::gluUnProject( winX, winY, winZ, modelview, projection, viewport, &posXh, &posYh, &posZh);

    tgt::vec3 returned = tgt::vec3(static_cast<float>(posXh), static_cast<float>(posYh), static_cast<float>(posZh));
    return returned;
}

tgt::vec3 GeometryRendererBase::getWindowPos(tgt::vec3 pos, tgt::mat4 pModelview, tgt::mat4 pProjection) const {
    GLint viewport[4];
    GLdouble modelview[16];
    GLdouble projection[16];

    tgt::mat4 projection_tgt = camera_.getProjectionMatrix(idManager_->getRenderTarget()->getSize()) * pProjection;
    tgt::mat4 modelview_tgt = camera_.getViewMatrix() * pModelview;
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
    viewport[2] = idManager_->getRenderTarget()->getSize().x;
    viewport[3] = idManager_->getRenderTarget()->getSize().y;

    GLdouble pointProjectedGL[3];
    tgt::gluProject(pos.x, pos.y, pos.z, modelview, projection, viewport,
               &pointProjectedGL[0], &pointProjectedGL[1], &pointProjectedGL[2]);

    return tgt::vec3(static_cast<float>(pointProjectedGL[0]),
                     static_cast<float>(pointProjectedGL[1]),
                     static_cast<float>(pointProjectedGL[2]));
}

void GeometryRendererBase::process() {
    // DO NOTHING
}

bool GeometryRendererBase::usesTransparency() const {
    return false;
}
void GeometryRendererBase::renderTransparent() {
    tgtAssert(false, "renderTransparent not implemented.");
}
void GeometryRendererBase::renderPicking() {
    // DO NOTHING
}

tgt::Bounds GeometryRendererBase::getBoundingBox() const {
    // Return an undefined box.
    return tgt::Bounds();
}

} // namespace voreen
