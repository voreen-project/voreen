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

#ifndef VRN_FLOWSLICERENDERER3D_H
#define VRN_FLOWSLICERENDERER3D_H

#include "flowslicerenderer.h"

namespace voreen {

class CameraInteractionHandler;

/**
 * A processor for rendering 2D flow images on orthogonal slices in 3D using
 * various 2D flow visualization techniques.
 */
class FlowSliceRenderer3D : public FlowSliceRenderer {
public:
    FlowSliceRenderer3D();
    virtual ~FlowSliceRenderer3D();

    virtual Processor* create() const { return new FlowSliceRenderer3D(); }
    virtual std::string getModule() const { return "flowreen"; }
    virtual std::string getCategory() const { return "Flow Visualization"; }
    virtual std::string getClassName() const { return "FlowSliceRenderer3D"; }
    virtual CodeState getCodeState() const { return CODE_STATE_STABLE; }
    virtual void process();

protected:
    virtual void setDescriptions() {
        setDescription("Renders 2D flow images from slices in 3D vector fields into the volume.");
    }

    virtual std::string generateShaderHeader();

private:
    inline void invalidateXYTexture() { rebuildTextureXY_ = true; }
    inline void invalidateXZTexture() { rebuildTextureXZ_ = true; }
    inline void invalidateZYTexture() { rebuildTextureZY_ = true; }

    static tgt::vec3 permuteComponents(const tgt::vec3& input, const tgt::ivec3& permutation);

    void renderSlice(const Flow3D& flow3D, const float sliceNo, const tgt::vec2& viewportSize,
        const SliceAlignment& sliceAlignment);

    void updateNumSlices();

private:
    BoolProperty useXYSliceProp_;
    BoolProperty useXZSliceProp_;
    BoolProperty useZYSliceProp_;
    IntProperty sliceNoXYProp_;
    IntProperty sliceNoXZProp_;
    IntProperty sliceNoZYProp_;
    CameraProperty camProp_;
    CameraInteractionHandler* cameraHandler_;

    tgt::Texture* textureXY_;
    tgt::Texture* textureXZ_;
    tgt::Texture* textureZY_;

    bool rebuildTextureXY_;
    bool rebuildTextureXZ_;
    bool rebuildTextureZY_;

    FlowSliceRenderer::SliceAlignment alignment_;
};

}   // namespace

#endif  // VRN_FLOWSLICERENDERER3D_H
