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

#ifndef VRN_ORTHOGONALSLICERENDERER_H
#define VRN_ORTHOGONALSLICERENDERER_H

#include "voreen/core/processors/volumerenderer.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/interaction/camerainteractionhandler.h"
#include "voreen/core/datastructures/volume/slice/slicehelper.h"
#include "tgt/camera.h"

#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/fontproperty.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/renderport.h"
#include "voreen/core/ports/genericcoprocessorport.h"

namespace voreen {

class FlowOrthogonalSliceRenderer : public VolumeRenderer {
public:
    FlowOrthogonalSliceRenderer();
    virtual ~FlowOrthogonalSliceRenderer();

    virtual Processor* create() const;
    virtual std::string getCategory() const { return "Slice Rendering"; }
    virtual std::string getClassName() const { return "FlowOrthogonalSliceRenderer"; }
    virtual CodeState getCodeState() const { return CODE_STATE_STABLE; }

    /**
     * The processor is ready when its inport and outports are connected. The
     * coprocessor ports are optional and can be ignored.
     */
    virtual bool isReady() const;
    virtual void process();

    const tgt::ivec3& getSlicePositions() const;

protected:
    enum SliceAlignment { SLICE_XY, SLICE_XZ, SLICE_ZY };

    enum TextureMode {
        TEXTURE_2D,
        TEXTURE_3D
    };

protected:
    virtual void initialize();
    virtual void deinitialize();

    /**
     * \see Processor::beforeProcess()
     */
    virtual void beforeProcess();

    virtual std::string buildShaderHeader();
    bool rebuildShader();

    static tgt::vec3 permuteComponents(const tgt::vec3& input, const tgt::ivec3& permutation);
    void onSlicePositionChange();
    virtual void renderSlice(const SliceAlignment& sliceAlign, const float sliceNo);
    void setSlicePropertiesVisible(const bool visible);
    void setupSliceProperties(const tgt::ivec3& positions);
    virtual void updateNumSlices();

protected:
    virtual void setDescriptions() {
        setDescription("Renders slices orthogonal to each other and aligned to the x-, y- or/and z-axis.");
    }

    RenderPort outport_;
    VolumePort inport_;

    TransFunc1DKeysProperty transferFunc_;

    OptionProperty<TextureMode> texMode_;     ///< use 2D slice textures or 3D volume texture?

    tgt::Shader* sliceShader_;

    tgt::ivec3 volumeDimensions_;
    tgt::ivec3 slicePositions_;

    BoolProperty useXYSliceProp_;
    BoolProperty useXZSliceProp_;
    BoolProperty useZYSliceProp_;
    IntProperty sliceNoXYProp_;
    IntProperty sliceNoXZProp_;
    IntProperty sliceNoZYProp_;
    CameraProperty camProp_;
    CameraInteractionHandler* cameraHandler_;

    GenericCoProcessorPort<FlowOrthogonalSliceRenderer> cpInport_;
    GenericCoProcessorPort<FlowOrthogonalSliceRenderer> cpOutport_;
};

}   // namespace

#endif
