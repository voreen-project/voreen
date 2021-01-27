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

#ifndef VRN_SLICEPOSITIONRENDERER_H
#define VRN_SLICEPOSITIONRENDERER_H

#include "voreen/core/processors/geometryrendererbase.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/shaderproperty.h"

#include "voreen/core/ports/volumeport.h"

namespace voreen {

class VRN_CORE_API SlicePositionRenderer : public GeometryRendererBase {
public:
    SlicePositionRenderer();

    Processor* create() const { return new SlicePositionRenderer; }
    std::string getCategory() const  { return "Geometry"; }
    std::string getClassName() const { return "SlicePositionRenderer"; }
    CodeState getCodeState() const   { return CODE_STATE_STABLE; }

    virtual tgt::Bounds getBoundingBox() const;

    virtual bool usesTransparency() const;
    virtual void renderTransparent();
    void render();

    void render(tgt::Shader* planeShader);

    void invalidate(int inv = INVALID_RESULT);

protected:
    virtual void setDescriptions() {
        setDescription("Indicates the position of axis-aligned slices in a 3D view. Can be used when linking multiple views, to show the position of the 2D views with respect to the 3D view.");
    }

    virtual void initialize();
    virtual void deinitialize();
    void process();

    virtual void adjustPropertiesToInput();
    void setVBOData(const std::array<tgt::vec3, 4>& vertices);
private:
    void togglePropertyVisibility();

    //port
    VolumePort inport_;

    BoolProperty enable_;
    //x slice
    BoolProperty renderXSlice_;
    IntProperty xSliceIndexProp_;
    ColorProperty xColor_;
    FloatProperty lineXWidth_;
    FloatProperty alphaFactorXPlane_;
    //y slice
    BoolProperty renderYSlice_;
    IntProperty ySliceIndexProp_;
    ColorProperty yColor_;
    FloatProperty lineYWidth_;
    FloatProperty alphaFactorYPlane_;
    //z slice
    BoolProperty renderZSlice_;
    IntProperty zSliceIndexProp_;
    ColorProperty zColor_;
    FloatProperty lineZWidth_;
    FloatProperty alphaFactorZPlane_;

    ShaderProperty planeShaderOpaque_;
    ShaderProperty planeShaderTransparent_;
    GLuint vao_;
    GLuint vbo_;
};

}

#endif

