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

#ifndef VRN_SLICEPOINTRENDERER3D_H
#define VRN_SLICEPOINTRENDERER3D_H

#include "voreen/core/processors/geometryrendererbase.h"

#include "voreen/core/ports/volumeport.h"

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/shaderproperty.h"

#include "voreen/core/datastructures/geometry/trianglemeshgeometryindexed.h"

namespace voreen {

    /**
     * Used to render points in 3D. Must be linked to a SlicePointRenderer2D.
     *
     * @see SliceViewer
     * @see SlicePointRenderer2D
     */
class VRN_CORE_API SlicePointRenderer3D : public GeometryRendererBase {
public:
    SlicePointRenderer3D();

    virtual Processor* create() const { return new SlicePointRenderer3D; }
    virtual std::string getCategory() const  { return "Geometry"; }
    virtual std::string getClassName() const { return "SlicePointRenderer3D"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_TESTING; }

    virtual void render();
    virtual void invalidate(int inv = INVALID_RESULT);

    virtual tgt::Bounds getBoundingBox() const;

protected:
    virtual void setDescriptions() {
        setDescription("Render Points at the position selected by the SlicePointRenderer2D. Use in combination with the SliceViewer.");
    }

    void process();

    virtual void initialize();
    virtual void deinitialize();

    void updateGeometry();
    void renderDiskHelper(tgt::Shader* prog, const tgt::ivec4 viewport, const tgt::mat4& pvMat, const tgt::vec3& pointPos, const tgt::vec4 col);
    void renderLinesHelper(const tgt::vec3& geomLlf, const tgt::vec3& geomUrb, const tgt::vec3& pointPos, const tgt::vec4& col);
private:
    //port
    VolumePort inport_;
    //general
    BoolProperty enable_;
    IntProperty pointRadius_;
    BoolProperty renderLines_;

    //points
    BoolProperty renderPoint0_;
    ColorProperty pointColor0_;
    IntVec3Property pointPos0_;
    BoolProperty renderPoint1_;
    ColorProperty pointColor1_;
    IntVec3Property pointPos1_;
    BoolProperty renderPoint2_;
    ColorProperty pointColor2_;
    IntVec3Property pointPos2_;
    BoolProperty renderPoint3_;
    ColorProperty pointColor3_;
    IntVec3Property pointPos3_;
    //shader
    ShaderProperty shaderProp_;
    //member
    TriangleMeshGeometryUInt16IndexedSimple* currentGeometry_;
};

}

#endif

