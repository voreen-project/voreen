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

#ifndef VRN_QUADRICRENDERER_H
#define VRN_QUADRICRENDERER_H

#include "voreen/core/processors/geometryrendererbase.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include "tgt/shadermanager.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/colorproperty.h"

namespace voreen {

/**
 * Renders a GLU quadric.
 */
class VRN_CORE_API QuadricRenderer : public GeometryRendererBase {
public:
    QuadricRenderer();
    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "QuadricRenderer";  }
    virtual std::string getCategory() const     { return "Geometry";         }
    virtual CodeState getCodeState() const      { return CODE_STATE_STABLE; }

    virtual tgt::Bounds getBoundingBox() const;

    virtual void render();

protected:
    virtual void setDescriptions() {
        setDescription("Renders a quadric geometry using a generated mesh. Currently supported: sphere and cylinder.");
    }

    virtual void initialize();
    virtual void deinitialize();

    void invalidateGeometry();

private:
    void adjustPropertyVisibilities();

    GlMeshGeometryUInt16Normal mesh_;
    tgt::Shader* shader_; //< mesh rendering shader with illumination

    BoolProperty enabled_;
    StringOptionProperty quadricType_;

    BoolProperty closeCylinder_;
    FloatVec3Property position_;
    FloatVec3Property start_;
    FloatVec3Property end_;
    FloatProperty radius_;
    ColorProperty color_;

    BoolProperty applyLighting_;
    FloatVec3Property lightPosition_;
    ColorProperty lightAmbient_;
    ColorProperty lightDiffuse_;
    ColorProperty lightSpecular_;
    FloatProperty materialShininess_;

};

}

#endif // VRN_QUADRICRENDERER_H

