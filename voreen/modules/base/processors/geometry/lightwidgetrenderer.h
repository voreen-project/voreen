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

#ifndef VRN_LIGHTWIDGETRENDERER_H
#define VRN_LIGHTWIDGETRENDERER_H

#include "voreen/core/processors/geometryrendererbase.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include "tgt/immediatemode/immediatemode.h"

namespace voreen {

/**
 * Draws a yellow sphere that indicates the position of the light source and can be moved.
 * Use linking to connect this widget to a light position of another processor.
 */
class VRN_CORE_API LightWidgetRenderer : public GeometryRendererBase {
public:
    LightWidgetRenderer();
    ~LightWidgetRenderer();
    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "LightWidgetRenderer"; }
    virtual std::string getCategory() const     { return "Geometry"; }
    virtual CodeState getCodeState() const      { return CODE_STATE_STABLE; }

    /// Renders the light widget
    virtual void render();
    virtual void renderPicking();
    virtual void setIDManager(IDManager* idm);
    virtual tgt::Bounds getBoundingBox() const;

private:
    virtual void setDescriptions() {
        setDescription("Draws a yellow sphere that indicates the position of the lightsource and can be moved. " \
                       "Use linking to connect this widget to a light source of another processor.");
    }

    void moveSphere(tgt::MouseEvent* e);

    BoolProperty showLightWidget_;
    EventProperty<LightWidgetRenderer>* moveSphereProp_;

    bool isClicked_;

    FloatVec4Property lightPosition_;
    FloatProperty sphereRadius_;
    ColorProperty sphereColor_;
    tgt::vec4 lightPositionAbs_;
    tgt::ivec2 startCoord_;

    GlMeshGeometryUInt16Normal mesh_;

    tgt::ImmediateMode::LightSource lightSource_;
    tgt::ImmediateMode::Material material_;
};

}

#endif

