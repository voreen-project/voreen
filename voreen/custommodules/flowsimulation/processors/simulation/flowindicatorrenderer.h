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

#ifndef VRN_FLOWINDICATORRENDERER_H
#define VRN_FLOWINDICATORRENDERER_H

#include "voreen/core/processors/geometryrendererbase.h"

#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include "../../ports/flowparametrizationport.h"

namespace voreen {

/**
 * Used to render in and out flow locations.
 */
class VRN_CORE_API FlowIndicatorRenderer : public GeometryRendererBase {
public:
    FlowIndicatorRenderer();

    virtual Processor* create() const { return new FlowIndicatorRenderer(); }
    virtual std::string getCategory() const  { return "Geometry"; }
    virtual std::string getClassName() const { return "FlowIndicatorRenderer"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL; }

    virtual void render();

    virtual tgt::Bounds getBoundingBox() const;

protected:
    virtual void setDescriptions() {
        setDescription("Used to render in and out flow locations.");
    }

    void process();

    virtual void initialize();
    virtual void deinitialize();

private:

    FlowParametrizationPort inport_;

    BoolProperty enable_;
    ColorProperty flowGeneratorColor_;
    ColorProperty pressureBoundaryColor_;
    ColorProperty measureFluxColor_;

    std::unique_ptr<GlMeshGeometryUInt16Simple> diskGeometry_;
    std::unique_ptr<GlMeshGeometryUInt16Simple> coneGeometry_;
};

}

#endif

