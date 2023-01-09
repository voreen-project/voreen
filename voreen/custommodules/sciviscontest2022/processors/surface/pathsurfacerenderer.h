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

#ifndef VRN_PATHSURFACERENDERER_H
#define VRN_PATHSURFACERENDERER_H
#include "modules/base/processors/geometry/geometryrenderer.h"
#include "voreen/core/properties/numeric/intervalproperty.h"

#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"

namespace voreen {

    class VRN_CORE_API PathSurfaceRenderer: public GeometryRenderer {

    public:
        Processor* create() const override { return new PathSurfaceRenderer(); }

        virtual std::string getCategory() const override { return "Path surface Processing"; }
        virtual std::string getClassName() const override { return "PathSurfaceRenderer"; }
        virtual Processor::CodeState getCodeState() const override { return CODE_STATE_EXPERIMENTAL; }
    protected:
        virtual void setDescriptions() override {
            setDescription(
                    "This processor is used to render path surfaces.");
        }

    public:
        PathSurfaceRenderer();

        virtual void render() override;
        virtual void renderTransparent() override;
        virtual void render(ShaderProperty& shaderProp) override;

        virtual bool isReady() const override;
    protected:
        virtual void adjustPropertiesToInput() override;
        void setupPathsurfaceShader();

        std::unique_ptr<voreen::VolumeBase> tfVolume_;

        BoolProperty enableTimelines_;
        BoolProperty enableTransparency_;
        BoolProperty useGradient_;
        FloatProperty texCoordYScale_;
        FloatProperty timelineWidth_;
        IntIntervalProperty timeInterval_;
        TransFunc1DKeysProperty gradient_;
        ShaderProperty pathSurfaceShaderProperty_;
    };

}

#endif //VRN_PATHSURFACERENDERER_H
