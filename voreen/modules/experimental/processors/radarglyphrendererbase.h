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

#ifndef VRN_RADARGLYPHRENDERERBASE_H
#define VRN_RADARGLYPHRENDERERBASE_H

//super class
#include "voreen/core/processors/renderprocessor.h"
//ports
#include "voreen/core/ports/volumeport.h"
//properties
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/vectorproperty.h"

namespace voreen {

class RadarGlyphRendererBase : public RenderProcessor
{
public:
    RadarGlyphRendererBase(const std::string& radarGlyphVertexFileName = "rg_radarglyph.vert",
                           const std::string& radarGlyphGeometryFileName = "rg_radarglyph.geom",
                           const std::string& radarGlyphFragmentFileName = "rg_radarglyph.frag",
                           const std::string& billBoardVertexFileName = "rg_billboard.vert",
                           const std::string& billBoardGeometryFileName = "rg_billboard.geom",
                           const std::string& billBoardFragmentFileName = "rg_billboard.frag");
    ~RadarGlyphRendererBase();

    virtual std::string getCategory() const { return "Radar Glyph"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }
    virtual bool usesExpensiveComputation() const { return true; }
    virtual void initialize();
    virtual void deinitialize();
protected:
    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0); /// adds geometry shader to header
    virtual void beforeProcess();
    virtual void compile();
    virtual void setDescriptions() {
        setDescription("Renders RadarGlyphs of a volume collection");
    }

    //ports
    VolumeListPort vcInport_;
    RenderPort imgOutport_;
    //properties
        //representation
    ShaderProperty rgShaderProp_;
    ShaderProperty bbShaderProp_;
    //ShaderProperty pdShaderProp_;
    OptionProperty<int> renderModeProp_;
    TransFunc1DKeysProperty tfProp_;
    FloatProperty lineWidthProp_;
        //performance
    IntProperty volumeStepProp_;
    IntProperty interpolationStepProp_;
    FloatProperty maxLengthProp_;

        //active index
    BoolProperty showIndexProp_;
    IntProperty activeIndexProp_;
    FloatProperty sphereSizeProp_;
    ColorProperty sphereColorProp_;

    //VBO
    void updateGlyphVBO();
    virtual void fillGlyphVBO() = 0;
    //void updateIndexVBO();
    //virtual void fillIndexVBO() = 0;

    //render
    virtual void renderRadarGlyphs();
    virtual void renderActiveIndex();
    //virtual void renderPointData();

    //onChange
    virtual void inportOnChange();
    virtual void renderModeOnChange() = 0;
    virtual void volumeStepOnChange();
    virtual void interpolationStepOnChange();
    virtual void showIndexOnChange();
    virtual void activeIndexOnChange();

    /*struct PointData {
        PointData(tgt::vec3 pos, float color) : pos_(pos), color_(color) {}
        tgt::vec3 pos_;
        float color_;
    };
    std::vector<PointData> pointData_;*/

    GLuint glyphVBO_;
    //GLuint indexVBO_;

    GLint verticesPerGlyph_;
    GLint numberOfGlyphs_;

    tgt::Shader* rgProgram_;        /// the radar shader
    tgt::Shader* bbProgram_;        /// the bill borad shader
    //tgt::Shader* pdProgram_;        /// the point data shader

    bool calcGlyphsNew_;//, calcIndexNew_;
    bool gpuGsSupport_;             /// true, if gs are supported
    std::string gpuGsErrorString_;  /// the error string
};

}   // namespace

#endif  // VRN_RADARGLYPHRENDERERBASE_H
