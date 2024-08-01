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

#include "radarglyphrendererbase.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/ports/genericport.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "tgt/gpucapabilities.h"
#include "tgt/glmath.h"

namespace voreen {

RadarGlyphRendererBase::RadarGlyphRendererBase(const std::string& radarGlyphVertexFileName,
                                               const std::string& radarGlyphGeometryFileName,
                                               const std::string& radarGlyphFragmentFileName,
                                               const std::string& billBoardVertexFileName,
                                               const std::string& billBoardGeometryFileName,
                                               const std::string& billBoardFragmentFileName)
    : RenderProcessor()
    //ports
    , vcInport_(Port::INPORT, "inport.volumecollection")
    , imgOutport_(Port::OUTPORT, "outport.image","Outport",false,Processor::INVALID_RESULT,RenderPort::RENDERSIZE_RECEIVER)
    //representation
    , rgShaderProp_("radarglyph.prg", "RG Shader:", radarGlyphFragmentFileName,radarGlyphVertexFileName,radarGlyphGeometryFileName)
    , bbShaderProp_("billboard.prg", "BB Shader:", billBoardFragmentFileName,billBoardVertexFileName,billBoardGeometryFileName)
    //, pdShaderProp_("pointdata.prg", "PD Shader:", "rg_pointdata.frag", "rg_pointdata.vert")
    , renderModeProp_("renderModeProp", "Render Mode:")
    , tfProp_("tfProp","Transfer Function:")
    , lineWidthProp_("lineWidthProp","Line Width:",1.f,0.01f,10.f)
    //performance
    , volumeStepProp_("volumeStepProp", "Volume Steps:",1,1,100)
    , interpolationStepProp_("interpolateStepProp","Shader Steps:",0,0,3)
    , maxLengthProp_("maxLengthProp","Max Length:",0.f,0.f,1000.f)
    //active index
    , showIndexProp_("showIndexProp","Show Index:",true)
    , activeIndexProp_("activeIndex","Active Index:",0,0,1000)
    , sphereSizeProp_("sphereSizeProp","Size:",0.5f,0.f,10.f)
    , sphereColorProp_("sphereColorProp","Color:",tgt::vec4(1.f))
    //members
    //, pointData_()
    , glyphVBO_(0), /*indexVBO_(0),*/ verticesPerGlyph_(0), numberOfGlyphs_(0)
    , rgProgram_(0), bbProgram_(0)/*, pdProgram_(0)*/
    , calcGlyphsNew_(true)//, calcIndexNew_(true)
    , gpuGsSupport_(true)
    , gpuGsErrorString_("")

{
    //ports
    addPort(vcInport_);
        //vcInport_.addCondition(new PortConditionVolumeChannelCount(3));
    addPort(imgOutport_);

    addProperty(rgShaderProp_);
        rgShaderProp_.setGroupID("Representation");
    addProperty(bbShaderProp_);
        bbShaderProp_.setGroupID("Representation");
    /*addProperty(pdShaderProp_);
        pdShaderProp_.setGroupID("Representation");*/
    addProperty(renderModeProp_);
        renderModeProp_.setGroupID("Representation");
        ON_PROPERTY_CHANGE(renderModeProp_,RadarGlyphRendererBase,renderModeOnChange);
    addProperty(tfProp_);
        tfProp_.setGroupID("Representation");
    addProperty(lineWidthProp_);
        lineWidthProp_.setGroupID("Representation");
    setPropertyGroupGuiName("Representation","Representation");

    addProperty(volumeStepProp_);
        volumeStepProp_.setGroupID("Performance");
        ON_PROPERTY_CHANGE(volumeStepProp_,RadarGlyphRendererBase,volumeStepOnChange);
    addProperty(interpolationStepProp_);
        interpolationStepProp_.setGroupID("Performance");
        ON_PROPERTY_CHANGE(interpolationStepProp_,RadarGlyphRendererBase,interpolationStepOnChange);
    addProperty(maxLengthProp_);
        maxLengthProp_.setGroupID("Performance");
    setPropertyGroupGuiName("Performance","Performance");

    addProperty(showIndexProp_);
        showIndexProp_.setGroupID("Active Index");
        ON_PROPERTY_CHANGE(showIndexProp_,RadarGlyphRendererBase,showIndexOnChange);
    addProperty(activeIndexProp_);
        activeIndexProp_.setGroupID("Active Index");
        ON_PROPERTY_CHANGE(activeIndexProp_,RadarGlyphRendererBase,activeIndexOnChange);
    addProperty(sphereSizeProp_);
        sphereSizeProp_.setGroupID("Active Index");
    addProperty(sphereColorProp_);
        sphereColorProp_.setGroupID("Active Index");
    setPropertyGroupGuiName("Active Index","Active Index");
}

RadarGlyphRendererBase::~RadarGlyphRendererBase() {
}

std::string RadarGlyphRendererBase::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string header = RenderProcessor::generateHeader(version);
    std::stringstream strstr;
    strstr << header;
    strstr << "#define ARRAY_SIZE_ " << static_cast<int>(std::pow(2.f,interpolationStepProp_.get()))+1 << "\n";
    strstr << "#extension GL_EXT_geometry_shader4 : enable" <<"\n";
    return strstr.str();
}

void RadarGlyphRendererBase::initialize() {
    RenderProcessor::initialize();
    //test if GSs are supported
    std::stringstream strstr;
    strstr.clear();
    strstr.str("");
    if(GpuCaps.getShaderModel() <= tgt::GpuCapabilities::SHADER_MODEL_4){
        gpuGsSupport_ = false;
        strstr << "No GS support. GPU supports ShaderModel "<< GpuCaps.getShaderModel() << ".0, but 4.0 is needed.";
    } else {
        if(GpuCaps.getMaxGeometryShaderVertices() < 10){
            gpuGsSupport_ = false;
            strstr << "The GS supports just " << GpuCaps.getMaxGeometryShaderVertices() << " Vertices, but 10 are needed";
        }
    }

    if(gpuGsSupport_)
        compile();
    else
        gpuGsErrorString_ = strstr.str();

    glGenBuffers(1, &glyphVBO_);
    //glGenBuffers(1, &indexVBO_);
    setProgress(0.f);
}

void RadarGlyphRendererBase::deinitialize() {
    RenderProcessor::deinitialize();
    glDeleteBuffers(1,&glyphVBO_);
    //glDeleteBuffers(1,&indexVBO_);
}

void RadarGlyphRendererBase::compile() {
    rgShaderProp_.setHeader(generateHeader());
    rgShaderProp_.rebuild();
    if(rgShaderProp_.hasValidShader()) {
        rgProgram_ = rgShaderProp_.getShader();
        rgProgram_->deactivate();
    } else {
        LERROR(gpuGsErrorString_);
        rgProgram_ = 0;
    }
    bbShaderProp_.setHeader(generateHeader());
    bbShaderProp_.rebuild();
    if(bbShaderProp_.hasValidShader()) {
        bbProgram_ = bbShaderProp_.getShader();
        bbProgram_->deactivate();
    } else {
        LERROR(gpuGsErrorString_);
        bbProgram_ = 0;
    }
    /*pdShaderProp_.setHeader(generateHeader());
    pdShaderProp_.rebuild();
    if(pdShaderProp_.hasValidShader()) {
        pdProgram_ = pdShaderProp_.getShader();
        pdProgram_->deactivate();
    } else {
        LERROR(gpuGsErrorString_);
        pdProgram_ = 0;
    }*/
}

void RadarGlyphRendererBase::beforeProcess() {
    RenderProcessor::beforeProcess();

    // compile program if needed
    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM && gpuGsSupport_) {
        PROFILING_BLOCK("compile");
        compile();
    }
    LGL_ERROR;
}

//-----------------------------------------------------------------------
//      VBO
//-----------------------------------------------------------------------
void RadarGlyphRendererBase::updateGlyphVBO() {
    glDeleteBuffers(1,&glyphVBO_);
    glGenBuffers(1, &glyphVBO_);
    glBindBuffer(GL_ARRAY_BUFFER, glyphVBO_); //dir.x,dir.y,dir.z,center.x,center.y,center.z,step
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*7*numberOfGlyphs_*verticesPerGlyph_, NULL, GL_STATIC_DRAW);
    fillGlyphVBO();
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    //updateIndex
    /*if(showIndexProp_.get()){
        updateIndexVBO();
    }*/
    calcGlyphsNew_ = false;
}

/*void RadarGlyphRendererBase::updateIndexVBO() {
    glDeleteBuffers(1,&indexVBO_);
    glGenBuffers(1, &indexVBO_);
    glBindBuffer(GL_ARRAY_BUFFER, indexVBO_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*3*numberOfGlyphs_, NULL, GL_STATIC_DRAW);
    fillIndexVBO();
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    calcIndexNew_ = false;
}*/

//-----------------------------------------------------------------------
//      render
//-----------------------------------------------------------------------
void RadarGlyphRendererBase::renderRadarGlyphs() {
    glBindBuffer(GL_ARRAY_BUFFER, glyphVBO_);
    GLint locCenter = rgProgram_->getAttributeLocation("center_");
    GLint locStep  = rgProgram_->getAttributeLocation("timestep_");

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableVertexAttribArray(locCenter);
    glEnableVertexAttribArray(locStep);

    glVertexPointer(3, GL_FLOAT, sizeof(GLfloat)*7 , NULL);
    glVertexAttribPointer(locCenter, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat)*7, (float*)(sizeof(GLfloat)*3));
    glVertexAttribPointer(locStep, 1, GL_FLOAT, GL_FALSE, sizeof(GLfloat)*7, (float*)(sizeof(GLfloat)*6));

    glLineWidth(lineWidthProp_.get());

    for(int i = 0; i < numberOfGlyphs_; i++)
        glDrawArrays(GL_LINE_STRIP, i*verticesPerGlyph_, verticesPerGlyph_);

    glLineWidth(1.f);

    glDisableVertexAttribArray(locStep);
    glDisableVertexAttribArray(locCenter);
    glDisableClientState(GL_VERTEX_ARRAY);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}


void RadarGlyphRendererBase::renderActiveIndex() {
    //glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glColor4fv(sphereColorProp_.get().elem);
    GLint locCenter = rgProgram_->getAttributeLocation("center_");
    int t = activeIndexProp_.get()+1;

    glBindBuffer(GL_ARRAY_BUFFER, glyphVBO_);

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableVertexAttribArray(locCenter);

    glVertexPointer(3, GL_FLOAT, sizeof(GLfloat)*7*verticesPerGlyph_ , (float*)(t*sizeof(GLfloat)*7));
    glVertexAttribPointer(locCenter, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat)*7*verticesPerGlyph_, (float*)(sizeof(GLfloat)*(t*7+3)));

    glDrawArrays(GL_POINTS, 0, numberOfGlyphs_); // TODO: 2D case are less glyphs

    glDisableVertexAttribArray(locCenter);
    glDisableClientState(GL_VERTEX_ARRAY);

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glColor4fv(tgt::vec4(1.f).elem);
    glBlendFunc(GL_ONE, GL_ZERO);
    glDisable(GL_BLEND);
    //glEnable(GL_DEPTH_TEST);
}

/*void RadarGlyphRendererBase::renderPointData() {
     GLint locColor = pdProgram_->getAttributeLocation("color_");

    glBegin(GL_POINTS);
    for(std::vector<PointData>::const_iterator it = pointData_.begin(); it != pointData_.end(); it++) {
        glVertexAttrib1f(locColor,it->color_);
        glVertex3fv(it->pos_.elem);
    }
    glEnd();
}*/

//-----------------------------------------------------------------------
//      onChange
//-----------------------------------------------------------------------
void RadarGlyphRendererBase::inportOnChange() {
    //set maximum
    float maximum = 0.f;
    for(size_t i = 0; i < vcInport_.getData()->size(); i++){
        float tmp = vcInport_.getData()->at(i)->getRepresentation<VolumeRAM_3xFloat>()->maxNormalizedMagnitude();
        if(tmp > maximum)
            maximum = tmp;
    }
    maxLengthProp_.set(maximum);

    verticesPerGlyph_ = static_cast<int>(vcInport_.getData()->size())/volumeStepProp_.get()+1;
    activeIndexProp_.setMaxValue(verticesPerGlyph_-2);
    tgt:: svec3 dim = vcInport_.getData()->at(0)->getDimensions();
    numberOfGlyphs_ = static_cast<GLint>(dim.x*dim.y*dim.z);
    calcGlyphsNew_ = true;
}

void RadarGlyphRendererBase::volumeStepOnChange() {
    activeIndexProp_.setStepping(volumeStepProp_.get());
    if(int tmp= (activeIndexProp_.get() % volumeStepProp_.get())){
        activeIndexProp_.set(activeIndexProp_.get() - tmp);
    }
    calcGlyphsNew_ = true;
    invalidate();
}

void RadarGlyphRendererBase::interpolationStepOnChange() {
    if(isInitialized())
        compile();
    invalidate();
}

void RadarGlyphRendererBase::showIndexOnChange() {
    activeIndexProp_.setReadOnlyFlag(!showIndexProp_.get());
    sphereSizeProp_.setReadOnlyFlag(!showIndexProp_.get());
    sphereColorProp_.setReadOnlyFlag(!showIndexProp_.get());
    invalidate();
}

void RadarGlyphRendererBase::activeIndexOnChange() {
    invalidate();
}

}   // namespace
