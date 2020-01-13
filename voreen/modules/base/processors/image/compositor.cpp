/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "compositor.h"

#include "tgt/textureunit.h"

using tgt::TextureUnit;

namespace voreen {

Compositor::Compositor()
    : ImageProcessor("image/compositor")
    , compositingMode_("blendMode", "Blend mode", Processor::INVALID_PROGRAM)
    , weightingFactor_("weightingFactor", "Weighting factor", 0.5f, 0.0f, 1.0f)
    , weightFirst_("weightFirst", "Weight (First Image)", 1.0f, 0.f, 1.0f)
    , weightSecond_("weightSecond", "Weight (Second Image)", 1.0f, 0.f, 1.0f)
    , addDepth_("addDepth", "Add Depth", false)
    , pipeThroughIfSecondNotReady_("pipethrough", "Pipe first image through if second is not ready", false)
    , inport0_(Port::INPORT, "image.inport0", "Image0 Input")
    , inport1_(Port::INPORT, "image.inport1", "Image1 Input")
    , outport_(Port::OUTPORT, "image.outport", "Image Output")
{

    compositingMode_.addOption("depth-test",                "Depth Test",                   "MODE_DEPTH_TEST");
    compositingMode_.addOption("alpha-compositing",         "Alpha Compositing",            "MODE_ALPHA_COMPOSITING");
    compositingMode_.addOption("alpha-blending",            "Alpha Blending",               "MODE_ALPHA_BLENDING");
    compositingMode_.addOption("alpha-blending-b-over-a",   "Alpha Blending (B over A)",    "MODE_ALPHA_BLENDING_B_OVER_A");
    compositingMode_.addOption("weighted-average",          "Weighted Average",             "MODE_WEIGHTED_AVERAGE");
    compositingMode_.addOption("take-first",                "Take First",                   "MODE_TAKE_FIRST");
    compositingMode_.addOption("take-second",               "Take Second",                  "MODE_TAKE_SECOND");
    compositingMode_.addOption("take-second-if-ready",      "Take Second if ready",         "MODE_TAKE_SECOND_IF_READY");
    compositingMode_.addOption("first-has-priority",        "First Has Priority",           "MODE_FIRST_HAS_PRIORITY");
    compositingMode_.addOption("second-has-priority",       "Second Has Priority",          "MODE_SECOND_HAS_PRIORITY");
    compositingMode_.addOption("maximum-alpha",             "Maximum Alpha",                "MODE_MAXIMUM_ALPHA");
    compositingMode_.addOption("difference",                "Difference",                   "MODE_DIFFERENCE");
    compositingMode_.addOption("difference-alpha-mul",      "Difference after alpha mul",   "MODE_DIFFERENCE_ALPHA_MUL");
    compositingMode_.addOption("add",                       "Add",                          "MODE_ADD");
    compositingMode_.set("alpha-compositing");
    compositingMode_.onChange(MemberFunctionCallback<Compositor>(this, &Compositor::compositingModeChanged));
    addProperty(compositingMode_);

    addProperty(weightingFactor_);
    addProperty(weightFirst_);
    addProperty(weightSecond_);
    addProperty(addDepth_);
    addProperty(pipeThroughIfSecondNotReady_);

    addPort(inport0_);
    addPort(inport1_);
    addPort(outport_);
}

Compositor::~Compositor() {}

void Compositor::initialize() {
    //call parent's initialize method
    ImageProcessor::initialize();
    //Hide properties that should be invisible at beginning
    compositingModeChanged();

    copyShader_ = ShdrMgr.loadSeparate("passthrough.vert", "copyimage.frag", ImageProcessor::generateHeader(), false);
}

void Compositor::deinitialize() {
    ImageProcessor::deinitialize();
    ShdrMgr.dispose(copyShader_);
    copyShader_ = 0;
}

bool Compositor::isReady() const {
    bool ready = false;
    if((compositingMode_.get() == "depth-test") || (compositingMode_.get() == "alpha-compositing") ||
       (compositingMode_.get() == "alpha-blending") || (compositingMode_.get() == "alpha-blending-b-over-a") ||
       (compositingMode_.get() == "weighted-average") || (compositingMode_.get() == "maximum-alpha") ||
       (compositingMode_.get() == "difference") || (compositingMode_.get() == "add") ||
       (compositingMode_.get() == "first-has-priority") || (compositingMode_.get() == "second-has-priority")||
       (compositingMode_.get() == "difference-alpha-mul")) {
        if(pipeThroughIfSecondNotReady_.get())
            ready = inport0_.isReady() || (inport0_.isReady() && inport1_.isReady());
        else
            ready = (inport0_.isReady() && inport1_.isReady());
    } else if(compositingMode_.get() == "take-first") {
        ready = inport0_.isReady();
    } else if(compositingMode_.get() == "take-second") {
        if(pipeThroughIfSecondNotReady_.get())
            ready = inport0_.isReady() || inport1_.isReady();
        else
            ready = inport1_.isReady();
    } else if(compositingMode_.get() == "take-second-if-ready") {
        ready = inport0_.isReady() || inport1_.isReady();
    }
    return ready;
}

void Compositor::process() {

    if(!inport1_.isReady() && pipeThroughIfSecondNotReady_.get()) {
        outport_.activateTarget();
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // bind input image to tex unit
        TextureUnit imageUnit, imageUnitDepth;
        inport0_.bindTextures(imageUnit.getEnum(), imageUnitDepth.getEnum());

        // copy input image to outport
        copyShader_->activate();
        setGlobalShaderParameters(copyShader_);
        inport0_.setTextureParameters(copyShader_, "texParams_");
        copyShader_->setUniform("useTexcoords", false);
        copyShader_->setUniform("colorTex_", imageUnit.getUnitNumber());
        copyShader_->setUniform("depthTex_", imageUnitDepth.getUnitNumber());
        renderQuad();
        copyShader_->deactivate();
        outport_.deactivateTarget();
        tgt::TextureUnit::setZeroUnit();
        LGL_ERROR;
        return;
    }

    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM)
        compile();

    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDepthFunc(GL_ALWAYS);

    TextureUnit colorUnit, depthUnit;
    TextureUnit colorUnit1, depthUnit1;
    inport0_.bindTextures(colorUnit.getEnum(), depthUnit.getEnum());
    inport1_.bindTextures(colorUnit1.getEnum(), depthUnit1.getEnum());

    // initialize shader
    program_->activate();
    setGlobalShaderParameters(program_);
    if (compositingMode_.get() != "take-second") {
        program_->setUniform("colorTex0_", colorUnit.getUnitNumber());
        program_->setUniform("depthTex0_", depthUnit.getUnitNumber());
        inport0_.setTextureParameters(program_, "textureParameters0_");
    }
    if (compositingMode_.get() != "take-first") {
        program_->setUniform("colorTex1_", colorUnit1.getUnitNumber());
        program_->setUniform("depthTex1_", depthUnit1.getUnitNumber());
        inport1_.setTextureParameters(program_, "textureParameters1_");
    }
    if (compositingMode_.get() == "take-second-if-ready") {
        if(inport1_.isReady()) {
            program_->setUniform("colorTex0_", colorUnit1.getUnitNumber());
            program_->setUniform("depthTex0_", depthUnit1.getUnitNumber());
            inport1_.setTextureParameters(program_, "textureParameters0_");
        }
        else if(inport0_.isReady()){
            program_->setUniform("colorTex0_", colorUnit.getUnitNumber());
            program_->setUniform("depthTex0_", depthUnit.getUnitNumber());
            inport0_.setTextureParameters(program_, "textureParameters0_");
        }
    }
    if (compositingMode_.get() == "weighted-average")
        program_->setUniform("weightingFactor_", weightingFactor_.get());

    if (compositingMode_.get() == "add") {
            program_->setUniform("weightFirst_", weightFirst_.get());
            program_->setUniform("weightSecond_", weightSecond_.get());
            program_->setUniform("addDepth_", addDepth_.get());
        }

    renderQuad();

    glDepthFunc(GL_LESS);

    program_->deactivate();
    outport_.deactivateTarget();
    TextureUnit::setZeroUnit();
    LGL_ERROR;
}

std::string Compositor::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string header = ImageProcessor::generateHeader(version);
    header += "#define " + compositingMode_.getValue() + "\n";
    return header;
}

void Compositor::compile() {
    if (program_)
        program_->setHeaders(generateHeader());
    ImageProcessor::compile();
}

void Compositor::compositingModeChanged() {
    weightingFactor_.setVisibleFlag(compositingMode_.get() == "weighted-average");
    weightFirst_.setVisibleFlag(compositingMode_.get() == "add");
    weightSecond_.setVisibleFlag(compositingMode_.get() == "add");
    addDepth_.setVisibleFlag(compositingMode_.get() == "add");
}

} // voreen namespace
