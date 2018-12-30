/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "imagesequenceloopfinalizer.h"

#include "voreen/core/datastructures/imagesequence.h"
#include "tgt/gpucapabilities.h"

namespace voreen {

const std::string ImageSequenceLoopFinalizer::loggerCat_("voreen.ImageSequenceLoopFinalizer");

ImageSequenceLoopFinalizer::ImageSequenceLoopFinalizer()
    : RenderProcessor(),
      inport_(Port::INPORT, "rendering.in", "Image Input"),
      outport_(Port::OUTPORT, "imagesequence.out", "ImageSequence Output"),
      loopOutport_(Port::OUTPORT, "loop.out", "Loop Outport"),
      textureFormat_("textureFormat", "Texture Format"),
      textureDataType_("textureDataType", "Texture Data Type"),
      keepPixelData_("keepPixelData", "Keep Pixel Data", true),
      imageSequence_(0)
{
    loopOutport_.setLoopPort(true);

    addPort(inport_);
    addPort(outport_);
    addPort(loopOutport_);

    textureFormat_.addOption("same-as-inport",      "Same as Inport",   0);
    textureFormat_.addOption("luminance",           "Luminance",        GL_RED);
    textureFormat_.addOption("rgb",                 "RGB",              GL_RGB);
    textureFormat_.addOption("rgba",                "RGBA",             GL_RGBA);

    textureDataType_.addOption("same-as-inport",    "Same as Inport",     (GLenum)0);
    textureDataType_.addOption("ubyte",             "GLubyte (8 bit)",    GL_UNSIGNED_BYTE);
    textureDataType_.addOption("ushort",            "GLushort (16 bit)",  GL_UNSIGNED_SHORT);
    textureDataType_.addOption("float",             "GLfloat",            GL_FLOAT);

    addProperty(textureFormat_);
    addProperty(textureDataType_);
    addProperty(keepPixelData_);
}

Processor* ImageSequenceLoopFinalizer::create() const {
    return new ImageSequenceLoopFinalizer();
}

bool ImageSequenceLoopFinalizer::isReady() const {
    return (outport_.isReady());
}

void ImageSequenceLoopFinalizer::process() {
    // clear current sequence, if no valid data at inport
    if (!inport_.isReady()) {
        outport_.setData(0);
        return;
    }

    // first iteration: clear previous image sequence
    if (loopOutport_.getLoopIteration() == 0) {
        delete imageSequence_;
        imageSequence_ = new ImageSequence();
    }

    // add current rendering to image sequence
    tgt::Texture* tex = createTextureFromInport();
    imageSequence_->add(tex);

    // last iteration: image sequence is complete => write it to outport
    if (loopOutport_.getLoopIteration() == loopOutport_.getNumLoopIterations()-1) {
        outport_.setData(imageSequence_, false);
    } else {
        loopOutport_.invalidatePort();
    }

    LGL_ERROR;

}

tgt::Texture* ImageSequenceLoopFinalizer::createTextureFromInport() {

    tgtAssert(inport_.getColorTexture(), "No color texture");

    GLint format;
    if (textureFormat_.isSelected("same-as-inport"))
        format = inport_.getColorTexture()->getGLFormat();
    else
        format = textureFormat_.getValue();

    GLenum dataType;
    if (textureDataType_.isSelected("same-as-inport"))
        dataType = inport_.getColorTexture()->getGLDataType();
    else
        dataType = textureDataType_.getValue();

    // read pixel data of render target and construct texture from it
    GLubyte* pixels = inport_.getColorTexture()->downloadTextureToBuffer(format, dataType);
    tgt::ivec3 dim(inport_.getSize().x, inport_.getSize().y, 1);
    tgt::Texture* texture = new tgt::Texture(dim,format,format,dataType,tgt::Texture::LINEAR,tgt::Texture::REPEAT, pixels,true);
    texture->uploadTexture();

    // swizzle so the red texture appears as luminance to the shader
    GLenum swizzleMask[] = {GL_RED, GL_RED, GL_RED, GL_ONE};
    texture->setSwizzle((GLint*)swizzleMask);

    // remove texture data from cpu memory
    if (!keepPixelData_.get()) {
        texture->setCpuTextureData(0,false);
    }

    LGL_ERROR;
    return texture;
}

} // voreen namespace
