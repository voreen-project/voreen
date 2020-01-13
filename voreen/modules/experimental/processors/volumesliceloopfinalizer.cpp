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

#include "volumesliceloopfinalizer.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/imagesequence.h"
#include "tgt/gpucapabilities.h"

namespace voreen {

using tgt::Texture;
using tgt::ivec3;

const std::string VolumeSliceLoopFinalizer::loggerCat_("voreen.experimental.VolumeSliceLoopFinalizer");

VolumeSliceLoopFinalizer::VolumeSliceLoopFinalizer()
    : RenderProcessor(),
      inport_(Port::INPORT, "rendering.in"),
      outport_(Port::OUTPORT, "imagesequence.out"),
      loopOutport_(Port::OUTPORT,"loop.out"),
      sliceAlignment_("sliceAlignmentProp", "Slice Alignment"),
      textureFormat_("textureFormat", "Texture Format"),
      textureDataType_("textureDataType", "Texture Data Type"),
      volume_(0)
{
    loopOutport_.setLoopPort(true);

    addPort(inport_);
    addPort(outport_);
    addPort(loopOutport_);

    sliceAlignment_.addOption("xy-plane", "XY-Plane (axial)", XY_PLANE);
    sliceAlignment_.addOption("xz-plane", "XZ-Plane (coronal)", XZ_PLANE);
    sliceAlignment_.addOption("yz-plane", "YZ-Plane (sagittal)", YZ_PLANE);
    addProperty(sliceAlignment_);

    textureFormat_.addOption("same-as-inport",      "Same as Inport",   0);
    textureFormat_.addOption("luminance",           "Luminance",        GL_LUMINANCE);
    textureFormat_.addOption("rgb",                 "RGB",              GL_RGB);
    textureFormat_.addOption("rgba",                "RGBA",             GL_RGBA);

    textureDataType_.addOption("same-as-inport",    "Same as Inport",     (GLenum)0);
    textureDataType_.addOption("ubyte",             "GLubyte (8 bit)",    GL_UNSIGNED_BYTE);
    textureDataType_.addOption("ushort",            "GLushort (16 bit)",  GL_UNSIGNED_SHORT);
    textureDataType_.addOption("float",             "GLfloat",            GL_FLOAT);

    addProperty(textureFormat_);
    addProperty(textureDataType_);
}

VolumeSliceLoopFinalizer::~VolumeSliceLoopFinalizer() {
    delete volume_;
}

Processor* VolumeSliceLoopFinalizer::create() const {
    return new VolumeSliceLoopFinalizer();
}

bool VolumeSliceLoopFinalizer::isReady() const {
    //return (outport_.isReady());
    return (inport_.isReady() && outport_.isReady());
}

void VolumeSliceLoopFinalizer::initialize() {
    RenderProcessor::initialize();
    //volume_ = new Volume();
}

void VolumeSliceLoopFinalizer::process() {

    // clear current sequence, if no valid data at inport
    if (!inport_.isReady()) {
        delete volume_;
        volume_ = 0;
        outport_.setData(0);
        return;
    }

    // first iteration: clear previous image sequence
    if (loopOutport_.getLoopIteration() == 0) {
        VolumeBase* vhb = const_cast<VolumeBase*>(outport_.getData());
        if(vhb) {
            Volume* temp = dynamic_cast<Volume*>(vhb);
            temp->deleteAllRepresentations();
            delete temp;
            outport_.setData(0);
        }
        //delete volume_; //the volume is one of the representations we just deleted
        volume_ = 0;

        ivec3 dims; //(loopOutport_.getNumLoopIterations()-1, inport_.getSize().x, inport_.getSize().y);

        switch(sliceAlignment_.getValue()) {
            case YZ_PLANE: dims = ivec3(loopOutport_.getNumLoopIterations(), inport_.getSize().x, inport_.getSize().y);
                           break;
            case XZ_PLANE: dims = ivec3(inport_.getSize().x, loopOutport_.getNumLoopIterations(), inport_.getSize().y);
                           break;
            case XY_PLANE: dims = ivec3(inport_.getSize().x, inport_.getSize().y, loopOutport_.getNumLoopIterations());
                           break;
            default: tgtAssert(false, "should not get here!");
        }

        volume_ = new VolumeRAM_UInt8(dims);
    }

    tgt::ivec3 dims = volume_->getDimensions();
    // add current rendering to image sequence
    switch(sliceAlignment_.getValue()) {
        case YZ_PLANE: {
                           int x = loopOutport_.getLoopIteration();
                           if((x < 0) || (x >= dims.x)) {
                                LERROR("Out of bounds! " << x);
                                return;
                           }

                           Texture* tex = inport_.getColorTexture();
                           tex->downloadTexture();

                           for(int y=0; y<dims.y; y++) {
                               for(int z=0; z<dims.z; z++) {
                                   float value = tex->texelAsFloat(y, z).r;
                                   volume_->setVoxelNormalized(value, ivec3(x,y,z));
                               }
                           }
                       }
                       break;
        case XZ_PLANE: {
                           int y = loopOutport_.getLoopIteration();
                           if((y < 0) || (y >= dims.y)) {
                                LERROR("Out of bounds! " << y);
                                return;
                           }

                           Texture* tex = inport_.getColorTexture();
                           tex->downloadTexture();

                           for(int x=0; x<dims.x; x++) {
                               for(int z=0; z<dims.z; z++) {
                                   float value = tex->texelAsFloat(x, z).r;
                                   volume_->setVoxelNormalized(value, ivec3(x,y,z));
                               }
                           }
                       }
                       break;
        case XY_PLANE: {
                           int z = loopOutport_.getLoopIteration();
                           if((z < 0) || (z >= dims.z)) {
                                LERROR("Out of bounds! " << z);
                                return;
                           }

                           Texture* tex = inport_.getColorTexture();
                           tex->downloadTexture();

                           for(int y=0; y<dims.y; y++) {
                               for(int x=0; x<dims.x; x++) {
                                   float value = tex->texelAsFloat(x, y).r;
                                   volume_->setVoxelNormalized(value, ivec3(x,y,z));
                               }
                           }
                       }
                       break;
        default: tgtAssert(false, "should not get here!");
    }

    // last iteration: image sequence is complete => write it to outport
    if (loopOutport_.getLoopIteration() == loopOutport_.getNumLoopIterations()-1) {
        //volume_->invalidate();//TODO
        Volume* vh = new Volume(volume_, tgt::vec3(1.0f), tgt::vec3(0.0f));//FIXME: metadata
        oldVolumePosition(vh); //FIXME
        outport_.setData(vh, false);
    } else {
        loopOutport_.invalidatePort();
    }

    LGL_ERROR;

}

} // voreen namespace
