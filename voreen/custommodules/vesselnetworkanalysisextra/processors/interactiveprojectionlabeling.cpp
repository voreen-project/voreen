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

#include "interactiveprojectionlabeling.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include <iostream>
#include "tgt/init.h"

namespace voreen {

const std::string InteractiveProjectionLabeling::loggerCat_("voreen.vesselnetworkanalysisextra.interactiveprojectionlabeling");

LabelProjection::LabelProjection()
    : LabelProjection(VolumeAtomic<float>(tgt::svec3(2)))
{
}
LabelProjection::LabelProjection(VolumeAtomic<float>&& projection)
    : projection_(std::move(projection))
    , labels_(projection_.getDimensions())
    , projectionTexture_(boost::none)
    , labelTexture_(boost::none)
{
    for(size_t i=0; i < tgt::hmul(labels_.getDimensions()); ++i) {
        labels_.voxel(i) = UNLABELED;
    }
}
void LabelProjection::ensureTexturesPresent() {
    if(!projectionTexture_) {
        projectionTexture_ = tgt::Texture(projection_.getDimensions(), GL_RED, GL_RED, GL_FLOAT, tgt::Texture::NEAREST, tgt::Texture::CLAMP_TO_EDGE, (GLubyte*) projection_.voxel(), false);
        projectionTexture_->uploadTexture();
    }
    if(!labelTexture_) {
        labelTexture_ = tgt::Texture(labels_.getDimensions(), GL_RED, GL_RED, GL_UNSIGNED_BYTE, tgt::Texture::NEAREST, tgt::Texture::CLAMP_TO_EDGE, (GLubyte*) labels_.voxel(), false);
        labelTexture_->uploadTexture();
    }
}
void LabelProjection::bindProjectionTexture() {
    ensureTexturesPresent();
    projectionTexture_->bind();
}
void LabelProjection::bindLabelTexture() {
    ensureTexturesPresent();
    labelTexture_->bind();
}
void LabelProjection::withLabels(std::function<void(VolumeAtomic<uint8_t>&)> fun) {
    fun(labels_);
    ensureTexturesPresent();
    labelTexture_->uploadTexture();
}
void InteractiveProjectionLabeling::drawEvent(tgt::MouseEvent* e, LabelProjection& p) {
    if(e->modifiers() != tgt::Event::CTRL || e->button() != tgt::MouseEvent::MOUSE_BUTTON_LEFT) {
        return;
    }
    auto event_type = e->getEventType();
    tgt::ivec2 coords = e->coord();
    tgt::ivec2 viewport = e->viewport();

    coords.y = viewport.y - coords.y - 1;
    auto norm_coords = tgt::vec2(coords)/tgt::vec2(viewport);
    tgt::svec2 labelcoords = tgt::round(norm_coords * tgt::vec2(p.labels().getDimensions().xy()));

    p.withLabels([&] (VolumeAtomic<uint8_t>& vol) {
            vol.voxel(labelcoords.x, labelcoords.y, 0) = LabelProjection::FOREGROUND;
    });
    if(event_type == tgt::MouseEvent::MOUSERELEASEEVENT) {
        //TODO sync labels => 3D
    }
    invalidate();
    e->accept();
}

void InteractiveProjectionLabeling::onPortEvent(tgt::Event* e, Port* port) {
    if(tgt::MouseEvent* me = dynamic_cast<tgt::MouseEvent*>(e)) {
        if(port == &xyProjectionOutport_) {
            drawEvent(me, xy_);
        }
        else if(port == &xzProjectionOutport_) {
            drawEvent(me, xz_);
        }
        else if(port == &yzProjectionOutport_) {
            drawEvent(me, yz_);
        }
    }
    if(!e->isAccepted()) {
        RenderProcessor::onPortEvent(e, port);
    }
}

InteractiveProjectionLabeling::InteractiveProjectionLabeling()
    : RenderProcessor()
    , inport_(Port::INPORT, "interactiveprojectionlabeling.inport", "Volume Input")
    , labelVolume_(Port::OUTPORT, "interactiveprojectionlabeling.labelVolume", "Labels Output")
    , xyProjectionOutport_(Port::OUTPORT, "interactiveprojectionlabeling.xyProjectionOutport", "Projection (XY)", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , xzProjectionOutport_(Port::OUTPORT, "interactiveprojectionlabeling.xzProjectionOutport", "Projection (XZ)", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , yzProjectionOutport_(Port::OUTPORT, "interactiveprojectionlabeling.yzProjectionOutport", "Projection (YZ)", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , shader_("shader", "Shader", "interactiveprojectionlabeling.frag", "oit_passthrough.vert")
    , xy_()
    , xz_()
    , yz_()
{
    addPort(inport_);
    addPort(labelVolume_);
    addPort(xyProjectionOutport_);
    addPort(xzProjectionOutport_);
    addPort(yzProjectionOutport_);

    addProperty(shader_);
}

void InteractiveProjectionLabeling::renderToPort(RenderPort& port, LabelProjection& p) {
    //TODO not sure why the texture from p does not work.
    //p.projectionTexture_ = tgt::Texture(p.projection().getDimensions(), GL_RED, GL_RED, GL_FLOAT, tgt::Texture::NEAREST, tgt::Texture::CLAMP_TO_EDGE, (GLubyte*) p.projection().voxel(), false);
    //p.projectionTexture_.uploadTexture();
    //p.labelTexture_ = tgt::Texture(p.labels().getDimensions(), GL_RED, GL_RED, GL_UNSIGNED_BYTE, tgt::Texture::NEAREST, tgt::Texture::CLAMP_TO_EDGE, (GLubyte*) p.labels().voxel(), false);
    //p.labelTexture_.uploadTexture();
    //VolumeAtomic<float> proj(p.projection().getDimensions());
    //for(size_t i=0; i < tgt::hmul(proj.getDimensions()); ++i) {
    //    proj.voxel(i) = p.projection().voxel(i);
    //}
    //LabelProjection o(std::move(proj));
    //o.withLabels([&] (VolumeAtomic<uint8_t>& vol) {
    //        for(size_t i=0; i < tgt::hmul(vol.getDimensions()); ++i) {
    //            vol.voxel(i) = p.labels().voxel(i);
    //        }
    //});

    port.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    tgt::TextureUnit projUnit;
    projUnit.activate();
    p.bindProjectionTexture();

    tgt::TextureUnit labelUnit;
    labelUnit.activate();
    p.bindLabelTexture();

    auto program = shader_.getShader();
    if(!program || !program->isLinked()) {
        LGL_ERROR;
        LERROR("Shader not compiled!");
        return;
    }
    program->activate();
    program->setUniform("projectionTex_", projUnit.getUnitNumber());
    program->setUniform("labelTex_", labelUnit.getUnitNumber());

    glDepthFunc(GL_ALWAYS);
    renderQuad();
    glDepthFunc(GL_LESS);

    program->deactivate();
    port.deactivateTarget();
    glActiveTexture(GL_TEXTURE0);
    LGL_ERROR;
}

InteractiveProjectionLabeling::~InteractiveProjectionLabeling() {
}

void InteractiveProjectionLabeling::process() {
    // initialize shader
    if (getInvalidationLevel() == INVALID_PROGRAM) {
        shader_.rebuild();
    }
    renderToPort(xyProjectionOutport_, xy_);
    renderToPort(xzProjectionOutport_, xz_);
    renderToPort(yzProjectionOutport_, yz_);


}
void InteractiveProjectionLabeling::adjustPropertiesToInput() {
    if(inport_.hasData()) {
        auto dim = inport_.getData()->getDimensions();

        VolumeAtomic<float> xyProjection(tgt::svec3(dim.x, dim.y, 1));
        VolumeAtomic<float> xzProjection(tgt::svec3(dim.x, dim.z, 1));
        VolumeAtomic<float> yzProjection(tgt::svec3(dim.y, dim.z, 1));

        for(size_t i=0; i < tgt::hmul(xyProjection.getDimensions()); ++i) {
            xyProjection.voxel(i) = -std::numeric_limits<float>::infinity();
        }
        for(size_t i=0; i < tgt::hmul(xzProjection.getDimensions()); ++i) {
            xzProjection.voxel(i) = -std::numeric_limits<float>::infinity();
        }
        for(size_t i=0; i < tgt::hmul(yzProjection.getDimensions()); ++i) {
            yzProjection.voxel(i) = -std::numeric_limits<float>::infinity();
        }

        for(size_t z=0; z < dim.z; ++z) {
            std::unique_ptr<VolumeRAM> slice(inport_.getData()->getSlice(z));
            for(size_t y=0; y < dim.y; ++y) {
                for(size_t x=0; x < dim.x; ++x) {
                    float volVoxel = slice->getVoxelNormalized(x, y, 0);
                    auto processVoxel = [volVoxel](float& current) {
                        current = std::max(volVoxel, current);
                    };
                    processVoxel(xyProjection.voxel(x, y, 0));
                    processVoxel(xzProjection.voxel(x, z, 0));
                    processVoxel(yzProjection.voxel(y, z, 0));
                }
            }
        }

        xy_ = LabelProjection(std::move(xyProjection));
        xz_ = LabelProjection(std::move(xzProjection));
        yz_ = LabelProjection(std::move(yzProjection));
    }
}
VoreenSerializableObject* InteractiveProjectionLabeling::create() const {
    return new InteractiveProjectionLabeling();
}

} // namespace voreen
