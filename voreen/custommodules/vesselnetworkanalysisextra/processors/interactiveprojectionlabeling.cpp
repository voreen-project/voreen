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

LabelGuard::LabelGuard(LabelProjection& labelProjection)
    : labelProjection_(labelProjection)
{
}
LabelGuard::~LabelGuard() {
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    labelProjection_.labelTexture_->uploadTexture();
}
uint8_t& LabelGuard::at(tgt::svec3 p) {
    return labelProjection_.labels_.voxel(p);
}
void LabelGuard::set(size_t x, size_t y, tgt::svec2 range, uint8_t val) {
    for(size_t z=range.x; z<=range.y; ++z) {
        tgt::vec3 projected(x, y, z);
        tgt::svec3 real = labelProjection_.projectedToReal() * projected;
        labelProjection_.labels_.voxel(real) = val;
    }
}

LabelProjection::LabelProjection()
    : LabelProjection(tgt::svec3(2), tgt::mat3::identity)
{
}
LabelProjection::LabelProjection(tgt::svec3 dimensions, tgt::mat3 realToProjectedMat)
    : labels_(dimensions)
    , labelTexture_(boost::none)
    , realToProjectedMat_(realToProjectedMat)
{
    for(size_t i=0; i < tgt::hmul(labels_.getDimensions()); ++i) {
        labels_.voxel(i) = UNLABELED;
    }
}
void LabelProjection::ensureTexturesPresent() {
    if(!labelTexture_) {
        labelTexture_ = tgt::Texture(labels_.getDimensions(), GL_RED, GL_RED, GL_UNSIGNED_BYTE, tgt::Texture::NEAREST, tgt::Texture::CLAMP_TO_EDGE, (GLubyte*) labels_.voxel(), false);
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        labelTexture_->uploadTexture();
    }
}
void LabelProjection::bindLabelTexture() {
    ensureTexturesPresent();
    labelTexture_->bind();
}
void InteractiveProjectionLabeling::drawEvent(tgt::MouseEvent* e, LabelProjection& p) {
    auto button = e->button();
    if(e->modifiers() != tgt::Event::CTRL
            || ((button & (tgt::MouseEvent::MOUSE_BUTTON_LEFT | tgt::MouseEvent::MOUSE_BUTTON_MIDDLE | tgt::MouseEvent::MOUSE_BUTTON_RIGHT)) == 0))
             {
        return;
    }
    auto event_type = e->getEventType();
    tgt::ivec2 coords = e->coord();
    tgt::ivec2 viewport = e->viewport();

    tgt::vec3 projectionDimensions = p.realToProjected() * p.labels().getDimensions();

    coords.y = viewport.y - coords.y - 1;
    auto norm_coords = tgt::vec2(coords)/tgt::vec2(viewport);
    tgt::svec2 labelcoords = tgt::floor(norm_coords * projectionDimensions.xy());

    uint8_t label=LabelProjection::UNLABELED;
    if((button & tgt::MouseEvent::MOUSE_BUTTON_LEFT) == button) {
        label = LabelProjection::FOREGROUND;
    }
    if((button & tgt::MouseEvent::MOUSE_BUTTON_MIDDLE) == button) {
        label = LabelProjection::UNLABELED;
    }
    if((button & tgt::MouseEvent::MOUSE_BUTTON_RIGHT) == button) {
        label = LabelProjection::BACKGROUND;
    }
    p.labels_mut().set(labelcoords.x, labelcoords.y, projectionRange(p), label);
    if(event_type == tgt::MouseEvent::MOUSERELEASEEVENT) {
        syncLabels();
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
    , projectionRegion_("projectionRegion", "Projection Region")
    , xy_()
    , xz_()
    , yz_()
    , outputVolume_(boost::none)
{
    addPort(inport_);
    addPort(labelVolume_);
    addPort(xyProjectionOutport_);
    addPort(xzProjectionOutport_);
    addPort(yzProjectionOutport_);

    addProperty(shader_);
    addProperty(projectionRegion_);
}

void InteractiveProjectionLabeling::renderToPort(RenderPort& port, LabelProjection& p) {
    if(!inport_.hasData()) {
        return;
    }
    const auto& vol = *inport_.getData();

    port.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    tgt::TextureUnit volUnit;
    volUnit.activate();
    auto volgl = vol.getRepresentation<VolumeGL>();
    auto voltex = volgl->getTexture();
    voltex->bind();

    tgt::TextureUnit labelUnit;
    labelUnit.activate();
    p.bindLabelTexture();

    auto program = shader_.getShader();
    if(!program || !program->isLinked()) {
        LGL_ERROR;
        LERROR("Shader not compiled!");
        return;
    }
    tgt::ivec2 projectionRange(this->projectionRange(p));
    program->activate();
    program->setUniform("dimensions_", tgt::ivec3(vol.getDimensions()));
    program->setUniform("realToProjectedMat_", p.realToProjected());
    program->setUniform("projectedToRealMat_", p.projectedToReal());
    program->setUniform("projectionRange_", projectionRange);
    program->setUniform("volumeTex_", volUnit.getUnitNumber());
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
        std::string header;
        header += "#define UNLABELED " + std::to_string(LabelProjection::UNLABELED) + "\n";
        header += "#define BACKGROUND " + std::to_string(LabelProjection::BACKGROUND) + "\n";
        header += "#define FOREGROUND " + std::to_string(LabelProjection::FOREGROUND) + "\n";
        header += "#define SUGGESTED_FOREGROUND " + std::to_string(LabelProjection::SUGGESTED_FOREGROUND) + "\n";
        header += "#define INCONSISTENT " + std::to_string(LabelProjection::INCONSISTENT) + "\n";
        shader_.setHeader(header);
        shader_.rebuild();
    }
    renderToPort(xyProjectionOutport_, xy_);
    renderToPort(xzProjectionOutport_, xz_);
    renderToPort(yzProjectionOutport_, yz_);

    //auto* vol = new VolumeAtomic<uint8_t>(xz_.labels().getDimensions());
    //for(int i=0; i<tgt::hmul(xz_.labels().getDimensions()); ++i) {
    //    vol->voxel(i) = xz_.labels().voxel(i);
    //}
    //labelVolume_.setData(new Volume(vol, inport_.getData()));
}

tgt::svec2 InteractiveProjectionLabeling::projectionRange(LabelProjection& p) {
    return tgt::svec2(
            (p.realToProjected() * projectionRegion_.get().getLLF()).z,
            (p.realToProjected() * projectionRegion_.get().getURB()).z
            );
}
void InteractiveProjectionLabeling::syncLabels() {
    if(!inport_.hasData()) {
        return;
    }
    const auto& vol = *inport_.getData();
    auto dim = vol.getDimensions();

    auto xy = xy_.labels_mut();
    auto xz = xz_.labels_mut();
    auto yz = yz_.labels_mut();

    auto mark_suggested = [] (tgt::svec3 p, LabelGuard& m, LabelGuard& o1, LabelGuard& o2) {
        if(m.at(p) == LabelProjection::UNLABELED && o1.at(p) == LabelProjection::FOREGROUND && o2.at(p) == LabelProjection::FOREGROUND) {
            m.at(p) = LabelProjection::SUGGESTED_FOREGROUND;
        }
    };

    withOutputVolume([&] (LZ4SliceVolume<uint8_t>& vol) {
        for(size_t z=0; z < dim.z; ++z) {
            auto slice = vol.getWriteableSlice(z);
            for(size_t y=0; y < dim.y; ++y) {
                for(size_t x=0; x < dim.x; ++x) {
                    tgt::svec3 p(x,y,z);
                    mark_suggested(p, xy, yz, xz);
                    mark_suggested(p, yz, xy, xz);
                    mark_suggested(p, xz, yz, xy);

                    if(xy.at(p) == LabelProjection::FOREGROUND && yz.at(p) == LabelProjection::FOREGROUND && xz.at(p) == LabelProjection::FOREGROUND) {
                        slice->voxel(p.x, p.y, 0) = LabelProjection::FOREGROUND;
                    } else if(xy.at(p) == LabelProjection::BACKGROUND || yz.at(p) == LabelProjection::BACKGROUND || xz.at(p) == LabelProjection::BACKGROUND) {
                        slice->voxel(p.x, p.y, 0) = LabelProjection::BACKGROUND;
                    }
                }
            }
        }
    });
}
void InteractiveProjectionLabeling::withOutputVolume(std::function<void(LZ4SliceVolume<uint8_t>&)> func) {
    if(outputVolume_) {
        labelVolume_.setData(nullptr);
        func(*outputVolume_);
        labelVolume_.setData(LZ4SliceVolume<uint8_t>::open(outputVolume_->getFilePath()).toVolume().release());
    }
}

void InteractiveProjectionLabeling::adjustPropertiesToInput() {
    labelVolume_.setData(nullptr);
    if(!inport_.hasData()) {
        return;
    }
    const auto& vol = *inport_.getData();
    auto dim = vol.getDimensions();

    LZ4SliceVolumeBuilder<uint8_t> builder(
            VoreenApplication::app()->getUniqueTmpFilePath("." + LZ4SliceVolumeBase::FILE_EXTENSION),
            LZ4SliceVolumeMetadata::fromVolume(vol).withRealWorldMapping(RealWorldMapping::createDenormalizingMapping<uint8_t>())
            );
    builder.fill(LabelProjection::UNLABELED);
    outputVolume_ = std::move(builder).finalize();
    withOutputVolume([&] (LZ4SliceVolume<uint8_t>& vol) { });

    projectionRegion_.setMaxValue(dim - tgt::svec3(1));

    xy_ = LabelProjection(dim, tgt::mat3(1,0,0, 0,1,0, 0,0,1));
    xz_ = LabelProjection(dim, tgt::mat3(1,0,0, 0,0,1, 0,1,0));
    yz_ = LabelProjection(dim, tgt::mat3(0,1,0, 0,0,1, 1,0,0));
}
VoreenSerializableObject* InteractiveProjectionLabeling::create() const {
    return new InteractiveProjectionLabeling();
}

} // namespace voreen
