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

#include "componentidselector.h"

#include "voreen/core/voreenapplication.h"

#include "tgt/textureunit.h"

#include <sstream>

using tgt::TextureUnit;

namespace voreen {


ComponentIdSelector::ComponentIdSelector()
    : ImageProcessorBypassable("image/distance"), mouseDown_(false), mousePosition2D_(0)
    , raycastingResultInport_(Port::INPORT, "raycasting.result", "Raycasting Result Input", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_DEFAULT, GL_RGBA16F)
    , fhpInport_(Port::INPORT, "fhp", "First-hit-points Input", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_DEFAULT, GL_RGBA16F)
    , refInport_(Port::INPORT, "refvol", "Reference Volume", false)
    , outport_(Port::OUTPORT, "image.output", "Image Output")
    , selectedID_("selected.id", "Selected ID", 0, 0, std::numeric_limits<int>::max())
    , maxID_("max.id", "Max ID", 0, 0, std::numeric_limits<int>::max())
    , components_("components", "Selected Components")
    , selectedPosition_("selected.position", "Selected Position", tgt::ivec3(0), tgt::ivec3(0), tgt::ivec3(std::numeric_limits<int>::max()))
    , idSearchRadius_("id.searchradius", "ID Search Radius", 1, 0, 10)
{
    addPort(raycastingResultInport_);
    addPort(fhpInport_);
    addPort(refInport_);
    addPort(outport_);

    addProperty(selectedID_);
    addProperty(maxID_);
    addProperty(components_);
    addProperty(selectedPosition_);
    addProperty(idSearchRadius_);

    ON_CHANGE(maxID_, ComponentIdSelector, fillList);
}

ComponentIdSelector::~ComponentIdSelector() {
}

Processor* ComponentIdSelector::create() const {
    return new ComponentIdSelector();
}

bool ComponentIdSelector::isReady() const {
    if (!isInitialized() || !fhpInport_.isReady() || !outport_.isReady())
        return false;

    return true;
}

void ComponentIdSelector::process() {
    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM)
        compile();

    const VolumeBase* refVolume = refInport_.getData();
    if(!refVolume) {
        LERROR("No reference volume");
        return;
    }

    if (mouseDown_) {
        tgt::ivec2 viewport = outport_.getSize();

        tgt::ivec2 clampedMousePos2D = clampToViewport(tgt::ivec2(mousePosition2D_.x, viewport.y-mousePosition2D_.y));

        tgt::vec4 fhp;
        fhpInport_.activateTarget();
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        glReadPixels(clampedMousePos2D.x, clampedMousePos2D.y, 1, 1, GL_RGBA, GL_FLOAT, &fhp);

        if (length(fhp) > 0.0f) {
            tgt::vec4 voxelPos = refVolume->getTextureToVoxelMatrix() * fhp;

            selectedPosition_.set(voxelPos.xyz());

            auto rwm = refVolume->getRealWorldMapping();

            // I'm assuming the ID is in the red channel of the texture
            auto volumeData = VolumeRAMRepresentationLock(refVolume);
            int id = 0;
            int radius = idSearchRadius_.get();
            for (int z = -radius; z <= radius; ++z) {
                for (int y = -radius; y <= radius; ++y) {
                    for (int x = -radius; x <= radius; ++x) {
                        float sample = volumeData->getVoxelNormalized(voxelPos.xyz() + tgt::vec3(x, y, z));
                        int tempId = static_cast<int>(std::round(rwm.normalizedToRealWorld(sample)));
                        id = std::max(id, tempId);
                    }
                }
            }

            selectedID_.set(id);

            // Only one component is selected
            auto iter = selection_.find(selectedID_.get()-1);
            if(iter == selection_.end()) {
                selection_.insert(selectedID_.get() - 1);
            }
            else {
                selection_.erase(iter);
            }

            components_.setSelectedRowIndices(std::vector<int>(selection_.begin(), selection_.end()));

            invalidate();
        }

        mouseDown_ = false;
    }

    // dummy outport
    bypass(&raycastingResultInport_, &outport_);
}

void ComponentIdSelector::onEvent(tgt::Event* e) {
    tgt::MouseEvent* event = dynamic_cast<tgt::MouseEvent*>(e);
    if (event) mouseEvent(event);

    ImageProcessorBypassable::onEvent(e);
}

void ComponentIdSelector::mouseEvent(tgt::MouseEvent* e) {
    if (e->button() == tgt::MouseEvent::MOUSE_BUTTON_LEFT && e->action() == tgt::MouseEvent::RELEASED && !mouseDown_) {

        if(e->modifiers() != tgt::MouseEvent::SHIFT) {
            selection_.clear();
        }

        int x = e->x();
        int y = e->y();

        mouseDown_ = true;
        mousePosition2D_.x = x;
        mousePosition2D_.y = y;
    }

    // We always accept the event and force a redraw.
    e->accept();
    invalidate();
}

tgt::ivec2 ComponentIdSelector::clampToViewport(tgt::ivec2 mousePos) {
    tgt::ivec2 result = mousePos;
    tgt::ivec2 size = outport_.getSize();
    if (result.x < 0) result.x = 0;
    else if (result.x > size.x-1) result.x = size.x-1;
    if (result.y < 0) result.y = 0;
    else if (result.y > size.y-1) result.y = size.y-1;
    return result;
}

void ComponentIdSelector::fillList() {
    components_.reset();

    // Add all possible IDs up to max id
    int lastId = maxID_.get();
    for(int id = 0; id <= lastId; id++) {
        components_.addRow(std::to_string(id));
    }
}

void ComponentIdSelector::serialize(Serializer& s) const {
    RenderProcessor::serialize(s);
    //s.serialize("selection", selection_);
}
void ComponentIdSelector::deserialize(Deserializer& s) {
    RenderProcessor::deserialize(s);
    //s.optionalDeserialize("selection", selection_, std::set<int>{});
}

} // namespace voreen
