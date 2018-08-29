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

#include "multiview.h"

namespace voreen {

MultiView::MultiView()
    : RenderProcessor()
    , showGrid_("showGrid", "Show grid", true)
    , gridColor_("gridColor", "Grid color", tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f))
    , maximized_("maximized", "Maximized sub-view", 0, 0, 100000)
    , maximizeOnDoubleClick_("maximizeOnDoubleClick", "Maximize on double click", true)
    , maximizeEventProp_("mouseEvent.maximize", "Maximize Event", this, &MultiView::toggleMaximization,
        tgt::MouseEvent::MOUSE_BUTTON_LEFT, tgt::MouseEvent::DOUBLECLICK, tgt::MouseEvent::MODIFIER_NONE)
    , mouseMoveEventProp_("mouseEvent.move", "Move Event", this, &MultiView::mouseMove,
    tgt::MouseEvent::MOUSE_BUTTON_NONE, tgt::MouseEvent::MOTION | tgt::MouseEvent::CLICK | tgt::MouseEvent::ENTER_EXIT, tgt::MouseEvent::MODIFIER_NONE)
    , outport_(Port::OUTPORT, "outport", "Image Output")
    , inport_(Port::INPORT, "inport", "Image Input", true)
    , renderPorts_(std::vector<RenderPort*>())
    , gridResolution_(tgt::ivec2(0))
    , currentPort_(-1)
    , isDragging_(false)
{
    addProperty(showGrid_);
    addProperty(gridColor_);
    addProperty(maximized_);
    maximized_.setVisibleFlag(false);
    addProperty(maximizeOnDoubleClick_);
    addEventProperty(maximizeEventProp_);
    addEventProperty(mouseMoveEventProp_);

    addPort(outport_);
    addPort(inport_);
}

MultiView::~MultiView() {
}

Processor* MultiView::create() const {
    return new MultiView();
}

bool MultiView::isReady() const {
    if (!outport_.isReady())
        return false;

    if (!inport_.isReady())
        return false;

    if(maximized_.get() != 0 && !inport_.getConnected().at(maximized_.get() - 1)->isReady())
        return false;

    return true;
}

void MultiView::updateGridResolution() {
    renderPorts_.clear();

    for(size_t i = 0; i < inport_.getConnected().size(); i++) {
        const RenderPort* rp = dynamic_cast<const RenderPort*>(inport_.getConnected().at(i));
        if(!rp)
            continue;
        renderPorts_.push_back(const_cast<RenderPort*>(rp));
    }

    if(renderPorts_.empty())
        return;

    if(!procNameOrder_.empty() && portNameOrder_.size() == procNameOrder_.size() && renderPorts_.size() == procNameOrder_.size())
        establishConnectionOrder();
    else {
        procNameOrder_.clear();
        portNameOrder_.clear();
    }

    float root = std::sqrt(static_cast<float>(renderPorts_.size()));
    gridResolution_ = tgt::ivec2(tgt::ifloor(root));
    if(root > std::floor(root))
        gridResolution_.x += 1;
    if(std::floor(root + 0.5f) > std::floor(root))
        gridResolution_.y += 1;
}

void MultiView::establishConnectionOrder() {
    std::vector<RenderPort*> tmpPorts(renderPorts_.size(), 0);
    for(size_t i = 0; i < procNameOrder_.size(); i++) {
        for(size_t j = 0; j < renderPorts_.size(); j++) {
            if(procNameOrder_.at(i) == renderPorts_.at(j)->getProcessor()->getID() && portNameOrder_.at(i) == renderPorts_.at(j)->getID()) {
                tmpPorts[i] = renderPorts_.at(j);
                break;
            }
        }
    }

    for(size_t i = 0; i < tmpPorts.size(); i++) {
        if(tmpPorts.at(i) == 0) {
            procNameOrder_.clear();
            portNameOrder_.clear();
            return;
        }
    }

    renderPorts_ = tmpPorts;
}

void MultiView::beforeProcess() {
    if(inport_.hasChanged() && inport_.getConnected().size() != renderPorts_.size())
        updateSizes();
}

void MultiView::process() {
    if(renderPorts_.empty())
        return;

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.loadIdentity();
    MatStack.multMatrix(tgt::mat4::createOrtho(0.f, 1.f, 0.f, 1.f, 0.f, 1.f));

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();
    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (maximized_.get() == 0) {
        //tgt::vec2 outSize = tgt::vec2(outport_.getSize());

        tgt::vec2 invRes = tgt::vec2(1.f) / tgt::vec2(gridResolution_);

        for(size_t i = 0; i < renderPorts_.size(); i++) {
            RenderPort* target = renderPorts_.at(i);

            if(target->isReady()) {
                MatStack.pushMatrix();
                float locTransX = static_cast<float>(i - ((i / gridResolution_.x) * gridResolution_.x));
                float locTransY = static_cast<float>(gridResolution_.y - (i / gridResolution_.x) - 1);

                tgt::TextureUnit colUnit;
                colUnit.activate();
                target->getColorTexture()->enable();
                target->bindColorTexture();

                MatStack.scale(invRes.x, invRes.y, 1.0f);
                MatStack.translate(locTransX, locTransY, 0.0f);

                renderQuad();

                MatStack.loadIdentity();
                target->getColorTexture()->disable();
                MatStack.popMatrix();
            }
        }

        if(showGrid_.get()) {
            glDepthFunc(GL_ALWAYS);
            glColor4f(gridColor_.get().r, gridColor_.get().g, gridColor_.get().b, gridColor_.get().a);
            for(int i = 1; i < gridResolution_.x; i++) {
                glBegin(GL_LINES);
                    glVertex2f(i * invRes.x, 0.0f);
                    glVertex2f(i * invRes.x, 1.0f);
                glEnd();
            }
            for(int i = 1; i < gridResolution_.y; i++) {
                glBegin(GL_LINES);
                    glVertex2f(0.f, i * invRes.y);
                    glVertex2f(1.f, i * invRes.y);
                glEnd();
            }
            glDepthFunc(GL_LESS);
        }
    }
    else {
        //maximized:
        size_t targetIndex = maximized_.get() - 1;
        tgtAssert(targetIndex < renderPorts_.size(), "Index of maximized target too large");
        RenderPort * target = renderPorts_.at(targetIndex);

        if(target->isReady()) {
            tgt::TextureUnit colUnit;
            colUnit.activate();
            target->getColorTexture()->enable();
            target->bindColorTexture();

            renderQuad();

            target->getColorTexture()->disable();
        }
    }

    outport_.deactivateTarget();
    tgt::TextureUnit::setZeroUnit();
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();
    LGL_ERROR;
}

void MultiView::renderQuad() {
    glDepthFunc(GL_ALWAYS);
    glBegin(GL_QUADS);
        glTexCoord2f(0.f, 0.f);
        glVertex2f(0.f, 0.f);

        glTexCoord2f(1.f, 0.f);
        glVertex2f(1.f, 0.f);

        glTexCoord2f(1.f, 1.f);
        glVertex2f(1.f, 1.f);

        glTexCoord2f(0.f, 1.f);
        glVertex2f(0.f, 1.f);
    glEnd();
    glDepthFunc(GL_LESS);
}

void MultiView::initialize() {
    RenderProcessor::initialize();
}

/* TODO: size linking adaptation
void MultiView::portResized(RenderPort* p, tgt::ivec2 newsize) {
    outport_.resize(newsize);

    updateSizes();
    invalidate();
} */

void MultiView::updateSizes() {
    if (outport_.getSize() == tgt::ivec2(0))
        return;

    updateGridResolution();
    if(renderPorts_.empty())
        return;

    /* TODO: sizelinking adaptation
    if(maximized_.get() == 0) {
        tgt::ivec2 subsize = outport_.getSize() / gridResolution_;
        inport_.resize(subsize);
    } else {
        // TODO:  inport->resize(...) would additionally do "size_ = newsize;", does that matter?
        RenderPort* rp = renderPorts_.at(maximized_.get() - 1);
        static_cast<RenderProcessor*>(rp->getProcessor())->portResized(rp, outport_.getSize());
    } */
}

void MultiView::mouseMove(tgt::MouseEvent* e) {
    if(renderPorts_.empty())
        return;

    e->accept();
    int prevCurrentPort = currentPort_;

    if (maximized_.get() == 0) {
        if ((e->action() & tgt::MouseEvent::PRESSED) == tgt::MouseEvent::PRESSED)
            isDragging_ = true;
        if ((e->action() & tgt::MouseEvent::RELEASED) == tgt::MouseEvent::RELEASED)
            isDragging_ = false;

        if (!isDragging_) {
            int indexX = e->x() / (e->viewport().x / gridResolution_.x);
            int indexY = e->y() / (e->viewport().y / gridResolution_.y);
            currentPort_ = indexY * gridResolution_.x + indexX;
            if((size_t)currentPort_ >= renderPorts_.size() || (e->action() & tgt::MouseEvent::EXIT) == tgt::MouseEvent::EXIT)
                currentPort_ = -1;
        }

        RenderPort* curRP = (currentPort_ == -1) ? 0 : renderPorts_.at(currentPort_);

        if (currentPort_ != prevCurrentPort) {
            tgt::MouseEvent leaveEvent(1, 1, tgt::MouseEvent::EXIT, e->modifiers(), e->button(), e->viewport() / gridResolution_);
            tgt::MouseEvent enterEvent(1, 1, tgt::MouseEvent::ENTER, e->modifiers(), e->button(), e->viewport() / gridResolution_);
            leaveEvent.ignore();
            enterEvent.ignore();
            RenderPort* prevRP = (prevCurrentPort == -1) ? 0 : renderPorts_.at(prevCurrentPort);
            if(prevRP && prevRP->isReady())
                prevRP->distributeEvent(&leaveEvent);
            if(curRP && curRP->isReady())
                curRP->distributeEvent(&enterEvent);
        }

        if(curRP && curRP->isReady()) {
            tgt::MouseEvent moveEvent(e->x() % (e->viewport().x / gridResolution_.x),
                                      e->y() % (e->viewport().y / gridResolution_.y),
                                      tgt::MouseEvent::MOTION, e->modifiers(), e->button(), e->viewport() / gridResolution_);
            moveEvent.ignore();
            curRP->distributeEvent(&moveEvent);
        }

    } else {
        RenderPort* rp = renderPorts_.at(maximized_.get() - 1);
        if(rp && rp->isReady())
            rp->distributeEvent(e);
    }
}

void MultiView::invalidate(int inv) {
    RenderProcessor::invalidate(inv);
}

void MultiView::onEvent(tgt::Event* e) {

    tgt::MouseEvent* me = dynamic_cast<tgt::MouseEvent*>(e);

    if (!me || mouseMoveEventProp_.accepts(me) || (maximizeEventProp_.accepts(me) && maximizeOnDoubleClick_.get())) {
        RenderProcessor::onEvent(e);
        return;
    }

    if (maximized_.get() == 0) {
        int indexX = me->x() / (me->viewport().x / gridResolution_.x);
        int indexY = me->y() / (me->viewport().y / gridResolution_.y);
        size_t targetIndex = indexY * gridResolution_.x + indexX;
        if(targetIndex >= renderPorts_.size()) {
            RenderProcessor::onEvent(e);
            return;
        }

        RenderPort* rp = renderPorts_.at(targetIndex);
        if(rp && rp->isReady()) {
            tgt::MouseEvent newme(me->x() % (me->viewport().x / gridResolution_.x),
                                  me->y() % (me->viewport().y / gridResolution_.y),
                                  me->action(), me->modifiers(), me->button(), me->viewport() / gridResolution_);
            newme.ignore();  // accepted is set to true by default
            rp->distributeEvent(&newme);
            if (newme.isAccepted())
                me->accept();
        }
    }
    else {
        RenderPort* rp = renderPorts_.at(maximized_.get() - 1);
        if(rp && rp->isReady())
            rp->distributeEvent(me);
    }
}

void MultiView::toggleMaximization(tgt::MouseEvent* me) {

    if (maximizeOnDoubleClick_.get()) {
        if (maximized_.get() == 0) {
            int indexX = me->x() / (me->viewport().x / gridResolution_.x);
            int indexY = me->y() / (me->viewport().y / gridResolution_.y);
            size_t targetIndex = indexY * gridResolution_.x + indexX;
            if(targetIndex >= renderPorts_.size()) {
                me->ignore();
                return;
            }
            maximized_.set(static_cast<int>(targetIndex) + 1);
        } else {
            maximized_.set(0);
        }
        updateSizes();

        me->accept();
    }
}

void MultiView::serialize(Serializer& s) const {
    std::vector<std::string> procNames;
    std::vector<std::string> portNames;
    for(size_t i = 0; i < renderPorts_.size(); i++) {
        portNames.push_back(renderPorts_.at(i)->getID());
        procNames.push_back(renderPorts_.at(i)->getProcessor()->getID());
    }

    s.serialize("PortNames", portNames, "PortName");
    s.serialize("ProcNames", procNames, "ProcName");

    Processor::serialize(s);
}

void MultiView::deserialize(Deserializer& s) {

    try {
        s.deserialize("PortNames", portNameOrder_, "PortName");
        s.deserialize("ProcNames", procNameOrder_, "ProcName");
    } catch (SerializationNoSuchDataException) {
        s.removeLastError();
    }

    Processor::deserialize(s);
}

} // namespace voreen
