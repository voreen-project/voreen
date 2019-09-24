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

#include "splitter.h"
#include "tgt/immediatemode/immediatemode.h"


namespace voreen {

const int Splitter::HANDLE_GRAB_TOLERANCE = 5;

Splitter::Splitter()
    : RenderProcessor()
    , showGrid_("showGrid", "Show grid", true)
    , gridColor_("gridColor", "Grid color", tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f))
    , lineWidth_("lineWidth", "Line Width", 1.0f, 0.0f, 10.0f)
    , vertical_("vertical", "Orientation", true)
    , overlay_("overlay", "Overlay", false)
    , fixSplitPosition_("fixSplitPosition", "Fix Split Position", false)
    , maximized_("maximized", "Maximized sub-view", 0, 0, 4)
    , maximizeOnDoubleClick_("maximizeOnDoubleClick", "Maximize on double click", true)
    , maximizeEventProp_("mouseEvent.maximize", "Maximize Event", this, &Splitter::toggleMaximization,
    tgt::MouseEvent::MOUSE_BUTTON_LEFT, tgt::MouseEvent::DOUBLECLICK, tgt::MouseEvent::MODIFIER_NONE)
    , position_("position", "Position", 0.5f, 0.0f, 1.0f)
    , outport_(Port::OUTPORT, "outport", "Image Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , inport1_(Port::INPORT, "inport1", "Image1 Input", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN)
    , inport2_(Port::INPORT, "inport2", "Image2 Input", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN)
    , currentPort_(-1)
    , isDragging_(false)
    , lastWindow_(1)
{
    vertical_.addOption("vertical", "Vertical", true);
    vertical_.addOption("horizontal", "Horizontal", false);

    addProperty(showGrid_);
    addProperty(gridColor_);
    addProperty(lineWidth_);
    addProperty(vertical_);
    addProperty(position_);
    addProperty(overlay_);
    addProperty(fixSplitPosition_);
    ON_CHANGE(overlay_, Splitter, updateSizes);
    position_.onChange( MemberFunctionCallback<Splitter>(this, &Splitter::updateSizes));
    vertical_.onChange( MemberFunctionCallback<Splitter>(this, &Splitter::updateSizes));
    maximized_.onChange( MemberFunctionCallback<Splitter>(this, &Splitter::updateSizes));
    overlay_.onChange( MemberFunctionCallback<Splitter>(this, &Splitter::updateSizes));

    addProperty(maximized_);
    maximized_.setVisibleFlag(false);
    addProperty(maximizeOnDoubleClick_);
    addEventProperty(maximizeEventProp_);
    addPort(outport_);
    addPort(inport1_);
    addPort(inport2_);

    outport_.onSizeReceiveChange<Splitter>(this, &Splitter::updateSizes);
}

Splitter::~Splitter() {
}

Processor* Splitter::create() const {
    return new Splitter();
}

bool Splitter::isReady() const {
    return (inport1_.isReady() || inport2_.isReady()) && outport_.isReady();
}


void Splitter::process() {
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    outport_.activateTarget();

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    WindowConfiguration config = getWindowConfiguration();
    float position = position_.get();
    tgt::vec2 viewport = tgt::vec2(outport_.getSize());
    if (inport1_.isReady()) {
        tgt::vec2 minPx, maxPx, textureMin = tgt::vec2::zero, textureMax = tgt::vec2::one;
        switch (config)
        {
        case voreen::Splitter::FIRST:
            minPx = tgt::vec2::zero;
            maxPx = tgt::vec2(viewport);
            break;
        case voreen::Splitter::SECOND:
            minPx = tgt::vec2::zero;
            maxPx = tgt::vec2::zero;
            break;
        case voreen::Splitter::HORIZONTAL:
            minPx = tgt::vec2(0, viewport.y*position);
            maxPx = tgt::vec2(viewport);
            break;
        case voreen::Splitter::VERTICAL:
            minPx = tgt::vec2::zero;
            maxPx = tgt::vec2(viewport.x*position, viewport.y);
            break;
        case voreen::Splitter::HORIZONTAL_OVERLAY:
            textureMin = tgt::vec2(0, position);
            minPx = tgt::vec2(0, viewport.y*position);
            maxPx = tgt::vec2(viewport);
            break;
        case voreen::Splitter::VERTICAL_OVERLAY:
            textureMax = tgt::vec2(position, 1);
            minPx = tgt::vec2::zero;
            maxPx = tgt::vec2(viewport.x*position, viewport.y);
            break;
        default:
            break;
        }
        renderSubview(&inport1_, minPx, maxPx, textureMin, textureMax);
    }

    if (inport2_.isReady()) {
        tgt::vec2 minPx, maxPx, textureMin = tgt::vec2::zero, textureMax = tgt::vec2::one;
        switch (config)
        {
        case voreen::Splitter::FIRST:
            minPx = tgt::vec2::zero;
            maxPx = tgt::vec2::zero;
            break;
        case voreen::Splitter::SECOND:
            minPx = tgt::vec2::zero;
            maxPx = tgt::vec2(viewport);
            break;
        case voreen::Splitter::HORIZONTAL:
            minPx = tgt::vec2::zero;
            maxPx = tgt::vec2(viewport.x, viewport.y*position);
            break;
        case voreen::Splitter::VERTICAL:
            minPx = tgt::vec2(viewport.x*position, 0);
            maxPx = tgt::vec2(viewport);

            break;
        case voreen::Splitter::HORIZONTAL_OVERLAY:
            textureMax = tgt::vec2(1, position);
            minPx = tgt::vec2::zero;
            maxPx = tgt::vec2(viewport.x, viewport.y*position);
            break;
        case voreen::Splitter::VERTICAL_OVERLAY:
            textureMin = tgt::vec2(position, 0);
            minPx = tgt::vec2(viewport.x*position, 0);
            maxPx = tgt::vec2(viewport);
            break;
        default:
            break;
        }
        renderSubview(&inport2_, minPx, maxPx, textureMin, textureMax);
    }
    //shader_->deactivate();

    if (maximized_.get() == 0){
        glActiveTexture(GL_TEXTURE0);
        if(showGrid_.get()) {
            glDepthFunc(GL_ALWAYS);
            IMode.color(gridColor_.get().r, gridColor_.get().g, gridColor_.get().b, gridColor_.get().a);
            glLineWidth(lineWidth_.get());
            IMode.begin(tgt::ImmediateMode::LINES);
            if(vertical_.getValue()) {
                IMode.vertex((position_.get() * 2.0f) - 1.0f, -1.0f);
                IMode.vertex((position_.get() * 2.0f) - 1.0f, 1.0f);
            }
            else {
                IMode.vertex(1.0f, (position_.get() * 2.0f) - 1.0f);
                IMode.vertex(-1.0f, (position_.get() * 2.0f) - 1.0f);
            }
            IMode.end();
            IMode.color(1.f,1.f,1.f,1.f);
            glLineWidth(1.0f);
            glDepthFunc(GL_LESS);
        }
    }

    outport_.deactivateTarget();

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();
    LGL_ERROR;
}

void Splitter::initialize() {
    RenderProcessor::initialize();
    updateSizes();
    //shader_ = ShdrMgr.loadSeparate("passthrough.vert", "copyimage.frag", generateHeader(), false);
}

void Splitter::deinitialize() {

    //ShdrMgr.dispose(shader_);

    RenderProcessor::deinitialize();
}

void Splitter::updateSizes() {
    if( (outport_.getReceivedSize().x == 0) || (outport_.getReceivedSize().y == 0) )
        return;

    // Retrieve the min and max size for render targets.
    tgt::ivec2 minValue = outport_.getSizeReceiveProperty()->getMinValue();
    tgt::ivec2 maxValue = outport_.getSizeReceiveProperty()->getMaxValue();

    tgt::ivec2 size = outport_.getReceivedSize();
    tgt::ivec2 subsize1 = size;
    tgt::ivec2 subsize2 = size;
    float position = position_.get();
    switch (getWindowConfiguration())
    {
    case voreen::Splitter::HORIZONTAL:
        subsize2.y = tgt::clamp(static_cast<int>(subsize1.y * position), minValue.y, std::min(size.y, maxValue.y));
        subsize1.y = tgt::clamp(size.y - subsize2.y, minValue.y, std::min(size.y, maxValue.y));
        break;
    case voreen::Splitter::VERTICAL:
        subsize1.x = tgt::clamp(static_cast<int>(subsize1.x * position), minValue.x, std::min(size.x, maxValue.x));
        subsize2.x = tgt::clamp(size.x - subsize1.x, minValue.x, std::min(size.x, maxValue.x));
        break;
    default:
        break;
    }

    inport1_.requestSize(subsize1);
    inport2_.requestSize(subsize2);
}

void Splitter::mouseEvent(tgt::MouseEvent* e) {

    if (maximizeEventProp_.accepts(e)) {
        RenderProcessor::onEvent(e);
        return;
    }

    tgt::MouseEvent::MouseAction action = e->action();
    tgt::MouseEvent::Modifier modifier = e->modifiers();
    tgt::ivec2 viewport = e->viewport();
    int x = e->x();
    int y = viewport.y - e->y();
    float position = position_.get();

    if (maximized_.get() == 0 && !fixSplitPosition_.get()) {
        if (!vertical_.getValue()) {
            if (action == tgt::MouseEvent::PRESSED && modifier == tgt::MouseEvent::MODIFIER_NONE && std::abs(y - viewport.y * position) < HANDLE_GRAB_TOLERANCE) {
                isDragging_ = true;
                e->accept();
                return;
            }
            if (isDragging_ && (action == tgt::MouseEvent::RELEASED || action == tgt::MouseEvent::EXIT)) {
                isDragging_ = false;
                e->accept();
                return;
            }

            if (isDragging_) {
                position_.set(tgt::clamp(1.0f*y / viewport.y, 0.0f, 1.0f));
                return;
            }
        }
        else {

            if (action == tgt::MouseEvent::PRESSED && modifier == tgt::MouseEvent::MODIFIER_NONE && std::abs(x - viewport.x * position) < HANDLE_GRAB_TOLERANCE) {
                isDragging_ = true;
                e->accept();
                return;
            }
            if (isDragging_ && (action == tgt::MouseEvent::RELEASED || action == tgt::MouseEvent::EXIT)) {
                isDragging_ = false;
                e->accept();
                return;
            }

            if (isDragging_) {
                position_.set(tgt::clamp(1.0f*x / viewport.x, 0.0f, 1.0f));
                return;
            }
        }
    }

    tgt::MouseEvent newme(0, 0, tgt::MouseEvent::ACTION_NONE, tgt::Event::Modifier::MODIFIER_NONE); // no empty constructor
    int window = getWindowForEvent(*e, &newme);
    if (lastWindow_ != window && !overlay_.get()){
        tgt::MouseEvent leaveEvent(1, 1, tgt::MouseEvent::EXIT, e->modifiers(), e->button(), getWindowViewport(getWindowConfiguration(), lastWindow_));
        tgt::MouseEvent enterEvent(1, 1, tgt::MouseEvent::ENTER, e->modifiers(), e->button(), getWindowViewport(getWindowConfiguration(), window));
        leaveEvent.ignore();
        enterEvent.ignore();
        distributeMouseEvent(lastWindow_, &leaveEvent);
        distributeMouseEvent(window, &enterEvent);
    }
    lastWindow_ = window;
    distributeMouseEvent(window, &newme);
    if (newme.isAccepted())
        e->accept();
}

void Splitter::onEvent(tgt::Event* e) {
    tgt::MouseEvent* me = dynamic_cast<tgt::MouseEvent*>(e);

    if (!me)
        RenderProcessor::onEvent(e);
    else
        mouseEvent(me);
}

Splitter::WindowConfiguration Splitter::getWindowConfiguration() const{
    if (maximized_.get() != 0)
    {
        switch (maximized_.get())
        {
        case 1:
            return FIRST;
        case 2:
            return SECOND;
        default:
            tgtAssert(false, "Splitter: Should not happen!")
            return FIRST;
        }
    }
    bool vertical = vertical_.getValue();
    bool overlay = overlay_.get();
    if (overlay){
        if (vertical)
            return VERTICAL_OVERLAY;
        else
            return HORIZONTAL_OVERLAY;
    }else{
        if (vertical)
            return VERTICAL;
        else
            return HORIZONTAL;
    }
}

int Splitter::getWindowForEvent(tgt::MouseEvent me, tgt::MouseEvent *translatedme){
    tgt::ivec2          viewport  = me.viewport();
    float               position  = position_.get();
    WindowConfiguration config    = getWindowConfiguration();
    bool                translate = (config == VERTICAL) || (config == HORIZONTAL);
    int                 vsplit    = static_cast<int>(position*viewport.x);
    int                 hsplit    = static_cast<int>((1.0f - position /* Different coordinate systems: origin in bottom or top left */)*viewport.y);
    switch (getWindowConfiguration())
    {
    case FIRST:
        if (translatedme)
            *translatedme = me;
        return 1;
    case SECOND:
        if (translatedme)
            *translatedme = me;
        return 2;
    case VERTICAL:
    case VERTICAL_OVERLAY:
        if (vsplit > me.x()) {
            // left side
            if (translatedme) {
                if (translate)
                    *translatedme = tgt::MouseEvent(me.x(), me.y(), me.action(), me.modifiers(), me.button(), tgt::ivec2(vsplit, viewport.y));
                else
                    *translatedme = me;
            }
            return 1;
        }
        else {
            // right side
            if (translatedme) {
                if (translate)
                    *translatedme = tgt::MouseEvent(me.x()-vsplit, me.y(), me.action(), me.modifiers(), me.button(), tgt::ivec2(viewport.x-vsplit, viewport.y));
                else
                    *translatedme = me;
            }
            return 2;
        }
    case HORIZONTAL:
    case HORIZONTAL_OVERLAY:
        if (hsplit > me.y()) {
            // top side
            if (translatedme) {
                if (translate)
                    *translatedme = tgt::MouseEvent(me.x(), me.y(), me.action(), me.modifiers(), me.button(), tgt::ivec2(viewport.x, hsplit));
                else
                    *translatedme = me;
            }
            return 1;
        }
        else {
            // bottom side
            if (translatedme) {
                if (translate)
                    *translatedme = tgt::MouseEvent(me.x(), me.y()-hsplit, me.action(), me.modifiers(), me.button(), tgt::ivec2(viewport.x, viewport.y-hsplit));
                else
                    *translatedme = me;
            }
            return 2;
        }
    default:
        break;
    }
    tgtAssert(false, "Splitter: Should not happen");
    return FIRST;
}

void Splitter::renderSubview(RenderPort* image, tgt::vec2 minPx, tgt::vec2 maxPx, tgt::vec2 textureMin, tgt::vec2 textureMax){
    if (!image->isReady())
        return;

    image->bindColorTexture(GL_TEXTURE0);
    image->getColorTexture()->enable();

    tgt::ivec2 viewport = outport_.getSize();

    MatStack.ortho(0, viewport.x, 0, viewport.y, -1, 1);

    glDepthFunc(GL_ALWAYS);
    IMode.setTextureMode(tgt::ImmediateMode::TEX2D);
    IMode.begin(tgt::ImmediateMode::TRIANGLE_STRIP);
        IMode.texcoord(textureMax.x, textureMin.y); IMode.vertex(maxPx.x, minPx.y);
        IMode.texcoord(textureMin); IMode.vertex(minPx);
        IMode.texcoord(textureMax); IMode.vertex(maxPx);
        IMode.texcoord(textureMin.x, textureMax.y); IMode.vertex(minPx.x, maxPx.y);
    IMode.end();
    glDepthFunc(GL_LESS);
    IMode.setTextureMode(tgt::ImmediateMode::TEXNONE);
    MatStack.loadIdentity();
    image->getColorTexture()->disable();
}

void Splitter::toggleMaximization(tgt::MouseEvent* me){
    if (!maximizeOnDoubleClick_.get())
        return;
    if (maximized_.get() != 0){
        maximized_.set(0);
        return;
    }
    int window = getWindowForEvent(*me, nullptr);
    maximized_.set(window);
}

void Splitter::distributeMouseEvent( int window, tgt::MouseEvent *newme ){
    RenderPort* distributionPort = 0;
    switch (window)
    {
    case 1:
        distributionPort = &inport1_;
        break;
    case 2:
        distributionPort = &inport2_;
        break;
    default:
        tgtAssert(false, "Splitter: Should not happen!")
        break;
    }
    newme->ignore();
    distributionPort->distributeEvent(newme);
}

tgt::ivec2 Splitter::getWindowViewport( WindowConfiguration windowConfiguration, int window )
{
    tgt::ivec2 size = outport_.getReceivedSize();
    float position = position_.get();
    switch (getWindowConfiguration())
    {
    case  FIRST:
        return (window==1)?size:tgt::ivec2::zero;
    case SECOND:
        return (window==2)?size:tgt::ivec2::zero;
    case VERTICAL:
        return tgt::ivec2(((window==1)?position:(1-position))*size.x, size.y);
    case VERTICAL_OVERLAY:
        return size;
    case HORIZONTAL:
        return tgt::ivec2(size.x, ((window==1)?position:(1-position))*size.y);
    case HORIZONTAL_OVERLAY:
        return size;
    default:
        return tgt::ivec2::zero;
    }
}

} // namespace voreen
