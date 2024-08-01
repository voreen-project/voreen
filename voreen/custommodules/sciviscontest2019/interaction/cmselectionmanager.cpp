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

#include "cmselectionmanager.h"
namespace voreen {
CMSelectionManager::CMSelectionManager()
    : InteractionHandler("dummy", "dummy")
    , noSelectionID_(-1)
    , mousePressedPos_(-1, -1)
{}
CMSelectionManager::CMSelectionManager(const std::string& id, const std::string& guiName, tgt::ivec2 size, IntProperty* selectionIDProp, IntProperty* mouseOverIDProp, const int noSelectionID)
    : InteractionHandler(id, guiName)
    , noSelectionID_(noSelectionID)
    , selectionIDProp_(selectionIDProp)
    , mouseOverIDProp_(mouseOverIDProp)
{
    moveEvent_ = new EventProperty<CMSelectionManager>(id + ".move", "Move", this,
        &CMSelectionManager::mouseMoveEvent,
        tgt::MouseEvent::MOUSE_BUTTON_NONE,
        tgt::MouseEvent::ACTION_ALL,
        tgt::Event::MODIFIER_NONE, true, true);
    addEventProperty(moveEvent_);

    clickEvent_ = new EventProperty<CMSelectionManager>(id + ".click", "Click", this,
        &CMSelectionManager::mousePressEvent,
        tgt::MouseEvent::MOUSE_BUTTON_LEFT,
        tgt::MouseEvent::CLICK,
        tgt::Event::MODIFIER_NONE, true, true);
    addEventProperty(clickEvent_);

    renderTarget_.resize(size);
}
CMSelectionManager::~CMSelectionManager() {
}
void CMSelectionManager::resize(tgt::ivec2 newsize) {
    renderTarget_.resize(newsize);
}
void CMSelectionManager::initialize() {
    renderTarget_.initialize(GL_R32UI, GL_DEPTH_COMPONENT32);
}
void CMSelectionManager::deinitialize() {
    renderTarget_.deinitialize();
}
void CMSelectionManager::activate() {
    renderTarget_.activateTarget("cmselectionmanager");
}
void CMSelectionManager::deactivate() {
    renderTarget_.deactivateTarget();
}
void CMSelectionManager::clear() {
    //glClearColor(???);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    LGL_ERROR;
}
int CMSelectionManager::idAt(tgt::ivec2 pos) {
    int id = noSelectionID_;
    //activate();
    //glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    //glReadPixels(pos.x, pos.y, 1, 1, GL_RED, GL_UNSIGNED_INT, &id);
    //deactivate();
    //tgt::vec4 col = renderTarget_.getColorAtPos(pos);
    //std::cout << col << std::endl;
    tgt::Texture* t = renderTarget_.getColorTexture();
    if(t == nullptr) {
        return id;
    }
    pos.y = t->getHeight() - pos.y;
    if(pos.x >= 0 && pos.y>=0 && pos.y < t->getHeight() && pos.x < t->getWidth()) {
        t->downloadTexture();
        //the fbo will be cleared to 0, but we want to return
        //NO_SELECTION_ID if there is no item at pos
        id = t->texel<GLuint>(pos) + noSelectionID_;
        LGL_ERROR;
    }
    return id;
}
void CMSelectionManager::mouseMoveEvent(tgt::MouseEvent* mouseEve) {
    if(mouseOverIDProp_) {
        tgt::ivec2 mousePos = tgt::ivec2(mouseEve->x(), mouseEve->y());
        int id = idAt(mousePos);
        if(id != mouseOverIDProp_->get()) {
            mouseOverIDProp_->set(id);
        }
    }
    mouseEve->accept();
}
void CMSelectionManager::mousePressEvent(tgt::MouseEvent* mouseEve) {
    if(selectionIDProp_) {
        tgt::ivec2 mousePos = tgt::ivec2(mouseEve->x(), mouseEve->y());
        if(mouseEve->action() == tgt::MouseEvent::PRESSED) {
            mousePressedPos_ = mousePos;
        } else if(mouseEve->action() == tgt::MouseEvent::RELEASED) {
            if(mousePressedPos_ == mousePos) {
                int id = idAt(mousePos);
                if(id != selectionIDProp_->get() && id != noSelectionID_) {
                    selectionIDProp_->set(id);
                }
            }
        }
    }
    mouseEve->accept();
}
void CMSelectionManager::onEvent(tgt::Event* eve) {
    // invalidate processor and update camera prop widgets, if event has been accepted
    if (eve->isAccepted()) {
        //center_->invalidate();
    }
}
} //namespace voreen
