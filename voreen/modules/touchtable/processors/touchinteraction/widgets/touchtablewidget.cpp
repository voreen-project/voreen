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

#include "touchtablewidget.h"

namespace voreen {

TouchTableWidget::TouchTableWidget()
    : Processor()
    , port_(Port::OUTPORT, "port", "Port", false, Processor::VALID)
    , isActive_("isActive", "is Active", false)
    , position_("position", "Position", tgt::vec2(200, 200), tgt::vec2(0,0), tgt::vec2(1920, 1080))
    , oldPosition_("oldposition", "Old Position", tgt::vec2(200, 200), tgt::vec2(0,0), tgt::vec2(1920, 1080))
    , inMotion_("motion", "In Motion", false)
    , invalidPosition_("invalidPosition", "Position invalid", false)
    , isOnTable_("ontable", "Is on Table", false)
    , symbolFile_("symbolfile", "Symbol Image File", "Image File for Widget Symbol", VoreenApplication::app()->getModulePath("touchtable") + "/textures", "png-Files (*.png)")
    , symbolTex_(0)
    , overlay_(0)
{
    addPort(port_);

    addProperty(isActive_);
    addProperty(position_);

    addProperty(oldPosition_);
    oldPosition_.setVisibleFlag(false);

    addProperty(isOnTable_);
    isOnTable_.setVisibleFlag(false);

    addProperty(inMotion_);
    inMotion_.setVisibleFlag(false);

    addProperty(invalidPosition_);
    invalidPosition_.setVisibleFlag(false);

    addProperty(symbolFile_);
    // Module path only available once Application is initialized.
    // However object factories need an instance beforehand.
    if(VoreenApplication::app()->isInitialized())
        symbolFile_.set(VoreenApplication::app()->getModulePath("touchtable") + "/textures" + "/nosymbol.png");
}

void TouchTableWidget::initialize() {
    Processor::initialize();

    //set default texture
    TexMgr.addPath(VoreenApplication::app()->getModulePath("touchtable") + "/textures");
    symbolTex_ = TexMgr.load("nosymbol.png");
}

void TouchTableWidget::deinitialize() {

    TexMgr.dispose(symbolTex_);
    symbolTex_ = 0;

    Processor::deinitialize();
}

std::string TouchTableWidget::getClassName() const {
    return "TouchTableWidget";
}

std::string TouchTableWidget::getCategory() const {
    return "Touch Table";
}

void TouchTableWidget::setDescriptions() {
    setDescription("Abstract base class of widgets for use with TouchScreenOverlay.");
}

bool TouchTableWidget::isReady() const {
    return port_.isConnected();
}

void TouchTableWidget::process() {
    if (!port_.isConnected())
       overlay_ = 0;
}

bool TouchTableWidget::isActive() const {
    return isActive_.get();
}

void TouchTableWidget::setActive(bool active) {
    isActive_.set(active);
}

tgt::vec2 TouchTableWidget::getPosition() const {
    return position_.get();
}

void TouchTableWidget::setPosition(tgt::vec2 pos) {
    position_.set(pos);
}

tgt::vec2 TouchTableWidget::getOldPosition() const {
    return oldPosition_.get();
}

void TouchTableWidget::setOldPosition(tgt::vec2 pos) {
    oldPosition_.set(pos);
}

tgt::Texture* TouchTableWidget::getSymbolTexture() const {
    return symbolTex_;
}

bool TouchTableWidget::isInMotion() const{
    return inMotion_.get();
}

void TouchTableWidget::setInMotion(bool inMotion){
    inMotion_.set(inMotion);
}

bool TouchTableWidget::isPositionInvalid() const{
    return invalidPosition_.get();
}

void TouchTableWidget::setPositionInvalid(bool invalid){
    invalidPosition_.set(invalid);
}

void TouchTableWidget::setOnTable(bool t) {
    isOnTable_.set(t);
}
bool TouchTableWidget::isOnTable() const {
    return isOnTable_.get();
}

void TouchTableWidget::invalidate(int inv) {

    if (!isInitialized())
        return;

    if (!symbolTex_ || (symbolFile_.get() != symbolTex_->getOptionalName())) {
        TexMgr.dispose(symbolTex_);
        symbolTex_ = 0;
        symbolTex_ = TexMgr.load(symbolFile_.get());
    }

    Processor::invalidate(inv);
}

void TouchTableWidget::setOverlay(TouchTableOverlay *overlay) {
    overlay_ = overlay;
}


} // namespace
