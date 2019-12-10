/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "voreen/core/interaction/idmanager.h"

#include "voreen/core/processors/processor.h"

#include "tgt/logmanager.h"
#include "tgt/immediatemode/immediatemode.h"

namespace voreen {

using tgt::col4;
using tgt::vec4;
using tgt::ivec2;

IDManager::IDManager() :
      rt_(0)
       , currentID_(0,0,0,255)
{
}

void IDManager::initializeTarget() {
    if (!rt_)
        return;

    rt_->initialize(GL_RGBA);
}

col4 IDManager::registerObject(const void* obj) {
    increaseID();
    colorToID_[currentID_] = obj;
    IDToColor_[obj] = currentID_;
    return currentID_;
}

void IDManager::deregisterObject(const void* obj) {
    col4 col = getColorFromObject(obj);
    colorToID_.erase(col);
    IDToColor_.erase(obj);
}

void IDManager::deregisterObject(const tgt::col4& col) {
    const void* obj = getObjectFromColor(col);
    colorToID_.erase(col);
    IDToColor_.erase(obj);
}

void IDManager::clearRegisteredObjects() {
    colorToID_.clear();
    IDToColor_.clear();
    currentID_ = col4(0,0,0,255);
}

col4 IDManager::getColorFromId(int id) {
    if(id < 0)
        return col4(0,0,0,0);

    if(id > 16777216) {
        LERRORC("voreen.IDManager", "id to big!");
        return col4(0,0,0,0);
    }

    col4 c;
    c.b = id & 255;
    c.g = (id >> 8) & 255;
    c.r = (id >> 16) & 255;
    c.a = 255;
    return c;
}

col4 IDManager::getColorFromObject(const void* obj) {
    if (obj == 0)
        return col4(0,0,0,0);

    if (isRegistered(obj))
        return IDToColor_[obj];
    else
        return col4(0,0,0,0);
}

vec4 IDManager::getColorFromObjectFloat(const void *obj) {
    col4 res = getColorFromObject(obj);
    vec4 f = vec4(IMode.normalize(res.x), IMode.normalize(res.y), IMode.normalize(res.z), IMode.normalize(res.w));
    return f;
}

void IDManager::setGLColor(const void* obj) {
    col4 col = getColorFromObject(obj);
    IMode.color(IMode.normalize(col));
}

void IDManager::setGLColor(int id) {
    col4 col = getColorFromId(id);
    IMode.color(IMode.normalize(col));
}

int IDManager::getIdFromColor(tgt::col4 col) {
    if (col.a == 0)
        return -1;
    else
        return (col.r * 256 * 256) + (col.g * 256) + col.b;
}

const void* IDManager::getObjectFromColor(tgt::col4 col) {
    if (col.a == 0)
        return 0;

    if (isRegistered(col))
        return colorToID_[col];
    else
        return 0;
}

int IDManager::getIdAtPos(tgt::ivec2 pos) {
    return getIdFromColor(getColorAtPos(pos));
}

const void* IDManager::getObjectAtPos(tgt::ivec2 pos) {
    return getObjectFromColor(getColorAtPos(pos));
}

col4 IDManager::getColorAtPos(tgt::ivec2 pos) {
    rt_->activateTarget();
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    col4 pixels;
    glReadPixels(pos.x, pos.y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, &pixels);
    rt_->deactivateTarget();
    return pixels;
}

bool IDManager::isHit(ivec2 pos, const void* obj) {
    return (getObjectAtPos(pos) == obj);
}

void IDManager::activateTarget(std::string debugLabel) {
    if (rt_) {
        rt_->activateTarget();
        rt_->increaseNumUpdates();
        rt_->setDebugLabel("ID target" + (debugLabel.empty() ? "" : " (" + debugLabel + ")"));
    }
    else
        LERRORC("voreen.idmanager", "No RenderTarget set!");
}

void IDManager::deactivateTarget() {
    if (rt_)
        rt_->deactivateTarget();
    else
        LERRORC("voreen.idmanager", "No RenderTarget set!");
}

void IDManager::clearTarget() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void IDManager::setRenderTarget(RenderTarget* rt) {
    rt_ = rt;
}

RenderTarget* IDManager::getRenderTarget() {
    return rt_;
}

void IDManager::increaseID() {
    if (currentID_.b == 255) {
        currentID_.b = 0;
        if (currentID_.g == 255) {
            currentID_.g = 0;
            if (currentID_.r == 255) {
                LERRORC("voreen.idmanager", "Out of ids...");
            }
            else
                currentID_.r++;

        }
        else
            currentID_.g++;
    }
    else
        currentID_.b++;
}

bool IDManager::isRegistered(const void* obj) {
    if(IDToColor_.find(obj) != IDToColor_.end())
        return true;
    else
        return false;
}

bool IDManager::isRegistered(tgt::col4 col) {
    if (colorToID_.find(col) != colorToID_.end())
        return true;
    else
        return false;
}

} // namespace voreen
