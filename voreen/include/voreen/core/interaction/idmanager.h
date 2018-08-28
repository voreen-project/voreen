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

#ifndef VRN_IDMANAGER_H
#define VRN_IDMANAGER_H

#include "tgt/vector.h"
#include "voreen/core/datastructures/rendertarget/rendertarget.h"

namespace voreen {

struct colComp {
    bool operator()(tgt::col4 c1, tgt::col4 c2) const {
        if (c1.r < c2.r)
            return true;
        else if (c1.r == c2.r) {
            if (c1.g < c2.g)
                return true;
            else if (c1.g == c2.g) {
                if (c1.b < c2.b)
                    return true;
                else
                    return false;
            }
            else
                return false;
        }
        else
            return false;
    }
};

///IDManager class for picking
class VRN_CORE_API IDManager {
public:
    IDManager();
    void initializeTarget();

    tgt::col4 registerObject(const void* obj);
    void deregisterObject(const void* obj);
    void deregisterObject(const tgt::col4& col);
    void clearRegisteredObjects();

    tgt::col4 getColorFromObject(const void* obj);
    tgt::vec4 getColorFromObjectFloat(const void* obj);
    tgt::col4 getColorFromId(int id);

    void setGLColor(const void* obj);
    void setGLColor(int obj);

    int getIdFromColor(tgt::col4 col);
    const void* getObjectFromColor(tgt::col4 col);

    int getIdAtPos(tgt::ivec2 pos);
    const void* getObjectAtPos(tgt::ivec2 pos);
    tgt::col4 getColorAtPos(tgt::ivec2 pos);

    bool isHit(tgt::ivec2 pos, const void* obj);

    void activateTarget(std::string debugLabel = "");
    void deactivateTarget();
    void clearTarget();
    void setRenderTarget(RenderTarget* rt);
    RenderTarget* getRenderTarget();

    bool isRegistered(const void* obj);
    bool isRegistered(tgt::col4 col);
private:
    void increaseID();

    std::map<tgt::col4, const void*, colComp> colorToID_;
    std::map<const void*, tgt::col4> IDToColor_;
    RenderTarget* rt_;
    tgt::col4 currentID_;
};


} //namespace voreen

#endif
