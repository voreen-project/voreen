/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2020 University of Muenster, Germany,           *
 * Department of Computer Science.                                    *
 *                                                                    *
 * This file is part of the tgt library. This library is free         *
 * software; you can redistribute it and/or modify it under the terms *
 * of the GNU Lesser General Public License version 2.1 as published  *
 * by the Free Software Foundation.                                   *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       *
 * GNU Lesser General Public License for more details.                *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License in the file "LICENSE.txt" along with this library.         *
 * If not, see <http://www.gnu.org/licenses/>.                        *
 *                                                                    *
 **********************************************************************/

#ifndef TGT_TEXTUREUNIT_H
#define TGT_TEXTUREUNIT_H

#include <string>
#include "tgt/tgt_gl.h"
#include "tgt/shadermanager.h"

namespace tgt {

/**
 * OpenGL Texture
 */
class TGT_API TextureUnit {
public:
    TextureUnit(bool keep = false);
    ~TextureUnit();

    void activate() const;

    GLint getEnum() const;
    GLint getUnitNumber() const;

    static void setZeroUnit();
    static void cleanup();

    static unsigned short numLocalActive();

    /**
     * Returns true if no texture units are currently assigned.
     */
    static bool unused();

protected:
    void assignUnit() const;
    static void init();

    mutable GLint number_;
    mutable GLint glEnum_;
    mutable bool assigned_;
    bool keep_;

    static bool initialized_;
    static unsigned short totalActive_;
    static unsigned short maxTexUnits_;
    static unsigned short numKeptUnits_;
    static std::vector<bool> busyUnits_;
};

} // namespace tgt

#endif // TGT_TEXTUREUNIT_H
