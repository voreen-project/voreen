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

#ifndef VRN_TEXTURE2DARRAY_H
#define VRN_TEXTURE2DARRAY_H

#include "tgt/texture.h"
#include "tgt/types.h" //for api macro

namespace tgt {

/**
 * This class is basically a Texture, only the 2D_Array type is set in the constructor.
  */
class TGT_API Texture2DArray : public Texture {
public:
    /**
     * Calls the normal constructor and sets the type to TEXTURE_2D_ARRAY
     */
    Texture2DArray(const tgt::ivec3& dimensions, GLint format = GL_RGBA, GLint internalformat = GL_RGBA,
                   GLenum dataType = GL_UNSIGNED_BYTE, Filter filter = LINEAR, Wrapping wrapping = REPEAT,
                   GLubyte* data = 0, bool ownership = true);
};

} // namespace tgt

#endif // VRN_TEXTURE2DARRAY_H
