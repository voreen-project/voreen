/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2021 University of Muenster, Germany,           *
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

#ifndef TGT_GLCONTEXTBASE_H
#define TGT_GLCONTEXTBASE_H

#include "tgt/types.h" //for TGT_API
#include <string>

namespace tgt {

    /**
     * This is an interface for OpenGL context objects.
     * They are managed by tgt::GLContextManager.
     */
class TGT_API GLContextBase {
    friend class GLContextManager;
protected:
    /** Constructor */
    GLContextBase(const std::string& debugName);
    /** Destructor */
    virtual ~GLContextBase();

    /**
     * Makes "this" the current OpenGL context.
     * Called by the GLContextManager.
     * @note It's strongly recommended to avoid calling elsewhere, except you know what you are doing!
     */

    virtual void activate() = 0;

    /**
    * Determines, whether "this" is the current OpenGL context.
    * Called by the GLContextManager.
    */
    virtual bool isActive() = 0;
};

} //namespace

#endif //TGT_GLCONTEXTBASE_H
