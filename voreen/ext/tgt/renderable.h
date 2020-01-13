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

#ifndef TGT_RENDERABLE_H
#define TGT_RENDERABLE_H

#include "tgt/bounds.h"
#include "tgt/types.h"

namespace tgt {

/**
 * Inherit from this abstract class in order to create a renderable object
*/
class TGT_API Renderable : public HasBounds {
protected:

    bool    static_;
    bool    visible_;

public:

    /**
     * The constructor sets basic properties of the renderable object.
     * Static objects are not expected to move
     * and groups of them can be cached in
     * displaylists or vertexarrays.
     * @param _static  static objects do not move
     * @param visible  Invisible objects are not rendered
    */
    Renderable(const Bounds& bounds, bool _static, bool visible = true)
      : HasBounds(bounds),
        static_(_static),
        visible_(visible)
    {}

    /**
     * virtual destructor to prevent memory leaks
    */
    virtual ~Renderable() {}

    /**
     * This method should provide the functionality
     * to init the rendering process.
     * You can init display lists here for example.
    */
    virtual void init() = 0;

    /**
     * This method should provide the functionality
     * to deinit the rendering procress.
     * You can deinit display lists here for example.
    */
    virtual void deinit() = 0;

    /**
     * This method should provide the functionality
     * to render this object.
    */
    virtual void render() = 0;

    /**
    *   Set visibility of the object.
    *   Invisible objects are not rendered.
    *   @param v Visibilty of the object
    */
    virtual void setVisible(bool v) {
        visible_ = v;
    }

    /**
    *   Set visibility of the object.
    *   Invisible objects are not rendered.
    *   @return Visibilty of the object
    */
    virtual bool getVisible() const {
        return visible_;
    }

    /**
    *   Static objects are not expected to move
    *   and groups of them can be cached in
    *   displaylists or vertexarrays.
    */
    virtual void setStatic(bool _static) {
        static_ = _static;
    }

    /**
    *   Static objects are not expected to move
    *   and groups of them can be cached in
    *   displaylists or vertexarrays.
    */
    virtual bool getStatic() const {
        return static_;
    }
};

}

#endif //TGT_RENDERABLE_H
