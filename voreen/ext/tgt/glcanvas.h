/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2018 University of Muenster, Germany,           *
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

#ifndef TGT_GLCANVAS_H
#define TGT_GLCANVAS_H

#include "tgt/types.h"
#include "tgt/vector.h"
#include "tgt/event/event.h"
#include "tgt/event/timeevent.h"
#include "tgt/event/mouseevent.h"
#include "tgt/event/keyevent.h"
#include "tgt/event/eventhandler.h"
#include "tgt/glcontextbase.h"

namespace tgt {

class Camera;
class Painter;

/**
 * This class is the base class for all tgt-Canvases. It provides the functionality of using
 * Painter-objects to actually render things to the canvas. The methods one has to implement in a
 * subclass are, most importantly, those taking care of turning incoming events into tgt-events.
 * This makes sure that all kinds of APIs can be used with tgt without specializing
 * a lot of code.
 */
class TGT_API GLCanvas : public GLContextBase {

    friend class Painter;

public:
    enum {
        DEFAULT_WINDOW_WIDTH  = 512,
        DEFAULT_WINDOW_HEIGHT = 512
    };

    enum Buffers {
        RGB_BUFFER          = 1 << 0,
        ALPHA_BUFFER        = 1 << 1,
        DEPTH_BUFFER        = 1 << 2,
        DOUBLE_BUFFER       = 1 << 3,
        STENCIL_BUFFER      = 1 << 4,
        ACCUM_BUFFER        = 1 << 5,
        QUAD_BUFFER         = 1 << 6, //activates stereo viewing in qt. must be used with double_buffer
        MULTISAMPLING       = 1 << 7,
        RGBA_BUFFER         = RGB_BUFFER | ALPHA_BUFFER,

        // frequently used settings
        RGBD    = RGB_BUFFER | DEPTH_BUFFER,
        RGBDS   = RGB_BUFFER | DEPTH_BUFFER | STENCIL_BUFFER,
        RGBDD   = RGB_BUFFER | DEPTH_BUFFER | DOUBLE_BUFFER,
        RGBDDS  = RGB_BUFFER | DEPTH_BUFFER | DOUBLE_BUFFER | STENCIL_BUFFER,
        RGBA    = RGB_BUFFER | ALPHA_BUFFER,
        RGBAD   = RGB_BUFFER | ALPHA_BUFFER | DEPTH_BUFFER,
        RGBADS  = RGB_BUFFER | ALPHA_BUFFER | DEPTH_BUFFER | STENCIL_BUFFER,
        RGBADD  = RGB_BUFFER | ALPHA_BUFFER | DEPTH_BUFFER | DOUBLE_BUFFER,
        RGBADDS = RGB_BUFFER | ALPHA_BUFFER | DEPTH_BUFFER | DOUBLE_BUFFER | STENCIL_BUFFER,
        RGBADDQ  = RGB_BUFFER | ALPHA_BUFFER | DEPTH_BUFFER | DOUBLE_BUFFER | QUAD_BUFFER,
        RGBADDQS = RGB_BUFFER | ALPHA_BUFFER | DEPTH_BUFFER | DOUBLE_BUFFER | QUAD_BUFFER | STENCIL_BUFFER
    };

    /// A Constructor
    /// @param title     window title if canvas is standalone window
    /// @param size      size of canvas in pixels
    /// @param buffers   which buffer setting to use, default is to use double-buffered
    ///                  RGB buffer with alpha and depth buffer
    ///                  @see GLCanvas::Buffers
    GLCanvas(const std::string& title = "",
             const ivec2& size = ivec2(DEFAULT_WINDOW_WIDTH, DEFAULT_WINDOW_HEIGHT),
             const Buffers buffers = RGBADD );

    /**
     * The GLCanvas and the Painter which owns the GLCanvas should be destroyed
     * manually otherwise endless recursions can occur.
     */
    virtual ~GLCanvas();

public:
    /// This is the method that swaps the front- and backbuffer. To be overridden by derived
    /// Canvas subclasses.
    virtual void swap() = 0;

    /**
     * Toggles fullscreen mode.
     *
     * To be overridden by derived Canvas classes.
     */
    virtual void toggleFullScreen() = 0;

    /**
     * Should be called when the canvas was resized.
     * To be called by the according toolkit (Qt, GLUT, SDL, ...).
     */
    void sizeChanged(const ivec2& size);

    /// to be called by application to cause (re)painting on the canvas
    virtual void repaint() = 0;

    /// to be called by application to cause the canvas to repaint when idle next time
    virtual void update() = 0;

    /// initialize the canvas -- e.g. create window, set it's size and title
    virtual void init();

    EventHandler* getEventHandler() const;
    void setEventHandler(EventHandler* handler);

    /// make the canvas call glFlush/swap automatically or not
    void setAutoFlush(bool a) { autoFlush_ = a; };

    /// wheater canvas automatically calls glFlush or swap
    bool getAutoFlush();

    /// force screen update using either swap() or glFlush()
    void finalizeRendering();

    ///Take a screenshot and save it as TGA file.
    ///@param fname Target filename
    bool takeScreenshot(std::string fname);

    const ivec4& getRgbaSize() const;
    int getDepthSize() const;
    int getStencilSize() const;
    ivec2 getSize() const;
    int getWidth() const;
    int getHeight() const;
    virtual ivec2 getPhysicalSize() const; // Size in phyical pixels on screen. Important for Hidpi displays
    int getPhysicalWidth() const;
    int getPhysicalHeight() const;
    Buffers getBuffers() const;
    bool isDoubleBuffered() const { return doubleBuffered_; }
    bool isQuadBuffered() const { return quadBuffered_; }
    bool isInitialized() const { return initialized_; }

protected:
    /// Set the painter the Canvas will use to draw it's content. Called from painter
    /// @param p the Painter
    void setPainter(Painter* p);

    /// Use the painter_ to actually paint something on the canvas
    /// For internal use, not to be called by application directly. Call repaint instead.
    virtual void paint();

    std::string title_; /// window title if canvas is represented by a window
    ivec2   size_;    /// the size of the canvas
    Buffers buffers_; /// the kind of buffers that will be used

    ivec4   rgbaSize_;
    int     depthSize_;
    int     stencilSize_;
    bool    doubleBuffered_;
    bool    quadBuffered_;
    bool    fullscreen_;
    bool    autoFlush_; ///< whether to call glFlush or swap automatically
    bool    initialized_;

    Painter* painter_;  ///< the painter that will be used for rendering

    EventHandler* eventHandler_;  ///< the eventHandler that will distribute incoming events
};

}

#endif // TGT_GLCANVAS_H
