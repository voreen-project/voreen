/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2024 University of Muenster, Germany,           *
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

#ifndef TGT_QTCANVAS_H
#define TGT_QTCANVAS_H

#include "tgt/glcanvas.h"
#include "voreen/qt/voreenqtapi.h"

#ifdef __APPLE__
#define __gl3_h_
#endif

#include <QWidget>
#include <QSurfaceFormat>

#include <memory>

class QInputEvent;

namespace tgt {

class CanvasBackend;
class CanvasBackendOwner {
public:
    virtual ~CanvasBackendOwner() {}
    virtual void onInitializeGL() = 0;
    virtual void onPaintGL() = 0;
    virtual void onResizeGL(int w, int h) = 0;
};

/**
 * Qt implementation of GLCanvas. Inherits QOpenGLWidget and combines the Qt methods and tgt methods.
 */
class VRN_QT_API QtCanvas : public QWidget, public GLCanvas, public CanvasBackendOwner {
public:
    /**
     * The constructor.
     *
     * @param parent The parent widget of this canvas.
     * @param f Qt::Wflags can be passed to this constructor to control the qt features, like stereo-buffering.
     */
    QtCanvas(const std::string& title = "",
             const ivec2& size = ivec2(DEFAULT_WINDOW_WIDTH, DEFAULT_WINDOW_HEIGHT),
             const Buffers buffers = RGBADD,
             QWidget* parent = 0,
             Qt::WindowFlags f = Qt::WindowFlags());

    /**
     * Destructor. Closes window (if canvas is a window).
     */
    virtual ~QtCanvas();

    /**
     * @return an image of the current framebuffer state.
     */
    QImage grabFramebuffer();

public:
    // Inherited from GLCanvas

    /**
     * If you manually want to cause a paint-event, use this function.
     * This will cause immediate repainting.
     */
    virtual void repaint();

    /**
     * If you manually want to cause a paint-event, use this function. It will call QWidget::update()
     * and repaint when entering main loop next time.
     */
    virtual void update();

    virtual void swap();
    virtual void toggleFullScreen();
    virtual ivec2 getPhysicalSize() const;

protected:

    // Inherited from CanvasBackendOwner
    virtual void onInitializeGL();
    virtual void onPaintGL();
    virtual void onResizeGL(int w, int h);

protected:

    // Inherited from QWidget
    virtual void enterEvent(QEvent* e);
    virtual void leaveEvent(QEvent* e);
    virtual void mousePressEvent(QMouseEvent* e);
    virtual void mouseReleaseEvent(QMouseEvent* e);
    virtual void mouseMoveEvent(QMouseEvent*  e);
    virtual void mouseDoubleClickEvent(QMouseEvent* e);
    virtual void wheelEvent(QWheelEvent* e);
    virtual void timerEvent(QTimerEvent* e);
    virtual void keyPressEvent(QKeyEvent* event);
    virtual void keyReleaseEvent(QKeyEvent* event);

    // necessary to handle touch events
    virtual bool event(QEvent *event);

    void broadcastEvent(tgt::Event* event);

protected:
    // Inherited from GLContextBase
    virtual void activate();
    virtual bool isActive();

private:

    std::unique_ptr<CanvasBackend> backend_;
    bool initializedGL_;

public:
    ///
    /// Helpers used to generate tgt-Events out of qt-Events
    ///

    // map one Qt-mousebutton to one tgt-mousebutton
    static tgt::MouseEvent::MouseButtons getButton(QMouseEvent* e);
    // map a set of Qt-mousebuttons to a set of tgt-mousebuttons
    static tgt::MouseEvent::MouseButtons getButtons(QMouseEvent* e);

    static tgt::Event::Modifier getModifier(QInputEvent* e);
    static KeyEvent::KeyCode getKey(int key);
    static QSurfaceFormat getSurfaceFormat(const Buffers buffers);
};

} // namespace

#endif // TGT_QTCANVAS_H
