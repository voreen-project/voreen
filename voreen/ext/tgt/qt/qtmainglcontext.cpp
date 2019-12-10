/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2019 University of Muenster, Germany,           *
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

#include <QOpenGLContext>
#include <QOffscreenSurface>
#include <QSurfaceFormat>

#include "qtmainglcontext.h"
#include "qtcanvas.h"

//#define VRN_DEBUG_QT_MAIN_GL_CONTEXT

#ifdef VRN_DEBUG_QT_MAIN_GL_CONTEXT
#include <QOpenGLDebugLogger>

class LogHandler : public QObject {
public:

    QOpenGLDebugLogger logger_;

public slots:
    void onMessageLogged(QOpenGLDebugMessage message) {
        //Set breakpoint here
        qDebug() << message;
    }

} handler;

#endif

namespace tgt {

QtMainGLContext::QtMainGLContext() 
    : GLContextBase("main")
{
    QSurfaceFormat format = QtCanvas::getSurfaceFormat(tgt::GLCanvas::RGBADD);
#ifdef VRN_DEBUG_QT_MAIN_GL_CONTEXT
    format.setOption(QSurfaceFormat::DebugContext);
#endif

    context_ = new QOpenGLContext();
    context_->setFormat(format);
    context_->setShareContext(QOpenGLContext::globalShareContext());
    context_->create();

    surface_ = new QOffscreenSurface();
    surface_->setFormat(format);
    surface_->create();

    // Activate directly because further initialization needs an active gl context.
    activate();

#ifdef VRN_DEBUG_QT_MAIN_GL_CONTEXT
    // Attach debugger.
    handler.logger_.initialize();
    QObject::connect(&handler.logger_, &QOpenGLDebugLogger::messageLogged, &handler, &LogHandler::onMessageLogged);
    handler.logger_.startLogging(QOpenGLDebugLogger::SynchronousLogging);
#endif
}
QtMainGLContext::~QtMainGLContext() {
    // Release active context before deletion.
    context_->doneCurrent();

    delete context_;
    delete surface_;
}

void QtMainGLContext::activate() {
    context_->makeCurrent(surface_);
}

bool QtMainGLContext::isActive() {
    return QOpenGLContext::currentContext() == context_;
}

} // namespace
