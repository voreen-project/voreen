/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2018 Visualization and Computer Graphics Group, *
 * Department of Computer Science, University of Muenster, Germany.   *
 * <http://viscg.uni-muenster.de>                                     *
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

#ifndef TGT_GUIAPPLICATION_H
#define TGT_GUIAPPLICATION_H

#include <list>

#include "tgt/init.h"
#include "tgt/vector.h"
#include "tgt/glcanvas.h"
#include "tgt/timer.h"

namespace tgt {

class GUIApplication {

protected:
    /// used to hold the canvases registered in our app.
    std::list<GLCanvas*> canvases_;

    // to store command line argumentes
    int* argc_;
    char** argv_;

    /// walks through canvas list to call init-method of each canvas
    void initCanvases();

    bool hasRun;

public:
    /// Constructor
    /// @param argc number of command line arguments
    /// @param vector of command line arguments
    GUIApplication(int& argc, char** argv);

    /// Destructor
    virtual ~GUIApplication();

    /// register an existing Canvas at the application
    virtual void addCanvas(GLCanvas* canvas);

    /// unregister Canvas
    virtual void removeCanvas( GLCanvas* );

    /// Initialize GUI toolkit, canvases and some more tgt stuff.
    /// At least one canvas must exist (must have been constructed) when calling this method.
    /// This method should be called as early as possible in programm.
    /// There should be no canvases initialized before calling this method as this might
    /// cause problems when using certain tgt features on this canvas.
    /// Method will initialize all canvases that exist at the time of method call.
    /// No OpenGL context shoud exist before calling this method, at least one will exist afterwards.
    virtual int init(InitFeature::Features featureset = InitFeature::ALL);

    /// enter GUI toolkit's event handling mainloop
    virtual int run() = 0;

    /// to leave GUI toolkit's mainloop
    virtual void quit(int status = 0) = 0;
};

} // namespace tgt

#endif // TGT_GUIAPPLICATION_H
