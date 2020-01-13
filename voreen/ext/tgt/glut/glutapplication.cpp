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

#include <list>

#include "tgt/glut/glutapplication.h"
#include "tgt/glut/tgt_glut.h"
#include "tgt/glut/glutcanvas.h"

namespace tgt {

GLUTApplication::GLUTApplication(int& argc, char** argv)
  : GUIApplication(argc, argv)
{}


GLUTApplication::~GLUTApplication(){
    deinitGL();
    deinit();
}

int GLUTApplication::init(InitFeature::Features featureset) {
    tgt::init(featureset);
    glutInit(argc_, argv_);
    initCanvases();
    tgt::initGL(featureset);
    return 0;
}

int GLUTApplication::run() {
    if (hasRun) {
        std::cout << "Failure: GLUTApplication::run() should be called only once.\n";
        return -1;
    }
    hasRun = true;
    glutMainLoop();
    return 0;
}

void GLUTApplication::quit(int status) {
    exit(status); // FIXME: exiting a glut application is not done correct. no good concept how to use destructors, atexit, canvas-application, and so on. glut will never return form mainloop. one has to use atexit, which differs from usage of other toolkits in factory. freeglut offers the ability to return from mainloop. perhaps we should limit tgt to use freeglut only, no glut...
}

}
