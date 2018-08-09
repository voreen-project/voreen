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

#ifndef TGT_GLUTAPPLICATION_H
#define TGT_GLUTAPPLICATION_H

#include "tgt/guiapplication.h"
#include "tgt/glut/glutcanvas.h"
#include "tgt/glut/gluttimer.h"

#include <GL/glut.h>

namespace tgt {

class GLUTApplication : public GUIApplication {
public:
    /// Constructor
    /// @param argc number of command line arguments
    /// @param vector of command line arguments
    GLUTApplication(int& argc, char** argv);

    /**
     * Destructor.
     *
     * Deletes all GLCanvases and deinits tgt. This will never be called in a normal way
     * because of GLUT's restrictions
     * @see onExitFunction
     */
    virtual ~GLUTApplication();

    /// initialize tgt and glut, create a window for every GLUTCanvas and stuff ...
    virtual int init(InitFeature::Features featureset = InitFeature::ALL);

    /// enter glut mainloop (glutMainLoop)
    virtual int run();

    /// tells glut to quit mainloop and exit program
    virtual void quit(int status = 0);

}; // class GLUTApplication

} // namespace tgt

#endif // TGT_GLUTAPPLICATION_H
