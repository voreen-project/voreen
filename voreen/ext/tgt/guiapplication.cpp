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

#include "tgt/guiapplication.h"
#include "tgt/tgt_gl.h"

namespace tgt {

GUIApplication::GUIApplication( int& argc, char** argv ) {
    argc_ = &argc;
    argv_ = argv;
    hasRun = false;
}

GUIApplication::~GUIApplication() {
}

int GUIApplication::init(InitFeature::Features featureset) {
    tgt::init(featureset);
    initCanvases();
    tgt::initGL(featureset);
    return 0;
}

void GUIApplication::addCanvas(GLCanvas* canvas) {
    canvases_.push_back(canvas);
}

void GUIApplication::removeCanvas(GLCanvas* canvas) {
    canvases_.remove(canvas);
}

void GUIApplication::initCanvases() {
    std::list<GLCanvas*>::iterator iter;
    for ( iter = canvases_.begin(); iter != canvases_.end(); iter++ ) {
        if (!(*iter)->isInitialized()) (*iter)->init();
    }
}

} // namespace tgt
