/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
 * Department of Computer Science.                                                 *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

#ifndef TGT_DummyCanvas_H
#define TGT_DummyCanvas_H

#include "tgt/glcanvas.h"
#include "tgt/event/keyevent.h"

namespace voreen {

class VtkDummyCanvas : public tgt::GLCanvas {
public:

    VtkDummyCanvas(const std::string& title = "",
        const tgt::ivec2& size = tgt::ivec2(DEFAULT_WINDOW_WIDTH, DEFAULT_WINDOW_HEIGHT),
        const GLCanvas::Buffers buffers = RGBADD );

    ~VtkDummyCanvas();

    void init();

    void swap();

    void getGLFocus();

    void toggleFullScreen();
 
    bool isFullScreen();
    
    void repaint();
 
    void update();

    inline int getWindowID();
};

}
#endif // TGT_DummyCanvas_H
