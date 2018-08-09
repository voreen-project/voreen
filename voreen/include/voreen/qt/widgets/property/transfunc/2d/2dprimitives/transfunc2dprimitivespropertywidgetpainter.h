/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_TRANSFUNC2DPRIMITIVESPROPERTYWIDGETPAINTER_H
#define VRN_TRANSFUNC2DPRIMITIVESPROPERTYWIDGETPAINTER_H

#include "voreen/qt/widgets/property/transfunc/2d/transfunc2dpropertywidgetpainter.h"

namespace tgt {
    class GLCanvas;
}

#include <QObject>

namespace voreen {

/**
 * Painter implementation for drawing onto a tgt::QtCanvas for editing 2D transfer functions.
 * There are 2 types of primitives that can be used to specify the transfer function, a banana
 * and a quad. Every primitive can be colored and moved. Furthermore you can change the extent
 *  and shape of the primitive by dragging the control points.
 */
class TransFunc2DPrimitivesPropertyWidgetPainter : public TransFunc2DPropertyWidgetPainter {
    Q_OBJECT
public:
    /**
     * Constructor
     *
     * @param canvas the canvas that uses this painter for painting
     */
    TransFunc2DPrimitivesPropertyWidgetPainter(tgt::GLCanvas* canvas, QColor backgroundColor);

    /**
     * Destructor
     */
    ~TransFunc2DPrimitivesPropertyWidgetPainter();

};

} // namespace voreen

#endif // VRN_TRANSFUNC2DPRIMITIVESPROPERTYWIDGETPAINTER_H
