/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_TRANSFUNC1DGAUSSIANPROPERTYWIDGETPAINTER_H
#define VRN_TRANSFUNC1DGAUSSIANPROPERTYWIDGETPAINTER_H

#include "voreen/qt/widgets/property/transfunc/1d/transfunc1dpropertywidgetpainter.h"

namespace voreen {

/**
 * The painter used to draw the tf texture into the associated texture canvas of the tf property widget.
 * @note this class has no functionality yet. 1dgaussian tfs use the default tf 1d painter.
 */
class TransFunc1DGaussianPropertyWidgetPainter : public TransFunc1DPropertyWidgetPainter {
public:
    /** Constructor */
    TransFunc1DGaussianPropertyWidgetPainter(tgt::GLCanvas* canvas);
    /** Destructor */
    virtual ~TransFunc1DGaussianPropertyWidgetPainter();
};

} // namespace voreen

#endif // VRN_TRANSFUNC1DGAUSSIANPROPERTYWIDGETPAINTER_H
