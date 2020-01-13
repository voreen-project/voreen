/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_NORMALESTIMATION_H
#define VRN_NORMALESTIMATION_H

#include "voreen/core/processors/imageprocessor.h"
#include "voreen/core/properties/cameraproperty.h"

namespace voreen {

class CameraProperty;

/**
 * Image-based normal estimation used for deformation
 */
class NormalEstimation : public ImageProcessor {
public:
    NormalEstimation();
    virtual std::string getCategory() const { return "Image Processing"; }
    virtual std::string getClassName() const { return "NormalEstimation"; }
    virtual Processor* create() const;
    void process();

protected:
    virtual void setDescriptions() {
        setDescription("Image-based normal estimation used for deformation.");
    }

    CameraProperty camera_;

    RenderPort inport_;
    RenderPort outport_;
};

} // namespace

#endif // VRN_NORMALESTIMATION_H
