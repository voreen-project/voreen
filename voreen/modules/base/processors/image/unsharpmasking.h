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

#ifndef VRN_UNSHARPMASKING_H
#define VRN_UNSHARPMASKING_H

#include "voreen/core/processors/imageprocessorbypassable.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"

namespace voreen {

/**
 * Performs an unsharp masking operation by subtracting a blurred version from the input image
 */
class VRN_CORE_API UnsharpMasking : public ImageProcessorBypassable {
public:
    UnsharpMasking();
    virtual Processor* create() const { return new UnsharpMasking(); }

    virtual std::string getCategory() const { return "Image Processing"; }
    virtual std::string getClassName() const { return "UnsharpMasking"; }
    virtual CodeState getCodeState() const { return CODE_STATE_STABLE; }

protected:
    virtual void setDescriptions() {
        setDescription("Performs an unsharped masking operation by subtracting a (Gaussian) blurred version of the image.");
        alpha_.setDescription("Controls, if the unsharp masking is applied to alpha values as well.");
        amount_.setDescription("Amount controls the amount contrast is increased.");
        radius_.setDescription("Radius controls the radius of the blur, which is co-dependent with Sigma.");
        sigma_.setDescription("Sigma is the parameter for the Gaussian blur applied to get an unsharp version of the image, which is later subtracted. \
\
Co-dependent with Radius");
        threshold_.setDescription("By using Threshold, masking is only applied when a minimal value is reached.");
    }

    void process();
    virtual void initialize();
    virtual void deinitialize();

    //methods for synchronizing radius and sigma slider
    void radiusChanged();
    void sigmaChanged();

    IntProperty radius_;
    FloatProperty sigma_;
    FloatProperty amount_;
    FloatProperty threshold_;
    BoolProperty alpha_;

    RenderPort inport_;
    RenderPort outport_;

    RenderPort privatePort_;

    tgt::Shader* secondProgram_; //used for the actual masking after blurring

    //flags for the linking of radius and sigma
    bool sigmaFlag_;
    bool radiusFlag_;
};


} // namespace voreen

#endif //VRN_UNSHARPMASKING_H
