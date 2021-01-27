/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_IMAGEPROCESSOR_H
#define VRN_IMAGEPROCESSOR_H

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/shaderproperty.h"

#include "voreen/core/processors/renderprocessor.h"

namespace voreen {

/**
 * Base class for fragment shader based post processors.
 *
 * Normally in a derived class you have only to call the ctor with the
 * appropriate shader name and overwrite the render method.
 */
class VRN_CORE_API ImageProcessor : public RenderProcessor {
public:
    /**
     * Constructor.
     *
     * @param shaderFilename The filename of the shader that will be used.
     */
    ImageProcessor(const std::string& shaderFilename = "", bool enableCoarsenessAdaptation = false);
    virtual ~ImageProcessor();

    virtual std::string getCategory() const  { return "Image Processing"; }

    virtual void beforeProcess();

protected:
    virtual void initialize();
    //virtual void deinitialize();

    /// Load the needed shader.
    virtual void compile();

    virtual tgt::vec2 computeDepthRange(RenderPort* port);

    ShaderProperty shaderProp_;
    BoolProperty interactionAdapt_;
    IntProperty interactionFactor_;

    tgt::Shader* program_;
    std::string shaderFilename_;

    static const std::string loggerCat_; ///< category used in logging
};

} // namespace

#endif // VRN_IMAGEPROCESSOR_H
