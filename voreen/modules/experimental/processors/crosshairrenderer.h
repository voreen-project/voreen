/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_CROSSHAIRRENDERER_H
#define VRN_CROSSHAIRRENDERER_H

#include "voreen/core/processors/imageprocessor.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "tgt/texturemanager.h"
#include "tgt/gpucapabilities.h"

namespace voreen {

class CrosshairRenderer : public ImageProcessor {
public:
    CrosshairRenderer();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "CrosshairRenderer";     }
    virtual std::string getCategory() const   { return "Image Processing"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL;  }

protected:
    virtual void setDescriptions() {
        setDescription("");
    }

    virtual void process();

    virtual void initialize();
    virtual void deinitialize();

    /**
     * Get screen position from 3D vector (project),
     * using the camera property's view and projection matrices.
     *
     * @param pos the point to project
     * @param modelview the modelview matrix to apply in addition to the camera's matrix.
     *      The camera's matrix is multiplied by the passed one, i.e., the passed one is applied first.
     * @param projection the projection matrix to apply in addition to the camera's matrix.
     *      The camera's matrix is multiplied by the passed one, i.e., the passed one is applied first.
     */
    tgt::vec3 getWindowPos(tgt::vec3 pos) const;
private:
    RenderPort inport_;
    RenderPort outport_;

    BoolProperty render_;
    FloatProperty width_;
    FloatProperty opacity_;
    ColorProperty color_;

    FloatVec3Property position_;
    CameraProperty camera_;

    tgt::Shader* copyShader_;
};

} // namespace

#endif
