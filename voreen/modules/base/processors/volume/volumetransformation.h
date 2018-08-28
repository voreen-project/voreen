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

#ifndef VRN_VOLUMETRANSFORMATION_H
#define VRN_VOLUMETRANSFORMATION_H

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/matrixproperty.h"

#include <string>

namespace voreen {

class Volume;

class VRN_CORE_API VolumeTransformation : public VolumeProcessor {
public:
    VolumeTransformation();
    virtual ~VolumeTransformation();
    virtual Processor* create() const;

    virtual std::string getCategory() const   { return "Volume Processing";    }
    virtual std::string getClassName() const  { return "VolumeTransformation"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_STABLE;      }

protected:
    virtual void setDescriptions() {
        setDescription("Transforms the input volume's position and orientation by modifying its transformation matrix.");
        sourceCoordinateSystem_.setDescription("Coordinates system from which the specified transformation should be applied (only in replace mode)");
        transformationMode_.setDescription("Concatenate: the volume's current transformation matrix is multiplied by the specified transformation matrix\
Replace: the volume's transformation matrix is replaced by the specified one");
        transformMatrix_.setDescription("The matrix that is applied for transforming the volume.");
    }

    virtual void process();
    virtual void initialize();

private:
    void adaptToTransformationMode();
    void adaptToTransformationType();

    void setTransformationMatrix(const tgt::mat4& mat);
    void setRotationMatrix();
    void setTranslationMatrix();

    VolumePort inport_;
    VolumePort outport_;

    BoolProperty enableProcessing_;
    StringOptionProperty transformationMode_;
    StringOptionProperty sourceCoordinateSystem_;

    enum TransformationType {
        ROTATION,
        TRANSLATION,
        ARBITRARY,
    };
    /// The selected type of transformation (see above)
    OptionProperty<TransformationType> transformationType_;
    /*
     * The transformation matrix itself.
     * Do not manipulate directly:
     * Use setTransformationMatrix() instead.
     */
    FloatMat4Property transformMatrix_;
    /*
     * Is used to distinguish between internal (e.g. via changing
     * of rotation properties) and external changes.
     */
    bool changingTransformationInternally_;

    /// Properties used to specify a rotation as a transformation.
    FloatProperty rotationAngle_;
    FloatVec3Property rotationAxis_;
    FloatVec3Property rotationReferencePoint_;

    /// Properties used to specify a translation as a transformation.
    FloatVec3Property translationAmount_;

    static const std::string loggerCat_; ///< category used in logging
};

}   //namespace

#endif // VRN_VOLUMETRANSFORMATION_H
