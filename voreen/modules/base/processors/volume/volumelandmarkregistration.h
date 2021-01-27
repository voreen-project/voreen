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

#ifndef VRN_VOLUMELANDMARKREGISTRATION_H
#define VRN_VOLUMELANDMARKREGISTRATION_H

#include <string>
#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/matrixproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "tgt/vector.h"

namespace voreen {

class Volume;
class Volume;

class VRN_CORE_API VolumeLandmarkRegistration : public VolumeProcessor {
public:
    VolumeLandmarkRegistration();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "VolumeLandmarkRegistration"; }
    virtual std::string getCategory() const  { return "Volume Processing";          }
    virtual CodeState getCodeState() const   { return CODE_STATE_TESTING;           }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Registers a volume to a reference volume based on landmarks. The landmarks are expected in voxel coordinates and are either provided as PointListGeometries through the geometry inports (see VolumePicking) or entered via properties. The computed registration matrix is written to the output volume and additionally provided through a matrix property.");
    }

    virtual void process();
    virtual void initialize();

private:
    tgt::mat4 computeRigidBodyTransformation(const VolumeBase* inputHandle, const VolumeBase* refHandle,
        tgt::vec3 p0, tgt::vec3 p0ref, tgt::vec3 p1, tgt::vec3 p1ref, tgt::vec3 p2, tgt::vec3 p2ref);

    tgt::mat4 computeAffineTransformation(const VolumeBase* inputHandle, const VolumeBase* refHandle,
        std::vector<tgt::vec3> inputPoints, std::vector<tgt::vec3> refPoints);

    tgt::vec3 transformFromVoxelToVolumeCoords(tgt::vec3 point, const VolumeBase* vol);

    void updateGuiState();

    VolumePort inportVolume_;
    VolumePort inportReference_;
    GeometryPort inportLandmarks_;
    GeometryPort inportReferenceLandmarks_;
    VolumePort outport_;

    BoolProperty enableProcessing_;
    StringOptionProperty registrationMode_;
    StringOptionProperty landmarkSpecification_;
    FloatVec3Property p1_;
    FloatVec3Property ref_p1_;
    FloatVec3Property p2_;
    FloatVec3Property ref_p2_;
    FloatVec3Property p3_;
    FloatVec3Property ref_p3_;
    FloatVec3Property p4_;
    FloatVec3Property ref_p4_;
    FloatMat4Property registrationMatrix_;

    static const std::string loggerCat_;
};

} //namespace

#endif // VRN_VOLUMELANDMARKREGISTRATION_H
