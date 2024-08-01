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

#ifndef VRN_CLASSIFICATIONMODES_H
#define VRN_CLASSIFICATIONMODES_H

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/transfunc/transfuncbase.h"
#include "tgt/shadermanager.h"
#include "tgt/camera.h"
#include "tgt/textureunit.h"

namespace voreen {

class StringOptionProperty;

class VRN_CORE_API ClassificationModes {
public:
    static std::string getShaderDefineSamplerType(const std::string mode, const TransFuncBase* tf, const std::string& defineName = "TF_SAMPLER_TYPE");
    static std::string getShaderDefineFunction(const std::string mode, const std::string& defineName = "RC_APPLY_CLASSIFICATION");
    static void bindTexture(const std::string mode, TransFuncBase* tf, float samplingStepSize);
    static bool usesTransferFunction(const std::string mode);
    static void fillProperty(StringOptionProperty* prop);
};

}  // namespace voreen

#endif
