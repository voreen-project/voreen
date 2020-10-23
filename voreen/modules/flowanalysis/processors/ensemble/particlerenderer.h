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

#ifndef VRN_PARTICLERENDERER_H
#define VRN_PARTICLERENDERER_H

#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include "voreen/core/interaction/camerainteractionhandler.h"

#include "modules/ensembleanalysis/ports/ensembledatasetport.h"

#include <memory>

namespace voreen {
class ParticleRenderer : public RenderProcessor {
public:
    ParticleRenderer();

    Processor* create() const override;
    std::string getCategory() const override;
    std::string getClassName() const override;

private:
    void initialize() override;
    void deinitialize() override;
    void process() override;

    void updateMesh();

    EnsembleDatasetPort _inportEnsemble;
    RenderPort _outportImage;

    IntOptionProperty _propertySelectedMember;
    IntOptionProperty _propertyFlowComponentX, _propertyFlowComponentY, _propertyFlowComponentZ;
    IntOptionProperty _propertySelectedField;
    IntProperty _propertySeedPoints;
    ButtonProperty _propertyUpdateMesh;

    TransFunc1DKeysProperty _propertyTransferFunction;
    ButtonProperty _propertyResetDomain;
    FloatProperty _propertyLineWidth;
    IntIntervalProperty _propertyTimesteps;
    ShaderProperty _propertyShader;
    CameraProperty _propertyCamera;

    std::unique_ptr<CameraInteractionHandler> _interactionHandlerCamera;

    tgt::vec2 _domain;
    std::unique_ptr<GlMeshGeometryBase> _mesh;
};
}

#endif // VRN_PARTICLERENDERER_H