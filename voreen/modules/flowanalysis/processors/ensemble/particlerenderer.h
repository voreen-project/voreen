#pragma once
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

namespace voreen
{
	class ParticleRenderer : public RenderProcessor
	{
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