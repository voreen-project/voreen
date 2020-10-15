#pragma once
#ifndef VRN_VORTEXPROCESSOR_H
#define VRN_VORTEXPROCESSOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen
{
	class VortexProcessor : public Processor
	{
	public:
		VortexProcessor();

		Processor* create() const override
		{
			return new VortexProcessor();
		}
		std::string getClassName() const override
		{
			return "VortexProcessor";
		}
		std::string getCategory() const override
		{
			return "Vortex Extraction";
		}
		bool isReady() const override
		{
			return _inportVelocity.isReady() && _inportMask.isReady();
		}

		static void Process( const VolumeRAM& velocity, const VolumeRAM& mask, VolumeRAM_Mat3Float& outJacobi, VolumeRAM_Float* outDelta, VolumeRAM_Float* outLambda2, VolumeRAM_Float* outQ );

	private:
		void process() override;

		VolumePort _inportVelocity, _inportMask;
		VolumePort _outportJacobi, _outportDelta, _outportLambda2, _outportQ;
	};
}

#endif // VRN_VORTEXPROCESSOR_H