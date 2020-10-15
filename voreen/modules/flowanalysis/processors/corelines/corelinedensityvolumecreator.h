#pragma once

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen
{
	class CorelineDensityVolumeCreator : public Processor
	{
	public:
		CorelineDensityVolumeCreator();
		virtual Processor *create() const { return new CorelineDensityVolumeCreator(); }
		virtual std::string getClassName() const { return "CorelineDensityVolumeCreator"; }
		virtual std::string getCategory() const { return "Volume Processing"; }

		static void Process(const std::vector<std::vector<tgt::vec3>>& corelines, VolumeRAM_Float& outBinaryVolume);

	protected:
		virtual void process();

	private:
		GeometryPort _inCorelines;
		VolumePort _inVolume;
		VolumePort _out;
	};
} // namespace voreen
